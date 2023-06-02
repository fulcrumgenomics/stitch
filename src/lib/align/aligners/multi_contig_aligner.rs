use crate::align::aligners::constants::DEFAULT_ALIGNER_CAPACITY;
use crate::align::aligners::single_contig_aligner::SingleContigAligner;
use crate::align::alignment::Alignment;
use crate::align::scoring::Scoring;
use crate::align::traceback::traceback;
use bio::alignment::pairwise::MatchFunc;
use bio::utils::TextSlice;
use itertools::Itertools;

struct ContigAligner<'a, F: MatchFunc> {
    pub name: String,
    pub is_forward: bool,
    pub aligner: SingleContigAligner<F>,
    pub seq: &'a [u8],
}

impl<'a, F: MatchFunc> ContigAligner<'a, F> {
    pub fn new(
        name: String,
        is_forward: bool,
        scoring: Scoring<F>,
        seq: TextSlice<'a>,
        contig_idx: usize,
        circular: bool,
    ) -> ContigAligner<'a, F> {
        let mut aligner = SingleContigAligner::with_capacity_and_scoring(
            DEFAULT_ALIGNER_CAPACITY,
            DEFAULT_ALIGNER_CAPACITY,
            scoring,
        );
        aligner.set_contig_idx(contig_idx);
        aligner.set_circular(circular);
        Self {
            name,
            is_forward,
            aligner,
            seq,
        }
    }

    pub fn len(&self) -> usize {
        self.seq.len()
    }
}

pub struct MultiContigAligner<'a, F: MatchFunc> {
    contigs: Vec<ContigAligner<'a, F>>,
}

impl<'a, F: MatchFunc> MultiContigAligner<'a, F> {
    /// Create new aligner instance with given scorer.
    ///
    /// # Arguments
    ///
    /// * `scoring_fn` - function that returns an alignment scorer
    ///    (see also [`bio::alignment::pairwise::Scoring`](struct.Scoring.html))
    pub fn new() -> Self {
        MultiContigAligner {
            contigs: Vec::new(),
        }
    }

    /// Adds a new aligner for the given contig and strand.
    pub fn add_contig(
        &mut self,
        name: &str,
        is_forward: bool,
        seq: TextSlice<'a>,
        circular: bool,
        scoring: Scoring<F>,
    ) {
        let alread_exists = self
            .contigs
            .iter()
            .any(|a| a.name == name && a.is_forward == is_forward);
        assert!(
            !alread_exists,
            "Contig already added! name: {name} is_forward: {is_forward}"
        );

        let contig_idx = self.contigs.len();
        let contig_aligner = ContigAligner::new(
            name.to_string(),
            is_forward,
            scoring,
            seq,
            contig_idx,
            circular,
        );
        self.contigs.push(contig_aligner);
    }

    /// The core function to compute the alignment
    ///
    /// # Arguments
    ///
    /// * `x` - Textslice
    /// * `y` - Textslice
    pub fn custom(&mut self, y: TextSlice<'_>) -> Alignment {
        let n = y.len();

        // Set the initial conditions
        // We are repeating some work, but that's okay!
        for contig in &mut self.contigs {
            contig.aligner.init_matrices(contig.len(), n);
        }

        for j in 1..=n {
            let curr = j % 2;
            let prev = 1 - curr;

            // Initialize the column
            for contig in &mut self.contigs {
                contig.aligner.init_column(j, curr, contig.len(), n);
            }

            // Initialize the jump buffers
            // TODO: jump score across contigs, and jump score within contigs.
            // TODO: jump score across the same strand vs. opposite strand of the same contig
            let best_jump_info = self
                .contigs
                .iter()
                .map(|contig| contig.aligner.get_jump_info(contig.len(), j - 1))
                .max_by_key(|contig| (contig.score, contig.len))
                .unwrap();

            // Fill the column
            for contig in &mut self.contigs {
                // Prefer a jump to the same contig
                let jump_info: crate::align::aligners::JumpInfo =
                    contig.aligner.get_jump_info(contig.len(), j - 1);
                let best_jump_info = if jump_info.score == best_jump_info.score
                    && jump_info.len == best_jump_info.len
                {
                    jump_info
                } else {
                    best_jump_info
                };

                contig.aligner.fill_column(
                    contig.seq,
                    y,
                    contig.len(),
                    n,
                    j,
                    prev,
                    curr,
                    best_jump_info,
                );
            }
        }

        for contig in &mut self.contigs {
            contig
                .aligner
                .fill_last_column_and_end_clipping(contig.len(), n);
        }

        let aligners = self
            .contigs
            .iter()
            .map(|contig| &contig.aligner)
            .collect_vec();
        traceback(&aligners, n)
    }
}

// Tests
#[cfg(test)]
pub mod tests {
    use bio::alignment::pairwise::MatchParams;
    use itertools::Itertools;
    use rstest::rstest;

    use crate::{
        align::{aligners::constants::MIN_SCORE, scoring::Scoring},
        util::dna::reverse_complement,
    };

    use super::{Alignment, MultiContigAligner};

    /// Upper-cases and remove display-related characters from a string.
    fn s(bases: &str) -> Vec<u8> {
        bases
            .chars()
            .filter(|base| *base != '-' && *base != ' ' && *base != '_')
            .map(|base| base.to_ascii_uppercase() as u8)
            .collect_vec()
    }

    fn assert_alignment(
        alignment: &Alignment,
        xstart: usize,
        xend: usize,
        ystart: usize,
        yend: usize,
        score: i32,
        contig_idx: usize,
        cigar: &str,
        length: usize,
    ) {
        assert_eq!(alignment.xstart, xstart, "xstart {alignment}");
        assert_eq!(alignment.xend, xend, "xend {alignment}");
        assert_eq!(alignment.ystart, ystart, "ystart {alignment}");
        assert_eq!(alignment.yend, yend, "yend {alignment}");
        assert_eq!(alignment.score, score, "score {alignment}");
        assert_eq!(alignment.contig_idx, contig_idx, "contig_idx {alignment}");
        assert_eq!(alignment.cigar(), cigar, "cigar {alignment}");
        assert_eq!(alignment.length, length, "length {alignment}");
    }

    fn scoring_global_custom(
        mismatch_score: i32,
        gap_open: i32,
        gap_extend: i32,
        jump_score: i32,
    ) -> Scoring<MatchParams> {
        let match_fn = MatchParams::new(1, mismatch_score);
        let mut scoring = Scoring::new(gap_open, gap_extend, jump_score, match_fn);
        scoring.xclip_prefix = MIN_SCORE;
        scoring.xclip_suffix = MIN_SCORE;
        scoring.yclip_prefix = MIN_SCORE;
        scoring.yclip_suffix = MIN_SCORE;
        scoring
    }

    fn scoring_global() -> Scoring<MatchParams> {
        scoring_global_custom(-1, -5, -1, -10)
    }

    /// Identical sequences, all matches
    #[rstest]
    fn test_identical() {
        let x = s("ACGTAACC");
        let x_revcomp = reverse_complement(&x);
        let y = s("ACGTAACC");
        let mut aligner = MultiContigAligner::new();
        aligner.add_contig("fwd", true, &x, false, scoring_global());
        aligner.add_contig("revcomp", false, &x_revcomp, false, scoring_global());
        let alignment = aligner.custom(&y);
        assert_alignment(&alignment, 0, 8, 0, 8, 8, 0, "8=", 8);
    }

    /// Identical sequences, all matches, reverse complemented
    #[rstest]
    fn test_identical_revcomp() {
        let x = s("ACGTAACC");
        let x_revcomp = reverse_complement(&x);
        let y = reverse_complement(s("ACGTAACC"));
        let mut aligner = MultiContigAligner::new();
        aligner.add_contig("fwd", true, &x, false, scoring_global());
        aligner.add_contig("revcomp", false, &x_revcomp, false, scoring_global());
        let alignment = aligner.custom(&y);
        assert_alignment(&alignment, 0, 8, 0, 8, 8, 1, "8=", 8);
    }

    #[rstest]
    fn test_fwd_to_fwd_jump() {
        let x = s("AAGGCCTT");
        let x_revcomp = reverse_complement(&x);
        let y = s("AACCGGTT");
        let mut aligner = MultiContigAligner::new();
        aligner.add_contig(
            "fwd",
            true,
            &x,
            false,
            scoring_global_custom(-1, -100_000, -100_000, -1),
        );
        aligner.add_contig(
            "revcomp",
            false,
            &x_revcomp,
            false,
            scoring_global_custom(-1, -100_000, -100_000, -1),
        );
        let alignment = aligner.custom(&y);
        assert_alignment(
            &alignment,
            0,
            8,
            0,
            8,
            8 - 1 - 1 - 1,
            0,
            "2=2J2=4j2=2J2=",
            8,
        );
    }

    #[rstest]
    fn test_fwd_to_rev_jump() {
        let x = s("AACCTTGG");
        let x_revcomp = reverse_complement(&x); // CCAAGGTT
        let y = s("AACCGGTT");
        let mut aligner = MultiContigAligner::new();
        aligner.add_contig(
            "fwd",
            true,
            &x,
            false,
            scoring_global_custom(-100_000, -100_000, -100_000, -1),
        );
        aligner.add_contig(
            "revcomp",
            false,
            &x_revcomp,
            false,
            scoring_global_custom(-100_000, -100_000, -100_000, -1),
        );
        let alignment = aligner.custom(&y);
        assert_alignment(&alignment, 0, 8, 0, 8, 8 - 1, 0, "4=1C0J4=", 8);
    }

    #[rstest]
    fn test_rev_to_fwd_jump() {
        let x = s("CCAAGGTT");
        let x_revcomp = reverse_complement(&x);
        let y = s("AACCGGTT");
        let mut aligner = MultiContigAligner::new();
        aligner.add_contig(
            "fwd",
            true,
            &x,
            false,
            scoring_global_custom(-100_000, -100_000, -100_000, -1),
        );
        aligner.add_contig(
            "revcomp",
            false,
            &x_revcomp,
            false,
            scoring_global_custom(-100_000, -100_000, -100_000, -1),
        );
        let alignment = aligner.custom(&y);
        assert_alignment(&alignment, 0, 8, 0, 8, 8 - 1, 1, "4=1c0J4=", 8);
    }

    #[rstest]
    fn test_fwd_to_rev_long_jump() {
        // x fwd: AACCAAAATTGG
        //        ||||
        // y    : AACCNNNNGGTT
        //                ||||
        // x rev: CCAA____GGTT
        let x = s("AACCAAAATTGG");
        let x_revcomp = reverse_complement(&x); // CCAATTTTGGTT
        let y = s("AACCGGTT");
        let mut aligner = MultiContigAligner::new();
        aligner.add_contig(
            "fwd",
            true,
            &x,
            false,
            scoring_global_custom(-100_000, -100_000, -100_000, -1),
        );
        aligner.add_contig(
            "revcomp",
            false,
            &x_revcomp,
            false,
            scoring_global_custom(-100_000, -100_000, -100_000, -1),
        );
        let alignment = aligner.custom(&y);
        assert_alignment(&alignment, 0, 12, 0, 8, 8 - 1, 0, "4=1C4J4=", 8);
    }

    #[rstest]
    fn test_rev_to_fwd_long_jump() {
        let x = s("CCAANNNNGGTT");
        let x_revcomp = reverse_complement(&x);
        let y = s("AACCGGTT");
        let mut aligner = MultiContigAligner::new();
        aligner.add_contig(
            "fwd",
            true,
            &x,
            false,
            scoring_global_custom(-100_000, -100_000, -100_000, -1),
        );
        aligner.add_contig(
            "revcomp",
            false,
            &x_revcomp,
            false,
            scoring_global_custom(-100_000, -100_000, -100_000, -1),
        );
        let alignment = aligner.custom(&y);
        assert_alignment(&alignment, 0, 12, 0, 8, 8 - 1, 1, "4=1c4J4=", 8);
    }

    #[rstest]
    fn test_many_contigs() {
        let x1 = s("TATATCCCCCTATATATATATATATATA");
        let x2 = s("ATATATTATATATATATATATATGGGGG");
        let x3 = s("AAAAA");
        let x4 = s("TTTTTTTTTTTTTTTT");
        let y1 = s("AAAAACCCCCGGGGGAAAAATTTTTTTTTTTTTTTT");
        let mut aligner = MultiContigAligner::new();
        let xs = vec![x1, x2, x3, x4];
        for (i, x) in xs.iter().enumerate() {
            aligner.add_contig(
                &format!("contig-{i}").to_string(),
                true,
                x,
                false,
                scoring_global_custom(-1, -100_000, -100_000, -1),
            );
        }
        let alignment = aligner.custom(&y1);
        assert_alignment(
            &alignment,
            0,
            16,
            0,
            36,
            36 - 1 - 1 - 1 - 1,
            2,
            "5=2c0J5=1C13J5=1C28j5=1C5j16=",
            36,
        );
    }
}
