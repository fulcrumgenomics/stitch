// Original Code was copied from:
//     https://github.com/rust-bio/rust-bio/blob/master/src/alignment/pairwise/mod.rs
// Copyright 2014-2015 Johannes Köster, Vadim Nazarov, Patrick Marks
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use crate::align::aligners::constants::AlignmentOperation;
use crate::align::aligners::constants::DEFAULT_ALIGNER_CAPACITY;
use crate::align::aligners::single_strand::SingleStrandAligner;
use crate::align::alignment::Alignment;
use crate::align::traceback::traceback_double_stranded;
use bio::alignment::pairwise::MatchFunc;
use bio::alignment::pairwise::MatchParams;
use bio::utils::TextSlice;
use std::i32;

use crate::align::aligners::constants::AlignmentMode;
use crate::align::aligners::constants::MIN_SCORE;
use crate::align::scoring::Scoring;

/// A generalized Smith-Waterman aligner, allowing for the alignment to jump forward
/// (or backward) in `x` (on either same strand).
///
/// `M(i,j,k)` is the best score such that the alignment jumps from where `x[k]` and `y[j]` ends
/// in a match (or substitution), to where `x[i]` and `y[j]` ends in a match (or substitution).
/// `M(i,j,i-1)` is the equivalent to Smith-Waterman, whereas `M(i,j,k)` where `k != i - 1`
/// includes a jump from a position in `x` (can be before or after `i - 1`).
/// ```ignore
///              .... A   G  x_i
///              .... C   G  y_j
/// ```
///
/// `I(i,j)` is the best score such that `x[i]` is aligned with a gap
/// ```ignore
///              .... A   G  x_i
///              .... G  y_j  -
/// ```
/// This is interpreted as an insertion into `x` w.r.t reference `y`
///
/// `D(i,j)` is the best score such that `y[j]` is aligned with a gap
/// ```ignore
///              .... A  x_i  -
///              .... G   G  y_j
/// ```
/// This is interpreted as a deletion from `x` w.r.t reference `y`
///
/// `S(i,j)` is the best score for prefixes `x[0..i]`, `y[0..j]`
///
/// To save space, only two columns of these matrices are stored at
/// any point - the current column and the previous one. Moreover
/// `M(i,j)` is not explicitly stored
///
/// `Lx` is the optimal x suffix clipping lengths from each position of the
/// sequence y
///
/// `Ly` is the optimal y suffix clipping lengths from each position of the
/// sequence x
///
/// `Sn` is the last column of the matrix. This is needed to keep track of
/// suffix clipping scores
///
/// `traceback` - see [`bio::alignment::pairwise::TracebackCell`](struct.TracebackCell.html)
///
/// `scoring` - see [`bio::alignment::pairwise::Scoring`](struct.Scoring.html)
pub struct DoubleStrandAligner<F: MatchFunc> {
    pub forward: SingleStrandAligner<F>,
    pub reverse: SingleStrandAligner<F>,
}

impl<F: MatchFunc> DoubleStrandAligner<F> {
    /// Create new aligner instance with given gap open and gap extend penalties
    /// and the score function.
    ///
    /// # Arguments
    ///
    /// * `gap_open` - the score for opening a gap (should be negative)
    /// * `gap_extend` - the score for extending a gap (should be negative)
    /// * `jump_score` - the score for jumping back in the query (should not be positive)
    /// * `match_fn` - function that returns the score for substitutions
    ///    (see also [`bio::alignment::pairwise::Scoring`](struct.Scoring.html))
    pub fn new(gap_open: i32, gap_extend: i32, jump_score: i32, match_fn: fn() -> F) -> Self {
        DoubleStrandAligner::with_capacity(
            DEFAULT_ALIGNER_CAPACITY,
            DEFAULT_ALIGNER_CAPACITY,
            gap_open,
            gap_extend,
            jump_score,
            match_fn,
        )
    }

    /// Create new aligner instance. The size hints help to
    /// avoid unnecessary memory allocations.
    ///
    /// # Arguments
    ///
    /// * `m` - the expected size of x
    /// * `n` - the expected size of y
    /// * `gap_open` - the score for opening a gap (should be negative)
    /// * `gap_extend` - the score for extending a gap (should be negative)
    /// * `jump_score` - the score for jumping back in the query (should not be positive)
    /// * `match_fn` - function that returns the score for substitutions
    ///    (see also [`bio::alignment::pairwise::Scoring`](struct.Scoring.html))
    pub fn with_capacity(
        m: usize,
        n: usize,
        gap_open: i32,
        gap_extend: i32,
        jump_score: i32,
        match_fn: fn() -> F,
    ) -> Self {
        assert!(gap_open <= 0, "gap_open can't be positive");
        assert!(gap_extend <= 0, "gap_extend can't be positive");

        DoubleStrandAligner {
            forward: SingleStrandAligner::with_capacity(
                m,
                n,
                gap_open,
                gap_extend,
                jump_score,
                match_fn(),
            ),
            reverse: SingleStrandAligner::with_capacity(
                m,
                n,
                gap_open,
                gap_extend,
                jump_score,
                match_fn(),
            ),
        }
    }

    /// Create new aligner instance with given the scoring struct
    ///
    /// # Arguments
    ///
    /// * `scoring` - the scoring struct (see bio::alignment::pairwise::Scoring)
    #[allow(dead_code)]
    pub fn with_scoring(scoring_fwd: Scoring<F>, scoring_rev: Scoring<F>) -> Self {
        DoubleStrandAligner::with_capacity_and_scoring(
            DEFAULT_ALIGNER_CAPACITY,
            DEFAULT_ALIGNER_CAPACITY,
            scoring_fwd,
            scoring_rev,
        )
    }

    /// Create new aligner instance with scoring and size hint. The size hints help to
    /// avoid unnecessary memory allocations.
    ///
    /// # Arguments
    ///
    /// * `m` - the expected size of x
    /// * `n` - the expected size of y
    /// * `scoring` - the scoring struct
    pub fn with_capacity_and_scoring(
        m: usize,
        n: usize,
        scoring_fwd: Scoring<F>,
        scoring_rev: Scoring<F>,
    ) -> Self {
        DoubleStrandAligner {
            forward: SingleStrandAligner::with_capacity_and_scoring(m, n, scoring_fwd),
            reverse: SingleStrandAligner::with_capacity_and_scoring(m, n, scoring_rev),
        }
    }

    fn fill_x_buffer_stranded(&mut self, m: usize) {
        for i in 0..=m {
            let fwd = self.forward.x_buffer.get(i);
            let rev = self.reverse.x_buffer.get(i);
            assert!(!fwd.flip_strand, "Bug: fwd strand");
            assert!(!rev.flip_strand, "Bug: rev strand");

            match fwd.score.cmp(&rev.score) {
                std::cmp::Ordering::Less => self.forward.x_buffer.set(i, rev.score, rev.from, true),
                std::cmp::Ordering::Greater => {
                    self.reverse.x_buffer.set(i, fwd.score, fwd.from, true)
                }
                std::cmp::Ordering::Equal => (),
            };
        }
    }

    /// The core function to compute the alignment
    ///
    /// # Arguments
    ///
    /// * `x` - Textslice
    /// * `y` - Textslice
    pub fn custom(
        &mut self,
        x_forward: TextSlice<'_>,
        x_revcomp: TextSlice<'_>,
        y: TextSlice<'_>,
    ) -> Alignment {
        let (m, n) = (x_forward.len(), y.len());
        assert!(x_forward.len() == x_revcomp.len(), "X length mismatch");

        // Set the initial conditions
        // We are repeating some work, but that's okay!
        self.forward.init_matrices(m, n);
        self.reverse.init_matrices(m, n);

        for j in 1..=n {
            let curr = j % 2;
            let prev = 1 - curr;

            // Initialize the column
            self.forward.init_column(j, curr, m, n);
            self.reverse.init_column(j, curr, m, n);

            // Initiliaze the jump buffers
            self.forward
                .x_buffer
                .fill(m, j - 1, &self.forward.S, &self.forward.scoring);
            self.reverse
                .x_buffer
                .fill(m, j - 1, &self.reverse.S, &self.reverse.scoring);
            self.fill_x_buffer_stranded(m);

            // Fill the column
            self.forward.fill_column(x_forward, y, m, n, j, prev, curr);
            self.reverse.fill_column(x_revcomp, y, m, n, j, prev, curr);
        }

        self.forward.fill_last_column_and_end_clipping(m, n);
        self.reverse.fill_last_column_and_end_clipping(m, n);

        // Traceback...
        traceback_double_stranded(&self.forward, &self.reverse, m, n)
    }

    /// Calculate global alignment of x against y.
    #[allow(dead_code)]
    pub fn global(
        &mut self,
        x_forward: TextSlice<'_>,
        x_revcomp: TextSlice<'_>,
        y: TextSlice<'_>,
    ) -> Alignment {
        // Store the current clip penalties
        let forward_clip_penalties = [
            self.forward.scoring.xclip_prefix,
            self.forward.scoring.xclip_suffix,
            self.forward.scoring.yclip_prefix,
            self.forward.scoring.yclip_suffix,
        ];
        let reverse_clip_penalties = [
            self.reverse.scoring.xclip_prefix,
            self.reverse.scoring.xclip_suffix,
            self.reverse.scoring.yclip_prefix,
            self.reverse.scoring.yclip_suffix,
        ];

        // Temporarily Over-write the clip penalties
        self.forward.scoring.xclip_prefix = MIN_SCORE;
        self.forward.scoring.xclip_suffix = MIN_SCORE;
        self.forward.scoring.yclip_prefix = MIN_SCORE;
        self.forward.scoring.yclip_suffix = MIN_SCORE;
        self.reverse.scoring.xclip_prefix = MIN_SCORE;
        self.reverse.scoring.xclip_suffix = MIN_SCORE;
        self.reverse.scoring.yclip_prefix = MIN_SCORE;
        self.reverse.scoring.yclip_suffix = MIN_SCORE;

        // Compute the alignment
        let mut alignment = self.custom(x_forward, x_revcomp, y);
        alignment.mode = AlignmentMode::Global;

        // Set the clip penalties to the original values
        self.forward.scoring.xclip_prefix = forward_clip_penalties[0];
        self.forward.scoring.xclip_suffix = forward_clip_penalties[1];
        self.forward.scoring.yclip_prefix = forward_clip_penalties[2];
        self.forward.scoring.yclip_suffix = forward_clip_penalties[3];
        self.reverse.scoring.xclip_prefix = reverse_clip_penalties[0];
        self.reverse.scoring.xclip_suffix = reverse_clip_penalties[1];
        self.reverse.scoring.yclip_prefix = reverse_clip_penalties[2];
        self.reverse.scoring.yclip_suffix = reverse_clip_penalties[3];
        alignment
    }

    /// Calculate semiglobal alignment of x against y (x is global, y is local).
    #[allow(dead_code)]
    pub fn semiglobal(
        &mut self,
        x_forward: TextSlice<'_>,
        x_revcomp: TextSlice<'_>,
        y: TextSlice<'_>,
    ) -> Alignment {
        // Store the current clip penalties
        let forward_clip_penalties = [
            self.forward.scoring.xclip_prefix,
            self.forward.scoring.xclip_suffix,
            self.forward.scoring.yclip_prefix,
            self.forward.scoring.yclip_suffix,
        ];
        let reverse_clip_penalties = [
            self.reverse.scoring.xclip_prefix,
            self.reverse.scoring.xclip_suffix,
            self.reverse.scoring.yclip_prefix,
            self.reverse.scoring.yclip_suffix,
        ];

        // Temporarily Over-write the clip penalties
        self.forward.scoring.xclip_prefix = MIN_SCORE;
        self.forward.scoring.xclip_suffix = MIN_SCORE;
        self.forward.scoring.yclip_prefix = 0;
        self.forward.scoring.yclip_suffix = 0;
        self.reverse.scoring.xclip_prefix = MIN_SCORE;
        self.reverse.scoring.xclip_suffix = MIN_SCORE;
        self.reverse.scoring.yclip_prefix = 0;
        self.reverse.scoring.yclip_suffix = 0;

        // Compute the alignment
        let mut alignment = self.custom(x_forward, x_revcomp, y);
        alignment.mode = AlignmentMode::Semiglobal;

        // Filter out Xclip and Yclip from alignment.operations
        {
            use self::AlignmentOperation::{Del, Ins, Match, Subst, Xclip, Xflip, Xskip};
            alignment.operations.retain(|x| {
                *x == Match
                    || *x == Subst
                    || *x == Ins
                    || *x == Del
                    || matches!(*x, Xclip(_))
                    || matches!(*x, Xskip(_))
                    || matches!(*x, Xflip(_))
            });
        }

        // Set the clip penalties to the original values

        self.forward.scoring.xclip_prefix = forward_clip_penalties[0];
        self.forward.scoring.xclip_suffix = forward_clip_penalties[1];
        self.forward.scoring.yclip_prefix = forward_clip_penalties[2];
        self.forward.scoring.yclip_suffix = forward_clip_penalties[3];
        self.reverse.scoring.xclip_prefix = reverse_clip_penalties[0];
        self.reverse.scoring.xclip_suffix = reverse_clip_penalties[1];
        self.reverse.scoring.yclip_prefix = reverse_clip_penalties[2];
        self.reverse.scoring.yclip_suffix = reverse_clip_penalties[3];

        alignment
    }

    /// Calculate local alignment of x against y.
    #[allow(dead_code)]
    pub fn local(
        &mut self,
        x_forward: TextSlice<'_>,
        x_revcomp: TextSlice<'_>,
        y: TextSlice<'_>,
    ) -> Alignment {
        // Store the current clip penalties
        let forward_clip_penalties = [
            self.forward.scoring.xclip_prefix,
            self.forward.scoring.xclip_suffix,
            self.forward.scoring.yclip_prefix,
            self.forward.scoring.yclip_suffix,
        ];
        let reverse_clip_penalties = [
            self.reverse.scoring.xclip_prefix,
            self.reverse.scoring.xclip_suffix,
            self.reverse.scoring.yclip_prefix,
            self.reverse.scoring.yclip_suffix,
        ];

        // Temporarily Over-write the clip penalties
        self.forward.scoring.xclip_prefix = 0;
        self.forward.scoring.xclip_suffix = 0;
        self.forward.scoring.yclip_prefix = 0;
        self.forward.scoring.yclip_suffix = 0;
        self.reverse.scoring.xclip_prefix = 0;
        self.reverse.scoring.xclip_suffix = 0;
        self.reverse.scoring.yclip_prefix = 0;
        self.reverse.scoring.yclip_suffix = 0;

        // Compute the alignment
        let mut alignment = self.custom(x_forward, x_revcomp, y);
        alignment.mode = AlignmentMode::Local;

        // Filter out Xclip and Yclip from alignment.operations
        {
            use self::AlignmentOperation::{Del, Ins, Match, Subst, Xclip, Xflip, Xskip};
            alignment.operations.retain(|x| {
                *x == Match
                    || *x == Subst
                    || *x == Ins
                    || *x == Del
                    || matches!(*x, Xclip(_))
                    || matches!(*x, Xskip(_))
                    || matches!(*x, Xflip(_))
            });
        }

        // Set the clip penalties to the original values

        self.forward.scoring.xclip_prefix = forward_clip_penalties[0];
        self.forward.scoring.xclip_suffix = forward_clip_penalties[1];
        self.forward.scoring.yclip_prefix = forward_clip_penalties[2];
        self.forward.scoring.yclip_suffix = forward_clip_penalties[3];
        self.reverse.scoring.xclip_prefix = reverse_clip_penalties[0];
        self.reverse.scoring.xclip_suffix = reverse_clip_penalties[1];
        self.reverse.scoring.yclip_prefix = reverse_clip_penalties[2];
        self.reverse.scoring.yclip_suffix = reverse_clip_penalties[3];

        alignment
    }
}

impl Default for DoubleStrandAligner<MatchParams> {
    fn default() -> Self {
        let match_fn = || MatchParams::new(1, -1);
        DoubleStrandAligner::new(-5, -1, -10, match_fn)
    }
}

// Tests
#[cfg(test)]
pub mod tests {
    use bio::alignment::pairwise::MatchParams;
    use itertools::Itertools;
    use rstest::rstest;

    use crate::util::dna::reverse_complement;

    use super::{Alignment, DoubleStrandAligner};

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
        is_forward: bool,
        cigar: &str,
        length: usize,
    ) {
        assert_eq!(alignment.xstart, xstart, "xstart {alignment}");
        assert_eq!(alignment.xend, xend, "xend {alignment}");
        assert_eq!(alignment.ystart, ystart, "ystart {alignment}");
        assert_eq!(alignment.yend, yend, "yend {alignment}");
        assert_eq!(alignment.score, score, "score {alignment}");
        assert_eq!(alignment.is_forward, is_forward, "strand {alignment}");
        assert_eq!(alignment.cigar(), cigar, "cigar {alignment}");
        assert_eq!(alignment.length, length, "length {alignment}");
    }

    /// Identical sequences, all matches
    #[rstest]
    fn test_identical() {
        let x = s("ACGTAACC");
        let x_revcomp = reverse_complement(&x);
        let y = s("ACGTAACC");
        let mut aligner = DoubleStrandAligner::default();
        let alignment = aligner.global(&x, &x_revcomp, &y);
        assert_alignment(&alignment, 0, 8, 0, 8, 8, true, "8=", 8);
    }

    /// Identical sequences, all matches, reverse complemented
    #[rstest]
    fn test_identical_revcomp() {
        let x = s("ACGTAACC");
        let x_revcomp = reverse_complement(&x);
        let y = reverse_complement(s("ACGTAACC"));
        let mut aligner = DoubleStrandAligner::default();
        let alignment = aligner.global(&x, &x_revcomp, &y);
        assert_alignment(&alignment, 0, 8, 0, 8, 8, false, "8=", 8);
    }

    #[rstest]
    fn test_fwd_to_fwd_jump() {
        let x = s("AAGGCCTT");
        let x_revcomp = reverse_complement(&x);
        let y = s("AACCGGTT");
        let match_fn = || MatchParams::new(1, -100_000);
        let mut aligner = DoubleStrandAligner::new(-100_000, -100_000, -1, match_fn);
        let alignment = aligner.global(&x, &x_revcomp, &y);
        assert_alignment(
            &alignment,
            0,
            8,
            0,
            8,
            8 - 1 - 1 - 1,
            true,
            "2=2J2=4j2=2J2=",
            8,
        );
    }

    #[rstest]
    fn test_fwd_to_rev_jump() {
        let x = s("AACCTTGG");
        let x_revcomp = reverse_complement(&x);
        let y = s("AACCGGTT");
        let match_fn = || MatchParams::new(1, -100_000);
        let mut aligner: DoubleStrandAligner<MatchParams> =
            DoubleStrandAligner::new(-100_000, -100_000, -1, match_fn);
        let alignment = aligner.global(&x, &x_revcomp, &y);
        assert_alignment(&alignment, 0, 8, 0, 8, 8 - 1, true, "4=0f4=", 8);
    }

    #[rstest]
    fn test_rev_to_fwd_jump() {
        let x = s("CCAAGGTT");
        let x_revcomp = reverse_complement(&x);
        let y = s("AACCGGTT");
        let match_fn = || MatchParams::new(1, -100_000);
        let mut aligner: DoubleStrandAligner<MatchParams> =
            DoubleStrandAligner::new(-100_000, -100_000, -1, match_fn);
        let alignment = aligner.global(&x, &x_revcomp, &y);
        assert_alignment(&alignment, 0, 8, 0, 8, 8 - 1, false, "4=0f4=", 8);
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
        let match_fn = || MatchParams::new(1, -100_000);
        let mut aligner: DoubleStrandAligner<MatchParams> =
            DoubleStrandAligner::new(-100_000, -100_000, -1, match_fn);
        let alignment = aligner.global(&x, &x_revcomp, &y);
        assert_alignment(&alignment, 0, 12, 0, 8, 8 - 1, true, "4=4F4=", 8);
    }

    #[rstest]
    fn test_rev_to_fwd_long_jump() {
        let x = s("CCAANNNNGGTT");
        let x_revcomp = reverse_complement(&x);
        let y = s("AACCGGTT");
        let match_fn = || MatchParams::new(1, -100_000);
        let mut aligner: DoubleStrandAligner<MatchParams> =
            DoubleStrandAligner::new(-100_000, -100_000, -1, match_fn);
        let alignment = aligner.global(&x, &x_revcomp, &y);
        assert_alignment(&alignment, 0, 12, 0, 8, 8 - 1, false, "4=4F4=", 8);
    }
}
