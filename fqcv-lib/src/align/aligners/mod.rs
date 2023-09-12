// Original Code was copied from:
//     https://github.com/rust-bio/rust-bio/blob/master/src/alignment/pairwise/mod.rs
// Copyright 2014-2015 Johannes Köster, Vadim Nazarov, Patrick Marks
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.
pub(crate) mod constants;
pub(crate) mod multi_contig_aligner;
pub(crate) mod single_contig_aligner;

pub use constants::AlignmentMode;

use derive_builder::Builder;

use crate::{
    align::{
        aligners::{
            constants::{
                AlignmentOperation::{Del, Ins, Match, Subst, Xjump},
                MIN_SCORE,
            },
            multi_contig_aligner::MultiContigAligner,
        },
        alignment::Alignment,
        scoring::Scoring,
        sub_alignment::SubAlignmentBuilder,
        PrimaryPickingStrategy,
    },
    util::{
        dna::reverse_complement,
        index_map::IndexMap,
        target_seq::{TargetHash, TargetSeq},
    },
};
use anyhow::{ensure, Context, Result};
use bio::alignment::{
    pairwise::{banded::Aligner as BandedAligner, MatchFunc, MatchParams, Scoring as BioScoring},
    sparse::HashMapFx as BandedHashMapFx,
};
use bit_set::BitSet;
use constants::DEFAULT_ALIGNER_CAPACITY;
use itertools::Itertools;
use noodles::{
    core::Position,
    sam::{
        alignment::Record as SamRecord,
        record::{
            cigar::op::{Kind, Op},
            data::field::tag::{ALIGNMENT_HIT_COUNT, ALIGNMENT_SCORE, HIT_INDEX, TOTAL_HIT_COUNT},
            Cigar, Data, Flags, MappingQuality, QualityScores, ReadName as SamReadName, Sequence,
        },
    },
};
use seq_io::fastq::{OwnedRecord as FastqOwnedRecord, Record as FastqRecord};

#[derive(Default, Debug, PartialEq, Eq, Copy, Clone)]
pub struct JumpInfo {
    score: i32,
    len: u32,
    idx: u32,
    from: u32,
}

// TODO: impl Default with same values from CLI
#[derive(Copy, Clone, Debug, Builder)]
#[builder(name = "Builder", build_fn(name = "build_options"))]
pub struct Options {
    #[builder(default)]
    mode: AlignmentMode,
    #[builder(default = "1")]
    match_score: i32,
    #[builder(default = "-4")]
    mismatch_score: i32,
    #[builder(default = "-6")]
    gap_open: i32,
    #[builder(default = "-2")]
    gap_extend: i32,
    #[builder(default = "-10")]
    default_jump_score: i32,
    jump_score_same_contig_and_strand: Option<i32>,
    jump_score_same_contig_opposite_strand: Option<i32>,
    jump_score_inter_contig: Option<i32>,
    #[builder(default = "12")]
    kmer_size: usize,
    #[builder(default = "50")]
    band_width: usize,
    #[builder(default = "false")]
    double_strand: bool,
    #[builder(default = "false")]
    circular: bool,
    #[builder(default = "20")]
    circular_slop: usize,
    #[builder(default = "false")]
    pre_align: bool,
    #[builder(default = "100")]
    pre_align_min_score: i32,
    #[builder(default = "true")]
    pre_align_subset_contigs: bool,
    #[builder(default = "false")]
    suboptimal: bool,
    #[builder(default = "20.0")]
    suboptimal_pct: f32,
    #[builder(default = "false")]
    soft_clip: bool,
    #[builder(default = "false")]
    use_eq_and_x: bool,
    #[builder(default = "PrimaryPickingStrategy::default()")]
    pick_primary: PrimaryPickingStrategy,
    #[builder(default = "false")]
    filter_secondary: bool,
    #[builder(default = "10.0")]
    filter_secondary_pct: f32,
}

impl Options {
    fn match_params(&self) -> MatchParams {
        MatchParams::new(self.match_score, self.mismatch_score)
    }

    fn clipping(&self) -> (i32, i32, i32, i32) {
        match self.mode {
            AlignmentMode::Local => (0, 0, 0, 0),
            AlignmentMode::QueryLocal => (MIN_SCORE, MIN_SCORE, 0, 0),
            AlignmentMode::TargetLocal => (0, 0, MIN_SCORE, MIN_SCORE),
            AlignmentMode::Global => (MIN_SCORE, MIN_SCORE, MIN_SCORE, MIN_SCORE),
            AlignmentMode::Custom => panic!("Custom alignment mode not supported"), // TODO: move to main run method
        }
    }

    fn banded_scoring(&self) -> BioScoring<MatchParams> {
        let match_params = self.match_params();
        let (xclip_prefix, xclip_suffix, yclip_prefix, yclip_suffix) = self.clipping();
        BioScoring::new(self.gap_open, self.gap_extend, match_params)
            .xclip_prefix(xclip_prefix)
            .xclip_suffix(xclip_suffix)
            .yclip_prefix(yclip_prefix)
            .yclip_suffix(yclip_suffix)
    }

    fn contig_scoring(&self) -> Scoring<MatchParams> {
        let jump_score_same_contig_and_strand = self
            .jump_score_same_contig_and_strand
            .unwrap_or(self.default_jump_score);
        let jump_score_same_contig_opposite_strand = self
            .jump_score_same_contig_opposite_strand
            .unwrap_or(self.default_jump_score);
        let jump_score_inter_contig = self
            .jump_score_inter_contig
            .unwrap_or(self.default_jump_score);
        let match_params = self.match_params();
        let (xclip_prefix, xclip_suffix, yclip_prefix, yclip_suffix) = self.clipping();
        Scoring::with_jump_scores(
            self.gap_open,
            self.gap_extend,
            jump_score_same_contig_and_strand,
            jump_score_same_contig_opposite_strand,
            jump_score_inter_contig,
            match_params.clone(),
        )
        .set_xclip_prefix(xclip_prefix)
        .set_xclip_suffix(xclip_suffix)
        .set_yclip_prefix(yclip_prefix)
        .set_yclip_suffix(yclip_suffix)
    }
}

impl Builder {
    pub fn build_aligners<'a>(&self, target_seqs: &'a [TargetSeq]) -> Aligners<'a, MatchParams> {
        let opts = self.build_options().unwrap();
        // Banded alignment is always local since the goal is to find at least some minimal scoring
        // local alignment.
        let banded = BandedAligner::with_capacity_and_scoring(
            DEFAULT_ALIGNER_CAPACITY,
            DEFAULT_ALIGNER_CAPACITY,
            opts.banded_scoring(),
            opts.kmer_size,
            opts.band_width,
        );
        let capacity = target_seqs.len() * (if opts.double_strand { 2 } else { 1 });
        let mut multi_contig: MultiContigAligner<'a, MatchParams> =
            MultiContigAligner::with_capacity(capacity);
        let multi_contig_scoring = opts.contig_scoring();
        for target_seq in target_seqs {
            multi_contig.add_contig(
                &target_seq.name,
                true,
                &target_seq.fwd,
                opts.circular,
                multi_contig_scoring.clone(),
            );
        }
        if opts.double_strand {
            for target_seq in target_seqs {
                multi_contig.add_contig(
                    &target_seq.name,
                    false,
                    &target_seq.revcomp,
                    opts.circular,
                    multi_contig_scoring.clone(),
                );
            }
        }
        Aligners {
            banded,
            multi_contig,
            opts,
        }
    }

    pub fn build_sam_record_formatter<'a>(
        &self,
        target_seqs: &'a [TargetSeq],
    ) -> SamRecordFormatter<'a, MatchParams> {
        let opts = self.build_options().unwrap();
        let scoring = opts.contig_scoring();
        SamRecordFormatter {
            target_seqs,
            scoring,
            opts,
        }
    }
}

pub struct Aligners<'a, F: MatchFunc> {
    // Aligner used to quickly determine if there are ANY high-quality local alignments.
    banded: BandedAligner<MatchParams>,
    // Aligner used when there are more than one contig (or double strand, or both)
    multi_contig: MultiContigAligner<'a, F>,
    // The alignment mode
    opts: Options,
}

impl Aligners<'_, MatchParams> {
    pub fn align(
        &mut self,
        record: &FastqOwnedRecord,
        target_seqs: &[TargetSeq],
        target_hashes: &[TargetHash],
    ) -> (Vec<Alignment>, Option<i32>) {
        let query = record
            .seq()
            .iter()
            .map(u8::to_ascii_uppercase)
            .collect_vec();
        let mut contig_idx_to_prealign_score: IndexMap<i32> =
            IndexMap::new(self.multi_contig.len());
        if self.opts.pre_align {
            // Align the record to all the targets in a local banded alignment. If there is at least one
            // alignment with minimum score, then align uses the full aligner.
            for (index, target_seq) in target_seqs.iter().enumerate() {
                let target_hash = &target_hashes[index];
                let (score_fwd, score_revcomp) = prealign_local_banded(
                    &query,
                    target_seq,
                    target_hash,
                    &mut self.banded,
                    self.opts.double_strand,
                    self.opts.pre_align_min_score,
                );
                if let Some(score) = score_fwd {
                    contig_idx_to_prealign_score.put(
                        self.multi_contig
                            .contig_index_for_strand(true, &target_seq.name)
                            .unwrap(),
                        score,
                    );
                }
                if let Some(score) = score_revcomp {
                    contig_idx_to_prealign_score.put(
                        self.multi_contig
                            .contig_index_for_strand(false, &target_seq.name)
                            .unwrap(),
                        score,
                    );
                }
                // If we are going to align to all the contigs anyhow, then we can stop here.
                if !self.opts.pre_align_subset_contigs && !contig_idx_to_prealign_score.is_empty() {
                    break;
                }
            }
            // If there was no contig with a good enough alignment, return None now.
            if contig_idx_to_prealign_score.is_empty() {
                return (Vec::new(), None);
            }
        }

        // Get the contigs to align based on if we pre-aligned or not
        let contigs_to_align: Option<BitSet<u32>> =
            if self.opts.pre_align && self.opts.pre_align_subset_contigs {
                let indexes = contig_idx_to_prealign_score.keys().collect::<BitSet<u32>>();
                assert!(!indexes.is_empty(), "Bug: should have returned above");
                Some(indexes)
            } else {
                // Use all the contigs!
                None
            };

        // Align to all the contigs! (or those that had a "good enough" pre-align score)
        // This populates the traceback matrices too for suboptimal alignments.
        let original_alignment = self.multi_contig_align(&query, contigs_to_align.as_ref());

        // Get all alignments if we want to keep sub-optimal alignments, or just process this one
        let mut alignments = Vec::new();
        if self.opts.suboptimal {
            // Use the traceback from the original alignment to get sub-optimal alignments
            let new_alignments = self
                .multi_contig
                .traceback_all(query.len(), contigs_to_align.as_ref());
            // Re-align around the origin if the contig is circular or we force circular
            // NB: must do this after tracing back all the above since we modify the traceback matrix below.
            for alignment in new_alignments {
                // remove leading/trailing clipping, needed to for origin re-alignment
                let alignment = self.remove_clipping(alignment);
                let alignment =
                    self.realign_origin(&query, alignment, self.opts.circular_slop, false);
                alignments.push(alignment);
            }

            // Filter out sub-optimal alignments
            if alignments.len() > 1 {
                alignments.sort_by_key(|a| -a.score); // descending
                let min_score = alignments[0].score as f32 * self.opts.suboptimal_pct / 100.0;
                let mut new_alignments = Vec::new();
                for alignment in alignments {
                    if alignment.score as f32 >= min_score {
                        new_alignments.push(alignment);
                    }
                }
                alignments = new_alignments;
            }
        } else {
            // Re-align around the origin if the contig is circular or we force circular
            let alignment =
                self.realign_origin(&query, original_alignment, self.opts.circular_slop, false);
            alignments.push(alignment);
        }

        // Get the maximum pre-align score to return
        let prealign_score: Option<i32> = contig_idx_to_prealign_score.values().copied().max();
        (alignments, prealign_score)
    }

    /// Removes leading and trailing clipping
    fn remove_clipping(&self, mut aln: Alignment) -> Alignment {
        match self.opts.mode {
            AlignmentMode::Local | AlignmentMode::QueryLocal | AlignmentMode::TargetLocal => {
                aln.operations
                    .retain(|x| matches!(*x, Match | Subst | Ins | Del | Xjump(_, _)));
            }
            AlignmentMode::Global => (), // do nothing, there can be no clipping!
            AlignmentMode::Custom => unreachable!(),
        }
        aln
    }

    fn multi_contig_align(
        &mut self,
        query: &[u8],
        contig_indexes: Option<&BitSet<u32>>,
    ) -> Alignment {
        // Check if we should subset the contig indexes
        let aln = self.multi_contig.custom_with_subset(query, contig_indexes);
        self.remove_clipping(aln)
    }

    fn get_start_and_end_contig_indexes_for_realignment(
        &self,
        alignment: &Alignment,
        slop: usize,
    ) -> (Option<usize>, Option<usize>) {
        let contig_at_start: Option<usize> = if alignment.xstart <= slop
            && self.multi_contig.is_circular(alignment.start_contig_idx)
        {
            Some(alignment.start_contig_idx)
        } else {
            None
        };

        // Get the contig at the end of the read
        let contig_at_end: Option<usize> = if alignment.xlen <= alignment.xend + slop
            && self.multi_contig.is_circular(alignment.end_contig_idx)
        {
            Some(alignment.end_contig_idx)
        } else {
            None
        };

        // If the alignment starts and ends on the same contig, do not re-align.  If it doesn't
        // start or end close to the start or end of the contig respectively, also don't re-align.
        match (contig_at_start, contig_at_end) {
            (Some(start), Some(end)) if start == end => {
                return (None, None);
            }
            (None, None) => return (None, None),
            _ => (),
        }

        // Ensure that there are additional `y` bases to align!
        let contig_at_start = if contig_at_start.is_none() || alignment.yend == alignment.ylen {
            None
        } else {
            contig_at_start
        };
        let contig_at_end = if contig_at_end.is_none() || 0 == alignment.ystart {
            None
        } else {
            contig_at_end
        };

        (contig_at_start, contig_at_end)
    }

    fn realign_and_split_at_y(
        &mut self,
        query: &[u8],
        best_alignment: &Alignment,
        contig_indexes: &Option<BitSet<u32>>,
        contig_index: usize,
        y_pivot: usize,
    ) -> Option<Alignment> {
        self.multi_contig_align(query, contig_indexes.as_ref()); // to get the traceback matrix
        let new_alignment = self.multi_contig.traceback_from(query.len(), contig_index);
        if let Some(new_alignment) = new_alignment {
            if new_alignment.score > best_alignment.score
                && new_alignment.start_contig_idx == contig_index
                && best_alignment.end_contig_idx == contig_index
            {
                return Some(self.remove_clipping(new_alignment).split_at_y(y_pivot));
            }
        }
        None
    }

    /// Realign alignments where `y` may align across the origin.
    ///
    /// If the alignment is within `slop` from the begging of the alignment start contig, split
    /// the `y` into two, and append the prefix to the suffix, then realign the new `y`.  If the
    /// prefix aligns to the end of the contig, then we should have have an alignment of the "new"
    /// `y` where there's a Xjump.
    ///
    /// The same is true for if the original alignment is within `slop` from the end of the
    /// alignment start contig...
    fn realign_origin(
        &mut self,
        query: &[u8],
        alignment: Alignment,
        slop: usize,
        all_contigs: bool,
    ) -> Alignment {
        let (contig_at_start, contig_at_end) =
            self.get_start_and_end_contig_indexes_for_realignment(&alignment, slop);
        if contig_at_start.is_none() && contig_at_end.is_none() {
            return alignment;
        }

        // Build the contigs to which to align
        let contig_indexes: Option<BitSet<u32>> = if all_contigs {
            Some((0..self.multi_contig.len()).collect::<BitSet<_>>())
        } else {
            // Use the contigs in the existing alignment
            let mut indexes = BitSet::new();
            indexes.insert(alignment.start_contig_idx);
            indexes.insert(alignment.end_contig_idx);
            for op in &alignment.operations {
                if let Xjump(idx, _) = op {
                    indexes.insert(*idx);
                }
            }
            Some(indexes)
        };

        // Set the best alignment to the current alignment
        let mut best_alignment = alignment.clone();

        // The case where the current alignment aligns from the start of the contig
        if let Some(start_contig_idx) = contig_at_start {
            // The first new query to realign is all bases up to the end of the current alignment
            // and append it to the unaligned suffix.
            let first_query: Vec<u8> =
                [&query[alignment.yend..], &query[..alignment.yend]].concat();
            let first_query_and_yend = (first_query, alignment.yend);

            // The second new query ignores any bases in the alignment from other contigs.
            // Therefore the pivot point is last base not aligned to the given contig starting
            // from the beginning of the alignment.
            let mut yend = alignment.ystart;
            for op in &alignment.operations {
                if let Xjump(idx, _) = op {
                    if *idx != start_contig_idx {
                        break;
                    }
                }
                yend += op.length_on_y();
            }
            let second_query: Vec<u8> = [&query[yend..], &query[..yend]].concat();
            let second_query_and_yend = (second_query, yend);

            // Align!
            for (query, yend) in [first_query_and_yend, second_query_and_yend] {
                best_alignment = self
                    .realign_and_split_at_y(
                        &query,
                        &best_alignment,
                        &contig_indexes,
                        start_contig_idx,
                        alignment.ylen - yend,
                    )
                    .unwrap_or(best_alignment);
            }
        }

        // The case where the current alignment aligns at the end of the contig
        if let Some(end_contig_idx) = contig_at_end {
            // The first new query to realign is all bases up to the start of the current alignment
            // and append it to the unaligned prefix.
            let first_query: Vec<u8> =
                [&query[alignment.ystart..], &query[..alignment.ystart]].concat();
            let first_query_and_ystart = (first_query, alignment.ystart);

            // The second new query ignores any bases in the alignment from other contigs.
            // Therefore the pivot point is first base not aligned to the given contig starting
            // from the end of the alignment.
            let mut ystart = alignment.ystart;
            let mut ycur = alignment.ystart;
            let mut xidx = alignment.start_contig_idx;
            for op in &alignment.operations {
                if let Xjump(idx, _) = op {
                    // update ystart when we jump from another contig to the end contig
                    if *idx == end_contig_idx && xidx != end_contig_idx {
                        ystart = ycur;
                    }
                    xidx = *idx;
                }
                ycur += op.length_on_y();
            }
            let second_query: Vec<u8> = [&query[ystart..], &query[..ystart]].concat();
            let second_query_and_ystart = (second_query, ystart);

            // Align!
            for (query, ystart) in [first_query_and_ystart, second_query_and_ystart] {
                best_alignment = self
                    .realign_and_split_at_y(
                        &query,
                        &best_alignment,
                        &contig_indexes,
                        end_contig_idx,
                        alignment.ylen - ystart,
                    )
                    .unwrap_or(best_alignment);
            }
        }

        best_alignment
    }
}

fn align_local_banded<F: MatchFunc>(
    query: &[u8],
    target: &[u8],
    aligner: &mut BandedAligner<F>,
    target_kmer_hash: &BandedHashMapFx<&[u8], Vec<u32>>,
) -> i32 {
    // Compute the alignment
    aligner
        .custom_with_prehash(query, target, target_kmer_hash)
        .score
}

fn prealign_local_banded<F: MatchFunc>(
    query: &[u8],
    target_seq: &TargetSeq,
    target_hash: &TargetHash,
    banded_aligner: &mut BandedAligner<F>,
    double_strand: bool,
    pre_align_min_score: i32,
) -> (Option<i32>, Option<i32>) {
    // Try to the forward strand
    let banded_fwd = align_local_banded(
        query,
        &target_seq.fwd,
        banded_aligner,
        &target_hash.fwd_hash,
    );
    let fwd = if banded_fwd < pre_align_min_score {
        None
    } else {
        Some(banded_fwd)
    };
    let revcomp = if double_strand {
        let banded_revcomp = align_local_banded(
            query,
            &target_seq.revcomp,
            banded_aligner,
            &target_hash.revcomp_hash,
        );
        if banded_revcomp < pre_align_min_score {
            None
        } else {
            Some(banded_revcomp)
        }
    } else {
        None
    };
    (fwd, revcomp)
}

pub struct SamRecordFormatter<'a, F: MatchFunc> {
    target_seqs: &'a [TargetSeq],
    scoring: Scoring<F>,
    opts: Options,
}

fn header_to_name(header: &[u8]) -> Result<String> {
    let header: std::borrow::Cow<str> = String::from_utf8_lossy(header);
    header
        .split_whitespace()
        .next()
        .map(std::string::ToString::to_string)
        .context("empty read name")
}

impl<'a, F: MatchFunc> SamRecordFormatter<'a, F> {
    pub fn format(
        &self,
        fastq: &FastqOwnedRecord,
        alignments: &[Alignment],
        alt_score: Option<i32>,
    ) -> Result<Vec<SamRecord>> {
        let name = header_to_name(fastq.head())?;
        let read_name: SamReadName = name.parse()?;
        let bases = fastq.seq();
        let quals = fastq.qual();

        if alignments.is_empty() {
            let mut record = SamRecord::default();

            // read name
            *record.read_name_mut() = Some(read_name);

            // flags
            *record.flags_mut() = Flags::UNMAPPED;

            // bases
            *record.sequence_mut() = Sequence::try_from(bases.to_vec()).unwrap();

            // qualities
            *record.quality_scores_mut() = QualityScores::try_from(quals.to_vec()).unwrap();

            // cigar
            *record.cigar_mut() = Cigar::default();

            // mapping quality
            *record.mapping_quality_mut() = MappingQuality::new(0);

            if let Some(score) = alt_score {
                let mut data = Data::default();
                data.insert(
                    "xs".parse().unwrap(),
                    noodles::sam::record::data::field::Value::from(score),
                );
                *record.data_mut() = data;
            }

            return Ok([record].to_vec());
        }

        let mut records = Vec::new();
        let mut is_first = true;
        let mut primary_score = MIN_SCORE;
        for alignment in alignments {
            let hard_clip = !self.opts.soft_clip;
            let mut builder: SubAlignmentBuilder = SubAlignmentBuilder::new(self.opts.use_eq_and_x);
            let mut subs = builder.build(alignment, true, &self.scoring);
            ensure!(!subs.is_empty());

            let mut primary_index = if is_first {
                let idx = match self.opts.pick_primary {
                    PrimaryPickingStrategy::QueryLength => subs
                        .iter()
                        .enumerate()
                        .max_by_key(|(_, alignment)| {
                            (alignment.query_end - alignment.query_start, alignment.score)
                        })
                        .map_or(0, |(index, _)| index),
                    PrimaryPickingStrategy::Score => subs
                        .iter()
                        .enumerate()
                        .max_by_key(|(_, alignment)| {
                            (alignment.score, alignment.query_end - alignment.query_start)
                        })
                        .map_or(0, |(index, _)| index),
                };
                primary_score = subs[idx].score;
                idx
            } else {
                subs.len()
            };

            // Filter out sub-alignments that have score worse than X% of the primary
            // if opts.  This changes the primary index!
            if self.opts.filter_secondary {
                let min_score = primary_score as f32 * self.opts.filter_secondary_pct / 100.0;
                let mut new_subs = Vec::new();
                let mut old_index = 0;
                while !subs.is_empty() {
                    let sub = subs.remove(0);
                    if old_index == primary_index {
                        primary_index = new_subs.len();
                    }
                    if sub.score as f32 >= min_score {
                        new_subs.push(sub);
                    }
                    old_index += 1;
                }
                subs = new_subs;
            }

            for (index, sub) in subs.iter().enumerate() {
                let is_secondary = index != primary_index;
                let mut record = SamRecord::default();
                assert!(sub.contig_idx < 2 * self.target_seqs.len());
                let is_forward: bool = sub.contig_idx < self.target_seqs.len();

                // read name
                *record.read_name_mut() = Some(read_name.clone());

                // flags
                let mut new_flags = Flags::default();
                if !is_forward {
                    new_flags.insert(Flags::REVERSE_COMPLEMENTED);
                }
                if is_secondary {
                    new_flags.insert(Flags::SECONDARY);
                }
                *record.flags_mut() = new_flags;

                let (bases_vec, quals_vec, cigar) = match (is_forward, hard_clip && is_secondary) {
                    (true, false) => (bases.to_vec(), quals.to_vec(), sub.cigar.clone()),
                    (true, true) => (
                        bases[sub.query_start..sub.query_end].to_vec(),
                        quals[sub.query_start..sub.query_end].to_vec(),
                        Cigar::try_from(sub.cigar.iter().rev().copied().collect::<Vec<Op>>())
                            .unwrap(),
                    ),
                    (false, false) => (
                        reverse_complement(bases),
                        quals.iter().copied().rev().collect(),
                        Cigar::try_from(sub.cigar.iter().rev().copied().collect::<Vec<Op>>())
                            .unwrap(),
                    ),
                    (false, true) => (
                        reverse_complement(bases[sub.query_start..sub.query_end].to_vec()),
                        quals[sub.query_start..sub.query_end]
                            .iter()
                            .copied()
                            .rev()
                            .collect(),
                        Cigar::try_from(sub.cigar.iter().rev().copied().collect::<Vec<Op>>())
                            .unwrap(),
                    ),
                };

                // bases
                *record.sequence_mut() = Sequence::try_from(bases_vec).unwrap();

                // qualities
                *record.quality_scores_mut() = QualityScores::try_from(quals_vec).unwrap();

                // cigar
                let clip_op = if hard_clip && is_secondary {
                    Kind::HardClip
                } else {
                    Kind::SoftClip
                };
                let mut cigar_ops = Vec::new();
                // Clip the start of the alignment
                let clip_prefix_len = if is_forward {
                    sub.query_start
                } else {
                    bases.len() - sub.query_end
                };
                if clip_prefix_len > 0 {
                    cigar_ops.push(Op::new(clip_op, clip_prefix_len));
                }
                // Add the CIGAR from the alignment
                cigar_ops.extend(cigar.iter());
                // Clip the end of the alignment
                let clip_suffix_len = if is_forward {
                    bases.len() - sub.query_end
                } else {
                    sub.query_start
                };
                if clip_suffix_len > 0 {
                    cigar_ops.push(Op::new(clip_op, clip_suffix_len));
                }
                *record.cigar_mut() = Cigar::try_from(cigar_ops).unwrap();

                // target id
                *record.reference_sequence_id_mut() = Some(sub.contig_idx % self.target_seqs.len());

                // target start
                if is_forward {
                    *record.alignment_start_mut() = Position::new(sub.target_start + 1);
                } else {
                    let target_len =
                        self.target_seqs[sub.contig_idx % self.target_seqs.len()].len();
                    *record.alignment_start_mut() = Position::new(target_len - sub.target_end + 1);
                }

                // mapping quality
                // TODO: base this on the alt_score
                if is_first {
                    *record.mapping_quality_mut() = MappingQuality::new(60);
                } else {
                    *record.mapping_quality_mut() = MappingQuality::new(0);
                }

                // TODO: tags (e.g. XS, NM, MD)
                let mut data = Data::default();
                data.insert(
                    "qs".parse().unwrap(),
                    noodles::sam::record::data::field::Value::from(sub.query_start as u32),
                );
                data.insert(
                    "qe".parse().unwrap(),
                    noodles::sam::record::data::field::Value::from(sub.query_end as u32),
                );
                data.insert(
                    "ts".parse().unwrap(),
                    noodles::sam::record::data::field::Value::from(sub.target_start as u32),
                );
                data.insert(
                    "te".parse().unwrap(),
                    noodles::sam::record::data::field::Value::from(sub.target_end as u32),
                );
                data.insert(
                    "as".parse().unwrap(),
                    noodles::sam::record::data::field::Value::from(alignment.score),
                );
                if let Some(score) = alt_score {
                    data.insert(
                        "xs".parse().unwrap(),
                        noodles::sam::record::data::field::Value::from(score),
                    );
                }
                data.insert(
                    ALIGNMENT_SCORE,
                    noodles::sam::record::data::field::Value::from(sub.score),
                );
                data.insert(
                    HIT_INDEX,
                    noodles::sam::record::data::field::Value::from(index as i32 + 1),
                );
                data.insert(
                    ALIGNMENT_HIT_COUNT,
                    noodles::sam::record::data::field::Value::from(subs.len() as i32),
                );
                data.insert(
                    TOTAL_HIT_COUNT,
                    noodles::sam::record::data::field::Value::from(subs.len() as i32),
                );
                *record.data_mut() = data;

                records.push(record);
            }

            is_first = false;
        }

        Ok(records)
    }
}

#[cfg(test)]
pub mod tests {
    use super::Builder;
    use crate::util::target_seq::{self, TargetHash};
    use seq_io::fastq::OwnedRecord as FastqOwnedRecord;

    #[test]
    fn test_case_insensitive() {
        let seq = b"ACGGACAGATCGAATACGACAGGAC".to_vec();
        let target_seqs = [target_seq::TargetSeq::new("test-contig", &seq, false)];
        let mut aligners = Builder::default().build_aligners(&target_seqs);
        let record = FastqOwnedRecord {
            head: b"test-record".to_vec(),
            seq: seq.clone(),
            qual: vec![b'#'; seq.len()],
        };
        let k = 7;
        let target_hashes: Vec<TargetHash> = target_seqs
            .iter()
            .map(|target_seq| target_seq.build_target_hash(k))
            .collect();
        let (alignment, _) = aligners.align(&record, &target_seqs, &target_hashes);
        assert_eq!(alignment.len(), 1);
        assert_eq!(alignment[0].length, seq.len());
        assert_eq!(alignment[0].cigar(), format!("{}=", seq.len()));
    }
}
