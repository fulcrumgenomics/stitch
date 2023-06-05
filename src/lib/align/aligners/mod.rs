// Original Code was copied from:
//     https://github.com/rust-bio/rust-bio/blob/master/src/alignment/pairwise/mod.rs
// Copyright 2014-2015 Johannes Köster, Vadim Nazarov, Patrick Marks
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use anyhow::{ensure, Context, Result};
use bio::alignment::pairwise::MatchFunc;

use crate::align::aligners::constants::AlignmentOperation::{Del, Ins, Match, Subst, Xjump};
use crate::align::PrimaryPickingStrategy;
use crate::commands::align::Align;
use crate::util::dna::reverse_complement;
use crate::util::target_seq::{TargetHash, TargetSeq};
use itertools::{self, Itertools};

use crate::align::aligners::constants::{AlignmentMode, MIN_SCORE};
use crate::align::aligners::multi_contig_aligner::MultiContigAligner;
use crate::align::alignment::Alignment;
use crate::align::scoring::Scoring;
use bio::alignment::pairwise::banded::Aligner as BandedAligner;
use bio::alignment::pairwise::MatchParams;
use bio::alignment::pairwise::Scoring as BioScoring;
use bio::alignment::sparse::HashMapFx as BandedHashMapFx;
use noodles::core::Position;
use noodles::sam::alignment::Record as SamRecord;
use noodles::sam::record::cigar::op::Kind;
use noodles::sam::record::cigar::op::Op;
use noodles::sam::record::data::field::tag::ALIGNMENT_HIT_COUNT;
use noodles::sam::record::data::field::tag::ALIGNMENT_SCORE;
use noodles::sam::record::data::field::tag::HIT_INDEX;
use noodles::sam::record::data::field::tag::TOTAL_HIT_COUNT;
use noodles::sam::record::Cigar;
use noodles::sam::record::Data;
use noodles::sam::record::Flags;
use noodles::sam::record::MappingQuality;
use noodles::sam::record::QualityScores;
use noodles::sam::record::ReadName as SamReadName;
use noodles::sam::record::Sequence;
use seq_io::fastq::OwnedRecord as FastqOwnedRecord;
use seq_io::fastq::Record as FastqRecord;

use self::constants::DEFAULT_ALIGNER_CAPACITY;

use super::sub_alignment::SubAlignmentBuilder;

pub(crate) mod constants;
pub(crate) mod multi_contig_aligner;
pub(crate) mod single_contig_aligner;

#[derive(Default, Debug, PartialEq, Eq, Copy, Clone)]
pub struct JumpInfo {
    score: i32,
    len: u32,
    idx: u32,
    from: u32,
}

struct ScoringBuilder {
    pub match_score: i32,
    pub mismatch_score: i32,
    pub gap_open: i32,
    pub gap_extend: i32,
    pub jump_score_same_contig_and_strand: i32,
    pub jump_score_same_contig_opposite_strand: i32,
    pub jump_score_inter_contig: i32,
    pub xclip_prefix: i32,
    pub xclip_suffix: i32,
    pub yclip_prefix: i32,
    pub yclip_suffix: i32,
}

impl ScoringBuilder {
    fn build_match_fn(match_score: i32, mismatch_score: i32) -> MatchParams {
        MatchParams::new(match_score, mismatch_score)
    }

    pub fn build_scoring(&self) -> Scoring<MatchParams> {
        let mut scoring = Scoring::with_jump_scores(
            self.gap_open,
            self.gap_extend,
            self.jump_score_same_contig_and_strand,
            self.jump_score_same_contig_opposite_strand,
            self.jump_score_inter_contig,
            Self::build_match_fn(self.match_score, self.mismatch_score),
        );
        scoring.xclip_prefix = self.xclip_prefix;
        scoring.xclip_suffix = self.xclip_suffix;
        scoring.yclip_prefix = self.yclip_prefix;
        scoring.yclip_suffix = self.yclip_suffix;
        scoring
    }

    pub fn build_bio_scoring(&self) -> BioScoring<MatchParams> {
        let mut scoring = BioScoring::new(
            self.gap_open,
            self.gap_extend,
            Self::build_match_fn(self.match_score, self.mismatch_score),
        );
        scoring.xclip_prefix = self.xclip_prefix;
        scoring.xclip_suffix = self.xclip_suffix;
        scoring.yclip_prefix = self.yclip_prefix;
        scoring.yclip_suffix = self.yclip_suffix;
        scoring
    }
}

pub struct Aligners<'a, F: MatchFunc> {
    // Aligner used to quickly determine if there are ANY high-quality local alignments.
    banded: BandedAligner<MatchParams>,
    // Aligner used when there are more than one contig (or double strand, or both)
    multi_contig: MultiContigAligner<'a, F>,
    // The alignment mode
    mode: AlignmentMode,
}

pub fn build_aligners<'a>(opts: &Align, target_seqs: &'a [TargetSeq]) -> Aligners<'a, MatchParams> {
    let (xclip_prefix, xclip_suffix, yclip_prefix, yclip_suffix) = match opts.mode {
        AlignmentMode::Local => (0, 0, 0, 0),
        AlignmentMode::QueryLocal => (MIN_SCORE, MIN_SCORE, 0, 0),
        AlignmentMode::TargetLocal => (0, 0, MIN_SCORE, MIN_SCORE),
        AlignmentMode::Global => (MIN_SCORE, MIN_SCORE, MIN_SCORE, MIN_SCORE),
        AlignmentMode::Custom => panic!("Custom alignment mode not supported"), // TODO: move to main run method
    };
    let scoring_builder = ScoringBuilder {
        match_score: opts.match_score,
        mismatch_score: opts.mismatch_score,
        gap_open: opts.gap_open,
        gap_extend: opts.gap_extend,
        jump_score_same_contig_and_strand: opts
            .jump_score_same_contig_and_strand
            .unwrap_or(opts.jump_score),
        jump_score_same_contig_opposite_strand: opts
            .jump_score_same_contig_opposite_strand
            .unwrap_or(opts.jump_score),
        jump_score_inter_contig: opts.jump_score_inter_contig.unwrap_or(opts.jump_score),
        xclip_prefix,
        xclip_suffix,
        yclip_prefix,
        yclip_suffix,
    };

    // Banded alignment is always local since the goal is to find at leaset some minimal scoring
    // local alignment.
    let banded = BandedAligner::with_capacity_and_scoring(
        DEFAULT_ALIGNER_CAPACITY,
        DEFAULT_ALIGNER_CAPACITY,
        scoring_builder.build_bio_scoring(),
        opts.k,
        opts.w,
    );

    let mut multi_contig: MultiContigAligner<'a, MatchParams> = MultiContigAligner::new();
    for target_seq in target_seqs {
        let opts = &opts.clone();
        multi_contig.add_contig(
            &target_seq.name,
            true,
            &target_seq.fwd,
            opts.circular,
            scoring_builder.build_scoring(),
        );
    }
    if opts.double_strand {
        for target_seq in target_seqs {
            multi_contig.add_contig(
                &target_seq.name,
                false,
                &target_seq.revcomp,
                opts.circular,
                scoring_builder.build_scoring(),
            );
        }
    }
    Aligners {
        banded,
        multi_contig,
        mode: opts.mode,
    }
}

impl Aligners<'_, MatchParams> {
    pub fn align(
        &mut self,
        record: &FastqOwnedRecord,
        target_seqs: &[TargetSeq],
        target_hashes: &[TargetHash],
        pre_align: bool,
        pre_align_min_score: i32,
        circular_slop: usize,
    ) -> (Option<Alignment>, Option<i32>) {
        let query = record.seq();
        let mut prealign_score: Option<i32> = None;
        if pre_align {
            // Align the record to all the targets in a local banded alignment. If there is at least one
            // alignment with minimum score, then align uses the full aligner.
            for (index, target_seq) in target_seqs.iter().enumerate() {
                let target_hash = &target_hashes[index];
                prealign_score = prealign_local_banded(
                    query,
                    target_seq,
                    target_hash,
                    &mut self.banded,
                    pre_align_min_score,
                );
                if prealign_score.is_some() {
                    break;
                }
            }
        }
        if pre_align && prealign_score.is_none() {
            (None, None)
        } else {
            let original_alignment = self.multi_contig_align(query);
            let alignment = self
                .realign_origin(query, &original_alignment, circular_slop)
                .or(Some(original_alignment));
            (alignment, prealign_score)
        }
    }

    fn multi_contig_align(&mut self, query: &[u8]) -> Alignment {
        let mut aln = self.multi_contig.custom(query);
        match self.mode {
            AlignmentMode::Local | AlignmentMode::QueryLocal | AlignmentMode::TargetLocal => {
                aln.operations.retain(|x| {
                    *x == Match
                        || *x == Subst
                        || *x == Ins
                        || *x == Del
                        || matches!(*x, Xjump(_, _))
                });
            }
            AlignmentMode::Global => (), // do nothing
            AlignmentMode::Custom => unreachable!(),
        }
        aln
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
        alignment: &Alignment,
        slop: usize,
    ) -> Option<Alignment> {
        // self.multi_contig
        //     .contigs
        //     .iter()
        //     .find(|contig| contig.aligner.contig_idx);

        //alignment.start_contig_idx

        // Get the contig at the start of the read
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
        // start or end close to the start or end of the contig respecitvely, also don't re-align.
        match (contig_at_start, contig_at_end) {
            (Some(start), Some(end)) => {
                if start == end {
                    return None;
                }
            }
            (None, None) => return None,
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
            contig_at_start
        };

        let start_alignment: Option<Alignment> = if let Some(start_contig_idx) = contig_at_start {
            // The current alignment aligns to the start of the contig, so take all bases up to
            // the end of the current alignment and append it to the unaligned suffix.
            let new_query: Vec<u8> = [&query[alignment.yend..], &query[..alignment.yend]].concat();
            let new_alignment = self.multi_contig_align(&new_query);
            // TODO: ensure that the new alignment crosses the origin
            if new_alignment.score > alignment.score
                && new_alignment.start_contig_idx == start_contig_idx
                && alignment.end_contig_idx == start_contig_idx
            {
                Some(new_alignment.split_at_y(alignment.ylen - alignment.yend))
            } else {
                None
            }
        } else {
            None
        };

        let end_alignment: Option<Alignment> = if let Some(end_contig_idx) = contig_at_end {
            // The current alignment aligns to the end of the contig, so take all bases up to
            // the start of the current alignment and append it to the unaligned prefix.
            let new_query: Vec<u8> =
                [&query[alignment.ystart..], &query[..alignment.ystart]].concat();
            let new_alignment = self.multi_contig_align(&new_query);
            // TODO: ensure that the new alignment crosses the origin
            if new_alignment.score > alignment.score
                && new_alignment.end_contig_idx == end_contig_idx
                && alignment.start_contig_idx == end_contig_idx
            {
                Some(new_alignment.split_at_y(alignment.ylen - alignment.ystart))
            } else {
                None
            }
        } else {
            None
        };

        // Pick between the alignments
        match (start_alignment, end_alignment) {
            (Some(start), Some(end)) => {
                if start.score >= end.score {
                    Some(start)
                } else {
                    Some(end)
                }
            }
            (Some(start), _) => Some(start),
            (_, Some(end)) => Some(end),
            _ => None,
        }
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
    pre_align_min_score: i32,
) -> Option<i32> {
    // Try to the forward strand
    let banded_fwd = align_local_banded(
        query,
        &target_seq.fwd,
        banded_aligner,
        &target_hash.fwd_hash,
    );
    if banded_fwd >= pre_align_min_score {
        return Some(banded_fwd);
    }
    let banded_revcomp = align_local_banded(
        query,
        &target_seq.revcomp,
        banded_aligner,
        &target_hash.revcomp_hash,
    );
    if banded_revcomp >= pre_align_min_score {
        return Some(banded_revcomp);
    }
    None
}

fn header_to_name(header: &[u8]) -> Result<String> {
    let header: std::borrow::Cow<str> = String::from_utf8_lossy(header);
    header
        .split_whitespace()
        .next()
        .map(std::string::ToString::to_string)
        .context("empty read name")
}

pub fn to_records<F: MatchFunc>(
    fastq: &FastqOwnedRecord,
    result: Option<Alignment>,
    hard_clip: bool,
    use_eq_and_x: bool,
    alt_score: Option<i32>,
    scoring: &Scoring<F>,
    pick_primary: PrimaryPickingStrategy,
    target_seqs: &[TargetSeq],
) -> Result<Vec<SamRecord>> {
    let name = header_to_name(fastq.head())?;
    let read_name: SamReadName = name.parse()?;
    let bases = fastq.seq();
    let quals = fastq.qual();

    let mut builder = SubAlignmentBuilder::new(use_eq_and_x);

    if let Some(alignment) = result {
        let subs = builder.build(&alignment, true, scoring);
        ensure!(!subs.is_empty());

        let primary_index = match pick_primary {
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

        let records = subs
            .iter()
            .enumerate()
            .map(|(index, sub)| {
                let is_secondary = index != primary_index;
                let mut record = SamRecord::default();
                assert!(sub.contig_idx < 2 * target_seqs.len());
                let is_forward: bool = sub.contig_idx < target_seqs.len();

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
                *record.reference_sequence_id_mut() = Some(sub.contig_idx % target_seqs.len());

                // target start
                if is_forward {
                    *record.alignment_start_mut() = Position::new(sub.target_start + 1);
                } else {
                    let target_len = target_seqs[sub.contig_idx % target_seqs.len()].len();
                    *record.alignment_start_mut() = Position::new(target_len - sub.target_end + 1);
                }

                // mapping quality
                // TODO: base this on the alt_score
                *record.mapping_quality_mut() = MappingQuality::new(60);

                // TODO: tags (e.g. AS, XS, NM, MD)
                let mut data = Data::default();
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

                record
            })
            .collect_vec();

        Ok(records)
    } else {
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

        Ok([record].to_vec())
    }
}
