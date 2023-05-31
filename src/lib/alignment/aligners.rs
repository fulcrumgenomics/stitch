use anyhow::{ensure, Context, Result};
use bio::alignment::pairwise::MatchFunc;

use crate::alignment::sub_alignment::Cigar;
use crate::opts::{Opts, PrimaryPickingStrategy};
use crate::target_seq::{TargetHash, TargetSeq};
use crate::util::reverse_complement;
use itertools::{self, Itertools};
use std::str::FromStr;

use super::constants::{AlignmentMode, MIN_SCORE};
use super::pairwise::PairwiseAlignment;
use super::scoring::Scoring;
use super::sub_alignment::SubAlignmentBuilder;
use super::{double_strand::DoubleStrandAligner, single_strand::SingleStrandAligner};
use bio::alignment::pairwise::banded::Aligner as BandedAligner;
use bio::alignment::pairwise::MatchParams;
use bio::alignment::pairwise::Scoring as BioScoring;
use bio::alignment::sparse::HashMapFx as BandedHashMapFx;
use noodles::core::Position;
use noodles::sam::alignment::Record as SamRecord;
use noodles::sam::record::data::field::tag::ALIGNMENT_HIT_COUNT;
use noodles::sam::record::data::field::tag::ALIGNMENT_SCORE;
use noodles::sam::record::data::field::tag::HIT_INDEX;
use noodles::sam::record::data::field::tag::TOTAL_HIT_COUNT;
use noodles::sam::record::Cigar as SamCigar;
use noodles::sam::record::Data;
use noodles::sam::record::Flags;
use noodles::sam::record::MappingQuality;
use noodles::sam::record::QualityScores;
use noodles::sam::record::ReadName as SamReadName;
use noodles::sam::record::Sequence;
use seq_io::fastq::OwnedRecord as FastqOwnedRecord;
use seq_io::fastq::Record as FastqRecord;

pub struct Aligners<F: MatchFunc> {
    banded: BandedAligner<F>,
    single_strand: SingleStrandAligner<F>,
    double_strand: DoubleStrandAligner<F>,
}

impl Aligners<MatchParams> {
    fn build_match_fn(match_score: i32, mismatch_score: i32) -> MatchParams {
        MatchParams::new(match_score, mismatch_score)
    }

    fn build_stranded_scoring(
        opts: &Opts,
        xclip_prefix: i32,
        xclip_suffix: i32,
        yclip_prefix: i32,
        yclip_suffix: i32,
    ) -> Scoring<MatchParams> {
        let mut scoring = Scoring::new(
            opts.gap_open,
            opts.gap_extend,
            opts.jump_score,
            Self::build_match_fn(opts.match_score, opts.mismatch_score),
        );
        scoring.xclip_prefix = xclip_prefix;
        scoring.xclip_suffix = xclip_suffix;
        scoring.yclip_prefix = yclip_prefix;
        scoring.yclip_suffix = yclip_suffix;
        scoring
    }

    pub fn new(opts: &Opts, target_seq_len: usize) -> Aligners<MatchParams> {
        let (xclip_prefix, xclip_suffix, yclip_prefix, yclip_suffix) = match opts.mode {
            AlignmentMode::Local => (0, 0, 0, 0),
            AlignmentMode::Semiglobal => (MIN_SCORE, MIN_SCORE, 0, 0),
            AlignmentMode::Global => (MIN_SCORE, MIN_SCORE, MIN_SCORE, MIN_SCORE),
            AlignmentMode::Custom => panic!("Custom alignment mode not supported"), // TODO: move to main run method
        };

        // Banded alignment is always local since the goal is to find at leaset some minimal scoring
        // local alignment.
        let banded_scoring = {
            let mut scoring = BioScoring::new(
                opts.gap_open,
                opts.gap_extend,
                Self::build_match_fn(opts.match_score, opts.mismatch_score),
            );
            scoring.xclip_prefix = 0;
            scoring.xclip_suffix = 0;
            scoring.yclip_prefix = 0;
            scoring.yclip_suffix = 0;
            scoring
        };
        let banded = BandedAligner::with_capacity_and_scoring(
            10000,
            target_seq_len,
            banded_scoring,
            opts.k,
            opts.w,
        );

        let single_strand = SingleStrandAligner::with_capacity_and_scoring(
            10000,
            target_seq_len,
            Self::build_stranded_scoring(
                opts,
                xclip_prefix,
                xclip_suffix,
                yclip_prefix,
                yclip_suffix,
            ),
        );
        let double_strand: DoubleStrandAligner<MatchParams> = {
            let scoring_fwd = Self::build_stranded_scoring(
                &opts.clone(),
                xclip_prefix,
                xclip_suffix,
                yclip_prefix,
                yclip_suffix,
            );
            let scoring_rev = Self::build_stranded_scoring(
                &opts.clone(),
                xclip_prefix,
                xclip_suffix,
                yclip_prefix,
                yclip_suffix,
            );

            DoubleStrandAligner::with_capacity_and_scoring(
                10000,
                target_seq_len,
                scoring_fwd,
                scoring_rev,
            )
        };

        Self {
            banded,
            single_strand,
            double_strand,
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

fn maybe_prealign<F: MatchFunc>(
    query: &[u8],
    target_seq: &TargetSeq,
    target_hash: &TargetHash,
    banded_aligner: &mut BandedAligner<F>,
    pre_align: bool,
    pre_align_min_score: Option<i32>,
) -> (Option<i32>, Option<i32>, Option<i32>) {
    let min_score = pre_align_min_score.unwrap_or(MIN_SCORE);
    let (banded_fwd, banded_revcomp, prealign_score) = {
        if pre_align {
            let banded_fwd = align_local_banded(
                query,
                &target_seq.fwd,
                banded_aligner,
                &target_hash.fwd_hash,
            );
            if banded_fwd >= min_score {
                (Some(banded_fwd), None, Some(banded_fwd))
            } else {
                let banded_revcomp = align_local_banded(
                    query,
                    &target_seq.revcomp,
                    banded_aligner,
                    &target_hash.revcomp_hash,
                );
                let prealign_score = std::cmp::max(banded_fwd, banded_revcomp);
                (Some(banded_fwd), Some(banded_revcomp), Some(prealign_score))
            }
        } else {
            (None, None, None)
        }
    };

    (banded_fwd, banded_revcomp, prealign_score)
}

pub fn align_double_strand<F: MatchFunc>(
    record: FastqOwnedRecord,
    target_seq: &TargetSeq,
    target_hash: &TargetHash,
    aligners: &mut Aligners<F>,
    pre_align: bool,
    pre_align_min_score: i32,
) -> (
    FastqOwnedRecord,
    Option<(PairwiseAlignment, bool)>,
    Option<i32>,
) {
    let query = record.seq();

    let (banded_fwd, banded_revcomp, prealign_score) = maybe_prealign(
        query,
        target_seq,
        target_hash,
        &mut aligners.banded,
        pre_align,
        Some(pre_align_min_score),
    );

    let alignment = {
        if banded_fwd.map_or(true, |score| score >= pre_align_min_score)
            || banded_revcomp.map_or(true, |score| score >= pre_align_min_score)
        {
            Some(
                aligners
                    .double_strand
                    .custom(&target_seq.fwd, &target_seq.revcomp, query),
            )
        } else {
            None
        }
    };

    match alignment {
        // NB: we always align to the "forward strand"
        Some(result) => (record, Some((result, true)), prealign_score),
        None => (record, None, prealign_score),
    }
}

pub fn align_single_strand<F: MatchFunc>(
    record: FastqOwnedRecord,
    target_seq: &TargetSeq,
    target_hash: &TargetHash,
    aligners: &mut Aligners<F>,
    pre_align: bool,
    pre_align_min_score: i32,
) -> (
    FastqOwnedRecord,
    Option<(PairwiseAlignment, bool)>,
    Option<i32>,
) {
    let query = record.seq();
    let (banded_fwd, banded_revcomp, prealign_score) = maybe_prealign(
        query,
        target_seq,
        target_hash,
        &mut aligners.banded,
        pre_align,
        None,
    );

    let fwd: Option<PairwiseAlignment> =
        if banded_fwd.map_or(true, |score| score >= pre_align_min_score) {
            Some(aligners.single_strand.custom(&target_seq.fwd, query))
        } else {
            None
        };
    let revcomp = if banded_revcomp.map_or(true, |score| score >= pre_align_min_score) {
        Some(aligners.single_strand.custom(&target_seq.revcomp, query))
    } else {
        None
    };

    match (fwd, revcomp) {
        (None, None) => (record, None, prealign_score),
        (None, Some(r)) => (record, Some((r, false)), prealign_score),
        (Some(f), None) => (record, Some((f, true)), prealign_score),
        (Some(f), Some(r)) => {
            if f.score >= r.score {
                (record, Some((f, true)), prealign_score)
            } else {
                (record, Some((r, false)), prealign_score)
            }
        }
    }
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
    result: Option<(PairwiseAlignment, bool)>,
    hard_clip: bool,
    use_eq_and_x: bool,
    alt_score: Option<i32>,
    scoring: &Scoring<F>,
    by_primary: PrimaryPickingStrategy,
    target_len: usize,
) -> Result<Vec<SamRecord>> {
    let name = header_to_name(fastq.head())?;
    let read_name: SamReadName = name.parse()?;
    let bases = fastq.seq();
    let quals = fastq.qual();

    let mut builder = SubAlignmentBuilder::new(use_eq_and_x);

    if let Some((alignment, is_fwd)) = result {
        let subs = builder.build(&alignment, true, scoring);
        ensure!(!subs.is_empty());

        let primary_index = match by_primary {
            PrimaryPickingStrategy::QueryLength => subs
                .iter()
                .enumerate()
                .max_by_key(|(_, alignment)| alignment.score)
                .map_or(0, |(index, _)| index),
            PrimaryPickingStrategy::Score => subs
                .iter()
                .enumerate()
                .max_by_key(|(_, alignment)| alignment.query_end - alignment.query_start)
                .map_or(0, |(index, _)| index),
        };

        let records = subs
            .iter()
            .enumerate()
            .map(|(index, sub)| {
                let is_secondary = index != primary_index;
                let mut record = SamRecord::default();
                let is_forward: bool = is_fwd == sub.is_forward;

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
                        Cigar {
                            elements: sub.cigar.elements.iter().rev().copied().collect(),
                        },
                    ),
                    (false, false) => (
                        reverse_complement(bases),
                        quals.iter().copied().rev().collect(),
                        Cigar {
                            elements: sub.cigar.elements.iter().rev().copied().collect(),
                        },
                    ),
                    (false, true) => (
                        reverse_complement(bases[sub.query_start..sub.query_end].to_vec()),
                        quals[sub.query_start..sub.query_end]
                            .iter()
                            .copied()
                            .rev()
                            .collect(),
                        Cigar {
                            elements: sub.cigar.elements.iter().rev().copied().collect(),
                        },
                    ),
                };

                // bases
                *record.sequence_mut() = Sequence::try_from(bases_vec).unwrap();

                // qualities
                *record.quality_scores_mut() = QualityScores::try_from(quals_vec).unwrap();

                // cigar
                let clip_str: &str = if hard_clip && is_secondary { "H" } else { "S" };
                let mut cigar_str = String::new();
                if is_forward && sub.query_start > 0 {
                    cigar_str.push_str(&format!("{}{}", sub.query_start, clip_str));
                } else if !is_forward && sub.query_end < bases.len() {
                    cigar_str.push_str(&format!("{}{}", bases.len() - sub.query_end, clip_str));
                }
                for elem in &cigar.elements {
                    cigar_str.push_str(&format!("{}{}", elem.len, elem.op));
                }
                if is_forward && sub.query_end < bases.len() {
                    cigar_str.push_str(&format!("{}{}", bases.len() - sub.query_end, clip_str));
                } else if !is_forward && sub.query_start > 0 {
                    cigar_str.push_str(&format!("{}{}", sub.query_start, clip_str));
                }
                *record.cigar_mut() = SamCigar::from_str(&cigar_str).unwrap();

                // target id
                *record.reference_sequence_id_mut() = Some(0usize);

                // target start
                if is_forward {
                    *record.alignment_start_mut() = Position::new(sub.target_start + 1);
                } else {
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
        *record.cigar_mut() = SamCigar::default();

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
