use crate::alignment::PairwiseAlignment;
use crate::alignment::SubAlignment;
use crate::io::FastqThreadReader;
use crate::io::OutputMessage;
use crate::io::OutputResult;
use crate::io::BUFFER_SIZE;
use crate::opts::Opts;
use crate::pairwise::Aligner as PairwiseAligner;
use crate::pairwise::Scoring;
use crate::util::built_info::VERSION;
use crate::util::reverse_complement;
use anyhow::Context;
use anyhow::{ensure, Result};
use bio::alignment::pairwise::banded::Aligner as BandedAligner;
use bio::alignment::pairwise::Scoring as BioScoring;
use bio::alignment::pairwise::MIN_SCORE;
use bio::alignment::sparse::hash_kmers;
use bio::alignment::sparse::HashMapFx as BandedHashMapFx;
use fgoxide::io::Io;
use flume::unbounded;
use itertools::{self, Itertools};
use noodles::bam::Writer as BamWriter;
use noodles::bgzf;
use noodles::bgzf::writer::CompressionLevel;
use noodles::core::Position;
use noodles::sam::alignment::Record as SamRecord;
use noodles::sam::header::record::value::map::Program;
use noodles::sam::header::record::value::map::ReferenceSequence;
use noodles::sam::header::record::value::Map;
use noodles::sam::header::Header as SamHeader;
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
use proglog::CountFormatterKind;
use proglog::ProgLogBuilder;
use seq_io::fasta::Reader as FastaReader;
use seq_io::fasta::Record as FastaRecord;
use seq_io::fastq::OwnedRecord as FastqOwnedRecord;
use seq_io::fastq::Record as FastqRecord;
use std::env;
use std::io;
use std::io::BufRead;
use std::io::Write;
use std::num::NonZeroUsize;
use std::path::PathBuf;
use std::str::FromStr;
use std::thread::JoinHandle;
use std::time::Duration;

fn align_local_banded(
    query: &[u8],
    target: &[u8],
    aligner: &mut BandedAligner<impl Fn(u8, u8) -> i32>,
    target_kmer_hash: &BandedHashMapFx<&[u8], Vec<u32>>,
) -> i32 {
    // Compute the alignment
    aligner
        .custom_with_prehash(query, target, target_kmer_hash)
        .score
}

fn align(
    record: FastqOwnedRecord,
    target_seq: &[u8],
    pairwise_aligner: &mut PairwiseAligner<impl Fn(u8, u8) -> i32>,
    banded_aligner: &mut BandedAligner<impl Fn(u8, u8) -> i32>,
    target_kmer_hash: &BandedHashMapFx<&[u8], Vec<u32>>,
    pre_align: bool,
    pre_align_min_score: i32,
) -> (
    FastqOwnedRecord,
    Option<(PairwiseAlignment, bool)>,
    Option<i32>,
) {
    let query = record.seq();
    let query_revcomp = reverse_complement(query);

    let banded_fwd = if pre_align {
        let score = align_local_banded(query, target_seq, banded_aligner, target_kmer_hash);
        Some(score)
    } else {
        None
    };

    let banded_revcomp = if pre_align {
        let score =
            align_local_banded(&query_revcomp, target_seq, banded_aligner, target_kmer_hash);
        Some(score)
    } else {
        None
    };

    let fwd = if banded_fwd.map_or(true, |score| score >= pre_align_min_score) {
        // Some(pairwise_aligner.local(target_seq, query))
        Some(pairwise_aligner.semiglobal2(target_seq, query))
    } else {
        None
    };
    // let fwd: Option<PairwiseAlignment> = None;
    let revcomp = if banded_revcomp.map_or(true, |score| score >= pre_align_min_score) {
        // Some(pairwise_aligner.local(target_seq, &query_revcomp))
        Some(pairwise_aligner.semiglobal2(target_seq, &query_revcomp))
    } else {
        None
    };
    // let revcomp: Option<PairwiseAlignment> = None;

    let prealign_score = match (banded_fwd, banded_revcomp) {
        (None, None) => None,
        (Some(f), None) => {
            if f > MIN_SCORE {
                Some(f)
            } else {
                None
            }
        }
        (None, Some(r)) => {
            if r > MIN_SCORE {
                Some(r)
            } else {
                None
            }
        }
        (Some(f), Some(r)) => {
            let score = std::cmp::max(f, r);
            if score > MIN_SCORE {
                Some(score)
            } else {
                None
            }
        }
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

fn read_target(file: &PathBuf) -> Result<(Vec<u8>, String)> {
    let fg_io: Io = Io::new(5, BUFFER_SIZE);
    let source: FastaReader<Box<dyn BufRead + Send>> =
        FastaReader::with_capacity(fg_io.new_reader(file)?, BUFFER_SIZE);

    let sequences = source
        .into_records()
        .map(|r| r.unwrap_or_else(|_| panic!("Error reading FASTA")))
        .collect_vec();
    ensure!(!sequences.is_empty(), "Found no sequences in the FASTA");
    ensure!(
        sequences.len() == 1,
        "Found multiple sequences in the FASTA"
    );

    let record = sequences
        .first()
        .context("No sequences found in the FASTA")?;

    let sequence = record
        .seq()
        .to_owned()
        .iter()
        .map(u8::to_ascii_uppercase)
        .collect_vec();

    let name = header_to_name(record.head())?;
    Ok((sequence, name))
}

fn to_records(
    fastq: &FastqOwnedRecord,
    result: Option<(PairwiseAlignment, bool)>,
    hard_clip: bool,
    use_eq_and_x: bool,
    alt_score: Option<i32>,
    scoring: Scoring<impl Fn(u8, u8) -> i32>,
) -> Result<Vec<SamRecord>> {
    let name = header_to_name(fastq.head())?;
    let read_name: SamReadName = name.parse()?;
    let bases = fastq.seq();
    let quals = fastq.qual();

    if let Some((alignment, is_fwd)) = result {
        let subs = SubAlignment::build(&alignment, scoring, use_eq_and_x, true);
        let bases = if is_fwd {
            bases.to_vec()
        } else {
            let x = reverse_complement(bases);
            assert!(x != bases);
            reverse_complement(bases)
        };
        let quals = if is_fwd {
            quals.to_vec()
        } else {
            let mut q = quals.to_vec();
            q.reverse();
            q
        };

        ensure!(!subs.is_empty());

        let primary_index = subs
            .iter()
            .enumerate()
            .max_by_key(|(_, alignment)| alignment.query_end - alignment.query_start)
            .map_or(0, |(index, _)| index);

        let records = subs
            .iter()
            .enumerate()
            .map(|(index, sub)| {
                let is_secondary = index != primary_index;
                let mut record = SamRecord::default();

                // read name
                *record.read_name_mut() = Some(read_name.clone());

                // flags
                let mut new_flags = Flags::default();
                if !is_fwd {
                    new_flags.insert(Flags::REVERSE_COMPLEMENTED);
                }
                if is_secondary {
                    new_flags.insert(Flags::SECONDARY);
                }
                *record.flags_mut() = new_flags;

                if hard_clip && is_secondary {
                    // bases
                    let bases_vec = bases[sub.query_start..sub.query_end].to_vec();
                    *record.sequence_mut() = Sequence::try_from(bases_vec).unwrap();

                    // qualities
                    let quals_vec = quals[sub.query_start..sub.query_end].to_vec();
                    *record.quality_scores_mut() = QualityScores::try_from(quals_vec).unwrap();
                } else {
                    // bases
                    *record.sequence_mut() = Sequence::try_from(bases.clone()).unwrap();

                    // qualities
                    *record.quality_scores_mut() = QualityScores::try_from(quals.clone()).unwrap();
                }

                // cigar
                let clip_str: &str = if hard_clip && is_secondary { "H" } else { "S" };
                let mut cigar = String::new();
                if sub.query_start > 0 {
                    cigar.push_str(&format!("{}{}", sub.query_start, clip_str));
                }
                cigar.push_str(&sub.cigar);
                if sub.query_end < bases.len() {
                    cigar.push_str(&format!("{}{}", bases.len() - sub.query_end, clip_str));
                }
                *record.cigar_mut() = Cigar::from_str(&cigar).unwrap();

                // target id
                *record.reference_sequence_id_mut() = Some(0usize);

                // target start
                *record.alignment_start_mut() = Position::new(sub.target_start + 1);

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

pub fn run(opts: &Opts) -> Result<()> {
    // ensure!(opts.threads > 1, "Must specify at least two threads");

    let progress_logger = ProgLogBuilder::new()
        .name("fqcv-progress")
        .noun("reads")
        .verb("Aligned")
        .unit((125 * opts.threads).try_into().unwrap())
        .count_formatter(CountFormatterKind::Comma)
        .build();

    // Read in the refefence/target FASTA
    let (target_seq, target_name) = read_target(&opts.ref_fasta)?;

    // Create the thread to read in the FASTQ records
    let reader = FastqThreadReader::new(opts.reads_fastq.clone(), opts.decompress, opts.threads);

    // Create the channel to gracefully signal a shutdown of the aligner threads
    let (shutdown_tx, shutdown_rx) = unbounded::<()>();

    let sleep_delay = Duration::from_millis(25);

    let match_fn = {
        let match_score: i32 = opts.match_score;
        let mismatch_score: i32 = opts.mismatch_score;
        move |a: u8, b: u8| {
            if a == b {
                match_score
            } else {
                mismatch_score
            }
        }
    };
    let scoring = Scoring::new(opts.gap_open, opts.gap_extend, opts.jump_score, match_fn);
    let scoring_bio = BioScoring::new(opts.gap_open, opts.gap_extend, match_fn);

    let thread_handles: Vec<JoinHandle<Result<()>>> = (0..opts.threads)
        .map(|_| {
            let to_align_rx = reader.to_align_rx.clone();
            let shutdown_rx = shutdown_rx.clone();
            let target_seq = target_seq.clone();
            let mut pairwise_aligner =
                PairwiseAligner::with_capacity_and_scoring(10000, target_seq.len(), scoring);
            let opts = opts.clone();
            let mut banded_aligner = {
                let mut aligner = BandedAligner::with_capacity_and_scoring(
                    10000,
                    target_seq.len(),
                    scoring_bio,
                    opts.k,
                    opts.w,
                );
                let scoring: &mut bio::alignment::pairwise::Scoring<_> = aligner.get_mut_scoring();

                // Temporarily Over-write the clip penalties
                scoring.xclip_prefix = 0;
                scoring.xclip_suffix = 0;
                scoring.yclip_prefix = 0;
                scoring.yclip_suffix = 0;
                aligner
            };

            std::thread::spawn(move || {
                let target_kmers_hash = hash_kmers(&target_seq, opts.k);
                loop {
                    // Try to process one chunk of alignments
                    if let Ok(msg) = to_align_rx.try_recv() {
                        let results: Vec<OutputResult> = msg
                            .records
                            .into_iter()
                            .map(|record| {
                                align(
                                    record,
                                    &target_seq,
                                    &mut pairwise_aligner,
                                    &mut banded_aligner,
                                    &target_kmers_hash,
                                    opts.pre_align,
                                    opts.pre_align_min_score,
                                )
                            })
                            .collect_vec();
                        msg.oneshot
                            .send(OutputMessage { results })
                            .expect("Send failed");
                    } else {
                        if shutdown_rx.is_disconnected() && to_align_rx.is_empty() {
                            break;
                        }
                        std::thread::sleep(sleep_delay);
                    }
                }
                Ok(())
            })
        })
        .collect();

    let command_line = env::args_os().map(|s| s.into_string().unwrap()).join(" ");
    let stdout = io::stdout().lock();
    let encoder = bgzf::writer::Builder::default()
        .set_compression_level(CompressionLevel::try_from(opts.compression)?)
        .build_with_writer(stdout);
    let mut writer = BamWriter::from(encoder);
    let header = SamHeader::builder()
        .set_header(Map::default())
        .add_program(
            "fqcv",
            Map::<Program>::builder()
                .set_name("fqcv")
                .set_version(VERSION.clone())
                .set_command_line(command_line)
                .build()?,
        )
        .add_reference_sequence(
            target_name.parse()?,
            Map::<ReferenceSequence>::new(NonZeroUsize::try_from(target_seq.len())?),
        )
        .build();
    writer.write_header(&header)?;

    loop {
        // Get a receiver for the alignment of a record
        if let Ok(receiver) = reader.to_output_rx.try_recv() {
            let msg = receiver.recv()?;
            for (fastq, result, alt_score) in msg.results {
                progress_logger.record();
                let records = to_records(
                    &fastq,
                    result,
                    !opts.soft_clip,
                    opts.use_eq_and_x,
                    alt_score,
                    scoring,
                )?;
                for record in &records {
                    writer.write_record(&header, record)?;
                }
                io::stdout().flush()?;
            }
        } else {
            if reader.handle.is_finished()
                && reader.to_align_rx.is_empty()
                && reader.to_output_rx.is_empty()
            {
                break;
            }
            std::thread::sleep(sleep_delay);
        }
    }

    // All done, shut down the reader and alignment threads
    match reader.handle.join() {
        Ok(_) => (),
        Err(e) => std::panic::resume_unwind(e),
    };
    drop(shutdown_tx); // to signal the alignment threads
    thread_handles
        .into_iter()
        .try_for_each(|handle| match handle.join() {
            Ok(result) => result,
            Err(e) => std::panic::resume_unwind(e),
        })?;

    Ok(())
}
