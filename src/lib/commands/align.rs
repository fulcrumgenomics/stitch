use crate::align::aligners::align_double_strand;
use crate::align::aligners::align_single_strand;
use crate::align::aligners::constants::AlignmentMode;
use crate::align::aligners::to_records;
use crate::align::aligners::Aligners;
use crate::align::io::FastqGroupingIterator;
use crate::align::io::FastqThreadReader;
use crate::align::io::OutputMessage;
use crate::align::io::OutputResult;
use crate::align::io::BUFFER_SIZE;
use crate::align::io::READER_CHANNEL_NUM_CHUNKS;
use crate::align::scoring::Scoring;
use crate::align::PrimaryPickingStrategy;
use crate::util::target_seq::TargetHash;
use crate::util::target_seq::TargetSeq;
use crate::util::version::built_info;
use crate::util::version::built_info::VERSION;
use anyhow::Context;
use anyhow::{ensure, Result};
use bio::alignment::pairwise::MatchParams;
use clap::Parser;
use fgoxide::io::Io;
use flume::unbounded;
use itertools::{self, Itertools};
use log::info;
use noodles::bam::Writer as BamWriter;
use noodles::bgzf;
use noodles::bgzf::writer::CompressionLevel;
use noodles::sam::header::record::value::map::Program;
use noodles::sam::header::record::value::map::ReferenceSequence;
use noodles::sam::header::record::value::Map;
use noodles::sam::header::Header as SamHeader;
use proglog::CountFormatterKind;
use proglog::ProgLogBuilder;
use seq_io::fasta::Reader as FastaReader;
use seq_io::fasta::Record as FastaRecord;
use seq_io::fastq::OwnedRecord as FastqOwnedRecord;
use std::env;
use std::io;
use std::io::BufRead;
use std::io::Write;
use std::num::NonZeroUsize;
use std::path::PathBuf;
use std::thread::JoinHandle;
use std::time::Duration;

/// Converts the FASTQ header (which may contain whitespaces) to a QNAME for the SAM format.
fn header_to_name(header: &[u8]) -> Result<String> {
    let header: std::borrow::Cow<str> = String::from_utf8_lossy(header);
    header
        .split_whitespace()
        .next()
        .map(std::string::ToString::to_string)
        .context("empty read name")
}

/// Reads a FASTA containing a single target contig for the vector/plasmid/construct against which
/// to align.
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

////////////////////////////////////////////////////////////////////////////////
// Align (main class) and it's impls
////////////////////////////////////////////////////////////////////////////////

/// Perfoms alignment of a long reads against a reference/expected vector/plasmid/construct.
///
/// The alignment extends the traditional alignment algorithms by introducing a "jump"
/// move/operator, whereby the alignment is able to jump anywhere in the reference sequence for a
/// fixed cost.  This models the case where the read sequence is composed of multiple pieces of the
/// expected reference, in shuffled or truncated order.  In the default mode, the alignment is
/// allowed to jump on the same strand either before or after the current reference position,
/// where as with `--double-strand` the alignment is also able to jump and continue on the opposite
/// strand and continue on the opposite strand.  Thus, in both modes it is possible to align to the
/// same reference sequence multiple times.
///
/// The output is in SAM/BAM format.  If the alignment has `N` jumps, then the output will contain
/// `N-1` records for the input read.  One record is marked as primary (see `--pick-primary`),
/// while the remaining records are marked as secondary.  The HI/HN SAM tags are used to denote
/// the order in which the alignments occur.  Furthermore, the order of the records output by this
/// tool are in the order in which they align the query/read sequence.
///
/// If the input FASTQ may be sorted by read sequence (bases), then the aligner will be sped up,
/// since it need only read align the first read in a run of consecutive reads that have the same
/// read sequence.
#[derive(Parser, Debug, Clone)]
#[clap(version = built_info::VERSION.as_str(), term_width=0)]
pub struct Align {
    /// The path to the input FASTQ with sequenced vector/plasmid/construct long reads.
    #[clap(long, short = 'f', display_order = 1)]
    pub reads_fastq: PathBuf,

    /// The path to the referece vector/plasmid/construct FASTA sequence.
    #[clap(long, short = 'r', display_order = 2)]
    pub ref_fasta: PathBuf,

    /// Align to both strands of the reference simulataneously.
    #[clap(long, short = 'd', default_value = "false", display_order = 3)]
    pub double_strand: bool,

    /// The number of threads to use.
    #[clap(long, short = 't', default_value = "2", display_order = 4)]
    pub threads: usize,

    /// Assume an unrecognized input (based on file extension) is GZIP compressed.
    #[clap(long, short = 'z', default_value = "false", display_order = 5)]
    pub decompress: bool,

    /// Pre-align with banded local alignment.
    #[clap(long, short = 'p', default_value = "false", display_order = 6)]
    pub pre_align: bool,

    /// K-mer size for banded pre-alignment.
    #[clap(long, short = 'k', default_value = "12", display_order = 7)]
    pub k: usize,

    /// Band size for banded pre-alignment.
    #[clap(long, short = 'w', default_value = "50", display_order = 8)]
    pub w: usize,

    /// Pre-align with banded local alignment.
    #[clap(
        long,
        short = 's',
        default_value = "100",
        allow_hyphen_values = true,
        display_order = 9
    )]
    pub pre_align_min_score: i32,

    /// Use soft-clipping for all alignments, otherwise secondary alignemnts will use hard-clipping
    #[clap(long, short = 'S', default_value = "false", display_order = 10)]
    pub soft_clip: bool,

    /// Use =/X CIGAR operators, otherwise use M
    #[clap(long, short = 'X', default_value = "false", display_order = 11)]
    pub use_eq_and_x: bool,

    /// Score for a sequence match (must be positive)
    #[clap(long, short = 'A', default_value = "1", display_order = 12)]
    pub match_score: i32,

    /// Score for a sequence mismatch (must be negative)
    #[clap(
        long,
        short = 'B',
        default_value = "-1",
        allow_hyphen_values = true,
        display_order = 13
    )]
    pub mismatch_score: i32,

    /// Score for a gap open (must be negative)
    #[clap(
        long,
        short = 'O',
        default_value = "-5",
        allow_hyphen_values = true,
        display_order = 14
    )]
    pub gap_open: i32,

    /// Score for a gap extend (must be negative); a gap of size k cost '{-O} + {-E}*k'
    #[clap(
        long,
        short = 'E',
        default_value = "-1",
        allow_hyphen_values = true,
        display_order = 15
    )]
    pub gap_extend: i32,

    /// Score for a target jump (must be negative)
    #[clap(
        long,
        short = 'J',
        default_value = "-10",
        allow_hyphen_values = true,
        display_order = 16
    )]
    pub jump_score: i32,

    /// The alignment mode:
    /// - Local: aligns a sub-sequence of the read versus a sub-sequence of the reference.
    /// - QueryLocal: aligns a sub-sequence of the read versus the full reference.
    /// - TargetLocal: aligns the full read versus a sub-sequence of the reference.
    /// - Global: aligns the full read versus the full reference.
    #[clap(
        long,
        short = 'm',
        default_value = "local",
        display_order = 17,
        verbatim_doc_comment
    )]
    pub mode: AlignmentMode,

    /// Determines how to pick the primary alignment:
    /// - query-length: the longest aligned query length, then by `score`
    /// - score: the best scoring alignment, then by `query-length`
    #[clap(
        long,
        short = 'P',
        default_value = "query-length",
        display_order = 18,
        verbatim_doc_comment
    )]
    pub pick_primary: PrimaryPickingStrategy,

    /// True to treat the input target as circular.  This allows the alignment to jump back to the
    /// start of the target/reference at no cost.
    #[clap(long, short = 'C', default_value = "false", display_order = 19)]
    pub circular: bool,

    /// The compression level of the output BAM
    #[clap(long, short = 'c', default_value = "0", display_order = 20)]
    pub compression: u8,
}

impl Align {
    /// Aligns a chunk of records
    fn align(
        records: Vec<FastqOwnedRecord>,
        opts: &Align,
        target_seq: &TargetSeq,
        target_hash: &TargetHash,
        aligners: &mut Aligners<MatchParams>,
    ) -> Vec<OutputResult> {
        let iter = FastqGroupingIterator::new(records.into_iter());
        let mut results: Vec<OutputResult> = Vec::new();

        for group in iter {
            let first: &FastqOwnedRecord = group.first().unwrap();
            let (maybe_alignment_and_is_forward, maybe_score) = {
                if opts.double_strand {
                    align_double_strand(
                        first,
                        target_seq,
                        target_hash,
                        aligners,
                        opts.pre_align,
                        opts.pre_align_min_score,
                    )
                } else {
                    align_single_strand(
                        first,
                        target_seq,
                        target_hash,
                        aligners,
                        opts.pre_align,
                        opts.pre_align_min_score,
                    )
                }
            };
            match (maybe_alignment_and_is_forward, maybe_score) {
                (Some((alignment, is_forward)), Some(score)) => {
                    for record in group {
                        let alignment = alignment.clone();
                        results.push((record, Some((alignment, is_forward)), Some(score)));
                    }
                }
                (None, None) => {
                    for record in group {
                        results.push((record, None, None));
                    }
                }
                (None, Some(score)) => {
                    for record in group {
                        results.push((record, None, Some(score)));
                    }
                }
                _ => panic!("Bug: should not reach here"),
            };
        }

        results
    }

    /// Executes the align command
    pub fn execute(&self) -> anyhow::Result<()> {
        info!("Starting alignment...");
        info!("Reading reference FASTA from {}", self.ref_fasta.display());
        info!("Reading reads FASTQ from {}", self.reads_fastq.display());
        // ensure!(self.threads > 1, "Must specify at least two threads");
        let progress_logger = ProgLogBuilder::new()
            .name("fqcv-progress")
            .noun("reads")
            .verb("Aligned")
            .unit(
                (READER_CHANNEL_NUM_CHUNKS * self.threads)
                    .try_into()
                    .unwrap(),
            )
            .count_formatter(CountFormatterKind::Comma)
            .build();

        // Read in the refefence/target FASTA
        let (target_seq, target_name) = read_target(&self.ref_fasta)?;

        // Create the thread to read in the FASTQ records
        let reader =
            FastqThreadReader::new(self.reads_fastq.clone(), self.decompress, self.threads);

        // Create the channel to gracefully signal a shutdown of the aligner threads
        let (shutdown_tx, shutdown_rx) = unbounded::<()>();

        // Create and start the aligner threads
        let sleep_delay = Duration::from_millis(100);
        let thread_handles: Vec<JoinHandle<Result<()>>> = (0..self.threads)
            .map(|_| {
                let to_align_rx = reader.to_align_rx.clone();
                let shutdown_rx = shutdown_rx.clone();
                let target_name = target_name.clone();
                let target_seq = target_seq.clone();
                let mut aligners = Aligners::new(self, target_seq.len());
                let opts = self.clone();

                std::thread::spawn(move || {
                    let target_seq = TargetSeq::new(&target_name, &target_seq);
                    let target_hash = target_seq.build_target_hash(opts.k);
                    loop {
                        // Try to process one chunk of alignments
                        if let Ok(msg) = to_align_rx.try_recv() {
                            let results = Self::align(
                                msg.records,
                                &opts,
                                &target_seq,
                                &target_hash,
                                &mut aligners,
                            );
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

        // Setup and write the SAM header
        let command_line = env::args_os().map(|s| s.into_string().unwrap()).join(" ");
        let stdout = io::stdout().lock();
        let encoder = bgzf::writer::Builder::default()
            .set_compression_level(CompressionLevel::try_from(self.compression)?)
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

        // Convert the alignments to SAM records
        loop {
            let match_fn: MatchParams = MatchParams::new(self.match_score, self.mismatch_score);
            let scoring = Scoring::new(self.gap_open, self.gap_extend, self.jump_score, match_fn);
            // Get a receiver for the alignment of a record
            if let Ok(receiver) = reader.to_output_rx.try_recv() {
                let msg = receiver.recv()?;
                for (fastq, result, alt_score) in msg.results {
                    progress_logger.record();
                    let records = to_records(
                        &fastq,
                        result,
                        !self.soft_clip,
                        self.use_eq_and_x,
                        alt_score,
                        &scoring,
                        self.pick_primary,
                        target_seq.len(),
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
}
