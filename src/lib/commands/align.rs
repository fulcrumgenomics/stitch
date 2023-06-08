use crate::align::aligners::constants::AlignmentMode;
use crate::align::aligners::to_records;
use crate::align::io::FastqGroupingIterator;
use crate::align::io::FastqThreadReader;
use crate::align::io::OutputMessage;
use crate::align::io::OutputResult;
use crate::align::io::READER_CHANNEL_NUM_CHUNKS;
use crate::align::scoring::Scoring;
use crate::align::PrimaryPickingStrategy;
use crate::util::target_seq::from_fasta;
use crate::util::target_seq::TargetHash;
use crate::util::version::built_info;
use crate::util::version::built_info::VERSION;
use anyhow::Result;
use bio::alignment::pairwise::MatchParams;
use clap::Parser;
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
use seq_io::fastq::OwnedRecord as FastqOwnedRecord;
use std::env;
use std::io;
use std::io::Write;
use std::num::NonZeroUsize;
use std::path::PathBuf;
use std::thread::JoinHandle;
use std::time::Duration;

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
/// If the input FASTQ may be sorted by read sequence (bases), in which case the aligner will be
/// sped up, since it need only read align the first read in a run of consecutive reads that have
/// the same read sequence.
///
/// Multiple contigs in the input FASTA are supported.
///
/// ## Jump scores
///
/// The jump score can be specified with `--jump-score`.
///
/// The jump score may also be specified specific to the jump being within the same contig
/// and the same strand (`--jump-score-same-contig-and-strand`), the same contig but opposite
/// strand (`--jump-score-same-contig-opposite-strand`), and the across different contigs.
/// (`--jump-score-inter-contig`).  If any of these options are not specified, then they will
/// default to the the value specified by `--jump-score`.
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

    /// The minimum score required for pre-alignment.
    #[clap(
        long,
        short = 's',
        default_value = "100",
        allow_hyphen_values = true,
        display_order = 9
    )]
    pub pre_align_min_score: i32,

    /// Only align to contigs that had an alignment score greater than the `--pre-align-min-score`
    #[clap(long, short = 'x', default_value = "true", display_order = 8)]
    pub pre_align_subset_contigs: bool,

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
        default_value = "-4",
        allow_hyphen_values = true,
        display_order = 13
    )]
    pub mismatch_score: i32,

    /// Score for a gap open (must be negative)
    #[clap(
        long,
        short = 'O',
        default_value = "-6",
        allow_hyphen_values = true,
        display_order = 14
    )]
    pub gap_open: i32,

    /// Score for a gap extend (must be negative); a gap of size k cost '{-O} + {-E}*k'
    #[clap(
        long,
        short = 'E',
        default_value = "-2",
        allow_hyphen_values = true,
        display_order = 15
    )]
    pub gap_extend: i32,

    /// Score for a target jump (must be negative).
    #[clap(
        long,
        short = 'J',
        default_value = "-10",
        allow_hyphen_values = true,
        display_order = 16
    )]
    pub jump_score: i32,

    /// Score for a target jump within the same contig and strand (must be negative)
    #[clap(long, allow_hyphen_values = true, display_order = 16)]
    pub jump_score_same_contig_and_strand: Option<i32>,

    /// Score for a target jump within the same contig and oppoosite strand (must be negative)
    #[clap(long, allow_hyphen_values = true, display_order = 16)]
    pub jump_score_same_contig_opposite_strand: Option<i32>,

    /// Score for a target jump across different contigs (must be negative)
    #[clap(long, allow_hyphen_values = true, display_order = 16)]
    pub jump_score_inter_contig: Option<i32>,

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

    /// When using `--circular`, re-align reads that may cross the target origin.  This may occur
    /// when a prefix (or suffix) of the read is unaligned, and the alignment starts within the
    /// slop of the target start (or target end).  In this case, the prefix is appended to the
    /// aligned suffix and the new read is re-aligned.
    #[clap(long, short = 'C', default_value = "20", display_order = 19)]
    pub circular_slop: usize,

    /// The compression level of the output BAM
    #[clap(long, short = 'c', default_value = "0", display_order = 21)]
    pub compression: u8,
}

impl Align {
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
        let target_seqs = from_fasta(&self.ref_fasta, self.circular)?;

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
                let opts = self.clone();
                let target_seqs = target_seqs.clone();

                std::thread::spawn(move || {
                    let target_hashes: Vec<TargetHash> = target_seqs
                        .iter()
                        .map(|target_seq| target_seq.build_target_hash(opts.k))
                        .collect();
                    let mut aligners = crate::align::aligners::build_aligners(&opts, &target_seqs);
                    loop {
                        // Try to process one chunk of alignments
                        if let Ok(msg) = to_align_rx.try_recv() {
                            let iter = FastqGroupingIterator::new(msg.records.into_iter());
                            let mut results: Vec<OutputResult> = Vec::new();
                            for group in iter {
                                let first: &FastqOwnedRecord = group.first().unwrap();
                                let (maybe_alignment, maybe_score) =
                                    aligners.align(first, &target_seqs, &target_hashes, &opts);
                                match maybe_alignment {
                                    Some(alignment) => {
                                        for record in group {
                                            let alignment = alignment.clone();
                                            results.push((record, Some(alignment), maybe_score));
                                        }
                                    }
                                    None => {
                                        for record in group {
                                            results.push((record, None, maybe_score));
                                        }
                                    }
                                };
                            }

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
        let header = {
            let mut builder = SamHeader::builder().set_header(Map::default()).add_program(
                "fqcv",
                Map::<Program>::builder()
                    .set_name("fqcv")
                    .set_version(VERSION.clone())
                    .set_command_line(command_line)
                    .build()?,
            );
            for target_seq in &target_seqs {
                builder = builder.add_reference_sequence(
                    target_seq.name.parse()?,
                    Map::<ReferenceSequence>::new(NonZeroUsize::try_from(target_seq.len())?),
                );
            }
            builder.build()
        };
        writer.write_header(&header)?;

        // Convert the alignments to SAM records
        loop {
            let match_fn: MatchParams = MatchParams::new(self.match_score, self.mismatch_score);
            let scoring = Scoring::with_jump_scores(
                self.gap_open,
                self.gap_extend,
                self.jump_score_same_contig_and_strand
                    .unwrap_or(self.jump_score),
                self.jump_score_same_contig_opposite_strand
                    .unwrap_or(self.jump_score),
                self.jump_score_inter_contig.unwrap_or(self.jump_score),
                match_fn,
            );
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
                        &target_seqs,
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
