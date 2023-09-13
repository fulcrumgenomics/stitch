use super::command::{Command, ValueEnum};
use anyhow::Result;
use clap::{
    builder::{PossibleValuesParser, TypedValueParser as _},
    Parser,
};
use flume::unbounded;
use itertools::{self, Itertools};
use log::info;
use noodles::{
    bam::Writer as BamWriter,
    bgzf,
    bgzf::writer::CompressionLevel,
    sam::header::{
        record::value::{
            map::{Program, ReferenceSequence},
            Map,
        },
        Header as SamHeader,
    },
};
use proglog::{CountFormatterKind, ProgLogBuilder};
use std::{
    env, io, io::Write, num::NonZeroUsize, path::PathBuf, sync::Arc, thread::JoinHandle,
    time::Duration,
};
use stitch::{
    align::{
        io::{
            FastqGroupingIterator, FastqThreadReader, OutputMessage, OutputResult,
            READER_CHANNEL_NUM_CHUNKS,
        },
        AlignmentMode, Builder, PrimaryPickingStrategy,
    },
    util::{
        target_seq,
        version::{built_info, built_info::VERSION},
    },
};

////////////////////////////////////////////////////////////////////////////////
// Align (main class) and it's impls
////////////////////////////////////////////////////////////////////////////////

impl ValueEnum for AlignmentMode {
    fn variants<'a>() -> &'a [Self] {
        &[
            Self::Local,
            Self::QueryLocal,
            Self::TargetLocal,
            Self::Global,
        ]
    }
}

impl ValueEnum for PrimaryPickingStrategy {
    fn variants<'a>() -> &'a [Self] {
        &[Self::QueryLength, Self::Score]
    }
}

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
    reads_fastq: PathBuf,

    /// The path to the referece vector/plasmid/construct FASTA sequence.
    #[clap(long, short = 'r', display_order = 2)]
    ref_fasta: PathBuf,

    /// Align to both strands of the reference simulataneously.
    #[clap(long, short = 'd', default_value = "false", display_order = 3)]
    double_strand: bool,

    /// The number of threads to use.
    #[clap(long, short = 't', default_value = "2", display_order = 4)]
    threads: usize,

    /// Assume an unrecognized input (based on file extension) is GZIP compressed.
    #[clap(long, short = 'z', default_value = "false", display_order = 5)]
    decompress: bool,

    /// Pre-align with banded local alignment.
    #[clap(long, short = 'p', default_value = "false", display_order = 6)]
    pre_align: bool,

    /// K-mer size for banded pre-alignment.
    #[clap(long, short = 'k', default_value = "12", display_order = 7)]
    k: usize,

    /// Band size for banded pre-alignment.
    #[clap(long, short = 'w', default_value = "50", display_order = 8)]
    w: usize,

    /// The minimum score required for pre-alignment.
    #[clap(
        long,
        short = 's',
        default_value = "100",
        allow_hyphen_values = true,
        display_order = 9
    )]
    pre_align_min_score: i32,

    /// Only align to contigs that had an alignment score greater than the `--pre-align-min-score`
    #[clap(long, short = 'x', default_value = "true", display_order = 8)]
    pre_align_subset_contigs: bool,

    /// Use soft-clipping for all alignments, otherwise secondary alignemnts will use hard-clipping
    #[clap(long, short = 'S', default_value = "false", display_order = 10)]
    soft_clip: bool,

    /// Use =/X CIGAR operators, otherwise use M
    #[clap(long, short = 'X', default_value = "false", display_order = 11)]
    use_eq_and_x: bool,

    /// Score for a sequence match (must be positive)
    #[clap(long, short = 'A', default_value = "1", display_order = 12)]
    match_score: i32,

    /// Score for a sequence mismatch (must be negative)
    #[clap(
        long,
        short = 'B',
        default_value = "-4",
        allow_hyphen_values = true,
        display_order = 13
    )]
    mismatch_score: i32,

    /// Score for a gap open (must be negative)
    #[clap(
        long,
        short = 'O',
        default_value = "-6",
        allow_hyphen_values = true,
        display_order = 14
    )]
    gap_open: i32,

    /// Score for a gap extend (must be negative); a gap of size k costs '{-O} + {-E}*k'
    #[clap(
        long,
        short = 'E',
        default_value = "-2",
        allow_hyphen_values = true,
        display_order = 15
    )]
    gap_extend: i32,

    /// Score for a target jump (must be negative).
    #[clap(
        long,
        short = 'J',
        default_value = "-10",
        allow_hyphen_values = true,
        display_order = 16
    )]
    jump_score: i32,

    /// Score for a target jump within the same contig and strand (must be negative)
    #[clap(long, allow_hyphen_values = true, display_order = 16)]
    jump_score_same_contig_and_strand: Option<i32>,

    /// Score for a target jump within the same contig and oppoosite strand (must be negative)
    #[clap(long, allow_hyphen_values = true, display_order = 16)]
    jump_score_same_contig_opposite_strand: Option<i32>,

    /// Score for a target jump across different contigs (must be negative)
    #[clap(long, allow_hyphen_values = true, display_order = 16)]
    jump_score_inter_contig: Option<i32>,

    /// The alignment mode:
    /// - Local: aligns a sub-sequence of the read versus a sub-sequence of the reference.
    /// - QueryLocal: aligns a sub-sequence of the read versus the full reference.
    /// - TargetLocal: aligns the full read versus a sub-sequence of the reference.
    /// - Global: aligns the full read versus the full reference.
    #[clap(
        long,
        short = 'm',
        value_parser = PossibleValuesParser::new(AlignmentMode::possible_values())
            .map(|s| s.parse::<AlignmentMode>().unwrap()),
        default_value_t = AlignmentMode::Local,
        ignore_case = true,
        display_order = 17,
        verbatim_doc_comment
    )]
    mode: AlignmentMode,

    /// Determines how to pick the primary alignment:
    /// - QueryLength: the longest aligned query length, then by `Score`
    /// - Score: the best scoring alignment, then by `QueryLength`
    #[clap(
        long,
        short = 'P',
        value_parser = PossibleValuesParser::new(PrimaryPickingStrategy::possible_values())
            .map(|s| s.parse::<PrimaryPickingStrategy>().unwrap()),
        default_value_t = PrimaryPickingStrategy::QueryLength,
        ignore_case = true,
        display_order = 18,
        verbatim_doc_comment
    )]
    pick_primary: PrimaryPickingStrategy,

    /// True to treat the input target as circular.  This allows the alignment to jump back to the
    /// start of the target/reference at no cost.
    #[clap(long, short = 'C', default_value = "false", display_order = 19)]
    circular: bool,

    /// When using `--circular`, re-align reads that may cross the target origin.  This may occur
    /// when a prefix (or suffix) of the read is unaligned, and the alignment starts within the
    /// slop of the target start (or target end).  In this case, the prefix is appended to the
    /// aligned suffix and the new read is re-aligned.
    #[clap(long, default_value = "20", display_order = 20)]
    circular_slop: usize,

    /// Filter out secondary alignments with score X% worse than the primary alignment.
    #[clap(long, default_value = "false", display_order = 21)]
    filter_secondary: bool,

    /// Filter out secondary alignments with score X% worse than the primary alignment.
    #[clap(long, default_value = "10", display_order = 22)]
    filter_secondary_pct: f32,

    /// Generate sub-optimal alignments.
    #[clap(long, default_value = "false", display_order = 23)]
    suboptimal: bool,

    /// Generate sub-optimal alignments with score X% worse than the optimal alignment.
    #[clap(long, default_value = "20", display_order = 24)]
    suboptimal_pct: f32,

    /// The compression level of the output BAM
    #[clap(long, short = 'c', default_value = "0", display_order = 25)]
    compression: u8,
}

impl Align {
    /// Executes the align command
    pub fn execute(&self) -> anyhow::Result<()> {
        info!("Starting alignment...");
        info!("Reading reference FASTA from {}", self.ref_fasta.display());
        info!("Reading reads FASTQ from {}", self.reads_fastq.display());
        // ensure!(self. > 1, "Must specify at least two threads");
        let progress_logger = ProgLogBuilder::new()
            .name("stitch-progress")
            .noun("reads")
            .verb("Processed")
            .unit(
                (READER_CHANNEL_NUM_CHUNKS * self.threads)
                    .try_into()
                    .unwrap(),
            )
            .count_formatter(CountFormatterKind::Comma)
            .build();

        // Create the Builder - we will use this when initializing each thread to create the
        // Aligners object, and to create the SamRecordFormatter for converting alignments to
        // SAM format for writing.
        let mut builder = Builder::default();
        builder
            .mode(self.mode)
            .match_score(self.match_score)
            .mismatch_score(self.mismatch_score)
            .gap_open(self.gap_open)
            .gap_extend(self.gap_extend)
            .default_jump_score(self.jump_score)
            .jump_score_same_contig_and_strand(self.jump_score_same_contig_and_strand)
            .jump_score_same_contig_opposite_strand(self.jump_score_same_contig_opposite_strand)
            .jump_score_inter_contig(self.jump_score_inter_contig)
            .kmer_size(self.k)
            .band_width(self.w)
            .double_strand(self.double_strand)
            .circular(self.circular)
            .circular_slop(self.circular_slop)
            .pre_align(self.pre_align)
            .pre_align_min_score(self.pre_align_min_score)
            .pre_align_subset_contigs(self.pre_align_subset_contigs)
            .suboptimal(self.suboptimal)
            .suboptimal_pct(self.suboptimal_pct)
            .soft_clip(self.soft_clip)
            .use_eq_and_x(self.use_eq_and_x)
            .pick_primary(self.pick_primary)
            .filter_secondary(self.filter_secondary)
            .filter_secondary_pct(self.filter_secondary_pct);
        // share the builder across threads
        let builder = Arc::new(builder);

        // Read in the refefence/target FASTA records - these are shared across threads
        let target_seqs = Arc::new(target_seq::from_fasta(&self.ref_fasta, self.circular)?);

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
                let target_seqs = Arc::clone(&target_seqs);
                let builder = Arc::clone(&builder);
                let k = self.k;

                std::thread::spawn(move || {
                    // Build target hashes to use during alignment
                    // TODO: there should be a way to do this once and share it across threads
                    let target_hashes = target_seqs
                        .iter()
                        .map(|target_seq| target_seq.build_target_hash(k))
                        .collect::<Vec<_>>();
                    let mut aligners = builder.build_aligners(&target_seqs);
                    loop {
                        // Try to process one chunk of alignments
                        if let Ok(msg) = to_align_rx.try_recv() {
                            let iter = FastqGroupingIterator::new(msg.records.into_iter());
                            let mut results: Vec<OutputResult> = Vec::new();
                            for group in iter {
                                let first = group.first().unwrap();
                                let (alignments, maybe_score) =
                                    aligners.align(first, &target_seqs, &target_hashes);

                                for record in group {
                                    let alignments = alignments.clone();
                                    results.push((record, alignments, maybe_score));
                                }
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
                "stitch",
                Map::<Program>::builder()
                    .set_name("stitch")
                    .set_version(VERSION.clone())
                    .set_command_line(command_line)
                    .build()?,
            );
            for target_seq in target_seqs.iter() {
                builder = builder.add_reference_sequence(
                    target_seq.name.parse()?,
                    Map::<ReferenceSequence>::new(NonZeroUsize::try_from(target_seq.len())?),
                );
            }
            builder.build()
        };
        writer.write_header(&header)?;

        // Convert the alignments to SAM records
        let record_formatter = builder.build_sam_record_formatter(&target_seqs);
        loop {
            // Get a receiver for the alignment of a record
            if let Ok(receiver) = reader.to_output_rx.try_recv() {
                let msg = receiver.recv()?;
                for (fastq, alignments, alt_score) in msg.results {
                    progress_logger.record();
                    let records = record_formatter.format(&fastq, &alignments, alt_score)?;
                    for record in records {
                        writer.write_record(&header, &record)?;
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

impl Command for Align {
    fn execute(&self) -> anyhow::Result<()> {
        Align::execute(self)
    }
}

#[cfg(test)]
mod tests {
    use clap::Parser;

    use super::Align;

    /// Check that the argument parser works
    #[test]
    fn test_parse() {
        Align::parse_from(["align", "-f", ".", "-r", "."]);
    }
}
