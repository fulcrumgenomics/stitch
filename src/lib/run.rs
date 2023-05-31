use crate::alignment::aligners::align_double_strand;
use crate::alignment::aligners::align_single_strand;
use crate::alignment::aligners::to_records;
use crate::alignment::aligners::Aligners;
use crate::alignment::scoring::Scoring;
use crate::io::FastqThreadReader;
use crate::io::OutputMessage;
use crate::io::OutputResult;
use crate::io::BUFFER_SIZE;
use crate::io::READER_CHANNEL_NUM_CHUNKS;
use crate::opts::Opts;
use crate::target_seq::TargetSeq;
use crate::util::built_info::VERSION;
use anyhow::Context;
use anyhow::{ensure, Result};
use bio::alignment::pairwise::MatchParams;
use fgoxide::io::Io;
use flume::unbounded;
use itertools::{self, Itertools};
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
use std::env;
use std::io;
use std::io::BufRead;
use std::io::Write;
use std::num::NonZeroUsize;
use std::path::PathBuf;
use std::thread::JoinHandle;
use std::time::Duration;

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

pub fn run(opts: &Opts) -> Result<()> {
    // ensure!(opts.threads > 1, "Must specify at least two threads");
    let progress_logger = ProgLogBuilder::new()
        .name("fqcv-progress")
        .noun("reads")
        .verb("Aligned")
        .unit(
            (READER_CHANNEL_NUM_CHUNKS * opts.threads)
                .try_into()
                .unwrap(),
        )
        .count_formatter(CountFormatterKind::Comma)
        .build();

    // Read in the refefence/target FASTA
    let (target_seq, target_name) = read_target(&opts.ref_fasta)?;

    // Create the thread to read in the FASTQ records
    let reader = FastqThreadReader::new(opts.reads_fastq.clone(), opts.decompress, opts.threads);

    // Create the channel to gracefully signal a shutdown of the aligner threads
    let (shutdown_tx, shutdown_rx) = unbounded::<()>();

    let sleep_delay = Duration::from_millis(25);

    let thread_handles: Vec<JoinHandle<Result<()>>> = (0..opts.threads)
        .map(|_| {
            let to_align_rx = reader.to_align_rx.clone();
            let shutdown_rx = shutdown_rx.clone();
            let target_name = target_name.clone();
            let target_seq = target_seq.clone();
            let mut aligners = Aligners::new(opts, target_seq.len());
            let opts = opts.clone();

            std::thread::spawn(move || {
                let target_seq = TargetSeq::new(&target_name, &target_seq);
                let target_hash = target_seq.build_target_hash(opts.k);
                loop {
                    // Try to process one chunk of alignments
                    if let Ok(msg) = to_align_rx.try_recv() {
                        let results: Vec<OutputResult> = msg
                            .records
                            .into_iter()
                            .map(|record| {
                                if opts.double_strand {
                                    align_double_strand(
                                        record,
                                        &target_seq,
                                        &target_hash,
                                        &mut aligners,
                                        opts.pre_align,
                                        opts.pre_align_min_score,
                                    )
                                } else {
                                    align_single_strand(
                                        record,
                                        &target_seq,
                                        &target_hash,
                                        &mut aligners,
                                        opts.pre_align,
                                        opts.pre_align_min_score,
                                    )
                                }
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
        let match_fn: MatchParams = MatchParams::new(opts.match_score, opts.mismatch_score);
        let scoring = Scoring::new(opts.gap_open, opts.gap_extend, opts.jump_score, match_fn);
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
                    &scoring,
                    opts.by_primary,
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
