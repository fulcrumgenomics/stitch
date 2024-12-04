use super::alignment::Alignment;
use crate::util::io::{is_fastq_path, is_gzip_path};
use anyhow::{Context, Result};
use derive_getters::Getters;
use flate2::bufread::MultiGzDecoder;
use flume::{bounded, Receiver, Sender};
use itertools::Itertools;
use seq_io::{
    fasta::{Reader as FastaReader, RefRecord as FastaRefRecord},
    fastq::{Reader as FastqReader, RefRecord as FastqRefRecord},
};
use std::{
    fs::File,
    io::{BufReader, Read},
    iter::Peekable,
    path::PathBuf,
    thread::JoinHandle,
};

/// 128 KB default buffer size, same as pigz.
pub const GZ_BUFSIZE: usize = 64 * (1 << 10) * 2;

/// The number of FASTQ records to includ per chunk, scaled but the # of CPUS.
pub const RECORDS_PER_CHUNK_PER_THREAD: usize = 10;

/// The number of chunks allowed in a channel
pub const READER_CHANNEL_NUM_CHUNKS: usize = 100;

/// Enumeration of supported input file formats
#[derive(Copy, Clone, Debug, Default, PartialEq)]
pub enum Format {
    #[default]
    FASTQ,
    FASTA,
}

#[derive(Clone, Debug, Getters)]
/// Common record struct that supports both FASTA and FASTQ
pub struct FastxOwnedRecord {
    pub head: Vec<u8>,
    pub seq: Vec<u8>,
    pub qual: Option<Vec<u8>>,
}

impl FastxOwnedRecord {
    pub fn from_fastq(record: &FastqRefRecord) -> Self {
        let owned_record = record.to_owned_record();
        Self {
            head: owned_record.head,
            seq: owned_record.seq,
            qual: Some(owned_record.qual),
        }
    }

    pub fn from_fasta(record: &FastaRefRecord) -> Self {
        let owned_record = record.to_owned_record();
        Self {
            head: owned_record.head,
            seq: owned_record.seq,
            qual: None,
        }
    }

    pub fn seq_upper_case(&self) -> Vec<u8> {
        self.seq.iter().map(u8::to_ascii_uppercase).collect_vec()
    }
}

pub struct FastaToFastxIterator(FastaReader<Box<dyn Read>>);

impl Iterator for FastaToFastxIterator {
    type Item = FastxOwnedRecord;

    fn next(&mut self) -> Option<Self::Item> {
        self.0
            .next()
            .map(|record| record.expect("Error reading fasta record"))
            .map(|record| FastxOwnedRecord::from_fasta(&record))
    }
}

pub struct FastqToFastxIterator(FastqReader<Box<dyn Read>>);

impl Iterator for FastqToFastxIterator {
    type Item = FastxOwnedRecord;

    fn next(&mut self) -> Option<Self::Item> {
        self.0
            .next()
            .map(|record| record.expect("Error reading fasta record"))
            .map(|record| FastxOwnedRecord::from_fastq(&record))
    }
}

/// A message that is sent from a [`FastqThreadReader`] to the aligner threadpool to align a chunk
/// of FASTQ records.
#[derive(Debug)]
pub struct InputMessage {
    /// The FASTQ records to align
    pub records: Vec<FastxOwnedRecord>,

    /// Where the records will be sent after alignment
    pub oneshot: Sender<OutputMessage>,
}

/// The output result of a single pairwise alignment as a triple:
/// 1. the FASTQ record that was aligned
/// 2. None if no alignment was found, otherwise some tuple of pairwise alignment and the target
///    strand to which the alignment was made.
/// 3. The alignment score, if aligned.
pub type OutputResult = (FastxOwnedRecord, Vec<Alignment>, Option<i32>);

/// The container for a chunk of pairwise alignments, one per input FASTQ record.
pub struct OutputMessage {
    pub results: Vec<OutputResult>,
}

pub struct FastxGroupingIterator<I: Iterator<Item = FastxOwnedRecord>>(Peekable<I>);

impl<I: Iterator<Item = FastxOwnedRecord>> FastxGroupingIterator<I> {
    pub fn new(iter: I) -> Self {
        Self(iter.peekable())
    }
}

impl<I: Iterator<Item = FastxOwnedRecord>> Iterator for FastxGroupingIterator<I> {
    type Item = Vec<FastxOwnedRecord>;

    #[inline]
    fn next(&mut self) -> Option<Vec<FastxOwnedRecord>> {
        match self.0.next() {
            None => None,
            Some(record) => {
                let mut items: Vec<FastxOwnedRecord> = vec![record];
                while let Some(record) = self.0.peek() {
                    if record.seq == items[0].seq {
                        items.push(self.0.next().expect("Expected record"));
                    } else {
                        break;
                    }
                }
                Some(items)
            }
        }
    }
}

/// A FASTQ reader that runs in its own thread and chunks reads to send to a pool of aligners.
pub struct FastxThreadReader {
    /// The [`JoinHandle`] for the thread that is reading.
    pub handle: JoinHandle<Result<()>>,
    /// The channel that will be receiving [`AlignmentsMessage`]s.
    pub to_align_rx: Receiver<InputMessage>,
    /// The channel that will be receiving oneshot recivers of chunks of alignments
    pub to_output_rx: Receiver<Receiver<OutputMessage>>,
}

impl FastxThreadReader {
    /// Writes the chunk of records to the alignment channel, as well as a receiver to the output
    /// channel.
    fn write_records_to_txs(
        records: Vec<FastxOwnedRecord>,
        to_align_tx: &Sender<InputMessage>,
        to_output_tx: &Sender<Receiver<OutputMessage>>,
    ) {
        let (records_tx, records_rx) = flume::unbounded(); // oneshot channel
        let input_msg = InputMessage {
            records,
            oneshot: records_tx,
        };
        to_align_tx.send(input_msg).expect("Error sending message");
        to_output_tx
            .send(records_rx)
            .expect("Error sending receiver");
    }

    /// Creates a new `FastqThreadReader` in a new thread.
    pub fn new(file: PathBuf, format: Format, decompress: bool, threads: usize) -> Self {
        // Channel to send chunks of records to align
        let (to_align_tx, to_align_rx): (Sender<InputMessage>, Receiver<InputMessage>) =
            bounded(READER_CHANNEL_NUM_CHUNKS * threads);

        // Channel to send receivers for each aligned chunk of records. The receivers maintain the
        // order of the records.
        let (to_output_tx, to_output_rx): (
            Sender<Receiver<OutputMessage>>,
            Receiver<Receiver<OutputMessage>>,
        ) = bounded(READER_CHANNEL_NUM_CHUNKS * threads);

        let handle = std::thread::spawn(move || {
            // Open the file or standad input
            let raw_handle = if file.as_os_str() == "-" {
                Box::new(std::io::stdin()) as Box<dyn Read>
            } else {
                let handle = File::open(&file)
                    .with_context(|| format!("Error opening input: {}", file.display()))
                    .unwrap();
                Box::new(handle) as Box<dyn Read>
            };
            // Wrap it in a buffer
            let buf_handle = BufReader::with_capacity(GZ_BUFSIZE, raw_handle);
            // Maybe wrap it in a decompressor
            let maybe_decoder_handle = {
                let is_gzip = is_gzip_path(&file) || (!is_fastq_path(&file) && decompress);
                if is_gzip {
                    Box::new(MultiGzDecoder::new(buf_handle)) as Box<dyn Read>
                } else {
                    Box::new(buf_handle) as Box<dyn Read>
                }
            };

            // Open an input file reader, group reads that have the same read sequence, then chunk
            // the reads to send over the output channel, keeping reads with the same read sequence
            // grouped together.
            let fastq_iter: Box<dyn Iterator<Item = FastxOwnedRecord>> = match format {
                Format::FASTQ => Box::new(FastqToFastxIterator(FastqReader::with_capacity(
                    maybe_decoder_handle,
                    GZ_BUFSIZE,
                ))),
                Format::FASTA => Box::new(FastaToFastxIterator(FastaReader::with_capacity(
                    maybe_decoder_handle,
                    GZ_BUFSIZE,
                ))),
            };
            let fastq_grouping_iter = FastxGroupingIterator::new(fastq_iter);
            let mut records = Vec::with_capacity(RECORDS_PER_CHUNK_PER_THREAD);
            for chunk in fastq_grouping_iter {
                records.extend(chunk);
                if records.len() >= RECORDS_PER_CHUNK_PER_THREAD {
                    Self::write_records_to_txs(records, &to_align_tx, &to_output_tx);
                    records = Vec::with_capacity(RECORDS_PER_CHUNK_PER_THREAD);
                }
            }
            if !records.is_empty() {
                Self::write_records_to_txs(records, &to_align_tx, &to_output_tx);
            }

            Ok(())
        });
        Self {
            handle,
            to_align_rx,
            to_output_rx,
        }
    }
}
