use anyhow::{Context, Result};
use flate2::bufread::MultiGzDecoder;
use flume::{bounded, Receiver, Sender};
use seq_io::fastq::{OwnedRecord as FastqOwnedRecord, Reader as FastqReader};
use std::path::PathBuf;
use std::thread::JoinHandle;
use std::{
    fs::File,
    io::{BufReader, Read},
};

use crate::util::io::{is_fastq_path, is_gzip_path};

use super::alignment::Alignment;

/// 128 KB default buffer size, same as pigz.
pub const GZ_BUFSIZE: usize = 64 * (1 << 10) * 2;

/// The buffer size for reading FASTAs
pub const BUFFER_SIZE: usize = 1024 * 1024;

/// The number of FASTQ records to includ per chunk, scaled but the # of CPUS.
pub const RECORDS_PER_CHUNK_PER_THREAD: usize = 10;

/// The number of chunks allowed in a channel
pub const READER_CHANNEL_NUM_CHUNKS: usize = 100;

/// A message that is sent from a [`FastqThreadReader`] to the aligner threadpool to align a chunk
/// of FASTQ records.
#[derive(Debug)]
pub struct InputMessage {
    /// The FASTQ records to align
    pub records: Vec<FastqOwnedRecord>,

    /// Where the records will be sent after alignment
    pub oneshot: Sender<OutputMessage>,
}

/// The output result of a single pairwise alignment as a triple:
/// 1. the FASTQ record that was aligned
/// 2. None if no alignment was found, otherwise some tuple of pairwise alignment and the target
///   strand to which the alignment was made.
/// 3. The alignment score, if aligned.
pub type OutputResult = (FastqOwnedRecord, Option<Alignment>, Option<i32>);

/// The container for a chunk of pairwise alignments, one per input FASTQ record.
pub struct OutputMessage {
    pub results: Vec<OutputResult>,
}

pub struct FastqGroupingIterator<I: Iterator<Item = FastqOwnedRecord>> {
    pub record: Option<FastqOwnedRecord>,
    pub iter: I,
}

impl<I: Iterator<Item = FastqOwnedRecord>> FastqGroupingIterator<I> {
    pub fn new(iter: I) -> Self {
        Self { record: None, iter }
    }
}

impl<I: Iterator<Item = FastqOwnedRecord>> Iterator for FastqGroupingIterator<I> {
    type Item = Vec<FastqOwnedRecord>;

    #[inline]
    fn next(&mut self) -> Option<Vec<FastqOwnedRecord>> {
        if self.record.is_none() {
            self.record = self.iter.next();
        }
        match self.record {
            None => None,
            Some(ref record) => {
                let mut items: Vec<FastqOwnedRecord> = Vec::new();
                let sequence = record.seq.clone();
                items.push(record.clone());
                self.record = None;
                for record in &mut self.iter {
                    if record.seq == sequence {
                        items.push(record);
                    } else {
                        self.record = Some(record);
                        break;
                    }
                }
                Some(items)
            }
        }
    }
}

/// A FASTQ reader that runs in its own thread and chunks reads to send to a pool of aligners.
pub struct FastqThreadReader {
    /// The [`JoinHandle`] for the thread that is reading.
    pub handle: JoinHandle<Result<()>>,
    /// The channel that will be receiving [`AlignmentsMessage`]s.
    pub to_align_rx: Receiver<InputMessage>,
    /// The channel that will be receiving oneshot recivers of chunks of alignments
    pub to_output_rx: Receiver<Receiver<OutputMessage>>,
}

impl FastqThreadReader {
    /// Writes the chunk of records to the alignment channel, as well as a receiver to the output
    /// channel.
    fn write_records_to_txs(
        records: &[FastqOwnedRecord],
        to_align_tx: &Sender<InputMessage>,
        to_output_tx: &Sender<Receiver<OutputMessage>>,
    ) {
        let (records_tx, records_rx) = flume::unbounded(); // oneshot channel
        let input_msg = InputMessage {
            records: records.to_vec(),
            oneshot: records_tx,
        };
        to_align_tx.send(input_msg).expect("Error sending message");
        to_output_tx
            .send(records_rx)
            .expect("Error sending receiver");
    }

    /// Creates a new `FastqThreadReader` in a new thread.
    pub fn new(file: PathBuf, decompress: bool, threads: usize) -> Self {
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

            // Open a FASTQ reader, then group reads that have the same read sequence, then chunk
            // chunk the reads to send over the output channel, keeping reads with the same read
            // sequence grouped together.
            let fastq_iter = FastqReader::with_capacity(maybe_decoder_handle, GZ_BUFSIZE)
                .into_records()
                .map(|r| r.expect("Error reading"));
            let fastq_grouping_iter = FastqGroupingIterator::new(fastq_iter);
            let mut records = Vec::new();
            for chunk in fastq_grouping_iter {
                records.extend(chunk);
                if records.len() >= RECORDS_PER_CHUNK_PER_THREAD {
                    Self::write_records_to_txs(&records, &to_align_tx, &to_output_tx);
                    records.clear();
                }
            }
            if !records.is_empty() {
                Self::write_records_to_txs(&records, &to_align_tx, &to_output_tx);
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
