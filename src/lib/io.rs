use crate::alignment::PairwiseAlignment;
use crate::{is_fastq_path, is_gzip_path};
use anyhow::{Context, Result};
use flate2::bufread::MultiGzDecoder;
use flume::{bounded, Receiver, Sender};
use itertools::{self, Itertools};
use seq_io::fastq::{OwnedRecord as FastqOwnedRecord, Reader as FastqReader};
use std::path::PathBuf;
use std::thread::JoinHandle;
use std::{
    fs::File,
    io::{BufReader, Read},
};

/// 128 KB default buffer size, same as pigz.
pub const GZ_BUFSIZE: usize = 64 * (1 << 10) * 2;
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

pub type OutputResult = (
    FastqOwnedRecord,
    Option<(PairwiseAlignment, bool)>,
    Option<i32>,
);

pub struct OutputMessage {
    pub results: Vec<OutputResult>,
}

pub struct FastqThreadReader {
    /// The [`JoinHandle`] for the thread that is reading.
    pub handle: JoinHandle<Result<()>>,
    /// The channel that will be receiving [`AlignmentsMessage`]s.
    pub to_align_rx: Receiver<InputMessage>,
    /// The channel that will be receiving oneshot recivers of chunks of alignments
    pub to_output_rx: Receiver<Receiver<OutputMessage>>,
}

impl FastqThreadReader {
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
            // Open a FASTQ reader, get an iterator over the records, and chunk them
            let fastq_reader = FastqReader::with_capacity(maybe_decoder_handle, GZ_BUFSIZE)
                .into_records()
                .chunks(RECORDS_PER_CHUNK_PER_THREAD * threads);
            for chunk in &fastq_reader {
                let records: Vec<FastqOwnedRecord> =
                    chunk.map(|r| r.expect("Error reading")).collect();
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

            Ok(())
        });
        Self {
            handle,
            to_align_rx,
            to_output_rx,
        }
    }
}
