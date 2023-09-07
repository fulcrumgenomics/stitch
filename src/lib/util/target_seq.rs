use std::path::PathBuf;

use crate::{align::io::BUFFER_SIZE, util::dna::reverse_complement};
use anyhow::Context;
use anyhow::{ensure, Result};
use bio::alignment::sparse::{hash_kmers, HashMapFx};
use fgoxide::io::Io;
use itertools::{self, Itertools};
use seq_io::fasta::Reader as FastaReader;
use seq_io::fasta::Record as FastaRecord;
use std::collections::HashMap;
use std::io::BufRead;

/// Contains the forward and reverse complement of a DNA sequence for a single contig.
#[derive(Default, Debug, PartialEq, Eq, Clone)]
pub struct TargetSeq {
    pub name: String,
    pub fwd: Vec<u8>,
    pub revcomp: Vec<u8>,
    pub circular: bool,
}

/// Contains the k-mer hash for the forward and reverse complement of a DNA sequence for a single
/// contig.
pub struct TargetHash<'a> {
    pub name: String,
    pub fwd_hash: HashMapFx<&'a [u8], Vec<u32>>,
    pub revcomp_hash: HashMapFx<&'a [u8], Vec<u32>>,
}

impl TargetSeq {
    /// Creates a new `TargetSeq` for the contig with the given name and DNA sequence.
    pub fn new(name: &str, seq: &Vec<u8>, circular: bool) -> Self {
        Self {
            name: name.to_string(),
            fwd: seq.clone(),
            revcomp: reverse_complement(seq),
            circular,
        }
    }

    pub fn len(&self) -> usize {
        self.fwd.len()
    }

    /// Creates a new `TargetHash` with the given k-mer size.
    pub fn build_target_hash(&self, k: usize) -> TargetHash {
        TargetHash {
            name: self.name.clone(),
            fwd_hash: hash_kmers(&self.fwd, k),
            revcomp_hash: hash_kmers(&self.revcomp, k),
        }
    }
}

/// Converts the FASTQ header (which may contain whitespaces) to a QNAME for the SAM format.
fn header_to_name(header: &[u8]) -> Result<String> {
    let header: std::borrow::Cow<str> = String::from_utf8_lossy(header);
    header
        .split_whitespace()
        .next()
        .map(std::string::ToString::to_string)
        .context("empty read name")
}

pub fn from_fasta(file: &PathBuf, circular: bool) -> Result<Vec<TargetSeq>> {
    let fg_io: Io = Io::new(5, BUFFER_SIZE);

    // Check if the .dict file exist, and if so, check if which contigs have a circular topology.
    let dict = file.as_path().with_extension(".dict");
    let circular_contigs = if dict.exists() {
        let mut circular_contigs = HashMap::new();
        for line in fg_io.read_lines(&dict)? {
            if !line.starts_with("@SQ") {
                continue;
            }
            let fields = line.split_ascii_whitespace().collect_vec();
            let contig_is_circular = fields
                .iter()
                .filter(|field| field.starts_with("TP"))
                .filter_map(|field| field.split_terminator(':').last())
                .any(|field| field == "circular");
            let contig_name = fields
                .iter()
                .filter(|field| field.starts_with("SN"))
                .find_map(|field| field.split_terminator(':').last())
                .unwrap();
            circular_contigs.insert(contig_name.to_owned(), contig_is_circular);
        }
        circular_contigs
    } else {
        HashMap::new()
    };

    let source: FastaReader<Box<dyn BufRead + Send>> =
        FastaReader::with_capacity(fg_io.new_reader(file)?, BUFFER_SIZE);

    let sequences = source
        .into_records()
        .map(|r| r.unwrap_or_else(|_| panic!("Error reading FASTA")))
        .collect_vec();

    ensure!(!sequences.is_empty(), "Found no sequences in the FASTA");

    sequences
        .iter()
        .map(|record| {
            let sequence = record
                .seq()
                .iter()
                .map(u8::to_ascii_uppercase)
                .collect_vec();
            let name = header_to_name(record.head())?;
            let contig_is_circular = circular_contigs
                .get(&name)
                .map_or(circular, |circular| *circular);
            Ok(TargetSeq::new(&name, &sequence, contig_is_circular))
        })
        .collect()
}
