use bio::alignment::sparse::{hash_kmers, HashMapFx};

use crate::util::dna::reverse_complement;

/// Contains the forward and reverse complement of a DNA sequence for a single contig.
#[derive(Default, Debug, PartialEq, Eq, Clone)]
pub struct TargetSeq {
    pub name: String,
    pub fwd: Vec<u8>,
    pub revcomp: Vec<u8>,
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
    pub fn new(name: &str, seq: &Vec<u8>) -> Self {
        Self {
            name: name.to_string(),
            fwd: seq.clone(),
            revcomp: reverse_complement(seq),
        }
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
