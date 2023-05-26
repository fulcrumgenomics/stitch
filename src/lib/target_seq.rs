use bio::alignment::sparse::{hash_kmers, HashMapFx};

use crate::util::reverse_complement;

pub struct TargetHash<'a> {
    pub fwd_hash: HashMapFx<&'a [u8], Vec<u32>>,
    pub revcomp_hash: HashMapFx<&'a [u8], Vec<u32>>,
}

#[derive(Default, Debug, PartialEq, Eq, Clone)]
pub struct TargetSeq {
    pub name: String,
    pub fwd: Vec<u8>,
    pub revcomp: Vec<u8>,
}

impl TargetSeq {
    pub fn new(name: &str, seq: &Vec<u8>) -> Self {
        Self {
            name: name.to_string(),
            fwd: seq.clone(),
            revcomp: reverse_complement(seq),
        }
    }

    pub fn build_target_hash(&self, k: usize) -> TargetHash {
        TargetHash {
            fwd_hash: hash_kmers(&self.fwd, k),
            revcomp_hash: hash_kmers(&self.revcomp, k),
        }
    }
}
