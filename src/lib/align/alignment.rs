use std::fmt;

use super::aligners::constants::{AlignmentMode, AlignmentOperation};

/// We consider alignment between two sequences x and  y. x is the query or read sequence
/// and y is the reference or template sequence. An alignment, consisting of a score,
/// the start and end position of the alignment on sequence x and sequence y, the
/// lengths of sequences x and y, and the alignment edit operations. The start position
/// and end position of the alignment does not include the clipped regions. The length
/// of clipped regions are already encapsulated in the Alignment Operation.
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(Debug, Eq, PartialEq, Clone, Default)]
pub struct Alignment {
    // FIXME: rename to Alignment
    /// Smith-Waterman alignment score
    pub score: i32,

    /// Start position of alignment in reference (0-based)
    pub ystart: usize,

    /// Start position of alignment in query (0-based)
    pub xstart: usize,

    /// End position of alignment in reference (0-based exclusive)
    pub yend: usize,

    /// End position of alignment in query (0-based exclusive)
    pub xend: usize,

    /// Length of the reference sequence (not aligned ylen, the original length of y!)
    pub ylen: usize,

    /// Length of the query sequence (not aligned xlen, the original length of x!)
    pub xlen: usize,

    /// The contig index (0-based), corresponding to the aligner for the aligned contig.
    pub contig_idx: usize,

    /// Vector of alignment operations
    pub operations: Vec<AlignmentOperation>,
    pub mode: AlignmentMode,

    /// Alignment length, excluding any clipped (x or y) bases and jumps.
    pub length: usize,
}

#[allow(dead_code)]
impl Alignment {
    pub fn cigar(&self) -> String {
        let mut cigar = String::new();
        if self.operations.is_empty() {
            return cigar;
        }
        let mut contig_idx = self.contig_idx;
        let mut x_index: i32 = self.xstart as i32;
        let mut last_op = self.operations.first().unwrap();
        let mut last_len = 0;
        for op in &self.operations {
            // Should we add the previous operation?
            if (op.is_special() || op != last_op) && last_len > 0 {
                cigar.push_str(&format!(
                    "{}{}",
                    last_len,
                    last_op.as_string(contig_idx, x_index as usize)
                ));
            }
            // Update
            if op.is_special() {
                cigar.push_str(op.as_string(contig_idx, x_index as usize).as_str());
                x_index += op.length_on_x(x_index);
                last_op = op;
                last_len = 0;
                if let AlignmentOperation::Xjump(new_contig_index, _) = op {
                    contig_idx = *new_contig_index;
                }
            } else if op == last_op {
                x_index += op.length_on_x(x_index);
                last_len += 1;
            } else {
                x_index += op.length_on_x(x_index);
                last_op = op;
                last_len = 1;
            }
        }
        if last_len > 0 {
            cigar.push_str(&format!(
                "{}{}",
                last_len,
                last_op.as_string(contig_idx, x_index as usize)
            ));
        }
        cigar
    }
}

impl fmt::Display for Alignment {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "contig-idx: {} x-span: {}-{} y-span: {}-{} score: {} cigar: {} aln-len: {}",
            self.contig_idx,
            self.xstart,
            self.xend,
            self.ystart,
            self.yend,
            self.score,
            self.cigar(),
            self.length
        )
    }
}
