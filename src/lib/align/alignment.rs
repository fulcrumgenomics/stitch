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

    /// If the aligmnent starts on the forward strand (used for double strand alignment)
    pub is_forward: bool,

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
        let mut x_index: i32 = self.xstart as i32;
        let mut last_op = self.operations.first().unwrap();
        let mut last_len = 0;
        for op in &self.operations {
            // Should we add the previous operation?
            if (op.is_special() || op != last_op) && last_len > 0 {
                cigar.push_str(&format!(
                    "{}{}",
                    last_len,
                    last_op.as_string(x_index as usize)
                ));
            }
            // Update
            if op.is_special() {
                cigar.push_str(op.as_string(x_index as usize).as_str());
                x_index += op.length_on_x(x_index);
                last_op = op;
                last_len = 0;
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
                last_op.as_string(x_index as usize)
            ));
        }
        cigar
    }

    /// Generates a padded text representation of the alignment for visualization. The returned
    /// sequence will consist of three lines as follows (minus the labels on the left):
    ///
    /// query : ACGTGAACTGACT-ACTGTATGCG
    /// align : |||||  |||||| ||||||||.|
    /// target: ACGTG--CTGACTGACTGTATGGG
    #[allow(dead_code)]
    pub fn padded_string(&self, query: &[u8], target: &[u8]) -> (String, String, String) {
        let mut query_buf: String = String::new();
        let mut align_buf: String = String::new();
        let mut target_buf: String = String::new();
        // TODO: need to revcomp query/target if Xflip is contained...

        let mut query_offset = self.xstart;
        let mut target_offset = self.ystart;

        self.operations.iter().for_each(|op| match *op {
            AlignmentOperation::Match | AlignmentOperation::Subst => {
                query_buf.push(query[query_offset] as char);
                if *op == AlignmentOperation::Match {
                    align_buf.push('|');
                } else {
                    align_buf.push('.');
                }
                target_buf.push(target[target_offset] as char);
                query_offset += 1;
                target_offset += 1;
            }
            AlignmentOperation::Ins => {
                query_buf.push(query[query_offset] as char);
                align_buf.push(' ');
                target_buf.push('-');
                query_offset += 1;
            }
            AlignmentOperation::Del => {
                query_buf.push('-');
                align_buf.push(' ');
                target_buf.push(target[target_offset] as char);
                target_offset += 1;
            }
            AlignmentOperation::Xflip(from_i) => {
                if query_buf.chars().nth(query_buf.len() - 1).unwrap() != ' '
                    && target_buf.chars().nth(target_buf.len() - 1).unwrap() != ' '
                {
                    query_buf.push('X');
                    align_buf.push('X');
                    target_buf.push('X');
                }
                query_offset = from_i;
            }
            AlignmentOperation::Xskip(from_i) => {
                if query_buf.chars().nth(query_buf.len() - 1).unwrap() != ' '
                    && target_buf.chars().nth(target_buf.len() - 1).unwrap() != ' '
                {
                    query_buf.push(' ');
                    align_buf.push(' ');
                    target_buf.push(' ');
                }
                query_offset = from_i;
            }
            AlignmentOperation::Xclip(from_i) => {
                query_offset = from_i;
            }
            _ => panic!("Unknown operator: {op:?}"),
        });
        (query_buf, align_buf, target_buf)
    }
}

impl fmt::Display for Alignment {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "x-span: {}-{} y-span: {}-{} score: {} strand: {} cigar: {} aln-len: {}",
            self.xstart,
            self.xend,
            self.ystart,
            self.yend,
            self.score,
            if self.is_forward { "+" } else { "-" },
            self.cigar(),
            self.length
        )
    }
}
