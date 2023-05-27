use crate::alignment::constants::AlignmentMode;
use crate::alignment::constants::AlignmentOperation;

/// We consider alignment between two sequences x and  y. x is the query or read sequence
/// and y is the reference or template sequence. An alignment, consisting of a score,
/// the start and end position of the alignment on sequence x and sequence y, the
/// lengths of sequences x and y, and the alignment edit operations. The start position
/// and end position of the alignment does not include the clipped regions. The length
/// of clipped regions are already encapsulated in the Alignment Operation.
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(Debug, Eq, PartialEq, Clone, Default)]
pub struct PairwiseAlignment {
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
}

/// Generates a padded text representation of the alignment for visualization. The returned
/// sequence will consist of three lines as follows (minus the labels on the left):
///
/// query : ACGTGAACTGACT-ACTGTATGCG
/// align : |||||  |||||| ||||||||.|
/// target: ACGTG--CTGACTGACTGTATGGG
pub fn padded_string(
    alignment: &PairwiseAlignment,
    query: &[u8],
    target: &[u8],
) -> (String, String, String) {
    let mut query_buf: String = String::new();
    let mut align_buf: String = String::new();
    let mut target_buf: String = String::new();
    // TODO: need to revcomp query/target if Xflip is contained...

    let mut query_offset = alignment.xstart;
    let mut target_offset = alignment.ystart;

    alignment.operations.iter().for_each(|op| match *op {
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
