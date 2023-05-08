use itertools::Itertools;

use crate::pairwise::{MatchFunc, Scoring};
use clap::ValueEnum;

/// Alignment operations supported are match, substitution, insertion, deletion
/// and clipping. Clipping is a special boundary condition where you are allowed
/// to clip off the beginning/end of the sequence for a fixed clip penalty. The
/// clip penalty could be different for the two sequences x and y, and the
/// clipping operations on both are distinguishable (Xclip and Yclip). The usize
/// value associated with the clipping operations are the lengths clipped. In case
/// of standard modes like Global, Semi-Global and Local alignment, the clip operations
/// are filtered out
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(Eq, PartialEq, Debug, Copy, Clone, Hash)]
pub enum AlignmentOperation {
    Match,
    Subst,
    Del,
    Ins,
    Xclip(usize),
    Yclip(usize),
    Xskip(usize),
}

/// The modes of alignment supported by the aligner include standard modes such as
/// Global, Semi-Global and Local alignment. In addition to this, user can also invoke
/// the custom mode. In the custom mode, users can explicitly specify the clipping penalties
/// for prefix and suffix of strings 'x' and 'y' independently. Under the hood the standard
/// modes are implemented as special cases of the custom mode with the clipping penalties
/// appropriately set.
///
/// The default alignment mode is Global.
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(Default, Debug, PartialEq, Eq, Copy, Clone, ValueEnum)]
pub enum AlignmentMode {
    #[default]
    Local,
    Semiglobal,
    Global,
    Custom,
}

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

    /// Start position of alignment in reference
    pub ystart: usize,

    /// Start position of alignment in query
    pub xstart: usize,

    /// End position of alignment in reference
    pub yend: usize,

    /// End position of alignment in query
    pub xend: usize,

    /// Length of the reference sequence
    pub ylen: usize,

    /// Length of the query sequence
    pub xlen: usize,

    /// Vector of alignment operations
    pub operations: Vec<AlignmentOperation>,
    pub mode: AlignmentMode,
}

#[derive(Debug, Eq, PartialEq, Clone, Default)]
pub struct SubAlignment {
    pub query_start: usize,
    pub query_end: usize,
    pub target_start: usize,
    pub target_end: usize,
    pub cigar: String,
    pub score: i32,
}

pub fn padded_string(
    alignment: &PairwiseAlignment,
    query: &[u8],
    target: &[u8],
) -> (String, String, String) {
    let mut query_buf: String = String::new();
    let mut align_buf: String = String::new();
    let mut target_buf: String = String::new();

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

fn cmp_op(last: AlignmentOperation, cur: AlignmentOperation, use_eq_and_x: bool) -> bool {
    if use_eq_and_x {
        last == cur
    } else {
        last == cur
            || (last == AlignmentOperation::Subst && cur == AlignmentOperation::Match)
            || (last == AlignmentOperation::Match && cur == AlignmentOperation::Subst)
    }
}

impl SubAlignment {
    // TODO: pass the query and target sequences to use a custom scoring function
    pub fn build(
        alignment: &PairwiseAlignment,
        scoring: Scoring<impl Fn(u8, u8) -> i32>,
        use_eq_and_x: bool,
        swap: bool,
    ) -> Vec<SubAlignment> {
        let match_str: &str = if use_eq_and_x { "=" } else { "M" };
        let mismatch_str: &str = if use_eq_and_x { "X" } else { "M" };
        let del_str: &str = if swap { "I" } else { "D" };
        let ins_str: &str = if swap { "D" } else { "I" };
        let mut cigar: String = String::new();
        let mut query_start = alignment.xstart;
        let mut target_start = alignment.ystart;
        let mut query_offset = query_start;
        let mut target_offset = target_start;
        let mut op_len: usize = 0;
        let mut alignments: Vec<SubAlignment> = Vec::new();
        let mut score: i32 = 0;

        let mut add_op = |op: AlignmentOperation, op_len| {
            match op {
                AlignmentOperation::Match => {
                    score += scoring.match_fn.score(b'A', b'A') * (op_len as i32);
                    query_offset += op_len;
                    target_offset += op_len;
                    cigar.push_str(&format!("{op_len}{match_str}"));
                }
                AlignmentOperation::Subst => {
                    score += scoring.match_fn.score(b'A', b'C') * (op_len as i32);
                    query_offset += op_len;
                    target_offset += op_len;
                    cigar.push_str(&format!("{op_len}{mismatch_str}"));
                }
                AlignmentOperation::Del => {
                    score += scoring.gap_open + (scoring.gap_extend * op_len as i32);
                    target_offset += op_len;
                    cigar.push_str(&format!("{op_len}{del_str}"));
                }
                AlignmentOperation::Ins => {
                    score += scoring.gap_open + (scoring.gap_extend * op_len as i32);
                    query_offset += op_len;
                    cigar.push_str(&format!("{op_len}{ins_str}"));
                }
                AlignmentOperation::Xskip(new_query_start) => {
                    let alignment = SubAlignment {
                        query_start,
                        query_end: query_offset,
                        target_start,
                        target_end: target_offset,
                        cigar: cigar.clone(),
                        score,
                    };
                    alignments.push(alignment);

                    // reset
                    cigar.clear();
                    target_start = target_offset;
                    query_start = new_query_start;
                    query_offset = new_query_start;
                    score = 0;
                }
                AlignmentOperation::Xclip(_) => (), // Ignore
                _ => panic!("Unknown operator: {op:?}"),
            };
        };

        let mut last = alignment.operations[0];
        for i in 0..alignment.operations.len() {
            let op = alignment.operations[i];
            if cmp_op(last, op, use_eq_and_x) {
                op_len += 1;
            } else {
                add_op(last, op_len);
                op_len = 1;
            }
            last = op;
        }
        add_op(last, op_len);

        let alignment = SubAlignment {
            query_start,
            query_end: query_offset,
            target_start,
            target_end: target_offset,
            cigar: cigar.clone(),
            score,
        };
        alignments.push(alignment);

        if swap {
            alignments
                .iter()
                .map(|a| SubAlignment {
                    query_start: a.target_start,
                    query_end: a.target_end,
                    target_start: a.query_start,
                    target_end: a.query_end,
                    cigar: a.cigar.clone(),
                    score: a.score,
                })
                .collect_vec()
        } else {
            alignments
        }
    }
}
