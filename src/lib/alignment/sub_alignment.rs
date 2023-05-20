use itertools::Itertools;

use crate::alignment::constants::AlignmentOperation;
use crate::alignment::pairwise::PairwiseAlignment;
use crate::alignment::scoring::Scoring;

use super::scoring::MatchFunc;

#[derive(Debug, Eq, PartialEq, Clone, Default)]
pub struct SubAlignment {
    pub query_start: usize,
    pub query_end: usize,
    pub target_start: usize,
    pub target_end: usize,
    pub cigar: String,
    pub score: i32,
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
