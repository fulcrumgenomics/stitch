use itertools::Itertools;

use crate::alignment::constants::AlignmentOperation;
use crate::alignment::pairwise::PairwiseAlignment;
use crate::alignment::scoring::Scoring;
use bio::alignment::pairwise::MatchFunc;

#[derive(Debug, Eq, PartialEq, Clone, Default)]
pub struct SubAlignment {
    pub query_start: usize,
    pub query_end: usize,
    pub target_start: usize,
    pub target_end: usize,
    pub is_forward: bool,
    pub cigar: String,
    pub score: i32,
}

pub struct SubAlignmentBuilder {
    use_eq_and_x: bool,
    match_char: char,
    mismatch_char: char,
    del_char: char,
    ins_char: char,
    cigar: String,
    query_start: usize,
    target_start: usize,
    query_offset: usize,
    target_offset: usize,
    score: i32,
    is_forward: bool,
}

impl SubAlignmentBuilder {
    fn cmp_op(&self, last: AlignmentOperation, cur: AlignmentOperation) -> bool {
        if self.use_eq_and_x {
            last == cur
        } else {
            last == cur
                || (last == AlignmentOperation::Subst && cur == AlignmentOperation::Match)
                || (last == AlignmentOperation::Match && cur == AlignmentOperation::Subst)
        }
    }

    fn add_op<F: MatchFunc>(
        &mut self,
        op: AlignmentOperation,
        op_len: usize,
        scoring: &Scoring<F>,
    ) -> Option<SubAlignment> {
        match op {
            AlignmentOperation::Match => {
                self.score += scoring.match_fn.score(b'A', b'A') * (op_len as i32);
                self.query_offset += op_len;
                self.target_offset += op_len;
                self.cigar.push_str(&format!("{op_len}{}", self.match_char));
                None
            }
            AlignmentOperation::Subst => {
                self.score += scoring.match_fn.score(b'A', b'C') * (op_len as i32);
                self.query_offset += op_len;
                self.target_offset += op_len;
                self.cigar
                    .push_str(&format!("{op_len}{}", self.mismatch_char));
                None
            }
            AlignmentOperation::Del => {
                self.score += scoring.gap_open + (scoring.gap_extend * op_len as i32);
                self.target_offset += op_len;
                self.cigar.push_str(&format!("{op_len}{}", self.del_char));
                None
            }
            AlignmentOperation::Ins => {
                self.score += scoring.gap_open + (scoring.gap_extend * op_len as i32);
                self.query_offset += op_len;
                self.cigar.push_str(&format!("{op_len}{}", self.ins_char));
                None
            }
            AlignmentOperation::Xskip(new_query_start) => {
                let alignment = SubAlignment {
                    query_start: self.query_start,
                    query_end: self.query_offset,
                    target_start: self.target_start,
                    target_end: self.target_offset,
                    is_forward: self.is_forward,
                    cigar: self.cigar.clone(),
                    score: self.score,
                };

                // reset
                self.cigar.clear();
                self.target_start = self.target_offset;
                self.query_start = new_query_start;
                self.query_offset = new_query_start;
                self.score = 0;

                Some(alignment)
            }
            AlignmentOperation::Xflip(new_query_start) => {
                let alignment = SubAlignment {
                    query_start: self.query_start,
                    query_end: self.query_offset,
                    target_start: self.target_start,
                    target_end: self.target_offset,
                    is_forward: self.is_forward,
                    cigar: self.cigar.clone(),
                    score: self.score,
                };

                // reset
                self.cigar.clear();
                self.target_start = self.target_offset;
                self.query_start = new_query_start;
                self.query_offset = new_query_start;
                self.score = 0;
                self.is_forward = !self.is_forward;

                Some(alignment)
            }
            AlignmentOperation::Xclip(_) => None, // Ignore
            _ => panic!("Unsupported operator: {op:?}"),
        }
    }

    pub fn new(use_eq_and_x: bool) -> Self {
        Self {
            use_eq_and_x,
            match_char: if use_eq_and_x { '=' } else { 'M' },
            mismatch_char: if use_eq_and_x { 'X' } else { 'M' },
            del_char: 'D',
            ins_char: 'I',
            cigar: String::new(),
            query_start: 0,
            target_start: 0,
            query_offset: 0,
            target_offset: 0,
            score: 0,
            is_forward: true,
        }
    }

    pub fn swap_cigar(cigar: &str) -> String {
        cigar
            .chars()
            .map(|c| {
                if c == 'I' {
                    'D'
                } else if c == 'D' {
                    'I'
                } else {
                    c
                }
            })
            .collect()
    }

    pub fn build<F: MatchFunc>(
        &mut self,
        alignment: &PairwiseAlignment,
        swap: bool,
        scoring: &Scoring<F>,
    ) -> Vec<SubAlignment> {
        self.cigar.clear();
        self.query_start = alignment.xstart;
        self.target_start = alignment.ystart;
        self.query_offset = self.query_start;
        self.target_offset = self.target_start;
        self.score = 0;
        self.is_forward = alignment.is_forward;

        let mut alignments = Vec::new();
        let mut last = alignment.operations[0];
        let mut op_len = 0;
        for i in 0..alignment.operations.len() {
            let op = alignment.operations[i];
            if self.cmp_op(last, op) {
                op_len += 1;
            } else {
                if let Some(alignment) = self.add_op(last, op_len, scoring) {
                    alignments.push(alignment);
                }
                op_len = 1;
            }
            last = op;
        }
        if let Some(alignment) = self.add_op(last, op_len, scoring) {
            alignments.push(alignment);
        }

        let alignment = SubAlignment {
            query_start: self.query_start,
            query_end: self.query_offset,
            target_start: self.target_start,
            target_end: self.target_offset,
            is_forward: self.is_forward,
            cigar: self.cigar.clone(),
            score: self.score,
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
                    is_forward: a.is_forward,
                    cigar: if swap {
                        Self::swap_cigar(&a.cigar)
                    } else {
                        a.cigar.clone()
                    },
                    score: a.score,
                })
                .collect_vec()
        } else {
            alignments
        }
    }
}
