use itertools::Itertools;
use noodles::sam::record::cigar::Op;
use noodles::sam::record::Cigar;

use super::aligners::constants::AlignmentOperation;
use super::alignment::Alignment;
use super::scoring::Scoring;
use bio::alignment::pairwise::MatchFunc;
use noodles::sam::record::cigar::op::Kind;

/// A pairwise alignment with no jumps allowed.
#[derive(Debug, Eq, PartialEq, Clone, Default)]
pub struct SubAlignment {
    pub query_start: usize,
    pub query_end: usize,
    pub target_start: usize,
    pub target_end: usize,
    pub is_forward: bool,
    pub cigar: Cigar,
    pub score: i32,
}

/// A builder for [`SubAlignment`]s.
pub struct SubAlignmentBuilder {
    use_eq_and_x: bool,
    match_kind: Kind,
    mismatch_kind: Kind,
    elements: Vec<Op>,
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
                self.elements.push(Op::new(self.match_kind, op_len));
                None
            }
            AlignmentOperation::Subst => {
                self.score += scoring.match_fn.score(b'A', b'C') * (op_len as i32);
                self.query_offset += op_len;
                self.target_offset += op_len;
                self.elements.push(Op::new(self.mismatch_kind, op_len));
                None
            }
            AlignmentOperation::Del => {
                self.score += scoring.gap_open + (scoring.gap_extend * op_len as i32);
                self.target_offset += op_len;
                self.elements.push(Op::new(Kind::Deletion, op_len));
                None
            }
            AlignmentOperation::Ins => {
                self.score += scoring.gap_open + (scoring.gap_extend * op_len as i32);
                self.query_offset += op_len;
                self.elements.push(Op::new(Kind::Insertion, op_len));
                None
            }
            AlignmentOperation::Xskip(new_query_start) => {
                let alignment = SubAlignment {
                    query_start: self.query_start,
                    query_end: self.query_offset,
                    target_start: self.target_start,
                    target_end: self.target_offset,
                    is_forward: self.is_forward,
                    cigar: Cigar::try_from(self.elements.clone()).unwrap(),
                    score: self.score,
                };

                // reset
                self.elements.clear();
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
                    cigar: Cigar::try_from(self.elements.clone()).unwrap(),
                    score: self.score,
                };

                // reset
                self.elements.clear();
                self.target_start = self.target_offset;
                self.query_start = new_query_start;
                self.query_offset = new_query_start;
                self.score = 0;
                self.is_forward = !self.is_forward;

                Some(alignment)
            }
            AlignmentOperation::Yclip(_) => {
                assert!(op_len == 1);
                None
            }
            AlignmentOperation::Xclip(_) => {
                assert!(op_len == 1);
                None
            }
        }
    }

    pub fn new(use_eq_and_x: bool) -> Self {
        Self {
            use_eq_and_x,
            match_kind: if use_eq_and_x {
                Kind::SequenceMatch
            } else {
                Kind::Match
            },
            mismatch_kind: if use_eq_and_x {
                Kind::SequenceMismatch
            } else {
                Kind::Match
            },
            elements: Vec::new(),
            query_start: 0,
            target_start: 0,
            query_offset: 0,
            target_offset: 0,
            score: 0,
            is_forward: true,
        }
    }

    pub fn swap_cigar(cigar: &Cigar) -> Cigar {
        let ops: Vec<Op> = cigar
            .iter()
            .map(|op| match op.kind() {
                Kind::Deletion => Op::new(Kind::Insertion, op.len()),
                Kind::Insertion => Op::new(Kind::Deletion, op.len()),
                _ => *op,
            })
            .collect();
        Cigar::try_from(ops).unwrap()
    }

    pub fn build<F: MatchFunc>(
        &mut self,
        alignment: &Alignment,
        swap: bool,
        scoring: &Scoring<F>,
    ) -> Vec<SubAlignment> {
        self.elements.clear();
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
        } else {
            let alignment = SubAlignment {
                query_start: self.query_start,
                query_end: self.query_offset,
                target_start: self.target_start,
                target_end: self.target_offset,
                is_forward: self.is_forward,
                cigar: Cigar::try_from(self.elements.clone()).unwrap(),
                score: self.score,
            };
            alignments.push(alignment);
        }

        if swap {
            alignments
                .iter()
                .map(|a| SubAlignment {
                    query_start: a.target_start,
                    query_end: a.target_end,
                    target_start: a.query_start,
                    target_end: a.query_end,
                    is_forward: a.is_forward,
                    cigar: Self::swap_cigar(&a.cigar),
                    score: a.score,
                })
                .collect_vec()
        } else {
            alignments
        }
    }
}