use std::fmt;

use super::aligners::constants::{AlignmentMode, AlignmentOperation};
use crate::align::aligners::constants::{
    AlignmentMode::{Global, QueryLocal, TargetLocal},
    AlignmentOperation::{Del, Ins, Match, Subst, Xclip, Xjump, Yclip, Yjump},
};

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

    /// The start contig index (0-based), corresponding to the aligner for the aligned contig.
    pub start_contig_idx: usize,

    /// The end contig index (0-based), corresponding to the aligner for the aligned contig.
    pub end_contig_idx: usize,

    /// Vector of alignment operations
    pub operations: Vec<AlignmentOperation>,
    pub mode: AlignmentMode,

    /// Alignment length, excluding any clipped (x or y) bases and jumps.
    pub length: usize,
}

#[allow(dead_code)]
impl Alignment {
    // Validate that the xend, yend, end_contig_idx, and length values are correct given the operations.
    pub fn validate(&self) {
        match self.mode {
            AlignmentMode::Global => {
                assert_eq!(self.xstart, 0);
                assert_eq!(self.xend, self.xlen);
                assert_eq!(self.ystart, 0);
                assert_eq!(self.yend, self.ylen);
            }
            AlignmentMode::TargetLocal => {
                assert!(self.xend <= self.xlen);
                assert_eq!(self.ystart, 0);
                assert_eq!(self.yend, self.ylen);
            }
            AlignmentMode::QueryLocal => {
                assert_eq!(self.xstart, 0);
                assert_eq!(self.xend, self.xlen);
                assert!(self.yend <= self.ylen);
            }
            AlignmentMode::Local => {
                assert!(self.xend <= self.xlen);
                assert!(self.yend <= self.ylen);
            }
            _ => (),
        }
        let mut xend: i32 = self.xstart as i32;
        let mut yend: usize = self.ystart;
        let mut end_contig_idx = self.end_contig_idx;
        let mut length = 0;
        for op in &self.operations {
            xend += op.length_on_x(xend as usize);
            yend += op.length_on_y();
            if let Xjump(new_contig_index, _) = op {
                end_contig_idx = *new_contig_index;
            }
            match op {
                Match | Subst | Del | Ins => {
                    length += 1;
                }
                _ => (),
            }
            assert!(xend <= self.xlen as i32);
            assert!(yend <= self.ylen);
        }
        assert_eq!(self.xend, xend as usize, "xend");
        assert_eq!(self.yend, yend, "yend");
        assert_eq!(self.end_contig_idx, end_contig_idx, "end_contig_idx");
        assert_eq!(self.length, length, "length");
    }

    pub fn cigar(&self) -> String {
        let mut cigar: String = String::new();
        if self.operations.is_empty() {
            return cigar;
        }
        let mut contig_idx = self.start_contig_idx;
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
                x_index += op.length_on_x(x_index as usize);
                last_op = op;
                last_len = 0;
                if let Xjump(new_contig_index, _) = op {
                    contig_idx = *new_contig_index;
                }
            } else if op == last_op {
                x_index += op.length_on_x(x_index as usize);
                last_len += 1;
            } else {
                x_index += op.length_on_x(x_index as usize);
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

    /// Returns the 0-based index in x of the earliest base in y that is aligned to the contig with
    /// the given index.
    pub fn earliest_x_base_for(&self, contig_idx: usize) -> Option<usize> {
        if self.operations.is_empty() {
            return None;
        }
        if self.start_contig_idx == contig_idx {
            return Some(self.xstart);
        }
        let mut x_contig_idx = self.start_contig_idx;
        let mut x_index: i32 = self.xstart as i32;
        for op in &self.operations {
            if x_contig_idx == contig_idx {
                return Some(x_index as usize);
            }
            // Update
            if let Xjump(new_contig_index, _) = op {
                x_contig_idx = *new_contig_index;
            }
            x_index += op.length_on_x(x_index as usize);
        }
        None
    }

    /// Returns the 0-based index in x of the latest base in y that is aligned to the contig with
    /// the given index.
    pub fn latest_x_base_for(&self, contig_idx: usize) -> Option<usize> {
        if self.operations.is_empty() {
            return None;
        }
        let mut x_contig_idx = self.start_contig_idx;
        let mut x_index: i32 = self.xstart as i32;
        let mut latest_x_base: Option<usize> = if x_contig_idx == contig_idx {
            Some(self.xstart)
        } else {
            None
        };
        for op in &self.operations {
            // Update the current contig index before checking if the latest x_base should be
            // updated.
            if let Xjump(new_contig_index, _) = op {
                x_contig_idx = *new_contig_index;
            }
            if x_contig_idx == contig_idx {
                latest_x_base = Some(x_index as usize);
            }
            x_index += op.length_on_x(x_index as usize);
        }
        latest_x_base
    }

    /// Splits the alignment into two halves, one half aligned up to `y_pivot` point, and the other
    /// half after and including the `y_pivot` point, then swaps their order, and joins them.
    ///
    /// Arguments:
    /// - `y_pivot` - The offset of the pivot point relative to the y used in this alignment.
    pub fn split_at_y(&self, y_pivot: usize) -> Alignment {
        if self.operations.is_empty() {
            return self.clone();
        }

        // There should be no leading or trailing clips
        assert!(!matches!(self.operations[0], Xclip(_) | Yclip(_)));
        assert!(!matches!(
            self.operations[self.operations.len() - 1],
            Xclip(_) | Yclip(_)
        ));

        let mut x_index: usize = self.xstart;
        let mut y_index: usize = self.ystart;
        let mut contig_index: usize = self.start_contig_idx;
        let mut op_index = 0;

        // skip over any X/Y clip/jumps at the start of the alignment
        for op in &self.operations {
            match op {
                Match | Subst | Del | Ins => break,
                _ => (),
            }

            if let Xjump(idx, _) = op {
                contig_index = *idx;
            }
            y_index += op.length_on_y();
            x_index = (x_index as i32 + op.length_on_x(x_index)) as usize;
            op_index += 1;
        }

        // build the alignment up to the pivot point
        for op in &self.operations[op_index..] {
            if y_index + op.length_on_y() >= y_pivot {
                break;
            }
            if let Xjump(idx, _) = op {
                contig_index = *idx;
            }
            y_index += op.length_on_y();
            x_index = (x_index as i32 + op.length_on_x(x_index)) as usize;
            op_index += 1;
        }
        let pre_pivot_aln = Alignment {
            xstart: self.xstart,
            xend: x_index + 1,
            ystart: self.ystart,
            yend: y_index + 1,
            xlen: 0,
            ylen: 0,
            start_contig_idx: self.start_contig_idx,
            end_contig_idx: contig_index,
            operations: self.operations[..=op_index].to_vec(),
            mode: self.mode,
            score: 0,
            length: 0,
        };
        assert!(y_pivot >= pre_pivot_aln.yend);

        // skip over any X/Y clip/jumps at the pivot point
        for op in &self.operations[op_index..] {
            if y_index >= y_pivot {
                match op {
                    Match | Subst | Del | Ins => break,
                    _ => (),
                }
            }
            if let Xjump(idx, _) = op {
                contig_index = *idx;
            }
            y_index += op.length_on_y();
            x_index = (x_index as i32 + op.length_on_x(x_index)) as usize;
            op_index += 1;
        }

        // build the alignment after the pivot point
        let post_pivot_aln = Alignment {
            xstart: x_index,
            xend: self.xend,
            ystart: y_index,
            yend: self.yend,
            xlen: 0,
            ylen: 0,
            start_contig_idx: contig_index,
            end_contig_idx: self.end_contig_idx,
            operations: self.operations[op_index..].to_vec(),
            mode: self.mode,
            score: 0,
            length: 0,
        };

        // join the two alignments
        let mut aln = Alignment {
            start_contig_idx: post_pivot_aln.start_contig_idx,
            end_contig_idx: pre_pivot_aln.end_contig_idx,
            xstart: post_pivot_aln.xstart,
            ystart: post_pivot_aln.ystart - y_pivot,
            xend: pre_pivot_aln.xend,
            yend: pre_pivot_aln.yend + self.ylen - y_pivot,
            ylen: self.ylen,
            xlen: self.xlen,
            score: self.score,
            operations: Vec::new(),
            mode: self.mode,
            length: self.length,
        };

        // True if we are to add prefix/suffix clipping to x/y respectively.
        let x_clip = matches!(aln.mode, Global | QueryLocal);
        let y_clip = matches!(aln.mode, Global | TargetLocal);

        // Add any X/Y prefix clipping
        if x_clip && aln.xstart > 0 {
            aln.operations.push(Xclip(aln.xstart));
            aln.xstart = 0;
        }
        if y_clip && aln.ystart > 0 {
            aln.operations.push(Yclip(aln.ystart));
            aln.ystart = 0;
        }

        // Add the post pivot alignment
        aln.operations.extend_from_slice(&post_pivot_aln.operations);

        // Add any Xjump
        if pre_pivot_aln.start_contig_idx != post_pivot_aln.end_contig_idx
            || pre_pivot_aln.xstart != post_pivot_aln.xend
        {
            aln.operations
                .push(Xjump(pre_pivot_aln.start_contig_idx, pre_pivot_aln.xstart));
        }

        // Add any Yjump
        let yjump_len = aln.ylen + pre_pivot_aln.ystart - post_pivot_aln.yend;
        if yjump_len > 0 {
            aln.operations.push(Yjump(yjump_len));
        }

        // Add the pre-pivot alignments
        aln.operations.extend_from_slice(&pre_pivot_aln.operations);

        // Add any X/Y suffix clipping
        if x_clip && aln.xend < aln.xlen {
            aln.operations.push(Xclip(aln.xlen - aln.xend));
            aln.xend = aln.xlen;
        }
        if y_clip && aln.yend < aln.ylen {
            aln.operations.push(Xclip(aln.ylen - aln.yend));
            aln.yend = aln.ylen;
        }

        aln
    }
}

impl fmt::Display for Alignment {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "contig-idx: {}-{} x-span: {}-{}/{} y-span: {}-{}/{} score: {} cigar: {} aln-len: {}",
            self.start_contig_idx,
            self.end_contig_idx,
            self.xstart,
            self.xend,
            self.xlen,
            self.ystart,
            self.yend,
            self.ylen,
            self.score,
            self.cigar(),
            self.length
        )
    }
}

#[cfg(test)]
pub mod tests {
    use rstest::rstest;

    use crate::align::aligners::constants::{
        AlignmentMode,
        AlignmentMode::Local,
        AlignmentOperation::{Del, Ins, Match, Subst, Xjump, Yclip, Yjump},
    };

    use super::Alignment;

    fn empty_alignment() -> Alignment {
        Alignment {
            score: 0,
            xstart: 0,
            xend: 0,
            ystart: 0,
            yend: 0,
            ylen: 0,
            xlen: 0,
            start_contig_idx: 0,
            end_contig_idx: 0,
            operations: Vec::new(),
            mode: AlignmentMode::Global,
            length: 0,
        }
    }

    fn non_empty_alignment() -> Alignment {
        Alignment {
            score: 0,
            xstart: 10,
            xend: 110,
            xlen: 110,
            ystart: 11,
            yend: 111,
            ylen: 111,
            start_contig_idx: 0,
            end_contig_idx: 0,
            operations: vec![Match; 100],
            mode: AlignmentMode::Local,
            length: 100,
        }
    }

    fn single_jump_alignment() -> Alignment {
        Alignment {
            score: 0,
            xstart: 10,
            xend: 4,
            xlen: 12,
            ystart: 11,
            yend: 15,
            ylen: 15,
            start_contig_idx: 0,
            end_contig_idx: 1,
            operations: [Match, Match, Xjump(1, 2), Match, Match].to_vec(),
            mode: AlignmentMode::Local,
            length: 4,
        }
    }

    fn double_jump_alignment() -> Alignment {
        Alignment {
            score: 0,
            xstart: 10,
            xend: 10,
            xlen: 12,
            ystart: 11,
            yend: 17,
            ylen: 17,
            start_contig_idx: 0,
            end_contig_idx: 0,
            operations: [
                Match,
                Match,
                Xjump(1, 2),
                Match,
                Match,
                Xjump(0, 8),
                Match,
                Match,
            ]
            .to_vec(),
            mode: AlignmentMode::Local,
            length: 6,
        }
    }

    fn jump_backwards() -> Alignment {
        Alignment {
            score: 0,
            xstart: 2,
            xend: 2,
            xlen: 4,
            ystart: 0,
            yend: 4,
            ylen: 4,
            start_contig_idx: 0,
            end_contig_idx: 0,
            operations: [Match, Match, Xjump(0, 0), Match, Match].to_vec(),
            mode: AlignmentMode::Local,
            length: 4,
        }
    }

    fn all_ops_alignmnent() -> Alignment {
        Alignment {
            score: 0,
            xstart: 10,
            xend: 7,
            xlen: 16,
            ystart: 11,
            yend: 28,
            ylen: 28,
            start_contig_idx: 0,
            end_contig_idx: 3,
            operations: [
                Match,       // [contig: 0, x: 10, y: 11] -> [contig: 0, x: 11, y: 12]
                Match,       // [contig: 0, x: 11, y: 12] -> [contig: 0, x: 12, y: 13]
                Xjump(1, 2), // [contig: 0, x: 12, y: 13] -> [contig: 1, x: 2,  y: 13]
                Match,       // [contig: 1, x: 2,  y: 13] -> [contig: 1, x: 3,  y: 14]
                Match,       // [contig: 1, x: 3,  y: 14] -> [contig: 1, x: 4,  y: 15]
                Xjump(0, 8), // [contig: 1, x: 4,  y: 15] -> [contig: 0, x: 8,  y: 15]
                Match,       // [contig: 1, x: 4,  y: 15] -> [contig: 0, x: 9,  y: 16]
                Match,       // [contig: 0, x: 9,  y: 16] -> [contig: 0, x: 10, y: 17]
                Subst,       // [contig: 0, x: 10, y: 17] -> [contig: 0, x: 11, y: 18]
                Yjump(3),    // [contig: 0, x: 11, y: 18] -> [contig: 0, x: 11, y: 21]
                Match,       // [contig: 0, x: 11, y: 21] -> [contig: 0, x: 12, y: 22]
                Ins,         // [contig: 0, x: 12, y: 22] -> [contig: 0, x: 13, y: 22]
                Ins,         // [contig: 0, x: 13, y: 22] -> [contig: 0, x: 14, y: 22]
                Ins,         // [contig: 0, x: 14, y: 22] -> [contig: 0, x: 15, y: 22]
                Match,       // [contig: 0, x: 15, y: 22] -> [contig: 0, x: 16, y: 23]
                Xjump(3, 4), // [contig: 0, x: 16, y: 23] -> [contig: 3, x:  4, y: 23]
                Subst,       // [contig: 3, x:  4, y: 23] -> [contig: 3, x:  5, y: 24]
                Match,       // [contig: 3, x:  5, y: 24] -> [contig: 3, x:  6, y: 25]
                Del,         // [contig: 3, x:  6, y: 25] -> [contig: 3, x:  6, y: 26]
                Del,         // [contig: 3, x:  6, y: 26] -> [contig: 3, x:  6, y: 27]
                Match,       // [contig: 3, x:  6, y: 27] -> [contig: 3, x:  7, y: 28]
            ]
            .to_vec(),
            mode: AlignmentMode::Local,
            length: 17,
        }
    }

    #[rstest]
    #[case(&empty_alignment())]
    #[case(&non_empty_alignment())]
    #[case(&single_jump_alignment())]
    #[case(&double_jump_alignment())]
    #[case(&jump_backwards())]
    #[case(&all_ops_alignmnent())]
    #[case(&test_no_y_jump())]
    #[case(&test_slop_5_on_x())]
    #[case(&test_slop_5_on_x_with_y_clipping(Local))]
    fn test_valid_alignments(#[case] alignment: &Alignment) {
        alignment.validate();
    }

    #[rstest]
    #[case(&empty_alignment(), 0, None)]
    #[case(&non_empty_alignment(), 0, Some(10))]
    #[case(&non_empty_alignment(), 1, None)]
    #[case(&single_jump_alignment(), 0, Some(10))]
    #[case(&single_jump_alignment(), 1, Some(2))]
    #[case(&double_jump_alignment(), 0, Some(10))]
    #[case(&double_jump_alignment(), 1, Some(2))]
    #[case(&jump_backwards(), 0, Some(2))]
    #[case(&all_ops_alignmnent(), 0, Some(10))]
    #[case(&all_ops_alignmnent(), 1, Some(2))]
    #[case(&all_ops_alignmnent(), 2, None)]
    #[case(&all_ops_alignmnent(), 3, Some(4))]
    fn test_earliest_x_base(
        #[case] alignment: &Alignment,
        #[case] contig_idx: usize,
        #[case] x: Option<usize>,
    ) {
        assert_eq!(alignment.earliest_x_base_for(contig_idx), x);
    }

    #[rstest]
    #[case(&empty_alignment(),0,  None)]
    #[case(&non_empty_alignment(), 0, Some(109))]
    #[case(&non_empty_alignment(), 1, None)]
    #[case(&single_jump_alignment(), 0, Some(11))]
    #[case(&single_jump_alignment(), 1, Some(3))]
    #[case(&double_jump_alignment(), 0, Some(9))]
    #[case(&double_jump_alignment(), 1, Some(3))]
    #[case(&jump_backwards(), 0, Some(1))]
    #[case(&all_ops_alignmnent(), 0, Some(15))]
    #[case(&all_ops_alignmnent(), 1, Some(3))]
    #[case(&all_ops_alignmnent(), 2, None)]
    #[case(&all_ops_alignmnent(), 3, Some(6))]
    fn test_latest_x_base_for(
        #[case] alignment: &Alignment,
        #[case] contig_idx: usize,
        #[case] x: Option<usize>,
    ) {
        assert_eq!(alignment.latest_x_base_for(contig_idx), x);
    }

    fn test_no_y_jump() -> Alignment {
        // aligns to the end of the target, then jumps to the start of the target, and aligns again
        Alignment {
            score: 0,
            xstart: 45,
            xend: 5,
            xlen: 50,
            ystart: 0,
            yend: 10,
            ylen: 10,
            start_contig_idx: 0,
            end_contig_idx: 0,
            operations: [
                Match,
                Match,
                Match,
                Match,
                Match,
                Xjump(0, 0),
                Match,
                Match,
                Match,
                Match,
                Match,
            ]
            .to_vec(),
            mode: AlignmentMode::Local,
            length: 10,
        }
    }

    fn test_slop_5_on_x() -> Alignment {
        Alignment {
            score: 0,
            xstart: 40,
            xend: 10,
            xlen: 50,
            ystart: 0,
            yend: 10,
            ylen: 10,
            start_contig_idx: 0,
            end_contig_idx: 0,
            operations: [
                Match,
                Match,
                Match,
                Match,
                Match,
                Xjump(0, 5),
                Match,
                Match,
                Match,
                Match,
                Match,
            ]
            .to_vec(),
            mode: AlignmentMode::Local,
            length: 10,
        }
    }

    fn test_slop_5_on_x_with_y_clipping(mode: AlignmentMode) -> Alignment {
        Alignment {
            score: 0,
            xstart: 40,
            xend: 10,
            xlen: 50,
            ystart: 0,
            yend: 20,
            ylen: 20,
            start_contig_idx: 0,
            end_contig_idx: 0,
            operations: [
                Match,
                Match,
                Match,
                Match,
                Match,
                Yclip(5),
                Xjump(0, 5),
                Yclip(5),
                Match,
                Match,
                Match,
                Match,
                Match,
            ]
            .to_vec(),
            mode,
            length: 10,
        }
    }

    #[rstest]
    #[case(&empty_alignment(), 0, 0, 0, 0, 0, 0, &String::new(), 0)]
    #[case(&test_no_y_jump(), 5, 0, 50, 0, 10, 0, "5=40J5=", 10)]
    #[case(&test_slop_5_on_x(), 5, 5, 45, 0, 10, 0, "5=30J5=", 10)]
    #[case(&test_slop_5_on_x_with_y_clipping(AlignmentMode::Global), 5, 0, 50, 0, 20, 0, "5A10B5=30J5=5A", 10)]
    #[case(&test_slop_5_on_x_with_y_clipping(AlignmentMode::Local), 5, 5, 45, 10, 20, 0, "5=30J5=", 10)]
    #[case(&test_slop_5_on_x_with_y_clipping(AlignmentMode::TargetLocal), 5, 5, 45, 0, 20, 0, "10B5=30J5=", 10)]
    #[case(&test_slop_5_on_x_with_y_clipping(AlignmentMode::QueryLocal), 5, 0, 50, 10, 20, 0, "5A5=30J5=5A", 10)]
    fn test_split_at_y(
        #[case] alignment: &Alignment,
        #[case] y_pivot: usize,
        #[case] xstart: usize,
        #[case] xend: usize,
        #[case] ystart: usize,
        #[case] yend: usize,
        #[case] score: i32,
        #[case] cigar: &str,
        #[case] length: usize,
    ) {
        let alignment = alignment.split_at_y(y_pivot);
        assert_eq!(alignment.xstart, xstart, "xstart {alignment}");
        assert_eq!(alignment.xend, xend, "xend {alignment}");
        assert_eq!(alignment.ystart, ystart, "ystart {alignment}");
        assert_eq!(alignment.yend, yend, "yend {alignment}");
        assert_eq!(alignment.score, score, "score {alignment}");
        assert_eq!(alignment.start_contig_idx, 0, "strand {alignment}");
        assert_eq!(alignment.cigar(), cigar, "cigar {alignment}");
        assert_eq!(alignment.length, length, "length {alignment}");
    }
}
