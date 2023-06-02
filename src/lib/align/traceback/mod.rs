use super::{
    aligners::{
        constants::{AlignmentMode, AlignmentOperation},
        single_strand::SingleStrandAligner,
    },
    alignment::Alignment,
};
use bio::alignment::pairwise::MatchFunc;
use serde::Deserialize;
use serde::Serialize;

#[derive(Default, Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug)]
pub struct SValue {
    pub tb: u8,
    pub len: u32,
    pub from: u32,
    pub flip_strand: bool,
}

pub trait TracebackCell: Clone {
    fn set_i(&mut self, tb: u8, len: u32);
    fn set_d(&mut self, tb: u8, len: u32);
    fn set_s(&mut self, tb: u8, len: u32);
    fn set_s_all(&mut self, tb: u8, len: u32, from: u32, strand: bool);
    fn get_i(self) -> (u8, u32);
    fn get_d(self) -> (u8, u32);
    fn get_s(self) -> SValue;
    fn set_all(&mut self, tb: u8, len: u32);

    // fn get_s_tb(self) -> u8;
    fn get_s_len(self) -> u32;
    fn get_i_len(self) -> u32;
    fn get_d_len(self) -> u32;
}

// Traceback moves
pub const TB_START: u8 = 0b0000;
pub const TB_INS: u8 = 0b0001; // 1
pub const TB_DEL: u8 = 0b0010; // 2
pub const TB_SUBST: u8 = 0b0011; // 3
pub const TB_MATCH: u8 = 0b0100; // 4
pub const TB_XCLIP_PREFIX: u8 = 0b0101; // prefix clip of x (5)
pub const TB_XCLIP_SUFFIX: u8 = 0b0110; // suffix clip of x (6)
pub const TB_YCLIP_PREFIX: u8 = 0b0111; // prefix clip of y (7)
pub const TB_YCLIP_SUFFIX: u8 = 0b1000; // suffix clip of y (8)
pub const TB_XJUMP: u8 = 0b1001; // jump (9)
pub const TB_MAX: u8 = 0b1001; // Useful in checking that the TB value we got is a valid one

cfg_if::cfg_if! {
    if #[cfg(low_mem)] {
        pub mod simple_cell;
        pub type Cell = simple_cell::SimpleCell;
    } else {
        pub mod packed_length_cell;
        pub type Cell = packed_length_cell::PackedLengthCell;
    }
}

pub fn default() -> Cell {
    Cell::default()
}

/// Internal traceback.
#[derive(Default, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize)]
pub struct Traceback {
    rows: usize,
    cols: usize,
    matrix: Vec<Cell>,
}

impl Traceback {
    pub fn with_capacity(m: usize, n: usize) -> Self {
        let rows = m + 1;
        let cols = n + 1;
        Traceback {
            rows,
            cols,
            matrix: Vec::with_capacity(rows * cols),
        }
    }

    pub fn init(&mut self, m: usize, n: usize) {
        self.matrix.clear();
        let mut start = crate::align::traceback::default();
        start.set_all(TB_START, 0);
        start.set_s_all(TB_START, 0, 0, false);
        // set every cell to start
        self.resize(m, n, start);
    }

    #[inline(always)]
    pub fn set(&mut self, i: usize, j: usize, v: Cell) {
        debug_assert!(i < self.rows);
        debug_assert!(j < self.cols);
        self.matrix[i * self.cols + j] = v;
    }

    #[inline(always)]
    pub fn get(&self, i: usize, j: usize) -> &Cell {
        debug_assert!(i < self.rows);
        debug_assert!(j < self.cols);
        &self.matrix[i * self.cols + j]
    }

    pub fn get_mut(&mut self, i: usize, j: usize) -> &mut Cell {
        debug_assert!(i < self.rows);
        debug_assert!(j < self.cols);
        &mut self.matrix[i * self.cols + j]
    }

    pub fn resize(&mut self, m: usize, n: usize, v: Cell) {
        self.rows = m + 1;
        self.cols = n + 1;
        self.matrix.resize(self.rows * self.cols, v);
    }
}

pub fn traceback_single_stranded<F: MatchFunc>(
    forward: &SingleStrandAligner<F>,
    m: usize,
    n: usize,
) -> Alignment {
    traceback_double_stranded(forward, forward, m, n)
}

pub fn traceback_double_stranded<F: MatchFunc>(
    forward: &SingleStrandAligner<F>,
    reverse: &SingleStrandAligner<F>,
    m: usize,
    n: usize,
) -> Alignment {
    let mut i = m;
    let mut j = n;
    let mut operations: Vec<AlignmentOperation> = Vec::with_capacity(m);
    let mut xstart: usize = 0usize;
    let mut ystart: usize = 0usize;
    let mut xend = m;
    let mut yend = n;

    // If the scores equal, pick the one with the longer alignment length
    let mut cur_is_forward = match forward.S[n % 2][m].cmp(&reverse.S[n % 2][m]) {
        std::cmp::Ordering::Less => false,
        std::cmp::Ordering::Greater => true,
        std::cmp::Ordering::Equal => {
            let fwd_len = forward.traceback.get(i, j).get_s_len();
            let rev_len = reverse.traceback.get(i, j).get_s_len();
            fwd_len >= rev_len
        }
    };
    let mut cur_aligner = if cur_is_forward { &forward } else { &reverse };
    let alignment_length = cur_aligner.traceback.get(i, j).get_s_len();
    let score = cur_aligner.S[n % 2][m];

    let mut last_layer = cur_aligner.traceback.get(i, j).get_s().tb;
    loop {
        cur_aligner = if cur_is_forward { &forward } else { &reverse };
        let next_layer: u8;
        match last_layer {
            TB_START => break,
            TB_INS => {
                operations.push(AlignmentOperation::Ins);
                next_layer = cur_aligner.traceback.get(i, j).get_i().0;
                i -= 1;
            }
            TB_DEL => {
                operations.push(AlignmentOperation::Del);
                next_layer = cur_aligner.traceback.get(i, j).get_d().0;
                j -= 1;
            }
            TB_MATCH | TB_SUBST => {
                if last_layer == TB_MATCH {
                    operations.push(AlignmentOperation::Match);
                } else {
                    operations.push(AlignmentOperation::Subst);
                }
                let s_value = cur_aligner.traceback.get(i, j).get_s();
                let s_from = s_value.from as usize;
                if s_value.flip_strand {
                    operations.push(AlignmentOperation::Xflip(i - 1));
                    cur_is_forward = !cur_is_forward;
                    cur_aligner = if cur_is_forward { &forward } else { &reverse };
                } else if s_from != i - 1 {
                    operations.push(AlignmentOperation::Xskip(i - 1));
                }
                next_layer = cur_aligner.traceback.get(s_from, j - 1).get_s().tb;
                i = s_from;
                j -= 1;
            }
            TB_XCLIP_PREFIX => {
                next_layer = cur_aligner.traceback.get(0, j).get_s().tb;
                // only add Xclip if there are only clip moves left, since we may have jumped!
                if next_layer == TB_START || next_layer == TB_YCLIP_PREFIX {
                    operations.push(AlignmentOperation::Xclip(i));
                    xstart = i;
                }
                i = 0;
            }
            TB_XCLIP_SUFFIX => {
                if operations.is_empty()
                    || matches!(operations.first().unwrap(), AlignmentOperation::Yclip(_))
                {
                    operations.push(AlignmentOperation::Xclip(cur_aligner.Lx[j]));
                    xend = i - cur_aligner.Lx[j];
                }
                i -= cur_aligner.Lx[j];
                next_layer = cur_aligner.traceback.get(i, j).get_s().tb;
            }
            TB_YCLIP_PREFIX => {
                operations.push(AlignmentOperation::Yclip(j));
                ystart = j;
                j = 0;
                next_layer = cur_aligner.traceback.get(i, 0).get_s().tb;
            }
            TB_YCLIP_SUFFIX => {
                operations.push(AlignmentOperation::Yclip(cur_aligner.Ly[i]));
                let s_from = cur_aligner.traceback.get(i, j).get_s().from as usize;
                j -= cur_aligner.Ly[i];
                if s_from != i {
                    operations.push(AlignmentOperation::Xskip(i));
                    i = s_from;
                }
                yend = j;
                next_layer = cur_aligner.traceback.get(i, j).get_s().tb;
            }
            TB_XJUMP => {
                let s_value = cur_aligner.traceback.get(i, j).get_s();
                let s_from;
                if s_value.flip_strand {
                    operations.push(AlignmentOperation::Xflip(i));
                    cur_is_forward = !cur_is_forward;
                    cur_aligner = if cur_is_forward { &forward } else { &reverse };
                    s_from = s_value.from;
                } else {
                    operations.push(AlignmentOperation::Xskip(i));
                    s_from = s_value.from;
                }
                next_layer = cur_aligner
                    .traceback
                    .get(s_value.from as usize, j)
                    .get_s()
                    .tb;
                i = s_from as usize;
            }
            _ => panic!("Dint expect this!"),
        }
        last_layer = next_layer;
    }

    operations.reverse();
    {
        use AlignmentOperation::{Xclip, Xskip, Yclip};
        if operations
            .iter()
            .all(|op| matches!(op, Xclip(_) | Yclip(_) | Xskip(_)))
        {
            xstart = 0;
            xend = 0;
            ystart = 0;
            yend = 0;
        }
    }
    Alignment {
        score,
        ystart,
        xstart,
        yend,
        xend,
        xlen: m,
        ylen: n,
        is_forward: cur_is_forward,
        operations,
        mode: AlignmentMode::Custom,
        length: alignment_length as usize,
    }
}
