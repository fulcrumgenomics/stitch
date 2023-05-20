use serde::Deserialize;
use serde::Serialize;

/// Packed representation of one cell of a Smith-Waterman traceback matrix.
/// Stores the I, D and S traceback matrix values in two bytes.
/// Possible traceback moves include : start, insert, delete, match, substitute,
/// prefix clip and suffix clip for x & y. So we need 4 bits each for matrices I, D, S
/// to keep track of these 9 moves.
#[derive(
    Default, Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize,
)]
pub struct TracebackCell {
    v: u32,
}

// Traceback bit positions (LSB)
pub const I_POS: u8 = 0; // Meaning bits 0,1,2,3 corresponds to I and so on
pub const D_POS: u8 = 4;
pub const S_POS: u8 = 8;
pub const L_POS: u8 = 12;

// Traceback moves
pub const TB_START: u32 = 0b0000;
pub const TB_INS: u32 = 0b0001;
pub const TB_DEL: u32 = 0b0010;
pub const TB_SUBST: u32 = 0b0011;
pub const TB_MATCH: u32 = 0b0100;

pub const TB_XCLIP_PREFIX: u32 = 0b0101; // prefix clip of x
pub const TB_XCLIP_SUFFIX: u32 = 0b0110; // suffix clip of x
pub const TB_YCLIP_PREFIX: u32 = 0b0111; // prefix clip of y
pub const TB_YCLIP_SUFFIX: u32 = 0b1000; // suffix clip of y

pub const TB_MAX: u32 = 0b1000; // Useful in checking that the
                                // TB value we got is a valid one

impl TracebackCell {
    /// Initialize a blank traceback cell
    #[inline(always)]
    pub fn new() -> TracebackCell {
        Default::default()
    }

    /// Sets 4 bits [pos, pos+4) with the 4 LSBs of value
    #[inline(always)]
    fn set_bits(&mut self, pos: u8, value: u32) {
        let bits: u32 = (0b1111) << pos;
        assert!(
            value <= TB_MAX,
            "Expected a value <= TB_MAX while setting traceback bits"
        );
        self.v = (self.v & !bits) // First clear the bits
            | (value << pos) // And set the bits
    }

    #[inline(always)]
    pub fn set_i_bits(&mut self, value: u32) {
        // Traceback corresponding to matrix I
        self.set_bits(I_POS, value);
    }

    #[inline(always)]
    pub fn set_d_bits(&mut self, value: u32) {
        // Traceback corresponding to matrix D
        self.set_bits(D_POS, value);
    }

    #[inline(always)]
    pub fn set_s_bits(&mut self, value: u32) {
        // Traceback corresponding to matrix S
        self.set_bits(S_POS, value);
    }

    #[inline(always)]
    pub fn set_l_bits(&mut self, value: u32) {
        // Traceback corresponding to matrix S
        let bits: u32 = (0b1111) << L_POS;
        self.v = (self.v & !bits) // First clear the bits
            | (value << L_POS) // And set the bits
    }

    // Gets 4 bits [pos, pos+4) of v
    #[inline(always)]
    fn get_bits(self, pos: u8) -> u32 {
        (self.v >> pos) & (0b1111)
    }

    #[inline(always)]
    pub fn get_i_bits(self) -> u32 {
        self.get_bits(I_POS)
    }

    #[inline(always)]
    pub fn get_d_bits(self) -> u32 {
        self.get_bits(D_POS)
    }

    #[inline(always)]
    pub fn get_s_bits(self) -> u32 {
        self.get_bits(S_POS)
    }

    #[inline(always)]
    pub fn get_l_bits(self) -> u32 {
        self.v >> L_POS
    }

    /// Set all matrices to the same value.
    pub fn set_all(&mut self, value: u32) {
        self.set_i_bits(value);
        self.set_d_bits(value);
        self.set_s_bits(value);
    }
}

/// Internal traceback.
#[derive(Default, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize)]
pub struct Traceback {
    rows: usize,
    cols: usize,
    matrix: Vec<TracebackCell>,
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
        let mut start = TracebackCell::new();
        start.set_all(TB_START);
        start.set_l_bits(0);
        // set every cell to start
        self.resize(m, n, start);
    }

    #[inline(always)]
    pub fn set(&mut self, i: usize, j: usize, v: TracebackCell) {
        debug_assert!(i < self.rows);
        debug_assert!(j < self.cols);
        self.matrix[i * self.cols + j] = v;
    }

    #[inline(always)]
    pub fn get(&self, i: usize, j: usize) -> &TracebackCell {
        debug_assert!(i < self.rows);
        debug_assert!(j < self.cols);
        &self.matrix[i * self.cols + j]
    }

    pub fn get_mut(&mut self, i: usize, j: usize) -> &mut TracebackCell {
        debug_assert!(i < self.rows);
        debug_assert!(j < self.cols);
        &mut self.matrix[i * self.cols + j]
    }

    pub fn resize(&mut self, m: usize, n: usize, v: TracebackCell) {
        self.rows = m + 1;
        self.cols = n + 1;
        self.matrix.resize(self.rows * self.cols, v);
    }
}
