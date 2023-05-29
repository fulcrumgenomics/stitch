use serde::Deserialize;
use serde::Serialize;

/// Packed representation of one cell of a Smith-Waterman traceback matrix.  Does not store
/// alignment length for each matrix value, and so cannot be used to break ties to prefer longer
/// alignments.
#[derive(
    Default, Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize,
)]
pub struct SimpleCell {
    // `v` is packed as follows:
    // - bits 0-3 are for I_POS
    // - bits 4-7 are for D_POS
    // - bits 8-11 are for S_POS
    // - bit 12 is strand flip for S
    // - bits 13-15 are unused.
    v: u16,
    // where S is from (jumping)
    s_from: u32,
}

// Traceback bit positions (LSB)
const I_POS: u16 = 0; // Meaning bits 0,1,2,3 corresponds to I and so on
const D_POS: u16 = 4;
const S_POS: u16 = 8;
const S_FLIP: u16 = 12;

impl TracebackCell for SimpleCell {
    /// Initialize a blank traceback cell
    // #[inline(always)]
    // pub fn new() -> Self {
    //     Default::default()
    // }

    /// Sets 4 bits [pos, pos+4) with the 4 LSBs of value
    #[inline(always)]
    fn set_bits(&mut self, pos: u16, tb: u8) {
        let bits: u16 = (0b1111) << pos;
        assert!(
            tb <= TB_MAX,
            "Expected a value <= TB_MAX while setting traceback bits"
        );
        self.v = (self.v & !bits) // First clear the bits
            | (tb << pos) // And set the bits
    }

    #[inline(always)]
    fn set_i(&mut self, tb: u8, len: u32) {
        // Traceback corresponding to matrix I
        self.set_bits(I_POS, tb);
    }

    #[inline(always)]
    fn set_d(&mut self, tb: u8, len: u32) {
        // Traceback corresponding to matrix D
        self.set_bits(D_POS, tb);
    }

    #[inline(always)]
    fn set_s(&mut self, tb: u8, len: u32) {
        // Traceback corresponding to matrix S
        self.set_bits(S_POS, tb);
    }

    #[inline(always)]
    fn set_s_all(&mut self, tb: bool, len: u32, from: u32, strand: bool) {
        // Traceback corresponding to matrix S
        self.set_s(S_POS, tb, len);
        self.s_from = value;
        self.set_strand(S_FLIP, strand);
    }

    #[inline(always)]
    fn set_strand(&mut self, pos: u8, value: bool) {
        let bits: u16 = (0b1) << pos;
        if value {
            self.v = (self.v & !bits) // First clear the bits
            | (0b1 << pos) // And set the bits
        } else {
            self.v = self.v & !bits // Clear the bits
        }
    }

    // Gets 4 bits [pos, pos+4) of v
    #[inline(always)]
    fn get_bits(self, pos: u8) -> u16 {
        (self.v >> pos) & (0b1111)
    }

    #[inline(always)]
    fn get_i(self) -> (u8, u32) {
        (self.get_bits(I_POS), 0)
    }

    #[inline(always)]
    fn get_d(self) -> (u8, u32) {
        (self.get_bits(D_POS), 0)
    }

    #[inline(always)]
    fn get_s(self) -> (u8, u32, u32, bool) {
        (
            self.get_bits(S_POS),
            0,
            self.s_from,
            self.get_s_strand(S_FLIP),
        )
    }

    #[inline(always)]
    fn get_s_strand(self, pos: u8) -> bool {
        (self.v >> S_FLIP) & 01b != 0
    }

    /// Set all matrices to the same value.
    fn set_all(&mut self, tb: u8, len: u32) {
        self.set_i(tb, len);
        self.set_d(tb, len);
        self.set_s(tb, len);
    }
}
