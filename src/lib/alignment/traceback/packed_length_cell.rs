use serde::Deserialize;
use serde::Serialize;

use crate::alignment::traceback::TB_MAX;

use super::SValue;
use super::TracebackCell;

/// Packed representation of one cell of a Smith-Waterman traceback matrix. Stores the
/// alignment length for each matrix value, and so can be used to break ties to prefer longer
/// alignments.
///
/// each matrix value (u32) is packed as follows:
/// - bits 0-3 are reserved for *_POS
/// - bit 4 is reserved for strand flip (used for S_POS only, so 2 unused)
/// - bits 5-32 are reserved for alignment length (27 bits each)
/// Also contains s_from (u32), so we have a totoal of u32 * 4 = 128 bits
#[derive(
    Default, Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize,
)]
pub struct PackedLengthCell {
    s: u32,
    i: u32,
    d: u32,
    s_from: u32,
}

// Traceback bit positions (LSB)
const TRACE_BIT_OFFSET: u8 = 0;
const STRAND_BIT_OFFSET: u8 = 4;
const ALIGN_LEN_BIT_OFFSET: u8 = 5;

impl PackedLengthCell {
    /// Sets 4 bits [pos, pos+4) with the 4 LSBs of value
    #[inline(always)]
    fn set_tb(&mut self, matrix: u32, tb: u8) -> u32 {
        let bits: u32 = (0b1111) << TRACE_BIT_OFFSET;
        assert!(
            tb <= TB_MAX,
            "Expected a value <= TB_MAX while setting traceback bits"
        );
        (matrix & !bits) | ((tb as u32) << TRACE_BIT_OFFSET) // And set the bits
    }

    #[inline(always)]
    fn set_strand(&mut self, matrix: u32, value: bool) -> u32 {
        let bits: u32 = (0b1) << STRAND_BIT_OFFSET;
        if value {
            (matrix & !bits) | (0b1 << STRAND_BIT_OFFSET) // Clear the bit and then set it
        } else {
            matrix & !bits // Clear the bit
        }
    }

    #[inline(always)]
    fn set_len(&mut self, matrix: u32, len: u32) -> u32 {
        let bits: u32 = (0b111_1111_1111_1111_1111_1111_1111) << ALIGN_LEN_BIT_OFFSET;
        (matrix & !bits) | (len << ALIGN_LEN_BIT_OFFSET) // And set the bits
    }

    // Gets 4 bits [pos, pos+4) of v
    #[inline(always)]
    fn get_tb(self, matrix: u32) -> u8 {
        ((matrix >> TRACE_BIT_OFFSET) & (0b1111)) as u8
    }

    #[inline(always)]
    fn get_strand(self, matrix: u32) -> bool {
        (matrix >> STRAND_BIT_OFFSET) & (0b1) != 0
    }

    #[inline(always)]
    fn get_len(self, matrix: u32) -> u32 {
        matrix >> ALIGN_LEN_BIT_OFFSET
    }
}

impl TracebackCell for PackedLengthCell {
    #[inline(always)]
    fn set_i(&mut self, tb: u8, len: u32) {
        // Traceback corresponding to matrix I
        self.i = self.set_tb(self.i, tb);
        self.i = self.set_len(self.i, len);
    }

    #[inline(always)]
    fn set_d(&mut self, tb: u8, len: u32) {
        // Traceback corresponding to matrix D
        self.d = self.set_tb(self.d, tb);
        self.d = self.set_len(self.d, len);
    }

    #[inline(always)]
    fn set_s(&mut self, tb: u8, len: u32) {
        // Traceback corresponding to matrix S
        self.s = self.set_tb(self.s, tb);
        self.s = self.set_len(self.s, len);
    }

    #[inline(always)]
    fn set_s_all(&mut self, tb: u8, len: u32, from: u32, strand: bool) {
        // Traceback corresponding to matrix S
        self.set_s(tb, len);
        self.s = self.set_strand(self.s, strand);
        self.s_from = from;
    }

    #[inline(always)]
    fn get_i_len(self) -> u32 {
        self.get_len(self.i)
    }

    #[inline(always)]
    fn get_i(self) -> (u8, u32) {
        (self.get_tb(self.i), self.get_len(self.i))
    }

    #[inline(always)]
    fn get_d_len(self) -> u32 {
        self.get_len(self.d)
    }

    #[inline(always)]
    fn get_d(self) -> (u8, u32) {
        (self.get_tb(self.d), self.get_len(self.d))
    }

    #[inline(always)]
    fn get_s_len(self) -> u32 {
        self.get_len(self.s)
    }

    #[inline(always)]
    fn get_s(self) -> SValue {
        SValue {
            tb: self.get_tb(self.s),
            len: self.get_s_len(),
            from: self.s_from,
            strand: self.get_strand(self.s),
        }
    }

    /// Set all matrices to the same value.
    fn set_all(&mut self, tb: u8, len: u32) {
        self.set_i(tb, len);
        self.set_d(tb, len);
        self.set_s(tb, len);
    }
}
