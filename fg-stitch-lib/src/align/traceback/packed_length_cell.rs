use serde::{Deserialize, Serialize};

use crate::align::traceback::TB_MAX;

use super::{SValue, TracebackCell};

/// Packed representation of one cell of a Smith-Waterman traceback matrix. Stores the
/// alignment length for each matrix value, and so can be used to break ties to prefer longer
/// alignments.
///
/// The maximum target (x) length is 2^27 - 1 = 134,217,727bp.  The maximum supported number of
/// contigs is 2^8 - 1 = 255 contigs.
///
/// Each matrix value (s/i/d as a u32) is packed as follows:
/// - bits 0-3 are reserved for *_POS
/// - bits 4-31 are reserved for alignment length (27 bits each)
/// - bit 31 for each of s/i/d is reserved for the upper 3-bits of the jump contig index (i.e. which contig did we jump from)
/// The `aux` field is packed as follows:
/// - bits 0-4 are for the lower 5-bits of the jump contig index (i.e. which contig did we jump from)
/// - bits 5-31 are for the "from" jump index in the contig (i.e. where in the given contig did we jump from)
/// Also contains s_from (u32), so we have a totoal of u32 * 4 = 128 bits
#[derive(
    Default, Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize,
)]
pub struct PackedLengthCell {
    s: u32,
    i: u32,
    d: u32,
    aux: u32,
}

// Traceback bit positions (LSB)
const TRACE_BIT_POS: u32 = 0;
const ALIGN_LEN_BIT_POS: u32 = 4;
const ALIGN_LEN_BIT_MASK: u32 = 0b0111_1111_1111_1111_1111_1111_1111;
const CONTIG_IDX_POS: u32 = 31;
const AUX_CONTIG_IDX_MASK: u32 = 0b1_1111;
const AUX_CONTIG_FROM_POS: u32 = 5;

impl PackedLengthCell {
    /// Sets 4 bits [pos, pos+4) with the 4 LSBs of value
    #[inline(always)]
    fn set_tb(&mut self, matrix: u32, tb: u16) -> u32 {
        let bits: u32 = (0b1111) << TRACE_BIT_POS;
        assert!(
            tb <= TB_MAX,
            "Expected a value <= TB_MAX while setting traceback bits"
        );
        (matrix & !bits) | ((tb as u32) << TRACE_BIT_POS) // And set the bits
    }

    #[inline(always)]
    fn set_len(&mut self, matrix: u32, len: u32) -> u32 {
        let bits: u32 = ALIGN_LEN_BIT_MASK << ALIGN_LEN_BIT_POS;
        (matrix & !bits) | (len << ALIGN_LEN_BIT_POS) // And set the bits
    }

    #[inline(always)]
    fn set_idx(&mut self, idx: u32) {
        // upper bit of idx goes to upper bit of self.s
        self.s = (self.s & !(1 << CONTIG_IDX_POS)) | (((idx >> 7) & 0b1) << CONTIG_IDX_POS);
        // second-most upper bit of idx goes to upper bit of self.i
        self.i = (self.i & !(1 << CONTIG_IDX_POS)) | (((idx >> 6) & 0b1) << CONTIG_IDX_POS);
        // third-most upper bit of idx goes to upper bit of self.d
        self.d = (self.d & !(1 << CONTIG_IDX_POS)) | (((idx >> 5) & 0b1) << CONTIG_IDX_POS);
        // all but the three-most upper bits go to self.aux
        self.aux = (self.aux & !AUX_CONTIG_IDX_MASK) | (idx & AUX_CONTIG_IDX_MASK);
    }

    #[inline(always)]
    fn set_from(&mut self, from: u32) {
        self.aux = (self.aux & AUX_CONTIG_IDX_MASK) | (from << AUX_CONTIG_FROM_POS)
        // And set the bits
    }

    // Gets 4 bits [pos, pos+4) of v
    #[inline(always)]
    fn get_tb(self, matrix: u32) -> u16 {
        ((matrix >> TRACE_BIT_POS) & (0b1111)) as u16
    }

    #[inline(always)]
    fn get_len(self, matrix: u32) -> u32 {
        (matrix >> ALIGN_LEN_BIT_POS) & ALIGN_LEN_BIT_MASK
    }

    #[inline(always)]
    fn get_idx(self) -> u32 {
        let mut value = 0;
        // upper bit of idx comes from upper bit of self.s
        value |= (self.s >> 31) << 7;
        // second-most upper bit of idx comes from upper bit of self.i
        value |= (self.i >> 31) << 6;
        // third-most upper bit of idx comes from upper bit of self.d
        value |= (self.d >> 31) << 5;
        // all but the three-most upper bits comes from self.aux
        value |= self.aux & AUX_CONTIG_IDX_MASK;
        value
    }

    #[inline(always)]
    fn get_from(self) -> u32 {
        self.aux >> AUX_CONTIG_FROM_POS
    }
}

impl TracebackCell for PackedLengthCell {
    fn max_target_len() -> u32 {
        134_217_727
    }

    fn max_num_contigs() -> u32 {
        255
    }

    #[inline(always)]
    fn set_i(&mut self, tb: u16, len: u32) {
        // Traceback corresponding to matrix I
        self.i = self.set_tb(self.i, tb);
        self.i = self.set_len(self.i, len);
    }

    #[inline(always)]
    fn set_d(&mut self, tb: u16, len: u32) {
        // Traceback corresponding to matrix D
        self.d = self.set_tb(self.d, tb);
        self.d = self.set_len(self.d, len);
    }

    #[inline(always)]
    fn set_s(&mut self, tb: u16, len: u32) {
        // Traceback corresponding to matrix S
        self.s = self.set_tb(self.s, tb);
        self.s = self.set_len(self.s, len);
    }

    #[inline(always)]
    fn set_s_all(&mut self, tb: u16, len: u32, idx: u32, from: u32) {
        assert!(idx <= Self::max_num_contigs());
        assert!(from <= Self::max_target_len());
        // Traceback corresponding to matrix S
        self.s = self.set_tb(self.s, tb);
        self.s = self.set_len(self.s, len);
        self.set_idx(idx);
        self.set_from(from);
    }

    #[inline(always)]
    fn get_i_len(self) -> u32 {
        self.get_len(self.i)
    }

    #[inline(always)]
    fn get_i(self) -> (u16, u32) {
        (self.get_tb(self.i), self.get_len(self.i))
    }

    #[inline(always)]
    fn get_d_len(self) -> u32 {
        self.get_len(self.d)
    }

    #[inline(always)]
    fn get_d(self) -> (u16, u32) {
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
            idx: self.get_idx(),
            from: self.get_from(),
        }
    }
}

// Tests
#[cfg(test)]
pub mod tests {
    use rstest::rstest;

    use crate::align::traceback::{
        packed_length_cell::PackedLengthCell, SValue, TracebackCell, TB_MAX,
    };

    #[rstest]
    fn test_set_and_get_i() {
        let mut cell = PackedLengthCell::default();
        for tb in 0..=TB_MAX {
            cell.set_i(0, 0);
            assert_eq!(cell.get_i(), (0, 0));
            assert_eq!(cell.get_i_len(), 0);
            cell.set_i(tb, 13);
            assert_eq!(cell.get_i(), (tb, 13));
            assert_eq!(cell.get_i_len(), 13);
            cell.set_i(0, 0);
        }
    }

    #[rstest]
    fn test_set_and_get_d() {
        let mut cell = PackedLengthCell::default();
        for tb in 0..=TB_MAX {
            cell.set_d(0, 0);
            assert_eq!(cell.get_d(), (0, 0));
            assert_eq!(cell.get_d_len(), 0);
            cell.set_d(tb, 13);
            assert_eq!(cell.get_d(), (tb, 13));
            assert_eq!(cell.get_d_len(), 13);
            cell.set_d(0, 0);
        }
    }

    #[rstest]
    fn test_set_and_get_s() {
        let mut cell = PackedLengthCell::default();
        for tb in 0..=TB_MAX {
            cell.set_s_all(0, 0, 0, 0);
            assert_eq!(
                cell.get_s(),
                SValue {
                    tb: 0,
                    len: 0,
                    idx: 0,
                    from: 0
                }
            );
            assert_eq!(cell.get_s_len(), 0);
            cell.set_s(tb, 13);
            assert_eq!(
                cell.get_s(),
                SValue {
                    tb,
                    len: 13,
                    idx: 0,
                    from: 0
                }
            );
            assert_eq!(cell.get_s_len(), 13);
            cell.set_s_all(tb, 14, 8, 22);
            assert_eq!(
                cell.get_s(),
                SValue {
                    tb,
                    len: 14,
                    idx: 8,
                    from: 22
                }
            );
            assert_eq!(cell.get_s_len(), 14);
        }
    }
}
