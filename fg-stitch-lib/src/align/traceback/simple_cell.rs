use crate::align::traceback::TB_MAX;

use super::{SValue, TracebackCell};
use serde::{Deserialize, Serialize};

/// Packed representation of one cell of a Smith-Waterman traceback matrix.  Does not store
/// alignment length for each matrix value, and so cannot be used to break ties to prefer longer
/// alignments.
///
/// The maximum target (x) length is 2^28 - 1 = 268,435,456bp.  The maximum supported number of
/// contigs is 2^8 - 1 = 255 contigs.
///
/// The `tb` field (traceback) is packed as follows:
/// - bits 0-3 are for insertion traceback (I_POS)
/// - bits 4-7 are for deletion traceback (D_POS)
/// - bits 8-11 are for S matrix traceback (S_POS)
/// - bit 12-15 are for the upper 4-bits of the jump contig index (i.e. which contig did we jump from)
/// The `aux` field is packed as follows:
/// - bits 0-3 are for the lower 4-bits of the jump contig index (i.e. which contig did we jump from)
/// - bits 4-31 are for the "from" jump index in the contig (i.e. where in the given contig did we jump from)
#[derive(
    Default, Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize,
)]
pub struct SimpleCell {
    tb: u16,
    aux: u32,
}

// Traceback bit positions (LSB)
const TB_I_POS: u16 = 0; // Meaning bits 0,1,2,3 corresponds to I and so on
const TB_D_POS: u16 = 4;
const TB_S_POS: u16 = 8;
const TB_IDX_POS: u16 = 12;
const TB_MASK: u16 = 0b1111;
const AUX_IDX_MASK: u32 = 0b1111;
const AUX_FROM_MASK: u32 = 0b1111_1111_1111_1111_1111_1111_1111;

impl SimpleCell {
    /// Sets 4 bits [pos, pos+4) with the 4 LSBs of value
    #[inline(always)]
    fn set_tb(&mut self, pos: u16, value: u16) {
        assert!(
            pos == TB_IDX_POS || value <= TB_MAX,
            "Expected a value <= TB_MAX while setting traceback bits"
        );
        let bits: u16 = (TB_MASK) << pos;
        self.tb = (self.tb & !bits) // First clear the bits
              | (value << pos) // And set the bits
    }

    // Gets 4 bits [pos, pos+4) of v
    #[inline(always)]
    fn get_tb(self, pos: u16) -> u16 {
        (self.tb >> pos) & (0b1111)
    }
}

impl TracebackCell for SimpleCell {
    fn max_target_len() -> u32 {
        268_435_456
    }

    fn max_num_contigs() -> u32 {
        255
    }

    #[inline(always)]
    fn set_i(&mut self, tb: u16, _len: u32) {
        // Traceback corresponding to matrix I
        self.set_tb(TB_I_POS, tb);
    }

    #[inline(always)]
    fn set_d(&mut self, tb: u16, _len: u32) {
        // Traceback corresponding to matrix D
        self.set_tb(TB_D_POS, tb);
    }

    #[inline(always)]
    fn set_s(&mut self, tb: u16, _len: u32) {
        // Traceback corresponding to matrix S
        self.set_tb(TB_S_POS, tb);
    }

    /// Sets the S traceback with contig index and from information.
    #[inline(always)]
    fn set_s_all(&mut self, tb: u16, len: u32, idx: u32, from: u32) {
        assert!(idx <= Self::max_num_contigs());
        assert!(from <= Self::max_target_len());
        // Traceback corresponding to matrix S
        self.set_s(tb, len);
        // contig index (upper 4 bits, then lower 4bits)
        self.set_tb(TB_IDX_POS, (idx >> 4) as u16);
        self.aux = (self.aux & !AUX_IDX_MASK) | (idx & AUX_IDX_MASK);
        // contig from
        self.aux = (self.aux & !(AUX_FROM_MASK << 4)) | (from << 4); // contig
    }

    #[inline(always)]
    fn get_i(self) -> (u16, u32) {
        (self.get_tb(TB_I_POS), 0)
    }

    #[inline(always)]
    fn get_d(self) -> (u16, u32) {
        (self.get_tb(TB_D_POS), 0)
    }

    #[inline(always)]
    fn get_s(self) -> SValue {
        let idx: u32 = (self.get_tb(TB_IDX_POS) << 4) as u32 | self.aux & AUX_IDX_MASK;
        let from: u32 = self.aux >> 4;

        SValue {
            tb: self.get_tb(TB_S_POS),
            len: 0, // alignment length is never stored, so always zero
            idx,
            from,
        }
    }

    #[inline(always)]
    fn get_i_len(self) -> u32 {
        // alignment length is never stored, so always zero
        0
    }

    #[inline(always)]
    fn get_d_len(self) -> u32 {
        // alignment length is never stored, so always zero
        0
    }

    #[inline(always)]
    fn get_s_len(self) -> u32 {
        // alignment length is never stored, so always zero
        0
    }
}

// Tests
#[cfg(test)]
pub mod tests {
    use rstest::rstest;

    use crate::align::traceback::{simple_cell::SimpleCell, SValue, TracebackCell, TB_MAX};

    #[rstest]
    fn test_set_and_get_i() {
        let mut cell = SimpleCell::default();
        for tb in 0..=TB_MAX {
            assert_eq!(cell.get_i(), (0, 0));
            assert_eq!(cell.get_i_len(), 0);
            cell.set_i(tb, 13);
            assert_eq!(cell.get_i(), (tb, 0));
            assert_eq!(cell.get_i_len(), 0);
            cell.set_i(0, 0);
        }
    }

    #[rstest]
    fn test_set_and_get_d() {
        let mut cell = SimpleCell::default();
        for tb in 0..=TB_MAX {
            assert_eq!(cell.get_d(), (0, 0));
            assert_eq!(cell.get_d_len(), 0);
            cell.set_d(tb, 13);
            assert_eq!(cell.get_d(), (tb, 0));
            assert_eq!(cell.get_d_len(), 0);
            cell.set_d(0, 0);
        }
    }

    #[rstest]
    fn test_set_and_get_s() {
        let mut cell = SimpleCell::default();
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
                    len: 0,
                    idx: 0,
                    from: 0
                }
            );
            assert_eq!(cell.get_s_len(), 0);
            cell.set_s_all(tb, 14, 8, 22);
            assert_eq!(
                cell.get_s(),
                SValue {
                    tb,
                    len: 0,
                    idx: 8,
                    from: 22
                }
            );
            assert_eq!(cell.get_s_len(), 0);
        }
    }
}
