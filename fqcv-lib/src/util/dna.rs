use lazy_static::lazy_static;
use std::borrow::Borrow;

/// Valid IUPAC DNA bases
pub const IUPAC_BASES: [u8; 15] = *b"AGCTYRWSKMDVHBN";

/// The complement of IUPAC DNA bases.
pub const IUPAC_BASES_COMPLEMENT: [u8; 15] = *b"TCGARYWSMKHBDVN";

lazy_static! {
    /// An array-based look up of the DNA complement for each IUPAC bases
    pub static ref COMPLEMENT: [u8; 256] = {
        let mut comp = [0; 256];
        for (v, a) in comp.iter_mut().enumerate() {
            *a = v as u8;
        }
        for (&a, &b) in IUPAC_BASES.iter().zip(IUPAC_BASES_COMPLEMENT.iter()) {
            comp[a as usize] = b;
            comp[a as usize + 32] = b + 32;  // lowercase variants
        }
        comp
    };
}

/// Complements a given DNA IUPAC base.
fn complement(a: u8) -> u8 {
    COMPLEMENT[a as usize]
}

/// Reverse complements a DNA IUPAC base.
pub fn reverse_complement<C, T>(text: T) -> Vec<u8>
where
    C: Borrow<u8>,
    T: IntoIterator<Item = C>,
    T::IntoIter: DoubleEndedIterator,
{
    text.into_iter()
        .rev()
        .map(|a| complement(*a.borrow()))
        .collect()
}
