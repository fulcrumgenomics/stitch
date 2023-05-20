use clap::ValueEnum;

/// Value to use as a 'negative infinity' score. Should be close to `i32::MIN`,
/// but avoid underflow when used with reasonable scoring parameters or even
/// adding two negative infinities. Use ~ `0.4 * i32::MIN`
pub const MIN_SCORE: i32 = -858_993_459;

/// Alignment operations supported are match, substitution, insertion, deletion
/// and clipping. Clipping is a special boundary condition where you are allowed
/// to clip off the beginning/end of the sequence for a fixed clip penalty. The
/// clip penalty could be different for the two sequences x and y, and the
/// clipping operations on both are distinguishable (Xclip and Yclip). The usize
/// value associated with the clipping operations are the lengths clipped. In case
/// of standard modes like Global, Semi-Global and Local alignment, the clip operations
/// are filtered out.
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(Eq, PartialEq, Debug, Copy, Clone, Hash)]
pub enum AlignmentOperation {
    Match,        // Consumes one query and one target base
    Subst,        // Consumes one query and one target base
    Del,          // Consumes a single target base
    Ins,          // Consumes a single query base
    Xclip(usize), // Consumes N query bases at the start or end of the query
    Yclip(usize), // Consumes N query bases at the start or end of the target
    Xskip(usize), // Consumes N query bases (can be negative)
    Xflip(usize), // Consumes N query bases (can be negative) and switches strand
}

/// The modes of alignment supported by the aligner include standard modes such as
/// Global, Semi-Global and Local alignment. In addition to this, user can also invoke
/// the custom mode. In the custom mode, users can explicitly specify the clipping penalties
/// for prefix and suffix of strings 'x' and 'y' independently. Under the hood the standard
/// modes are implemented as special cases of the custom mode with the clipping penalties
/// appropriately set.
///
/// The default alignment mode is Global.
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(Default, Debug, PartialEq, Eq, Copy, Clone, ValueEnum)]
pub enum AlignmentMode {
    #[default]
    Local,
    Semiglobal,
    Global,
    Custom,
}
