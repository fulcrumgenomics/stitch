use clap::ValueEnum;

/// Value to use as a 'negative infinity' score. Should be close to `i32::MIN`,
/// but avoid underflow when used with reasonable scoring parameters or even
/// adding two negative infinities. Use ~ `0.4 * i32::MIN`
pub const MIN_SCORE: i32 = -858_993_459;

pub const DEFAULT_ALIGNER_CAPACITY: usize = 200;

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
    Xskip(usize), // Consumes N query bases
    Xflip(usize), // Consumes N query bases and switches strand
}

impl AlignmentOperation {
    pub fn is_special(&self) -> bool {
        use crate::align::aligners::constants::AlignmentOperation::{Xclip, Xflip, Xskip, Yclip};
        matches!(self, Xclip(_) | Yclip(_) | Xskip(_) | Xflip(_))
    }

    pub fn as_string(&self, x_index: usize) -> String {
        match *self {
            AlignmentOperation::Match => "=".to_string(),
            AlignmentOperation::Subst => "X".to_string(),
            AlignmentOperation::Del => "D".to_string(),
            AlignmentOperation::Ins => "I".to_string(),
            AlignmentOperation::Xclip(l) => format!("{}A", l),
            AlignmentOperation::Yclip(l) => format!("{}B", l),
            AlignmentOperation::Xskip(new_x_index) => {
                if new_x_index > x_index {
                    format!("{}J", new_x_index - x_index)
                } else {
                    format!("{}j", x_index - new_x_index)
                }
            }
            AlignmentOperation::Xflip(new_x_index) => {
                if new_x_index > x_index {
                    format!("{}F", new_x_index - x_index)
                } else {
                    format!("{}f", x_index - new_x_index)
                }
            }
        }
    }

    pub fn length_on_x(&self, x_index: i32) -> i32 {
        use crate::align::aligners::constants::AlignmentOperation::{
            Del, Ins, Match, Subst, Xclip, Xflip, Xskip, Yclip,
        };
        match *self {
            Match | Subst | Ins => 1,
            Del | Yclip(_) => 0,
            Xclip(len) => len as i32,
            Xskip(to_x_index) | Xflip(to_x_index) => to_x_index as i32 - x_index,
        }
    }

    #[allow(dead_code)]
    pub fn length_on_y(&self) -> i32 {
        use crate::align::aligners::constants::AlignmentOperation::{
            Del, Ins, Match, Subst, Xclip, Xflip, Xskip, Yclip,
        };
        match *self {
            Match | Subst | Del => 1,
            Yclip(len) => len as i32,
            Ins | Xclip(_) | Xskip(_) | Xflip(_) => 0,
        }
    }
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
    /// Aligns a sub-sequence of the read versus a sub-sequence of the reference
    #[default]
    Local,
    /// Aligns a sub-sequence of the read versus the full reference.
    QueryLocal,
    /// Aligns the full read versus a sub-sequence of the reference
    TargetLocal,
    /// Aligns the full read versus the full reference.
    Global,
    Custom,
}
