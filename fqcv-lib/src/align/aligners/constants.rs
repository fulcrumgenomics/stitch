use anyhow::{anyhow, Error};
use std::{fmt::Display, str::FromStr};

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
    Match,               // Consumes one x and one y base
    Subst,               // Consumes one x and one y base
    Del,                 // Consumes a single x base
    Ins,                 // Consumes a single y base
    Xclip(usize),        // Consumes N x bases at the start or end of the x
    Yclip(usize),        // Consumes N x bases at the start or end of the y
    Xjump(usize, usize), // Consumes N x bases (contig_idx, from_idx)
    Yjump(usize),        // Consumes N y bases (from_idx)
}

impl AlignmentOperation {
    pub fn is_special(&self) -> bool {
        use crate::align::aligners::constants::AlignmentOperation::{Xclip, Xjump, Yclip};
        matches!(self, Xclip(_) | Yclip(_) | Xjump(_, _))
    }

    pub fn as_string(&self, contig_idx: usize, x_index: usize) -> String {
        match *self {
            AlignmentOperation::Match => "=".to_string(),
            AlignmentOperation::Subst => "X".to_string(),
            AlignmentOperation::Del => "D".to_string(),
            AlignmentOperation::Ins => "I".to_string(),
            AlignmentOperation::Xclip(l) => format!("{l}A"),
            AlignmentOperation::Yclip(l) => format!("{l}B"),
            AlignmentOperation::Xjump(new_contig_idx, new_x_index) => {
                let contig_jump_str = match new_contig_idx.cmp(&contig_idx) {
                    std::cmp::Ordering::Greater => format!("{}C", new_contig_idx - contig_idx),
                    std::cmp::Ordering::Less => format!("{}c", contig_idx - new_contig_idx),
                    std::cmp::Ordering::Equal => String::new(),
                };
                if new_x_index >= x_index {
                    format!("{contig_jump_str}{}J", new_x_index - x_index)
                } else {
                    format!("{contig_jump_str}{}j", x_index - new_x_index)
                }
            }
            AlignmentOperation::Yjump(y_jump_len) => format!("{}S", y_jump_len),
        }
    }

    pub fn length_on_x(&self, x_index: usize) -> i32 {
        use crate::align::aligners::constants::AlignmentOperation::{
            Del, Ins, Match, Subst, Xclip, Xjump, Yclip, Yjump,
        };
        match *self {
            Match | Subst | Ins => 1,
            Del | Yclip(_) | Yjump(_) => 0,
            Xclip(len) => len as i32,
            Xjump(_, to_x_index) => to_x_index as i32 - x_index as i32,
        }
    }

    #[allow(dead_code)]
    pub fn length_on_y(&self) -> usize {
        use crate::align::aligners::constants::AlignmentOperation::{
            Del, Ins, Match, Subst, Xclip, Xjump, Yclip, Yjump,
        };
        match *self {
            Match | Subst | Del => 1,
            Yclip(len) => len,
            Yjump(len) => len,
            Ins | Xclip(_) | Xjump(_, _) => 0,
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
#[derive(Default, Debug, PartialEq, Eq, Copy, Clone)]
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

impl Display for AlignmentMode {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Local => write!(f, "local"),
            Self::QueryLocal => write!(f, "query-local"),
            Self::TargetLocal => write!(f, "target-local"),
            Self::Global => write!(f, "global"),
            Self::Custom => write!(f, "custom"),
        }
    }
}

impl FromStr for AlignmentMode {
    type Err = Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_ascii_lowercase().as_str() {
            "local" => Ok(AlignmentMode::Local),
            "query-local" | "query_local" | "querylocal" | "query" => Ok(AlignmentMode::QueryLocal),
            "target-local" | "target_local" | "targetlocal" | "target" => {
                Ok(AlignmentMode::TargetLocal)
            }
            "global" => Ok(AlignmentMode::Global),
            "custom" => Ok(AlignmentMode::Custom),
            _ => Err(anyhow!("Invalid alignment mode: {}", s)),
        }
    }
}
