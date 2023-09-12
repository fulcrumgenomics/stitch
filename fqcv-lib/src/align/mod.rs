mod aligners;
mod alignment;
pub mod io;
mod scoring;
mod sub_alignment;
mod traceback;

pub use aligners::{AlignmentMode, Builder};

use anyhow::{anyhow, Error};
use std::{fmt::Display, str::FromStr};

/// The various strategies to pick the primary alignment amongst multiple sub-alignments.
#[derive(Default, Debug, PartialEq, Eq, Copy, Clone)]
pub enum PrimaryPickingStrategy {
    #[default]
    QueryLength,
    Score,
}

impl Display for PrimaryPickingStrategy {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::QueryLength => write!(f, "query-length"),
            Self::Score => write!(f, "score"),
        }
    }
}

impl FromStr for PrimaryPickingStrategy {
    type Err = Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "query-length" | "query_length" | "querylength" | "length" => {
                Ok(PrimaryPickingStrategy::QueryLength)
            }
            "score" => Ok(PrimaryPickingStrategy::Score),
            _ => Err(anyhow!("Invalid primary picking strategy: {}", s)),
        }
    }
}
