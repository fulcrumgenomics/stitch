use clap::ValueEnum;

pub(crate) mod aligners;
mod alignment;
pub(crate) mod io;
pub(crate) mod scoring;
pub(crate) mod sub_alignment;
mod traceback;

/// The various strategies to pick the primary alignment amonst multiple sub-alignments.
#[derive(Default, Debug, PartialEq, Eq, Copy, Clone, ValueEnum)]
pub enum PrimaryPickingStrategy {
    #[default]
    QueryLength,
    Score,
}
