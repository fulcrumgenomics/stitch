use crate::alignment::constants::AlignmentMode;
use crate::util::built_info;
use clap::Parser;
use env_logger::Env;
use std::path::PathBuf;

pub static TOOL_NAME: &str = "fqvqc";

/// Aligns stuff
#[derive(Parser, Debug, Clone)]
#[clap(name = TOOL_NAME, version = built_info::VERSION.as_str(), term_width=0)]
pub struct Opts {
    /// The path to the input FASTQ with sequenced vector/plasmid/construct long reads.
    #[clap(long, short = 'f', display_order = 1)]
    pub reads_fastq: PathBuf,

    /// The path to the referece vector/plasmid/construct FASTA sequence.
    #[clap(long, short = 'r', display_order = 2)]
    pub ref_fasta: PathBuf,

    /// The path to the referece vector/plasmid/construct FASTA sequence.
    #[clap(long, short = 'd', default_value = "false", display_order = 2)]
    pub double_strand: bool,

    /// The number of threads to use.
    #[clap(long, short = 't', default_value = "2", display_order = 3)]
    pub threads: usize,

    /// Assume an unrecognized input (based on file extension) is GZIP compressed.
    #[clap(long, short = 'z', default_value = "false", display_order = 4)]
    pub decompress: bool,

    /// Pre-align with banded local alignment.
    #[clap(long, short = 'p', default_value = "false", display_order = 5)]
    pub pre_align: bool,

    /// K-mer size for banded pre-alignment.
    #[clap(long, short = 'k', default_value = "12", display_order = 6)]
    pub k: usize,

    /// Band size for banded pre-alignment.
    #[clap(long, short = 'w', default_value = "50", display_order = 7)]
    pub w: usize,

    /// Pre-align with banded local alignment.
    #[clap(long, short = 's', default_value = "100", display_order = 8)]
    pub pre_align_min_score: i32,

    /// Use soft-clipping for all alignments, otherwise secondary alignemnts will use hard-clipping
    #[clap(long, short = 'S', default_value = "false", display_order = 9)]
    pub soft_clip: bool,

    /// Use =/X CIGAR operators, otherwise use M
    #[clap(long, short = 'X', default_value = "false", display_order = 10)]
    pub use_eq_and_x: bool,

    /// Score for a sequence match (must be positive)
    #[clap(long, short = 'A', default_value = "1", display_order = 11)]
    pub match_score: i32,

    /// Score for a sequence mismatch (must be negative)
    #[clap(long, short = 'B', default_value = "-1", display_order = 12)]
    pub mismatch_score: i32,

    /// Score for a gap open (must be negative)
    #[clap(long, short = 'O', default_value = "-5", display_order = 13)]
    pub gap_open: i32,

    /// Score for a gap extend (must be negative); a gap of size k cost '{-O} + {-E}*k'
    #[clap(long, short = 'E', default_value = "-1", display_order = 14)]
    pub gap_extend: i32,

    /// Score for a target jump (must be negative)
    #[clap(long, short = 'J', default_value = "-10", display_order = 15)]
    pub jump_score: i32,

    /// FIXME: not wired in
    /// The alignment mode
    #[clap(long, short = 'm', default_value = "local", display_order = 16)]
    pub mode: AlignmentMode,

    /// The compression level of the output BAM
    #[clap(long, short = 'c', default_value = "0", display_order = 17)]
    pub compression: u8,
}

/// Parse args and set up logging / tracing
pub fn setup() -> Opts {
    if std::env::var("RUST_LOG").is_err() {
        std::env::set_var("RUST_LOG", "info");
    }
    env_logger::Builder::from_env(Env::default().default_filter_or("info")).init();

    Opts::parse()
}
