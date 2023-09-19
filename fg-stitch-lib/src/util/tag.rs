//! Custom SAM tags for stitch.  See the README.md file for the description of each tag.

/// (`qs`)
pub const QUERY_START: &str = "qs";

/// (`qe`)
pub const QUERY_END: &str = "qe";

/// (`ts`)
pub const TARGET_START: &str = "ts";

/// (`te`)
pub const TARGET_END: &str = "te";

/// (`as`)
pub const CHAIN_ALIGNMENT_SCORE: &str = "as";

/// (`xs`)
pub const SUBOPTIMAL_SCORE: &str = "xs";

/// (`si`)
pub const SUB_ALIGNMENT_INDEX: &str = "si";

/// (`sc`)
pub const SUB_ALIGNMENT_CIGAR: &str = "sc";

/// (`cl`)
pub const CHAIN_LENGTH: &str = "cl";

/// (`ci`)
pub const CHAIN_INDEX: &str = "ci";

/// (`cn`)
pub const NUMBER_OF_CHAINS: &str = "cn";
