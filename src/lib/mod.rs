#![deny(unsafe_code)]
#![allow(
    clippy::must_use_candidate,
    clippy::missing_panics_doc,
    clippy::missing_errors_doc,
    clippy::module_name_repetitions,
    clippy::similar_names,
    clippy::too_many_lines,
    clippy::too_many_arguments,
    clippy::struct_excessive_bools
)]

mod align;
pub mod commands;
pub mod util;
