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
use std::path::Path;

/// Returns true if the path ends with a recognized GZIP file extension
fn is_path_with_extension<P: AsRef<Path>>(p: &P, extensions: [&str; 2]) -> bool {
    if let Some(ext) = p.as_ref().extension() {
        match ext.to_str() {
            Some(x) => extensions.contains(&x),
            None => false,
        }
    } else {
        false
    }
}

/// The set of file extensions to treat as GZIPPED
const GZIP_EXTENSIONS: [&str; 2] = ["gz", "bgz"];

/// Returns true if the path ends with a recognized GZIP file extension
pub fn is_gzip_path<P: AsRef<Path>>(p: &P) -> bool {
    is_path_with_extension(p, GZIP_EXTENSIONS)
}

/// The set of file extensions to treat as FASTQ
const FASTQ_EXTENSIONS: [&str; 2] = ["fastq", "fq"];

/// Returns true if the path ends with a recognized FASTQ file extension
pub fn is_fastq_path<P: AsRef<Path>>(p: &P) -> bool {
    is_path_with_extension(p, FASTQ_EXTENSIONS)
}
