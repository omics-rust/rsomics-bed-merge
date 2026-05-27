//! Merge overlapping BED intervals — `bedtools merge` equivalent.
//!
//! Requires the input to be sorted by chromosome then start. Adjacent or
//! overlapping intervals on the same chromosome are collapsed into a single
//! merged interval. The implementation wraps `rsomics_intervals::bed::merge_sorted`.

use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};
use rsomics_intervals::bed;

/// Merge sorted BED from `input` to `output`.
pub fn merge(input: &Path, output: &mut dyn Write) -> Result<()> {
    let w = BufWriter::new(output);
    let f = File::open(input).map_err(RsomicsError::Io)?;
    bed::merge_sorted(f, w)
}

/// Merge sorted BED from stdin to `output`.
pub fn merge_stdin(output: &mut dyn Write) -> Result<()> {
    let w = BufWriter::new(output);
    bed::merge_sorted(io::stdin().lock(), w)
}
