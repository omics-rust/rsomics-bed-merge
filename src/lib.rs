//! Merge overlapping BED intervals — `bedtools merge` equivalent.
//!
//! Input must be sorted by chromosome then start.
//!
//! A zero-length feature (`start == end`, a half-open insertion point) is
//! handled the way `bedtools merge` does: for overlap testing and for the
//! reported span it is virtually widened to `[start-1, start+1]`. So a
//! zero-length feature next to another interval merges with that virtual
//! footprint (`chr1 250 260` + `chr1 260 260` → `chr1 250 261`), while a
//! zero-length feature that merges nothing is emitted at its original
//! coordinates (`chr1 10 10` → `chr1 10 10`). Widening the low edge of a
//! feature at coordinate 0 yields -1, which bedtools emits verbatim, so span
//! bounds are signed.

use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Read, Write};
use std::path::Path;

use rsomics_common::{Context, Result, RsomicsError};
use rsomics_intervals::IntervalError;

pub fn merge(input: &Path, output: &mut dyn Write) -> Result<()> {
    let f = File::open(input).map_err(RsomicsError::Io)?;
    merge_reader(f, output)
}

pub fn merge_stdin(output: &mut dyn Write) -> Result<()> {
    merge_reader(io::stdin().lock(), output)
}

/// Sweep-merge a pre-sorted (chrom, start) BED3 stream. Out-of-order starts or
/// a chromosome that reappears after closing fail loud; chrom is opaque bytes,
/// so a non-UTF8 name passes through untouched.
fn merge_reader<R: Read, W: Write>(r: R, w: W) -> Result<()> {
    let mut rdr = BufReader::new(r);
    let mut w = BufWriter::new(w);
    let mut line: Vec<u8> = Vec::with_capacity(256);
    let mut chrom: Vec<u8> = Vec::with_capacity(32);
    let mut closed: Vec<Vec<u8>> = Vec::new();
    let mut have = false;
    // cluster start is fixed by the first feature's widened low edge; the end is
    // the running max of widened high edges. `cfirst` is the first feature's raw
    // start, used only for the sorted-order check.
    let (mut cfirst, mut clo, mut chi) = (0_u64, 0_i128, 0_i128);
    let mut csingle_zero: Option<(u64, u64)> = None;
    let mut lineno = 0_usize;

    loop {
        line.clear();
        if rdr.read_until(b'\n', &mut line).map_err(RsomicsError::Io)? == 0 {
            break;
        }
        lineno += 1;
        while matches!(line.last(), Some(b'\n' | b'\r')) {
            line.pop();
        }
        if line.is_empty()
            || line[0] == b'#'
            || line.starts_with(b"track")
            || line.starts_with(b"browser")
        {
            continue;
        }
        let (lc, ls, le) = parse_bed3_bytes(&line)
            .map_err(|e| RsomicsError::InvalidInput(format!("BED line {lineno}: {e}")))?;
        let (lo, hi) = widened(ls, le);

        if have && lc == chrom.as_slice() {
            if ls < cfirst {
                return Err(RsomicsError::InvalidInput(format!(
                    "BED line {lineno}: input not sorted (start {ls} < {cfirst} on same chrom) — \
                     sort with `sort -k1,1 -k2,2n` first"
                )));
            }
            if lo <= chi {
                chi = chi.max(hi);
                csingle_zero = None;
                continue;
            }
            emit(&mut w, &chrom, clo, chi, csingle_zero)?;
        } else {
            if have {
                emit(&mut w, &chrom, clo, chi, csingle_zero)?;
                closed.push(chrom.clone());
            }
            if closed.iter().any(|c| c.as_slice() == lc) {
                return Err(RsomicsError::InvalidInput(format!(
                    "BED line {lineno}: input not sorted (chromosome {} reappears after it \
                     closed) — sort with `sort -k1,1 -k2,2n` first",
                    String::from_utf8_lossy(lc)
                )));
            }
            chrom.clear();
            chrom.extend_from_slice(lc);
        }
        cfirst = ls;
        clo = lo;
        chi = hi;
        csingle_zero = (ls == le).then_some((ls, le));
        have = true;
    }
    if have {
        emit(&mut w, &chrom, clo, chi, csingle_zero)?;
    }
    w.flush().map_err(RsomicsError::Io)?;
    Ok(())
}

/// Virtual footprint used for overlap and span: a zero-length feature widens to
/// `[start-1, start+1]`; every other feature keeps its own bounds.
fn widened(start: u64, end: u64) -> (i128, i128) {
    if start == end {
        (i128::from(start) - 1, i128::from(start) + 1)
    } else {
        (i128::from(start), i128::from(end))
    }
}

/// A cluster of a single zero-length feature reports its original coordinates,
/// not the widened footprint; every other cluster reports the widened span.
fn emit<W: Write>(
    w: &mut W,
    chrom: &[u8],
    lo: i128,
    hi: i128,
    single_zero: Option<(u64, u64)>,
) -> Result<()> {
    w.write_all(chrom).rs_context("writing merged BED")?;
    match single_zero {
        Some((start, end)) => writeln!(w, "\t{start}\t{end}"),
        None => writeln!(w, "\t{lo}\t{hi}"),
    }
    .rs_context("writing merged BED")?;
    Ok(())
}

fn parse_bed3_bytes(s: &[u8]) -> std::result::Result<(&[u8], u64, u64), String> {
    let mut it = s.split(|&c| c == b'\t');
    let chrom = it.next().ok_or("missing chrom")?;
    let start = parse_u64(it.next().ok_or("missing start")?)?;
    let end = parse_u64(it.next().ok_or("missing end")?)?;
    if start > end {
        return Err(IntervalError::Inverted { start, end }.to_string());
    }
    Ok((chrom, start, end))
}

fn parse_u64(b: &[u8]) -> std::result::Result<u64, String> {
    if b.is_empty() {
        return Err("empty integer field".into());
    }
    let mut n: u64 = 0;
    for &c in b {
        let d = c.wrapping_sub(b'0');
        if d > 9 {
            return Err(format!("bad integer {:?}", String::from_utf8_lossy(b)));
        }
        n = n
            .checked_mul(10)
            .and_then(|n| n.checked_add(u64::from(d)))
            .ok_or_else(|| format!("integer overflows u64: {:?}", String::from_utf8_lossy(b)))?;
    }
    Ok(n)
}
