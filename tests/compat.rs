use std::io::Write;
use std::path::Path;
use std::process::{Command, Stdio};

use rsomics_bed_merge::merge;

fn golden(name: &str) -> std::path::PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

#[test]
fn basic_merge_correctness() {
    let input = golden("sorted.bed");
    let mut out = Vec::new();
    merge(&input, &mut out).unwrap();
    let result = String::from_utf8(out).unwrap();
    let lines: Vec<&str> = result.lines().filter(|l| !l.is_empty()).collect();
    // chr1: [100,200) + [150,300) merge to [100,300); [400,500) separate → 2 lines
    // chr2: [10,50) + [45,100) merge to [10,100); [200,300) separate → 2 lines
    assert_eq!(lines.len(), 4, "expected 4 merged intervals: {result}");
    assert!(
        result.contains("chr1\t100\t300"),
        "chr1 merge wrong: {result}"
    );
    assert!(
        result.contains("chr1\t400\t500"),
        "chr1 separate wrong: {result}"
    );
    assert!(
        result.contains("chr2\t10\t100"),
        "chr2 merge wrong: {result}"
    );
    assert!(
        result.contains("chr2\t200\t300"),
        "chr2 separate wrong: {result}"
    );
}

#[test]
fn golden_matches_committed_upstream() {
    let input = golden("sorted.bed");
    let expected = std::fs::read_to_string(golden("merge.upstream.expected")).unwrap();

    let mut ours = Vec::new();
    merge(&input, &mut ours).unwrap();
    let ours_str = String::from_utf8(ours).unwrap();

    assert_eq!(
        ours_str, expected,
        "output differs from bedtools merge golden"
    );
}

// Zero-length features (start == end) are virtually widened to [start-1, start+1]
// for overlap and span, exactly as bedtools merge does: a standalone one passes
// through verbatim, an adjacent one drags the merged span out by one, and a
// widened low edge at coordinate 0 reaches -1.
#[test]
fn zero_length_golden_matches_committed_upstream() {
    let input = golden("zero_length.bed");
    let expected = std::fs::read_to_string(golden("zero_length.upstream.expected")).unwrap();

    let mut ours = Vec::new();
    merge(&input, &mut ours).unwrap();
    let ours_str = String::from_utf8(ours).unwrap();

    assert_eq!(
        ours_str, expected,
        "zero-length output differs from bedtools merge golden"
    );
}

#[test]
fn zero_length_key_cases() {
    let mut ours = Vec::new();
    merge(&golden("zero_length.bed"), &mut ours).unwrap();
    let out = String::from_utf8(ours).unwrap();
    // standalone zero-length is emitted verbatim
    assert!(out.contains("chr1\t10\t10\n"), "{out}");
    assert!(out.contains("chr4\t40\t40\n"), "{out}");
    // adjacent zero-length widens the merged span by one
    assert!(out.contains("chr1\t250\t261\n"), "{out}");
    // two coincident zero-length features merge to their shared widened footprint
    assert!(out.contains("chr4\t9\t11\n"), "{out}");
    // widening the low edge at coordinate 0 reaches -1
    assert!(out.contains("chr1\t-1\t5\n"), "{out}");
}

#[test]
fn inverted_interval_fails_loud() {
    let bin = env!("CARGO_BIN_EXE_rsomics-bed-merge");
    let mut child = Command::new(bin)
        .stdin(Stdio::piped())
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .spawn()
        .unwrap();
    child
        .stdin
        .take()
        .unwrap()
        .write_all(b"chr1\t200\t100\n")
        .unwrap();
    let status = child.wait().unwrap();
    assert!(
        !status.success(),
        "start > end must exit non-zero, got {status}"
    );
}

#[test]
fn bedtools_compat() {
    let bedtools = Command::new("bedtools").arg("--version").output();
    if bedtools.is_err() || !bedtools.unwrap().status.success() {
        eprintln!("bedtools not available — skipping compat test");
        return;
    }

    for name in ["sorted.bed", "zero_length.bed"] {
        let input = golden(name);
        let mut ours = Vec::new();
        merge(&input, &mut ours).unwrap();
        let ours_str = String::from_utf8(ours).unwrap();

        let bt = Command::new("bedtools")
            .args(["merge", "-i"])
            .arg(&input)
            .output()
            .expect("bedtools merge failed");
        let bt_str = String::from_utf8(bt.stdout).unwrap();

        assert_eq!(
            ours_str, bt_str,
            "output differs from bedtools merge on {name}"
        );
    }
}
