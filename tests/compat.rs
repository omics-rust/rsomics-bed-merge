use std::path::Path;
use std::process::Command;

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
fn bedtools_compat() {
    let bedtools = Command::new("bedtools").arg("--version").output();
    if bedtools.is_err() || !bedtools.unwrap().status.success() {
        eprintln!("bedtools not available — skipping compat test");
        return;
    }

    let input = golden("sorted.bed");
    let mut ours = Vec::new();
    merge(&input, &mut ours).unwrap();
    let ours_str = String::from_utf8(ours).unwrap();

    let bt = Command::new("bedtools")
        .args(["merge", "-i"])
        .arg(&input)
        .output()
        .expect("bedtools merge failed");
    let bt_str = String::from_utf8(bt.stdout).unwrap();

    let mut ours_lines: Vec<&str> = ours_str.lines().filter(|l| !l.is_empty()).collect();
    let mut bt_lines: Vec<&str> = bt_str.lines().filter(|l| !l.is_empty()).collect();
    ours_lines.sort_unstable();
    bt_lines.sort_unstable();

    assert_eq!(ours_lines, bt_lines, "output differs from bedtools merge");
}
