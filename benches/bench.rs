use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;
use std::path::PathBuf;
use std::process::Command;

fn bench_bed_merge(c: &mut Criterion) {
    let bin = env!("CARGO_BIN_EXE_rsomics-bed-merge");
    let manifest = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let bed = manifest.join("tests/golden/sorted.bed");
    c.bench_function("rsomics-bed-merge golden", |b| {
        b.iter(|| {
            let out = Command::new(black_box(bin))
                .arg(bed.to_str().unwrap())
                .output()
                .unwrap();
            assert!(out.status.success());
        });
    });
}

criterion_group!(benches, bench_bed_merge);
criterion_main!(benches);
