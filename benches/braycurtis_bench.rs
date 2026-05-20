use braycurtis::BrayCurtisPlugin;
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use pluma_plugin_trait::PluMAPlugin;
use std::hint::black_box;

fn generate_matrix(n_samples: usize, n_otus: usize) -> (Vec<String>, Vec<String>, Vec<Vec<f64>>) {
    let samples: Vec<String> = (0..n_samples).map(|i| format!("S{i}")).collect();
    let otus: Vec<String> = (0..n_otus).map(|j| format!("OTU{j}")).collect();
    let counts: Vec<Vec<f64>> = (0..n_samples)
        .map(|i| {
            (0..n_otus)
                .map(|j| {
                    let v = ((i * 31 + j * 17) % 100) as f64;
                    if (i + j) % 5 == 0 { 0.0 } else { v }
                })
                .collect()
        })
        .collect();
    (samples, otus, counts)
}

fn bench_braycurtis(c: &mut Criterion) {
    let mut group = c.benchmark_group("braycurtis");
    for (n_samples, n_otus) in [(25, 500), (100, 1000), (250, 2000), (500, 2000)] {
        let (samples, otus, counts) = generate_matrix(n_samples, n_otus);
        let plugin = BrayCurtisPlugin::from_matrix(samples, otus, counts);
        group.bench_with_input(
            BenchmarkId::new("compute", format!("{n_samples}x{n_otus}")),
            &plugin,
            |b, p| {
                b.iter(|| {
                    let mut p = p.clone();
                    p.run().unwrap();
                    black_box(p.dissimilarity().len())
                })
            },
        );
    }
    group.finish();
}

criterion_group!(benches, bench_braycurtis);
criterion_main!(benches);
