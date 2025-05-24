use RustConstantTermSequences::dfao::DFAO;
use RustConstantTermSequences::laurent_poly::LaurentPoly;
use RustConstantTermSequences::mod_int::ModInt;
use RustConstantTermSequences::sequences::{constant_term, constant_term_slow};
use RustConstantTermSequences::*;
use criterion::{BenchmarkId, Criterion, black_box, criterion_group, criterion_main};
use std::time::Duration;

pub fn rowland_zeilberger(c: &mut Criterion) {
    let p = 29;
    let P = LaurentPoly::from_string("x + 1 + x^-1", p);
    let Q = LaurentPoly::from_string("1 - x^2", p);
    let mut group = c.benchmark_group("Rowland-Zeilberger");
    for n in 1u64..10 {
        group.bench_with_input(BenchmarkId::new("By Reduction", n), &n, |b, n| {
            b.iter(|| constant_term(&P, &Q, &(1000 * n)));
        });
        group.bench_with_input(BenchmarkId::new("Direct", n), &n, |b, n| {
            b.iter(|| constant_term_slow(&P, &Q, &(10 * n)));
        });
    }
    group.finish();
}

pub fn first_zero_motzkin(c: &mut Criterion) {
    let mut group = c.benchmark_group("first_zero_motzkin");
    for p in [2, 11, 17, 23, 29].iter() {
        group.bench_with_input(BenchmarkId::new("DFAO-based", p), p, |b, p| {
            b.iter(|| {
                let P = LaurentPoly::from_string("x + 1 + x^-1", *p);
                let Q = LaurentPoly::from_string("1 - x^2", *p);
                let _ = DFAO::compute_shortest_zero(&P, &Q, 10000);
            });
        });
        group.bench_with_input(BenchmarkId::new("Direct", p), p, |b, p| {
            b.iter(|| {
                let P = LaurentPoly::from_string("x + 1 + x^-1", *p);
                let Q = LaurentPoly::from_string("1 - x^2", *p);
                let mut n = 0;
                while constant_term(&P, &Q, &n) != ModInt::zero(*p) {
                    n += 1;
                }
            });
        });
    }
    group.finish();
}

pub fn first_zero_big(c: &mut Criterion) {
    let mut group = c.benchmark_group("first_zero_big");
    let Q = LaurentPoly::one(2);
    for P in [
        LaurentPoly::from_string("x^2 + 1 + x^-2", 2),
        LaurentPoly::from_string("x^5 + x^3 + x^2 + 1 + x^-2 + x^-3 + x^-5", 2),
        LaurentPoly::from_string("x^7 + x^4 + x^2 + x + 1 + x^-1 + x^-2 + x^-4 + x^-7", 2),
        LaurentPoly::from_string("x^8 + x^3 + x + 1 + x^-3 + x^-8", 2),
        LaurentPoly::from_string(
            "x^9 + x^8 + x^6 + x^4 + x^2 + x + 1 + x^-1 + x^-2 + x^-4 + x^-6 + x^-8 + x^-9",
            2,
        ),
    ]
    .iter()
    {
        group.bench_with_input(BenchmarkId::new("DFAO-based", P.degree()), P, |b, P| {
            b.iter(|| {
                let _ = DFAO::compute_shortest_zero(&P, &Q, 10000);
            });
        });
        group.bench_with_input(BenchmarkId::new("Direct", P.degree()), P, |b, P| {
            b.iter(|| {
                let mut n = 0;
                while constant_term(&P, &Q, &n) != ModInt::zero(2) && n < 2^(2^(2*P.degree() + 1)) {
                    n += 1;
                }
            });
        });
    }
    group.finish();
}

criterion_group! {
    name = benches;
    config = Criterion::default().sample_size(20).warm_up_time(Duration::from_secs(1));
    targets = rowland_zeilberger, first_zero_motzkin, first_zero_big
}
criterion_main!(benches);
