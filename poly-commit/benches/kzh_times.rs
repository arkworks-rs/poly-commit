//! End-to-end KZH setup, commitment, opening, and verification benchmarks.
//!
//! The complete KZH-2/3/4/6 × n=12/16/20 matrix is intentionally expensive.
//! Set `KZH_BENCH_K` and/or `KZH_BENCH_N` to benchmark one family member or
//! evaluation-table size while tuning an implementation.

use ark_bls12_381::{Bls12_381, Fr};
use ark_crypto_primitives::sponge::{
    poseidon::{PoseidonConfig, PoseidonSponge},
    CryptographicSponge,
};
use ark_ff::{One, Zero};
use ark_pcs_bench_templates::*;
use ark_poly::{DenseMultilinearExtension, MultilinearExtension};
use ark_poly_commit::{
    kzh::{PreparedVerifierKey, KZH},
    LabeledPolynomial, PolynomialCommitment,
};
use ark_serialize::{CanonicalSerialize, Compress};
use ark_std::UniformRand;
use rand_chacha::{rand_core::SeedableRng, ChaCha20Rng};
use std::time::Duration;

type Mle = DenseMultilinearExtension<Fr>;
type Kzh<const K: usize> = KZH<Bls12_381, Mle, K>;

const NUM_VARS: [usize; 3] = [12, 16, 20];

fn matches_requested_value(requested: Option<usize>, value: usize) -> bool {
    match requested {
        Some(requested) => requested == value,
        None => true,
    }
}

/// Returns deterministic benchmark-only sponge parameters.
fn benchmark_sponge_config() -> PoseidonConfig<Fr> {
    let full_rounds = 8;
    let partial_rounds = 31;
    let alpha = 17;
    let mds = vec![
        vec![Fr::one(), Fr::zero(), Fr::one()],
        vec![Fr::one(), Fr::one(), Fr::zero()],
        vec![Fr::zero(), Fr::one(), Fr::one()],
    ];

    let mut rng = ChaCha20Rng::from_seed([91u8; 32]);
    let ark = (0..full_rounds + partial_rounds)
        .map(|_| (0..3).map(|_| Fr::rand(&mut rng)).collect())
        .collect();
    PoseidonConfig::new(full_rounds, partial_rounds, alpha, mds, ark, 2, 1)
}

fn bench_family<const K: usize>(criterion: &mut Criterion, num_vars: usize) {
    let num_evaluations = 1usize << num_vars;
    let seed = (K as u8).wrapping_mul(31).wrapping_add(num_vars as u8);
    let mut fixture_rng = ChaCha20Rng::from_seed([seed; 32]);

    // KZH setup is tied to one exact number of variables. It is non-hiding, so
    // both trim and commitment use a zero hiding bound and no commitment RNG.
    let params = Kzh::<K>::setup(1, Some(num_vars), &mut fixture_rng).unwrap();
    let params_size = params.serialized_size(Compress::Yes);
    let block_sizes = params.num_vars_per_block().to_vec();
    let (committer_key, verifier_key) = Kzh::<K>::trim(&params, 1, 0, None).unwrap();
    // At n=20 the SRS contains roughly one million G1 points. `trim` clones
    // the committer material, so release the universal copy before measuring.
    drop(params);
    let polynomial = LabeledPolynomial::new(
        format!("kzh-{K}-n-{num_vars}"),
        Mle::rand(num_vars, &mut fixture_rng),
        None,
        None,
    );
    let point = (0..num_vars)
        .map(|_| Fr::rand(&mut fixture_rng))
        .collect::<Vec<_>>();
    let value = polynomial.evaluate(&point);
    let (commitments, states) = Kzh::<K>::commit(&committer_key, [&polynomial], None).unwrap();

    let sponge_config = benchmark_sponge_config();
    let mut opening_sponge = PoseidonSponge::new(&sponge_config);
    let proof = Kzh::<K>::open(
        &committer_key,
        [&polynomial],
        &commitments,
        &point,
        &mut opening_sponge,
        &states,
        None,
    )
    .unwrap();
    let mut checking_sponge = PoseidonSponge::new(&sponge_config);
    assert!(Kzh::<K>::check(
        &verifier_key,
        &commitments,
        &point,
        [value],
        &proof,
        &mut checking_sponge,
        None,
    )
    .unwrap());
    let prepared_verifier_key = PreparedVerifierKey::prepare(&verifier_key);
    let mut prepared_checking_sponge = PoseidonSponge::new(&sponge_config);
    assert!(Kzh::<K>::check_prepared(
        &prepared_verifier_key,
        &commitments,
        &point,
        [value],
        &proof,
        &mut prepared_checking_sponge,
        None,
    )
    .unwrap());

    println!(
        "KZH-{K}, n={num_vars}, N={num_evaluations}, blocks={block_sizes:?}: params={} B, ck={} B, vk={} B, commitment={} B, state={} B, proof={} B (compressed)",
        params_size,
        committer_key.serialized_size(Compress::Yes),
        verifier_key.serialized_size(Compress::Yes),
        commitments[0]
            .commitment()
            .serialized_size(Compress::Yes),
        states[0].serialized_size(Compress::Yes),
        proof.serialized_size(Compress::Yes),
    );

    let mut group = criterion.benchmark_group(format!("KZH-{K}/n={num_vars}"));
    if num_vars == 20 {
        group
            .sample_size(10)
            .warm_up_time(Duration::from_millis(100))
            .measurement_time(Duration::from_secs(1))
            .sampling_mode(SamplingMode::Flat);
    }

    // Setup is useful at the smaller sizes, but ten repeated million-point
    // SRS generations obscure the commit/open/check scaling of interest.
    if num_vars < 20 {
        group.bench_function("setup", |bencher| {
            bencher.iter_batched(
                || ChaCha20Rng::from_seed([seed; 32]),
                |mut rng| {
                    black_box(Kzh::<K>::setup(1, Some(num_vars), &mut rng).unwrap());
                },
                BatchSize::SmallInput,
            );
        });
    }

    group.bench_function("commit", |bencher| {
        bencher.iter(|| {
            black_box(Kzh::<K>::commit(&committer_key, [&polynomial], None).unwrap());
        });
    });

    group.bench_function("open", |bencher| {
        bencher.iter_batched(
            || PoseidonSponge::new(&sponge_config),
            |mut sponge| {
                black_box(
                    Kzh::<K>::open(
                        &committer_key,
                        [&polynomial],
                        &commitments,
                        &point,
                        &mut sponge,
                        &states,
                        None,
                    )
                    .unwrap(),
                );
            },
            BatchSize::SmallInput,
        );
    });

    group.bench_function("check", |bencher| {
        bencher.iter_batched(
            || PoseidonSponge::new(&sponge_config),
            |mut sponge| {
                assert!(black_box(
                    Kzh::<K>::check(
                        &verifier_key,
                        &commitments,
                        &point,
                        [value],
                        &proof,
                        &mut sponge,
                        None,
                    )
                    .unwrap()
                ));
            },
            BatchSize::SmallInput,
        );
    });

    // Preparation is intentionally measured separately: applications that
    // reuse a verifier key should prepare it once and benchmark the amortized
    // `check_prepared` path independently from this one-time cost.
    group.bench_function("prepare_vk", |bencher| {
        bencher.iter(|| {
            black_box(PreparedVerifierKey::prepare(black_box(&verifier_key)));
        });
    });

    group.bench_function("check_prepared", |bencher| {
        bencher.iter_batched(
            || PoseidonSponge::new(&sponge_config),
            |mut sponge| {
                assert!(black_box(
                    Kzh::<K>::check_prepared(
                        &prepared_verifier_key,
                        &commitments,
                        &point,
                        [value],
                        &proof,
                        &mut sponge,
                        None,
                    )
                    .unwrap()
                ));
            },
            BatchSize::SmallInput,
        );
    });

    group.bench_function("prepare_and_check", |bencher| {
        bencher.iter_batched(
            || PoseidonSponge::new(&sponge_config),
            |mut sponge| {
                let one_shot_key = PreparedVerifierKey::prepare(&verifier_key);
                assert!(black_box(
                    Kzh::<K>::check_prepared(
                        &one_shot_key,
                        &commitments,
                        &point,
                        [value],
                        &proof,
                        &mut sponge,
                        None,
                    )
                    .unwrap()
                ));
            },
            BatchSize::SmallInput,
        );
    });

    group.finish();
}

fn bench_kzh(criterion: &mut Criterion) {
    let requested_k = std::env::var("KZH_BENCH_K").ok().map(|value| {
        value
            .parse::<usize>()
            .expect("KZH_BENCH_K must be an integer")
    });
    let requested_num_vars = std::env::var("KZH_BENCH_N").ok().map(|value| {
        value
            .parse::<usize>()
            .expect("KZH_BENCH_N must be an integer")
    });
    if let Some(k) = requested_k {
        assert!(
            [2, 3, 4, 6].contains(&k),
            "KZH_BENCH_K must be 2, 3, 4, or 6"
        );
    }
    if let Some(num_vars) = requested_num_vars {
        assert!(
            NUM_VARS.contains(&num_vars),
            "KZH_BENCH_N must be one of {:?}",
            NUM_VARS
        );
    }

    // K=3 is an odd arity; K=2, 4, and 6 are even. n increases by four so that
    // each step multiplies the evaluation table by 16.
    for num_vars in NUM_VARS
        .iter()
        .copied()
        .filter(|num_vars| matches_requested_value(requested_num_vars, *num_vars))
    {
        if matches_requested_value(requested_k, 2) {
            bench_family::<2>(criterion, num_vars);
        }
        if matches_requested_value(requested_k, 3) {
            bench_family::<3>(criterion, num_vars);
        }
        if matches_requested_value(requested_k, 4) {
            bench_family::<4>(criterion, num_vars);
        }
        if matches_requested_value(requested_k, 6) {
            bench_family::<6>(criterion, num_vars);
        }
    }
}

criterion_group! {
    name = benches;
    config = Criterion::default()
        .sample_size(10)
        .warm_up_time(Duration::from_secs(1))
        .measurement_time(Duration::from_secs(3));
    targets = bench_kzh
}
criterion_main!(benches);
