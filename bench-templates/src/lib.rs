use ark_crypto_primitives::{
    crh::{sha256::digest::Digest, CRHScheme},
    sponge::{
        poseidon::{PoseidonConfig, PoseidonSponge},
        CryptographicSponge,
    },
};
use ark_ff::PrimeField;
use ark_poly::Polynomial;
use ark_serialize::{CanonicalSerialize, Compress};
use ark_std::{test_rng, UniformRand};
use rand_chacha::{
    rand_core::{RngCore, SeedableRng},
    ChaCha20Rng,
};

use core::time::Duration;
use std::{borrow::Borrow, marker::PhantomData, time::Instant};

use ark_poly_commit::{to_bytes, LabeledPolynomial, PolynomialCommitment};

pub use criterion::*;
pub use paste::paste;

/// Measure the time cost of `method` (i.e., commit/open/verify) of a
/// multilinear PCS for all `num_vars` specified in `nv_list`.
/// `rand_poly` is a function that outputs a random multilinear polynomial.
/// `rand_point` is a function that outputs a random point in the domain of polynomial.
pub fn bench_pcs_method<F: PrimeField, P: Polynomial<F>, PCS: PolynomialCommitment<F, P>>(
    c: &mut Criterion,
    nv_list: Vec<usize>,
    msg: &str,
    method: impl Fn(
        &PCS::CommitterKey,
        &PCS::VerifierKey,
        usize,
        fn(usize, &mut ChaCha20Rng) -> P,
        fn(usize, &mut ChaCha20Rng) -> P::Point,
    ) -> Duration,
    rand_poly: fn(usize, &mut ChaCha20Rng) -> P,
    rand_point: fn(usize, &mut ChaCha20Rng) -> P::Point,
) {
    let mut group = c.benchmark_group(msg);
    let rng = &mut ChaCha20Rng::from_rng(test_rng()).unwrap();

    for num_vars in nv_list {
        let pp = PCS::setup(num_vars, Some(num_vars), rng).unwrap();
        let (ck, vk) = PCS::trim(&pp, num_vars, num_vars, None).unwrap();

        group.bench_with_input(
            BenchmarkId::from_parameter(num_vars),
            &num_vars,
            |b, num_vars| {
                b.iter_custom(|i| {
                    let mut time = Duration::from_nanos(0);
                    for _ in 0..i {
                        time += method(&ck, &vk, *num_vars, rand_poly, rand_point);
                    }
                    time
                });
            },
        );
    }

    group.finish();
}

/// Report the time cost of a commitment
pub fn commit<F: PrimeField, P: Polynomial<F>, PCS: PolynomialCommitment<F, P>>(
    ck: &PCS::CommitterKey,
    _vk: &PCS::VerifierKey,
    num_vars: usize,
    rand_poly: fn(usize, &mut ChaCha20Rng) -> P,
    _rand_point: fn(usize, &mut ChaCha20Rng) -> P::Point,
) -> Duration {
    let rng = &mut ChaCha20Rng::from_rng(test_rng()).unwrap();

    let labeled_poly =
        LabeledPolynomial::new("test".to_string(), rand_poly(num_vars, rng), None, None);

    let start = Instant::now();
    let (_, _) = PCS::commit(&ck, [&labeled_poly], Some(rng)).unwrap();
    start.elapsed()
}

/// Report the size of a commitment
pub fn commitment_size<F: PrimeField, P: Polynomial<F>, PCS: PolynomialCommitment<F, P>>(
    num_vars: usize,
    rand_poly: fn(usize, &mut ChaCha20Rng) -> P,
) -> usize {
    let rng = &mut ChaCha20Rng::from_rng(test_rng()).unwrap();

    let pp = PCS::setup(num_vars, Some(num_vars), rng).unwrap();

    let (ck, _) = PCS::trim(&pp, num_vars, num_vars, None).unwrap();

    let labeled_poly =
        LabeledPolynomial::new("test".to_string(), rand_poly(num_vars, rng), None, None);

    let (coms, _) = PCS::commit(&ck, [&labeled_poly], Some(rng)).unwrap();

    coms[0].commitment().serialized_size(Compress::No)
}

/// Report the time cost of an opening
pub fn open<F, P, PCS>(
    ck: &PCS::CommitterKey,
    _vk: &PCS::VerifierKey,
    num_vars: usize,
    rand_poly: fn(usize, &mut ChaCha20Rng) -> P,
    rand_point: fn(usize, &mut ChaCha20Rng) -> P::Point,
) -> Duration
where
    F: PrimeField,
    P: Polynomial<F>,
    PCS: PolynomialCommitment<F, P>,
{
    let rng = &mut ChaCha20Rng::from_rng(test_rng()).unwrap();

    let labeled_poly =
        LabeledPolynomial::new("test".to_string(), rand_poly(num_vars, rng), None, None);

    let (coms, states) = PCS::commit(&ck, [&labeled_poly], Some(rng)).unwrap();
    let point = rand_point(num_vars, rng);

    let start = Instant::now();
    let _ = PCS::open(
        &ck,
        [&labeled_poly],
        &coms,
        &point,
        &mut test_sponge::<F>(),
        &states,
        Some(rng),
    )
    .unwrap();
    start.elapsed()
}

/// Report the size of a proof
pub fn proof_size<F, P, PCS>(num_vars: usize, rand_poly: fn(usize, &mut ChaCha20Rng) -> P) -> usize
where
    F: PrimeField,
    P: Polynomial<F>,
    PCS: PolynomialCommitment<F, P>,
    P::Point: UniformRand,
{
    let rng = &mut ChaCha20Rng::from_rng(test_rng()).unwrap();

    let pp = PCS::setup(num_vars, Some(num_vars), rng).unwrap();

    let (ck, _) = PCS::trim(&pp, num_vars, num_vars, None).unwrap();
    let labeled_poly =
        LabeledPolynomial::new("test".to_string(), rand_poly(num_vars, rng), None, None);

    let (coms, states) = PCS::commit(&ck, [&labeled_poly], Some(rng)).unwrap();
    let point = P::Point::rand(rng);

    let proofs = PCS::open(
        &ck,
        [&labeled_poly],
        &coms,
        &point,
        &mut test_sponge::<F>(),
        &states,
        Some(rng),
    )
    .unwrap();

    let bproof: PCS::BatchProof = vec![proofs].into();

    bproof.serialized_size(Compress::No)
}

/// Report the time cost of a verification
pub fn verify<F, P, PCS>(
    ck: &PCS::CommitterKey,
    vk: &PCS::VerifierKey,
    num_vars: usize,
    rand_poly: fn(usize, &mut ChaCha20Rng) -> P,
    rand_point: fn(usize, &mut ChaCha20Rng) -> P::Point,
) -> Duration
where
    F: PrimeField,
    P: Polynomial<F>,
    PCS: PolynomialCommitment<F, P>,
{
    let rng = &mut ChaCha20Rng::from_rng(test_rng()).unwrap();

    let labeled_poly =
        LabeledPolynomial::new("test".to_string(), rand_poly(num_vars, rng), None, None);

    let (coms, states) = PCS::commit(&ck, [&labeled_poly], Some(rng)).unwrap();
    let point = rand_point(num_vars, rng);
    let claimed_eval = labeled_poly.evaluate(&point);
    let proof = PCS::open(
        &ck,
        [&labeled_poly],
        &coms,
        &point,
        &mut test_sponge::<F>(),
        &states,
        Some(rng),
    )
    .unwrap();

    let start = Instant::now();
    PCS::check(
        &vk,
        &coms,
        &point,
        [claimed_eval],
        &proof,
        &mut test_sponge::<F>(),
        None,
    )
    .unwrap();
    start.elapsed()
}

/*************** Auxiliary functions ***************/

fn test_sponge<F: PrimeField>() -> PoseidonSponge<F> {
    let full_rounds = 8;
    let partial_rounds = 31;
    let alpha = 17;

    let mds = vec![
        vec![F::one(), F::zero(), F::one()],
        vec![F::one(), F::one(), F::zero()],
        vec![F::zero(), F::one(), F::one()],
    ];

    let mut v = Vec::new();
    let mut ark_rng = test_rng();

    for _ in 0..(full_rounds + partial_rounds) {
        let mut res = Vec::new();

        for _ in 0..3 {
            res.push(F::rand(&mut ark_rng));
        }
        v.push(res);
    }
    let config = PoseidonConfig::new(full_rounds, partial_rounds, alpha, mds, v, 2, 1);
    PoseidonSponge::new(&config)
}

#[macro_export]
macro_rules! bench_method {
    ($c:expr, $method:ident, $scheme_type:ty, $rand_poly:ident, $rand_point:ident) => {
        let scheme_type_str = stringify!($scheme_type);
        let bench_name = format!("{} {}", stringify!($method), scheme_type_str);
        bench_pcs_method::<_, _, $scheme_type>(
            $c,
            (MIN_NUM_VARS..MAX_NUM_VARS).step_by(2).collect(),
            &bench_name,
            $method::<_, _, $scheme_type>,
            $rand_poly::<_>,
            $rand_point::<_>,
        );
    };
}

#[macro_export]
macro_rules! bench {
    (
        $scheme_type:ty, $rand_poly:ident, $rand_point:ident
    ) => {
        fn bench_pcs(c: &mut Criterion) {
            bench_method!(c, commit, $scheme_type, $rand_poly, $rand_point);
            bench_method!(c, open, $scheme_type, $rand_poly, $rand_point);
            bench_method!(c, verify, $scheme_type, $rand_poly, $rand_point);
        }

        criterion_group!(benches, bench_pcs);

        paste! {
            criterion_main!(
                benches
            );
        }
    };
}

/**** Auxiliary methods for linear-code-based PCSs ****/

/// Needed for benches and tests.
pub struct LeafIdentityHasher;

impl CRHScheme for LeafIdentityHasher {
    type Input = Vec<u8>;
    type Output = Vec<u8>;
    type Parameters = ();

    fn setup<R: RngCore>(_: &mut R) -> Result<Self::Parameters, ark_crypto_primitives::Error> {
        Ok(())
    }

    fn evaluate<T: Borrow<Self::Input>>(
        _: &Self::Parameters,
        input: T,
    ) -> Result<Self::Output, ark_crypto_primitives::Error> {
        Ok(input.borrow().to_vec().into())
    }
}

/// Needed for benches and tests.
pub struct FieldToBytesColHasher<F, D>
where
    F: PrimeField + CanonicalSerialize,
    D: Digest,
{
    _phantom: PhantomData<(F, D)>,
}

impl<F, D> CRHScheme for FieldToBytesColHasher<F, D>
where
    F: PrimeField + CanonicalSerialize,
    D: Digest,
{
    type Input = Vec<F>;
    type Output = Vec<u8>;
    type Parameters = ();

    fn setup<R: RngCore>(_rng: &mut R) -> Result<Self::Parameters, ark_crypto_primitives::Error> {
        Ok(())
    }

    fn evaluate<T: Borrow<Self::Input>>(
        _parameters: &Self::Parameters,
        input: T,
    ) -> Result<Self::Output, ark_crypto_primitives::Error> {
        let mut dig = D::new();
        dig.update(to_bytes!(input.borrow()).unwrap());
        Ok(dig.finalize().to_vec())
    }
}
