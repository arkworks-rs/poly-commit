use crate::hyrax::HyraxPC;
use crate::tests::*;
use crate::utils::test_sponge;
use crate::{LabeledPolynomial, PolynomialCommitment};
use ark_bls12_377::G1Affine;
use ark_crypto_primitives::sponge::poseidon::PoseidonSponge;
use ark_ec::AffineRepr;
use ark_ed_on_bls12_381::EdwardsAffine;
use ark_ff::PrimeField;
use ark_poly::{DenseMultilinearExtension, MultilinearExtension};
use ark_std::test_rng;
use rand_chacha::{rand_core::SeedableRng, ChaCha20Rng};

// The test structure is largely taken from the multilinear_ligero module
// inside this crate

// ****************** types ******************

type Fq = <G1Affine as AffineRepr>::ScalarField;
type Hyrax377 = HyraxPC<G1Affine, DenseMultilinearExtension<Fq>, PoseidonSponge<Fq>>;

type Fr = <EdwardsAffine as AffineRepr>::ScalarField;
type Hyrax381 = HyraxPC<EdwardsAffine, DenseMultilinearExtension<Fr>, PoseidonSponge<Fr>>;

// ******** auxiliary test functions ********

fn rand_poly<Fr: PrimeField>(
    _: usize, // degree: unused
    num_vars: Option<usize>,
    rng: &mut ChaCha20Rng,
) -> DenseMultilinearExtension<Fr> {
    match num_vars {
        Some(n) => DenseMultilinearExtension::rand(n, rng),
        None => panic!("Must specify the number of variables"),
    }
}

fn constant_poly<Fr: PrimeField>(
    _: usize, // degree: unused
    num_vars: Option<usize>,
    rng: &mut ChaCha20Rng,
) -> DenseMultilinearExtension<Fr> {
    match num_vars {
        Some(0) => DenseMultilinearExtension::rand(0, rng),
        _ => panic!("Must specify the number of variables: 0"),
    }
}

fn rand_point<F: PrimeField>(num_vars: Option<usize>, rng: &mut ChaCha20Rng) -> Vec<F> {
    match num_vars {
        Some(n) => (0..n).map(|_| F::rand(rng)).collect(),
        None => panic!("Must specify the number of variables"),
    }
}

// ****************** tests ******************

#[test]
fn test_hyrax_construction() {
    // Desired number of variables (must be even!)
    let n = 8;

    let chacha = &mut ChaCha20Rng::from_rng(test_rng()).unwrap();

    let pp = Hyrax381::setup(1, Some(n), chacha).unwrap();

    let (ck, vk) = Hyrax381::trim(&pp, 1, 1, None).unwrap();

    let l_poly = LabeledPolynomial::new(
        "test_poly".to_string(),
        rand_poly::<Fr>(0, Some(n), chacha),
        None,
        None,
    );

    let (c, rands) = Hyrax381::commit(&ck, &[l_poly.clone()], Some(chacha)).unwrap();

    let point: Vec<Fr> = rand_point(Some(n), chacha);
    let value = l_poly.evaluate(&point);

    // Dummy argument
    let mut test_sponge = test_sponge::<Fr>();

    let proof = Hyrax381::open(
        &ck,
        &[l_poly],
        &c,
        &point,
        &mut (test_sponge.clone()),
        &rands,
        Some(chacha),
    )
    .unwrap();

    assert!(Hyrax381::check(
        &vk,
        &c,
        &point,
        [value],
        &proof,
        &mut test_sponge,
        Some(chacha),
    )
    .unwrap());
}

#[test]
fn hyrax_single_poly_test() {
    single_poly_test::<_, _, Hyrax377, _>(
        Some(10),
        rand_poly,
        rand_point,
        poseidon_sponge_for_test,
    )
    .expect("test failed for bls12-377");
    single_poly_test::<_, _, Hyrax381, _>(
        Some(10),
        rand_poly,
        rand_point,
        poseidon_sponge_for_test,
    )
    .expect("test failed for bls12-381");
}

#[test]
fn hyrax_constant_poly_test() {
    single_poly_test::<_, _, Hyrax377, _>(
        Some(0),
        constant_poly,
        rand_point,
        poseidon_sponge_for_test,
    )
    .expect("test failed for bls12-377");
    single_poly_test::<_, _, Hyrax381, _>(
        Some(0),
        constant_poly,
        rand_point,
        poseidon_sponge_for_test,
    )
    .expect("test failed for bls12-381");
}

#[test]
fn hyrax_full_end_to_end_test() {
    full_end_to_end_test::<_, _, Hyrax377, _>(
        Some(8),
        rand_poly,
        rand_point,
        poseidon_sponge_for_test,
    )
    .expect("test failed for bls12-377");
    full_end_to_end_test::<_, _, Hyrax381, _>(
        Some(10),
        rand_poly,
        rand_point,
        poseidon_sponge_for_test,
    )
    .expect("test failed for bls12-381");
}

#[test]
fn hyrax_single_equation_test() {
    single_equation_test::<_, _, Hyrax377, _>(
        Some(6),
        rand_poly,
        rand_point,
        poseidon_sponge_for_test,
    )
    .expect("test failed for bls12-377");
    single_equation_test::<_, _, Hyrax381, _>(
        Some(6),
        rand_poly,
        rand_point,
        poseidon_sponge_for_test,
    )
    .expect("test failed for bls12-381");
}

#[test]
fn hyrax_two_equation_test() {
    two_equation_test::<_, _, Hyrax377, _>(
        Some(10),
        rand_poly,
        rand_point,
        poseidon_sponge_for_test,
    )
    .expect("test failed for bls12-377");
    two_equation_test::<_, _, Hyrax381, _>(
        Some(10),
        rand_poly,
        rand_point,
        poseidon_sponge_for_test,
    )
    .expect("test failed for bls12-381");
}

#[test]
fn hyrax_full_end_to_end_equation_test() {
    full_end_to_end_equation_test::<_, _, Hyrax377, _>(
        Some(8),
        rand_poly,
        rand_point,
        poseidon_sponge_for_test,
    )
    .expect("test failed for bls12-377");
    full_end_to_end_equation_test::<_, _, Hyrax381, _>(
        Some(8),
        rand_poly,
        rand_point,
        poseidon_sponge_for_test,
    )
    .expect("test failed for bls12-381");
}
