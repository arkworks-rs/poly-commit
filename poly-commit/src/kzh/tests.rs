#[cfg(not(feature = "std"))]
use ark_std::{string::ToString, vec::Vec};

use crate::{tests::poseidon_sponge_for_test, LabeledPolynomial, PolynomialCommitment};
use ark_bls12_381::{Bls12_381, Fr};
use ark_ff::One;
use ark_poly::{DenseMultilinearExtension, MultilinearExtension, SparseMultilinearExtension};
use ark_std::{test_rng, UniformRand};
use rand_chacha::{rand_core::SeedableRng, ChaCha20Rng};

use super::KZH;

type DenseKZH<const K: usize> = KZH<Bls12_381, DenseMultilinearExtension<Fr>, K>;

fn commit_open_check<const K: usize, P: MultilinearExtension<Fr>>(
    num_vars: usize,
    polynomial: &LabeledPolynomial<Fr, P>,
    rng: &mut ChaCha20Rng,
) {
    let params = KZH::<Bls12_381, P, K>::setup(1, Some(num_vars), rng).unwrap();
    let (ck, vk) = KZH::<Bls12_381, P, K>::trim(&params, 1, 0, None).unwrap();
    let point: Vec<Fr> = (0..num_vars).map(|_| Fr::rand(rng)).collect();
    let value = polynomial.evaluate(&point);
    let (commitments, states) = KZH::<Bls12_381, P, K>::commit(&ck, [polynomial], None).unwrap();
    let base_sponge = poseidon_sponge_for_test::<Fr>();
    let proof = KZH::<Bls12_381, P, K>::open(
        &ck,
        [polynomial],
        &commitments,
        &point,
        &mut base_sponge.clone(),
        &states,
        None,
    )
    .unwrap();
    assert!(KZH::<Bls12_381, P, K>::check(
        &vk,
        &commitments,
        &point,
        [value],
        &proof,
        &mut base_sponge.clone(),
        None,
    )
    .unwrap());
}

#[test]
fn commit_open_check_supports_dense_and_sparse_polynomials() {
    let mut rng = ChaCha20Rng::from_rng(test_rng()).unwrap();
    const NUM_VARS: usize = 8;

    let dense = LabeledPolynomial::new(
        "dense".to_string(),
        DenseMultilinearExtension::rand(NUM_VARS, &mut rng),
        None,
        None,
    );
    commit_open_check::<4, _>(NUM_VARS, &dense, &mut rng);

    let sparse = LabeledPolynomial::new(
        "sparse".to_string(),
        SparseMultilinearExtension::rand_with_config(NUM_VARS, 1 << 4, &mut rng),
        None,
        None,
    );
    commit_open_check::<4, _>(NUM_VARS, &sparse, &mut rng);
}

#[test]
fn commit_open_check_supports_an_uneven_tensor_partition() {
    let mut rng = ChaCha20Rng::from_rng(test_rng()).unwrap();
    const NUM_VARS: usize = 5;

    let dense = LabeledPolynomial::new(
        "uneven".to_string(),
        DenseMultilinearExtension::rand(NUM_VARS, &mut rng),
        None,
        None,
    );
    commit_open_check::<3, _>(NUM_VARS, &dense, &mut rng);
}

#[test]
fn check_rejects_an_incorrect_evaluation() {
    let mut rng = ChaCha20Rng::from_rng(test_rng()).unwrap();
    const NUM_VARS: usize = 8;

    let params = DenseKZH::<4>::setup(1, Some(NUM_VARS), &mut rng).unwrap();
    let (ck, vk) = DenseKZH::<4>::trim(&params, 1, 0, None).unwrap();
    let polynomial = LabeledPolynomial::new(
        "wrong-value".to_string(),
        DenseMultilinearExtension::rand(NUM_VARS, &mut rng),
        None,
        None,
    );
    let point: Vec<Fr> = (0..NUM_VARS).map(|_| Fr::rand(&mut rng)).collect();
    let value = polynomial.evaluate(&point);
    let (commitments, states) = DenseKZH::<4>::commit(&ck, [&polynomial], None).unwrap();
    let base_sponge = poseidon_sponge_for_test::<Fr>();
    let proof = DenseKZH::<4>::open(
        &ck,
        [&polynomial],
        &commitments,
        &point,
        &mut base_sponge.clone(),
        &states,
        None,
    )
    .unwrap();

    assert!(!DenseKZH::<4>::check(
        &vk,
        &commitments,
        &point,
        [value + Fr::one()],
        &proof,
        &mut base_sponge.clone(),
        None,
    )
    .unwrap());
}
