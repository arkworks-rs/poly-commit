//! KZH correctness, compatibility, cost, and API tests.
//!
//! The tensor-oracle tests intentionally reconstruct commitments and opening
//! layers without calling KZH helpers. Keeping that arithmetic explicit makes
//! variable-order and indexing regressions observable. Higher-level API tests
//! keep the complete `commit`/`open`/`check` flow visible for reviewers.

use super::utils::{
    arkworks_lagrange_evaluations, auxiliary_prefix_lengths, balanced_block_sizes,
    block_dimensions, opening_work_plan, product, OpeningLayerSource,
};
use crate::{
    kzh::{Commitment, CommitterKey, PreparedVerifierKey, UniversalParams, VerifierKey, KZH},
    tests::poseidon_sponge_for_test,
    Error, Evaluations, LCTerm, LabeledCommitment, LabeledPolynomial, LinearCombination,
    PCCommitment, PolynomialCommitment, QuerySet,
};
use ark_bls12_381::{Bls12_381, Fr, G1Affine, G1Projective};
use ark_crypto_primitives::sponge::CryptographicSponge;
use ark_ec::{pairing::Pairing, CurveGroup, PrimeGroup, VariableBaseMSM};
use ark_ff::{Field, One, Zero};
use ark_poly::{
    DenseMultilinearExtension, MultilinearExtension, Polynomial, SparseMultilinearExtension,
};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, Compress, Validate};
use core::fmt::Debug;
use rand_chacha::{rand_core::SeedableRng, ChaCha20Rng};

type TestPolynomial = DenseMultilinearExtension<Fr>;
type TestKZH<const K: usize> = KZH<Bls12_381, TestPolynomial, K>;
type SparseTestPolynomial = SparseMultilinearExtension<Fr>;
type SparseTestKZH<const K: usize> = KZH<Bls12_381, SparseTestPolynomial, K>;

// Shared deterministic fixtures and independent oracles.

fn seeded_rng(seed: u8) -> ChaCha20Rng {
    ChaCha20Rng::from_seed([seed; 32])
}

fn boolean_hypercube_point(num_vars: usize, index: usize) -> Vec<Fr> {
    (0..num_vars)
        .map(|bit| Fr::from(((index >> bit) & 1) as u64))
        .collect()
}

fn non_boolean_point(num_vars: usize, offset: u64) -> Vec<Fr> {
    (0..num_vars)
        .map(|index| Fr::from(offset + 3 * index as u64 + 2))
        .collect()
}

fn explicit_equality_weight(point: &[Fr], assignment: usize) -> Fr {
    point
        .iter()
        .enumerate()
        .fold(Fr::one(), |weight, (bit, coordinate)| {
            if assignment & (1usize << bit) == 0 {
                weight * (Fr::one() - coordinate)
            } else {
                weight * coordinate
            }
        })
}

fn explicit_prefix_weight(
    point_blocks: &[&[Fr]],
    dimensions: &[usize],
    mut flat_assignment: usize,
) -> Fr {
    let mut weight = Fr::one();
    for block in 0..point_blocks.len() {
        let assignment = flat_assignment % dimensions[block];
        flat_assignment /= dimensions[block];
        weight *= explicit_equality_weight(point_blocks[block], assignment);
    }
    assert_eq!(flat_assignment, 0);
    weight
}

fn explicit_auxiliary_tables<const K: usize>(
    committer_key: &CommitterKey<Bls12_381, K>,
    evaluations: &[Fr],
    dimensions: &[usize],
    level_count: usize,
) -> Vec<Vec<G1Affine>> {
    (0..level_count)
        .map(|level| {
            let prefix_dimension = dimensions[..level].iter().product::<usize>();
            let current_dimension = dimensions[level];
            let suffix_dimension = dimensions[level + 1..].iter().product::<usize>();
            assert_eq!(committer_key.h[level + 1].len(), suffix_dimension);

            let projective = (0..current_dimension)
                .flat_map(|current_assignment| {
                    (0..prefix_dimension).map(move |prefix_assignment| {
                        let coefficients = (0..suffix_dimension)
                            .map(|suffix_assignment| {
                                let index = prefix_assignment
                                    + prefix_dimension
                                        * (current_assignment
                                            + current_dimension * suffix_assignment);
                                evaluations[index]
                            })
                            .collect::<Vec<_>>();
                        <G1Projective as VariableBaseMSM>::msm(
                            &committer_key.h[level + 1],
                            &coefficients,
                        )
                        .unwrap()
                    })
                })
                .collect::<Vec<_>>();
            G1Projective::normalize_batch(&projective)
        })
        .collect()
}

fn explicit_opening_layers<const K: usize>(
    committer_key: &CommitterKey<Bls12_381, K>,
    evaluations: &[Fr],
    point_blocks: &[&[Fr]],
    dimensions: &[usize],
) -> Vec<Vec<G1Affine>> {
    assert_eq!(dimensions.len(), K);
    assert_eq!(point_blocks.len(), K);

    (0..K - 1)
        .map(|level| {
            let prefix_dimension = dimensions[..level].iter().product::<usize>();
            let current_dimension = dimensions[level];
            let suffix_dimension = dimensions[level + 1..].iter().product::<usize>();

            (0..current_dimension)
                .map(|current_assignment| {
                    let mut coefficients = vec![Fr::zero(); suffix_dimension];
                    for prefix_assignment in 0..prefix_dimension {
                        let prefix_weight = explicit_prefix_weight(
                            &point_blocks[..level],
                            &dimensions[..level],
                            prefix_assignment,
                        );
                        for (suffix_assignment, coefficient) in coefficients.iter_mut().enumerate()
                        {
                            let index = prefix_assignment
                                + prefix_dimension
                                    * (current_assignment + current_dimension * suffix_assignment);
                            *coefficient += prefix_weight * evaluations[index];
                        }
                    }
                    <G1Projective as VariableBaseMSM>::msm(
                        &committer_key.h[level + 1],
                        &coefficients,
                    )
                    .unwrap()
                    .into_affine()
                })
                .collect()
        })
        .collect()
}

fn explicit_final_evaluations(
    evaluations: &[Fr],
    point_blocks: &[&[Fr]],
    dimensions: &[usize],
) -> Vec<Fr> {
    let final_level = dimensions.len() - 1;
    let prefix_dimension = dimensions[..final_level].iter().product::<usize>();
    (0..dimensions[final_level])
        .map(|final_assignment| {
            (0..prefix_dimension).fold(Fr::zero(), |sum, prefix_assignment| {
                let prefix_weight = explicit_prefix_weight(
                    &point_blocks[..final_level],
                    &dimensions[..final_level],
                    prefix_assignment,
                );
                sum + prefix_weight
                    * evaluations[prefix_assignment + prefix_dimension * final_assignment]
            })
        })
        .collect()
}

fn add_generator(point: G1Affine) -> G1Affine {
    (G1Projective::from(point) + G1Projective::generator()).into_affine()
}

/// Evaluate a deliberately nontrivial multilinear polynomial directly from
/// its closed-form expression. This is independent of both ark-poly's MLE
/// evaluator and KZH's tensor folding code.
fn structured_evaluation(point: &[Fr], salt: u64) -> Fr {
    let mut result = Fr::from(17 + salt);

    for (index, coordinate) in point.iter().enumerate() {
        let coefficient = Fr::from(salt + 3 * index as u64 + 5);
        result += coefficient * coordinate;
    }

    for index in 0..point.len().saturating_sub(1) {
        let coefficient = Fr::from(2 * salt + 5 * index as u64 + 11);
        result += coefficient * point[index] * point[index + 1];
    }

    if point.len() >= 3 {
        result +=
            Fr::from(3 * salt + 29) * point[0] * point[point.len() / 2] * point[point.len() - 1];
    }

    result
}

fn structured_polynomial(
    label: &str,
    num_vars: usize,
    salt: u64,
) -> LabeledPolynomial<Fr, TestPolynomial> {
    let evaluations = (0..(1usize << num_vars))
        .map(|index| structured_evaluation(&boolean_hypercube_point(num_vars, index), salt))
        .collect();
    let polynomial = TestPolynomial::from_evaluations_vec(num_vars, evaluations);
    LabeledPolynomial::new(label.to_string(), polynomial, None, None)
}

fn setup_and_trim<const K: usize>(
    num_vars: usize,
    seed: u8,
) -> (
    UniversalParams<Bls12_381, K>,
    CommitterKey<Bls12_381, K>,
    VerifierKey<Bls12_381, K>,
) {
    let mut rng = seeded_rng(seed);
    let params = TestKZH::<K>::setup(1, Some(num_vars), &mut rng).unwrap();
    let (committer_key, verifier_key) = TestKZH::<K>::trim(&params, 1, 0, None).unwrap();
    (params, committer_key, verifier_key)
}

mod batching;
mod costs;
mod layout;
mod serialization;
mod validation;
