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
    kzh::{
        KZHCommitment, KZHCommitterKey, KZHPreparedVerifierKey, KZHUniversalParams, KZHVerifierKey,
        KZH,
    },
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
    committer_key: &KZHCommitterKey<Bls12_381, K>,
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
    committer_key: &KZHCommitterKey<Bls12_381, K>,
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
    KZHUniversalParams<Bls12_381, K>,
    KZHCommitterKey<Bls12_381, K>,
    KZHVerifierKey<Bls12_381, K>,
) {
    let mut rng = seeded_rng(seed);
    let params = TestKZH::<K>::setup(1, Some(num_vars), &mut rng).unwrap();
    let (committer_key, verifier_key) = TestKZH::<K>::trim(&params, 1, 0, None).unwrap();
    (params, committer_key, verifier_key)
}

// Arkworks representation and execution compatibility.

#[test]
fn lagrange_evaluations_match_arkworks_unit_basis_polynomials() {
    let point = [Fr::from(2u64), Fr::from(5u64), Fr::from(11u64)];
    assert!(point
        .iter()
        .all(|coordinate| !coordinate.is_zero() && !coordinate.is_one()));

    let dimension = 1usize << point.len();
    let expected: Vec<_> = (0..dimension)
        .map(|assignment| {
            let mut basis_evaluations = vec![Fr::zero(); dimension];
            basis_evaluations[assignment] = Fr::one();
            TestPolynomial::from_evaluations_vec(point.len(), basis_evaluations)
                .evaluate(&point.to_vec())
        })
        .collect();

    assert_eq!(arkworks_lagrange_evaluations(&point).unwrap(), expected);

    // Arkworks collapses a dense MLE's arity when scalar multiplication makes
    // it identically zero. The adapter must nevertheless retain both halves
    // of the two-variable equality table at Boolean coordinates.
    assert_eq!(
        arkworks_lagrange_evaluations(&[Fr::zero(), Fr::one()]).unwrap(),
        [Fr::zero(), Fr::zero(), Fr::one(), Fr::zero()]
    );
}

#[test]
fn dense_and_sparse_arkworks_mles_produce_identical_kzh_artifacts() {
    const NUM_VARS: usize = 6;
    const K: usize = 3;
    let (_, committer_key, verifier_key) = setup_and_trim::<K>(NUM_VARS, 71);

    // Use isolated, nonsymmetric entries whose indices exercise both low and
    // high variable bits. Arkworks must expose the same canonical evaluation
    // vector for its dense and sparse MLE representations.
    let sparse_entries = [
        (0usize, Fr::from(3u64)),
        (1, Fr::from(5u64)),
        (6, Fr::from(11u64)),
        (19, Fr::from(17u64)),
        (40, Fr::from(23u64)),
        (63, Fr::from(31u64)),
    ];
    let mut evaluations = vec![Fr::zero(); 1usize << NUM_VARS];
    for &(index, value) in &sparse_entries {
        evaluations[index] = value;
    }

    let dense = LabeledPolynomial::new(
        "arkworks-representation".to_string(),
        TestPolynomial::from_evaluations_vec(NUM_VARS, evaluations.clone()),
        None,
        None,
    );
    let sparse = LabeledPolynomial::new(
        "arkworks-representation".to_string(),
        SparseTestPolynomial::from_evaluations(NUM_VARS, &sparse_entries),
        None,
        None,
    );
    assert_eq!(dense.polynomial().to_evaluations(), evaluations);
    assert_eq!(sparse.polynomial().to_evaluations(), evaluations);

    let point = non_boolean_point(NUM_VARS, 101);
    let expected_value = dense.evaluate(&point);
    assert_eq!(sparse.evaluate(&point), expected_value);

    let (dense_commitments, dense_states) =
        TestKZH::<K>::commit(&committer_key, [&dense], None).unwrap();
    let (sparse_commitments, sparse_states) =
        SparseTestKZH::<K>::commit(&committer_key, [&sparse], None).unwrap();
    assert_eq!(dense_commitments.len(), sparse_commitments.len());
    for (dense_commitment, sparse_commitment) in dense_commitments.iter().zip(&sparse_commitments) {
        assert_eq!(dense_commitment.label(), sparse_commitment.label());
        assert_eq!(
            dense_commitment.commitment(),
            sparse_commitment.commitment()
        );
        assert_eq!(
            dense_commitment.degree_bound(),
            sparse_commitment.degree_bound()
        );
    }
    assert_eq!(dense_states, sparse_states);

    let base_sponge = poseidon_sponge_for_test::<Fr>();
    let dense_proof = TestKZH::<K>::open(
        &committer_key,
        [&dense],
        &dense_commitments,
        &point,
        &mut base_sponge.clone(),
        &dense_states,
        None,
    )
    .unwrap();
    let sparse_proof = SparseTestKZH::<K>::open(
        &committer_key,
        [&sparse],
        &sparse_commitments,
        &point,
        &mut base_sponge.clone(),
        &sparse_states,
        None,
    )
    .unwrap();
    assert_eq!(dense_proof, sparse_proof);

    assert!(TestKZH::<K>::check(
        &verifier_key,
        &sparse_commitments,
        &point,
        [expected_value],
        &sparse_proof,
        &mut base_sponge.clone(),
        None,
    )
    .unwrap());
    assert!(SparseTestKZH::<K>::check(
        &verifier_key,
        &dense_commitments,
        &point,
        [expected_value],
        &dense_proof,
        &mut base_sponge.clone(),
        None,
    )
    .unwrap());
}

#[cfg(feature = "parallel")]
#[test]
fn parallel_auxiliary_construction_matches_the_serial_path() {
    const NUM_VARS: usize = 8;
    const K: usize = 6;
    let (_, committer_key, _) = setup_and_trim::<K>(NUM_VARS, 52);
    let polynomial = structured_polynomial("parallel-auxiliary", NUM_VARS, 29);

    let commit_with_threads = |num_threads| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(num_threads)
            .build()
            .unwrap()
            .install(|| TestKZH::<K>::commit(&committer_key, [&polynomial], None).unwrap())
    };

    // A one-thread local pool exercises the serial auxiliary-table path; four
    // threads exercise the bounded row-parallel path. Both must produce the
    // exact same commitments and cached state.
    let (serial_commitments, serial_states) = commit_with_threads(1);
    let (parallel_commitments, parallel_states) = commit_with_threads(4);
    assert_eq!(serial_commitments.len(), parallel_commitments.len());
    for (serial, parallel) in serial_commitments.iter().zip(&parallel_commitments) {
        assert_eq!(serial.label(), parallel.label());
        assert_eq!(serial.commitment(), parallel.commitment());
        assert_eq!(serial.degree_bound(), parallel.degree_bound());
    }
    assert_eq!(serial_states, parallel_states);
}

// Tensor decomposition, stored-key shape, and asymptotic work bounds.

fn assert_independent_end_to_end<const K: usize>(
    num_vars: usize,
    expected_block_sizes: &[usize],
    seed: u8,
) {
    let (params, committer_key, verifier_key) = setup_and_trim::<K>(num_vars, seed);
    assert_eq!(params.num_vars_per_block, expected_block_sizes);
    assert_eq!(committer_key.num_vars_per_block, expected_block_sizes);
    assert_eq!(verifier_key.num_vars_per_block, expected_block_sizes);

    let salt = seed as u64 + 7;
    let polynomial = structured_polynomial("structured", num_vars, salt);
    let point = non_boolean_point(num_vars, salt + 3);
    assert!(point
        .iter()
        .all(|coordinate| !coordinate.is_zero() && !coordinate.is_one()));
    let expected = structured_evaluation(&point, salt);

    // This assertion checks the independently computed expression against
    // ark-poly before either value is passed to KZH.
    assert_eq!(polynomial.evaluate(&point), expected);

    let (commitments, states) = TestKZH::<K>::commit(&committer_key, [&polynomial], None).unwrap();
    let base_sponge = poseidon_sponge_for_test::<Fr>();
    let proof = TestKZH::<K>::open(
        &committer_key,
        [&polynomial],
        &commitments,
        &point,
        &mut base_sponge.clone(),
        &states,
        None,
    )
    .unwrap();

    let block_sizes = balanced_block_sizes(num_vars, K).unwrap();
    assert_eq!(block_sizes, expected_block_sizes);
    assert_eq!(block_sizes.len(), K);
    let dimensions = block_dimensions(&block_sizes).unwrap();
    assert_eq!(proof.layers().len(), K - 1);
    for (layer, expected_dimension) in proof.layers().iter().zip(&dimensions[..K - 1]) {
        assert_eq!(layer.len(), *expected_dimension);
    }
    assert_eq!(proof.final_evaluations.len(), dimensions[K - 1]);

    // KZH consumes the lowest Arkworks variables first. After the first K - 1
    // blocks, its disclosed field table must therefore be exactly the result
    // of ark-poly's canonical low-prefix binding. Uneven decompositions make
    // this assertion sensitive to both block direction and block boundaries.
    let fixed_prefix_len = block_sizes[..K - 1].iter().sum::<usize>();
    let arkworks_final = polynomial
        .polynomial()
        .fix_variables(&point[..fixed_prefix_len]);
    assert_eq!(proof.final_evaluations, arkworks_final.evaluations);
    assert_eq!(
        arkworks_final.evaluate(&point[fixed_prefix_len..].to_vec()),
        expected
    );

    assert!(TestKZH::<K>::check(
        &verifier_key,
        &commitments,
        &point,
        [expected],
        &proof,
        &mut base_sponge.clone(),
        None,
    )
    .unwrap());
}

#[test]
fn kzh2_supports_an_uneven_tensor_partition() {
    assert_independent_end_to_end::<2>(5, &[3, 2], 2);
}

#[test]
fn kzh2_supports_the_minimal_valid_tensor_partition() {
    assert_independent_end_to_end::<2>(2, &[1, 1], 1);
}

#[test]
fn kzh3_supports_an_uneven_tensor_partition() {
    assert_independent_end_to_end::<3>(7, &[3, 2, 2], 3);
}

#[test]
fn kzh4_supports_an_uneven_tensor_partition() {
    assert_independent_end_to_end::<4>(6, &[2, 2, 1, 1], 4);
}

#[test]
fn kzh5_supports_an_uneven_tensor_partition() {
    assert_independent_end_to_end::<5>(7, &[2, 2, 1, 1, 1], 5);
}

fn assert_v_tau_has_only_proof_layers<const K: usize>(num_vars: usize, seed: u8) {
    let (params, _, verifier_key) = setup_and_trim::<K>(num_vars, seed);
    let block_sizes = balanced_block_sizes(num_vars, K).unwrap();
    assert_eq!(block_sizes.len(), K);
    assert_eq!(params.v_tau.len(), K - 1);
    assert_eq!(verifier_key.v_tau.len(), K - 1);

    for ((setup_layer, verifier_layer), num_block_vars) in params
        .v_tau
        .iter()
        .zip(&verifier_key.v_tau)
        .zip(&block_sizes[..K - 1])
    {
        let expected_dimension = 1usize << num_block_vars;
        assert_eq!(setup_layer.len(), expected_dimension);
        assert_eq!(verifier_layer.len(), expected_dimension);
        assert_eq!(setup_layer, verifier_layer);
    }
}

#[test]
fn setup_and_trim_store_exactly_k_minus_one_verifier_layers() {
    assert_v_tau_has_only_proof_layers::<2>(5, 6);
    assert_v_tau_has_only_proof_layers::<3>(6, 7);
    assert_v_tau_has_only_proof_layers::<4>(8, 7);
    assert_v_tau_has_only_proof_layers::<5>(7, 8);
    assert_v_tau_has_only_proof_layers::<5>(10, 9);
}

fn assert_symbolic_cost_bounds<const K: usize>(num_vars: usize) {
    let block_sizes = balanced_block_sizes(num_vars, K).unwrap();
    assert_eq!(block_sizes.len(), K);
    assert_eq!(block_sizes.iter().sum::<usize>(), num_vars);

    let dimensions = block_dimensions(&block_sizes).unwrap();
    assert_eq!(dimensions.len(), K);
    let num_evaluations = product(&dimensions).unwrap();
    assert_eq!(num_evaluations, 1usize << num_vars);

    let plan = opening_work_plan(&dimensions, 1).unwrap();
    assert_eq!(plan.len(), K - 1);
    let auxiliary_lengths = auxiliary_prefix_lengths(&dimensions).unwrap();
    let cached_count = plan
        .iter()
        .take_while(|layer| layer.source == OpeningLayerSource::Cached)
        .count();
    assert_eq!(cached_count, auxiliary_lengths.len());
    assert!(plan[cached_count..]
        .iter()
        .all(|layer| layer.source == OpeningLayerSource::Direct));

    let mut prefix_dimension = 1usize;
    let mut suffix_dimension = num_evaluations;
    for (level, (&dimension, layer)) in dimensions.iter().zip(&plan).enumerate() {
        prefix_dimension *= dimension;
        if level < cached_count {
            assert_eq!(auxiliary_lengths[level], prefix_dimension);
            let expected_terms = if level == 0 { 0 } else { prefix_dimension };
            assert_eq!(layer.scalar_terms, expected_terms);
        } else {
            assert_eq!(layer.scalar_terms, suffix_dimension);
        }
        suffix_dimension /= dimension;
    }

    // Every retained auxiliary table is one N-term preprocessing pass. The
    // omitted suffix consists precisely of tables unused by a generic opening.
    let commitment_scalar_terms = num_evaluations * (1 + auxiliary_lengths.len());
    assert!(commitment_scalar_terms <= K * num_evaluations);

    // Write q = ceil(n / K). Every balanced dimension is at most 2^q, and the
    // cheaper cached/direct construction at a transition spans at most
    // ceil(K / 2) axes. Hence each layer has at most
    // (2^q)^ceil(K/2) <= 2^ceil(K/2) * N^(ceil(K/2)/K) group-scalar terms.
    let q = num_vars.div_ceil(K);
    let kth_root_bound = 1usize << q;
    let half_axis_count = K.div_ceil(2);
    let opening_layer_bound = kth_root_bound.checked_pow(half_axis_count as u32).unwrap();
    assert!(plan
        .iter()
        .all(|layer| layer.scalar_terms <= opening_layer_bound));

    // Cached terms grow geometrically by at least two and direct terms shrink
    // geometrically by at least two. Thus the two sides together cost less
    // than four times either valid maximum-layer bound. The crossing-block
    // argument gives the sharper shape-aware maximum
    // M = 2^floor((n + ceil(n/K)) / 2), including uneven partitions.
    let crossing_layer_bound = 1usize << ((num_vars + q) / 2);
    assert!(plan
        .iter()
        .all(|layer| layer.scalar_terms <= crossing_layer_bound));
    let total_opening_terms = plan.iter().map(|layer| layer.scalar_terms).sum::<usize>();
    assert!(total_opening_terms < 4 * crossing_layer_bound);
    // This second form records the requested fixed-K asymptotic directly.
    assert!(total_opening_terms < 4 * opening_layer_bound);

    // The plan counts online group-scalar terms. Dense tensor folding is the
    // separately documented O(N) field work and is intentionally not inferred
    // from this group-work schedule.

    // There are exactly K verifier axes, each of dimension at most 2^q.
    let verifier_msm_terms = dimensions.iter().sum::<usize>();
    assert!(verifier_msm_terms <= K * kth_root_bound);
    let verifier_pairing_inputs = dimensions[..K - 1].iter().sum::<usize>() + K - 1;
    assert!(verifier_pairing_inputs <= K * (kth_root_bound + 1));
}

fn assert_family_cost_sweep<const K: usize>() {
    for num_vars in K..=K + 12 {
        assert_symbolic_cost_bounds::<K>(num_vars);
    }
}

#[test]
fn generic_cost_bounds_hold_for_kzh2_through_kzh8() {
    assert_family_cost_sweep::<2>();
    assert_family_cost_sweep::<3>();
    assert_family_cost_sweep::<4>();
    assert_family_cost_sweep::<5>();
    assert_family_cost_sweep::<6>();
    assert_family_cost_sweep::<7>();
    assert_family_cost_sweep::<8>();
}

fn assert_exact_generic_costs<const K: usize>(
    num_vars: usize,
    expected_axes: &[usize],
    expected_dimensions: &[usize],
    expected_auxiliary_lengths: &[usize],
    expected_opening_terms: &[usize],
    seed: u8,
) {
    assert_symbolic_cost_bounds::<K>(num_vars);
    let (params, committer_key, verifier_key) = setup_and_trim::<K>(num_vars, seed);
    assert_eq!(expected_axes.len(), K);
    assert_eq!(balanced_block_sizes(num_vars, K).unwrap(), expected_axes);
    assert_eq!(params.num_vars_per_block, expected_axes);
    assert_eq!(committer_key.num_vars_per_block, expected_axes);
    assert_eq!(verifier_key.num_vars_per_block, expected_axes);
    assert_eq!(params.h.len(), K);
    assert_eq!(params.v_tau.len(), K - 1);
    assert_eq!(verifier_key.v_tau.len(), K - 1);

    let dimensions = block_dimensions(expected_axes).unwrap();
    assert_eq!(dimensions, expected_dimensions);
    let plan = opening_work_plan(&dimensions, 1).unwrap();
    assert_eq!(
        plan.iter()
            .map(|layer| layer.scalar_terms)
            .collect::<Vec<_>>(),
        expected_opening_terms
    );
    assert!(plan[..expected_auxiliary_lengths.len()]
        .iter()
        .all(|layer| layer.source == OpeningLayerSource::Cached));
    assert!(plan[expected_auxiliary_lengths.len()..]
        .iter()
        .all(|layer| layer.source == OpeningLayerSource::Direct));

    let salt = seed as u64 + 17;
    let polynomial = structured_polynomial("generic-cost", num_vars, salt);
    let point = non_boolean_point(num_vars, salt + 31);
    assert!(point
        .iter()
        .all(|coordinate| !coordinate.is_zero() && !coordinate.is_one()));
    let expected_value = structured_evaluation(&point, salt);
    let (commitments, states) = TestKZH::<K>::commit(&committer_key, [&polynomial], None).unwrap();
    assert_eq!(states.len(), 1);
    assert_eq!(
        states[0]
            .auxiliary_tables
            .iter()
            .map(Vec::len)
            .collect::<Vec<_>>(),
        expected_auxiliary_lengths
    );

    let num_evaluations = 1usize << num_vars;
    let commitment_scalar_terms = num_evaluations * (1 + states[0].auxiliary_tables.len());
    assert!(commitment_scalar_terms <= K * num_evaluations);

    let base_sponge = poseidon_sponge_for_test::<Fr>();
    let proof = TestKZH::<K>::open(
        &committer_key,
        [&polynomial],
        &commitments,
        &point,
        &mut base_sponge.clone(),
        &states,
        None,
    )
    .unwrap();
    assert_eq!(proof.layers().len(), K - 1);
    assert_eq!(
        proof.layers().iter().map(Vec::len).collect::<Vec<_>>(),
        dimensions[..K - 1]
    );
    assert_eq!(proof.final_evaluations.len(), dimensions[K - 1]);

    let verifier_msm_terms =
        proof.layers().iter().map(Vec::len).sum::<usize>() + proof.final_evaluations.len();
    let kth_root_bound = 1usize << num_vars.div_ceil(K);
    assert!(verifier_msm_terms <= K * kth_root_bound);
    let verifier_pairing_inputs = proof
        .layers()
        .iter()
        .map(|layer| layer.len() + 1)
        .sum::<usize>();
    assert!(verifier_pairing_inputs <= K * (kth_root_bound + 1));

    assert!(TestKZH::<K>::check(
        &verifier_key,
        &commitments,
        &point,
        [expected_value],
        &proof,
        &mut base_sponge.clone(),
        None,
    )
    .unwrap());
}

#[test]
fn odd_k_openings_have_expected_block_shapes_and_work_counts() {
    assert_exact_generic_costs::<3>(6, &[2, 2, 2], &[4, 4, 4], &[4], &[0, 16], 9);
    assert_exact_generic_costs::<5>(
        10,
        &[2, 2, 2, 2, 2],
        &[4, 4, 4, 4, 4],
        &[4, 16],
        &[0, 16, 64, 16],
        10,
    );
}

#[test]
fn kzh6_opening_matches_arkworks_low_variable_order() {
    const NUM_VARS: usize = 6;
    const K: usize = 6;
    let (_, committer_key, verifier_key) = setup_and_trim::<K>(NUM_VARS, 11);
    assert_eq!(committer_key.num_vars_per_block, [1; K]);

    let dimensions = block_dimensions(&committer_key.num_vars_per_block).unwrap();
    assert_eq!(dimensions, [2; K]);
    assert_eq!(
        opening_work_plan(&dimensions, 1)
            .unwrap()
            .iter()
            .map(|layer| layer.source)
            .collect::<Vec<_>>(),
        [
            OpeningLayerSource::Cached,
            OpeningLayerSource::Cached,
            OpeningLayerSource::Cached,
            OpeningLayerSource::Direct,
            OpeningLayerSource::Direct,
        ]
    );

    let polynomial = structured_polynomial("kzh6-native", NUM_VARS, 43);
    let point = non_boolean_point(NUM_VARS, 59);
    let expected_value = structured_evaluation(&point, 43);
    let (commitments, states) = TestKZH::<K>::commit(&committer_key, [&polynomial], None).unwrap();
    assert_eq!(
        states[0]
            .auxiliary_tables
            .iter()
            .map(Vec::len)
            .collect::<Vec<_>>(),
        [2, 4, 8]
    );

    let base_sponge = poseidon_sponge_for_test::<Fr>();
    let proof = TestKZH::<K>::open(
        &committer_key,
        [&polynomial],
        &commitments,
        &point,
        &mut base_sponge.clone(),
        &states,
        None,
    )
    .unwrap();

    // Reconstruct every layer from ark-poly's native low-prefix table. In
    // particular, layer two contracts a cached table with the Kronecker
    // weights for two already-fixed blocks, so reversing the prefix-product
    // order changes this exact proof vector.
    let mut staged = polynomial.polynomial().clone();
    let mut fixed_low_num_vars = 0usize;
    let mut expected_layers = Vec::with_capacity(K - 1);
    for (level, &dimension) in dimensions.iter().take(K - 1).enumerate() {
        assert_eq!(
            staged,
            polynomial
                .polynomial()
                .fix_variables(&point[..fixed_low_num_vars])
        );
        let suffix_dimension = committer_key.h[level + 1].len();
        let projective = (0..dimension)
            .map(|current_assignment| {
                let coefficients = (0..suffix_dimension)
                    .map(|suffix_assignment| {
                        staged.evaluations[current_assignment + dimension * suffix_assignment]
                    })
                    .collect::<Vec<_>>();
                <G1Projective as VariableBaseMSM>::msm(&committer_key.h[level + 1], &coefficients)
                    .unwrap()
            })
            .collect::<Vec<_>>();
        expected_layers.push(G1Projective::normalize_batch(&projective));

        staged = staged.fix_variables(&point[fixed_low_num_vars..fixed_low_num_vars + 1]);
        fixed_low_num_vars += 1;
    }

    assert_eq!(proof.layers(), expected_layers);
    assert_eq!(proof.final_evaluations, staged.evaluations);
    assert_eq!(fixed_low_num_vars, K - 1);
    assert!(TestKZH::<K>::check(
        &verifier_key,
        &commitments,
        &point,
        [expected_value],
        &proof,
        &mut base_sponge.clone(),
        None,
    )
    .unwrap());
}

// Independent tensor-order and generic-opening oracles.

#[test]
fn unequal_dimensions_match_complete_native_order_oracle() {
    const NUM_VARS: usize = 7;
    const K: usize = 6;
    let (_, committer_key, verifier_key) = setup_and_trim::<K>(NUM_VARS, 14);
    assert_eq!(committer_key.num_vars_per_block, [2, 1, 1, 1, 1, 1]);

    let dimensions = block_dimensions(&committer_key.num_vars_per_block).unwrap();
    assert_eq!(dimensions, [4, 2, 2, 2, 2, 2]);
    assert_eq!(
        opening_work_plan(&dimensions, 1)
            .unwrap()
            .iter()
            .map(|layer| layer.source)
            .collect::<Vec<_>>(),
        [
            OpeningLayerSource::Cached,
            OpeningLayerSource::Cached,
            OpeningLayerSource::Direct,
            OpeningLayerSource::Direct,
            OpeningLayerSource::Direct,
        ]
    );

    let polynomial = structured_polynomial("unequal-native-oracle", NUM_VARS, 47);
    let evaluations = polynomial.polynomial().evaluations.clone();
    let point = non_boolean_point(NUM_VARS, 67);
    let point_blocks = [
        &point[0..2],
        &point[2..3],
        &point[3..4],
        &point[4..5],
        &point[5..6],
        &point[6..7],
    ];
    let expected_value = structured_evaluation(&point, 47);
    let (commitments, states) = TestKZH::<K>::commit(&committer_key, [&polynomial], None).unwrap();
    assert_eq!(
        states[0]
            .auxiliary_tables
            .iter()
            .map(Vec::len)
            .collect::<Vec<_>>(),
        [4, 8]
    );

    // Independently reconstruct every cached suffix commitment from the raw
    // arkworks evaluation vector. Here the second table is especially useful:
    // its four-entry prior prefix and two-entry current block make every axis
    // swap observable instead of being hidden by equal tensor dimensions.
    assert_eq!(
        states[0].auxiliary_tables,
        explicit_auxiliary_tables(
            &committer_key,
            &evaluations,
            &dimensions,
            states[0].auxiliary_tables.len(),
        )
    );

    let base_sponge = poseidon_sponge_for_test::<Fr>();
    let proof = TestKZH::<K>::open(
        &committer_key,
        [&polynomial],
        &commitments,
        &point,
        &mut base_sponge.clone(),
        &states,
        None,
    )
    .unwrap();

    // Derive every proof layer directly from the raw evaluation table and the
    // explicit Boolean-Lagrange formula. This covers both cached layers and all
    // three direct layers without calling KZH's folding or equality helpers.
    assert_eq!(
        proof.layers(),
        explicit_opening_layers(&committer_key, &evaluations, &point_blocks, &dimensions)
    );
    assert_eq!(
        proof.final_evaluations,
        explicit_final_evaluations(&evaluations, &point_blocks, &dimensions)
    );

    let fixed_prefix_len = committer_key.num_vars_per_block[..K - 1]
        .iter()
        .sum::<usize>();
    assert_eq!(fixed_prefix_len, NUM_VARS - 1);
    assert_eq!(
        proof.final_evaluations,
        polynomial
            .polynomial()
            .fix_variables(&point[..fixed_prefix_len])
            .evaluations
    );
    assert!(TestKZH::<K>::check(
        &verifier_key,
        &commitments,
        &point,
        [expected_value],
        &proof,
        &mut base_sponge.clone(),
        None,
    )
    .unwrap());
}

#[test]
fn generic_opening_layers_match_explicit_and_arkworks_oracles() {
    const NUM_VARS: usize = 8;
    const K: usize = 4;
    let (_, committer_key, verifier_key) = setup_and_trim::<K>(NUM_VARS, 12);
    let polynomial = structured_polynomial("tensor-oracle", NUM_VARS, 37);
    let evaluations = polynomial.polynomial().evaluations.clone();
    let point = non_boolean_point(NUM_VARS, 73);
    assert!(point
        .iter()
        .all(|coordinate| !coordinate.is_zero() && !coordinate.is_one()));

    let (commitments, states) = TestKZH::<K>::commit(&committer_key, [&polynomial], None).unwrap();
    assert_eq!(states[0].auxiliary_tables.len(), 2);
    assert_eq!(states[0].auxiliary_tables[0].len(), 4);
    assert_eq!(states[0].auxiliary_tables[1].len(), 16);
    let base_sponge = poseidon_sponge_for_test::<Fr>();
    let proof = TestKZH::<K>::open(
        &committer_key,
        [&polynomial],
        &commitments,
        &point,
        &mut base_sponge.clone(),
        &states,
        None,
    )
    .unwrap();

    // Tensor axes follow arkworks directly: the first block contains the
    // lowest variables, and the first axis is the fastest-varying one.
    let dimensions = [4usize, 4, 4, 4];
    let point_blocks = [&point[0..2], &point[2..4], &point[4..6], &point[6..8]];

    // Independently reconstruct every stored auxiliary entry. For a low
    // prefix assignment p, current-block assignment a, and high-suffix
    // assignment s, arkworks stores the original evaluation at
    // p + P * (a + d * s), where P is the prior-prefix dimension and d is the
    // current block dimension. Auxiliary tables are current-major, then
    // prior-prefix-minor, which is also the canonical combined-prefix order.
    assert_eq!(
        states[0].auxiliary_tables,
        explicit_auxiliary_tables(
            &committer_key,
            &evaluations,
            &dimensions,
            states[0].auxiliary_tables.len(),
        )
    );
    assert_eq!(
        proof.layers(),
        explicit_opening_layers(&committer_key, &evaluations, &point_blocks, &dimensions)
    );
    assert_eq!(
        proof.final_evaluations,
        explicit_final_evaluations(&evaluations, &point_blocks, &dimensions)
    );

    // Independently derive every partially evaluated table through ark-poly's
    // evaluator. This checks the cached and direct proof paths against the
    // exact evaluation-table convention used by other arkworks components.
    let mut arkworks_expected_layers = Vec::new();
    let mut staged = polynomial.polynomial().clone();
    let mut fixed_low_num_vars = 0usize;
    for (level, &dimension) in dimensions.iter().take(K - 1).enumerate() {
        let one_shot = polynomial
            .polynomial()
            .fix_variables(&point[..fixed_low_num_vars]);
        assert_eq!(staged, one_shot);

        let suffix_dimension = committer_key.h[level + 1].len();
        assert_eq!(staged.evaluations.len(), dimension * suffix_dimension);
        let projective = (0..dimension)
            .map(|current_assignment| {
                let coefficients = (0..suffix_dimension)
                    .map(|suffix_assignment| {
                        staged.evaluations[current_assignment + dimension * suffix_assignment]
                    })
                    .collect::<Vec<_>>();
                <G1Projective as VariableBaseMSM>::msm(&committer_key.h[level + 1], &coefficients)
                    .unwrap()
            })
            .collect::<Vec<_>>();
        arkworks_expected_layers.push(G1Projective::normalize_batch(&projective));

        let block_len = committer_key.num_vars_per_block[level];
        staged = staged.fix_variables(&point[fixed_low_num_vars..fixed_low_num_vars + block_len]);
        fixed_low_num_vars += block_len;
    }
    assert_eq!(proof.layers(), arkworks_expected_layers);
    assert_eq!(
        staged,
        polynomial
            .polynomial()
            .fix_variables(&point[..fixed_low_num_vars])
    );
    assert_eq!(proof.final_evaluations, staged.evaluations);
    assert_eq!(fixed_low_num_vars, NUM_VARS - point_blocks[K - 1].len());

    let expected_value = structured_evaluation(&point, 37);
    assert!(TestKZH::<K>::check(
        &verifier_key,
        &commitments,
        &point,
        [expected_value],
        &proof,
        &mut base_sponge.clone(),
        None,
    )
    .unwrap());

    // Exact proof equality alone would also pass if all layers were recomputed
    // directly. Perturb every stored table without changing its shape and
    // require the corresponding proof layer to change. The omitted third
    // table is covered above by the independent and Arkworks-derived
    // direct-path layer.
    for level in 0..states[0].auxiliary_tables.len() {
        let mut perturbed_state = states[0].clone();
        perturbed_state.auxiliary_tables[level][0] =
            add_generator(perturbed_state.auxiliary_tables[level][0]);
        let perturbed_proof = TestKZH::<K>::open(
            &committer_key,
            [&polynomial],
            &commitments,
            &point,
            &mut base_sponge.clone(),
            [&perturbed_state],
            None,
        )
        .unwrap();

        assert_ne!(perturbed_proof.layers()[level], proof.layers()[level]);
        assert!(!TestKZH::<K>::check(
            &verifier_key,
            &commitments,
            &point,
            [expected_value],
            &perturbed_proof,
            &mut base_sponge.clone(),
            None,
        )
        .unwrap());
    }
}

#[test]
fn commitment_state_omits_auxiliary_tables_unused_by_generic_opening() {
    const NUM_VARS: usize = 6;
    let (_, committer_key, _) = setup_and_trim::<3>(NUM_VARS, 13);
    let polynomial = structured_polynomial("aux-crossover", NUM_VARS, 38);
    let (_, states) = TestKZH::<3>::commit(&committer_key, [&polynomial], None).unwrap();

    let block_sizes = balanced_block_sizes(NUM_VARS, 3).unwrap();
    let dimensions = block_dimensions(&block_sizes).unwrap();
    let plan = opening_work_plan(&dimensions, 1).unwrap();
    let cached_count = plan
        .iter()
        .take_while(|layer| layer.source == OpeningLayerSource::Cached)
        .count();

    // Only the first table contracted by a generic opening is retained. The
    // second transition is direct, so no Boolean-selection-only table exists.
    assert_eq!(committer_key.num_vars_per_block, [2, 2, 2]);
    assert_eq!(block_sizes, [2, 2, 2]);
    assert_eq!(states[0].auxiliary_tables.len(), cached_count);
    assert_eq!(
        states[0]
            .auxiliary_tables
            .iter()
            .map(Vec::len)
            .collect::<Vec<_>>(),
        [4]
    );
    assert_eq!(
        plan.iter().map(|layer| layer.source).collect::<Vec<_>>(),
        [OpeningLayerSource::Cached, OpeningLayerSource::Direct,]
    );
}

// Same-point batching and prepared verification.

#[test]
fn same_point_opening_batches_multiple_polynomials_into_one_proof() {
    const NUM_VARS: usize = 8;
    let (_, committer_key, verifier_key) = setup_and_trim::<4>(NUM_VARS, 21);
    let polynomials = vec![
        structured_polynomial("first", NUM_VARS, 1),
        structured_polynomial("second", NUM_VARS, 9),
        structured_polynomial("third", NUM_VARS, 27),
    ];
    let point = non_boolean_point(NUM_VARS, 41);
    let honest_values: Vec<_> = [1u64, 9, 27]
        .iter()
        .map(|salt| structured_evaluation(&point, *salt))
        .collect();
    for (polynomial, expected) in polynomials.iter().zip(&honest_values) {
        assert_eq!(polynomial.evaluate(&point), *expected);
    }

    let (commitments, states) = TestKZH::<4>::commit(&committer_key, &polynomials, None).unwrap();
    assert_eq!(
        states
            .iter()
            .map(|state| state.auxiliary_tables.len())
            .collect::<Vec<_>>(),
        [2, 2, 2]
    );
    let base_sponge = poseidon_sponge_for_test::<Fr>();
    let proof = TestKZH::<4>::open(
        &committer_key,
        &polynomials,
        &commitments,
        &point,
        &mut base_sponge.clone(),
        &states,
        None,
    )
    .unwrap();

    // The proof has one transition layer per adjacent tensor-block pair, not
    // one proof per polynomial.
    assert_eq!(proof.layers().len(), 3);
    assert!(TestKZH::<4>::check(
        &verifier_key,
        &commitments,
        &point,
        honest_values.clone(),
        &proof,
        &mut base_sponge.clone(),
        None,
    )
    .unwrap());

    // The second cached layer must survive challenge aggregation and be used
    // by the generic opening path. The first polynomial has challenge one, so
    // this perturbation cannot disappear through a zero batching challenge.
    let mut perturbed_states = states.clone();
    perturbed_states[0].auxiliary_tables[1][0] =
        add_generator(perturbed_states[0].auxiliary_tables[1][0]);
    let perturbed_proof = TestKZH::<4>::open(
        &committer_key,
        &polynomials,
        &commitments,
        &point,
        &mut base_sponge.clone(),
        &perturbed_states,
        None,
    )
    .unwrap();
    assert_ne!(perturbed_proof.layers()[1], proof.layers()[1]);
    assert!(!TestKZH::<4>::check(
        &verifier_key,
        &commitments,
        &point,
        honest_values.clone(),
        &perturbed_proof,
        &mut base_sponge.clone(),
        None,
    )
    .unwrap());

    let mut incorrect_values = honest_values.clone();
    incorrect_values[1] += Fr::one();
    assert!(!TestKZH::<4>::check(
        &verifier_key,
        &commitments,
        &point,
        incorrect_values,
        &proof,
        &mut base_sponge.clone(),
        None,
    )
    .unwrap());

    // Preserve the honest aggregate under the *honest* challenge. This is the
    // cancellation attack that succeeds if batching challenges do not bind
    // the individual claimed values. Verification must derive a different
    // challenge from these false claims and reject them.
    let commitment_refs: Vec<_> = commitments.iter().collect();
    let honest_challenges = TestKZH::<4>::batch_challenges(
        &mut base_sponge.clone(),
        &commitment_refs,
        &point,
        &honest_values,
        NUM_VARS,
    )
    .unwrap();
    let mut cancelling_values = honest_values;
    let delta = Fr::from(123u64);
    cancelling_values[0] += delta;
    cancelling_values[1] -= delta * honest_challenges[1].inverse().unwrap();
    assert!(!TestKZH::<4>::check(
        &verifier_key,
        &commitments,
        &point,
        cancelling_values,
        &proof,
        &mut base_sponge.clone(),
        None,
    )
    .unwrap());
}

fn assert_prepared_key_matches_verifier<const K: usize>(
    prepared_verifier_key: &KZHPreparedVerifierKey<Bls12_381, K>,
    verifier_key: &KZHVerifierKey<Bls12_381, K>,
) {
    assert_eq!(prepared_verifier_key.num_vars(), verifier_key.num_vars());
    assert_eq!(
        prepared_verifier_key.num_vars_per_block(),
        verifier_key.num_vars_per_block()
    );
    assert_eq!(
        prepared_verifier_key.prepared_v(),
        &<Bls12_381 as Pairing>::G2Prepared::from(&verifier_key.v)
    );
    assert_eq!(
        prepared_verifier_key.prepared_v_tau().len(),
        verifier_key.v_tau.len()
    );
    for (prepared_layer, affine_layer) in prepared_verifier_key
        .prepared_v_tau()
        .iter()
        .zip(&verifier_key.v_tau)
    {
        assert_eq!(prepared_layer.len(), affine_layer.len());
        for (prepared, affine) in prepared_layer.iter().zip(affine_layer) {
            assert_eq!(prepared, &<Bls12_381 as Pairing>::G2Prepared::from(affine));
        }
    }
}

fn assert_prepared_verification_matches_standard<const K: usize>(num_vars: usize, seed: u8) {
    let (_, committer_key, verifier_key) = setup_and_trim::<K>(num_vars, seed);
    let prepared_verifier_key = KZHPreparedVerifierKey::prepare(&verifier_key);
    assert_prepared_key_matches_verifier(&prepared_verifier_key, &verifier_key);

    // Use two polynomials so this also checks that the prepared path derives
    // exactly the same batching challenges and aggregate commitment.
    let polynomials = vec![
        structured_polynomial("prepared-first", num_vars, seed as u64 + 3),
        structured_polynomial("prepared-second", num_vars, seed as u64 + 17),
    ];
    let point = non_boolean_point(num_vars, seed as u64 + 29);
    let values = polynomials
        .iter()
        .map(|polynomial| polynomial.evaluate(&point))
        .collect::<Vec<_>>();
    let (commitments, states) = TestKZH::<K>::commit(&committer_key, &polynomials, None).unwrap();
    let base_sponge = poseidon_sponge_for_test::<Fr>();
    let proof = TestKZH::<K>::open(
        &committer_key,
        &polynomials,
        &commitments,
        &point,
        &mut base_sponge.clone(),
        &states,
        None,
    )
    .unwrap();

    let mut standard_sponge = base_sponge.clone();
    let mut prepared_sponge = base_sponge.clone();
    assert!(TestKZH::<K>::check(
        &verifier_key,
        &commitments,
        &point,
        values.clone(),
        &proof,
        &mut standard_sponge,
        None,
    )
    .unwrap());
    assert!(TestKZH::<K>::check_prepared(
        &prepared_verifier_key,
        &commitments,
        &point,
        values.clone(),
        &proof,
        &mut prepared_sponge,
        None,
    )
    .unwrap());
    let standard_tail = standard_sponge.squeeze_field_elements::<Fr>(3);
    assert_eq!(
        prepared_sponge.squeeze_field_elements::<Fr>(3),
        standard_tail
    );

    // A prepared key is reusable. Its Clone implementation must produce an
    // independent cache with identical verification behavior.
    for reusable_key in [&prepared_verifier_key, &prepared_verifier_key.clone()] {
        let mut reusable_sponge = base_sponge.clone();
        assert!(TestKZH::<K>::check_prepared(
            reusable_key,
            &commitments,
            &point,
            values.clone(),
            &proof,
            &mut reusable_sponge,
            None,
        )
        .unwrap());
        assert_eq!(
            reusable_sponge.squeeze_field_elements::<Fr>(3),
            standard_tail
        );
    }
}

#[test]
fn prepared_check_matches_standard_for_balanced_and_uneven_families() {
    assert_prepared_verification_matches_standard::<2>(5, 23);
    assert_prepared_verification_matches_standard::<3>(7, 24);
    assert_prepared_verification_matches_standard::<4>(6, 25);
    assert_prepared_verification_matches_standard::<6>(7, 26);
}

#[test]
fn prepared_check_has_the_same_rejection_behavior() {
    const NUM_VARS: usize = 7;
    const K: usize = 6;
    let (_, committer_key, verifier_key) = setup_and_trim::<K>(NUM_VARS, 27);
    let prepared_verifier_key = KZHPreparedVerifierKey::prepare(&verifier_key);
    let polynomial = structured_polynomial("prepared-rejection", NUM_VARS, 31);
    let point = non_boolean_point(NUM_VARS, 43);
    let value = polynomial.evaluate(&point);
    let (commitments, states) = TestKZH::<K>::commit(&committer_key, [&polynomial], None).unwrap();
    let base_sponge = poseidon_sponge_for_test::<Fr>();
    let proof = TestKZH::<K>::open(
        &committer_key,
        [&polynomial],
        &commitments,
        &point,
        &mut base_sponge.clone(),
        &states,
        None,
    )
    .unwrap();

    let assert_both_reject = |candidate: &_, candidate_value| {
        let mut standard_sponge = base_sponge.clone();
        let mut prepared_sponge = base_sponge.clone();
        assert!(!TestKZH::<K>::check(
            &verifier_key,
            &commitments,
            &point,
            [candidate_value],
            candidate,
            &mut standard_sponge,
            None,
        )
        .unwrap());
        assert!(!TestKZH::<K>::check_prepared(
            &prepared_verifier_key,
            &commitments,
            &point,
            [candidate_value],
            candidate,
            &mut prepared_sponge,
            None,
        )
        .unwrap());
        assert_eq!(
            prepared_sponge.squeeze_field_elements::<Fr>(3),
            standard_sponge.squeeze_field_elements::<Fr>(3)
        );
    };

    assert_both_reject(&proof, value + Fr::one());

    let mut tampered_final = proof.clone();
    tampered_final.final_evaluations[0] += Fr::one();
    assert_both_reject(&tampered_final, value);

    for level in 0..proof.layers().len() {
        let mut tampered_layer = proof.clone();
        tampered_layer.layer_commitments[level][0] = G1Affine::identity();
        assert_both_reject(&tampered_layer, value);
    }

    // Malformed inputs are rejected before either path mutates its sponge.
    let mut malformed_layer = proof;
    malformed_layer.layer_commitments[0].pop();
    let mut prepared_sponge = base_sponge.clone();
    let mut untouched_sponge = base_sponge.clone();
    assert!(matches!(
        TestKZH::<K>::check_prepared(
            &prepared_verifier_key,
            &commitments,
            &point,
            [value],
            &malformed_layer,
            &mut prepared_sponge,
            None,
        ),
        Err(Error::IncorrectCommitmentSize { .. })
    ));
    assert_eq!(
        prepared_sponge.squeeze_field_elements::<Fr>(3),
        untouched_sponge.squeeze_field_elements::<Fr>(3)
    );
}

#[test]
fn four_polynomial_batch_uses_the_batch_aware_auxiliary_crossover() {
    const NUM_VARS: usize = 8;
    let (_, committer_key, verifier_key) = setup_and_trim::<4>(NUM_VARS, 22);
    let salts = [2u64, 5, 11, 23];
    let polynomials: Vec<_> = salts
        .iter()
        .enumerate()
        .map(|(index, salt)| structured_polynomial(&format!("batch-four-{index}"), NUM_VARS, *salt))
        .collect();
    let point = non_boolean_point(NUM_VARS, 47);
    let expected_values: Vec<_> = salts
        .iter()
        .map(|salt| structured_evaluation(&point, *salt))
        .collect();
    let (commitments, states) = TestKZH::<4>::commit(&committer_key, &polynomials, None).unwrap();
    assert!(states.iter().all(|state| state.auxiliary_tables.len() == 2));

    let base_sponge = poseidon_sponge_for_test::<Fr>();
    let proof = TestKZH::<4>::open(
        &committer_key,
        &polynomials,
        &commitments,
        &point,
        &mut base_sponge.clone(),
        &states,
        None,
    )
    .unwrap();
    assert!(TestKZH::<4>::check(
        &verifier_key,
        &commitments,
        &point,
        expected_values.clone(),
        &proof,
        &mut base_sponge.clone(),
        None,
    )
    .unwrap());

    // At level one, four prefix contractions cost 4 * 16 terms, exactly the
    // 64 terms required by the direct route. The strict crossover must ignore
    // this otherwise-valid cached table.
    let mut ignored_second_layer = states.clone();
    ignored_second_layer[0].auxiliary_tables[1][0] =
        add_generator(ignored_second_layer[0].auxiliary_tables[1][0]);
    let direct_proof = TestKZH::<4>::open(
        &committer_key,
        &polynomials,
        &commitments,
        &point,
        &mut base_sponge.clone(),
        &ignored_second_layer,
        None,
    )
    .unwrap();
    assert_eq!(direct_proof, proof);
    assert!(TestKZH::<4>::check(
        &verifier_key,
        &commitments,
        &point,
        expected_values.clone(),
        &direct_proof,
        &mut base_sponge.clone(),
        None,
    )
    .unwrap());

    // Level zero remains strictly beneficial: 4 * 4 < 256. Perturbing it must
    // affect the proof and invalidate it against the original commitment.
    let mut used_first_layer = states.clone();
    used_first_layer[0].auxiliary_tables[0][0] =
        add_generator(used_first_layer[0].auxiliary_tables[0][0]);
    let perturbed_proof = TestKZH::<4>::open(
        &committer_key,
        &polynomials,
        &commitments,
        &point,
        &mut base_sponge.clone(),
        &used_first_layer,
        None,
    )
    .unwrap();
    assert_ne!(perturbed_proof.layers()[0], proof.layers()[0]);
    assert!(!TestKZH::<4>::check(
        &verifier_key,
        &commitments,
        &point,
        expected_values,
        &perturbed_proof,
        &mut base_sponge.clone(),
        None,
    )
    .unwrap());
}

// Default multi-point and linear-combination PCS APIs.

#[test]
fn default_batch_api_opens_multiple_polynomials_at_multiple_points() {
    const NUM_VARS: usize = 6;
    let (_, committer_key, verifier_key) = setup_and_trim::<4>(NUM_VARS, 31);
    let polynomials = vec![
        structured_polynomial("alpha", NUM_VARS, 4),
        structured_polynomial("beta", NUM_VARS, 15),
    ];
    let (commitments, states) = TestKZH::<4>::commit(&committer_key, &polynomials, None).unwrap();

    let points = [
        non_boolean_point(NUM_VARS, 7),
        non_boolean_point(NUM_VARS, 29),
        non_boolean_point(NUM_VARS, 53),
    ];
    let mut query_set = QuerySet::new();
    let mut evaluations = Evaluations::new();
    for (point_index, point) in points.iter().enumerate() {
        for polynomial in &polynomials {
            query_set.insert((
                polynomial.label().clone(),
                (format!("point-{point_index}"), point.clone()),
            ));
            evaluations.insert(
                (polynomial.label().clone(), point.clone()),
                polynomial.evaluate(point),
            );
        }
    }

    let base_sponge = poseidon_sponge_for_test::<Fr>();
    let proof = TestKZH::<4>::batch_open(
        &committer_key,
        &polynomials,
        &commitments,
        &query_set,
        &mut base_sponge.clone(),
        &states,
        None,
    )
    .unwrap();
    assert_eq!(proof.len(), points.len());

    let mut verifier_rng = seeded_rng(32);
    assert!(TestKZH::<4>::batch_check(
        &verifier_key,
        &commitments,
        &query_set,
        &evaluations,
        &proof,
        &mut base_sponge.clone(),
        &mut verifier_rng,
    )
    .unwrap());

    let first_key = (polynomials[0].label().clone(), points[0].clone());
    let mut incorrect_evaluations = evaluations;
    *incorrect_evaluations.get_mut(&first_key).unwrap() += Fr::one();
    let mut verifier_rng = seeded_rng(33);
    assert!(!TestKZH::<4>::batch_check(
        &verifier_key,
        &commitments,
        &query_set,
        &incorrect_evaluations,
        &proof,
        &mut base_sponge.clone(),
        &mut verifier_rng,
    )
    .unwrap());
}

#[test]
fn default_linear_combination_api_opens_multiple_equations_at_multiple_points() {
    const NUM_VARS: usize = 6;
    let (_, committer_key, verifier_key) = setup_and_trim::<4>(NUM_VARS, 35);
    let polynomials = vec![
        structured_polynomial("alpha", NUM_VARS, 6),
        structured_polynomial("beta", NUM_VARS, 19),
    ];
    let (commitments, states) = TestKZH::<4>::commit(&committer_key, &polynomials, None).unwrap();

    let linear_combinations = vec![
        LinearCombination::new(
            "weighted-sum",
            vec![
                (
                    Fr::from(2u64),
                    LCTerm::PolyLabel(polynomials[0].label().clone()),
                ),
                (
                    Fr::from(3u64),
                    LCTerm::PolyLabel(polynomials[1].label().clone()),
                ),
                (Fr::from(5u64), LCTerm::One),
            ],
        ),
        LinearCombination::new(
            "affine-difference",
            vec![
                (
                    Fr::from(7u64),
                    LCTerm::PolyLabel(polynomials[0].label().clone()),
                ),
                (
                    -Fr::one(),
                    LCTerm::PolyLabel(polynomials[1].label().clone()),
                ),
                (Fr::from(11u64), LCTerm::One),
            ],
        ),
    ];
    let points = [
        non_boolean_point(NUM_VARS, 13),
        non_boolean_point(NUM_VARS, 41),
    ];

    let mut query_set = QuerySet::new();
    let mut evaluations = Evaluations::new();
    for (point_index, point) in points.iter().enumerate() {
        let alpha = polynomials[0].evaluate(point);
        let beta = polynomials[1].evaluate(point);
        let expected = [
            Fr::from(2u64) * alpha + Fr::from(3u64) * beta + Fr::from(5u64),
            Fr::from(7u64) * alpha - beta + Fr::from(11u64),
        ];

        for (linear_combination, value) in linear_combinations.iter().zip(expected) {
            query_set.insert((
                linear_combination.label().clone(),
                (format!("lc-point-{point_index}"), point.clone()),
            ));
            evaluations.insert((linear_combination.label().clone(), point.clone()), value);
        }
    }

    let base_sponge = poseidon_sponge_for_test::<Fr>();
    let proof = TestKZH::<4>::open_combinations(
        &committer_key,
        &linear_combinations,
        &polynomials,
        &commitments,
        &query_set,
        &mut base_sponge.clone(),
        &states,
        None,
    )
    .unwrap();
    assert_eq!(proof.proof.len(), points.len());
    assert_eq!(
        proof.evals.as_ref().map(Vec::len),
        Some(polynomials.len() * points.len())
    );

    let mut verifier_rng = seeded_rng(36);
    assert!(TestKZH::<4>::check_combinations(
        &verifier_key,
        &linear_combinations,
        &commitments,
        &query_set,
        &evaluations,
        &proof,
        &mut base_sponge.clone(),
        &mut verifier_rng,
    )
    .unwrap());

    let first_key = (linear_combinations[0].label().clone(), points[0].clone());
    let mut incorrect_evaluations = evaluations;
    *incorrect_evaluations.get_mut(&first_key).unwrap() += Fr::one();
    let mut verifier_rng = seeded_rng(37);
    assert!(!TestKZH::<4>::check_combinations(
        &verifier_key,
        &linear_combinations,
        &commitments,
        &query_set,
        &incorrect_evaluations,
        &proof,
        &mut base_sponge.clone(),
        &mut verifier_rng,
    )
    .unwrap());
}

// Malformed-input and negative security tests.

#[test]
fn tampered_claims_commitments_and_proof_elements_are_rejected() {
    const NUM_VARS: usize = 5;
    const K: usize = 3;
    let (_, committer_key, verifier_key) = setup_and_trim::<K>(NUM_VARS, 41);
    let polynomial = structured_polynomial("tamper", NUM_VARS, 8);
    let point = non_boolean_point(NUM_VARS, 19);
    let value = structured_evaluation(&point, 8);
    let (commitments, states) = TestKZH::<K>::commit(&committer_key, [&polynomial], None).unwrap();
    // With dimensions [4,4,2], only the first cached prefix is strictly
    // smaller than the corresponding partially evaluated tensor.
    assert_eq!(states[0].auxiliary_tables.len(), 1);
    assert_eq!(states[0].auxiliary_tables[0].len(), 4);
    let base_sponge = poseidon_sponge_for_test::<Fr>();
    let proof = TestKZH::<K>::open(
        &committer_key,
        [&polynomial],
        &commitments,
        &point,
        &mut base_sponge.clone(),
        &states,
        None,
    )
    .unwrap();

    assert!(!TestKZH::<K>::check(
        &verifier_key,
        &commitments,
        &point,
        [value + Fr::one()],
        &proof,
        &mut base_sponge.clone(),
        None,
    )
    .unwrap());

    let wrong_point = non_boolean_point(NUM_VARS, 20);
    assert!(!TestKZH::<K>::check(
        &verifier_key,
        &commitments,
        &wrong_point,
        [polynomial.evaluate(&wrong_point)],
        &proof,
        &mut base_sponge.clone(),
        None,
    )
    .unwrap());

    let mut wrong_commitment = *commitments[0].commitment();
    wrong_commitment.comm = G1Affine::identity();
    let wrong_commitment =
        LabeledCommitment::new(commitments[0].label().clone(), wrong_commitment, None);
    assert!(!TestKZH::<K>::check(
        &verifier_key,
        [&wrong_commitment],
        &point,
        [value],
        &proof,
        &mut base_sponge.clone(),
        None,
    )
    .unwrap());

    let mut tampered_scalar = proof.clone();
    tampered_scalar.final_evaluations[0] += Fr::one();
    assert!(!TestKZH::<K>::check(
        &verifier_key,
        &commitments,
        &point,
        [value],
        &tampered_scalar,
        &mut base_sponge.clone(),
        None,
    )
    .unwrap());

    for level in 0..proof.layers().len() {
        let mut tampered_group = proof.clone();
        tampered_group.layer_commitments[level][0] = G1Affine::identity();
        assert!(!TestKZH::<K>::check(
            &verifier_key,
            &commitments,
            &point,
            [value],
            &tampered_group,
            &mut base_sponge.clone(),
            None,
        )
        .unwrap());
    }
}

#[test]
fn malformed_proof_shapes_are_rejected_before_transcript_changes() {
    const NUM_VARS: usize = 5;
    const K: usize = 3;
    let (_, committer_key, verifier_key) = setup_and_trim::<K>(NUM_VARS, 41);
    let polynomial = structured_polynomial("malformed-proof", NUM_VARS, 8);
    let point = non_boolean_point(NUM_VARS, 19);
    let value = structured_evaluation(&point, 8);
    let (commitments, states) = TestKZH::<K>::commit(&committer_key, [&polynomial], None).unwrap();
    let base_sponge = poseidon_sponge_for_test::<Fr>();
    let proof = TestKZH::<K>::open(
        &committer_key,
        [&polynomial],
        &commitments,
        &point,
        &mut base_sponge.clone(),
        &states,
        None,
    )
    .unwrap();

    let mut malformed_layer = proof.clone();
    malformed_layer.layer_commitments[0].pop();
    let mut malformed_proof_sponge = base_sponge.clone();
    let mut untouched_proof_sponge = base_sponge.clone();
    assert!(matches!(
        TestKZH::<K>::check(
            &verifier_key,
            &commitments,
            &point,
            [value],
            &malformed_layer,
            &mut malformed_proof_sponge,
            None,
        ),
        Err(Error::IncorrectCommitmentSize { .. })
    ));
    assert_eq!(
        malformed_proof_sponge.squeeze_field_elements::<Fr>(3),
        untouched_proof_sponge.squeeze_field_elements::<Fr>(3),
    );

    let mut missing_layer = proof.clone();
    missing_layer.layer_commitments.pop();
    assert!(matches!(
        TestKZH::<K>::check(
            &verifier_key,
            &commitments,
            &point,
            [value],
            &missing_layer,
            &mut base_sponge.clone(),
            None,
        ),
        Err(Error::InvalidParameters(_))
    ));

    let mut malformed_final_vector = proof.clone();
    malformed_final_vector.final_evaluations.pop();
    assert!(matches!(
        TestKZH::<K>::check(
            &verifier_key,
            &commitments,
            &point,
            [value],
            &malformed_final_vector,
            &mut base_sponge.clone(),
            None,
        ),
        Err(Error::IncorrectCommitmentSize { .. })
    ));

    let mut wrong_family = proof.clone();
    wrong_family.k = 2;
    assert!(matches!(
        TestKZH::<K>::check(
            &verifier_key,
            &commitments,
            &point,
            [value],
            &wrong_family,
            &mut base_sponge.clone(),
            None,
        ),
        Err(Error::InvalidParameters(_))
    ));
}

#[test]
fn opening_rejects_mismatched_inputs_and_malformed_state() {
    const NUM_VARS: usize = 5;
    const K: usize = 3;
    let (_, committer_key, verifier_key) = setup_and_trim::<K>(NUM_VARS, 41);
    let polynomial = structured_polynomial("malformed-opening", NUM_VARS, 8);
    let point = non_boolean_point(NUM_VARS, 19);
    let value = structured_evaluation(&point, 8);
    let (commitments, states) = TestKZH::<K>::commit(&committer_key, [&polynomial], None).unwrap();
    let base_sponge = poseidon_sponge_for_test::<Fr>();
    let proof = TestKZH::<K>::open(
        &committer_key,
        [&polynomial],
        &commitments,
        &point,
        &mut base_sponge.clone(),
        &states,
        None,
    )
    .unwrap();

    let mismatched_label = LabeledCommitment::new(
        "different-label".to_string(),
        *commitments[0].commitment(),
        None,
    );
    assert!(matches!(
        TestKZH::<K>::open(
            &committer_key,
            [&polynomial],
            [&mismatched_label],
            &point,
            &mut base_sponge.clone(),
            &states,
            None,
        ),
        Err(Error::MismatchedLabels { .. })
    ));

    assert!(matches!(
        TestKZH::<K>::open(
            &committer_key,
            [&polynomial],
            &commitments[..0],
            &point,
            &mut base_sponge.clone(),
            &states,
            None,
        ),
        Err(Error::IncorrectInputLength(_))
    ));

    let malformed_table = LabeledPolynomial::new(
        polynomial.label().clone(),
        TestPolynomial {
            evaluations: vec![Fr::zero(); (1usize << NUM_VARS) - 1],
            num_vars: NUM_VARS,
        },
        None,
        None,
    );
    let mut malformed_table_sponge = base_sponge.clone();
    let mut untouched_table_sponge = base_sponge.clone();
    assert!(matches!(
        TestKZH::<K>::open(
            &committer_key,
            [&malformed_table],
            &commitments,
            &point,
            &mut malformed_table_sponge,
            &states,
            None,
        ),
        Err(Error::IncorrectInputLength(_))
    ));
    assert_eq!(
        malformed_table_sponge.squeeze_field_elements::<Fr>(3),
        untouched_table_sponge.squeeze_field_elements::<Fr>(3),
    );
    assert!(matches!(
        TestKZH::<K>::open(
            &committer_key,
            [&polynomial],
            &commitments,
            &point,
            &mut base_sponge.clone(),
            &states[..0],
            None,
        ),
        Err(Error::IncorrectInputLength(_))
    ));
    assert!(matches!(
        TestKZH::<K>::check(
            &verifier_key,
            &commitments,
            &point,
            Vec::<Fr>::new(),
            &proof,
            &mut base_sponge.clone(),
            None,
        ),
        Err(Error::IncorrectInputLength(_))
    ));

    for malformed_point in [point[..NUM_VARS - 1].to_vec(), {
        let mut point = point.clone();
        point.push(Fr::from(101u64));
        point
    }] {
        assert!(matches!(
            TestKZH::<K>::open(
                &committer_key,
                [&polynomial],
                &commitments,
                &malformed_point,
                &mut base_sponge.clone(),
                &states,
                None,
            ),
            Err(Error::MismatchedNumVars { .. })
        ));
        assert!(matches!(
            TestKZH::<K>::check(
                &verifier_key,
                &commitments,
                &malformed_point,
                [value],
                &proof,
                &mut base_sponge.clone(),
                None,
            ),
            Err(Error::MismatchedNumVars { .. })
        ));
    }

    let mut malformed_state = states[0].clone();
    malformed_state.k = 2;
    assert!(matches!(
        TestKZH::<K>::open(
            &committer_key,
            [&polynomial],
            &commitments,
            &point,
            &mut base_sponge.clone(),
            [&malformed_state],
            None,
        ),
        Err(Error::InvalidParameters(_))
    ));

    let mut malformed_auxiliary_state = states[0].clone();
    malformed_auxiliary_state.auxiliary_tables[0].pop();
    assert!(matches!(
        TestKZH::<K>::open(
            &committer_key,
            [&polynomial],
            &commitments,
            &point,
            &mut base_sponge.clone(),
            [&malformed_auxiliary_state],
            None,
        ),
        Err(Error::IncorrectCommitmentSize { .. })
    ));

    let mut wrong_size_state = states[0].clone();
    wrong_size_state.num_vars -= 1;
    assert!(TestKZH::<K>::open(
        &committer_key,
        [&polynomial],
        &commitments,
        &point,
        &mut base_sponge.clone(),
        [&wrong_size_state],
        None,
    )
    .is_err());
}

#[test]
fn malformed_keys_are_rejected_before_transcript_changes() {
    const NUM_VARS: usize = 5;
    const K: usize = 3;
    let (_, committer_key, verifier_key) = setup_and_trim::<K>(NUM_VARS, 41);
    let polynomial = structured_polynomial("malformed-keys", NUM_VARS, 8);
    let point = non_boolean_point(NUM_VARS, 19);
    let value = structured_evaluation(&point, 8);
    let (commitments, states) = TestKZH::<K>::commit(&committer_key, [&polynomial], None).unwrap();
    let base_sponge = poseidon_sponge_for_test::<Fr>();
    let proof = TestKZH::<K>::open(
        &committer_key,
        [&polynomial],
        &commitments,
        &point,
        &mut base_sponge.clone(),
        &states,
        None,
    )
    .unwrap();

    let mut malformed_committer_key = committer_key.clone();
    malformed_committer_key.h[1].pop();
    assert!(matches!(
        TestKZH::<K>::commit(&malformed_committer_key, [&polynomial], None),
        Err(Error::IncorrectCommitmentSize { .. })
    ));

    let mut malformed_verifier_key = verifier_key.clone();
    malformed_verifier_key.v_tau[0].pop();
    assert!(matches!(
        TestKZH::<K>::check(
            &malformed_verifier_key,
            &commitments,
            &point,
            [value],
            &proof,
            &mut base_sponge.clone(),
            None,
        ),
        Err(Error::IncorrectCommitmentSize { .. })
    ));

    let malformed_prepared_key = KZHPreparedVerifierKey::prepare(&malformed_verifier_key);
    let mut malformed_key_sponge = base_sponge.clone();
    let mut untouched_key_sponge = base_sponge.clone();
    assert!(matches!(
        TestKZH::<K>::check_prepared(
            &malformed_prepared_key,
            &commitments,
            &point,
            [value],
            &proof,
            &mut malformed_key_sponge,
            None,
        ),
        Err(Error::IncorrectCommitmentSize { .. })
    ));
    assert_eq!(
        malformed_key_sponge.squeeze_field_elements::<Fr>(3),
        untouched_key_sponge.squeeze_field_elements::<Fr>(3)
    );
}

#[test]
fn hiding_and_strict_degree_bounds_are_rejected() {
    const NUM_VARS: usize = 4;
    let mut rng = seeded_rng(51);
    let params = TestKZH::<2>::setup(1, Some(NUM_VARS), &mut rng).unwrap();

    assert!(matches!(
        TestKZH::<2>::trim(&params, 1, 1, None),
        Err(Error::InvalidParameters(_))
    ));
    assert!(matches!(
        TestKZH::<2>::trim(&params, 1, 0, Some(&[1])),
        Err(Error::UnsupportedDegreeBound(1))
    ));
    assert!(matches!(
        TestKZH::<2>::trim(&params, 1, 0, Some(&[])),
        Err(Error::EmptyDegreeBounds)
    ));

    let (committer_key, _) = TestKZH::<2>::trim(&params, 1, 0, None).unwrap();
    let plain = structured_polynomial("plain", NUM_VARS, 5);
    let (plain_commitments, plain_states) =
        TestKZH::<2>::commit(&committer_key, [&plain], None).unwrap();
    assert!(!plain_commitments[0].commitment().has_degree_bound());

    let hiding = LabeledPolynomial::new(
        "hiding".to_string(),
        plain.polynomial().clone(),
        None,
        Some(1),
    );
    assert!(matches!(
        TestKZH::<2>::commit(&committer_key, [&hiding], Some(&mut rng)),
        Err(Error::InvalidParameters(_))
    ));

    let bounded = LabeledPolynomial::new(
        "bounded".to_string(),
        plain.polynomial().clone(),
        Some(1),
        None,
    );
    assert!(matches!(
        TestKZH::<2>::commit(&committer_key, [&bounded], None),
        Err(Error::UnsupportedDegreeBound(1))
    ));

    let bounded_commitment = LabeledCommitment::new(
        plain.label().clone(),
        *plain_commitments[0].commitment(),
        Some(1),
    );
    let point = non_boolean_point(NUM_VARS, 17);
    let base_sponge = poseidon_sponge_for_test::<Fr>();
    let proof = TestKZH::<2>::open(
        &committer_key,
        [&plain],
        &plain_commitments,
        &point,
        &mut base_sponge.clone(),
        &plain_states,
        None,
    )
    .unwrap();
    assert!(matches!(
        TestKZH::<2>::open(
            &committer_key,
            [&plain],
            [&bounded_commitment],
            &point,
            &mut base_sponge.clone(),
            &plain_states,
            None,
        ),
        Err(Error::UnsupportedDegreeBound(1))
    ));
    let (_, verifier_key) = TestKZH::<2>::trim(&params, 1, 0, None).unwrap();
    assert!(matches!(
        TestKZH::<2>::check(
            &verifier_key,
            [&bounded_commitment],
            &point,
            [plain.evaluate(&point)],
            &proof,
            &mut base_sponge.clone(),
            None,
        ),
        Err(Error::UnsupportedDegreeBound(1))
    ));
}

// Canonical serialization and setup-boundary behavior.

fn canonical_round_trip<T>(value: &T) -> T
where
    T: CanonicalSerialize + CanonicalDeserialize + Debug + PartialEq,
{
    let mut compressed = Vec::new();
    value
        .serialize_with_mode(&mut compressed, Compress::Yes)
        .unwrap();
    let mut compressed_reader = compressed.as_slice();
    let decoded =
        T::deserialize_with_mode(&mut compressed_reader, Compress::Yes, Validate::Yes).unwrap();
    assert_eq!(decoded, *value);
    assert!(compressed_reader.is_empty());

    let mut uncompressed = Vec::new();
    value
        .serialize_with_mode(&mut uncompressed, Compress::No)
        .unwrap();
    let mut uncompressed_reader = uncompressed.as_slice();
    let uncompressed_decoded =
        T::deserialize_with_mode(&mut uncompressed_reader, Compress::No, Validate::Yes).unwrap();
    assert_eq!(uncompressed_decoded, *value);
    assert!(uncompressed_reader.is_empty());

    decoded
}

fn assert_family_serialization_round_trips<const K: usize>(num_vars: usize, seed: u8) {
    let (params, committer_key, verifier_key) = setup_and_trim::<K>(num_vars, seed);
    let polynomial = structured_polynomial("serialization", num_vars, seed as u64);
    let point = non_boolean_point(num_vars, seed as u64 + 9);
    let (commitments, states) = TestKZH::<K>::commit(&committer_key, [&polynomial], None).unwrap();
    let base_sponge = poseidon_sponge_for_test::<Fr>();
    let proof = TestKZH::<K>::open(
        &committer_key,
        [&polynomial],
        &commitments,
        &point,
        &mut base_sponge.clone(),
        &states,
        None,
    )
    .unwrap();

    let decoded_params = canonical_round_trip(&params);
    let decoded_committer_key = canonical_round_trip(&committer_key);
    let decoded_verifier_key = canonical_round_trip(&verifier_key);
    let decoded_commitment = canonical_round_trip(commitments[0].commitment());
    let decoded_state = canonical_round_trip(&states[0]);
    let decoded_proof = canonical_round_trip(&proof);
    let decoded_batch = canonical_round_trip(&vec![proof]);

    let (trimmed_committer_key, trimmed_verifier_key) =
        TestKZH::<K>::trim(&decoded_params, 1, 0, None).unwrap();
    assert_eq!(trimmed_committer_key, decoded_committer_key);
    assert_eq!(trimmed_verifier_key, decoded_verifier_key);

    let decoded_labeled_commitment =
        LabeledCommitment::new(polynomial.label().clone(), decoded_commitment, None);
    let value = structured_evaluation(&point, seed as u64);
    assert!(TestKZH::<K>::check(
        &decoded_verifier_key,
        [&decoded_labeled_commitment],
        &point,
        [value],
        &decoded_proof,
        &mut base_sponge.clone(),
        None,
    )
    .unwrap());
    let prepared_decoded_verifier_key = KZHPreparedVerifierKey::prepare(&decoded_verifier_key);
    assert!(TestKZH::<K>::check_prepared(
        &prepared_decoded_verifier_key,
        [&decoded_labeled_commitment],
        &point,
        [value],
        &decoded_proof,
        &mut base_sponge.clone(),
        None,
    )
    .unwrap());
    assert!(TestKZH::<K>::check(
        &decoded_verifier_key,
        [&decoded_labeled_commitment],
        &point,
        [value],
        &decoded_batch[0],
        &mut base_sponge.clone(),
        None,
    )
    .unwrap());

    let reopened = TestKZH::<K>::open(
        &decoded_committer_key,
        [&polynomial],
        [&decoded_labeled_commitment],
        &point,
        &mut base_sponge.clone(),
        [&decoded_state],
        None,
    )
    .unwrap();
    assert!(TestKZH::<K>::check(
        &decoded_verifier_key,
        [&decoded_labeled_commitment],
        &point,
        [value],
        &reopened,
        &mut base_sponge.clone(),
        None,
    )
    .unwrap());
}

#[test]
fn canonical_serialization_round_trips_for_k2_k3_and_k4() {
    assert_family_serialization_round_trips::<2>(3, 61);
    assert_family_serialization_round_trips::<3>(4, 62);
    assert_family_serialization_round_trips::<4>(5, 63);
}

#[test]
fn serialized_runtime_family_marker_is_checked_at_operation_boundaries() {
    let (params, _, _) = setup_and_trim::<2>(4, 71);
    let mut bytes = Vec::new();
    params.serialize_compressed(&mut bytes).unwrap();

    // Const generics do not alter the byte layout. Decoding under the wrong
    // Rust type therefore succeeds, but the serialized runtime marker and
    // block shape must cause the first KZH operation to reject the object.
    let mistyped =
        KZHUniversalParams::<Bls12_381, 3>::deserialize_compressed(bytes.as_slice()).unwrap();
    assert!(matches!(
        TestKZH::<3>::trim(&mistyped, 1, 0, None),
        Err(Error::InvalidParameters(_))
    ));
}

#[test]
fn setup_rejects_invalid_family_dimensions() {
    let mut rng = seeded_rng(81);
    assert!(matches!(
        TestKZH::<1>::setup(1, Some(4), &mut rng),
        Err(Error::InvalidNumberOfVariables)
    ));
    assert!(matches!(
        TestKZH::<4>::setup(1, Some(3), &mut rng),
        Err(Error::InvalidNumberOfVariables)
    ));
    assert!(matches!(
        TestKZH::<2>::setup(1, None, &mut rng),
        Err(Error::InvalidNumberOfVariables)
    ));
}

#[test]
fn setup_and_trim_enforce_multilinear_degree_one() {
    let mut rng = seeded_rng(82);
    assert!(matches!(
        TestKZH::<2>::setup(0, Some(4), &mut rng),
        Err(Error::DegreeIsZero)
    ));
    assert!(matches!(
        TestKZH::<2>::setup(2, Some(4), &mut rng),
        Err(Error::InvalidParameters(_))
    ));

    let params = TestKZH::<2>::setup(1, Some(4), &mut rng).unwrap();
    assert!(matches!(
        TestKZH::<2>::trim(&params, 0, 0, None),
        Err(Error::InvalidParameters(_))
    ));
    assert!(matches!(
        TestKZH::<2>::trim(&params, 2, 0, None),
        Err(Error::TrimmingDegreeTooLarge)
    ));
    assert!(TestKZH::<2>::trim(&params, 1, 0, None).is_ok());
}

#[test]
fn empty_commitment_is_a_usable_commitment_to_zero() {
    const NUM_VARS: usize = 4;
    let (_, committer_key, verifier_key) = setup_and_trim::<2>(NUM_VARS, 91);
    let zero_polynomial = LabeledPolynomial::new(
        "zero".to_string(),
        TestPolynomial::from_evaluations_vec(NUM_VARS, vec![Fr::from(0u64); 1usize << NUM_VARS]),
        None,
        None,
    );
    let (_, states) = TestKZH::<2>::commit(&committer_key, [&zero_polynomial], None).unwrap();
    let empty = LabeledCommitment::new(
        zero_polynomial.label().clone(),
        KZHCommitment::<Bls12_381, 2>::empty(),
        None,
    );
    let point = non_boolean_point(NUM_VARS, 23);
    let base_sponge = poseidon_sponge_for_test::<Fr>();
    let proof = TestKZH::<2>::open(
        &committer_key,
        [&zero_polynomial],
        [&empty],
        &point,
        &mut base_sponge.clone(),
        &states,
        None,
    )
    .unwrap();
    assert!(TestKZH::<2>::check(
        &verifier_key,
        [&empty],
        &point,
        [Fr::from(0u64)],
        &proof,
        &mut base_sponge.clone(),
        None,
    )
    .unwrap());
}
