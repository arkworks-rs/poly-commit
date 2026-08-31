use super::*;

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
