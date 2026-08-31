use super::*;

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

    // The plan counts online group-scalar terms, not the dense tensor folding
    // that walks the evaluation table.

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
