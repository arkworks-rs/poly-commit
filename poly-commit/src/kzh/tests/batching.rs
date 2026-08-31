use super::*;

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
    prepared_verifier_key: &PreparedVerifierKey<Bls12_381, K>,
    verifier_key: &VerifierKey<Bls12_381, K>,
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
    let prepared_verifier_key = PreparedVerifierKey::prepare(&verifier_key);
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
    let prepared_verifier_key = PreparedVerifierKey::prepare(&verifier_key);
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
