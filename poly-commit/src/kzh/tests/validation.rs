use super::*;

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

    let malformed_prepared_key = PreparedVerifierKey::prepare(&malformed_verifier_key);
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
        Commitment::<Bls12_381, 2>::empty(),
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
