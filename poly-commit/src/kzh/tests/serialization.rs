use super::*;

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
    let prepared_decoded_verifier_key = PreparedVerifierKey::prepare(&decoded_verifier_key);
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
        UniversalParams::<Bls12_381, 3>::deserialize_compressed(bytes.as_slice()).unwrap();
    assert!(matches!(
        TestKZH::<3>::trim(&mistyped, 1, 0, None),
        Err(Error::InvalidParameters(_))
    ));
}
