use ark_bls12_381::{Bls12_381, Fr};
use ark_poly::univariate::DensePolynomial;
use ark_poly::UVPolynomial;
use ark_std::test_rng;
use ark_std::vec::Vec;
use ark_std::UniformRand;

use crate::iterable::{Iterable, Reversed};
use crate::kzg::space::CommitterKeyStream;
use crate::kzg::time::CommitterKey;
use crate::kzg::VerifierKey;
use crate::misc::evaluate_le;

#[test]
fn test_commitment_consistency() {
    let rng = &mut ark_std::test_rng();
    let d = 15;
    let polynomial = DensePolynomial::<Fr>::rand(d, rng);
    let polynomial_stream = Reversed::new(polynomial.coeffs());
    let time_ck = CommitterKey::<Bls12_381>::new(d + 1, 3, rng);
    let space_ck = CommitterKeyStream::from(&time_ck);

    // compute the time commitment
    let time_commitment = time_ck.commit(&polynomial);
    let space_commitment = space_ck.commit(&polynomial_stream);

    assert_eq!(space_commitment, time_commitment);
}

#[test]
fn test_srs() {
    use ark_bls12_381::Bls12_381;

    let rng = &mut ark_std::test_rng();
    let time_ck = CommitterKey::<Bls12_381>::new(10, 3, rng);
    let space_ck = CommitterKeyStream::from(&time_ck);
    // Make sure that there are enough elements for the entire array.
    assert_eq!(time_ck.powers_of_g.len(), space_ck.powers_of_g.len());
}

#[test]
fn test_open_consistency() {
    let rng = &mut ark_std::test_rng();
    let d = 15;
    let max_msm_buffer = 1 << 20;
    let polynomials = DensePolynomial::<Fr>::rand(d, rng);
    let polynomial_stream = Reversed::new(polynomials.coeffs());
    let time_ck = CommitterKey::<Bls12_381>::new(d + 1, 3, rng);
    let space_ck = CommitterKeyStream::from(&time_ck);
    let alpha = Fr::rand(rng);

    // compute the time commitment
    let (time_evaluation, time_open) = time_ck.open(&polynomials, &alpha);
    let (space_evaluation, space_open) = space_ck.open(&polynomial_stream, &alpha, max_msm_buffer);
    // compute the space commitment
    assert_eq!(time_evaluation, space_evaluation);
    assert_eq!(time_open, space_open);
}

#[test]
fn test_open_multipoints_correctness() {
    let mut rng = &mut test_rng();
    let d = 100;

    let eval_points = (0..5).map(|_| Fr::rand(rng)).collect::<Vec<_>>();
    let polynomials = (0..15)
        .map(|_| DensePolynomial::<Fr>::rand(d, rng).coeffs)
        .collect::<Vec<_>>();
    let evals = polynomials
        .iter()
        .map(|p| {
            eval_points
                .iter()
                .map(|e| evaluate_le(p, e))
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();

    let time_ck = CommitterKey::<Bls12_381>::new(d + 1, eval_points.len(), rng);
    let time_vk = VerifierKey::from(&time_ck);

    let time_batched_commitments = time_ck.batch_commit(&polynomials);

    let eta: Fr = u128::rand(&mut rng).into();

    let proof = time_ck.batch_open_multi_points(
        &polynomials.iter().collect::<Vec<_>>()[..],
        &eval_points,
        &eta,
    );

    let verification_result = time_vk.verify_multi_points(
        &time_batched_commitments,
        &eval_points,
        &evals,
        &proof,
        &eta,
    );

    assert!(verification_result.is_ok());
}
