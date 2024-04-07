use ark_bls12_381::{Bls12_381, Fr};
use ark_poly::univariate::DensePolynomial;
use ark_poly::DenseUVPolynomial;
use ark_std::{UniformRand, Zero};

use crate::streaming_kzg::space::CommitterKeyStream;
use crate::streaming_kzg::time::CommitterKey;
use crate::streaming_kzg::{vanishing_polynomial, VerifierKey};
use ark_ff::Field;
use ark_std::borrow::Borrow;
use ark_std::iterable::{Iterable, Reverse};
#[cfg(not(feature = "std"))]
use ark_std::vec::Vec;

/// Polynomial evaluation, assuming that the
/// coefficients are in little-endian.
#[inline]
fn evaluate_le<F>(polynomial: &[F], x: &F) -> F
where
    F: Field,
{
    evaluate_be(polynomial.iter().rev(), x)
}

/// Polynomial evaluation, assuming that the
/// coeffients are in big-endian.
#[inline]
fn evaluate_be<I, F>(polynomial: I, x: &F) -> F
where
    F: Field,
    I: IntoIterator,
    I::Item: Borrow<F>,
{
    polynomial
        .into_iter()
        .fold(F::zero(), |previous, c| previous * x + c.borrow())
}

#[test]
fn test_commitment_consistency() {
    let rng = &mut ark_std::test_rng();
    let d = 15;
    let polynomial = DensePolynomial::<Fr>::rand(d, rng);
    let polynomial_stream = Reverse(polynomial.coeffs());
    let time_ck = CommitterKey::<Bls12_381>::new(d + 1, 3, rng);
    let space_ck = CommitterKeyStream::from(&time_ck);

    // compute the time commitment
    let time_commitment = time_ck.commit(&polynomial);
    let space_commitment = space_ck.commit(&polynomial_stream);

    assert_eq!(space_commitment, time_commitment);
}

#[test]
fn test_ck_consistency() {
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
    let polynomial_stream = Reverse(polynomials.coeffs());
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
    let mut rng = &mut ark_std::test_rng();
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

#[test]
fn test_vanishing_polynomial() {
    use ark_bls12_381::Fr as F;
    use ark_ff::Zero;

    let points = [F::from(10u64), F::from(5u64), F::from(13u64)];
    let zeros = vanishing_polynomial(&points);
    assert_eq!(evaluate_le(&zeros, &points[0]), F::zero());
    assert_eq!(evaluate_le(&zeros, &points[1]), F::zero());
    assert_eq!(evaluate_le(&zeros, &points[2]), F::zero());
}

#[test]
fn test_srs() {
    use ark_bls12_381::Bls12_381;

    let rng = &mut ark_std::test_rng();
    let ck = CommitterKey::<Bls12_381>::new(10, 3, rng);
    let vk = VerifierKey::from(&ck);
    // Make sure that there are enough elements for the entire array.
    assert_eq!(ck.powers_of_g.len(), 11);
    assert_eq!(ck.powers_of_g2, &vk.powers_of_g2[..]);
}

#[test]
fn test_trivial_commitment() {
    use ark_bls12_381::Bls12_381;
    use ark_bls12_381::Fr;
    use ark_poly::univariate::DensePolynomial;
    use ark_poly::DenseUVPolynomial;
    use ark_std::One;

    let rng = &mut ark_std::test_rng();
    let ck = CommitterKey::<Bls12_381>::new(10, 3, rng);
    let vk = VerifierKey::from(&ck);
    let polynomial = DensePolynomial::from_coefficients_slice(&[Fr::zero(), Fr::one(), Fr::one()]);
    let alpha = Fr::zero();

    let commitment = ck.commit(&polynomial);
    let (evaluation, proof) = ck.open(&polynomial, &alpha);
    assert_eq!(evaluation, Fr::zero());
    assert!(vk.verify(&commitment, &alpha, &evaluation, &proof).is_ok())
}

#[test]
fn test_commitment() {
    use ark_bls12_381::Bls12_381;
    use ark_bls12_381::Fr;
    use ark_poly::univariate::DensePolynomial;
    use ark_poly::DenseUVPolynomial;
    use ark_poly::Polynomial;

    let rng = &mut ark_std::test_rng();
    let ck = CommitterKey::<Bls12_381>::new(100, 3, rng);
    let vk = VerifierKey::from(&ck);
    let polynomial = DensePolynomial::rand(100, rng);
    let alpha = Fr::zero();

    let commitment = ck.commit(&polynomial);
    let (evaluation, proof) = ck.open(&polynomial, &alpha);
    let expected_evaluation = polynomial.evaluate(&alpha);
    assert_eq!(evaluation, expected_evaluation);
    assert!(vk.verify(&commitment, &alpha, &evaluation, &proof).is_ok())
}

#[test]
fn test_open_multi_points() {
    use crate::ark_std::UniformRand;
    use ark_bls12_381::{Bls12_381, Fr};
    use ark_ff::Field;
    use ark_poly::univariate::DensePolynomial;
    use ark_poly::DenseUVPolynomial;
    use ark_std::test_rng;

    let max_msm_buffer = 1 << 20;
    let rng = &mut test_rng();
    // f = 80*x^6 + 80*x^5 + 88*x^4 + 3*x^3 + 73*x^2 + 7*x + 24
    let polynomial = [
        Fr::from(80u64),
        Fr::from(80u64),
        Fr::from(88u64),
        Fr::from(3u64),
        Fr::from(73u64),
        Fr::from(7u64),
        Fr::from(24u64),
    ];
    let polynomial_stream = &polynomial[..];
    let beta = Fr::from(53u64);

    let time_ck = CommitterKey::<Bls12_381>::new(200, 3, rng);
    let space_ck = CommitterKeyStream::from(&time_ck);

    let (remainder, _commitment) = space_ck.open_multi_points(
        &polynomial_stream,
        &[beta.square(), beta, -beta],
        max_msm_buffer,
    );
    let evaluation_remainder = evaluate_be(&remainder, &beta);
    assert_eq!(evaluation_remainder, Fr::from(1807299544171u64));

    let (remainder, _commitment) =
        space_ck.open_multi_points(&polynomial_stream, &[beta], max_msm_buffer);
    assert_eq!(remainder.len(), 1);

    // get a random polynomial with random coefficient,
    let polynomial = DensePolynomial::rand(100, rng).coeffs().to_vec();
    let polynomial_stream = &polynomial[..];
    let beta = Fr::rand(rng);
    let (_, evaluation_proof_batch) =
        space_ck.open_multi_points(&polynomial_stream, &[beta], max_msm_buffer);
    let (_, evaluation_proof_single) = space_ck.open(&polynomial_stream, &beta, max_msm_buffer);
    assert_eq!(evaluation_proof_batch, evaluation_proof_single);

    let (remainder, _evaluation_poof) = space_ck.open_multi_points(
        &polynomial_stream,
        &[beta, -beta, beta.square()],
        max_msm_buffer,
    );
    let expected_evaluation = evaluate_be(&remainder, &beta);
    let obtained_evaluation = evaluate_be(&polynomial, &beta);
    assert_eq!(expected_evaluation, obtained_evaluation);
    let expected_evaluation = evaluate_be(&remainder, &beta.square());
    let obtained_evaluation = evaluate_be(&polynomial, &beta.square());
    assert_eq!(expected_evaluation, obtained_evaluation);
    // let expected_evaluation = evaluate_be(&remainder, &beta.square());
    // let obtained_evaluation = evaluate_be(&polynomial, &beta.square());
    // assert_eq!(expected_evaluation, obtained_evaluation);
    // let expected_evaluation = evaluate_be(&remainder, &beta.square());
    // let obtained_evaluation = evaluate_be(&polynomial, &beta.square());
    // assert_eq!(expected_evaluation, obtained_evaluation);
}
