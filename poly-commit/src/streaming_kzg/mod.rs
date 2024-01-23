//! The KZG (or Kate) polynomial commitment, space- and time-efficient.
//!
//! # Background
//! [[KZG](https://www.iacr.org/archive/asiacrypt2010/6477178/6477178.pdf)]
//! commitments are pretty simple:
//! - A [`CommitterKey`](self::CommitterKey) is consists of a sequence \\(\vec G \defeq (G, \tau G, \dots, \tau^DG)\\).
//! - A [`Commitment`](self::EvaluationProof) is a polynomial \\(f(x)\\) is \\(C \defeq \langle \vec f, \vec G \rangle \\).
//! - An [`EvaluationProof`](self::EvaluationProof)
//! for the polynomial \\(f\\)
//! in the evaluation point \\(\alpha\\)
//! is a commitment to the quotient of \\(f(x)\\) by \\((\tau - \alpha)\\).
//! The remainder is the evaluation \\(f(\alpha)\\).
//! When evaluation over points \\((\alpha_0, \dots, \alpha_m)\\),
//! we can consider at once the quotient of \\(f(x)\\) by \\(Z\\) (the polynomial whose roots are \\(\alpha_i\\)).
//! The remainder is a polynomial \\(r\\) such that \\(r(\alpha_i) = f(\alpha_i)\\).
//! We refer to the proof as \\(\pi\\).
//!
//! To verify a proof \\(\pi\\) proving that \\(f(\alpha) = \mu\\), one considers the pairing equation:
//! \\[
//! e(C, \tau H - \mu H) = e(f - \mu G, H)
//! \\]
//! To verify a proof \\(\pi\\) over a set of points \\(f(\alpha_i) = \mu_i\\),
//! consider the polynomial \\(\nu\\) such that \\(\nu(\alpha_i) = \mu_i \\), and check:
//! \\[
//! e(C, Z) = e(f - \nu, H).
//! \\]
//!
//! It is also possible to open multiple polynomials \\(f_0, \dots, f_n\\)
//!  _on the same set of evaluation points_
//! by asking the verifier a random challenge \\(\eta\\), and opening instead
//! \\(\sum_i \eta^i f_i \\).
//!
//! _Nota bene:_ despite it is also possible to open multiple polynomials
//! over different points [[BDFG20](https://eprint.iacr.org/2020/081.pdf)],
//! however this is not currently supported by our implementation.
//!
//!
//! # Examples
//!
//! When creating a new SRS, one must specify a degree bound `max_degree`
//! for the commitment polynomials, and a degree bound `max_evals` for
//! the maximum number of opening points.
//! From the SRS, it is possible to derive the verification key
//! [`VerifierKey`](self::VerifierKey).
//!
//! ```
//! use ark_poly_commit::streaming_kzg::CommitterKey;
//! use ark_bls12_381::{Fr, Bls12_381};
//!
//! let rng = &mut ark_std::test_rng();
//! let max_degree = 100;
//! let max_evals = 10;
//!
//! let ck = CommitterKey::<Bls12_381>::new(max_degree, max_evals, rng);
//! # // XXX. if you change the following lines,
//! # // please note that documentation below might break.
//! # let f = vec![Fr::from(1u64), Fr::from(2u64), Fr::from(4u64), Fr::from(8u64)];
//! # let commitment  = ck.commit(&f);
//! # let alpha = Fr::from(42u64);
//! # let (evaluation, proof) = ck.open(&f, &alpha);
//! # use ark_poly_commit::streaming_kzg::VerifierKey;
//! # let vk = VerifierKey::from(&ck);
//! # assert!(vk.verify(&commitment, &alpha, &evaluation, &proof).is_ok())
//! ```
//!
//! Then to commit to a polynomial `f`:
//! ```ignore
//! let f = vec![Fr::from(1u64), Fr::from(2u64), Fr::from(4u64), Fr::from(8u64)];
//! let commitment  = ck.commit(&f);
//! ```
//! To prove the evaluation of `f` in a point `alpha`:
//!
//! ```ignore
//! let alpha = Fr::from(42u64);
//! let (evaluation, proof) = ck.open(&f, &alpha);
//! ```
//! To veify that an opening is correct:
//! ```ignore
//! use gemini::kzg::VerifierKey;
//!
//! let vk = VerifierKey::from(&ck);
//! assert!(vk.verify(&commitment, &alpha, &evaluation, &proof).is_ok())
//! ```

mod data_structures;
mod space;
mod time;

use ark_ec::CurveGroup;
use ark_serialize::{CanonicalSerialize, Compress};
use ark_std::vec::Vec;
pub use data_structures::*;
pub use space::CommitterKeyStream;
pub use time::CommitterKey;

#[cfg(test)]
pub mod tests;

use ark_ff::{Field, One, PrimeField, Zero};
use ark_poly::{univariate::DensePolynomial, DenseUVPolynomial};
use ark_std::ops::{Add, Mul};

use ark_std::borrow::Borrow;
use ark_std::fmt;

use ark_ec::{pairing::Pairing, AffineRepr, VariableBaseMSM};

/// A Kate polynomial commitment over a bilinear group, represented as a single \\(\GG_1\\) element.
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub struct Commitment<E: Pairing>(pub(crate) E::G1Affine);

impl<E: Pairing> Commitment<E> {
    /// Return the size of Commitment in bytes.
    pub fn size_in_bytes(&self) -> usize {
        E::G1Affine::zero().serialized_size(Compress::Yes)
    }
}

#[inline]
fn msm<E: Pairing>(bases: &[E::G1Affine], scalars: &[E::ScalarField]) -> E::G1Affine {
    let scalars = scalars.iter().map(|x| x.into_bigint()).collect::<Vec<_>>();
    let sp = <E::G1 as VariableBaseMSM>::msm_bigint(bases, &scalars);
    sp.into_affine()
}

/// Polynomial evaluation proof, represented as a single \\(\GG_1\\) element.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct EvaluationProof<E: Pairing>(pub E::G1Affine);

impl<E: Pairing> Add for EvaluationProof<E> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        EvaluationProof((self.0 + rhs.0).into_affine())
    }
}

impl<E: Pairing> core::iter::Sum for EvaluationProof<E> {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        let zero = EvaluationProof(E::G1Affine::zero());
        iter.fold(zero, |x, y| x + y)
    }
}

/// Error type denoting an incorrect evaluation proof.
#[derive(Debug, Clone)]
pub struct VerificationError;

impl fmt::Display for VerificationError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Error in stream.")
    }
}

pub(crate) type VerificationResult = Result<(), VerificationError>;

// XXX.  add const generic argument for the size.
/// The verification key for the polynomial commitment scheme.
/// It also implements verification functions for the evaluation proof.
#[derive(Debug, PartialEq, Eq)]
pub struct VerifierKey<E: Pairing> {
    /// The generator of  \\(\GG_1\\)
    powers_of_g: Vec<E::G1Affine>,
    /// The generator og \\(\GG_2\\), together with its multiplication by the trapdoor.
    powers_of_g2: Vec<E::G2Affine>,
}

impl<E: Pairing> VerifierKey<E> {
    /// The verification procedure for the EvaluationProof with a single polynomial evaluated at a single evaluation point.
    /// The polynomial are evaluated at the point ``alpha`` and is committed as ``commitment``.
    /// The evaluation proof can be obtained either in a space-efficient or a time-efficient flavour.
    pub fn verify(
        &self,
        commitment: &Commitment<E>,
        &alpha: &E::ScalarField,
        evaluation: &E::ScalarField,
        proof: &EvaluationProof<E>,
    ) -> VerificationResult {
        let scalars = [(-alpha).into_bigint(), E::ScalarField::one().into_bigint()];
        let ep = <E::G2 as VariableBaseMSM>::msm_bigint(&self.powers_of_g2, &scalars);
        let lhs = commitment.0.into_group() - self.powers_of_g[0].mul(evaluation);
        let g2 = self.powers_of_g2[0];

        if E::pairing(lhs, g2) == E::pairing(proof.0, ep) {
            Ok(())
        } else {
            Err(VerificationError)
        }
    }

    /// The verification procedure for the EvaluationProof with a set of polynomials evaluated at a set of evaluation points.
    /// All the polynomials are evaluated at the set of points ``eval_points`` and are committed as ``commitments``.
    /// ``evaluations`` contains evaluations of each polynomial at each point in ``eval_points``.
    /// ``evaluations`` follows the same polynomial order as ``commitments`` and the same evaluation point order as ``eval_points``.
    /// The evaluation proof can be obtained either in a space-efficient or a time-efficient flavour.
    /// ``open_chal`` is a random challenge for batching evaluation proofs across different polynomials.
    pub fn verify_multi_points(
        &self,
        commitments: &[Commitment<E>],
        eval_points: &[E::ScalarField],
        evaluations: &[Vec<E::ScalarField>],
        proof: &EvaluationProof<E>,
        open_chal: &E::ScalarField,
    ) -> VerificationResult {
        // Computing the vanishing polynomial over eval_points
        let zeros = vanishing_polynomial(eval_points);
        let zeros_repr = zeros.iter().map(|x| x.into_bigint()).collect::<Vec<_>>();
        let zeros = <E::G2 as VariableBaseMSM>::msm_bigint(&self.powers_of_g2, &zeros_repr);

        // Computing the inverse for the interpolation
        let mut sca_inverse = Vec::new();
        for (j, x_j) in eval_points.iter().enumerate() {
            let mut sca = E::ScalarField::one();
            for (k, x_k) in eval_points.iter().enumerate() {
                if j == k {
                    continue;
                }
                sca *= *x_j - x_k;
            }
            sca = sca.inverse().unwrap();
            sca_inverse.push(sca);
        }

        // Computing the lagrange polynomial for the interpolation
        let mut lang = Vec::new();
        for (j, _x_j) in eval_points.iter().enumerate() {
            let mut l_poly = DensePolynomial::from_coefficients_vec(vec![E::ScalarField::one()]);
            for (k, x_k) in eval_points.iter().enumerate() {
                if j == k {
                    continue;
                }
                let tmp_poly =
                    DensePolynomial::from_coefficients_vec(vec![-(*x_k), E::ScalarField::one()]);
                l_poly = l_poly.mul(&tmp_poly);
            }
            lang.push(l_poly);
        }

        // Computing the commitment for the interpolated polynomials
        let etas = powers(*open_chal, evaluations.len());
        let interpolated_polynomials = evaluations
            .iter()
            .map(|e| interpolate_poly::<E>(eval_points, e, &sca_inverse, &lang).coeffs)
            .collect::<Vec<_>>();
        let i_poly = linear_combination(&interpolated_polynomials[..], &etas).unwrap();

        let i_comm = msm::<E>(&self.powers_of_g, &i_poly);

        // Gathering commitments
        let comm_vec = commitments.iter().map(|x| x.0).collect::<Vec<_>>();
        let etas_repr = etas.iter().map(|e| e.into_bigint()).collect::<Vec<_>>();
        let f_comm = <E::G1 as VariableBaseMSM>::msm_bigint(&comm_vec, &etas_repr);

        let g2 = self.powers_of_g2[0];

        if E::pairing(f_comm - i_comm.into_group(), g2) == E::pairing(proof.0, zeros) {
            Ok(())
        } else {
            Err(VerificationError)
        }
    }
}

fn interpolate_poly<E: Pairing>(
    eval_points: &[E::ScalarField],
    evals: &[E::ScalarField],
    sca_inverse: &[E::ScalarField],
    lang: &[DensePolynomial<E::ScalarField>],
) -> DensePolynomial<E::ScalarField> {
    let mut res = DensePolynomial::from_coefficients_vec(vec![E::ScalarField::zero()]);
    for (j, (_x_j, y_j)) in eval_points.iter().zip(evals.iter()).enumerate() {
        let l_poly = lang[j].mul(sca_inverse[j] * y_j);
        res = (&res).add(&l_poly);
    }
    res
}

/// The polynomial in \\(\FF\\) that vanishes in all the points `points`.
pub(crate) fn vanishing_polynomial<F: Field>(points: &[F]) -> DensePolynomial<F> {
    let one = DensePolynomial::from_coefficients_vec(vec![F::one()]);
    points
        .iter()
        .map(|&point| DensePolynomial::from_coefficients_vec(vec![-point, F::one()]))
        .fold(one, |x, y| x.naive_mul(&y))
}

/// Compute a linear combination of the polynomials `polynomials` with the given challenges.
pub(crate) fn linear_combination<F: Field, PP>(
    polynomials: &[PP],
    challenges: &[F],
) -> Option<Vec<F>>
where
    PP: Borrow<Vec<F>>,
{
    polynomials
        .iter()
        .zip(challenges.iter())
        .map(|(p, &c)| &DensePolynomial::from_coefficients_vec(p.borrow().to_vec()) * c)
        .reduce(|x, y| x + y)?
        .coeffs
        .into()
}

/// Return a vector of length `len` containing the consecutive powers of element.
pub(crate) fn powers<F: Field>(element: F, len: usize) -> Vec<F> {
    let mut powers = vec![F::one(); len];
    for i in 1..len {
        powers[i] = element * powers[i - 1];
    }
    powers
}
