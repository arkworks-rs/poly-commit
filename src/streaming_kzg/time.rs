//! An impementation of a time-efficient version of Kate et al's polynomial commitment,
//! with optimization from [\[BDFG20\]](https://eprint.iacr.org/2020/081.pdf).
use ark_ec::msm::FixedBase;
use ark_ec::PairingEngine;
use ark_ec::{AffineCurve, ProjectiveCurve};
use ark_ff::{PrimeField, Zero};
use ark_poly::{univariate::DensePolynomial, UVPolynomial};
use ark_std::borrow::Borrow;
use ark_std::ops::Div;
use ark_std::rand::RngCore;
use ark_std::vec::Vec;
use ark_std::UniformRand;

use crate::streaming_kzg::{
    linear_combination, msm, powers, Commitment, EvaluationProof, VerifierKey,
};

use super::vanishing_polynomial;

/// The SRS for the polynomial commitment scheme for a max
///
/// The SRS consists of the `max_degree` powers of \\(\tau\\) in \\(\GG_1\\)
/// plus the `max_eval_degree` powers over \\(\GG_2\\),
/// where `max_degree` is the max polynomial degree to commit to,
/// and `max_eval_degree` is the max number of different points to open simultaneously.
pub struct CommitterKey<E: PairingEngine> {
    pub(crate) powers_of_g: Vec<E::G1Affine>,
    pub(crate) powers_of_g2: Vec<E::G2Affine>,
}

impl<E: PairingEngine> From<&CommitterKey<E>> for VerifierKey<E> {
    fn from(ck: &CommitterKey<E>) -> VerifierKey<E> {
        let max_eval_points = ck.max_eval_points();
        let powers_of_g2 = ck.powers_of_g2[..max_eval_points + 1].to_vec();
        let powers_of_g = ck.powers_of_g[..max_eval_points].to_vec();

        VerifierKey {
            powers_of_g,
            powers_of_g2,
        }
    }
}

impl<E: PairingEngine> CommitterKey<E> {
    /// The setup algorithm for the commitment scheme.
    ///
    /// Given a degree bound `max_degree`,
    /// an evaluation point bound `max_eval_points`,
    /// and a cryptographically-secure random number generator `rng`,
    /// construct the committer key.
    pub fn new(max_degree: usize, max_eval_points: usize, rng: &mut impl RngCore) -> Self {
        // Compute the consecutive powers of an element.
        let tau = E::Fr::rand(rng);
        let powers_of_tau = powers(tau, max_degree + 1);

        let g = E::G1Projective::rand(rng);
        let window_size = FixedBase::get_mul_window_size(max_degree + 1);
        let scalar_bits = E::Fr::MODULUS_BIT_SIZE as usize;
        let g_table = FixedBase::get_window_table(scalar_bits, window_size, g);
        let powers_of_g_proj = FixedBase::msm(scalar_bits, window_size, &g_table, &powers_of_tau);
        let powers_of_g = E::G1Projective::batch_normalization_into_affine(&powers_of_g_proj);

        let g2 = E::G2Projective::rand(rng).into_affine();
        let powers_of_g2 = powers_of_tau
            .iter()
            .take(max_eval_points + 1)
            .map(|t| g2.mul(t.into_bigint()).into_affine())
            .collect::<Vec<_>>();

        CommitterKey {
            powers_of_g,
            powers_of_g2,
        }
    }

    /// Return the bound on evaluation points.
    #[inline]
    pub fn max_eval_points(&self) -> usize {
        self.powers_of_g2.len() - 1
    }

    /// Given a polynomial `polynomial` of degree less than `max_degree`, return a commitment to `polynomial`.
    pub fn commit(&self, polynomial: &[E::Fr]) -> Commitment<E> {
        Commitment(msm::<E>(&self.powers_of_g, polynomial))
    }

    /// Obtain a new preprocessed committer key defined by the indices `indices`.
    pub fn index_by(&self, indices: &[usize]) -> Self {
        let mut indexed_powers_of_g = vec![E::G1Affine::zero(); self.powers_of_g.len()];
        indices
            .iter()
            .zip(self.powers_of_g.iter())
            .for_each(|(&i, &g)| indexed_powers_of_g[i] = indexed_powers_of_g[i] + g);
        Self {
            powers_of_g2: self.powers_of_g2.clone(),
            powers_of_g: indexed_powers_of_g,
        }
    }

    /// Given an iterator over `polynomials`, expressed as vectors of coefficients, return a vector of commitmetns to all of them.
    pub fn batch_commit<J>(&self, polynomials: J) -> Vec<Commitment<E>>
    where
        J: IntoIterator,
        J::Item: Borrow<Vec<E::Fr>>,
    {
        polynomials
            .into_iter()
            .map(|p| self.commit(p.borrow()))
            .collect::<Vec<_>>()
    }

    /// Given a polynomial `polynomial` and an evaluation point `evaluation_point`,
    /// return the evaluation of `polynomial in `evaluation_point`,
    /// together with an evaluation proof.
    pub fn open(
        &self,
        polynomial: &[E::Fr],
        evalualtion_point: &E::Fr,
    ) -> (E::Fr, EvaluationProof<E>) {
        let mut quotient = Vec::new();

        let mut previous = E::Fr::zero();
        for &c in polynomial.iter().rev() {
            let coefficient = c + previous * evalualtion_point;
            quotient.insert(0, coefficient);
            previous = coefficient;
        }

        let (&evaluation, quotient) = quotient.split_first().unwrap_or((&E::Fr::zero(), &[]));
        let evaluation_proof = msm::<E>(&self.powers_of_g, quotient);
        (evaluation, EvaluationProof(evaluation_proof))
    }

    /// Evaluate a single polynomial at a set of points `eval_points`, and provide a single evaluation proof.
    pub fn open_multi_points(
        &self,
        polynomial: &[E::Fr],
        eval_points: &[E::Fr],
    ) -> EvaluationProof<E> {
        // Computing the vanishing polynomial over eval_points
        let z_poly = vanishing_polynomial(eval_points);

        let f_poly = DensePolynomial::from_coefficients_slice(polynomial);
        let q_poly = f_poly.div(&z_poly);
        EvaluationProof(self.commit(&q_poly.coeffs).0)
    }

    /// Evaluate a set of polynomials at a set of points `eval_points`, and provide a single batched evaluation proof.
    /// `eval_chal` is the random challenge for batching evaluation proofs across different polynomials.
    pub fn batch_open_multi_points(
        &self,
        polynomials: &[&Vec<E::Fr>],
        eval_points: &[E::Fr],
        eval_chal: &E::Fr,
    ) -> EvaluationProof<E> {
        assert!(eval_points.len() < self.powers_of_g2.len());
        let etas = powers(*eval_chal, polynomials.len());
        let batched_polynomial =
            linear_combination(polynomials, &etas).unwrap_or_else(|| vec![E::Fr::zero()]);
        self.open_multi_points(&batched_polynomial, eval_points)
    }
}
