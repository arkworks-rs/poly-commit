//! Space-efficient implementation of the polynomial commitment of Kate et al.
use ark_ec::{pairing::Pairing, CurveGroup};
use ark_ff::{PrimeField, Zero};
use ark_poly::Polynomial;
use ark_std::borrow::Borrow;
use ark_std::collections::VecDeque;
use ark_std::vec::Vec;

use crate::streaming_kzg::{ceil_div, vanishing_polynomial, FoldedPolynomialTree};
use ark_ec::scalar_mul::variable_base::{ChunkedPippenger, HashMapPippenger, VariableBaseMSM};
use ark_std::iterable::{Iterable, Reverse};

use super::{time::CommitterKey, VerifierKey};
use super::{Commitment, EvaluationProof};

const LENGTH_MISMATCH_MSG: &str = "Expecting at least one element in the committer key.";

/// The streaming SRS for the polynomial commitment scheme consists of a stream of consecutive powers of g.
/// It also implements functions for `setup`, `commit` and `open`.
#[derive(Clone)]
pub struct CommitterKeyStream<E, SG>
where
    E: Pairing,
    SG: Iterable,
    SG::Item: Borrow<E::G1Affine>,
{
    /// Stream of G1 elements.
    pub powers_of_g: SG,
    /// Two G2 elements needed for the committer.
    pub powers_of_g2: Vec<E::G2Affine>,
}

impl<E, SG> CommitterKeyStream<E, SG>
where
    E: Pairing,
    SG: Iterable,
    SG::Item: Borrow<E::G1Affine>,
{
    /// Turn a streaming SRS into a normal SRS.
    pub fn as_committer_key(&self, max_degree: usize) -> CommitterKey<E> {
        let offset = self.powers_of_g.len() - max_degree;
        let mut powers_of_g = self
            .powers_of_g
            .iter()
            .skip(offset)
            .map(|x| *x.borrow())
            .collect::<Vec<_>>();
        powers_of_g.reverse();
        let powers_of_g2 = self.powers_of_g2.clone().to_vec();
        CommitterKey {
            powers_of_g,
            powers_of_g2,
        }
    }

    /// Evaluate a single polynomial at the point `alpha`, and provide an evaluation proof along with the evaluation.
    pub fn open<SF>(
        &self,
        polynomial: &SF,
        alpha: &E::ScalarField,
        max_msm_buffer: usize,
    ) -> (E::ScalarField, EvaluationProof<E>)
    where
        SF: Iterable,
        SF::Item: Borrow<E::ScalarField>,
    {
        let mut quotient: ChunkedPippenger<E::G1> = ChunkedPippenger::new(max_msm_buffer);

        let bases_init = self.powers_of_g.iter();
        let scalars = polynomial.iter();

        // align the streams and remove one degree
        // TODO: change `skip` to `advance_by` once rust-lang/rust#7774 is fixed.
        // See <https://github.com/rust-lang/rust/issues/77404>
        let bases = bases_init.skip(self.powers_of_g.len() - polynomial.len());

        let mut previous = E::ScalarField::zero();
        for (scalar, base) in scalars.zip(bases) {
            quotient.add(base, previous.into_bigint());
            let coefficient = previous * alpha + scalar.borrow();
            previous = coefficient;
        }

        let evaluation = previous;
        let evaluation_proof = quotient.finalize().into_affine();
        (evaluation, EvaluationProof(evaluation_proof))
    }

    /// Evaluate a single polynomial at a set of points `points`, and provide an evaluation proof along with evaluations.
    pub fn open_multi_points<SF>(
        &self,
        polynomial: &SF,
        points: &[E::ScalarField],
        max_msm_buffer: usize,
    ) -> (Vec<E::ScalarField>, EvaluationProof<E>)
    where
        SF: Iterable,
        SF::Item: Borrow<E::ScalarField>,
    {
        let zeros = vanishing_polynomial(points);
        let mut quotient: ChunkedPippenger<E::G1> = ChunkedPippenger::new(max_msm_buffer);
        let bases_init = self.powers_of_g.iter();
        // TODO: change `skip` to `advance_by` once rust-lang/rust#7774 is fixed.
        // See <https://github.com/rust-lang/rust/issues/77404>
        let mut bases = bases_init.skip(self.powers_of_g.len() - polynomial.len() + zeros.degree());

        let mut state = VecDeque::<E::ScalarField>::with_capacity(points.len());

        let mut polynomial_iterator = polynomial.iter();

        (0..points.len()).for_each(|_| {
            state.push_back(*polynomial_iterator.next().unwrap().borrow());
        });

        for coefficient in polynomial_iterator {
            let coefficient = coefficient.borrow();
            let quotient_coefficient = state.pop_front().unwrap();
            state.push_back(*coefficient);
            (0..points.len()).for_each(|i| {
                state[i] -= zeros.coeffs[zeros.degree() - i - 1] * quotient_coefficient;
            });
            let base = bases.next().unwrap();
            quotient.add(base, quotient_coefficient.into_bigint());
        }
        let remainder = state.make_contiguous().to_vec();
        let commitment = EvaluationProof(quotient.finalize().into_affine());
        (remainder, commitment)
    }

    /// The commitment procedures, that takes as input a committer key and the streaming coefficients of polynomial, and produces the desired commitment.
    pub fn commit<SF: ?Sized>(&self, polynomial: &SF) -> Commitment<E>
    where
        SF: Iterable,
        SF::Item: Borrow<E::ScalarField>,
    {
        assert!(self.powers_of_g.len() >= polynomial.len());

        Commitment(
            <E::G1 as VariableBaseMSM>::msm_chunks(&self.powers_of_g, polynomial).into_affine(),
        )
    }

    /// The batch commitment procedures, that takes as input a committer key and the streaming coefficients of a list of polynomials, and produces the desired commitments.
    pub fn batch_commit<'a, F>(
        &self,
        polynomials: &[&'a dyn Iterable<Item = F, Iter = &mut dyn Iterator<Item = F>>],
    ) -> Vec<Commitment<E>>
    where
        F: Borrow<E::ScalarField>,
    {
        polynomials.iter().map(|&p| self.commit(p)).collect()
    }

    /// The commitment procedures for our tensor check protocol.
    /// The algorithm takes advantage of the tree structure of folding polynomials in our protocol. Please refer to our paper for more details.
    /// The function takes as input a committer key and the tree structure of all the folding polynomials, and produces the desired commitment for each polynomial.
    pub fn commit_folding<SF>(
        &self,
        polynomials: &FoldedPolynomialTree<E::ScalarField, SF>,
        max_msm_buffer: usize,
    ) -> Vec<Commitment<E>>
    where
        SF: Iterable,
        SF::Item: Borrow<E::ScalarField>,
    {
        let n = polynomials.depth();
        let mut pippengers: Vec<ChunkedPippenger<E::G1>> = Vec::new();
        let mut folded_bases = Vec::new();
        for i in 1..n + 1 {
            let pippenger: ChunkedPippenger<<E as Pairing>::G1> =
                ChunkedPippenger::with_size(max_msm_buffer / n);
            let bases_init = self.powers_of_g.iter();

            let delta = self.powers_of_g.len() - ceil_div(polynomials.len(), 1 << i);
            // TODO: change `skip` to `advance_by` once rust-lang/rust#7774 is fixed.
            // See <https://github.com/rust-lang/rust/issues/77404>
            let bases = bases_init.skip(delta);
            folded_bases.push(bases);
            pippengers.push(pippenger);
        }

        for (i, coefficient) in polynomials.iter() {
            let base = folded_bases[i - 1].next().unwrap();
            pippengers[i - 1].add(base.borrow(), coefficient.into_bigint());
        }

        pippengers
            .into_iter()
            .map(|p| Commitment(p.finalize().into_affine()))
            .collect::<Vec<_>>()
    }

    /// The commitment procedures for our tensor check protocol.
    /// The algorithm takes advantage of the tree structure of folding polynomials in our protocol. Please refer to our paper for more details.
    /// The function evaluates all the folding polynomials at a set of evaluation points `points` and produces a single batched evaluation proof.
    /// `eta` is the random challenge for batching folding polynomials.
    pub fn open_folding<'a, SF>(
        &self,
        polynomials: FoldedPolynomialTree<'a, E::ScalarField, SF>,
        points: &[E::ScalarField],
        etas: &[E::ScalarField],
        max_msm_buffer: usize,
    ) -> (Vec<Vec<E::ScalarField>>, EvaluationProof<E>)
    where
        SG: Iterable,
        SF: Iterable,
        E: Pairing,
        SG::Item: Borrow<E::G1Affine>,
        SF::Item: Borrow<E::ScalarField> + Copy,
    {
        let n = polynomials.depth();
        let mut pippenger = HashMapPippenger::<E::G1>::new(max_msm_buffer);
        let mut folded_bases = Vec::new();
        let zeros = vanishing_polynomial(points);
        let mut remainders = vec![VecDeque::new(); n];

        for i in 1..n + 1 {
            let bases_init = self.powers_of_g.iter();
            let delta = self.powers_of_g.len() - ceil_div(polynomials.len(), 1 << i);
            // TODO: change `skip` to `advance_by` once rust-lang/rust#7774 is fixed.
            // See <https://github.com/rust-lang/rust/issues/77404>
            let bases = bases_init.skip(delta);

            (0..points.len()).for_each(|_| {
                remainders[i - 1].push_back(E::ScalarField::zero());
            });

            folded_bases.push(bases);
        }

        for (i, coefficient) in polynomials.iter() {
            if i == 0 {
                continue;
            } // XXX. skip the 0th elements automatically

            let base = folded_bases[i - 1].next().unwrap();
            let quotient_coefficient = remainders[i - 1].pop_front().unwrap();
            remainders[i - 1].push_back(coefficient);
            (0..points.len()).for_each(|j| {
                remainders[i - 1][j] -= zeros.coeffs[zeros.degree() - j - 1] * quotient_coefficient;
            });

            let scalar = etas[i - 1] * quotient_coefficient;
            pippenger.add(base, scalar);
        }

        let evaluation_proof = pippenger.finalize().into_affine();
        let remainders = remainders
            .iter_mut()
            .map(|x| x.make_contiguous().to_vec())
            .collect::<Vec<_>>();

        (remainders, EvaluationProof(evaluation_proof))
    }
}

impl<'a, E: Pairing> From<&'a CommitterKey<E>>
    for CommitterKeyStream<E, Reverse<&'a [E::G1Affine]>>
{
    fn from(ck: &'a CommitterKey<E>) -> Self {
        CommitterKeyStream {
            powers_of_g: Reverse(ck.powers_of_g.as_slice()),
            powers_of_g2: ck.powers_of_g2.clone(),
        }
    }
}

impl<E, SG> From<&CommitterKeyStream<E, SG>> for VerifierKey<E>
where
    E: Pairing,
    SG: Iterable,
    SG::Item: Borrow<E::G1Affine>,
{
    fn from(ck: &CommitterKeyStream<E, SG>) -> Self {
        let powers_of_g2 = ck.powers_of_g2.to_vec();
        // take the first element from the stream
        let g = *ck
            .powers_of_g
            .iter()
            .last()
            .expect(LENGTH_MISMATCH_MSG)
            .borrow();
        Self {
            powers_of_g2,
            powers_of_g: vec![g],
        }
    }
}
