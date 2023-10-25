use ark_ec::AffineRepr;
use ark_ff::PrimeField;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::{rand::RngCore, vec::Vec};

use crate::{
    PCCommitment, PCCommitterKey, PCPreparedCommitment, PCPreparedVerifierKey, PCRandomness,
    PCUniversalParams, PCVerifierKey,
};

/// `UniversalParams` amounts to a Pederson commitment key of sufficient length
#[derive(Derivative, CanonicalSerialize, CanonicalDeserialize)]
#[derivative(Default(bound = ""), Clone(bound = ""), Debug(bound = ""))]
pub struct HyraxUniversalParams<G: AffineRepr> {
    /// A list of generators of the group.
    pub com_key: Vec<G>,
    /// A generator of the group.
    pub h: G,
}

impl<G: AffineRepr> PCUniversalParams for HyraxUniversalParams<G> {
    fn max_degree(&self) -> usize {
        // Only MLEs are supported
        1
    }
}

/// The committer key, which coincides with the universal parameters
pub type HyraxCommitterKey<G> = HyraxUniversalParams<G>;

/// The verifier key, which coincides with the committer key
pub type HyraxVerifierKey<G> = HyraxCommitterKey<G>;

impl<G: AffineRepr> PCCommitterKey for HyraxCommitterKey<G> {
    fn max_degree(&self) -> usize {
        // Only MLEs are supported
        1
    }
    fn supported_degree(&self) -> usize {
        // Only MLEs are supported
        1
    }
}

impl<G: AffineRepr> PCVerifierKey for HyraxVerifierKey<G> {
    // Only MLEs are supported
    fn max_degree(&self) -> usize {
        1
    }
    // Only MLEs are supported
    fn supported_degree(&self) -> usize {
        1
    }
}

/// Nothing to do to prepare this prover-verifier key.
pub type HyraxPreparedVerifierKey<G> = HyraxVerifierKey<G>;

impl<G: AffineRepr> PCPreparedVerifierKey<HyraxVerifierKey<G>> for HyraxPreparedVerifierKey<G> {
    /// Simply clone the prover-verifier key
    fn prepare(vk: &HyraxVerifierKey<G>) -> Self {
        vk.clone()
    }
}

/// Hyrax commitment to a polynomial consisting of one multi-commit per row of
/// the coefficient matrix
#[derive(Derivative, CanonicalSerialize, CanonicalDeserialize)]
#[derivative(Default(bound = ""), Clone(bound = ""), Debug(bound = ""))]
pub struct HyraxCommitment<G: AffineRepr> {
    /// A list of multi-commits to each row of the matrix representing the
    /// polynomial.
    pub row_coms: Vec<G>,
}

impl<G: AffineRepr> PCCommitment for HyraxCommitment<G> {
    #[inline]
    fn empty() -> Self {
        HyraxCommitment {
            row_coms: Vec::new(),
        }
    }

    // The degree bound is always 1, since only multilinear polynomials are
    // supported
    fn has_degree_bound(&self) -> bool {
        true
    }
}

/// No preparation is needed for Hyrax commitments
pub type HyraxPreparedCommitment<E> = HyraxCommitment<E>;

impl<G: AffineRepr> PCPreparedCommitment<HyraxCommitment<G>> for HyraxPreparedCommitment<G> {
    /// Simply clone the prover-verifier key
    fn prepare(vk: &HyraxCommitment<G>) -> Self {
        vk.clone()
    }
}

pub(crate) type HyraxRandomness<F> = Vec<F>;

/// A vector of scalars, each of which multiplies the distinguished group
/// element in the Pederson commitment key for a different commitment
impl<F: PrimeField> PCRandomness for HyraxRandomness<F> {
    fn empty() -> Self {
        unimplemented!()
    }

    fn rand<R: RngCore>(
        num_queries: usize,
        _has_degree_bound: bool,
        _num_vars: Option<usize>,
        rng: &mut R,
    ) -> Self {
        (0..num_queries).map(|_| F::rand(rng)).collect()
    }
}

/// Proof of a Hyrax opening, containing various commitments
/// and auxiliary values generated randomly during the opening
#[derive(Derivative, CanonicalSerialize, CanonicalDeserialize)]
#[derivative(Default(bound = ""), Clone(bound = ""), Debug(bound = ""))]
pub struct HyraxProof<G: AffineRepr> {
    /// Commitment to the evaluation of the polynomial at the requested point
    pub com_eval: G,
    /// Commitment to auxiliary random vector `d`
    pub com_d: G,
    /// Commitment to auxiliary random scalar `b`
    pub com_b: G,
    /// Auxiliary random vector
    pub z: Vec<G::ScalarField>,
    /// Auxiliary random scalar
    pub z_d: G::ScalarField,
    /// Auxiliary random scalar
    pub z_b: G::ScalarField,
    /// The hiding scalar r_eval is not part of a Hyrax PCS proof as described
    /// in the reference article. Cf. the "Modification note" at the beginning
    /// of `mod.rs`
    pub r_eval: G::ScalarField,
}
