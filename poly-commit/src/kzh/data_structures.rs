use crate::{
    PCCommitment, PCCommitmentState, PCCommitterKey, PCPreparedCommitment, PCPreparedVerifierKey,
    PCUniversalParams, PCVerifierKey,
};
use ark_ec::{pairing::Pairing, AffineRepr};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::{rand::RngCore, vec::Vec};

/// Universal parameters for the KZH-`K` commitment scheme.
///
/// Tensor blocks and `h` layers follow arkworks' little-endian variable order,
/// from the lowest variables (`x_0` first) to the highest. The `j`-th `h` layer
/// contains the flattened KZH tensor for the suffix of blocks beginning at
/// block `j`, with that suffix's first block varying fastest. The `v_tau`
/// vectors cover only the first `K - 1` blocks because the highest/final block
/// is checked from the disclosed final evaluation vector and never uses a
/// pairing layer.
#[derive(Derivative, CanonicalSerialize, CanonicalDeserialize)]
#[derivative(
    Clone(bound = ""),
    Debug(bound = ""),
    PartialEq(bound = ""),
    Eq(bound = "")
)]
pub struct UniversalParams<E: Pairing, const K: usize> {
    /// Runtime copy of the family parameter, included in the serialized form.
    pub(crate) k: usize,
    /// Number of variables supported by these parameters.
    pub(crate) num_vars: usize,
    /// Number of variables in each of the `K` balanced tensor blocks, ordered
    /// low-to-high in arkworks' variable order.
    pub(crate) num_vars_per_block: Vec<usize>,
    /// Flattened suffix commitment tensors, one for each tensor block.
    pub(crate) h: Vec<Vec<E::G1Affine>>,
    /// The verifier's base generator in `G2`.
    pub(crate) v: E::G2Affine,
    /// Trapdoor-scaled `G2` bases for every tensor block except the last.
    pub(crate) v_tau: Vec<Vec<E::G2Affine>>,
}

impl<E: Pairing, const K: usize> UniversalParams<E, K> {
    /// Returns the exact number of variables supported by these parameters.
    #[inline]
    pub const fn num_vars(&self) -> usize {
        self.num_vars
    }

    /// Returns the low-to-high variable counts of the balanced tensor blocks.
    #[inline]
    pub fn num_vars_per_block(&self) -> &[usize] {
        &self.num_vars_per_block
    }
}

impl<E: Pairing, const K: usize> PCUniversalParams for UniversalParams<E, K> {
    fn max_degree(&self) -> usize {
        // KZH commits to multilinear polynomials.
        1
    }
}

/// Committer key for the KZH-`K` commitment scheme.
#[derive(Derivative, CanonicalSerialize, CanonicalDeserialize)]
#[derivative(
    Clone(bound = ""),
    Debug(bound = ""),
    PartialEq(bound = ""),
    Eq(bound = "")
)]
pub struct CommitterKey<E: Pairing, const K: usize> {
    /// Runtime copy of the family parameter, included in the serialized form.
    pub(crate) k: usize,
    /// Number of variables supported by this key.
    pub(crate) num_vars: usize,
    /// Number of variables in each of the `K` balanced tensor blocks, ordered
    /// low-to-high in arkworks' variable order.
    pub(crate) num_vars_per_block: Vec<usize>,
    /// Flattened suffix commitment tensors, one for each tensor block.
    pub(crate) h: Vec<Vec<E::G1Affine>>,
}

impl<E: Pairing, const K: usize> CommitterKey<E, K> {
    /// Returns the exact number of variables supported by this key.
    #[inline]
    pub const fn num_vars(&self) -> usize {
        self.num_vars
    }

    /// Returns the low-to-high variable counts of the balanced tensor blocks.
    #[inline]
    pub fn num_vars_per_block(&self) -> &[usize] {
        &self.num_vars_per_block
    }
}

impl<E: Pairing, const K: usize> PCCommitterKey for CommitterKey<E, K> {
    fn max_degree(&self) -> usize {
        // KZH commits to multilinear polynomials.
        1
    }

    fn supported_degree(&self) -> usize {
        // KZH commits to multilinear polynomials.
        1
    }
}

/// Verifier key for the KZH-`K` commitment scheme.
#[derive(Derivative, CanonicalSerialize, CanonicalDeserialize)]
#[derivative(
    Clone(bound = ""),
    Debug(bound = ""),
    PartialEq(bound = ""),
    Eq(bound = "")
)]
pub struct VerifierKey<E: Pairing, const K: usize> {
    /// Runtime copy of the family parameter, included in the serialized form.
    pub(crate) k: usize,
    /// Number of variables supported by this key.
    pub(crate) num_vars: usize,
    /// Number of variables in each of the `K` balanced tensor blocks, ordered
    /// low-to-high in arkworks' variable order.
    pub(crate) num_vars_per_block: Vec<usize>,
    /// The commitment layer for the highest/final tensor block.
    pub(crate) h_last: Vec<E::G1Affine>,
    /// The verifier's base generator in `G2`.
    pub(crate) v: E::G2Affine,
    /// Trapdoor-scaled `G2` bases for every tensor block except the last.
    pub(crate) v_tau: Vec<Vec<E::G2Affine>>,
}

impl<E: Pairing, const K: usize> VerifierKey<E, K> {
    /// Returns the exact number of variables supported by this key.
    #[inline]
    pub const fn num_vars(&self) -> usize {
        self.num_vars
    }

    /// Returns the low-to-high variable counts of the balanced tensor blocks.
    #[inline]
    pub fn num_vars_per_block(&self) -> &[usize] {
        &self.num_vars_per_block
    }
}

impl<E: Pairing, const K: usize> PCVerifierKey for VerifierKey<E, K> {
    fn max_degree(&self) -> usize {
        // KZH commits to multilinear polynomials.
        1
    }

    fn supported_degree(&self) -> usize {
        // KZH commits to multilinear polynomials.
        1
    }
}

/// Prepared verifier key for repeated KZH verification.
///
/// This is an ephemeral runtime cache: every `G2` input used by the pairing
/// equations is converted to Arkworks' prepared representation once and then
/// reused by [`KZH::check_prepared`](super::KZH::check_prepared). Prepared line
/// coefficients are deliberately not canonically serializable; deserialize
/// and validate an ordinary [`VerifierKey`] and prepare it locally instead.
#[derive(Derivative)]
#[derivative(Clone(bound = ""), Debug(bound = ""))]
pub struct PreparedVerifierKey<E: Pairing, const K: usize> {
    verifier_key: VerifierKey<E, K>,
    prepared_v: E::G2Prepared,
    prepared_v_tau: Vec<Vec<E::G2Prepared>>,
}

impl<E: Pairing, const K: usize> PreparedVerifierKey<E, K> {
    /// Prepares all `G2` pairing inputs in `vk` for repeated verification.
    pub fn prepare(vk: &VerifierKey<E, K>) -> Self {
        Self::from(vk)
    }

    /// Returns the ordinary verifier key from which this cache was derived.
    #[inline]
    pub fn verifier_key(&self) -> &VerifierKey<E, K> {
        &self.verifier_key
    }

    /// Returns the exact number of variables supported by this key.
    #[inline]
    pub const fn num_vars(&self) -> usize {
        self.verifier_key.num_vars
    }

    /// Returns the low-to-high variable counts of the balanced tensor blocks.
    #[inline]
    pub fn num_vars_per_block(&self) -> &[usize] {
        &self.verifier_key.num_vars_per_block
    }

    #[inline]
    pub(crate) fn prepared_v(&self) -> &E::G2Prepared {
        &self.prepared_v
    }

    #[inline]
    pub(crate) fn prepared_v_tau(&self) -> &[Vec<E::G2Prepared>] {
        &self.prepared_v_tau
    }
}

impl<E: Pairing, const K: usize> From<&VerifierKey<E, K>> for PreparedVerifierKey<E, K> {
    fn from(vk: &VerifierKey<E, K>) -> Self {
        let prepared_v = E::G2Prepared::from(&vk.v);
        let prepared_v_tau = vk
            .v_tau
            .iter()
            .map(|layer| layer.iter().map(E::G2Prepared::from).collect())
            .collect();

        Self {
            verifier_key: vk.clone(),
            prepared_v,
            prepared_v_tau,
        }
    }
}

impl<E: Pairing, const K: usize> PCPreparedVerifierKey<VerifierKey<E, K>>
    for PreparedVerifierKey<E, K>
{
    fn prepare(vk: &VerifierKey<E, K>) -> Self {
        Self::from(vk)
    }
}

/// A constant-size KZH commitment.
#[derive(Derivative, CanonicalSerialize, CanonicalDeserialize)]
#[derivative(
    Clone(bound = ""),
    Copy(bound = ""),
    Debug(bound = ""),
    PartialEq(bound = ""),
    Eq(bound = "")
)]
pub struct Commitment<E: Pairing, const K: usize> {
    /// Runtime copy of the family parameter, included in the serialized form.
    pub(crate) k: usize,
    /// Number of variables in the committed polynomial.
    pub(crate) num_vars: usize,
    /// The commitment group element.
    pub(crate) comm: E::G1Affine,
}

impl<E: Pairing, const K: usize> Commitment<E, K> {
    /// Returns the number of variables in the committed polynomial.
    #[inline]
    pub const fn num_vars(&self) -> usize {
        self.num_vars
    }

    /// Returns the underlying commitment group element.
    #[inline]
    pub fn comm(&self) -> &E::G1Affine {
        &self.comm
    }
}

impl<E: Pairing, const K: usize> Default for Commitment<E, K> {
    fn default() -> Self {
        Self {
            k: K,
            num_vars: 0,
            comm: E::G1Affine::zero(),
        }
    }
}

impl<E: Pairing, const K: usize> PCCommitment for Commitment<E, K> {
    #[inline]
    fn empty() -> Self {
        Self::default()
    }

    fn has_degree_bound(&self) -> bool {
        // KZH enforces multilinearity through the polynomial type, but does
        // not authenticate strict degree-bound metadata on commitments.
        false
    }
}

/// Prepared KZH commitment.
///
/// KZH currently performs no additional commitment preparation.
pub type PreparedCommitment<E, const K: usize> = Commitment<E, K>;

impl<E: Pairing, const K: usize> PCPreparedCommitment<Commitment<E, K>>
    for PreparedCommitment<E, K>
{
    fn prepare(commitment: &Commitment<E, K>) -> Self {
        *commitment
    }
}

/// Private state cached while committing to a polynomial.
///
/// `auxiliary_tables[j]` contains higher-variable suffix commitments for all
/// assignments to low blocks `0..=j`, stored current-block-index first and
/// prior-prefix-index second. During a generic opening, early proof layers are
/// obtained by contracting these tables with the already fixed low point
/// blocks. Later layers are committed directly from the partially evaluated
/// polynomial; only the beneficial prefix of tables is stored. Consequently,
/// the number of tables is generally smaller than the number of proof
/// transition vectors. Batched openings may consume an even shorter prefix
/// when fusing many states would cost more than the direct route.
#[derive(Derivative, CanonicalSerialize, CanonicalDeserialize)]
#[derivative(
    Clone(bound = ""),
    Debug(bound = ""),
    PartialEq(bound = ""),
    Eq(bound = "")
)]
pub struct CommitmentState<E: Pairing, const K: usize> {
    /// Runtime copy of the family parameter, included in the serialized form.
    pub(crate) k: usize,
    /// Number of variables in the committed polynomial.
    pub(crate) num_vars: usize,
    /// Suffix-commitment tables for the beneficial low-block prefix.
    pub(crate) auxiliary_tables: Vec<Vec<E::G1Affine>>,
}

impl<E: Pairing, const K: usize> CommitmentState<E, K> {
    /// Returns the number of variables supported by this commitment state.
    #[inline]
    pub const fn num_vars(&self) -> usize {
        self.num_vars
    }
}

impl<E: Pairing, const K: usize> Default for CommitmentState<E, K> {
    fn default() -> Self {
        Self {
            k: K,
            num_vars: 0,
            auxiliary_tables: Vec::new(),
        }
    }
}

impl<E: Pairing, const K: usize> PCCommitmentState for CommitmentState<E, K> {
    type Randomness = ();

    fn empty() -> Self {
        Self::default()
    }

    fn rand<R: RngCore>(
        _num_queries: usize,
        _has_degree_bound: bool,
        _num_vars: Option<usize>,
        _rng: &mut R,
    ) -> Self::Randomness {
        // KZH is non-hiding, so its commitment randomness is the unit type.
    }
}

/// An opening proof for the KZH-`K` scheme.
///
/// The proof contains one vector of slice commitments for each of the lowest
/// `K - 1` tensor blocks, plus the final canonical arkworks evaluation table
/// for the highest block.
#[derive(Derivative, CanonicalSerialize, CanonicalDeserialize)]
#[derivative(
    Clone(bound = ""),
    Debug(bound = ""),
    PartialEq(bound = ""),
    Eq(bound = "")
)]
pub struct Proof<E: Pairing, const K: usize> {
    /// Runtime copy of the family parameter, included in the serialized form.
    pub(crate) k: usize,
    /// Number of variables in the opened polynomial.
    pub(crate) num_vars: usize,
    /// Slice commitments for every low-to-high tensor block except the last.
    pub(crate) layer_commitments: Vec<Vec<E::G1Affine>>,
    /// Canonical evaluations remaining for the highest/final block.
    pub(crate) final_evaluations: Vec<E::ScalarField>,
}

impl<E: Pairing, const K: usize> Proof<E, K> {
    /// Returns the number of variables opened by this proof.
    #[inline]
    pub const fn num_vars(&self) -> usize {
        self.num_vars
    }

    /// Returns the low-to-high slice-commitment layers in this proof.
    #[inline]
    pub fn layers(&self) -> &[Vec<E::G1Affine>] {
        &self.layer_commitments
    }

    /// Returns the disclosed evaluation table for the highest tensor block.
    #[inline]
    pub fn final_evaluations(&self) -> &[E::ScalarField] {
        &self.final_evaluations
    }
}

impl<E: Pairing, const K: usize> Default for Proof<E, K> {
    fn default() -> Self {
        Self {
            k: K,
            num_vars: 0,
            layer_commitments: Vec::new(),
            final_evaluations: Vec::new(),
        }
    }
}
