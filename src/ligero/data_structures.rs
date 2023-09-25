use ark_crypto_primitives::crh::{CRHScheme, TwoToOneCRHScheme};
use ark_crypto_primitives::{
    merkle_tree::{Config, LeafParam, Path, TwoToOneParam},
    sponge::CryptographicSponge,
};
use ark_ff::PrimeField;
use ark_poly::DenseUVPolynomial;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::fmt::Debug;
use core::marker::PhantomData;
use digest::Digest;

use crate::{
    PCCommitment, PCCommitterKey, PCPreparedCommitment, PCPreparedVerifierKey, PCRandomness,
    PCUniversalParams, PCVerifierKey,
};

use ark_std::rand::RngCore;

// TODO: Disclaimer: no hiding prop
/// The Ligero polynomial commitment scheme.
pub struct Ligero<
    F: PrimeField,
    C: Config,
    D: Digest,
    S: CryptographicSponge,
    P: DenseUVPolynomial<F>,
    // one over the rate rho
    const RHO_INV: usize,
    // security parameter, used in calculating t
    const SEC_PARAM: usize,
> {
    pub(crate) _field: PhantomData<F>,
    pub(crate) _config: PhantomData<C>,
    pub(crate) _digest: PhantomData<D>,
    pub(crate) _sponge: PhantomData<S>,
    pub(crate) _poly: PhantomData<P>,
}

// TODO come up with reasonable defaults
const DEFAULT_RHO_INV: usize = 2;
const DEFAULT_SEC_PARAM: usize = 128;

#[derive(Derivative, CanonicalSerialize, CanonicalDeserialize)]
#[derivative(Clone(bound = ""), Debug(bound = ""))]
pub struct LigeroPCUniversalParams {
    /// number of rows resp. columns of the square matrix containing the coefficients of the polynomial
    pub(crate) num_rows: usize,
    pub(crate) num_cols: usize,
    pub(crate) num_ext_cols: usize,
}

impl PCUniversalParams for LigeroPCUniversalParams {
    fn max_degree(&self) -> usize {
        todo!()
    }
}

/// Ligero commitment structure
#[derive(Derivative, CanonicalSerialize, CanonicalDeserialize)]
#[derivative(Clone(bound = ""), Debug(bound = ""))]
pub struct LigeroPCCommitterKey<C>
where
    C: Config,
    <<C as Config>::TwoToOneHash as TwoToOneCRHScheme>::Parameters: Debug,
    <<C as Config>::LeafHash as CRHScheme>::Parameters: Debug,
{
    #[derivative(Debug = "ignore")]
    /// Parameters for hash function of Merkle tree leaves
    pub leaf_hash_params: LeafParam<C>,
    #[derivative(Debug = "ignore")]
    /// Parameters for hash function of Merke tree combining two nodes into one
    pub two_to_one_params: TwoToOneParam<C>,
    /// This is a flag which determines if the random linear combination is done.
    pub check_well_formedness: bool,
}

impl<C> PCCommitterKey for LigeroPCCommitterKey<C>
where
    C: Config,
    <<C as Config>::TwoToOneHash as TwoToOneCRHScheme>::Parameters: Debug,
    <<C as Config>::LeafHash as CRHScheme>::Parameters: Debug,
{
    fn max_degree(&self) -> usize {
        todo!()
    }

    fn supported_degree(&self) -> usize {
        todo!()
    }
}
/// The verifier key which holds some scheme parameters
#[derive(Derivative, CanonicalSerialize, CanonicalDeserialize)]
#[derivative(Clone(bound = ""), Debug(bound = ""))]
pub struct LigeroPCVerifierKey<C>
where
    C: Config,
    <<C as Config>::TwoToOneHash as TwoToOneCRHScheme>::Parameters: Debug,
    <<C as Config>::LeafHash as CRHScheme>::Parameters: Debug,
{
    /// Parameters for hash function of Merkle tree leaves
    pub leaf_hash_params: LeafParam<C>,
    /// Parameters for hash function of Merke tree combining two nodes into one
    pub two_to_one_params: TwoToOneParam<C>,
    /// This is a flag which determines if the random linear combination is done.
    pub check_well_formedness: bool,
}

impl<C> PCVerifierKey for LigeroPCVerifierKey<C>
where
    C: Config,
    <<C as Config>::TwoToOneHash as TwoToOneCRHScheme>::Parameters: Debug,
    <<C as Config>::LeafHash as CRHScheme>::Parameters: Debug,
{
    fn max_degree(&self) -> usize {
        todo!()
    }

    fn supported_degree(&self) -> usize {
        todo!()
    }
}

pub(crate) type LigeroPCPreparedVerifierKey = ();

impl<Unprepared: PCVerifierKey> PCPreparedVerifierKey<Unprepared> for LigeroPCPreparedVerifierKey {
    fn prepare(_vk: &Unprepared) -> Self {
        todo!()
    }
}

/// The commitment to a polynomial is a root of the merkle tree,
/// where each node is a hash of the column of the encoded coefficient matrix U.
#[derive(Derivative, CanonicalSerialize, CanonicalDeserialize)]
#[derivative(Default(bound = ""), Clone(bound = ""), Debug(bound = ""))]
pub struct LigeroPCCommitment<C: Config> {
    // number of rows resp. columns of the square matrix containing the coefficients of the polynomial
    pub(crate) n_rows: usize,
    pub(crate) n_cols: usize,
    pub(crate) n_ext_cols: usize,
    pub(crate) root: C::InnerDigest,
}

impl<C: Config> PCCommitment for LigeroPCCommitment<C> {
    fn empty() -> Self {
        todo!()
    }

    fn has_degree_bound(&self) -> bool {
        todo!()
    }
}

pub(crate) type LigeroPCPreparedCommitment = ();

impl<Unprepared: PCCommitment> PCPreparedCommitment<Unprepared> for LigeroPCPreparedCommitment {
    fn prepare(_cm: &Unprepared) -> Self {
        todo!()
    }
}

pub(crate) type LigeroPCRandomness = ();

impl PCRandomness for LigeroPCRandomness {
    fn empty() -> Self {
        todo!()
    }

    fn rand<R: RngCore>(
        _num_queries: usize,
        _has_degree_bound: bool,
        _num_vars: Option<usize>,
        _rng: &mut R,
    ) -> Self {
        todo!()
    }
}

/// Proof of an individual Ligero well-formedness check or opening
#[derive(Derivative, CanonicalSerialize, CanonicalDeserialize)]
#[derivative(Default(bound = ""), Clone(bound = ""), Debug(bound = ""))]
pub(crate) struct LigeroPCProofSingle<F, C>
where
    F: PrimeField,
    C: Config,
{
    /// For each of the indices in q, `paths` contains the path from the root of the merkle tree to the leaf
    pub(crate) paths: Vec<Path<C>>,

    /// v, s.t. E(v) = w
    pub(crate) v: Vec<F>,

    pub(crate) columns: Vec<Vec<F>>,
}

/// The Proof type for Ligero, which amounts to an array of individual ligero proofs
#[derive(Derivative, CanonicalSerialize, CanonicalDeserialize)]
#[derivative(Default(bound = ""), Clone(bound = ""), Debug(bound = ""))]
pub struct LigeroPCProof<F, C>
where
    F: PrimeField,
    C: Config,
{
    pub(crate) opening: LigeroPCProofSingle<F, C>,
    pub(crate) well_formedness: Option<Vec<F>>,
}

// Multiple poly at one point
pub(crate) type LPCPArray<F, C> = Vec<LigeroPCProof<F, C>>;
