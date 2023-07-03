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
    const rho_inv: usize,
    // security parameter, used in calculating t
    const sec_param: usize,
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

pub(crate) type LigeroPCUniversalParams = ();

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
    fn prepare(vk: &Unprepared) -> Self {
        todo!()
    }
}

/// The commitment to a polynomial is a root of the merkle tree,
/// where each node is a hash of the column of the encoded coefficient matrix U.
#[derive(Derivative, CanonicalSerialize, CanonicalDeserialize)]
#[derivative(Default(bound = ""), Clone(bound = ""), Debug(bound = ""))]
pub struct LigeroPCCommitment<F: PrimeField, C: Config> {
    // number of rows of the square matrix containing the coefficients of the polynomial
    pub(crate) m: usize,
    // TODO is InnerDigest the right type?
    pub(crate) root: C::InnerDigest,
    pub(crate) proof: LigeroPCProof<F, C>,
}

impl<F: PrimeField, C: Config> PCCommitment for LigeroPCCommitment<F, C> {
    fn empty() -> Self {
        todo!()
    }

    fn has_degree_bound(&self) -> bool {
        todo!()
    }
}

pub(crate) type LigeroPCPreparedCommitment = ();

impl<Unprepared: PCCommitment> PCPreparedCommitment<Unprepared> for LigeroPCPreparedCommitment {
    fn prepare(cm: &Unprepared) -> Self {
        todo!()
    }
}

pub(crate) type LigeroPCRandomness = ();

impl PCRandomness for LigeroPCRandomness {
    fn empty() -> Self {
        todo!()
    }

    fn rand<R: RngCore>(
        num_queries: usize,
        has_degree_bound: bool,
        num_vars: Option<usize>,
        rng: &mut R,
    ) -> Self {
        todo!()
    }
}

/// Ligero proof
#[derive(Derivative, CanonicalSerialize, CanonicalDeserialize)]
#[derivative(Default(bound = ""), Clone(bound = ""), Debug(bound = ""))]
pub struct LigeroPCProof<F, C>
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
