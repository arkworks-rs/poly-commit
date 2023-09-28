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
> {
    pub(crate) _field: PhantomData<F>,
    pub(crate) _config: PhantomData<C>,
    pub(crate) _digest: PhantomData<D>,
    pub(crate) _sponge: PhantomData<S>,
    pub(crate) _poly: PhantomData<P>,
}

#[derive(Derivative, CanonicalSerialize, CanonicalDeserialize)]
#[derivative(Clone(bound = ""), Debug(bound = ""))]
pub struct LigeroPCUniversalParams<F: PrimeField, C: Config>
where
    C: Config,
    <<C as Config>::TwoToOneHash as TwoToOneCRHScheme>::Parameters: Debug,
    <<C as Config>::LeafHash as CRHScheme>::Parameters: Debug,
{
    pub _field: PhantomData<F>,
    /// The security parameter
    pub sec_param: usize,
    /// The inverse of the code rate.
    pub rho_inv: usize,
    /// This is a flag which determines if the random linear combination is done.
    pub check_well_formedness: bool,
    /// Parameters for hash function of Merkle tree leaves
    #[derivative(Debug = "ignore")]
    pub leaf_hash_params: LeafParam<C>,
    /// Parameters for hash function of Merke tree combining two nodes into one
    #[derivative(Debug = "ignore")]
    pub two_to_one_params: TwoToOneParam<C>,
}

impl<F, C> PCUniversalParams for LigeroPCUniversalParams<F, C>
where
    F: PrimeField,
    C: Config,
    <<C as Config>::TwoToOneHash as TwoToOneCRHScheme>::Parameters: Debug,
    <<C as Config>::LeafHash as CRHScheme>::Parameters: Debug,
{
    fn max_degree(&self) -> usize {
        if F::TWO_ADICITY < self.rho_inv as u32 {
            0
        } else if (F::TWO_ADICITY - self.rho_inv as u32) * 2 < 64 {
            2_usize.pow((F::TWO_ADICITY - self.rho_inv as u32) * 2)
        } else {
            usize::MAX
        }
    }
}

/// Ligero commitment structure
#[derive(Derivative, CanonicalSerialize, CanonicalDeserialize)]
#[derivative(Clone(bound = ""), Debug(bound = ""))]
pub struct LigeroPCCommitterKey<F, C>
where
    F: PrimeField,
    C: Config,
    <<C as Config>::TwoToOneHash as TwoToOneCRHScheme>::Parameters: Debug,
    <<C as Config>::LeafHash as CRHScheme>::Parameters: Debug,
{
    pub(crate) _field: PhantomData<F>,
    /// The security parameter
    pub(crate) sec_param: usize,
    /// The inverse of code rate
    pub(crate) rho_inv: usize,
    /// Parameters for hash function of Merkle tree leaves
    #[derivative(Debug = "ignore")]
    pub(crate) leaf_hash_params: LeafParam<C>,
    /// Parameters for hash function of Merke tree combining two nodes into one
    #[derivative(Debug = "ignore")]
    pub(crate) two_to_one_params: TwoToOneParam<C>,
    /// This is a flag which determines if the random linear combination is done.
    pub(crate) check_well_formedness: bool,
}

impl<F, C> PCCommitterKey for LigeroPCCommitterKey<F, C>
where
    F: PrimeField,
    C: Config,
    <<C as Config>::TwoToOneHash as TwoToOneCRHScheme>::Parameters: Debug,
    <<C as Config>::LeafHash as CRHScheme>::Parameters: Debug,
{
    fn max_degree(&self) -> usize {
        if (F::TWO_ADICITY - self.rho_inv as u32) * 2 < 64 {
            2_usize.pow((F::TWO_ADICITY - self.rho_inv as u32) * 2)
        } else {
            usize::MAX
        }
    }

    fn supported_degree(&self) -> usize {
        self.max_degree()
    }
}
/// The verifier key which holds some scheme parameters
#[derive(Derivative, CanonicalSerialize, CanonicalDeserialize)]
#[derivative(Clone(bound = ""), Debug(bound = ""))]
pub struct LigeroPCVerifierKey<F, C>
where
    F: PrimeField,
    C: Config,
    <<C as Config>::TwoToOneHash as TwoToOneCRHScheme>::Parameters: Debug,
    <<C as Config>::LeafHash as CRHScheme>::Parameters: Debug,
{
    pub(crate) _field: PhantomData<F>,
    /// The security parameter
    pub(crate) sec_param: usize,
    /// The inverse of code rate
    pub(crate) rho_inv: usize,
    /// Parameters for hash function of Merkle tree leaves
    pub(crate) leaf_hash_params: LeafParam<C>,
    /// Parameters for hash function of Merke tree combining two nodes into one
    pub(crate) two_to_one_params: TwoToOneParam<C>,
    /// This is a flag which determines if the random linear combination is done.
    pub(crate) check_well_formedness: bool,
}

impl<F, C> PCVerifierKey for LigeroPCVerifierKey<F, C>
where
    F: PrimeField,
    C: Config,
    <<C as Config>::TwoToOneHash as TwoToOneCRHScheme>::Parameters: Debug,
    <<C as Config>::LeafHash as CRHScheme>::Parameters: Debug,
{
    fn max_degree(&self) -> usize {
        if (F::TWO_ADICITY - self.rho_inv as u32) * 2 < 64 {
            2_usize.pow((F::TWO_ADICITY - self.rho_inv as u32) * 2)
        } else {
            usize::MAX
        }
    }

    fn supported_degree(&self) -> usize {
        self.max_degree()
    }
}

pub(crate) type LigeroPCPreparedVerifierKey = ();

impl<Unprepared: PCVerifierKey> PCPreparedVerifierKey<Unprepared> for LigeroPCPreparedVerifierKey {
    fn prepare(_vk: &Unprepared) -> Self {
        todo!()
    }
}
#[derive(Derivative, CanonicalSerialize, CanonicalDeserialize)]
#[derivative(Default(bound = ""), Clone(bound = ""), Debug(bound = ""))]
pub(crate) struct Metadata {
    pub(crate) n_rows: usize,
    pub(crate) n_cols: usize,
    pub(crate) n_ext_cols: usize,
}

/// The commitment to a polynomial is a root of the merkle tree,
/// where each node is a hash of the column of the encoded coefficient matrix U.
#[derive(Derivative, CanonicalSerialize, CanonicalDeserialize)]
#[derivative(Default(bound = ""), Clone(bound = ""), Debug(bound = ""))]
pub struct LigeroPCCommitment<C: Config> {
    // number of rows resp. columns of the square matrix containing the coefficients of the polynomial
    pub(crate) metadata: Metadata,
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
