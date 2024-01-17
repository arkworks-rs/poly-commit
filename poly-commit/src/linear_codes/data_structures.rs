use super::utils::SprsMat;
use crate::{utils::Matrix, PCCommitment, PCCommitmentState};
use ark_crypto_primitives::{
    crh::CRHScheme,
    merkle_tree::{Config, LeafParam, Path, TwoToOneParam},
};
use ark_ff::PrimeField;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::rand::RngCore;
use ark_std::vec::Vec;

#[derive(Derivative, CanonicalSerialize, CanonicalDeserialize)]
#[derivative(Clone(bound = ""), Debug(bound = ""))]
/// The public parameters for Brakedown PCS.
pub struct BrakedownPCParams<F: PrimeField, C: Config, H: CRHScheme> {
    /// The security parameter
    pub(crate) sec_param: usize,
    /// alpha in the paper
    pub(crate) alpha: (usize, usize),
    /// beta in the paper
    pub(crate) beta: (usize, usize),
    /// The inverse of the code rate.
    pub(crate) rho_inv: (usize, usize),
    /// Threshold of the base case to encode with RS
    pub(crate) base_len: usize,
    /// Length of each column in the matrix that represents the polynomials
    pub(crate) n: usize,
    /// Length of each row in the matrix that represents the polynomials
    pub(crate) m: usize,
    /// Length of each row in the matrix that represents the polynomials, **after encoding**
    pub(crate) m_ext: usize,
    /// Constarints on A matrices. `a_dims[i]` is `(n, m, c)`, where `n` is
    /// the number of rows, `m` is the number of columns, `c` is the number of
    /// non-zero elements in each row, for the matrix A in the `i`th step of
    /// the encoding.
    pub(crate) a_dims: Vec<(usize, usize, usize)>,
    /// Same as `a_dims`, but for B matrices.
    pub(crate) b_dims: Vec<(usize, usize, usize)>,
    /// By having `a_dims` and `b_dims`, we compute a vector of indices that
    /// specfies where is the beginning of the sub-chunk that we need to
    /// encode during the recursive encoding. Notice that we do not recurse
    /// in this implementation, instead we do it iteratively.
    pub(crate) start: Vec<usize>,
    /// Same as `start`, but stores the end index of those chunks.
    pub(crate) end: Vec<usize>,
    /// A vector of all A matrices we need for encoding.
    pub(crate) a_mats: Vec<SprsMat<F>>,
    /// A vector of all B matrices we need for encoding.
    pub(crate) b_mats: Vec<SprsMat<F>>,
    /// This is a flag which determines if the random linear combination is done.
    pub(crate) check_well_formedness: bool,
    /// Parameters for hash function of Merkle tree leaves
    #[derivative(Debug = "ignore")]
    pub(crate) leaf_hash_param: LeafParam<C>,
    /// Parameters for hash function of Merke tree combining two nodes into one
    #[derivative(Debug = "ignore")]
    pub(crate) two_to_one_hash_param: TwoToOneParam<C>,
    // Parameters for obtaining leaf digest from leaf value.
    #[derivative(Debug = "ignore")]
    pub(crate) col_hash_params: H::Parameters,
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
pub struct LinCodePCCommitment<C: Config> {
    // number of rows resp. columns of the square matrix containing the coefficients of the polynomial
    pub(crate) metadata: Metadata,
    pub(crate) root: C::InnerDigest,
}

impl<C: Config> PCCommitment for LinCodePCCommitment<C> {
    fn empty() -> Self {
        LinCodePCCommitment::default()
    }

    fn has_degree_bound(&self) -> bool {
        false
    }
}

#[derive(Derivative, CanonicalSerialize, CanonicalDeserialize)]
#[derivative(Default(bound = ""), Clone(bound = ""), Debug(bound = ""))]
pub struct LinCodePCCommitmentState<F, H>
where
    F: PrimeField,
    H: CRHScheme,
{
    pub(crate) mat: Matrix<F>,
    pub(crate) ext_mat: Matrix<F>,
    pub(crate) leaves: Vec<H::Output>,
}

impl<F, H> PCCommitmentState for LinCodePCCommitmentState<F, H>
where
    F: PrimeField,
    H: CRHScheme,
{
    type Randomness = ();
    fn empty() -> Self {
        unimplemented!()
    }

    fn rand<R: RngCore>(
        _num_queries: usize,
        _has_degree_bound: bool,
        _num_vars: Option<usize>,
        _rng: &mut R,
    ) -> Self::Randomness {
        unimplemented!()
    }
}

/// Proof of an individual linear code well-formedness check or opening
#[derive(Derivative, CanonicalSerialize, CanonicalDeserialize)]
#[derivative(Default(bound = ""), Clone(bound = ""), Debug(bound = ""))]
pub(crate) struct LinCodePCProofSingle<F, C>
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

/// The Proof type for linear code PCS, which amounts to an array of individual proofs
#[derive(Derivative, CanonicalSerialize, CanonicalDeserialize)]
#[derivative(Default(bound = ""), Clone(bound = ""), Debug(bound = ""))]
pub struct LinCodePCProof<F, C>
where
    F: PrimeField,
    C: Config,
{
    pub(crate) opening: LinCodePCProofSingle<F, C>,
    pub(crate) well_formedness: Option<Vec<F>>,
}

// Multiple poly at one point
pub(crate) type LPCPArray<F, C> = Vec<LinCodePCProof<F, C>>;
