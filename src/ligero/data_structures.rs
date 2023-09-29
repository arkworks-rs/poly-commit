use crate::{
    Error, PCCommitment, PCCommitterKey, PCPreparedCommitment, PCPreparedVerifierKey, PCRandomness,
    PCUniversalParams, PCVerifierKey,
};
use ark_crypto_primitives::crh::{CRHScheme, TwoToOneCRHScheme};
use ark_crypto_primitives::merkle_tree::MerkleTree;
use ark_crypto_primitives::sponge::Absorb;
use ark_crypto_primitives::{
    merkle_tree::{Config, LeafParam, Path, TwoToOneParam},
    sponge::CryptographicSponge,
};
use ark_ff::PrimeField;
use ark_poly::DenseUVPolynomial;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::borrow::Borrow;
use ark_std::fmt::Debug;
use ark_std::marker::PhantomData;
use ark_std::rand::RngCore;
use digest::Digest;

use super::utils::Matrix;
use super::utils::{
    calculate_t, compute_dimensions, get_indices_from_transcript, hash_column, reed_solomon,
    IOPTranscript,
};

/// The univariate Ligero polynomial commitment scheme based on [[Ligero]][ligero].
/// The scheme defaults to the naive batching strategy.
///
/// Note: The scheme currently does not support hiding.
///
/// [ligero]: https://eprint.iacr.org/2022/1608.pdf
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

impl<F, C, D, S, P> Ligero<F, C, D, S, P>
where
    F: PrimeField,
    C: Config,
    Vec<u8>: Borrow<C::Leaf>,
    C::InnerDigest: Absorb,
    D: Digest,
    S: CryptographicSponge,
    P: DenseUVPolynomial<F>,
{
    /// Create a new instance of Ligero.
    pub fn new() -> Self {
        Self {
            _config: PhantomData,
            _field: PhantomData,
            // TODO potentially can get rid of digest and sponge
            _digest: PhantomData,
            _sponge: PhantomData,
            _poly: PhantomData,
        }
    }

    pub(crate) fn compute_matrices(polynomial: &P, rho_inv: usize) -> (Matrix<F>, Matrix<F>) {
        let mut coeffs = polynomial.coeffs().to_vec();

        // 1. Computing parameters and initial matrix
        let (n_rows, n_cols) = compute_dimensions::<F>(polynomial.degree() + 1); // for 6 coefficients, this is returning 4 x 2 with a row of 0s: fix

        // padding the coefficient vector with zeroes
        // TODO is this the most efficient/safest way to do it?
        coeffs.resize(n_rows * n_cols, F::zero());

        let mat = Matrix::new_from_flat(n_rows, n_cols, &coeffs);

        // 2. Apply Reed-Solomon encoding row-wise
        let ext_mat = Matrix::new_from_rows(
            mat.rows()
                .iter()
                .map(|r| reed_solomon(r, rho_inv))
                .collect(),
        );

        (mat, ext_mat)
    }
    pub(crate) fn create_merkle_tree(
        ext_mat: &Matrix<F>,
        leaf_hash_params: &<<C as Config>::LeafHash as CRHScheme>::Parameters,
        two_to_one_params: &<<C as Config>::TwoToOneHash as TwoToOneCRHScheme>::Parameters,
    ) -> MerkleTree<C> {
        let mut col_hashes: Vec<Vec<u8>> = Vec::new();
        let ext_mat_cols = ext_mat.cols();

        for col in ext_mat_cols.iter() {
            col_hashes.push(hash_column::<D, F>(col));
        }

        // pad the column hashes with zeroes
        let next_pow_of_two = col_hashes.len().next_power_of_two();
        col_hashes.resize(next_pow_of_two, vec![0; <D as Digest>::output_size()]);

        MerkleTree::<C>::new(leaf_hash_params, two_to_one_params, col_hashes).unwrap()
    }

    pub(crate) fn generate_proof(
        sec_param: usize,
        rho_inv: usize,
        b: &[F],
        mat: &Matrix<F>,
        ext_mat: &Matrix<F>,
        col_tree: &MerkleTree<C>,
        transcript: &mut IOPTranscript<F>,
    ) -> Result<LigeroPCProofSingle<F, C>, Error> {
        let t = calculate_t::<F>(sec_param, rho_inv, ext_mat.n)?;

        // 1. left-multiply the matrix by `b`, where for a requested query point `z`,
        // `b = [1, z^m, z^(2m), ..., z^((m-1)m)]`
        let v = mat.row_mul(b);

        transcript
            .append_serializable_element(b"v", &v)
            .map_err(|_| Error::TranscriptError)?;

        // 2. Generate t column indices to test the linear combination on
        let indices = get_indices_from_transcript(ext_mat.m, t, transcript)?;

        // 3. Compute Merkle tree paths for the requested columns
        let mut queried_columns = Vec::with_capacity(t);
        let mut paths = Vec::with_capacity(t);

        let ext_mat_cols = ext_mat.cols();

        for i in indices {
            queried_columns.push(ext_mat_cols[i].clone());
            paths.push(
                col_tree
                    .generate_proof(i)
                    .map_err(|_| Error::TranscriptError)?,
            );
        }

        Ok(LigeroPCProofSingle {
            paths,
            v,
            columns: queried_columns,
        })
    }
}

/// The public parameters for the Ligero polynomial commitment scheme.
/// This is only a default setup with reasonable parameters.
/// To create your own public parameters, use:
/// # Example
/// ```rust
/// use ark_bls12_377::Fr;
/// use ark_crypto_primitives::{
///     crh::{pedersen, sha256::Sha256, CRHScheme, TwoToOneCRHScheme},
///     merkle_tree::{ByteDigestConverter, Config},
/// };
/// use ark_std::test_rng;
/// use ark_poly_commit::ligero::LigeroPCUniversalParams;
/// use core::marker::PhantomData;
///
/// type LeafH = Sha256;
/// type CompressH = Sha256;
/// struct MerkleTreeParams;
/// impl Config for MerkleTreeParams {
///     type Leaf = [u8];
///     type LeafDigest = <LeafH as CRHScheme>::Output;
///     type LeafInnerDigestConverter = ByteDigestConverter<Self::LeafDigest>;
///     type InnerDigest = <CompressH as TwoToOneCRHScheme>::Output;
///     type LeafHash = LeafH;
///     type TwoToOneHash = CompressH;
/// }
/// type MTConfig = MerkleTreeParams;
/// let mut rng = &mut test_rng();
/// let leaf_hash_params = <LeafH as CRHScheme>::setup(&mut rng).unwrap();
/// let two_to_one_params = <CompressH as TwoToOneCRHScheme>::setup(&mut rng)
///     .unwrap()
///     .clone();
/// let pp: LigeroPCUniversalParams<Fr, MTConfig> = LigeroPCUniversalParams::new(128, 2, true,
///     leaf_hash_params, two_to_one_params);
#[derive(Derivative, CanonicalSerialize, CanonicalDeserialize)]
#[derivative(Clone(bound = ""), Debug(bound = ""))]
pub struct LigeroPCUniversalParams<F: PrimeField, C: Config>
where
    C: Config,
    <<C as Config>::TwoToOneHash as TwoToOneCRHScheme>::Parameters: Debug,
    <<C as Config>::LeafHash as CRHScheme>::Parameters: Debug,
{
    _field: PhantomData<F>,
    /// The security parameter
    pub(crate) sec_param: usize,
    /// The inverse of the code rate.
    pub(crate) rho_inv: usize,
    /// This is a flag which determines if the random linear combination is done.
    pub(crate) check_well_formedness: bool,
    /// Parameters for hash function of Merkle tree leaves
    #[derivative(Debug = "ignore")]
    pub(crate) leaf_hash_params: LeafParam<C>,
    /// Parameters for hash function of Merke tree combining two nodes into one
    #[derivative(Debug = "ignore")]
    pub(crate) two_to_one_params: TwoToOneParam<C>,
}

impl<F, C> LigeroPCUniversalParams<F, C>
where
    F: PrimeField,
    C: Config,
    <<C as Config>::TwoToOneHash as TwoToOneCRHScheme>::Parameters: Debug,
    <<C as Config>::LeafHash as CRHScheme>::Parameters: Debug,
{
    /// Create new LigeroPCUniversalParams
    pub fn new(
        sec_param: usize,
        rho_inv: usize,
        check_well_formedness: bool,
        leaf_hash_params: LeafParam<C>,
        two_to_one_params: TwoToOneParam<C>,
    ) -> Self {
        Self {
            _field: PhantomData,
            sec_param,
            rho_inv,
            check_well_formedness,
            leaf_hash_params,
            two_to_one_params,
        }
    }
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
    fn prepare(_vk: &Unprepared) -> Self {}
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
        LigeroPCCommitment::default()
    }

    fn has_degree_bound(&self) -> bool {
        false
    }
}

pub(crate) type LigeroPCPreparedCommitment = ();

impl<Unprepared: PCCommitment> PCPreparedCommitment<Unprepared> for LigeroPCPreparedCommitment {
    fn prepare(_cm: &Unprepared) -> Self {}
}

pub(crate) type LigeroPCRandomness = ();

impl PCRandomness for LigeroPCRandomness {
    fn empty() -> Self {}

    fn rand<R: RngCore>(
        _num_queries: usize,
        _has_degree_bound: bool,
        _num_vars: Option<usize>,
        _rng: &mut R,
    ) -> Self {
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
