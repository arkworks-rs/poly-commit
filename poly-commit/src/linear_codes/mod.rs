use crate::utils::{inner_product, Matrix};
use crate::{
    to_bytes, Error, LabeledCommitment, LabeledPolynomial, PCCommitterKey, PCUniversalParams,
    PCVerifierKey, PolynomialCommitment,
};

use ark_crypto_primitives::crh::{CRHScheme, TwoToOneCRHScheme};
use ark_crypto_primitives::merkle_tree::MerkleTree;
use ark_crypto_primitives::{
    merkle_tree::Config,
    sponge::{Absorb, CryptographicSponge},
};
use ark_ff::PrimeField;
use ark_poly::Polynomial;
use ark_std::borrow::Borrow;
use ark_std::marker::PhantomData;
use ark_std::rand::RngCore;
use ark_std::string::ToString;
use ark_std::vec::Vec;

#[cfg(feature = "parallel")]
use rayon::iter::{IntoParallelIterator, IntoParallelRefIterator, ParallelIterator};

mod utils;

mod multilinear_brakedown;

pub use multilinear_brakedown::MultilinearBrakedown;

mod brakedown;
mod data_structures;
use data_structures::*;

pub use data_structures::LinCodePCProof;

use utils::{calculate_t, get_indices_from_sponge};

const FIELD_SIZE_ERROR: &str = "This field is not suitable for the proposed parameters";

/// For linear code PC schemes, the universal paramters, committer key
/// and verifier key are all the same. This trait abstracts the common
/// information contained in these.
pub trait LinCodeParametersInfo<C, H>
where
    C: Config,
    H: CRHScheme,
{
    /// Get the security parameter.
    fn sec_param(&self) -> usize;

    /// Get the distance of the code.
    fn distance(&self) -> (usize, usize);

    /// See whether there should be a well-formedness check.
    fn check_well_formedness(&self) -> bool;

    /// Compute the dimensions of the coefficient matrix.
    fn compute_dimensions(&self, n: usize) -> (usize, usize);

    /// Get the hash parameters for obtaining leaf digest from leaf value.
    fn leaf_hash_param(&self) -> &<<C as Config>::LeafHash as CRHScheme>::Parameters;

    /// Get the parameters for hashing nodes in the merkle tree.
    fn two_to_one_hash_param(
        &self,
    ) -> &<<C as Config>::TwoToOneHash as TwoToOneCRHScheme>::Parameters;

    /// Get the parameters for hashing a vector of values,
    /// representing a column of the coefficient matrix, into a leaf value.
    fn col_hash_params(&self) -> &H::Parameters;
}

/// A trait for linear codes.
pub trait LinearEncode<F, C, P, H>
where
    F: PrimeField,
    C: Config,
    H: CRHScheme,
    P: Polynomial<F>,
{
    /// For schemes like Brakedown and Ligero, PCCommiiterKey and
    /// PCVerifierKey and PCUniversalParams are all the same.
    type LinCodePCParams: PCUniversalParams
        + PCCommitterKey
        + PCVerifierKey
        + LinCodeParametersInfo<C, H>;

    /// Does a default setup for the PCS.
    fn setup<R: RngCore>(
        max_degree: usize,
        num_vars: Option<usize>,
        rng: &mut R,
        leaf_hash_param: <<C as Config>::LeafHash as CRHScheme>::Parameters,
        two_to_one_hash_param: <<C as Config>::TwoToOneHash as TwoToOneCRHScheme>::Parameters,
        col_hash_params: H::Parameters,
    ) -> Self::LinCodePCParams;

    /// Encode a message, which is interpreted as a vector of coefficients
    /// of a polynomial of degree m - 1.
    fn encode(msg: &[F], param: &Self::LinCodePCParams) -> Result<Vec<F>, Error>;

    /// Represent the polynomial as either coefficients,
    /// in the univariate case, or evaluations over
    /// the Boolean hypercube, in the multilinear case.
    fn poly_to_vec(polynomial: &P) -> Vec<F>;

    /// Represent the query point as a vector of Field elements.
    fn point_to_vec(point: P::Point) -> Vec<F>;

    /// Arrange the coefficients of the polynomial into a matrix,
    /// and apply encoding to each row.
    /// Returns the tuple (original_matrix, encoded_matrix).
    fn compute_matrices(polynomial: &P, param: &Self::LinCodePCParams) -> (Matrix<F>, Matrix<F>) {
        let mut coeffs = Self::poly_to_vec(polynomial);

        // 1. Computing the matrix dimensions.
        let (n_rows, n_cols) = param.compute_dimensions(coeffs.len());

        // padding the coefficient vector with zeroes
        coeffs.resize(n_rows * n_cols, F::zero());

        let mat = Matrix::new_from_flat(n_rows, n_cols, &coeffs);

        // 2. Apply encoding row-wise
        let rows = mat.rows();
        let ext_mat = Matrix::new_from_rows(
            cfg_iter!(rows)
                .map(|r| Self::encode(r, param).unwrap())
                .collect(),
        );

        (mat, ext_mat)
    }

    /// Tensor the query point z in the following sense:
    /// For a polynomial p(X) represented by a matrix M
    /// with n rows and m columns such that M_{i,j} = p_{i + n*j},
    /// we define the tensoring of `z`: (a, b) = tensor(z, n, m) such that:
    /// p(z) = b^T.M.a
    /// returns the evaluation of p at z.
    fn tensor(z: &P::Point, n: usize, m: usize) -> (Vec<F>, Vec<F>);
}

/// Any linear-code-based commitment scheme.
pub struct LinearCodePCS<L, F, P, S, C, H>
where
    F: PrimeField,
    C: Config,
    S: CryptographicSponge,
    P: Polynomial<F>,
    H: CRHScheme,
    L: LinearEncode<F, C, P, H>,
{
    _phantom: PhantomData<(L, F, P, S, C, H)>,
}

impl<L, F, P, S, C, H> PolynomialCommitment<F, P, S> for LinearCodePCS<L, F, P, S, C, H>
where
    L: LinearEncode<F, C, P, H>,
    F: PrimeField + Absorb,
    P: Polynomial<F>,
    S: CryptographicSponge,
    C: Config + 'static,
    Vec<F>: Borrow<<H as CRHScheme>::Input>,
    H::Output: Into<C::Leaf> + Send,
    C::Leaf: Sized + Clone + Default + Send + AsRef<C::Leaf>,
    H: CRHScheme + 'static,
{
    type UniversalParams = L::LinCodePCParams;

    type CommitterKey = L::LinCodePCParams;

    type VerifierKey = L::LinCodePCParams;

    type Commitment = LinCodePCCommitment<C>;

    type CommitmentState = LinCodePCCommitmentState<F, H>;

    type Proof = LPCPArray<F, C>;

    type BatchProof = Vec<Self::Proof>;

    type Error = Error;

    /// This is only a default setup with reasonable parameters.
    /// To create your own public parameters (from which vk/ck can be derived by `trim`),
    /// see the documentation for `BrakedownPCUniversalParams`.
    fn setup<R: RngCore>(
        max_degree: usize,
        num_vars: Option<usize>,
        rng: &mut R,
    ) -> Result<Self::UniversalParams, Self::Error> {
        let leaf_hash_param = <C::LeafHash as CRHScheme>::setup(rng).unwrap();
        let two_to_one_hash_param = <C::TwoToOneHash as TwoToOneCRHScheme>::setup(rng)
            .unwrap()
            .clone();
        let col_hash_params = <H as CRHScheme>::setup(rng).unwrap();
        let pp = L::setup::<R>(
            max_degree,
            num_vars,
            rng,
            leaf_hash_param,
            two_to_one_hash_param,
            col_hash_params,
        );
        let real_max_degree = <Self::UniversalParams as PCUniversalParams>::max_degree(&pp);
        if max_degree > real_max_degree || real_max_degree == 0 {
            return Err(Error::InvalidParameters(FIELD_SIZE_ERROR.to_string()));
        }
        Ok(pp)
    }

    fn trim(
        pp: &Self::UniversalParams,
        _supported_degree: usize,
        _supported_hiding_bound: usize,
        _enforced_degree_bounds: Option<&[usize]>,
    ) -> Result<(Self::CommitterKey, Self::VerifierKey), Self::Error> {
        if <Self::UniversalParams as PCUniversalParams>::max_degree(pp) == 0 {
            return Err(Error::InvalidParameters(FIELD_SIZE_ERROR.to_string()));
        }
        Ok((pp.clone(), pp.clone()))
    }

    fn commit<'a>(
        ck: &Self::CommitterKey,
        polynomials: impl IntoIterator<Item = &'a LabeledPolynomial<F, P>>,
        _rng: Option<&mut dyn RngCore>,
    ) -> Result<
        (
            Vec<LabeledCommitment<Self::Commitment>>,
            Vec<Self::CommitmentState>,
        ),
        Self::Error,
    >
    where
        P: 'a,
    {
        let mut commitments = Vec::new();
        let mut states = Vec::new();

        for labeled_polynomial in polynomials {
            let polynomial = labeled_polynomial.polynomial();

            // 1. Arrange the coefficients of the polynomial into a matrix,
            // and apply encoding to get `ext_mat`.
            let (mat, ext_mat) = L::compute_matrices(polynomial, ck);
            let n_rows = mat.n;
            let n_cols = mat.m;
            let n_ext_cols = ext_mat.m;

            // 2. Create the Merkle tree from the hashes of each column.
            let ext_mat_cols = ext_mat.cols();
            let leaves: Vec<H::Output> = cfg_into_iter!(ext_mat_cols)
                .map(|col| {
                    H::evaluate(ck.col_hash_params(), col)
                        .map_err(|_| Error::HashingError)
                        .unwrap()
                })
                .collect();
            let state = Self::CommitmentState {
                mat,
                ext_mat,
                leaves,
            };
            let mut leaves: Vec<C::Leaf> =
                state.leaves.clone().into_iter().map(|h| h.into()).collect();
            let col_tree = create_merkle_tree::<C>(
                &mut leaves,
                ck.leaf_hash_param(),
                ck.two_to_one_hash_param(),
            )?;

            // 3. Obtain the MT root
            let root = col_tree.root();

            // 4. The commitment is just the root, but since each commitment could be to a differently-sized polynomial, we also add some metadata.
            let commitment = LinCodePCCommitment {
                metadata: Metadata {
                    n_rows,
                    n_cols,
                    n_ext_cols,
                },
                root,
            };

            commitments.push(LabeledCommitment::new(
                labeled_polynomial.label().clone(),
                commitment,
                None,
            ));
            states.push(state);
        }
        Ok((commitments, states))
    }

    fn open<'a>(
        ck: &Self::CommitterKey,
        _labeled_polynomials: impl IntoIterator<Item = &'a LabeledPolynomial<F, P>>,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Self::Commitment>>,
        point: &'a P::Point,
        sponge: &mut S,
        states: impl IntoIterator<Item = &'a Self::CommitmentState>,
        _rng: Option<&mut dyn RngCore>,
    ) -> Result<Self::Proof, Self::Error>
    where
        P: 'a,
        Self::CommitmentState: 'a,
        Self::Commitment: 'a,
    {
        let mut proof_array = LPCPArray::default();

        for (labeled_commitment, state) in commitments.into_iter().zip(states) {
            let commitment = labeled_commitment.commitment();
            let n_rows = commitment.metadata.n_rows;
            let n_cols = commitment.metadata.n_cols;

            // 1. Arrange the coefficients of the polynomial into a matrix,
            // and apply encoding to get `ext_mat`.
            // 2. Create the Merkle tree from the hashes of each column.
            let Self::CommitmentState {
                mat,
                ext_mat,
                leaves: col_hashes,
            } = state;
            let mut col_hashes: Vec<C::Leaf> =
                col_hashes.clone().into_iter().map(|h| h.into()).collect();

            let col_tree = create_merkle_tree::<C>(
                &mut col_hashes,
                ck.leaf_hash_param(),
                ck.two_to_one_hash_param(),
            )?;

            // 3. Generate vector `b` to left-multiply the matrix.
            let (_, b) = L::tensor(point, n_cols, n_rows);

            sponge.absorb(&to_bytes!(&commitment.root).map_err(|_| Error::TranscriptError)?);

            // If we are checking well-formedness, we need to compute the well-formedness proof (which is just r.M) and append it to the transcript.
            let well_formedness = if ck.check_well_formedness() {
                let r = sponge.squeeze_field_elements::<F>(n_rows);
                let v = mat.row_mul(&r);

                sponge.absorb(&v);
                Some(v)
            } else {
                None
            };

            let point_vec = L::point_to_vec(point.clone());
            sponge.absorb(&point_vec);

            proof_array.push(LinCodePCProof {
                // Compute the opening proof and append b.M to the transcript.
                opening: generate_proof(
                    ck.sec_param(),
                    ck.distance(),
                    &b,
                    &mat,
                    &ext_mat,
                    &col_tree,
                    sponge,
                )?,
                well_formedness,
            });
        }

        Ok(proof_array)
    }

    fn check<'a>(
        vk: &Self::VerifierKey,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Self::Commitment>>,
        point: &'a P::Point,
        values: impl IntoIterator<Item = F>,
        proof_array: &Self::Proof,
        sponge: &mut S,
        _rng: Option<&mut dyn RngCore>,
    ) -> Result<bool, Self::Error>
    where
        Self::Commitment: 'a,
    {
        let leaf_hash_param: &<<C as Config>::LeafHash as CRHScheme>::Parameters =
            vk.leaf_hash_param();
        let two_to_one_hash_param: &<<C as Config>::TwoToOneHash as TwoToOneCRHScheme>::Parameters =
            vk.two_to_one_hash_param();

        for (i, (labeled_commitment, value)) in commitments.into_iter().zip(values).enumerate() {
            let proof = &proof_array[i];
            let commitment = labeled_commitment.commitment();
            let n_rows = commitment.metadata.n_rows;
            let n_cols = commitment.metadata.n_cols;
            let n_ext_cols = commitment.metadata.n_ext_cols;
            let root = &commitment.root;
            let t = calculate_t::<F>(vk.sec_param(), vk.distance(), n_ext_cols)?;

            sponge.absorb(&to_bytes!(&commitment.root).map_err(|_| Error::TranscriptError)?);

            let out = if vk.check_well_formedness() {
                if proof.well_formedness.is_none() {
                    return Err(Error::InvalidCommitment);
                }
                let tmp = &proof.well_formedness.as_ref();
                let v = tmp.unwrap();
                let r = sponge.squeeze_field_elements::<F>(n_rows);
                // Upon sending `v` to the Verifier, add it to the sponge. The claim is that v = r.M.
                sponge.absorb(&v);

                (Some(v), Some(r))
            } else {
                (None, None)
            };

            // 1. Seed the transcript with the point and the recieved vector
            // TODO Consider removing the evaluation point from the transcript.
            let point_vec = L::point_to_vec(point.clone());
            sponge.absorb(&point_vec);
            sponge.absorb(&proof.opening.v);

            // 2. Ask random oracle for the `t` indices where the checks happen.
            let indices = get_indices_from_sponge(n_ext_cols, t, sponge)?;

            // 3. Hash the received columns into leaf hashes.
            let col_hashes: Vec<C::Leaf> = proof
                .opening
                .columns
                .iter()
                .map(|c| {
                    H::evaluate(vk.col_hash_params(), c.clone())
                        .map_err(|_| Error::HashingError)
                        .unwrap()
                        .into()
                })
                .collect();

            // 4. Verify the paths for each of the leaf hashes - this is only run once,
            // even if we have a well-formedness check (i.e., we save sending and checking the columns).
            // See "Concrete optimizations to the commitment scheme", p.12 of [Brakedown](https://eprint.iacr.org/2021/1043.pdf).
            for (j, (leaf, q_j)) in col_hashes.iter().zip(indices.iter()).enumerate() {
                let path = &proof.opening.paths[j];
                if path.leaf_index != *q_j {
                    return Err(Error::InvalidCommitment);
                }

                path.verify(leaf_hash_param, two_to_one_hash_param, root, leaf.clone())
                    .map_err(|_| Error::InvalidCommitment)?;
            }

            // Helper closure: checks if a.b = c.
            let check_inner_product = |a, b, c| -> Result<(), Error> {
                if inner_product(a, b) != c {
                    return Err(Error::InvalidCommitment);
                }

                Ok(())
            };

            // 5. Compute the encoding w = E(v).
            let w = L::encode(&proof.opening.v, vk)?;

            // 6. Compute `a`, `b` to right- and left- multiply with the matrix `M`.
            let (a, b) = L::tensor(point, n_cols, n_rows);

            // 7. Probabilistic checks that whatever the prover sent,
            // matches with what the verifier computed for himself.
            // Note: we sacrifice some code repetition in order not to repeat execution.
            if let (Some(well_formedness), Some(r)) = out {
                let w_well_formedness = L::encode(well_formedness, vk)?;
                for (transcript_index, matrix_index) in indices.iter().enumerate() {
                    check_inner_product(
                        &r,
                        &proof.opening.columns[transcript_index],
                        w_well_formedness[*matrix_index],
                    )?;
                    check_inner_product(
                        &b,
                        &proof.opening.columns[transcript_index],
                        w[*matrix_index],
                    )?;
                }
            } else {
                for (transcript_index, matrix_index) in indices.iter().enumerate() {
                    check_inner_product(
                        &b,
                        &proof.opening.columns[transcript_index],
                        w[*matrix_index],
                    )?;
                }
            }

            if inner_product(&proof.opening.v, &a) != value {
                eprintln!("Function check: claimed value in position {i} does not match the evaluation of the committed polynomial in the same position");
                return Ok(false);
            }
        }

        Ok(true)
    }
}

// TODO maybe this can go to utils
fn create_merkle_tree<C>(
    leaves: &mut Vec<C::Leaf>,
    leaf_hash_param: &<<C as Config>::LeafHash as CRHScheme>::Parameters,
    two_to_one_hash_param: &<<C as Config>::TwoToOneHash as TwoToOneCRHScheme>::Parameters,
) -> Result<MerkleTree<C>, Error>
where
    C: Config,
    C::Leaf: Default + Clone + Send + AsRef<C::Leaf>,
{
    // pad the column hashes with zeroes
    let next_pow_of_two = leaves.len().next_power_of_two();
    leaves.resize(next_pow_of_two, <C::Leaf>::default());

    MerkleTree::<C>::new(leaf_hash_param, two_to_one_hash_param, leaves)
        .map_err(|_| Error::HashingError)
}

fn generate_proof<F, C, S>(
    sec_param: usize,
    distance: (usize, usize),
    b: &[F],
    mat: &Matrix<F>,
    ext_mat: &Matrix<F>,
    col_tree: &MerkleTree<C>,
    sponge: &mut S,
) -> Result<LinCodePCProofSingle<F, C>, Error>
where
    F: PrimeField + Absorb,
    C: Config,
    S: CryptographicSponge,
{
    let t = calculate_t::<F>(sec_param, distance, ext_mat.m)?;

    // 1. left-multiply the matrix by `b`.
    let v = mat.row_mul(b);
    sponge.absorb(&v);

    // 2. Generate t column indices to test the linear combination on.
    let indices = get_indices_from_sponge(ext_mat.m, t, sponge)?;

    // 3. Compute Merkle tree paths for the requested columns.
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

    Ok(LinCodePCProofSingle {
        paths,
        v,
        columns: queried_columns,
    })
}
