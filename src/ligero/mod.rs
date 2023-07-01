use ark_crypto_primitives::crh::{CRHScheme, TwoToOneCRHScheme};
use ark_crypto_primitives::{
    merkle_tree::{Config, LeafParam, MerkleTree, TwoToOneParam},
    sponge::{Absorb, CryptographicSponge},
};
use ark_ff::PrimeField;
use ark_poly::{DenseUVPolynomial, EvaluationDomain, GeneralEvaluationDomain};
use ark_std::fmt::Debug;
use core::marker::PhantomData;
use digest::Digest;
use jf_primitives::pcs::transcript::IOPTranscript;
use std::borrow::Borrow;

use crate::{LabeledPolynomial, PolynomialCommitment, Error, LabeledCommitment};
use crate::{
    ligero::utils::{inner_product, reed_solomon},
};

use ark_std::rand::RngCore;

mod utils;
use utils::Matrix;

mod data_structures;
use data_structures::*;

pub use data_structures::{Ligero, LigeroPCCommitterKey, LigeroPCVerifierKey};

use utils::{get_indices_from_transcript, calculate_t, hash_column};
mod tests;

impl<F, C, D, S, P, const rho_inv: usize, const sec_param: usize>
    Ligero<F, C, D, S, P, rho_inv, sec_param>
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
    /// If either or both parameters are None, their default values are used.
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

    /// The verifier can check the well-formedness of the commitment by taking random linear combinations.
    fn check_well_formedness(
        commitment: &LigeroPCCommitment<F, C>,
        leaf_hash_params: &LeafParam<C>,
        two_to_one_params: &TwoToOneParam<C>,
    ) -> Result<(), Error> {
        let t = calculate_t(rho_inv, sec_param);

        // TODO replace unwraps by proper error handling
        let mut transcript: IOPTranscript<F> = IOPTranscript::new(b"well_formedness_transcript");
        transcript
            .append_serializable_element(b"root", &commitment.root)
            .unwrap();

        // 2. Get the linear combination coefficients from the transcript
        let mut r = Vec::new();
        for _ in 0..t {
            r.push(transcript.get_and_append_challenge(b"r").unwrap());
        }
        // Upon sending `v` to the Verifier, add it to the sponge. Claim is that v = r.M
        transcript
            .append_serializable_element(b"v", &commitment.proof.v)
            .unwrap();

        Self::check_random_linear_combination(
            &r,
            commitment,
            t,
            &mut transcript,
            leaf_hash_params,
            two_to_one_params,
        )
    }
    fn check_random_linear_combination(
        coeffs: &[F],
        commitment: &LigeroPCCommitment<F, C>,
        t: usize,
        transcript: &mut IOPTranscript<F>,
        leaf_hash_params: &<<C as Config>::LeafHash as CRHScheme>::Parameters,
        two_to_one_params: &<<C as Config>::TwoToOneHash as TwoToOneCRHScheme>::Parameters,
    ) -> Result<(), Error> {
        // 1. Hash the received columns into leaf hashes
        let mut col_hashes: Vec<Vec<u8>> = Vec::new();
        for c in commitment.proof.columns.iter() {
            col_hashes.push(hash_column::<D, F>(c).into());
        }

        // 2. Compute t column indices to check the linear combination at
        let num_encoded_rows = commitment.m * rho_inv;
        let indices = get_indices_from_transcript::<F>(num_encoded_rows, t, transcript);

        // 3. Verify the paths for each of the leaf hashes
        for (leaf, i) in col_hashes.into_iter().zip(indices.iter()) {
            // TODO handle the error here
            let path = &commitment.proof.paths[*i];
            assert!(path.leaf_index == *i, "Path is for a different index!"); // TODO return an error

            path.verify(
                leaf_hash_params,
                two_to_one_params,
                &commitment.root,
                leaf.clone(),
            )
            .unwrap();
        }

        // 4. Compute the encoding w = E(v)
        let w = reed_solomon(&commitment.proof.v, rho_inv);

        // 5. Verify the random linear combinations
        for (transcript_index, matrix_index) in indices.into_iter().enumerate() {
            if inner_product(coeffs, &commitment.proof.columns[transcript_index]) != w[matrix_index]
            {
                // TODO return proper error
                return Err(Error::IncorrectInputLength(
                    "Incorrect linear combination".to_string(),
                ));
            }
        }

        Ok(())
    }
    fn compute_matrices(polynomial: &P) -> (Matrix<F>, Matrix<F>) {

        let mut coeffs = polynomial.coeffs().to_vec();

        // 1. Computing parameters and initial matrix
        // want: ceil(sqrt(f.degree() + 1)); need to deal with usize -> f64 conversion
        let num_elems = polynomial.degree() + 1;

        // TODO move this check to the constructor?
        assert_eq!(
            (num_elems as f64) as usize,
            num_elems,
            "Degree of polynomial + 1 cannot be converted to f64: aborting"
        );
        let m = (num_elems as f64).sqrt().ceil() as usize;
        // TODO: check if fft_domain.compute_size_of_domain(m) is large enough

        // padding the coefficient vector with zeroes
        // TODO is this the most efficient/safest way to do it?
        coeffs.resize(m * m, F::zero());

        let mat = Matrix::new_from_flat(m, m, &coeffs);

        // 2. Apply Reed-Solomon encoding row-wise
        let ext_mat = Matrix::new_from_rows(
            mat.rows()
                .iter()
                .map(|r| reed_solomon(r, rho_inv))
                .collect(),
        );

        (mat, ext_mat)
    }
    fn create_merkle_tree(
        ext_mat: &Matrix<F>,
        leaf_hash_params: &<<C as Config>::LeafHash as CRHScheme>::Parameters,
        two_to_one_params: &<<C as Config>::TwoToOneHash as TwoToOneCRHScheme>::Parameters,
    ) -> MerkleTree::<C> 
    {
        let mut col_hashes: Vec<Vec<u8>> = Vec::new();
        let ext_mat_cols = ext_mat.cols();

        for col in ext_mat_cols.iter() {
            col_hashes.push(hash_column::<D, F>(col).into());
        }

        MerkleTree::<C>::new(leaf_hash_params, two_to_one_params, col_hashes).unwrap()
    }
    fn generate_proof(
        coeffs: &[F],
        mat: &Matrix<F>,
        ext_mat: &Matrix<F>,
        col_tree: &MerkleTree<C>,
        transcript: &mut IOPTranscript<F>
     ) -> LigeroPCProof<F, C> {

        let m = mat.m;
        let t =  calculate_t(rho_inv, sec_param);

        // 1. Compute the linear combination using the random coefficients
        let v = mat.row_mul(coeffs);

        transcript.append_serializable_element(b"v", &v).unwrap();

        // 2. Generate t column indices to test the linear combination on
        let indices = get_indices_from_transcript(m * rho_inv, t, transcript);

        // 3. Compute Merkle tree paths for the columns
        let mut queried_columns = Vec::new();
        let mut paths = Vec::new();

        let ext_mat_cols = ext_mat.cols();

        for i in indices {
            queried_columns.push(ext_mat_cols[i].clone());
            paths.push(col_tree.generate_proof(i).unwrap());
        }

        LigeroPCProof {
            paths,
            v,
            columns: queried_columns,
        }
    }
}

impl<F, P, S, C, D, const rho_inv: usize, const sec_param: usize> PolynomialCommitment<F, P, S>
    for Ligero<F, C, D, S, P, rho_inv, sec_param>
where
    F: PrimeField,
    P: DenseUVPolynomial<F>,
    S: CryptographicSponge,
    C: Config + 'static,
    Vec<u8>: Borrow<C::Leaf>,
    <<C as Config>::TwoToOneHash as TwoToOneCRHScheme>::Parameters: Debug,
    <<C as Config>::LeafHash as CRHScheme>::Parameters: Debug,
    C::InnerDigest: Absorb,
    D: Digest,
{
    type UniversalParams = LigeroPCUniversalParams;

    type CommitterKey = LigeroPCCommitterKey<C>;

    type VerifierKey = LigeroPCVerifierKey<C>;

    type PreparedVerifierKey = LigeroPCPreparedVerifierKey;

    type Commitment = LigeroPCCommitment<F, C>;

    type PreparedCommitment = LigeroPCPreparedCommitment;

    type Randomness = LigeroPCRandomness;

    type Proof = LigeroPCProof<F, C>;

    type BatchProof = Vec<Self::Proof>;

    type Error = Error;

    fn setup<R: RngCore>(
        max_degree: usize,
        num_vars: Option<usize>,
        rng: &mut R,
    ) -> Result<Self::UniversalParams, Self::Error> {
        assert!(rho_inv >= 1, "rho_inv is the inverse of the rate and must be at least 1");
        // The domain will have size m * rho_inv, but we already have the first m elements
        GeneralEvaluationDomain::<F>::compute_size_of_domain(max_degree * (rho_inv - 1)).ok_or(Error::UnsupportedDegreeBound(max_degree))?;

        Ok(LigeroPCUniversalParams::default())
    }

    fn trim(
        pp: &Self::UniversalParams,
        supported_degree: usize,
        supported_hiding_bound: usize,
        enforced_degree_bounds: Option<&[usize]>,
    ) -> Result<(Self::CommitterKey, Self::VerifierKey), Self::Error> {
        todo!();
    }

    fn commit<'a>(
        ck: &Self::CommitterKey,
        polynomials: impl IntoIterator<Item = &'a LabeledPolynomial<F, P>>,
        rng: Option<&mut dyn RngCore>,
    ) -> Result<
        (
            Vec<LabeledCommitment<Self::Commitment>>,
            Vec<Self::Randomness>,
        ),
        Self::Error,
    >
    where
        P: 'a,
    {
        // TODO loop over all polynomials

        // TODO decide what to do with label and degree bound (these are private! but the commitment also has them)
        let labeled_polynomial = polynomials.into_iter().next().unwrap();

        let polynomial = labeled_polynomial.polynomial();

        // 1. Compute matrices
        let (mat, ext_mat) = Self::compute_matrices(polynomial);

        // 2. Create the Merkle tree from the hashes of the columns
        let col_tree = Self::create_merkle_tree(
            &ext_mat,
            &ck.leaf_hash_params,
            &ck.two_to_one_params
        );

        // 3. Add root to transcript and generate random linear combination with it
        let root = col_tree.root();

        let mut transcript: IOPTranscript<F> = IOPTranscript::new(b"well_formedness_transcript");
        transcript
            .append_serializable_element(b"root", &root)
            .unwrap();

        let m = mat.m;
        let mut r = Vec::new();
        for _ in 0..m {
            r.push(transcript.get_and_append_challenge(b"r").unwrap());
        }

        // 4. Generate the proof by choosing random columns and proving their paths in the tree
        let proof = Self::generate_proof(
            &r,
            &mat,
            &ext_mat,
            &col_tree,
            &mut transcript
        );

        let commitment = LigeroPCCommitment { m, root, proof };

        Ok((
            vec![LabeledCommitment::new(
                labeled_polynomial.label().clone(),
                commitment,
                None, // TODO think about this (degree_bound)
            )],
            Vec::new(),
        ))
        // TODO when should this return Err?
    }

    fn open<'a>(
        ck: &Self::CommitterKey,
        labeled_polynomials: impl IntoIterator<Item = &'a LabeledPolynomial<F, P>>,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Self::Commitment>>,
        point: &'a P::Point,
        challenge_generator: &mut crate::challenge::ChallengeGenerator<F, S>,
        rands: impl IntoIterator<Item = &'a Self::Randomness>,
        rng: Option<&mut dyn RngCore>,
    ) -> Result<Self::Proof, Self::Error>
    where
        P: 'a,
        Self::Randomness: 'a,
        Self::Commitment: 'a,
    {
        let labeled_polynomial = labeled_polynomials.into_iter().next().unwrap();
        let polynomial = labeled_polynomial.polynomial();

        // 1. Compute matrices
        let (mat, ext_mat) = Self::compute_matrices(polynomial);

        // 2. Create the Merkle tree from the hashes of the columns
        let col_tree = Self::create_merkle_tree(
            &ext_mat,
            &ck.leaf_hash_params,
            &ck.two_to_one_params
        );

        // 3. Generate vector b and add v = b·M to the transcript
        let m = mat.m;

        let mut b = Vec::new();
        let point_pow_m = point.pow([m as u64]);
        let mut acc_b = F::one();
        for _ in 0..m {
            b.push(acc_b);
            acc_b *= point_pow_m;
        }

        let mut transcript: IOPTranscript<F> = IOPTranscript::new(b"opening_transcript");

        let v = mat.row_mul(&b);
        transcript.append_serializable_element(b"point", point).unwrap();
        transcript.append_serializable_element(b"v", &v).unwrap();

        Ok(Self::generate_proof(
            &b,
            &mat,
            &ext_mat,
            &col_tree,
            &mut transcript
        ))
    }

    fn check<'a>(
        vk: &Self::VerifierKey,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Self::Commitment>>,
        point: &'a P::Point,
        values: impl IntoIterator<Item = F>,
        proof: &Self::Proof,
        challenge_generator: &mut crate::challenge::ChallengeGenerator<F, S>,
        rng: Option<&mut dyn RngCore>,
    ) -> Result<bool, Self::Error>
    where
        Self::Commitment: 'a,
    {
        let labeled_commitment = commitments.into_iter().next().unwrap();
        let commitment = labeled_commitment.commitment();

        let m = commitment.m;
        let t = calculate_t(rho_inv, sec_param);

        // check if we've seen this commitment before. If not, we should verify it.
        Self::check_well_formedness(commitment, &vk.leaf_hash_params, &vk.two_to_one_params)
            .unwrap();

        // 1. Compute a and b
        let mut a = Vec::new();
        let mut acc_a = F::one();
        for _ in 0..m {
            a.push(acc_a);
            acc_a *= point;
        }

        // by now acc_a = point^m
        let mut b = Vec::new();
        let mut acc_b = F::one();
        for _ in 0..m {
            b.push(acc_b);
            acc_b *= acc_a;
        }

        // 2. Seed the transcript with the point and generate t random indices
        // TODO replace unwraps by proper error handling
        let mut transcript: IOPTranscript<F> = IOPTranscript::new(b"opening_transcript");
        transcript
            .append_serializable_element(b"point", point)
            .unwrap();
        transcript.append_serializable_element(b"v", &proof.v).unwrap();

        // 3. Check the linear combination in the proof
        Self::check_random_linear_combination(
            &b,
            commitment,
            t,
            &mut transcript,
            &vk.leaf_hash_params,
            &vk.two_to_one_params,
        )?;

        Ok(inner_product(&commitment.proof.v, &a) == values.into_iter().next().unwrap())
    }
}

// TODO start considering degree bound
