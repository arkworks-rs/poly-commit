use ark_crypto_primitives::crh::CRHScheme;
use ark_crypto_primitives::crh::TwoToOneCRHScheme;
use ark_crypto_primitives::{
    merkle_tree::{Config, LeafParam, MerkleTree, Path, TwoToOneParam},
    sponge::{Absorb, CryptographicSponge},
};
use ark_ff::PrimeField;
use ark_poly::{DenseUVPolynomial, EvaluationDomain, GeneralEvaluationDomain};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use core::marker::PhantomData;
use digest::Digest;
use jf_primitives::pcs::transcript::IOPTranscript;

use crate::LabeledPolynomial;
use crate::{
    ligero::utils::{inner_product, reed_solomon},
    Error, LabeledCommitment, PCCommitment, PCCommitterKey, PCPreparedCommitment,
    PCPreparedVerifierKey, PCRandomness, PCUniversalParams, PCVerifierKey, PolynomialCommitment,
};

use ark_std::rand::RngCore;

mod utils;
use utils::Matrix;

use self::utils::get_indices_from_transcript;
use self::utils::hash_array;
mod tests;

// TODO: Disclaimer: no hiding prop
/// The Ligero polynomial commitment scheme.
pub struct Ligero<
    F: PrimeField,
    C: Config,
    D: Digest,
    S: CryptographicSponge,
    /// one over the rate rho
    const rho_inv: usize,
    /// security parameter, used in calculating t
    const sec_param: usize,
> {
    _field: PhantomData<F>,
    _config: PhantomData<C>,
    _digest: PhantomData<D>,
    _sponge: PhantomData<S>,
}

// TODO come up with reasonable defaults
const DEFAULT_RHO_INV: usize = 2;
const DEFAULT_SEC_PARAM: usize = 128;

fn calculate_t(rho_inv: usize, sec_param: usize) -> usize {
    // TODO calculate t somehow
    let t = 5;
    println!("WARNING: you are using dummy t = {t}");
    t
}

impl<F, C, D, S, const rho_inv: usize, const sec_param: usize>
    Ligero<F, C, D, S, rho_inv, sec_param>
where
    F: PrimeField,
    C: Config,
    C::Leaf: Sized + From<Vec<u8>>,
    C::InnerDigest: Absorb,
    D: Digest,
    S: CryptographicSponge,
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
        }
    }

    /// The verifier can check the well-formedness of the commitment by taking random linear combinations.
    fn check_well_formedness(
        commitment: &LigeroPCCommitment<F, C>,
        leaf_hash_param: &LeafParam<C>,
        two_to_one_param: &TwoToOneParam<C>,
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
            leaf_hash_param,
            two_to_one_param
        )
    }
    fn check_random_linear_combination(
        coeffs: &[F],
        commitment: &LigeroPCCommitment<F, C>,
        t: usize,
        transcript: &mut IOPTranscript<F>,
        leaf_hash_param: &<<C as Config>::LeafHash as CRHScheme>::Parameters,
        two_to_one_param: &<<C as Config>::TwoToOneHash as TwoToOneCRHScheme>::Parameters
    ) -> Result<(), Error>
    {
        // 1. Hash the received columns into leaf hashes
        let mut col_hashes: Vec<C::Leaf> = Vec::new();
        for c in commitment.proof.columns.iter() {
            col_hashes.push(hash_array::<D, F>(c).into());
        }
    
        // 2. Compute t column indices to check the linear combination at
        let num_encoded_rows = commitment.m * rho_inv;
        let indices = get_indices_from_transcript::<F>(num_encoded_rows, t, transcript);

        // 3. Verify the paths for each of the leaf hashes
        for (leaf, i) in col_hashes.into_iter().zip(indices.iter()) {
            // TODO handle the error here
            let path = &commitment.proof.paths[*i];
            assert!(path.leaf_index == *i, "Path is for a different index!"); // TODO return an error

            path.verify(leaf_hash_param, two_to_one_param, &commitment.root, leaf)
                .unwrap();
        }

        // 4. Compute the encoding w = E(v)
        let fft_domain = GeneralEvaluationDomain::<F>::new(commitment.m).unwrap();
        let mut domain_iter = fft_domain.elements();
        let w = reed_solomon(&commitment.proof.v, rho_inv, fft_domain, &mut domain_iter);

        // 5. Verify the random linear combinations
        for (transcript_index, matrix_index) in indices.into_iter().enumerate() {
            if inner_product(coeffs, &commitment.proof.columns[transcript_index]) != w[matrix_index] {
                // TODO return proper error
                return Err(Error::IncorrectInputLength(
                    "Incorrect linear combination".to_string(),
                ));
            }
        }

        Ok(())

    }
}

type LigeroPCUniversalParams = ();

impl PCUniversalParams for LigeroPCUniversalParams {
    fn max_degree(&self) -> usize {
        todo!()
    }
}

type LigeroPCCommitterKey = ();

impl PCCommitterKey for LigeroPCCommitterKey {
    fn max_degree(&self) -> usize {
        todo!()
    }

    fn supported_degree(&self) -> usize {
        todo!()
    }
}

type LigeroPCVerifierKey = ();

impl PCVerifierKey for LigeroPCVerifierKey {
    fn max_degree(&self) -> usize {
        todo!()
    }

    fn supported_degree(&self) -> usize {
        todo!()
    }
}

type LigeroPCPreparedVerifierKey = ();

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
    m: usize,
    // TODO is InnerDigest the right type?
    root: C::InnerDigest,
    proof: LigeroPCProof<F, C>,
}

impl<F: PrimeField, C: Config> PCCommitment for LigeroPCCommitment<F, C> {
    fn empty() -> Self {
        todo!()
    }

    fn has_degree_bound(&self) -> bool {
        todo!()
    }
}

type LigeroPCPreparedCommitment = ();

impl<Unprepared: PCCommitment> PCPreparedCommitment<Unprepared> for LigeroPCPreparedCommitment {
    fn prepare(cm: &Unprepared) -> Self {
        todo!()
    }
}

type LigeroPCRandomness = ();

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
    C: Config
{
    /// For each of the indices in q, `paths` contains the path from the root of the merkle tree to the leaf
    paths: Vec<Path<C>>,

    /// v, s.t. E(v) = w
    v: Vec<F>,

    columns: Vec<Vec<F>>,
}

impl<F, P, S, C, D, const rho_inv: usize, const sec_param: usize> PolynomialCommitment<F, P, S>
    for Ligero<F, C, D, S, rho_inv, sec_param>
where
    F: PrimeField,
    P: DenseUVPolynomial<F>,
    S: CryptographicSponge,
    C: Config + 'static,
    C::Leaf: Sized + From<Vec<u8>>,
    C::InnerDigest: Absorb,
    D: Digest,
{
    type UniversalParams = LigeroPCUniversalParams;

    type CommitterKey = LigeroPCCommitterKey;

    type VerifierKey = LigeroPCVerifierKey;

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
        todo!()
    }

    fn trim(
        pp: &Self::UniversalParams,
        supported_degree: usize,
        supported_hiding_bound: usize,
        enforced_degree_bounds: Option<&[usize]>,
    ) -> Result<(Self::CommitterKey, Self::VerifierKey), Self::Error> {
        todo!()
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
        // 0. Recovering parameters
        let t = calculate_t(rho_inv, sec_param);
        let mut optional = crate::optional_rng::OptionalRng(rng); // TODO taken from Marlin code; change in the future?
        let leaf_hash_param = C::LeafHash::setup(&mut optional).unwrap();
        let two_to_one_param = C::TwoToOneHash::setup(&mut optional).unwrap();

        // TODO loop over all polynomials

        // TODO decide what to do with label and degree bound (these are private! but the commitment also has them)
        let labeled_polynomial = polynomials.into_iter().next().unwrap();

        let polynomial = labeled_polynomial.polynomial();

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
        let fft_domain = GeneralEvaluationDomain::<F>::new(m).unwrap();
        let mut domain_iter = fft_domain.elements();

        let ext_mat = Matrix::new_from_rows(
            mat.rows()
                .iter()
                .map(|r| reed_solomon(r, rho_inv, fft_domain, &mut domain_iter))
                .collect(),
        );

        // 3. Create the Merkle tree from the hashes of the columns
        let mut col_hashes: Vec<C::Leaf> = Vec::new();
        let ext_mat_cols = ext_mat.cols();

        for col in ext_mat_cols.iter() {
            col_hashes.push(hash_array::<D, F>(col).into());
        }

        let col_tree =
            MerkleTree::<C>::new(&leaf_hash_param, &two_to_one_param, col_hashes).unwrap();

        let root = col_tree.root();

        // 4. Add root to transcript and generate random linear combination with it
        let mut transcript: IOPTranscript<F> = IOPTranscript::new(b"well_formedness_transcript");
        transcript
            .append_serializable_element(b"root", &root)
            .unwrap();

        // 5. Generate linear combination using the matrix and random coefficients
        let mut r = Vec::new();
        for _ in 0..m {
            r.push(transcript.get_and_append_challenge(b"r").unwrap());
        }

        let v = mat.row_mul(&r);

        transcript.append_serializable_element(b"v", &v).unwrap();

        // 6. Generate t column indices to test the linear combination on
        let indices = get_indices_from_transcript(m * rho_inv, t, &mut transcript);

        // 7. Compute Merkle tree paths for the columns
        let mut queried_columns = Vec::new();
        let mut paths = Vec::new();

        for i in indices {
            queried_columns.push(ext_mat_cols[i].clone());
            paths.push(col_tree.generate_proof(i).unwrap());
        }

        let proof = LigeroPCProof {
            paths,
            v,
            columns: queried_columns,
        };

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
        let num_elems = polynomial.degree() + 1;
        let m = (num_elems as f64).sqrt().ceil() as usize;
        let t = calculate_t(rho_inv, sec_param);
        let mut optional = crate::optional_rng::OptionalRng(rng);
        let leaf_hash_param = C::LeafHash::setup(&mut optional).unwrap();
        let two_to_one_param = C::TwoToOneHash::setup(&mut optional).unwrap();

        let mut coeffs = polynomial.coeffs().to_vec();
        coeffs.resize(m * m, F::zero());
        let mat = Matrix::new_from_flat(m, m, &coeffs);

        let mut transcript: IOPTranscript<F> = IOPTranscript::new(b"opening_transcript");

        // 1. Generate vector b
        let mut b = Vec::new();
        let point_pow_m = point.pow([m as u64]);
        let mut acc_b = F::one();
        for _ in 0..m {
            b.push(acc_b);
            acc_b *= point_pow_m;
        }

        let v = mat.row_mul(&b);
        transcript.append_serializable_element(b"v", &v).unwrap();

        // 2. Generate t column indices to test the linear combination on
        let indices = get_indices_from_transcript(m * rho_inv, t, &mut transcript);

        // 3. Apply Reed-Solomon encoding row-wise
        let fft_domain = GeneralEvaluationDomain::<F>::new(m).unwrap();
        let mut domain_iter = fft_domain.elements();

        let ext_mat = Matrix::new_from_rows(
            mat.rows()
                .iter()
                .map(|r| reed_solomon(r, rho_inv, fft_domain, &mut domain_iter))
                .collect(),
        );

        // 3. Create the Merkle tree from the hashes of the columns
        let mut col_hashes: Vec<C::Leaf> = Vec::new();
        let ext_mat_cols = ext_mat.cols();

        for col in ext_mat_cols.iter() {
            col_hashes.push(hash_array::<D, F>(col).into());
        }

        let col_tree =
            MerkleTree::<C>::new(&leaf_hash_param, &two_to_one_param, col_hashes).unwrap();

        let root = col_tree.root();

        // 7. Compute Merkle tree paths for the columns
        let mut queried_columns = Vec::new();
        let mut paths = Vec::new();

        for i in indices {
            queried_columns.push(ext_mat_cols[i].clone());
            paths.push(col_tree.generate_proof(i).unwrap());
        }
        
        todo!()
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
        let mut rng = rng.unwrap();
        let labeled_commitment = commitments.into_iter().next().unwrap();
        let commitment = labeled_commitment.commitment();
        
        let m = commitment.m;
        let t = calculate_t(rho_inv, sec_param);

        let leaf_hash_param = C::LeafHash::setup(&mut rng).unwrap();
        let two_to_one_param = C::TwoToOneHash::setup(&mut rng).unwrap();

        // check if we've seen this commitment before. If not, we should verify it.
        Self::check_well_formedness(
            commitment,
            &leaf_hash_param,
            &two_to_one_param,
        ).unwrap();

        // 1. Seed the transcript with the point and generate t random indices
        // TODO replace unwraps by proper error handling
        let mut transcript: IOPTranscript<F> = IOPTranscript::new(b"opening_transcript");
        transcript
            .append_serializable_element(b"point", point)
            .unwrap();

        let indices = get_indices_from_transcript(m * rho_inv, t, &mut transcript);;

        // 2. Compute a and b
        let mut a = Vec::new();
        let mut acc_a = F::one();
        for i in 0..m {
            a.push(acc_a);
            acc_a *= point;
        }

        // by now acc_a = point^m
        let mut b = Vec::new();
        let mut acc_b = F::one();
        for i in 0..m {
            b.push(acc_b);
            acc_b *= acc_a;
        }

        Self::check_random_linear_combination(
            &b,
            commitment,
            t,
            &mut transcript,
            &leaf_hash_param,
            &two_to_one_param
        )?;

        return Ok(inner_product(&commitment.proof.v, &a) == values.into_iter().next().unwrap());

        // TODO what is "values"? does it contain the actual evaluations?
        todo!()
    }
}

// TODO start considering degree bound
