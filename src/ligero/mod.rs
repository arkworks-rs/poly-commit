use core::{borrow::Borrow, marker::PhantomData};
use jf_primitives::pcs::transcript::IOPTranscript;

use ark_crypto_primitives::{
    merkle_tree::{Config, LeafParam, Path, TwoToOneParam},
    sponge::{Absorb, CryptographicSponge},
};
use ark_ff::PrimeField;
use ark_poly::Polynomial;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use digest::Digest;

use crate::{
    Error, PCCommitment, PCCommitterKey, PCPreparedCommitment, PCPreparedVerifierKey, PCRandomness,
    PCUniversalParams, PCVerifierKey, PolynomialCommitment,
};

use ark_std::rand::RngCore;

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
    t
}

impl<F, C, D, S, const rho_inv: usize, const sec_param: usize>
    Ligero<F, C, D, S, rho_inv, sec_param>
where
    F: PrimeField + Borrow<<C as Config>::Leaf> + Absorb,
    C: Config,
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
    fn verify(
        commitment: &LigeroPCCommitment<F, C>,
        leaf_hash_params: &LeafParam<C>,
        two_to_one_params: &TwoToOneParam<C>,
    ) -> Result<(), Error> {
        let t = calculate_t(rho_inv, sec_param);

        // 1. Hash the received columns, to get the leaf hashes
        let col_hashes: Vec<F> = vec![F::zero(); t];
        for c in commitment.transcript.columns.iter() {
            // TODO some hashing, with the digest?
        }

        // 2. Verify the paths for each of the leaf hashes
        // TODO need a way to relate the index to the leaf hash
        for (i, leaf) in col_hashes.iter().enumerate() {
            // TODO handle the error here
            commitment.transcript.paths[i]
                .verify(
                    leaf_hash_params,
                    two_to_one_params,
                    &commitment.root,
                    leaf.clone(),
                )
                .unwrap();
        }

        // 3. verify the random linear combinations

        // 4. Verify the Fiat-Shamir transformation
        // TODO replace unwraps by proper error handling
        let mut transcript: IOPTranscript<F> = IOPTranscript::new(b"test");
        transcript
            .append_serializable_element(b"root", &commitment.root)
            .unwrap();

        let mut r = Vec::new();
        for _ in 0..t {
            r.push(transcript.get_and_append_challenge(b"r").unwrap());
        }
        // Upon sending `v` to the Verifier, add it to the sponge
        transcript
            .append_serializable_element(b"v", &commitment.transcript.v)
            .unwrap();

        // we want to squeeze enough bytes to get the indices in the range [0, rho_inv * m)
        // TODO check whether this is right
        let bytes_to_squeeze = ((commitment.m * rho_inv) >> 8) + 1;
        let mut indices = Vec::with_capacity(t);
        for _ in 0..t {
            let mut bytes: Vec<u8> = vec![0; bytes_to_squeeze];
            let _ = transcript
                .get_and_append_byte_challenge(b"i", &mut bytes)
                .unwrap();

            // get the usize from Vec<u8>:
            let ind = bytes.iter().fold(0, |acc, &x| (acc << 8) + x as usize);
            // modulo the number of columns in the encoded matrix
            indices.push(ind % (rho_inv * commitment.m));
        }

        todo!()
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

#[derive(Derivative, CanonicalSerialize, CanonicalDeserialize)]
#[derivative(Default(bound = ""), Clone(bound = ""), Debug(bound = ""))]

struct CommitmentTranscript<F: PrimeField, C: Config> {
    /// the randomness generated by the verifier
    r: Vec<F>,

    /// set of indices at which the verifier queries the commitment matrtix
    q: Vec<usize>,

    /// For each of the indices in q, `paths` contains the path from the root of the merkle tree to the leaf
    paths: Vec<Path<C>>,

    /// v, s.t. E(v) = w
    v: Vec<F>,

    columns: Vec<Vec<F>>,
}

/// The commitment to a polynomial is a root of the merkle tree,
/// where each node is a hash of the column of the encoded coefficient matrix U.
#[derive(Derivative, CanonicalSerialize, CanonicalDeserialize)]
#[derivative(Default(bound = ""), Clone(bound = ""), Debug(bound = ""))]
pub struct LigeroPCCommitment<F: PrimeField, C: Config> {
    // TODO is InnerDigest the right type?
    root: C::InnerDigest,
    m: usize,
    transcript: CommitmentTranscript<F, C>,
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

type LigeroPCProof = ();

impl<F, P, S, C, D, const rho_inv: usize, const sec_param: usize> PolynomialCommitment<F, P, S>
    for Ligero<F, C, D, S, rho_inv, sec_param>
where
    F: PrimeField,
    P: Polynomial<F>,
    S: CryptographicSponge,
    C: Config + 'static,
    D: Digest,
{
    type UniversalParams = LigeroPCUniversalParams;

    type CommitterKey = LigeroPCCommitterKey;

    type VerifierKey = LigeroPCVerifierKey;

    type PreparedVerifierKey = LigeroPCPreparedVerifierKey;

    type Commitment = LigeroPCCommitment<F, C>;

    type PreparedCommitment = LigeroPCPreparedCommitment;

    type Randomness = LigeroPCRandomness;

    type Proof = LigeroPCProof;

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
        polynomials: impl IntoIterator<Item = &'a crate::LabeledPolynomial<F, P>>,
        rng: Option<&mut dyn RngCore>,
    ) -> Result<
        (
            Vec<crate::LabeledCommitment<Self::Commitment>>,
            Vec<Self::Randomness>,
        ),
        Self::Error,
    >
    where
        P: 'a,
    {
        todo!()
    }

    fn open<'a>(
        ck: &Self::CommitterKey,
        labeled_polynomials: impl IntoIterator<Item = &'a crate::LabeledPolynomial<F, P>>,
        commitments: impl IntoIterator<Item = &'a crate::LabeledCommitment<Self::Commitment>>,
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
        todo!()
    }

    fn check<'a>(
        vk: &Self::VerifierKey,
        commitments: impl IntoIterator<Item = &'a crate::LabeledCommitment<Self::Commitment>>,
        point: &'a P::Point,
        values: impl IntoIterator<Item = F>,
        proof: &Self::Proof,
        challenge_generator: &mut crate::challenge::ChallengeGenerator<F, S>,
        rng: Option<&mut dyn RngCore>,
    ) -> Result<bool, Self::Error>
    where
        Self::Commitment: 'a,
    {
        todo!()
    }
}
