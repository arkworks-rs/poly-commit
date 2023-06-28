use core::marker::PhantomData;

use ark_crypto_primitives::{crh::TwoToOneCRHScheme, sponge::CryptographicSponge};
use ark_ff::PrimeField;
use ark_poly::Polynomial;

use crate::{
    Error, PCCommitment, PCCommitterKey, PCPreparedCommitment, PCPreparedVerifierKey, PCRandomness,
    PCUniversalParams, PCVerifierKey, PolynomialCommitment,
};

use ark_std::rand::RngCore;

// TODO: Disclaimer: no hiding prop
/// The Ligero polynomial commitment scheme.
pub struct Ligero<F: PrimeField, H: TwoToOneCRHScheme> {
    /// one over the rate rho
    rho_inv: usize,

    /// number of columns that the verifier queries
    t: usize,

    _phantom: PhantomData<(F, H)>,
}

// TODO come up with reasonable defaults
const DEFAULT_RHO_INV: usize = 2;
const DEFAULT_SEC_PARAM: usize = 128;

fn calculate_t(rho_inv: usize, sec_param: usize) -> usize {
    // TODO calculate t somehow
    let t = 5;
    t
}

impl<F: PrimeField, H: TwoToOneCRHScheme> Ligero<F, H> {
    /// Create a new instance of Ligero.
    /// If either or both parameters are None, their default values are used.
    pub fn new(rho_inv: Option<usize>, sec_param: Option<usize>) -> Self {
        let rho_inv = rho_inv.unwrap_or(DEFAULT_RHO_INV);
        let t = calculate_t(rho_inv, sec_param.unwrap_or(DEFAULT_SEC_PARAM));

        Self {
            rho_inv,
            t,
            _phantom: PhantomData,
        }
    }
}

impl<F: PrimeField, H: TwoToOneCRHScheme> Default for Ligero<F, H> {
    /// Create an instance of Ligero with the default rho (inverse: DEFAULT_RHO_INV) and security parameter (DEFAULT_SEC_PARAM).
    fn default() -> Self {
        Self::new(Some(DEFAULT_RHO_INV), Some(DEFAULT_SEC_PARAM))
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

type LigeroPCCommitment = ();

impl PCCommitment for LigeroPCCommitment {
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

impl<F: PrimeField, P: Polynomial<F>, S: CryptographicSponge, H: TwoToOneCRHScheme>
    PolynomialCommitment<F, P, S> for Ligero<F, H>
{
    type UniversalParams = LigeroPCUniversalParams;

    type CommitterKey = LigeroPCCommitterKey;

    type VerifierKey = LigeroPCVerifierKey;

    type PreparedVerifierKey = LigeroPCPreparedVerifierKey;

    type Commitment = LigeroPCCommitment;

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
