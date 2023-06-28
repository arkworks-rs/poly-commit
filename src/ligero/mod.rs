use ark_crypto_primitives::sponge::CryptographicSponge;
use ark_ff::PrimeField;
use ark_poly::DenseUVPolynomial;

use crate::{
    Error, PCCommitment, PCCommitterKey, PCPreparedCommitment, PCPreparedVerifierKey, PCRandomness,
    PCUniversalParams, PCVerifierKey, PolynomialCommitment,
};

use ark_std::rand::RngCore;

mod utils;
use utils::Matrix;
mod tests;

// TODO: Disclaimer: no hiding prop
/// The Ligero polynomial commitment scheme.
pub struct Ligero {}

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

impl<F: PrimeField, P: DenseUVPolynomial<F>, S: CryptographicSponge> PolynomialCommitment<F, P, S>
    for Ligero
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
        let f = polynomials.into_iter().next().unwrap().polynomial();

        let mut coeffs = f.coeffs().to_vec();

        // want: ceil(sqrt(f.degree() + 1)); need to deal with usize -> f64 conversion 
        let num_elems = f.degree() + 1;
        // TODO move this check to the constructor?
        assert_eq!((num_elems as f64) as usize, num_elems, "Degree of polynomial + 1 cannot be converted to f64: aborting");
        let m = (num_elems as f64).sqrt().ceil() as usize;

        // padding the coefficient vector with zeroes
        // TODO is this the most efficient way to do it?
        coeffs.resize(m * m, F::zero()); 

        let mat = Matrix::new_from_flat( m, m, &coeffs);

        // applying Reed-Solomon code row-wise
        let ext_mat = Matrix::new_from_rows(
            mat.rows().map(|r| reed_solomon(row, self.rho_inv()))
        );

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
