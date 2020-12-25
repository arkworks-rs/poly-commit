use crate::{Error, LabeledCommitment, LabeledPolynomial, PolynomialCommitment};
use ark_ec::PairingEngine;
use ark_poly::polynomial::MultilinearPolynomialEvaluationForm;
use ark_poly::Polynomial;
use ark_std::marker::PhantomData;
use rand_core::RngCore;

mod data_structures;

pub struct MultilinearPC<E: PairingEngine, P: MultilinearPolynomialEvaluationForm<E::Fr>> {
    _engine: PhantomData<E>,
    _polynomial: PhantomData<P>,
}

impl<E, P> PolynomialCommitment<E::Fr, P> for MultilinearPC<E, P> {
    type UniversalParams = ();
    type CommitterKey = ();
    type VerifierKey = ();
    type PreparedVerifierKey = ();
    type Commitment = ();
    type PreparedCommitment = ();
    type Randomness = ();
    type Proof = ();
    type BatchProof = ();
    type Error = ();

    fn setup<R: RngCore>(
        _max_degree: usize,
        num_vars: Option<usize>,
        rng: &mut R,
    ) -> Result<Self::UniversalParams, Self::Error> {
        let num_vars = num_vars.expect("num_vars should not be None. ");
        todo!()
    }

    fn trim(
        pp: &Self::UniversalParams,
        supported_degree: usize,
        supported_hiding_bound: usize,
        enforced_degree_bounds: Option<&[usize]>,
    ) -> Result<(Self::CommitterKey, Self::VerifierKey), Self::Error> {
        unimplemented!("Multilinear commitment does not have degree bound")
    }

    fn commit<'a>(
        ck: &Self::CommitterKey,
        polynomials: impl IntoIterator<Item = &'a LabeledPolynomial<_, P>, IntoIter = _>,
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
        todo!()
    }

    fn open_individual_opening_challenges<'a>(
        ck: &Self::CommitterKey,
        labeled_polynomials: impl IntoIterator<Item = &'a LabeledPolynomial<_, P>, IntoIter = _>,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Self::Commitment>, IntoIter = _>,
        point: &'a <P as Polynomial<_>>::Point,
        opening_challenges: &dyn Fn(u64) -> _,
        rands: impl IntoIterator<Item = &'a Self::Randomness, IntoIter = _>,
        rng: Option<&mut dyn RngCore>,
    ) -> Result<Self::Proof, Self::Error>
    where
        P: 'a,
        Self::Randomness: 'a,
        Self::Commitment: 'a,
    {
        todo!()
    }

    fn check_individual_opening_challenges<'a>(
        vk: &Self::VerifierKey,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Self::Commitment>, IntoIter = _>,
        point: &'a <P as Polynomial<_>>::Point,
        values: impl IntoIterator<Item = _, IntoIter = _>,
        proof: &Self::Proof,
        opening_challenges: &dyn Fn(u64) -> _,
        rng: Option<&mut dyn RngCore>,
    ) -> Result<bool, Self::Error>
    where
        Self::Commitment: 'a,
    {
        todo!()
    }
}
