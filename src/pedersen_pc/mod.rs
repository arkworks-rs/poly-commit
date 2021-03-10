use crate::{BTreeMap, BTreeSet, PolynomialLabel, ToString, Vec};
use crate::{
    Error, Evaluations, LabeledCommitment, LabeledPolynomial, PCCommitterKey, PCUniversalParams,
    PolynomialCommitment, QuerySet,
};
use ark_ec::AffineCurve;
use ark_ff::{UniformRand, Zero};
use ark_poly::UVPolynomial;
use ark_std::vec;
use core::marker::PhantomData;
use rand_core::RngCore;

mod data_structures;
pub use data_structures::*;

mod pedersen;
pub use pedersen::*;

/// A simple polynomial commitment scheme that based on Pedersen commitments.
/// Hiding bounds are not supported.
pub struct PedersenPC<G: AffineCurve, P: UVPolynomial<G::ScalarField>> {
    _field: PhantomData<G>,
    _polynomial: PhantomData<P>,
}

impl<G: AffineCurve, P: UVPolynomial<G::ScalarField>> PedersenPC<G, P> {
    fn check_degrees(supported_degree: usize, p: &P) -> Result<(), Error> {
        if p.degree() < 1 {
            return Err(Error::DegreeIsZero);
        } else if p.degree() > supported_degree {
            return Err(Error::TooManyCoefficients {
                num_coefficients: p.degree() + 1,
                num_powers: supported_degree + 1,
            });
        }

        Ok(())
    }

    fn check_degrees_and_bounds(
        supported_degree: usize,
        p: &LabeledPolynomial<G::ScalarField, P>,
    ) -> Result<(), Error> {
        Self::check_degrees(supported_degree, p.polynomial())?;
        if p.hiding_bound().is_some() {
            return Err(Error::HidingBoundsUnsupported);
        }

        if p.degree_bound().is_some() {
            return Err(Error::DegreeBoundsUnsupported);
        }

        Ok(())
    }
}

impl<G: AffineCurve, P: UVPolynomial<G::ScalarField>> PolynomialCommitment<G::ScalarField, P>
    for PedersenPC<G, P>
{
    type UniversalParams = UniversalParams<G>;

    type CommitterKey = CommitterKey<G>;
    type VerifierKey = CommitterKey<G>;
    type PreparedVerifierKey = CommitterKey<G>;

    type Commitment = Commitment<G>;
    type PreparedCommitment = Commitment<G>;

    type Randomness = Randomness;
    type Proof = Proof<G::ScalarField, P>;
    type BatchProof = Vec<Proof<G::ScalarField, P>>;
    type Error = Error;

    fn setup<R: RngCore>(
        max_degree: usize,
        _: Option<usize>,
        _rng: &mut R,
    ) -> Result<Self::UniversalParams, Self::Error> {
        Ok(PedersenCommitment::setup(max_degree + 1))
    }

    fn trim(
        pp: &Self::UniversalParams,
        supported_degree: usize,
        _supported_hiding_bound: usize,
        _enforced_degree_bounds: Option<&[usize]>,
    ) -> Result<(Self::CommitterKey, Self::VerifierKey), Self::Error> {
        if supported_degree > pp.max_degree() {
            return Err(Error::TrimmingDegreeTooLarge);
        }

        let ck = PedersenCommitment::trim(pp, supported_degree + 1);
        let vk = ck.clone();

        Ok((ck, vk))
    }

    fn commit<'a>(
        ck: &Self::CommitterKey,
        polynomials: impl IntoIterator<Item = &'a LabeledPolynomial<G::ScalarField, P>>,
        _rng: Option<&mut dyn RngCore>,
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
        let supported_degree = ck.supported_degree();
        let mut commitments = Vec::new();
        for labeled_polynomial in polynomials {
            Self::check_degrees_and_bounds(supported_degree, labeled_polynomial)?;

            let polynomial = labeled_polynomial.polynomial();
            let mut coeffs = polynomial.coeffs().to_vec();
            while coeffs.len() < ck.supported_degree() + 1 {
                coeffs.push(G::ScalarField::zero());
            }

            let elem = PedersenCommitment::commit(&ck, coeffs.as_slice(), None);
            let comm = Commitment { elem };

            let labeled_comm =
                LabeledCommitment::new(labeled_polynomial.label().clone(), comm, None);

            commitments.push(labeled_comm);
        }

        let randomness = vec![Randomness; commitments.len()];
        Ok((commitments, randomness))
    }

    fn open_individual_opening_challenges<'a>(
        ck: &Self::CommitterKey,
        labeled_polynomials: impl IntoIterator<Item = &'a LabeledPolynomial<G::ScalarField, P>>,
        _commitments: impl IntoIterator<Item = &'a LabeledCommitment<Self::Commitment>>,
        _point: &G::ScalarField,
        opening_challenges: &dyn Fn(u64) -> G::ScalarField,
        _rands: impl IntoIterator<Item = &'a Self::Randomness>,
        _rng: Option<&mut dyn RngCore>,
    ) -> Result<Self::Proof, Self::Error>
    where
        Self::Randomness: 'a,
        Self::Commitment: 'a,
        P: 'a,
    {
        let supported_degree = ck.supported_degree();
        let mut combined_polynomial = P::zero();

        let mut i = 0;
        for labeled_polynomial in labeled_polynomials {
            Self::check_degrees_and_bounds(supported_degree, labeled_polynomial)?;
            combined_polynomial += (opening_challenges(i), labeled_polynomial.polynomial());
            i += 1;
        }

        let combined_polynomial =
            LabeledPolynomial::new(PolynomialLabel::new(), combined_polynomial, None, None);

        Ok(Proof {
            polynomial: combined_polynomial,
        })
    }

    fn check_individual_opening_challenges<'a>(
        vk: &Self::VerifierKey,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Self::Commitment>>,
        point: &G::ScalarField,
        values: impl IntoIterator<Item = G::ScalarField>,
        proof: &Self::Proof,
        opening_challenges: &dyn Fn(u64) -> G::ScalarField,
        _rng: Option<&mut dyn RngCore>,
    ) -> Result<bool, Self::Error>
    where
        Self::Commitment: 'a,
    {
        let supported_degree = vk.supported_degree();
        let check = Self::check_degrees(supported_degree, &proof.polynomial);
        if check.is_err() {
            return Ok(false);
        }

        let mut accumulated_value = G::ScalarField::zero();
        let mut scalar_commitment_pairs = Vec::new();

        let mut i = 0;
        for (labeled_commitment, value) in commitments.into_iter().zip(values) {
            if labeled_commitment.degree_bound().is_some() {
                return Ok(false);
            }

            let cur_challenge = opening_challenges(i);
            accumulated_value += &(value * &cur_challenge);
            scalar_commitment_pairs.push((cur_challenge, labeled_commitment.commitment().clone()));
            i += 1;
        }

        let expected_value = proof.polynomial.evaluate(point);
        if accumulated_value != expected_value {
            return Ok(false);
        }

        let mut coeffs = proof.polynomial.coeffs().to_vec();
        while coeffs.len() < supported_degree + 1 {
            coeffs.push(G::ScalarField::zero());
        }

        let accumulated_commitment: Commitment<G> = scalar_commitment_pairs.into_iter().sum();
        let expected_commitment: G = PedersenCommitment::commit(&vk, coeffs.as_slice(), None);

        Ok(expected_commitment.eq(&accumulated_commitment.elem))
    }

    fn batch_check_individual_opening_challenges<'a, R: RngCore>(
        vk: &Self::VerifierKey,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Self::Commitment>>,
        query_set: &QuerySet<G::ScalarField>,
        values: &Evaluations<G::ScalarField, P::Point>,
        proof: &Self::BatchProof,
        opening_challenges: &dyn Fn(u64) -> G::ScalarField,
        rng: &mut R,
    ) -> Result<bool, Self::Error>
    where
        Self::Commitment: 'a,
    {
        let supported_degree = vk.supported_degree();
        let commitments: BTreeMap<_, _> = commitments.into_iter().map(|c| (c.label(), c)).collect();
        let mut query_to_labels_map = BTreeMap::new();

        for (label, (point_label, point)) in query_set.iter() {
            let labels = query_to_labels_map
                .entry(point_label)
                .or_insert((point, BTreeSet::new()));
            labels.1.insert(label);
        }

        assert_eq!(proof.len(), query_to_labels_map.len());

        let mut expected_value = G::ScalarField::zero();
        let mut randomizers = Vec::new();
        let mut proof_commitments = Vec::new();

        let mut accumulated_value = G::ScalarField::zero();
        let mut scalar_commitment_pairs = Vec::new();

        for ((_query, (point, labels)), p) in query_to_labels_map.into_iter().zip(proof) {
            let query_challenge: G::ScalarField = u128::rand(rng).into();
            expected_value += &(p.polynomial.evaluate(&point) * &query_challenge);

            let mut coeffs = p.polynomial.coeffs().to_vec();
            while coeffs.len() < supported_degree + 1 {
                coeffs.push(G::ScalarField::zero());
            }

            let proof_commitment = Commitment {
                elem: PedersenCommitment::commit(&vk, coeffs.as_slice(), None),
            };

            randomizers.push(query_challenge);
            proof_commitments.push(proof_commitment);

            let mut i = 0;
            for label in labels.into_iter() {
                let cur_challenge = query_challenge * &opening_challenges(i);
                let commitment = commitments.get(label).ok_or(Error::MissingPolynomial {
                    label: label.to_string(),
                })?;

                let v_i =
                    values
                        .get(&(label.clone(), *point))
                        .ok_or(Error::MissingEvaluation {
                            label: label.to_string(),
                        })?;

                accumulated_value += &(cur_challenge * v_i);
                scalar_commitment_pairs.push((cur_challenge, commitment.commitment().clone()));
                i += 1;
            }
        }

        if expected_value != accumulated_value {
            return Ok(false);
        }

        let expected_scalar_commitment_pairs: Vec<(G::ScalarField, Commitment<G>)> =
            randomizers.into_iter().zip(proof_commitments).collect();

        let accumulated_commitment: Commitment<G> = scalar_commitment_pairs.into_iter().sum();
        let expected_commitment: Commitment<G> = expected_scalar_commitment_pairs.into_iter().sum();

        Ok(accumulated_commitment == expected_commitment)
    }
}

#[cfg(test)]
mod tests {
    #![allow(non_camel_case_types)]
    use super::PedersenPC;
    use ark_ed_on_bls12_381::EdwardsAffine;
    use ark_ed_on_bls12_381::Fr;
    use ark_ff::PrimeField;
    use ark_poly::{univariate::DensePolynomial, UVPolynomial};
    use crate::tests::TestInfo;
    use crate::tests::test_template;
    use crate::Error;

    type PC_PED = PedersenPC<EdwardsAffine, DensePolynomial<Fr>>;

    fn rand_poly<F: PrimeField>(
        degree: usize,
        _: Option<usize>,
        rng: &mut rand::prelude::StdRng,
    ) -> DensePolynomial<F> {
        DensePolynomial::rand(degree, rng)
    }

    fn rand_point<F: PrimeField>(_: Option<usize>, rng: &mut rand::prelude::StdRng) -> F {
        F::rand(rng)
    }

    #[test]
    fn full_end_to_end_test() -> Result<(), Error> {
        let info = TestInfo {
            num_iters: 100,
            max_degree: None,
            supported_degree: None,
            num_vars: None,
            num_polynomials: 10,
            enforce_degree_bounds: false,
            make_hiding: false,
            max_num_queries: 5,
            num_equations: None,
            rand_poly,
            rand_point,
        };

        test_template::<Fr, DensePolynomial<Fr>, PC_PED>(info)
    }
}
