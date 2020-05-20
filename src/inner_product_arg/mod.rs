use crate::{PCCommitterKey, PCCommitment, PCVerifierKey};
use crate::{BTreeMap, BTreeSet, ToString, Vec};
use crate::{BatchLCProof, Error, Evaluations, QuerySet};
use crate::{LabeledCommitment, LabeledPolynomial, LinearCombination};
use crate::{PCRandomness, PCUniversalParams, Polynomial, PolynomialCommitment};

use algebra_core::{One, PairingEngine, ProjectiveCurve, UniformRand, Zero, VariableBaseMSM, FixedBaseMSM, PrimeField, ToBytes, Field, AffineCurve, Group};
use core::{convert::TryInto, marker::PhantomData};
use rand_core::RngCore;

mod data_structures;
pub use data_structures::*;

#[cfg(feature = "parallel")]
use rayon::prelude::*;

use std::iter::Iterator;
use digest::Digest;
use ff_fft::DensePolynomial;

pub struct InnerProductArg<G: AffineCurve, D: Digest> {
    _projective: PhantomData<G>,
    _digest: PhantomData<D>
}

impl<G: AffineCurve, D: Digest> InnerProductArg<G, D> {
    pub const PROTOCOL_NAME: &'static [u8] = b"PC-DL-2020";

    fn cm_commit(
        comm_key: &[G],
        scalars: &[G::ScalarField]
    ) -> G::Projective {
        let scalars_bigint = ff_fft::cfg_iter!(scalars)
            .map(|s| s.into_repr())
            .collect::<Vec<<G::ScalarField as PrimeField>::BigInt>>();

        VariableBaseMSM::multi_scalar_mul(
            comm_key,
            scalars_bigint.as_slice(),
        )
    }

    // Takes three ToByte elements and converts them into a scalar field element
    fn rand_oracle (
        a: impl ToBytes,
        b: impl ToBytes,
        c: impl ToBytes
    ) -> G::ScalarField {
        let mut i = G::ScalarField::zero();
        let mut point = None;
        while point.is_none() {
            let bytes = algebra_core::to_bytes![a, b, c, i].unwrap();
            let hash = D::digest(bytes.as_slice());
            point = <G::ScalarField as Field>::from_random_bytes (&hash);

            i += &G::ScalarField::one();
        }

        point.unwrap()
    }

    fn inner_product(
        l: &[G::ScalarField],
        r: &[G::ScalarField]
    ) -> G::ScalarField {
        ff_fft::cfg_iter!(l)
            .zip(r)
            .map(|(li, ri)| *li * ri)
            .sum()
    }

    fn rep_poly (
        log_d: usize,
        challenges: &Vec<G::ScalarField>
    ) -> Polynomial<G::ScalarField> {
        // Representative polynomial is degree 2^log_d
        let mut coeffs = vec![G::ScalarField::one(); 1 << log_d];

        for (i, challenge) in challenges.iter().enumerate() {
            let i = i + 1;
            let elem_degree = 1 << (log_d - i);
            for start in (elem_degree..coeffs.len()).step_by(elem_degree * 2) {
                for offset in 0..elem_degree {
                    coeffs[start + offset] *= challenge;
                }
            }
        }

        Polynomial::from_coefficients_vec(coeffs)
    }

    fn succinct_check<'a> (
        vk: &VerifierKey<G>,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Commitment<G>>>,
        point: G::ScalarField,
        values: impl IntoIterator<Item = G::ScalarField>,
        proof: &Proof<G>,
        opening_challenge: G::ScalarField,
    ) -> Option <DensePolynomial<G::ScalarField>> {
        let d = vk.supported_degree();
        let log_d = ((d + 1) as f32).log2() as usize;

        let mut combined_commitment_proj = G::Projective::zero();
        let mut combined_v = G::ScalarField::zero();

        let mut curr_challenge = opening_challenge;
        let labeled_commitment = commitments.into_iter();
        let values = values.into_iter();
        for (labeled_commitment, value) in labeled_commitment.zip(values) {
            combined_v += &(curr_challenge * &value);
            combined_commitment_proj += &labeled_commitment.commitment().0.mul(curr_challenge);
            curr_challenge *= &opening_challenge;
        }

        let mut wrapper = [combined_commitment_proj];
        G::Projective::batch_normalization(wrapper.as_mut());
        let combined_commitment_proj = wrapper[0];

        // Challenge for each round
        let mut round_challenges = Vec::with_capacity(log_d);
        let mut round_challenge = Self::rand_oracle(combined_commitment_proj, point, combined_v);

        let h_prime = vk.h.mul(round_challenge);

        let mut round_commitment_proj = combined_commitment_proj + &h_prime.mul(combined_v);

        let l_iter = proof.l_vec.iter();
        let r_iter = proof.r_vec.iter();

        for (l, r) in l_iter.zip(r_iter) {
            round_challenge = Self::rand_oracle(round_challenge, l, r);
            round_challenges.push(round_challenge);
            round_commitment_proj += &(l.mul(round_challenge.inverse().unwrap()) + &r.mul(round_challenge));
        }

        let rep_poly = Self::rep_poly(log_d, &round_challenges);
        let v_prime = rep_poly.evaluate (point) * &proof.c;
        let h_prime = G::Projective::batch_normalization_into_affine(&[h_prime]).pop().unwrap();

        let check_commitment_elem: G::Projective =
            Self::cm_commit(
                &[proof.final_comm_key.clone(), h_prime],
                &[proof.c.clone(), v_prime]
            );

        if !(round_commitment_proj - &check_commitment_elem).is_zero() {
            return None;
        }

        Some (rep_poly)
    }

    fn check_degrees_and_bounds (
        supported_degree: usize,
        p: &LabeledPolynomial<G::ScalarField>,
    ) -> Result<(), Error> {
        if p.degree() < 1 {
            return Err(Error::DegreeIsZero);
        } else if p.degree() > supported_degree {
            return Err(Error::TooManyCoefficients {
                num_coefficients: p.degree() + 1,
                num_powers: supported_degree + 1
            });
        }

        if let Some(bound) = p.degree_bound() {
            if bound < p.degree() || bound > supported_degree {
                return Err(Error::IncorrectDegreeBound {
                    poly_degree: p.degree(),
                    degree_bound: bound,
                    supported_degree,
                    label: p.label().to_string(),
                });
            }
        }

        Ok(())
    }
}

impl<G: AffineCurve, D: Digest> PolynomialCommitment<G::ScalarField> for InnerProductArg<G, D> {
    type UniversalParams = UniversalParams<G>;
    type CommitterKey = CommitterKey<G>;
    type VerifierKey = VerifierKey<G>;
    type Commitment = Commitment<G>;
    type Randomness = Randomness;
    type Proof = Proof<G>;
    type BatchProof = Vec<Self::Proof>;
    type Error = Error;

    fn setup<R: RngCore>(
        max_degree: usize,
        rng: &mut R,
    ) -> Result<Self::UniversalParams, Self::Error> {
        let max_degree = (1 << ((max_degree + 1) as f32).log2().ceil() as usize) - 1;

        let mut affines = Vec::with_capacity(max_degree + 2);

        let mut i = G::ScalarField::zero();
        while affines.len() < max_degree + 2 {
            let bytes: Vec<u8> = algebra_core::to_bytes![&Self::PROTOCOL_NAME, i].unwrap();

            let hash = D::digest(bytes.as_slice());
            let affine = G::from_random_bytes (&hash);
            if affine.is_some() {
                affines.push (affine.unwrap().mul_by_cofactor());
            }

            i += &G::ScalarField::one();
        }

        let h: G = affines.pop().unwrap();

        let pp = UniversalParams {
            comm_key: affines,
            h
        };

        Ok(pp)
    }

    fn trim(
        pp: &Self::UniversalParams,
        supported_degree: usize,
        _enforced_degree_bounds: Option<&[usize]>,
    ) -> Result<(Self::CommitterKey, Self::VerifierKey), Self::Error> {
        let supported_degree = (1 << ((supported_degree + 1) as f32).log2().ceil() as usize) - 1;
        if supported_degree > pp.max_degree () {
            return Err(Error::TrimmingDegreeTooLarge);
        }

        let ck = CommitterKey {
            comm_key: pp.comm_key[0..(supported_degree + 1)].to_vec(),
            h: pp.h.clone(),
            max_degree: pp.max_degree()
        };

        let vk = VerifierKey {
            comm_key: pp.comm_key[0..(supported_degree + 1)].to_vec(),
            h: pp.h.clone(),
            max_degree: pp.max_degree()
        };

        Ok((ck, vk))
    }

    /// Outputs a commitment to `polynomial`.
    fn commit<'a>(
        ck: &Self::CommitterKey,
        polynomials: impl IntoIterator<Item = &'a LabeledPolynomial<'a, G::ScalarField>>,
        _rng: Option<&mut dyn RngCore>,
    ) -> Result<
        (
            Vec<LabeledCommitment<Self::Commitment>>,
            Vec<Self::Randomness>,
        ),
        Self::Error,
    > {
        let mut comms: Vec<LabeledCommitment<Commitment<G>>> = Vec::new();
        for labeled_polynomial in polynomials {
            Self::check_degrees_and_bounds (ck.supported_degree(), labeled_polynomial)?;

            let polynomial = labeled_polynomial.polynomial();
            let label = labeled_polynomial.label();
            let degree_bound = labeled_polynomial.degree_bound();

            let comm = Self::cm_commit(&ck.comm_key[..(polynomial.degree() + 1)], &polynomial.coeffs).into();
            let labeled_comm = LabeledCommitment::new(
                label.to_string(),
                Commitment(comm),
                degree_bound
            );

            comms.push(labeled_comm);
        }

        let num_comms = comms.len();
        Ok((comms, vec![Randomness(); num_comms]))
    }

    fn open<'a>(
        ck: &Self::CommitterKey,
        labeled_polynomials: impl IntoIterator<Item = &'a LabeledPolynomial<'a, G::ScalarField>>,
        point: G::ScalarField,
        opening_challenge: G::ScalarField,
        _rands: impl IntoIterator<Item = &'a Self::Randomness>,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Self::Commitment>>,
    ) -> Result<Self::Proof, Self::Error>
    where
        Self::Commitment: 'a,
        Self::Randomness: 'a,
    {
        let mut combined_polynomial = Polynomial::zero();
        let mut combined_commitment_proj = G::Projective::zero();

        let labeled_polynomials = labeled_polynomials.into_iter();
        let labeled_commitments = commitments.into_iter();

        let mut curr_challenge = opening_challenge;
        for (labeled_polynomial, labeled_commitment) in labeled_polynomials.zip(labeled_commitments) {
            Self::check_degrees_and_bounds (ck.supported_degree(), labeled_polynomial)?;
            combined_polynomial += (curr_challenge, labeled_polynomial.polynomial());
            combined_commitment_proj += &labeled_commitment.commitment().0.mul(curr_challenge);
            curr_challenge *= &opening_challenge;
        }

        let combined_v = combined_polynomial.evaluate(point);

        let mut wrapper = [combined_commitment_proj];
        G::Projective::batch_normalization(wrapper.as_mut());
        let combined_commitment_proj = wrapper[0];

        // ith challenge
        let mut round_challenge = Self::rand_oracle(combined_commitment_proj, point, combined_v);
        let h_prime = ck.h.mul(round_challenge);
        let h_prime = G::Projective::batch_normalization_into_affine(&[h_prime]).pop().unwrap();

        // Pad the coefficients to the appropriate vector size
        let d = ck.supported_degree();
        let log_d = ((d + 1) as f32).log2() as usize;

        let mut coeffs = combined_polynomial.coeffs;
        if coeffs.len() < d + 1 {
            for _ in coeffs.len()..(d+1) {
                coeffs.push(G::ScalarField::zero());
            }
        }
        // This extra index will be used for storing an inner product
        coeffs.push(G::ScalarField::zero());

        // Powers of z
        let mut z: Vec<G::ScalarField> = Vec::with_capacity(d + 1);
        let mut curr_z: G::ScalarField = G::ScalarField::one();
        for _ in 0..(d+1) {
            z.push(curr_z);
            curr_z *= &point;
        }

        // Key for MSM
        let mut key: Vec<G> = ck.comm_key.iter().map(|x| (*x)).collect();
        // This extra index will be used for storing h_prime
        key.push(G::zero());

        // This will be used for transforming the key in each step
        let mut key_proj =
            ck.comm_key.iter()
                .map(|x| (*x).into())
                .collect::<Vec<G::Projective>>();

        let mut l_vec = Vec::with_capacity(log_d);
        let mut r_vec = Vec::with_capacity(log_d);

        let mut n = d + 1;
        while n > 1 {
            // TODO: This implementation or previous one? This implementation merges things into the same cm commitment
            let rl_inner_prod = Self::inner_product (&coeffs[n/2..n], &z[..n/2]);
            let lr_inner_prod = Self::inner_product (&coeffs[..n/2], &z[n/2..n]);

            let tmp_key = key[n/2].clone();
            key[n/2] = h_prime.clone();
            coeffs[n] = rl_inner_prod;
            let l = Self::cm_commit (&key[..(n/2 + 1)], &coeffs[n/2..(n+1)]);
            key[n/2] = tmp_key;

            let tmp_coeff = coeffs[n/2].clone();
            key[n] = h_prime.clone();
            coeffs[n/2] = lr_inner_prod;
            let r = Self::cm_commit (&key[n/2..(n+1)], &coeffs[..(n/2+1)]);
            coeffs[n/2] = tmp_coeff;

            let lr = G::Projective::batch_normalization_into_affine(&[l, r]);
            l_vec.push(lr[0]);
            r_vec.push(lr[1]);

            round_challenge = Self::rand_oracle (round_challenge, lr[0], lr[1]);
            let round_challenge_inv = round_challenge.inverse().unwrap();

            for i in 0..n/2 {
                let coeffs_update = &(round_challenge_inv * &coeffs[n/2 + i]);
                let z_update = &(round_challenge * &z[n/2 + i]);
                let key_update = &(key_proj[n/2 + i].mul(round_challenge));

                coeffs[i] += coeffs_update;
                z[i] += z_update;
                key_proj[i] += key_update;
            }

            key = G::Projective::batch_normalization_into_affine(&key_proj[..n/2]);
            // This extra index will be used for storing h_prime
            key.push(G::zero());

            n /= 2;
        }

        Ok(Proof {
            l_vec,
            r_vec,
            final_comm_key: key[0].clone(),
            c: coeffs[0].clone()
        })
    }

    fn check<'a>(
        vk: &Self::VerifierKey,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Self::Commitment>>,
        point: G::ScalarField,
        values: impl IntoIterator<Item = G::ScalarField>,
        proof: &Self::Proof,
        opening_challenge: G::ScalarField,
    ) -> Result<bool, Self::Error>
    where
        Self::Commitment: 'a,
    {
        let d = vk.supported_degree();
        let log_d = ((d + 1) as f32).log2() as usize;

        if proof.l_vec.len() != proof.r_vec.len() || proof.l_vec.len() != log_d {
            return Err(Error::IncorrectInputLength(
                format!(
                    "Expected proof vectors to be {:}. Instead, l_vec size is {:} and r_vec size is {:}",
                    log_d,
                    proof.l_vec.len(),
                    proof.r_vec.len()
                )
            ));
        }

        let rep_poly = Self::succinct_check(vk, commitments, point, values, proof, opening_challenge);
        if rep_poly.is_none() {
            return Ok(false);
        }

        let rep_poly = rep_poly.unwrap();
        let final_key = Self::cm_commit(vk.comm_key.as_slice(), rep_poly.coeffs.as_slice());
        if !(final_key - &proof.final_comm_key.into()).is_zero() {
            return Ok(false);
        }

        Ok(true)
    }

    fn batch_check<'a, R: RngCore>(
        vk: &Self::VerifierKey,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Self::Commitment>>,
        query_set: &QuerySet<G::ScalarField>,
        values: &Evaluations<G::ScalarField>,
        proof: &Self::BatchProof,
        opening_challenge: G::ScalarField,
        rng: &mut R,
    ) -> Result<bool, Self::Error>
        where
            Self::Commitment: 'a,
    {
        let commitments: BTreeMap<_, _> = commitments.into_iter().map(|c| (c.label(), c)).collect();
        let mut query_to_labels_map = BTreeMap::new();

        for (label, point) in query_set.iter() {
            let labels = query_to_labels_map.entry(point).or_insert(BTreeSet::new());
            labels.insert(label);
        }

        assert_eq!(proof.len(), query_to_labels_map.len());

        let mut randomizer = G::ScalarField::one();

        let mut combined_rep_poly = Polynomial::zero();
        let mut combined_final_key = G::Projective::zero();

        for ((query, labels), p) in query_to_labels_map.into_iter().zip(proof) {
            let mut comms: Vec<&'_ LabeledCommitment<_>> = Vec::new();
            let mut vals = Vec::new();
            for label in labels.into_iter() {
                let commitment = commitments.get(label).ok_or(Error::MissingPolynomial {
                    label: label.to_string(),
                })?;

                let v_i = values
                    .get(&(label.clone(), *query))
                    .ok_or(Error::MissingEvaluation {
                        label: label.to_string(),
                    })?;

                comms.push(commitment);
                vals.push(*v_i);
            }

            let rep_poly= Self::succinct_check(
                vk,
                comms.into_iter(),
                *query,
                vals.into_iter(),
                p,
                opening_challenge
            );

            if (rep_poly.is_none()) {
                return Ok(false);
            }

            let rep_poly = rep_poly.unwrap();
            combined_rep_poly += (randomizer, &rep_poly);
            combined_final_key += &p.final_comm_key.into_projective().mul(randomizer);

            randomizer = u128::rand(rng).into();
        }

        let final_key = Self::cm_commit(vk.comm_key.as_slice(), combined_rep_poly.coeffs.as_slice());
        if !(final_key - &combined_final_key).is_zero() {
            return Ok(false);
        }

        Ok(true)

    }

    fn open_combinations<'a>(
        ck: &Self::CommitterKey,
        lc_s: impl IntoIterator<Item = &'a LinearCombination<G::ScalarField>>,
        polynomials: impl IntoIterator<Item = &'a LabeledPolynomial<'a, G::ScalarField>>,
        query_set: &QuerySet<G::ScalarField>,
        opening_challenge: G::ScalarField,
        _rands: impl IntoIterator<Item = &'a Self::Randomness>,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Self::Commitment>>
    ) -> Result<BatchLCProof<G::ScalarField, Self>, Self::Error>
        where
            Self::Randomness: 'a,
            Self::Commitment: 'a,
    {
        let label_poly_comm_map = polynomials
            .into_iter()
            .zip(commitments)
            .map(|(p, c)| (p.label(), (p, c)))
            .collect::<BTreeMap<_, _>>();

        let mut lc_polynomials = Vec::new();
        let mut lc_commitments = Vec::new();
        let mut lc_info = Vec::new();

        for lc in lc_s {
            let lc_label = lc.label().clone();
            let mut poly = Polynomial::zero();
            let mut degree_bound = None;
            let mut hiding_bound = None;
            let mut commitment = G::Projective::zero();

            let num_polys = lc.len();
            for (coeff, label) in lc.iter().filter(|(_, l)| !l.is_one()) {
                let label: &String = label.try_into().expect("cannot be one!");
                let &(cur_poly, cur_comm) =
                    label_poly_comm_map
                        .get(label)
                        .ok_or(Error::MissingPolynomial {
                            label: label.to_string(),
                        })?;

                if num_polys == 1 && cur_poly.degree_bound().is_some() {
                    assert!(
                        coeff.is_one(),
                        "Coefficient must be one for degree-bounded equations"
                    );
                    degree_bound = cur_poly.degree_bound();
                } else if cur_poly.degree_bound().is_some() {
                    eprintln!("Degree bound when number of equations is non-zero");
                    return Err(Self::Error::EquationHasDegreeBounds(lc_label));
                }

                // Some(_) > None, always.
                hiding_bound = core::cmp::max(hiding_bound, cur_poly.hiding_bound());
                poly += (*coeff, cur_poly.polynomial());
                commitment += &cur_comm.commitment().0.mul(*coeff);
            }

            let lc_poly = LabeledPolynomial::new_owned(lc_label.clone(), poly, degree_bound, hiding_bound);
            lc_polynomials.push(lc_poly);
            lc_commitments.push(commitment);
            lc_info.push((lc_label, degree_bound));
        }

        let comms: Vec<Self::Commitment> =
            G::Projective::batch_normalization_into_affine(&lc_commitments)
                .into_iter()
                .map(|c| Commitment(c))
                .collect();

        let lc_commitments = lc_info
            .into_iter()
            .zip(comms)
            .map(|((label, d), c)| LabeledCommitment::new(label, c, d))
            .collect::<Vec<_>>();

        let randomness = vec![Randomness(); lc_polynomials.len()];
        let proof = Self::batch_open(
            ck,
            lc_polynomials.iter(),
            &query_set,
            opening_challenge,
            randomness.iter(),
            lc_commitments.iter(),
        )?;
        Ok(BatchLCProof { proof, evals: None })
    }

    /// Checks that `values` are the true evaluations at `query_set` of the polynomials
    /// committed in `labeled_commitments`.
    fn check_combinations<'a, R: RngCore>(
        vk: &Self::VerifierKey,
        lc_s: impl IntoIterator<Item = &'a LinearCombination<G::ScalarField>>,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Self::Commitment>>,
        query_set: &QuerySet<G::ScalarField>,
        evaluations: &Evaluations<G::ScalarField>,
        proof: &BatchLCProof<G::ScalarField, Self>,
        opening_challenge: G::ScalarField,
        rng: &mut R,
    ) -> Result<bool, Self::Error>
        where
            Self::Commitment: 'a,
    {
        let BatchLCProof { proof, .. } = proof;
        let label_comm_map = commitments
            .into_iter()
            .map(|c| (c.label(), c))
            .collect::<BTreeMap<_, _>>();

        let mut lc_commitments = Vec::new();
        let mut lc_info = Vec::new();
        let mut evaluations = evaluations.clone();
        for lc in lc_s {
            let lc_label = lc.label().clone();
            let num_polys = lc.len();

            let mut degree_bound = None;
            let mut combined_comm = G::Projective::zero();

            for (coeff, label) in lc.iter() {
                if label.is_one() {
                    for (&(ref label, _), ref mut eval) in evaluations.iter_mut() {
                        if label == &lc_label {
                            **eval -= coeff;
                        }
                    }
                } else {
                    let label: &String = label.try_into().unwrap();
                    let &cur_comm = label_comm_map.get(label).ok_or(Error::MissingPolynomial {
                        label: label.to_string(),
                    })?;

                    if num_polys == 1 && cur_comm.degree_bound().is_some() {
                        assert!(
                            coeff.is_one(),
                            "Coefficient must be one for degree-bounded equations"
                        );
                        degree_bound = cur_comm.degree_bound();
                    } else if cur_comm.degree_bound().is_some() {
                        return Err(Self::Error::EquationHasDegreeBounds(lc_label));
                    }
                    combined_comm += &cur_comm.commitment().0.mul(*coeff);
                }
            }

            lc_commitments.push(combined_comm);
            lc_info.push((lc_label, degree_bound));
        }

        let comms: Vec<Self::Commitment> =
            G::Projective::batch_normalization_into_affine(&lc_commitments)
                .into_iter()
                .map(|c| Commitment(c))
                .collect();

        let lc_commitments = lc_info
            .into_iter()
            .zip(comms)
            .map(|((label, d), c)| LabeledCommitment::new(label, c, d))
            .collect::<Vec<_>>();

        Self::batch_check(
            vk,
            &lc_commitments,
            &query_set,
            &evaluations,
            proof,
            opening_challenge,
            rng,
        )
    }

}

#[cfg(test)]
mod tests {
    #![allow(non_camel_case_types)]

    use crate::inner_product_arg::InnerProductArg;

    use algebra::jubjub;
    use blake2::Blake2s;

    type PC<E, D> = InnerProductArg<E, D>;
    type PC_JJB2S = PC <jubjub::JubJubAffine, Blake2s>;

    #[test]
    fn single_poly_test() {
        use crate::tests::*;
        single_poly_test::<_, PC_JJB2S>().expect("test failed for jubjub-blake2s");
    }

    #[test]
    fn quadratic_poly_degree_bound_multiple_queries_test() {
        use crate::tests::*;
        quadratic_poly_degree_bound_multiple_queries_test::<_, PC_JJB2S>()
            .expect("test failed for jubjub-blake2s");
    }

    #[test]
    fn linear_poly_degree_bound_test() {
        use crate::tests::*;
        linear_poly_degree_bound_test::<_, PC_JJB2S>().expect("test failed for jubjub-blake2s");
    }

    #[test]
    fn single_poly_degree_bound_test() {
        use crate::tests::*;
        single_poly_degree_bound_test::<_, PC_JJB2S>().expect("test failed for jubjub-blake2s");
    }

    #[test]
    fn single_poly_degree_bound_multiple_queries_test() {
        use crate::tests::*;
        single_poly_degree_bound_multiple_queries_test::<_, PC_JJB2S>()
            .expect("test failed for jubjub-blake2s");
    }

    #[test]
    fn two_polys_degree_bound_single_query_test() {
        use crate::tests::*;
        two_polys_degree_bound_single_query_test::<_, PC_JJB2S>()
            .expect("test failed for jubjub-blake2s");
    }

    #[test]
    fn full_end_to_end_test() {
        use crate::tests::*;
        full_end_to_end_test::<_, PC_JJB2S>().expect("test failed for jubjub-blake2s");
        println!("Finished jubjub-blake2s");
    }

    #[test]
    fn single_equation_test() {
        use crate::tests::*;
        single_equation_test::<_, PC_JJB2S>().expect("test failed for jubjub-blake2s");
        println!("Finished jubjub-blake2s");
    }

    #[test]
    fn two_equation_test() {
        use crate::tests::*;
        two_equation_test::<_, PC_JJB2S>().expect("test failed for jubjub-blake2s");
        println!("Finished jubjub-blake2s");
    }

    #[test]
    fn two_equation_degree_bound_test() {
        use crate::tests::*;
        two_equation_degree_bound_test::<_, PC_JJB2S>().expect("test failed for jubjub-blake2s");
        println!("Finished jubjub-blake2s");
    }

    #[test]
    fn full_end_to_end_equation_test() {
        use crate::tests::*;
        full_end_to_end_equation_test::<_, PC_JJB2S>().expect("test failed for jubjub-blake2s");
        println!("Finished jubjub-blake2s");
    }
}
