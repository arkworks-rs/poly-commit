use crate::{kzg10, PCCommitterKey};
use crate::{BTreeMap, BTreeSet, ToString, Vec};
use crate::{Error, Evaluations, QuerySet};
use crate::{LabeledCommitment, LabeledPolynomial};
use crate::{PCRandomness, PCUniversalParams, Polynomial, PolynomialCommitment};

use algebra_core::{AffineCurve, One, PairingEngine, ProjectiveCurve, UniformRand, Zero};
use core::marker::PhantomData;
use rand_core::RngCore;

mod data_structures;
pub use data_structures::*;

/// Polynomial commitment based on [[KZG10]][kzg], with degree enforcement and
/// batching taken from [[MBKM19, “Sonic”]][sonic] (more precisely, their
/// counterparts in [[Gabizon19, “AuroraLight”]][al] that avoid negative G1 powers).
/// The (optional) hiding property of the commitment scheme follows the approach
/// described in [[CHMMVW20, “Marlin”]][marlin].
///
/// [kzg]: http://cacr.uwaterloo.ca/techreports/2010/cacr2010-10.pdf
/// [sonic]: https://eprint.iacr.org/2019/099
/// [al]: https://eprint.iacr.org/2019/601
/// [marlin]: https://eprint.iacr.org/2019/1047
pub struct SonicKZG10<E: PairingEngine> {
    _engine: PhantomData<E>,
}

impl<E: PairingEngine> SonicKZG10<E> {
    fn accumulate_elems<'a>(
        combined_comms: &mut BTreeMap<Option<usize>, E::G1Projective>,
        combined_witness: &mut E::G1Projective,
        combined_adjusted_witness: &mut E::G1Projective,
        vk: &VerifierKey<E>,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Commitment<E>>>,
        point: E::Fr,
        values: impl IntoIterator<Item = E::Fr>,
        proof: &kzg10::Proof<E>,
        opening_challenge: E::Fr,
        randomizer: Option<E::Fr>,
    ) {
        let acc_time = start_timer!(|| "Accumulating elements");
        let mut curr_challenge = opening_challenge;

        // Keeps track of running combination of values
        let mut combined_values = E::Fr::zero();

        // Iterates through all of the commitments and accumulates common degree_bound elements in a BTreeMap
        for (labeled_comm, value) in commitments.into_iter().zip(values) {
            combined_values += &(value * &curr_challenge);

            let comm = labeled_comm.commitment();
            let degree_bound = labeled_comm.degree_bound();

            // Applying opening challenge and randomness (used in batch_checking)
            let mut comm_with_challenge: E::G1Projective = comm.0.mul(curr_challenge);

            if let Some(randomizer) = randomizer {
                comm_with_challenge = comm_with_challenge.mul(randomizer);
            }

            // Accumulate values in the BTreeMap
            *combined_comms
                .entry(degree_bound)
                .or_insert(E::G1Projective::zero()) += &comm_with_challenge;
            curr_challenge *= &opening_challenge;
        }

        // Push expected results into list of elems. Power will be the negative of the expected power
        let mut witness: E::G1Projective = proof.w.into_projective();
        let mut adjusted_witness: E::G1Projective = vk.g.into_projective().mul(combined_values)
            + &vk.gamma_g.into_projective().mul(proof.random_v)
            - &proof.w.into_projective().mul(point);

        if let Some(randomizer) = randomizer {
            witness = witness.mul(randomizer);
            adjusted_witness = adjusted_witness.mul(randomizer);
        }

        *combined_witness += &witness;
        *combined_adjusted_witness += &adjusted_witness;
        end_timer!(acc_time);
    }

    fn check_elems(
        combined_comms: BTreeMap<Option<usize>, E::G1Projective>,
        combined_witness: E::G1Projective,
        combined_adjusted_witness: E::G1Projective,
        vk: &VerifierKey<E>,
    ) -> Result<bool, Error> {
        let check_time = start_timer!(|| "Checking elems");
        let mut g1_projective_elems: Vec<E::G1Projective> = Vec::new();
        let mut g2_prepared_elems: Vec<E::G2Prepared> = Vec::new();

        for (degree_bound, comm) in combined_comms.into_iter() {
            let shift_power = if let Some(degree_bound) = degree_bound {
                vk.get_shift_power(degree_bound)
                    .ok_or(Error::UnsupportedDegreeBound(degree_bound))?
            } else {
                vk.prepared_h.clone()
            };

            g1_projective_elems.push(comm);
            g2_prepared_elems.push(shift_power);
        }

        g1_projective_elems.push(-combined_adjusted_witness);
        g2_prepared_elems.push(vk.prepared_h.clone());

        g1_projective_elems.push(-combined_witness);
        g2_prepared_elems.push(vk.prepared_beta_h.clone());

        let g1_prepared_elems_iter =
            E::G1Projective::batch_normalization_into_affine(g1_projective_elems.as_slice())
                .into_iter()
                .map(|a| a.into());

        let g1_g2_prepared: Vec<(E::G1Prepared, E::G2Prepared)> =
            g1_prepared_elems_iter.zip(g2_prepared_elems).collect();
        let is_one: bool = E::product_of_pairings(g1_g2_prepared.iter()).is_one();
        end_timer!(check_time);
        Ok(is_one)
    }
}

impl<E: PairingEngine> PolynomialCommitment<E::Fr> for SonicKZG10<E> {
    type UniversalParams = UniversalParams<E>;
    type CommitterKey = CommitterKey<E>;
    type VerifierKey = VerifierKey<E>;
    type Commitment = Commitment<E>;
    type Randomness = Randomness<E>;
    type Proof = kzg10::Proof<E>;
    type BatchProof = Vec<Self::Proof>;
    type Error = Error;

    fn setup<R: RngCore>(
        max_degree: usize,
        rng: &mut R,
    ) -> Result<Self::UniversalParams, Self::Error> {
        kzg10::KZG10::setup(max_degree, true, rng).map_err(Into::into)
    }

    fn trim(
        pp: &Self::UniversalParams,
        supported_degree: usize,
        enforced_degree_bounds: Option<&[usize]>,
    ) -> Result<(Self::CommitterKey, Self::VerifierKey), Self::Error> {
        let trim_time = start_timer!(|| "Trimming public parameters");
        let prepared_neg_powers_of_h = pp.prepared_neg_powers_of_h.as_ref().unwrap();
        let max_degree = pp.max_degree();
        if supported_degree > max_degree {
            return Err(Error::TrimmingDegreeTooLarge);
        }

        let enforced_degree_bounds = enforced_degree_bounds.map(|bounds| {
            let mut v = bounds.to_vec();
            v.sort();
            v.dedup();
            v
        });

        let (
            shifted_powers_of_g,
            shifted_powers_of_gamma_g,
            degree_bounds_and_prepared_neg_powers_of_h,
        ) = if let Some(enforced_degree_bounds) = enforced_degree_bounds.as_ref() {
            if enforced_degree_bounds.is_empty() {
                (None, None, None)
            } else {
                let highest_enforced_degree_bound = *enforced_degree_bounds.last().unwrap();
                if highest_enforced_degree_bound > supported_degree {
                    return Err(Error::UnsupportedDegreeBound(highest_enforced_degree_bound));
                }

                let lowest_shift_degree = max_degree - highest_enforced_degree_bound;

                let shifted_ck_time = start_timer!(|| format!(
                    "Constructing `shifted_powers` of size {}",
                    max_degree - lowest_shift_degree + 1
                ));

                let shifted_powers_of_g = pp.powers_of_g[lowest_shift_degree..].to_vec();
                let shifted_powers_of_gamma_g =
                    pp.powers_of_gamma_g[lowest_shift_degree..].to_vec();

                end_timer!(shifted_ck_time);

                let prepared_neg_powers_of_h_time = start_timer!(|| format!(
                    "Constructing `prepared_neg_powers_of_h` of size {}",
                    enforced_degree_bounds.len()
                ));

                let degree_bounds_and_prepared_neg_powers_of_h = enforced_degree_bounds
                    .iter()
                    .map(|bound| {
                        (
                            *bound,
                            prepared_neg_powers_of_h[max_degree - *bound].clone(),
                        )
                    })
                    .collect();

                end_timer!(prepared_neg_powers_of_h_time);

                (
                    Some(shifted_powers_of_g),
                    Some(shifted_powers_of_gamma_g),
                    Some(degree_bounds_and_prepared_neg_powers_of_h),
                )
            }
        } else {
            (None, None, None)
        };

        let powers_of_g = pp.powers_of_g[..=supported_degree].to_vec();
        let powers_of_gamma_g = pp.powers_of_gamma_g[..=(supported_degree + 1)].to_vec();

        let ck = CommitterKey {
            powers_of_g,
            powers_of_gamma_g,
            shifted_powers_of_g,
            shifted_powers_of_gamma_g,
            enforced_degree_bounds,
            max_degree,
        };

        let g = pp.powers_of_g[0];
        let h = pp.h;
        let beta_h = pp.beta_h;
        let gamma_g = pp.powers_of_gamma_g[0];
        let prepared_h = (&pp.prepared_h).clone();
        let prepared_beta_h = (&pp.prepared_beta_h).clone();

        let vk = VerifierKey {
            g,
            gamma_g,
            h,
            beta_h,
            prepared_h,
            prepared_beta_h,
            degree_bounds_and_prepared_neg_powers_of_h,
            supported_degree,
            max_degree,
        };

        end_timer!(trim_time);
        Ok((ck, vk))
    }

    /// Outputs a commitment to `polynomial`.
    fn commit<'a>(
        ck: &Self::CommitterKey,
        polynomials: impl IntoIterator<Item = &'a LabeledPolynomial<'a, E::Fr>>,
        rng: Option<&mut dyn RngCore>,
    ) -> Result<
        (
            Vec<LabeledCommitment<Self::Commitment>>,
            Vec<Self::Randomness>,
        ),
        Self::Error,
    > {
        let commit_time = start_timer!(|| "Committing to polynomials");
        let mut labeled_comms: Vec<LabeledCommitment<Self::Commitment>> = Vec::new();
        let mut randomness: Vec<Self::Randomness> = Vec::new();
        let rng = &mut kzg10::optional_rng::OptionalRng(rng);

        for labeled_polynomial in polynomials {
            let enforced_degree_bounds: Option<&[usize]> = ck
                .enforced_degree_bounds
                .as_ref()
                .map(|bounds| bounds.as_slice());

            kzg10::KZG10::<E>::check_degrees_and_bounds(
                ck.supported_degree(),
                ck.max_degree,
                enforced_degree_bounds,
                &labeled_polynomial,
            )?;

            let polynomial = labeled_polynomial.polynomial();
            let degree_bound = labeled_polynomial.degree_bound();
            let hiding_bound = labeled_polynomial.hiding_bound();
            let label = labeled_polynomial.label();

            let commit_time = start_timer!(|| format!(
                "Polynomial {} of degree {}, degree bound {:?}, and hiding bound {:?}",
                label,
                polynomial.degree(),
                degree_bound,
                hiding_bound,
            ));

            let powers = if let Some(degree_bound) = degree_bound {
                ck.shifted_powers(degree_bound).unwrap()
            } else {
                ck.powers()
            };

            let (comm, rand) = kzg10::KZG10::commit(&powers, &polynomial, hiding_bound, Some(rng))?;

            labeled_comms.push(LabeledCommitment::new(
                label.to_string(),
                comm,
                degree_bound,
            ));
            randomness.push(rand);
            end_timer!(commit_time);
        }

        end_timer!(commit_time);
        Ok((labeled_comms, randomness))
    }

    fn open<'a>(
        ck: &Self::CommitterKey,
        labeled_polynomials: impl IntoIterator<Item = &'a LabeledPolynomial<'a, E::Fr>>,
        point: E::Fr,
        opening_challenge: E::Fr,
        rands: impl IntoIterator<Item = &'a Self::Randomness>,
    ) -> Result<Self::Proof, Self::Error>
    where
        Self::Randomness: 'a,
    {
        let mut combined_polynomial = Polynomial::zero();
        let mut combined_rand = kzg10::Randomness::empty();
        let mut curr_challenge = opening_challenge;

        for (polynomial, rand) in labeled_polynomials.into_iter().zip(rands) {
            let enforced_degree_bounds: Option<&[usize]> = ck
                .enforced_degree_bounds
                .as_ref()
                .map(|bounds| bounds.as_slice());

            kzg10::KZG10::<E>::check_degrees_and_bounds(
                ck.supported_degree(),
                ck.max_degree,
                enforced_degree_bounds,
                &polynomial,
            )?;

            combined_polynomial += (curr_challenge, polynomial.polynomial());
            combined_rand += (curr_challenge, rand);
            curr_challenge *= &opening_challenge;
        }

        let proof_time = start_timer!(|| "Creating proof for polynomials");
        let proof = kzg10::KZG10::open(&ck.powers(), &combined_polynomial, point, &combined_rand)?;
        end_timer!(proof_time);

        Ok(proof)
    }

    fn check<'a>(
        vk: &Self::VerifierKey,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Self::Commitment>>,
        point: E::Fr,
        values: impl IntoIterator<Item = E::Fr>,
        proof: &Self::Proof,
        opening_challenge: E::Fr,
    ) -> Result<bool, Self::Error>
    where
        Self::Commitment: 'a,
    {
        let check_time = start_timer!(|| "Checking evaluations");
        let mut combined_comms: BTreeMap<Option<usize>, E::G1Projective> = BTreeMap::new();
        let mut combined_witness: E::G1Projective = E::G1Projective::zero();
        let mut combined_adjusted_witness: E::G1Projective = E::G1Projective::zero();

        Self::accumulate_elems(
            &mut combined_comms,
            &mut combined_witness,
            &mut combined_adjusted_witness,
            vk,
            commitments,
            point,
            values,
            proof,
            opening_challenge,
            None,
        );

        let res = Self::check_elems(
            combined_comms,
            combined_witness,
            combined_adjusted_witness,
            vk,
        );
        end_timer!(check_time);
        res
    }

    fn batch_check<'a, R: RngCore>(
        vk: &Self::VerifierKey,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Self::Commitment>>,
        query_set: &QuerySet<E::Fr>,
        values: &Evaluations<E::Fr>,
        proof: &Self::BatchProof,
        opening_challenge: E::Fr,
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

        let mut randomizer = E::Fr::one();

        let mut combined_comms: BTreeMap<Option<usize>, E::G1Projective> = BTreeMap::new();
        let mut combined_witness: E::G1Projective = E::G1Projective::zero();
        let mut combined_adjusted_witness: E::G1Projective = E::G1Projective::zero();

        for ((query, labels), p) in query_to_labels_map.into_iter().zip(proof) {
            let mut comms_to_combine: Vec<&'_ LabeledCommitment<_>> = Vec::new();
            let mut values_to_combine = Vec::new();
            for label in labels.into_iter() {
                let commitment = commitments.get(label).ok_or(Error::MissingPolynomial {
                    label: label.to_string(),
                })?;

                let v_i = values
                    .get(&(label.clone(), *query))
                    .ok_or(Error::MissingEvaluation {
                        label: label.to_string(),
                    })?;

                comms_to_combine.push(commitment);
                values_to_combine.push(*v_i);
            }

            Self::accumulate_elems(
                &mut combined_comms,
                &mut combined_witness,
                &mut combined_adjusted_witness,
                vk,
                comms_to_combine.into_iter(),
                *query,
                values_to_combine.into_iter(),
                p,
                opening_challenge,
                Some(randomizer),
            );

            randomizer = u128::rand(rng).into();
        }

        Self::check_elems(
            combined_comms,
            combined_witness,
            combined_adjusted_witness,
            vk,
        )
    }
}

#[cfg(test)]
mod tests {
    #![allow(non_camel_case_types)]

    use crate::sonic_kzg10::SonicKZG10;
    use algebra::Bls12_377;
    use algebra::Bls12_381;
    use algebra::MNT6;
    use algebra::SW6;

    type PC<E> = SonicKZG10<E>;
    type PC_Bls12_381 = PC<Bls12_381>;
    type PC_Bls12_377 = PC<Bls12_377>;
    type PC_MNT6 = PC<MNT6>;
    type PC_SW6 = PC<SW6>;

    #[test]
    fn single_poly_test() {
        use crate::tests::*;
        single_poly_test::<_, PC_Bls12_377>().expect("test failed for bls12-377");
        single_poly_test::<_, PC_Bls12_381>().expect("test failed for bls12-381");
        single_poly_test::<_, PC_MNT6>().expect("test failed for MNT6");
        single_poly_test::<_, PC_SW6>().expect("test failed for SW6");
    }

    #[test]
    fn quadratic_poly_degree_bound_multiple_queries_test() {
        use crate::tests::*;
        quadratic_poly_degree_bound_multiple_queries_test::<_, PC_Bls12_377>()
            .expect("test failed for bls12-377");
        quadratic_poly_degree_bound_multiple_queries_test::<_, PC_Bls12_381>()
            .expect("test failed for bls12-381");
        quadratic_poly_degree_bound_multiple_queries_test::<_, PC_MNT6>()
            .expect("test failed for MNT6");
        quadratic_poly_degree_bound_multiple_queries_test::<_, PC_SW6>()
            .expect("test failed for SW6");
    }

    #[test]
    fn linear_poly_degree_bound_test() {
        use crate::tests::*;
        linear_poly_degree_bound_test::<_, PC_Bls12_377>().expect("test failed for bls12-377");
        linear_poly_degree_bound_test::<_, PC_Bls12_381>().expect("test failed for bls12-381");
        linear_poly_degree_bound_test::<_, PC_MNT6>().expect("test failed for MNT6");
        linear_poly_degree_bound_test::<_, PC_SW6>().expect("test failed for SW6");
    }

    #[test]
    fn single_poly_degree_bound_test() {
        use crate::tests::*;
        single_poly_degree_bound_test::<_, PC_Bls12_377>().expect("test failed for bls12-377");
        single_poly_degree_bound_test::<_, PC_Bls12_381>().expect("test failed for bls12-381");
        single_poly_degree_bound_test::<_, PC_MNT6>().expect("test failed for MNT6");
        single_poly_degree_bound_test::<_, PC_SW6>().expect("test failed for SW6");
    }

    #[test]
    fn single_poly_degree_bound_multiple_queries_test() {
        use crate::tests::*;
        single_poly_degree_bound_multiple_queries_test::<_, PC_Bls12_377>()
            .expect("test failed for bls12-377");
        single_poly_degree_bound_multiple_queries_test::<_, PC_Bls12_381>()
            .expect("test failed for bls12-381");
        single_poly_degree_bound_multiple_queries_test::<_, PC_MNT6>()
            .expect("test failed for MNT6");
        single_poly_degree_bound_multiple_queries_test::<_, PC_SW6>().expect("test failed for SW6");
    }

    #[test]
    fn two_polys_degree_bound_single_query_test() {
        use crate::tests::*;
        two_polys_degree_bound_single_query_test::<_, PC_Bls12_377>()
            .expect("test failed for bls12-377");
        two_polys_degree_bound_single_query_test::<_, PC_Bls12_381>()
            .expect("test failed for bls12-381");
        two_polys_degree_bound_single_query_test::<_, PC_MNT6>().expect("test failed for MNT6");
        two_polys_degree_bound_single_query_test::<_, PC_SW6>().expect("test failed for SW6");
    }

    #[test]
    fn full_end_to_end_test() {
        use crate::tests::*;
        full_end_to_end_test::<_, PC_Bls12_377>().expect("test failed for bls12-377");
        println!("Finished bls12-377");
        full_end_to_end_test::<_, PC_Bls12_381>().expect("test failed for bls12-381");
        println!("Finished bls12-381");
        full_end_to_end_test::<_, PC_MNT6>().expect("test failed for MNT6");
        println!("Finished mnt6");
        full_end_to_end_test::<_, PC_SW6>().expect("test failed for SW6");
        println!("Finished sw6");
    }

    #[test]
    fn single_equation_test() {
        use crate::tests::*;
        single_equation_test::<_, PC_Bls12_377>().expect("test failed for bls12-377");
        println!("Finished bls12-377");
        single_equation_test::<_, PC_Bls12_381>().expect("test failed for bls12-381");
        println!("Finished bls12-381");
        single_equation_test::<_, PC_MNT6>().expect("test failed for MNT6");
        println!("Finished mnt6");
        single_equation_test::<_, PC_SW6>().expect("test failed for SW6");
        println!("Finished sw6");
    }

    #[test]
    fn two_equation_test() {
        use crate::tests::*;
        two_equation_test::<_, PC_Bls12_377>().expect("test failed for bls12-377");
        println!("Finished bls12-377");
        two_equation_test::<_, PC_Bls12_381>().expect("test failed for bls12-381");
        println!("Finished bls12-381");
        two_equation_test::<_, PC_MNT6>().expect("test failed for MNT6");
        println!("Finished mnt6");
        two_equation_test::<_, PC_SW6>().expect("test failed for SW6");
        println!("Finished sw6");
    }

    #[test]
    fn two_equation_degree_bound_test() {
        use crate::tests::*;
        two_equation_degree_bound_test::<_, PC_Bls12_377>().expect("test failed for bls12-377");
        println!("Finished bls12-377");
        two_equation_degree_bound_test::<_, PC_Bls12_381>().expect("test failed for bls12-381");
        println!("Finished bls12-381");
        two_equation_degree_bound_test::<_, PC_MNT6>().expect("test failed for MNT6");
        println!("Finished mnt6");
        two_equation_degree_bound_test::<_, PC_SW6>().expect("test failed for SW6");
        println!("Finished sw6");
    }

    #[test]
    fn full_end_to_end_equation_test() {
        use crate::tests::*;
        full_end_to_end_equation_test::<_, PC_Bls12_377>().expect("test failed for bls12-377");
        println!("Finished bls12-377");
        full_end_to_end_equation_test::<_, PC_Bls12_381>().expect("test failed for bls12-381");
        println!("Finished bls12-381");
        full_end_to_end_equation_test::<_, PC_MNT6>().expect("test failed for MNT6");
        println!("Finished mnt6");
        full_end_to_end_equation_test::<_, PC_SW6>().expect("test failed for SW6");
        println!("Finished sw6");
    }
}
