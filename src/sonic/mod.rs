use crate::kzg10;
use crate::{BTreeMap, BTreeSet, ToString, Vec};
use crate::{BatchLCProof, Evaluations, QuerySet, QuerySetError};
use crate::{LabeledCommitment, LabeledPolynomial, LinearCombination};
use crate::{PCRandomness, PCUniversalParams, Polynomial, PolynomialCommitment};

use algebra_core::{AffineCurve, Field, One, PairingEngine, ProjectiveCurve, Zero, UniformRand};
use core::marker::PhantomData;
use rand_core::RngCore;
use std::collections::HashMap;

mod data_structures;
pub use data_structures::*;

mod error;
pub use error::*;

pub struct Sonic<E: PairingEngine> {
    _engine: PhantomData<E>,
}

impl<E: PairingEngine> Sonic<E> {
    fn accumulate_elems<'a>(
        combined_comms: &mut HashMap<usize, E::G1Projective>,
        g1_projective_elems: &mut Vec<E::G1Projective>,
        g2_prepared_elems: &mut Vec<E::G2Prepared>,
        vk: &VerifierKey<E>,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Commitment<E>>>,
        point: E::Fr,
        values: impl IntoIterator<Item = E::Fr>,
        proof: &kzg10::Proof<E>,
        opening_challenge: E::Fr,
        randomizer: Option<E::Fr>
    ){
        let mut curr_challenge = opening_challenge;

        // Keeps track of running combination of values
        let mut combined_values = E::Fr::one();

        // Iterates through all of the commitments and accumulates common degree_bound elements in a HashMap
        for (labeled_comm, value) in commitments.into_iter().zip(values) {
            combined_values += &(value * &curr_challenge);

            // Gets the degree_bound. If there is none, then use default bound of maximum degree
            let degree_bound =
                if let Some(degree_bound) = labeled_comm.degree_bound() {
                    degree_bound
                }else{
                    vk.max_degree
                };

            let comm = labeled_comm.commitment();

            // Applying opening challenge and randomness (used in batch_checking)
            let comm_with_challenge: E::G1Projective = comm.comm.0.mul(curr_challenge);

            if let Some(randomizer) = randomizer {
                comm_with_challenge.mul(randomizer);
            }

            // Accumulate values in the HashMap
            if let Some(combined_comm) = combined_comms.get_mut(&degree_bound){
                *combined_comm += &comm_with_challenge;
            }else{
                combined_comms.insert(degree_bound, comm_with_challenge);
            }

            curr_challenge *= &opening_challenge;
        }

        // Push expected results into list of elems. Power will be the negative of the expected power
        if let Some(randomizer) = randomizer {
            g1_projective_elems.push((proof.w.mul(point) - &vk.g.mul(combined_values) - &vk.gamma_g.mul(proof.random_v)).mul(randomizer));
            g1_projective_elems.push(proof.w.mul(-randomizer));

        }else{
            g1_projective_elems.push(proof.w.mul(point) - &vk.g.mul(combined_values) - &vk.gamma_g.mul(proof.random_v));
            g1_projective_elems.push(proof.w.mul(-E::Fr::one()));
        }

        g2_prepared_elems.push((&vk.prepared_h).clone());
        g2_prepared_elems.push((&vk.prepared_beta_h).clone());
    }

    fn check_elems (
        combined_comms: HashMap<usize, E::G1Projective>,
        mut g1_projective_elems: Vec<E::G1Projective>,
        mut g2_prepared_elems: Vec<E::G2Prepared>,
        vk: &VerifierKey<E>
    ) -> Result<bool, Error>{
        for (degree_bound, comm) in combined_comms.into_iter() {
            let shift_power =
                if degree_bound == vk.max_degree {
                    vk.prepared_h.clone()
                }else{
                    vk.get_shift_power(degree_bound)
                        .ok_or(Error::UnsupportedDegreeBound(degree_bound))?
                };

            g1_projective_elems.push(comm.clone());
            g2_prepared_elems.push(shift_power);
        }

        let g1_prepared_elems_iter =
            E::G1Projective::batch_normalization_into_affine(g1_projective_elems.as_slice())
                .into_iter().map(|a| a.into());

        let g1_g2_prepared: Vec<(E::G1Prepared, E::G2Prepared)> = g1_prepared_elems_iter.zip(g2_prepared_elems).collect();

        let eq: bool =
            E::product_of_pairings(g1_g2_prepared.iter())
            == E::Fqk::one();
        Ok(eq)
    }
}

impl<E: PairingEngine> PolynomialCommitment<E::Fr> for Sonic<E> {
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

        let (shifted_powers_of_g, 
             shifted_powers_of_gamma_g, 
             degree_bounds_and_prepared_neg_powers_of_h) =

            if let Some(enforced_degree_bounds) = enforced_degree_bounds.as_ref() {
                if enforced_degree_bounds.is_empty(){
                    return Err(Error::EmptyDegreeBounds);
                }else{
                    let highest_enforced_degree_bound = *enforced_degree_bounds.last().unwrap();
                    if highest_enforced_degree_bound > supported_degree {
                        return Err(Error::UnsupportedDegreeBound(highest_enforced_degree_bound));
                    }

                    let lowest_shift_degree = max_degree - highest_enforced_degree_bound;
                    let shifted_powers_of_g = pp.powers_of_g[lowest_shift_degree..].to_vec();
                    let shifted_powers_of_gamma_g = pp.powers_of_gamma_g[lowest_shift_degree..].to_vec();
                    
                    let degree_bounds_and_prepared_neg_powers_of_h = enforced_degree_bounds
                        .iter()
                        .map(|bound| (*bound, prepared_neg_powers_of_h[max_degree - *bound].clone()))
                        .collect();

                    (Some(shifted_powers_of_g), 
                     Some(shifted_powers_of_gamma_g), 
                     Some(degree_bounds_and_prepared_neg_powers_of_h))
                }
            } else{
                (None, None, None)
            };

        let powers_of_g = pp.powers_of_g[..=supported_degree].to_vec();
        let powers_of_gamma_g = pp.powers_of_gamma_g[..=supported_degree].to_vec();

        let ck = CommitterKey {
            powers_of_g,
            powers_of_gamma_g,
            shifted_powers_of_g,
            shifted_powers_of_gamma_g,
            enforced_degree_bounds,
            max_degree,
        };

        let g = pp.powers_of_g[0];
        let gamma_g = pp.powers_of_gamma_g[0];
        let prepared_h = (&pp.prepared_h).clone();
        let prepared_beta_h = (&pp.prepared_beta_h).clone();

        let vk = VerifierKey {
            g,
            gamma_g,
            prepared_h,
            prepared_beta_h,
            degree_bounds_and_prepared_neg_powers_of_h,
            supported_degree,
            max_degree,
        };
        
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
        let mut labeled_comms: Vec<LabeledCommitment<Self::Commitment>> = Vec::new();
        let mut randomness: Vec<Self::Randomness> = Vec::new();
        let rng = &mut kzg10::optional_rng::OptionalRng(rng);

        for labeled_polynomial in polynomials {
            // TODO: Check that error function properly checks
            Self::Error::check_degrees_and_bounds(&ck, &labeled_polynomial)?;

            let polynomial = labeled_polynomial.polynomial();
            let degree_bound = labeled_polynomial.degree_bound();
            let hiding_bound = labeled_polynomial.hiding_bound();
            let label = labeled_polynomial.label();

            let powers =
                if let Some(degree_bound) = degree_bound{
                    ck.shifted_powers(degree_bound).unwrap()
                }else{
                    ck.powers()
                };

            let (comm, rand) = kzg10::KZG10::commit(&powers, &polynomial, hiding_bound, Some(rng))?;

            let commitment = Commitment {
                comm,
                has_degree_bound: degree_bound.is_some(), 
            };

            labeled_comms.push(LabeledCommitment::new(label.to_string(), commitment, degree_bound));
            randomness.push(rand);
        }

        Ok((labeled_comms, randomness))
    }

    /// On input a polynomial `p` and a point `point`, outputs a proof for the same.
    fn open<'a>(
        ck: &Self::CommitterKey,
        labeled_polynomials: impl IntoIterator<Item = &'a LabeledPolynomial<'a, E::Fr>>,
        point: E::Fr,
        opening_challenge: E::Fr,
        rands: impl IntoIterator<Item = &'a Self::Randomness>,
    ) -> Result<Self::Proof, Self::Error>
    where
        Self::Randomness: 'a
    {
        let mut combined_polynomial = Polynomial::zero();
        let mut combined_rand = kzg10::Randomness::empty();
        let mut curr_challenge = opening_challenge;

        for (polynomial, rand) in labeled_polynomials.into_iter().zip(rands){
            Self::Error::check_degrees_and_bounds(&ck, polynomial)?;
            combined_polynomial += (curr_challenge, polynomial.polynomial());
            combined_rand += (curr_challenge, rand);
            curr_challenge *= &opening_challenge;
        }

        let proof = kzg10::KZG10::open(&ck.powers(), &combined_polynomial, point, &combined_rand)?;
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
        Self::Commitment: 'a
    {
        let mut combined_comms : HashMap<usize, E::G1Projective> = HashMap::new();
        let mut g1_projective_elems: Vec<E::G1Projective> = Vec::new();
        let mut g2_prepared_elems: Vec<E::G2Prepared>= Vec::new();

        Self::accumulate_elems(
            &mut combined_comms, &mut g1_projective_elems, &mut g2_prepared_elems,
            vk, commitments, point, values, proof, opening_challenge, None
        );

        Self::check_elems(combined_comms, g1_projective_elems, g2_prepared_elems, vk)
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
        Self::Commitment: 'a
    {
        let commitments: BTreeMap<_, _> = commitments.into_iter().map(|c| (c.label(), c)).collect();
        let mut query_to_labels_map = BTreeMap::new();

        for (label, point) in query_set.iter() {
            let labels = query_to_labels_map.entry(point).or_insert(BTreeSet::new());
            labels.insert(label);
        }

        assert_eq!(proof.len(), query_to_labels_map.len());

        let mut randomizer = E::Fr::one();

        let mut combined_comms : HashMap<usize, E::G1Projective> = HashMap::new();
        let mut g1_projective_elems: Vec<E::G1Projective> = Vec::new();
        let mut g2_prepared_elems: Vec<E::G2Prepared>= Vec::new();

        for ((query, labels), p) in query_to_labels_map.into_iter().zip(proof) {
            let mut comms_to_combine: Vec<&'_ LabeledCommitment<_>> = Vec::new();
            let mut values_to_combine = Vec::new();
            for label in labels.into_iter() {
                let commitment =
                    commitments
                        .get(label)
                        .ok_or(QuerySetError::MissingPolynomial {
                            label: label.to_string(),
                        })?;

                let v_i = values.get(&(label.clone(), *query)).ok_or(
                    QuerySetError::MissingEvaluation {
                        label: label.to_string(),
                    },
                )?;

                comms_to_combine.push(commitment);
                values_to_combine.push(*v_i);
            }

            Self::accumulate_elems(
                &mut combined_comms, &mut g1_projective_elems, &mut g2_prepared_elems,
                vk, comms_to_combine.into_iter(), *query, values_to_combine.into_iter(), p, opening_challenge,
                Some(randomizer)
            );

            randomizer = u128::rand(rng).into();
        }

        Self::check_elems(combined_comms, g1_projective_elems, g2_prepared_elems, vk)
    }
}

#[cfg(test)]
mod tests {
    #![allow(non_camel_case_types)]

    use crate::sonic::Sonic;
    use algebra::curves::bls12_377::Bls12_377;
    use algebra::curves::bls12_381::Bls12_381;
    use algebra::curves::mnt6::MNT6;
    use algebra::curves::sw6::SW6;

    type PC<E> = Sonic<E>;
    type PC_Bls12_381 = PC<Bls12_381>;
    type PC_Bls12_377 = PC<Bls12_377>;
    type PC_MNT6 = PC<MNT6>;
    type PC_SW6 = PC<SW6>;

    #[test]
    fn single_poly_test() {
        use crate::tests::*;
        single_poly_test::<_, PC_Bls12_377>().expect("test failed for bls12-377");
        /*
        single_poly_test::<_, PC_Bls12_381>().expect("test failed for bls12-381");
        single_poly_test::<_, PC_MNT6>().expect("test failed for MNT6");
        single_poly_test::<_, PC_SW6>().expect("test failed for SW6");
        */
    }


    #[test]
    fn quadratic_poly_degree_bound_multiple_queries_test() {
        use crate::tests::*;
        quadratic_poly_degree_bound_multiple_queries_test::<_, PC_Bls12_377>().expect("test failed for bls12-377");
        /*
        quadratic_poly_degree_bound_multiple_queries_test::<_, PC_Bls12_381>().expect("test failed for bls12-381");
        quadratic_poly_degree_bound_multiple_queries_test::<_, PC_MNT6>().expect("test failed for MNT6");
        quadratic_poly_degree_bound_multiple_queries_test::<_, PC_SW6>().expect("test failed for SW6");
        */
    }

    #[test]
    fn linear_poly_degree_bound_test() {
        use crate::tests::*;
        linear_poly_degree_bound_test::<_, PC_Bls12_377>().expect("test failed for bls12-377");
        /*
        linear_poly_degree_bound_test::<_, PC_Bls12_381>().expect("test failed for bls12-381");
        linear_poly_degree_bound_test::<_, PC_MNT6>().expect("test failed for MNT6");
        linear_poly_degree_bound_test::<_, PC_SW6>().expect("test failed for SW6");
        */
    }


    #[test]
    fn single_poly_degree_bound_test() {
        use crate::tests::*;
        single_poly_degree_bound_test::<_, PC_Bls12_377>().expect("test failed for bls12-377");
        /*
        single_poly_degree_bound_test::<_, PC_Bls12_381>().expect("test failed for bls12-381");
        single_poly_degree_bound_test::<_, PC_MNT6>().expect("test failed for MNT6");
        single_poly_degree_bound_test::<_, PC_SW6>().expect("test failed for SW6");
        */
    }

    #[test]
    fn single_poly_degree_bound_multiple_queries_test() {
        use crate::tests::*;
        single_poly_degree_bound_multiple_queries_test::<_, PC_Bls12_377>()
            .expect("test failed for bls12-377");
        /*
        single_poly_degree_bound_multiple_queries_test::<_, PC_Bls12_381>()
            .expect("test failed for bls12-381");
        single_poly_degree_bound_multiple_queries_test::<_, PC_MNT6>()
            .expect("test failed for MNT6");
        single_poly_degree_bound_multiple_queries_test::<_, PC_SW6>()
            .expect("test failed for SW6");
            */
    }

    #[test]
    fn two_polys_degree_bound_single_query_test() {
        use crate::tests::*;
        two_polys_degree_bound_single_query_test::<_, PC_Bls12_377>()
            .expect("test failed for bls12-377");
        /*
        two_polys_degree_bound_single_query_test::<_, PC_Bls12_381>()
            .expect("test failed for bls12-381");
        two_polys_degree_bound_single_query_test::<_, PC_MNT6>()
            .expect("test failed for MNT6");
        two_polys_degree_bound_single_query_test::<_, PC_SW6>().expect("test failed for SW6");
        */
    }

    #[test]
    fn full_end_to_end_test() {
        use crate::tests::*;
        full_end_to_end_test::<_, PC_Bls12_377>().expect("test failed for bls12-377");
        println!("Finished bls12-377");
        /*
        full_end_to_end_test::<_, PC_Bls12_381>().expect("test failed for bls12-381");
        println!("Finished bls12-381");
        full_end_to_end_test::<_, PC_MNT6>().expect("test failed for MNT6");
        println!("Finished mnt6");
        full_end_to_end_test::<_, PC_SW6>().expect("test failed for SW6");
        println!("Finished sw6");
        */
    }

    #[test]
    fn single_equation_test() {
        use crate::tests::*;
        single_equation_test::<_, PC_Bls12_377>().expect("test failed for bls12-377");
        println!("Finished bls12-377");
        /*
        single_equation_test::<_, PC_Bls12_381>().expect("test failed for bls12-381");
        println!("Finished bls12-381");
        single_equation_test::<_, PC_MNT6>().expect("test failed for MNT6");
        println!("Finished mnt6");
        single_equation_test::<_, PC_SW6>().expect("test failed for SW6");
        println!("Finished sw6");
        */
    }

    #[test]
    fn two_equation_test() {
        use crate::tests::*;
        two_equation_test::<_, PC_Bls12_377>().expect("test failed for bls12-377");
        println!("Finished bls12-377");
        /*
        two_equation_test::<_, PC_Bls12_381>().expect("test failed for bls12-381");
        println!("Finished bls12-381");
        two_equation_test::<_, PC_MNT6>().expect("test failed for MNT6");
        println!("Finished mnt6");
        two_equation_test::<_, PC_SW6>().expect("test failed for SW6");
        println!("Finished sw6");
        */
    }

    #[test]
    fn two_equation_degree_bound_test() {
        use crate::tests::*;
        two_equation_degree_bound_test::<_, PC_Bls12_377>().expect("test failed for bls12-377");
        println!("Finished bls12-377");
        /*
        two_equation_degree_bound_test::<_, PC_Bls12_381>().expect("test failed for bls12-381");
        println!("Finished bls12-381");
        two_equation_degree_bound_test::<_, PC_MNT6>().expect("test failed for MNT6");
        println!("Finished mnt6");
        two_equation_degree_bound_test::<_, PC_SW6>().expect("test failed for SW6");
        println!("Finished sw6");
        */
    }

    #[test]
    fn full_end_to_end_equation_test() {
        use crate::tests::*;
        full_end_to_end_equation_test::<_, PC_Bls12_377>().expect("test failed for bls12-377");
        println!("Finished bls12-377");
        /*
        full_end_to_end_equation_test::<_, PC_Bls12_381>().expect("test failed for bls12-381");
        println!("Finished bls12-381");
        full_end_to_end_equation_test::<_, PC_MNT6>().expect("test failed for MNT6");
        println!("Finished mnt6");
        full_end_to_end_equation_test::<_, PC_SW6>().expect("test failed for SW6");
        println!("Finished sw6");
        */
    }
}
