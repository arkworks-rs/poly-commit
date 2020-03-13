use crate::{PCUniversalParams, PCRandomness, Polynomial, PolynomialCommitment};
use crate::{QuerySetError, EquationError, QuerySet, Evaluations};
use crate::{LabeledPolynomial, LabeledCommitment, Equation};
use crate::kzg10;

use algebra::{AffineCurve, Field, PairingEngine, ProjectiveCurve, One, Zero};
use algebra::msm::{FixedBaseMSM, VariableBaseMSM};
use rand_core::RngCore;
use std::marker::PhantomData;
use std::collections::{BTreeMap, BTreeSet};

mod data_structures;
pub use data_structures::*;

mod error;
pub use error::*;

pub struct Sonic<E: PairingEngine> {
    _engine: PhantomData<E>,
}

impl<E: PairingEngine> Sonic<E> {

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
        if max_degree < 1 {
            return Err(Error::DegreeIsZero);
        }

        let g = E::G1Projective::rand(rng);
        let gamma_g = E::G1Projective::rand(rng);
        let h = E::G2Projective::rand(rng);

        let beta = E::Fr::rand(rng);
        let mut powers_of_beta = vec![E::Fr::one()];
        let mut neg_powers_of_beta = vec![E::Fr::one()];

        let mut curr_pos = E::Fr::one();
        let mut curr_neg = E::Fr::one();

        for _ in 0..=max_degree {
            curr_pos *= &beta;
            curr_neg /= &beta;
            powers_of_beta.push(curr_pos);
            neg_powers_of_beta.push(curr_neg);
        }

        let window_size = FixedBaseMSM::get_mul_window_size(max_degree + 1);
        let neg_window_size = FixedBaseMSM::get_mul_window_size(max_degree + 1);
        let scalar_bits = E::Fr::size_in_bits();

        let g_table = FixedBaseMSM::get_window_table(scalar_bits, window_size, g);
        let powers_of_g = FixedBaseMSM::multi_scalar_mul::<E::G1Projective>(
            scalar_bits,
            window_size,
            &g_table,
            &powers_of_beta,
        ).into_iter().map(|e| e.into_affine()).collect();
        
        let gamma_g_table = FixedBaseMSM::get_window_table(scalar_bits, window_size, gamma_g);
        let mut powers_of_gamma_g = FixedBaseMSM::multi_scalar_mul::<E::G1Projective>(
            scalar_bits,
            window_size,
            &gamma_g_table,
            &powers_of_beta,
        ).into_iter().map(|e| e.into_affine()).collect();

        let neg_h_table = FixedBaseMSM::get_window_table(scalar_bits, neg_window_size, h);
        let mut neg_powers_of_h = FixedBaseMSM::multi_scalar_mul::<E::G2Projective>(
            scalar_bits,
            neg_window_size,
            &neg_h_table,
            &neg_powers_of_beta,
        ).into_iter().map(|e| e.into_affine()).collect();

        let beta_h = h.mul(&beta).into_affine();
        let pp = UniversalParams {
            powers_of_g,
            powers_of_gamma_g,
            neg_powers_of_h,
            beta_h
        };
        
        Ok(pp);
    }

    fn trim(
        pp: &Self::UniversalParams,
        supported_degree: usize,
        enforced_degree_bounds: Option<&[usize]>,
    ) -> Result<(Self::CommitterKey, Self::VerifierKey), Self::Error> {
        let max_degree = pp.max_degree();
        if supported_degree > max_degree {
            return Err(Error::TrimmingDegreeTooLarge);
        }

        let enforced_degree_bounds = enforced_degree_bounds.map(let|bounds| {
            let mut v = bounds.to_vec();
            v.sort();
            v.dedup();
            v
        });

        let (shifted_powers_of_g, 
             shifted_powers_of_gamma_g, 
             degree_bounds_and_neg_powers_of_h) = 

            if let Some(enforced_degree_bounds) = enforced_degree_bounds.as_ref() {
                if enforced_degree_bounds.is_empty(){
                    return Err(Error::EmptyDegreeBounds);
                }else{
                    let highest_enforced_degree_bound = enforced_degree_bounds.last().unwrap();
                    if highest_enforced_degree_bound > supported_degree {
                        return Err(Error::UnsupportedDegreeBound(highest_enforced_degree_bound));
                    }

                    let lowest_shift_degree = max_degree - highest_enforced_degree_bound;
                    let shifted_powers_of_g = pp.powers_of_g[lowest_shift_degree..].to_vec();
                    let shifted_powers_of_gamma_g = pp.powers_of_gamma_g[lowest_shift_degree..].to_vec();
                    
                    let degree_bounds_and_neg_powers_of_h = enforced_degree_bounds 
                        .iter()
                        .map(|bound| (*bound, pp.neg_powers_of_h[max_degree - bound]))
                        .collect();

                    (Some(shifted_powers_of_g), 
                     Some(shifted_powers_of_gamma_g), 
                     Some(degree_bounds_and_neg_powers_of_h))
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
        let h = pp.neg_powers_of_h[0];

        let vk = VerifierKey {
            g,
            gamma_g,
            h,
            degree_bounds_and_neg_powers_of_h,
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
        let commitments : Vec<LabeledCommitment<Self::Commitment>> = Vec::new();
        let randomness : Vec<Self::Randomness> = Vec::new();
        let rng = &mut kzg10::optional_rng::OptionalRng(rng);

        for p in polynomials {
            Self::Error::check_degrees_and_bounds(&ck, &p)?;

            let polynomial = p.polynomial();
            let degree_bound = p.degree_bound();
            let hiding_bound = p.hiding_bound();
            let label = p.label();

            let powers =
                if let Some(degree_bound) = degree_bound.as_ref(){
                    &ck.shifted_powers(degree_bound)
                }else{
                    &ck.powers()
                }

            let (comm, rand) = kzg10::KZG10::commit(powers, &polynomial, hiding_bound, Some(rng))?;

            commitment = Commitment {
                comm,
                has_degree_bound: degree_bound.is_some(), 
            };

            random = Randomness {
                rand,
            };

            commitments.push(LabeledCommitment::new(label.to_string(), commitment, degree_bound));
            randomness.push(random);
        }

        Ok((commitments, randomness));
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
        let combined_p = Polynomial::zero();
        let combined_r = kzg10::Randomness::empty();

        let curr_challenge = opening_challenge;
        for (p, rand) in labeled_polynomials.into_iter().zip(rands){
            let degree_bound = p.degree_bound();
            Self::Error::check_degrees_and_bounds(&ck, &p)?;
            combined_p += (curr_challenge, p.polynomial());
            combined_r += (curr_challenge, &r.rand);
            curr_challenge *= &opening_challenge;
        }

        let powers =
            if let Some(degree_bound) = degree_bound.as_ref(){
                &ck.shifted_powers(degree_bound)
            }else{
                &ck.powers()
            }

        let proof = kzg10::KZG10::open(powers, &combined_p, point, &combined_r)?;
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

    }

    fn open_equations<'a>(
        ck: &Self::CommitterKey,
        equations: impl IntoIterator<Item = &'a Equation<E::Fr>>,
        labeled_polynomials: impl IntoIterator<Item = &'a LabeledPolynomial<'a, E::Fr>>,
        opening_challenge: E::Fr,
        rands: impl IntoIterator<Item = &'a Self::Randomness>,
    ) -> Result<Self::BatchProof, Self::Error>
    where
        Self::Randomness: 'a
    {

    }

    fn check_equations<'a, R: RngCore>(
        vk: &Self::VerifierKey,
        equations: impl IntoIterator<Item = &'a Equation<E::Fr>>,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Self::Commitment>>,
        _: Option<&Evaluations<E::Fr>>,
        proof: &Self::BatchProof,
        opening_challenge: E::Fr,
        rng: &mut R,
    ) -> Result<bool, Self::Error>
    where
        Self::Commitment: 'a
    {

    }
}

#[cfg(test)]
mod tests {
    #![allow(non_camel_case_types)]

    use crate::marlin_kzg10::MarlinKZG10;
    use algebra::curves::bls12_377::Bls12_377;
    use algebra::curves::bls12_381::Bls12_381;
    use algebra::curves::mnt6::MNT6;
    use algebra::curves::sw6::SW6;

    type PC<E> = MarlinKZG10<E>;
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
