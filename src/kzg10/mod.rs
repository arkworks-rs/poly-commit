//! Here we constuct a polynomial commitment that enables users to commit to a
//! single polynomial `p`, and then later provide an evaluation proof that
//! convinces verifiers that a claimed value `v` is the true evaluation of `p`
//! at a chosen point `x`. Our construction follows the template of the construction
//! proposed by Kate, Zaverucha, and Goldberg ([KZG10](http://cacr.uwaterloo.ca/techreports/2010/cacr2010-10.pdf)).
//! This construction achieves extractability in the algebraic group model (AGM).

use crate::{PCRandomness, Polynomial};
use algebra::msm::{FixedBaseMSM, VariableBaseMSM};
use algebra::{
    AffineCurve, Field, Group, PairingCurve, PairingEngine, PrimeField, ProjectiveCurve,
    UniformRand,
};
use rand_core::RngCore;
use rayon::prelude::*;
use std::marker::PhantomData;

mod data_structures;
pub use data_structures::*;

mod error;
pub use error::*;

pub(crate) mod optional_rng;

/// `KZG10` is an implementation of the polynomial commitment scheme of
/// [Kate, Zaverucha and Goldbgerg][kzg10]
///
/// [kzg10]: http://cacr.uwaterloo.ca/techreports/2010/cacr2010-10.pdf
pub struct KZG10<E: PairingEngine> {
    _engine: PhantomData<E>,
}

impl<E: PairingEngine> KZG10<E> {
    /// Constructs public parameters when given as input the maximum degree `degree`
    /// for the polynomial commitment scheme.
    pub fn setup<R: RngCore>(
        mut max_degree: usize,
        _produce_g2_powers: bool,
        rng: &mut R,
    ) -> Result<UniversalParams<E>, Error> {
        if max_degree < 1 {
            return Err(Error::DegreeIsZero);
        } else if max_degree == 1 {
            // FIXME: hack to support hiding for degree one polynomials.
            max_degree += 1;
        }
        let setup_time = start_timer!(|| format!("KZG10::Setup with degree {}", degree));
        let beta = E::Fr::rand(rng);
        let g = E::G1Projective::rand(rng);
        let gamma_g = E::G1Projective::rand(rng);
        let h = E::G2Projective::rand(rng);

        let mut powers_of_beta = vec![E::Fr::one()];
        let mut cur = beta;
        for _ in 0..max_degree {
            powers_of_beta.push(cur);
            cur *= &beta;
        }

        let window_size = FixedBaseMSM::get_mul_window_size(max_degree + 1);

        let scalar_bits = E::Fr::size_in_bits();
        let g_time = start_timer!(|| "Generating powers of G");
        let g_table = FixedBaseMSM::get_window_table(scalar_bits, window_size, g);
        let mut powers_of_g = FixedBaseMSM::multi_scalar_mul::<E::G1Projective>(
            scalar_bits,
            window_size,
            &g_table,
            &powers_of_beta,
        );
        end_timer!(g_time);
        let gamma_g_time = start_timer!(|| "Generating powers of gamma * G");
        let gamma_g_table = FixedBaseMSM::get_window_table(scalar_bits, window_size, gamma_g);
        let mut powers_of_gamma_g = FixedBaseMSM::multi_scalar_mul::<E::G1Projective>(
            scalar_bits,
            window_size,
            &gamma_g_table,
            &powers_of_beta,
        );
        end_timer!(gamma_g_time);
        E::G1Projective::batch_normalization(powers_of_g.as_mut_slice());
        E::G1Projective::batch_normalization(powers_of_gamma_g.as_mut_slice());

        let beta_h = h.mul(&beta).into_affine();
        let h = h.into_affine();
        let prepared_h = h.prepare();
        let prepared_beta_h = beta_h.prepare();

        let pp = UniversalParams {
            powers_of_g: powers_of_g.into_iter().map(|e| e.into_affine()).collect(),
            powers_of_gamma_g: powers_of_gamma_g
                .into_iter()
                .map(|e| e.into_affine())
                .collect(),
            h,
            beta_h,
            prepared_h,
            prepared_beta_h,
        };
        end_timer!(setup_time);
        Ok(pp)
    }

    /// Outputs a commitment to `polynomial`.
    pub fn commit(
        powers: &Powers<E>,
        polynomial: &Polynomial<E::Fr>,
        hiding_bound: Option<usize>,
        rng: Option<&mut dyn RngCore>,
    ) -> Result<(Commitment<E>, Randomness<E>), Error> {
        Error::check_degree_is_within_bounds(polynomial.degree(), powers.size())?;

        let commit_time = start_timer!(|| format!(
            "Committing to polynomial of degree {} with hiding_bound: {:?}",
            polynomial.degree(),
            hiding_bound,
        ));

        let (num_leading_zeros, plain_coeffs) = skip_leading_zeros_and_convert_to_bigints(&polynomial);

        let msm_time = start_timer!(|| "MSM to compute commitment to plaintext poly");
        let mut commitment = VariableBaseMSM::multi_scalar_mul(&powers.powers_of_g[num_leading_zeros..], &plain_coeffs);
        end_timer!(msm_time);

        let mut randomness = Randomness::empty();
        if let Some(hiding_degree) = hiding_bound {
            let mut rng = rng.ok_or(Error::MissingRng)?;
            let sample_random_poly_time = start_timer!(|| format!(
                "Sampling a random polynomial of degree {}",
                hiding_degree
            ));

            // TODO: decide on degree of polynomial, and then adjust check in 
            // Error::check_hiding_bound.
            randomness = Randomness::rand(hiding_degree, &mut rng);
            Error::check_hiding_bound(randomness.blinding_polynomial.degree(), powers.powers_of_gamma_g.len())?;
            end_timer!(sample_random_poly_time);
        }

        let random_ints = convert_to_bigints(&randomness.blinding_polynomial.coeffs);
        let msm_time = start_timer!(|| "MSM to compute commitment to random poly");
        let random_commitment =
            VariableBaseMSM::multi_scalar_mul(&powers.powers_of_gamma_g, random_ints.as_slice())
                .into_affine();
        end_timer!(msm_time);

        commitment.add_assign_mixed(&random_commitment);

        end_timer!(commit_time);
        Ok((Commitment(commitment.into()), randomness))
    }

    /// Compute witness polynomial.
    pub fn compute_witness_polynomial(
        p: &Polynomial<E::Fr>,
        point: E::Fr,
        randomness: &Randomness<E>,
    ) -> Result<(Polynomial<E::Fr>, Option<Polynomial<E::Fr>>), Error> {
        let eval_time = start_timer!(|| "Evaluating polynomial");
        let value = p.evaluate(point);
        end_timer!(eval_time);

        let divisor = Polynomial::from_coefficients_vec(vec![-point, E::Fr::one()]);

        let witness_time = start_timer!(|| "Computing witness polynomial");
        let witness_polynomial = &(p - &Polynomial::from_coefficients_vec(vec![value])) / &divisor;
        end_timer!(witness_time);

        let random_witness_polynomial = if randomness.is_hiding() {
            let random_p = &randomness.blinding_polynomial;

            let rand_eval_time = start_timer!(|| "Evaluating random polynomial");
            let random_value = random_p.evaluate(point);
            end_timer!(rand_eval_time);

            let witness_time = start_timer!(|| "Computing random witness polynomial");
            let random_witness_polynomial = &(random_p - &Polynomial::from_coefficients_vec(vec![random_value])) / &divisor;
            end_timer!(witness_time);
            Some(random_witness_polynomial)
        } else {
            None
        };

        Ok((witness_polynomial, random_witness_polynomial))
    }

    pub(crate) fn open_with_witness_polynomial<'a>(
        powers: &Powers<E>,
        point: E::Fr,
        randomness: &Randomness<E>,
        witness_polynomial: &Polynomial<E::Fr>,
        hiding_witness_polynomial: Option<&Polynomial<E::Fr>>,
    ) -> Result<Proof<E>, Error> {
        Error::check_degree_is_too_large(witness_polynomial.degree(), powers.size())?;
        let (num_leading_zeros, witness_coeffs) = skip_leading_zeros_and_convert_to_bigints(&witness_polynomial);

        let witness_comm_time = start_timer!(|| "Computing commitment to witness polynomial");
        let mut w = VariableBaseMSM::multi_scalar_mul(&powers.powers_of_g[num_leading_zeros..], &witness_coeffs);
        end_timer!(witness_comm_time);

        let mut random_v = E::Fr::zero();
        if let Some(hiding_witness_polynomial) = hiding_witness_polynomial {
            let blinding_p = &randomness.blinding_polynomial;
            let blinding_eval_time = start_timer!(|| "Evaluating random polynomial");
            let blinding_evaluation = blinding_p.evaluate(point);
            end_timer!(blinding_eval_time);

            let random_witness_coeffs = convert_to_bigints(&hiding_witness_polynomial.coeffs);
            let witness_comm_time =
                start_timer!(|| "Computing commitment to random witness polynomial");
            w += &VariableBaseMSM::multi_scalar_mul(&powers.powers_of_gamma_g, &random_witness_coeffs);
            end_timer!(witness_comm_time);
            random_v = blinding_evaluation;
        }

        Ok(Proof {
            w: w.into_affine(),
            random_v,
        })
    }

    /// On input a polynomial `p` and a point `point`, outputs a proof for the same.
    pub(crate) fn open<'a>(
        powers: &Powers<E>,
        p: &Polynomial<E::Fr>,
        point: E::Fr,
        rand: &Randomness<E>,
    ) -> Result<Proof<E>, Error> {
        Error::check_degree_is_within_bounds(p.degree(), powers.size())?;
        let open_time = start_timer!(|| format!("Opening polynomial of degree {}", p.degree()));

        let witness_time = start_timer!(|| "Computing witness polynomials");
        let (witness_poly, hiding_witness_poly) = Self::compute_witness_polynomial(p, point, rand)?;
        end_timer!(witness_time);

        let proof = Self::open_with_witness_polynomial(powers, point, rand, &witness_poly, hiding_witness_poly.as_ref());

        end_timer!(open_time);
        proof
    }

    /// Verifies that `value` is the evaluation at `point` of the polynomial
    /// committed inside `comm`.
    pub fn check(
        vk: &VerifierKey<E>,
        comm: &Commitment<E>,
        point: E::Fr,
        value: E::Fr,
        proof: &Proof<E>,
    ) -> Result<bool, Error> {
        let check_time = start_timer!(|| "Checking evaluation");
        let inner = comm.0.into_projective()
            - &vk.g.into_projective().mul(&value)
            - &vk.gamma_g.into_projective().mul(&proof.random_v);
        let lhs = E::pairing(inner, vk.h);

        let inner = vk.beta_h.into_projective() - &vk.h.into_projective().mul(&point);
        let rhs = E::pairing(proof.w, inner);

        end_timer!(check_time, || format!("Result: {}", lhs == rhs));
        Ok(lhs == rhs)
    }

    /// Check that each `proof_i` in `proofs` is a valid proof of evaluation for 
    /// `commitment_i` at `point_i`.
    pub fn batch_check<R: RngCore>(
        vk: &VerifierKey<E>,
        commitments: &[Commitment<E>],
        points: &[E::Fr],
        values: &[E::Fr],
        proofs: &[Proof<E>],
        rng: &mut R,
    ) -> Result<bool, Error> {
        let check_time =
            start_timer!(|| format!("Checking {} evaluation proofs", commitments.len()));
        let g = vk.g.into_projective();
        let gamma_g = vk.gamma_g.into_projective();

        let mut total_c = <E::G1Projective as ProjectiveCurve>::zero();
        let mut total_w = <E::G1Projective as ProjectiveCurve>::zero();

        let combination_time = start_timer!(|| "Combining commitments and proofs");
        let mut randomizer = E::Fr::one();
        // Instead of multiplying g and gamma_g in each turn, we simply accumulate
        // their coefficients and perform a final multiplication at the end.
        let mut g_multiplier = E::Fr::zero();
        let mut gamma_g_multiplier = E::Fr::zero();
        for (((c, z), v), proof) in commitments.iter().zip(points).zip(values).zip(proofs) {
            let mut c = c.0.into_projective();
            let w = proof.w.into_projective();
            c += &w.mul(z);
            g_multiplier += &(randomizer * &v);
            gamma_g_multiplier += &(randomizer * &proof.random_v);
            total_c += &c.mul(&randomizer);
            total_w += &w.mul(&randomizer);
            // We don't need to sample randomizers from the full field,
            // only from 128-bit strings.
            randomizer = u128::rand(rng).into();
        }
        total_c -= &g.mul(&g_multiplier);
        total_c -= &gamma_g.mul(&gamma_g_multiplier);
        end_timer!(combination_time);

        let to_affine_time = start_timer!(|| "Converting results to affine for pairing");
        let mut to_affine = [-total_w, total_c];
        E::G1Projective::batch_normalization(&mut to_affine);
        let [total_w, total_c] = to_affine;
        let total_w = total_w.into_affine();
        let total_c = total_c.into_affine();
        end_timer!(to_affine_time);

        let pairing_time = start_timer!(|| "Performing product of pairings");
        let result = E::product_of_pairings(&[
            (&total_w.prepare(), &vk.prepared_beta_h),
            (&total_c.prepare(), &vk.prepared_h),
        ]) == E::Fqk::one();
        end_timer!(pairing_time);
        end_timer!(check_time, || format!("Result: {}", result));
        Ok(result)
    }
}

fn skip_leading_zeros_and_convert_to_bigints<F: PrimeField>(p: &Polynomial<F>) -> (usize, Vec<F::BigInt>) {
    let mut num_leading_zeros = 0;
    while p.coeffs[num_leading_zeros].is_zero() && num_leading_zeros < p.coeffs.len() {
        num_leading_zeros += 1;
    }
    let coeffs = convert_to_bigints(&p.coeffs[num_leading_zeros..]);
    (num_leading_zeros, coeffs)
}

fn convert_to_bigints<F: PrimeField>(p: &[F]) -> Vec<F::BigInt> {
    let to_bigint_time = start_timer!(|| "Converting polynomial coeffs to bigints");
    let coeffs = p.par_iter().map(|s| s.into_repr()).collect::<Vec<_>>();
    end_timer!(to_bigint_time);
    coeffs
}

#[cfg(test)]
mod tests {
    #![allow(non_camel_case_types)]
    use crate::kzg10::*;
    use crate::*;
    use algebra::fields::bls12_381::Fr;

    use algebra::curves::bls12_377::Bls12_377;
    use algebra::curves::bls12_381::Bls12_381;
    use algebra::curves::mnt6::MNT6;
    use algebra::curves::sw6::SW6;

    use rand::thread_rng;

    impl<E: PairingEngine> KZG10<E> {
        /// Specializes the public parameters for a given maximum degree `d` for polynomials
        /// `d` should be less that `pp.max_degree()`.
        pub(crate) fn trim(
            pp: &UniversalParams<E>,
            mut supported_degree: usize,
        ) -> Result<(Powers<E>, VerifierKey<E>), Error> {
            if supported_degree == 1 {
                supported_degree += 1;
            }
            let powers_of_g = pp.powers_of_g[..=supported_degree].to_vec();
            let powers_of_gamma_g = pp.powers_of_gamma_g[..=supported_degree].to_vec();

            let powers = Powers {
                powers_of_g: std::borrow::Cow::Owned(powers_of_g),
                powers_of_gamma_g: std::borrow::Cow::Owned(powers_of_gamma_g),
            };
            let vk = VerifierKey {
                g: pp.powers_of_g[0],
                gamma_g: pp.powers_of_gamma_g[0],
                h: pp.h,
                beta_h: pp.beta_h,
                prepared_h: pp.prepared_h.clone(),
                prepared_beta_h: pp.prepared_beta_h.clone(),
            };
            Ok((powers, vk))
        }
    }

    #[test]
    fn add_commitments_test() {
        let rng = &mut thread_rng();
        let p = Polynomial::from_coefficients_slice(&[
            Fr::rand(rng),
            Fr::rand(rng),
            Fr::rand(rng),
            Fr::rand(rng),
            Fr::rand(rng),
        ]);
        let f = Fr::rand(rng);
        let mut f_p = Polynomial::zero();
        f_p += (f, &p);

        let degree = 4;
        let pp = KZG_Bls12_381::setup(degree, false, rng).unwrap();
        let (powers, _) = KZG_Bls12_381::trim(&pp, degree).unwrap();

        let hiding_bound = None;
        let (comm, _) = KZG10::commit(&powers, &p, hiding_bound, Some(rng)).unwrap();
        let (f_comm, _) = KZG10::commit(&powers, &f_p, hiding_bound, Some(rng)).unwrap();
        let mut f_comm_2 = Commitment::empty();
        f_comm_2 += (f, &comm);

        assert_eq!(f_comm, f_comm_2);
    }

    type KZG_Bls12_381 = KZG10<Bls12_381>;

    fn end_to_end_test_template<E: PairingEngine>() -> Result<(), Error> {
        let rng = &mut thread_rng();
        for _ in 0..100 {
            let mut degree = 0;
            while degree <= 1 {
                degree = usize::rand(rng) % 20;
            }
            let pp = KZG10::<E>::setup(degree, false, rng)?;
            let (ck, vk) = KZG10::trim(&pp, degree)?;
            let p = Polynomial::rand(degree, rng);
            let hiding_bound = Some(1);
            let (comm, rand) = KZG10::<E>::commit(&ck, &p, hiding_bound, Some(rng))?;
            let point = E::Fr::rand(rng);
            let value = p.evaluate(point);
            let proof = KZG10::<E>::open(&ck, &p, point, &rand)?;
            assert!(
                KZG10::<E>::check(&vk, &comm, point, value, &proof)?,
                "proof was incorrect for max_degree = {}, polynomial_degree = {}, hiding_bound = {:?}",
                degree,
                p.degree(),
                hiding_bound,
            );
        }
        Ok(())
    }

    fn linear_polynomial_test_template<E: PairingEngine>() -> Result<(), Error> {
        let rng = &mut thread_rng();
        for _ in 0..100 {
            let degree = 50;
            let pp = KZG10::<E>::setup(degree, false, rng)?;
            let (ck, vk) = KZG10::trim(&pp, 2)?;
            let p = Polynomial::rand(1, rng);;
            let hiding_bound = Some(1);
            let (comm, rand) = KZG10::<E>::commit(&ck, &p, hiding_bound, Some(rng))?;
            let point = E::Fr::rand(rng);
            let value = p.evaluate(point);
            let proof = KZG10::<E>::open(&ck, &p, point, &rand)?;
            assert!(
                KZG10::<E>::check(&vk, &comm, point, value, &proof)?,
                "proof was incorrect for max_degree = {}, polynomial_degree = {}, hiding_bound = {:?}",
                degree,
                p.degree(),
                hiding_bound,
            );
        }
        Ok(())
    }

    fn batch_check_test_template<E: PairingEngine>() -> Result<(), Error> {
        let rng = &mut thread_rng();
        for _ in 0..10 {
            let mut degree = 0;
            while degree <= 1 {
                degree = usize::rand(rng) % 20;
            }
            let pp = KZG10::<E>::setup(degree, false, rng)?;
            let (ck, vk) = KZG10::trim(&pp, degree)?;
            let mut comms = Vec::new();
            let mut values = Vec::new();
            let mut points = Vec::new();
            let mut proofs = Vec::new();
            for _ in 0..10 {
                let p = Polynomial::rand(degree, rng);;
                let hiding_bound = Some(1);
                let (comm, rand) = KZG10::<E>::commit(&ck, &p, hiding_bound, Some(rng))?;
                let point = E::Fr::rand(rng);
                let value = p.evaluate(point);
                let proof = KZG10::<E>::open(&ck, &p, point, &rand)?;

                assert!(KZG10::<E>::check(&vk, &comm, point, value, &proof)?);
                comms.push(comm);
                values.push(value);
                points.push(point);
                proofs.push(proof);
            }
            assert!(KZG10::<E>::batch_check(
                &vk, &comms, &points, &values, &proofs, rng
            )?);
        }
        Ok(())
    }

    #[test]
    fn end_to_end_test() {
        end_to_end_test_template::<Bls12_377>().expect("test failed for bls12-377");
        end_to_end_test_template::<Bls12_381>().expect("test failed for bls12-381");
        end_to_end_test_template::<MNT6>().expect("test failed for MNT6");
        end_to_end_test_template::<SW6>().expect("test failed for SW6");
    }

    #[test]
    fn linear_polynomial_test() {
        linear_polynomial_test_template::<Bls12_377>().expect("test failed for bls12-377");
        linear_polynomial_test_template::<Bls12_381>().expect("test failed for bls12-381");
        linear_polynomial_test_template::<MNT6>().expect("test failed for MNT6");
        linear_polynomial_test_template::<SW6>().expect("test failed for SW6");
    }
    #[test]
    fn batch_check_test() {
        batch_check_test_template::<Bls12_377>().expect("test failed for bls12-377");
        batch_check_test_template::<Bls12_381>().expect("test failed for bls12-381");
        batch_check_test_template::<MNT6>().expect("test failed for MNT6");
        batch_check_test_template::<SW6>().expect("test failed for SW6");
    }
}
