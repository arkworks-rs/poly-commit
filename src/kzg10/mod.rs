//! Here we constuct a polynomial commitment that enables users to commit to a
//! single polynomial `p`, and then later provide an evaluation proof that
//! convinces verifiers that a claimed value `v` is the true evaluation of `p`
//! at a chosen point `x`. Our construction follows the template of the construction
//! proposed by Kate, Zaverucha, and Goldberg ([KZG10](http://cacr.uwaterloo.ca/techreports/2010/cacr2010-10.pdf)).
//! This construction achieves extractability in the algebraic group model (AGM).

use crate::{
    PCUniversalParams, PCCommitment, PCCommitterKey, PCRandomness, PCVerifierKey,
    Polynomial, SinglePolynomialCommitment,
};
use algebra::bytes::*;
use algebra::msm::{FixedBaseMSM, VariableBaseMSM};
use algebra::{
    AffineCurve, Field, Group, PairingCurve, PairingEngine, PrimeField, ProjectiveCurve,
    UniformRand,
};
use rand_core::RngCore;
use rayon::prelude::*;
use std::marker::PhantomData;
use std::ops::{Range, AddAssign};

mod data_structures;
pub use data_structures::*;

mod error;
pub use error::*;

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
        max_degree: usize,
        has_g2: bool,
        rng: &mut R,
    ) -> Result<UniversalParams<E>, Error> {
        if max_degree < 1 {
            return Err(Error::UnsupportedDegree);
        }
        let setup_time = start_timer!(|| format!("Started KZG10::Setup with degree {}", degree));
        let beta = E::Fr::rand(rng);
        let g = E::G1Projective::rand(rng);
        let gamma_g = E::G1Projective::rand(rng);
        let h = E::G2Projective::rand(rng);

        let mut powers_of_beta = vec![E::Fr::one()];
        let mut cur = beta;
        // TODO: can optimize by computing
        // beta, beta^2, beta^4, ..., beta^d and using those to compute arbitrary
        // powers of beta faster via a bit decomposition of each scalar.
        // But that seems like it would be slower? log(D) additions instead of
        // one multiplication.
        // Anyway it shouldn't be a bottleneck?
        // If we use MPC to run the setup, we can consider this optimization.
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
            max_degree,
            h,
            beta_h,
            prepared_h,
            prepared_beta_h,
        };
        end_timer!(setup_time);
        Ok(pp)
    }

    pub fn trim(
        pp: &UniversalParams<E>,
        max_degree: usize,
    ) -> Result<(CommitterKey<E>, VerifierKey<E>), Error> {
        let mut powers_of_g = Vec::new();
        let mut powers_of_gamma_g = Vec::new();
        let mut max_degree = 0;
        let mut suppoer = Vec::new();

        let powers_of_g = pp.powers_of_g[..=max_degree].to_vec();
        let powers_of_gamma_g = pp.powers_of_gamma_g[..=max_degree].to_vec();

        let max_degree = pp.max_degree;

        let start = 0;
        let ck = CommitterKey {
            powers_of_g,
            powers_of_gamma_g,
            supported_degree_bounds: degree_bounds_to_support.to_vec(),
            max_degree,
        };
        let vk = VerifierKey {
            g: pp.powers_of_g[0],
            gamma_g: pp.powers_of_gamma_g[0],
            h: pp.h,
            beta_h: pp.beta_h,
            prepared_h: pp.prepared_h,
            prepared_beta_h: pp.prepared_beta_h,
            max_degree,
        };
        Ok((ck, vk))
    }

    /// Outputs a commitment to `polynomial`.
    pub fn commit(
        ck: &CommitterKey<E>,
        polynomial: &Polynomial<E::Fr>,
        hiding_bound: Option<usize>,
        rng: Option<&mut dyn RngCore>,
    ) -> Result<(Commitment<E>, Randomness<E>), Error> {
        Error::check_degree(polynomial.degree(), ck.max_degree())?;

        let commit_time = start_timer!(|| format!(
            "Committing to polynomial of degree {} with hiding_bound: {:?}",
            polynomial.degree(),
            hiding_bound,
        ));

        let mut skip_leading_zeros = 0;
        while polynomial.coeffs[skip_leading_zeros].is_zero() && skip_leading_zeros < polynomial.coeffs.len() {
            skip_leading_zeros += 1;
        }

        let from_mont_repr_time =
            start_timer!(|| "Converting plaintext polynomial from Montgomery repr");
        let plain_coeffs = polynomial.coeffs[skip_leading_zeros..]
            .par_iter()
            .map(|s| s.into_repr())
            .collect::<Vec<_>>();
        end_timer!(from_mont_repr_time);

        let msm_time = start_timer!(|| "MSM to compute commitment to plaintext poly");
        let mut commitment = VariableBaseMSM::multi_scalar_mul(&ck.powers_of_g[skip_leading_zeros..], &plain_coeffs);
        end_timer!(msm_time);

        let mut randomness = Randomness::empty();
        if let Some(hiding_degree) = hiding_bound {
            let mut rng = rng.ok_or(Error::MissingRng)?;
            let sample_random_poly_time = start_timer!(|| format!(
                "Sampling a random polynomial of degree {}",
                hiding_degree
            ));

            randomness = Randomness::rand(hiding_degree, &mut rng);
            if randomness.blinding_polynomial.degree() > ck.max_degree() {
                eprintln!("The hiding bound is too large for the commitment key.");
                Err(Error::PolynomialDegreeTooLarge {
                    poly_degree: randomness.blinding_polynomial.degree(),
                    max_degree: ck.max_degree(),
                })?;
            }
            end_timer!(sample_random_poly_time);
        }

        let from_mont_repr_time =
            start_timer!(|| "Converting random polynomial from Montgomery repr");
        let random_ints = randomness
            .blinding_polynomial
            .coeffs
            .par_iter()
            .map(|s: &E::Fr| s.into_repr())
            .collect::<Vec<_>>();
        end_timer!(from_mont_repr_time);

        let msm_time = start_timer!(|| "MSM to compute commitment to random poly");
        let random_commitment =
            VariableBaseMSM::multi_scalar_mul(&ck.powers_of_gamma_g, random_ints.as_slice())
                .into_affine();
        end_timer!(msm_time);

        commitment.add_assign_mixed(&random_commitment);

        end_timer!(commit_time);
        Ok((Commitment(commitment.into()), randomness))
    }

    /// Outputs a commitments to `polynomials`.
    pub fn multi_commit<'a>(
        ck: &CommitterKey<E>,
        polynomials: impl IntoIterator<Item = &'a LabeledPolynomial<'a, E::Fr>>,
        rng: Option<&mut dyn RngCore>,
    ) -> Result<
        (
            Vec<LabeledCommitment<Self::Commitment>>,
            Vec<Self::Randomness>,
        ),
        Error,
    > {
        let commit_time = start_timer!(|| "Committing to polynomials");

        let mut commitments = Vec::new();
        let mut randomness = Vec::new();
        let max_degree = ck.max_degree();

        let rng = &mut optional_rng::OptionalRng(rng);
        for polynomial in polynomials {
            let label = polynomial.label();
            let degree_bound = polynomial.degree_bound();
            let hiding_bound = polynomial.hiding_bound();
            let polynomial = polynomial.polynomial();

            Error::check_degrees(
                polynomial.degree(),
                degree_bound,
                max_degree,
                label.to_string(),
            )?;

            let commit_time = start_timer!(|| format!(
                "Polynomial {} of degree {}, degree bound {:?}, and hiding bound {:?}",
                label,
                polynomial.degree(),
                degree_bound,
                hiding_bound,
            ));
            let (comm, rand) = Self::commit(ck, polynomial, hiding_bound, Some(rng))?;
            let (shifted_comm, shifted_rand) = if let Some(degree_bound) = degree_bound {
                if degree_bound < max_degree {
                    let s_polynomial = shift_polynomial(polynomial, degree_bound, max_degree);
                    assert!(
                        polynomial.degree() <= s_polynomial.degree()
                            && s_polynomial.degree() <= max_degree
                            && s_polynomial.degree()
                                == polynomial.degree() + max_degree - degree_bound,
                        "polynomial.degree(): {}; s_polynomial.degree(): {}; max_degree: {}.",
                        polynomial.degree(),
                        s_polynomial.degree(),
                        max_degree
                    );
                    let (shifted_comm, shifted_rand) =
                        Self::commit(ck, &s_polynomial, hiding_bound, Some(rng))?;
                    (Some(shifted_comm), Some(shifted_rand))
                } else {
                    (None, None)
                }
            } else {
                (None, None)
            };

            let comm = Commitment { comm, shifted_comm };
            let rand = Randomness { rand, shifted_rand };
            commitments.push(LabeledCommitment::new(label.to_string(), comm, degree_bound));
            randomness.push(rand);
            end_timer!(commit_time);
        }
        end_timer!(commit_time);
        Ok((commitments, randomness))
    }

    fn compute_witness_polynomial(
        polynomial: &Polynomial<E::Fr>,
        point: E::Fr,
        randomness: &Randomness<E>,
    ) -> Result<(Polynomial<E::Fr>, Option<Polynomial<E::Fr>>), Error> {
        let eval_time = start_timer!(|| "Evaluating polynomial");
        let value = p.evaluate(point);
        end_timer!(eval_time);

        let witness_time = start_timer!(|| "Computing witness polynomial");
        let witness_polynomial = &(p - &Polynomial::from_coefficients_vec(vec![value]))
            / &Polynomial::from_coefficients_vec(vec![-point, E::Fr::one()]);
        timer_end!(witness_time)

        let random_witness_polynomial = if randomness.is_hiding() {
            let random_p = &randomness.blinding_polynomial;

            let rand_eval_time = start_timer!(|| "Evaluating random polynomial");
            let random_value = random_p.evaluate(point);
            end_timer!(rand_eval_time);

            let witness_time = start_timer!(|| "Computing random witness polynomial");
            let random_witness_polynomial = &(random_p
                - &Polynomial::from_coefficients_vec(vec![random_value]))
                / &Polynomial::from_coefficients_vec(vec![-point, E::Fr::one()]);
            end_timer!(witness_time);
            Some(random_witness_polynomial)
        } else {
            None
        };

        Ok((witness_time, random_witness_polynomial))
    }

    /// On input a polynomial `p` and a point `point`, outputs a proof for the same.
    pub(crate) fn open<'a>(
        ck: &Self::CommitterKey,
        polynomial: &Polynomial<E::Fr>,
        point: E::Fr,
        rand: &Randomness<E>,
    ) -> Result<Proof<E>, Error> {
        Error::check_degree(p.degree(), ck.max_degree())?;
        let open_time = start_timer!(|| format!("Opening polynomial of degree {}", p.degree()));

        let witness_time = start_timer!(|| "Computing witness polynomials");
        let (witness_poly, hiding_witness_poly) = Self::compute_witness_polynomial(polynomial, point, rand)?;
        end_timer!(witness_time);

        let proof = Self::open_with_witness_polynomial(ck, point, rand, &witness_poly, hiding_witness_poly.as_ref());

        end_timer!(open_time);
        proof
    }

    pub(crate) fn open_with_witness_polynomial<'a>(
        ck: &Self::CommitterKey,
        point: E::Fr,
        randomness: &Randomness<E>,
        witness_polynomial: &Polynomial<E::Fr>,
        hiding_witness_polynomial: Option<&Polynomial<E::Fr>>,
    ) -> Result<Proof, Error> {
        let convert_time = start_timer!(|| "Converting witness polynomial from Montgomery repr");
        let mut skip_leading_zeros = 0;
        while witness_polynomial.coeffs[skip_leading_zeros].is_zero() && skip_leading_zeros < polynomial.coeffs.len() {
            skip_leading_zeros += 1;
        }
        let witness_coeffs = witness_polynomial[skip_leading_zeros..]
            .coeffs
            .par_iter()
            .map(|s| s.into_repr())
            .collect::<Vec<_>>();
        end_timer!(convert_time);

        let witness_comm_time = start_timer!(|| "Computing commitment to witness polynomial");
        let mut w = VariableBaseMSM::multi_scalar_mul(&ck.powers_of_g[skip_leading_zeros..], &witness_coeffs);
        end_timer!(witness_comm_time);

        let mut random_v = E::Fr::zero();
        if let Some(hiding_witness_polynomial) = hiding_witness_poly {
            assert_eq!(randomness.is_some(), "hiding_witness_polynomial is some, but randomness is none");
            let blinding_p = &randomness.blinding_polynomial;
            let blinding_eval_time = start_timer!(|| "Evaluating random polynomial");
            let blinding_evaluation = blinding_p.evaluate(point);
            end_timer!(blinding_eval_time);

            let random_witness_coeffs = random_witness_polynomial
                .coeffs
                .into_par_iter()
                .map(|s| s.into_repr())
                .collect::<Vec<_>>();

            let witness_comm_time =
                start_timer!(|| "Computing commitment to random witness polynomial");
            w += &VariableBaseMSM::multi_scalar_mul(&ck.powers_of_gamma_g, &random_witness_coeffs);
            end_timer!(witness_comm_time);
            random_v = blinding_evaluation;
        }

        Ok(Proof {
            w: w.into_affine(),
            random_v,
        })
    }

    /// Verifies that `value` is the evaluation at `point` of the polynomial
    /// committed inside `comm`.
    pub fn check(
        vk: &Self::VerifierKey,
        comm: &Self::Commitment,
        point: E::Fr,
        value: E::Fr,
        proof: &Self::Proof,
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

    pub fn batch_check<R: RngCore>(
        vk: &Self::VerifierKey,
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

impl<E: PairingEngine> SinglePolynomialCommitment<E::Fr> for KZG10<E> {
    type UniversalParams = UniversalParams<E>;
    type CommitterKey = CommitterKey<E>;
    type VerifierKey = VerifierKey<E>;
    type Commitment = Commitment<E>;
    type Randomness = Randomness<E>;
    type Proof = Proof<E>;
    type Error = Error;

    /// Constructs public parameters when given as input the maximum degree `degree`
    /// for the polynomial commitment scheme.
    fn setup<R: RngCore>(
        max_degree: usize,
        rng: &mut R,
    ) -> Result<Self::UniversalParams, Self::Error> {
        Self::setup(max_degree, rng)
    }

    fn trim(
        pp: &Self::UniversalParams,
        degree: usize,
    ) -> Result<(Self::CommitterKey, Self::VerifierKey), Self::Error> {
        Self::trim(pp, degree)
    }

    /// Outputs a commitment to `polynomial`.
    fn commit<'a>(
        ck: &Self::CommitterKey,
        polynomials: impl IntoIterator<Item = &'a LabeledPolynomial<'a, F>>,
        rng: Option<&mut dyn RngCore>,
    ) -> Result<
        (
            Vec<LabeledCommitment<Self::Commitment>>,
            Vec<Self::Randomness>,
        ),
        Self::Error,
    > {
        Self::multi_commit(ck, polynomials, rng)
    }


}

#[cfg(test)]
mod tests {
    #![allow(non_camel_case_types)]
    use crate::single_pc::kzg10::*;
    use crate::*;
    use algebra::fields::bls12_381::Fr;

    use algebra::curves::bls12_377::Bls12_377;
    use algebra::curves::bls12_381::Bls12_381;
    use algebra::curves::mnt6::MNT6;
    use algebra::curves::sw6::SW6;

    use rand::thread_rng;

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
        let pp = KZG_Bls12_381::setup(degree, rng).unwrap();
        let support = CoefficientSupport::from_dual_interval(4, 4);
        let (ck, _) = KZG_Bls12_381::trim(degree).unwrap();

        let hiding_bound = None;
        let (comm, _) = KZG10::commit(&ck, &p, hiding_bound, Some(rng)).unwrap();
        let (f_comm, _) = KZG10::commit(&ck, &f_p, hiding_bound, Some(rng)).unwrap();
        let mut f_comm_2 = Commitment::empty();
        f_comm_2 += (f, &comm);

        assert_eq!(f_comm, f_comm_2);
    }
    type KZG_Bls12_381 = KZG10<Bls12_381>;
    type KZG_Bls12_377 = KZG10<Bls12_377>;
    type KZG_MNT6 = KZG10<MNT6>;
    type KZG_SW6 = KZG10<SW6>;

    #[test]
    fn end_to_end_test() {
        use crate::single_pc::tests::*;

        end_to_end_test::<_, KZG_Bls12_377>().expect("test failed for bls12-377");
        end_to_end_test::<_, KZG_Bls12_381>().expect("test failed for bls12-381");
        end_to_end_test::<_, KZG_MNT6>().expect("test failed for MNT6");
        end_to_end_test::<_, KZG_SW6>().expect("test failed for SW6");
    }

    #[test]
    fn linear_polynomial_test() {
        use crate::single_pc::tests::*;

        linear_polynomial_test::<_, KZG_Bls12_377>().expect("test failed for bls12-377");
        linear_polynomial_test::<_, KZG_Bls12_381>().expect("test failed for bls12-381");
        linear_polynomial_test::<_, KZG_MNT6>().expect("test failed for MNT6");
        linear_polynomial_test::<_, KZG_SW6>().expect("test failed for SW6");
    }
    #[test]
    fn batch_check_test() {
        use crate::single_pc::tests::*;

        batch_check_test::<_, KZG_Bls12_377>().expect("test failed for bls12-377");
        batch_check_test::<_, KZG_Bls12_381>().expect("test failed for bls12-381");
        batch_check_test::<_, KZG_MNT6>().expect("test failed for MNT6");
        batch_check_test::<_, KZG_SW6>().expect("test failed for SW6");
    }
}
