//! Here we construct a polynomial commitment that enables users to commit to a
//! single polynomial `p`, and then later provide an evaluation proof that
//! convinces verifiers that a claimed value `v` is the true evaluation of `p`
//! at a chosen point `x`. Our construction follows the template of the construction
//! proposed by Kate, Zaverucha, and Goldberg ([KZG10](http://cacr.uwaterloo.ca/techreports/2010/cacr2010-10.pdf)).
//! This construction achieves extractability in the algebraic group model (AGM).

use crate::{BTreeMap, Error, LabeledPolynomial, PCRandomness, ToString, Vec};
use ark_ec::AffineRepr;
use ark_ec::{pairing::Pairing, CurveGroup};
use ark_ec::{scalar_mul::fixed_base::FixedBase, VariableBaseMSM};
use ark_ff::{One, PrimeField, UniformRand, Zero};
use ark_poly::DenseUVPolynomial;
use ark_std::{format, marker::PhantomData, ops::Div, ops::Mul, vec};

use ark_std::rand::RngCore;
#[cfg(feature = "parallel")]
use rayon::prelude::*;

mod data_structures;
pub use data_structures::*;

/// `KZG10` is an implementation of the polynomial commitment scheme of
/// [Kate, Zaverucha and Goldbgerg][kzg10]
///
/// [kzg10]: http://cacr.uwaterloo.ca/techreports/2010/cacr2010-10.pdf
pub struct KZG10<E: Pairing, P: DenseUVPolynomial<E::ScalarField>> {
    _engine: PhantomData<E>,
    _poly: PhantomData<P>,
}

impl<E, P> KZG10<E, P>
where
    E: Pairing,
    P: DenseUVPolynomial<E::ScalarField, Point = E::ScalarField>,
    for<'a, 'b> &'a P: Div<&'b P, Output = P>,
{
    /// Constructs public parameters when given as input the maximum degree `degree`
    /// for the polynomial commitment scheme.
    ///
    /// # Examples
    ///
    /// ```
    /// use ark_poly_commit::kzg10::KZG10;
    /// use ark_bls12_381::Bls12_381;
    /// use ark_bls12_381::Fr;
    /// use ark_poly::univariate::DensePolynomial;
    /// use ark_ec::pairing::Pairing;
    /// use ark_std::test_rng;
    /// type UniPoly_381 = DensePolynomial<<Bls12_381 as Pairing>::ScalarField>;
    ///
    /// let rng = &mut test_rng();
    /// let params = KZG10::<Bls12_381, UniPoly_381>::setup(10, false, rng).expect("Setup failed");
    /// ```
    pub fn setup<R: RngCore>(
        max_degree: usize,
        produce_g2_powers: bool,
        rng: &mut R,
    ) -> Result<UniversalParams<E>, Error> {
        if max_degree < 1 {
            return Err(Error::DegreeIsZero);
        }
        let setup_time = start_timer!(|| format!("KZG10::Setup with degree {}", max_degree));
        let beta = E::ScalarField::rand(rng);
        let g = E::G1::rand(rng);
        let gamma_g = E::G1::rand(rng);
        let h = E::G2::rand(rng);

        let mut powers_of_beta = vec![E::ScalarField::one()];

        let mut cur = beta;
        for _ in 0..max_degree {
            powers_of_beta.push(cur);
            cur *= &beta;
        }

        let window_size = FixedBase::get_mul_window_size(max_degree + 1);

        let scalar_bits = E::ScalarField::MODULUS_BIT_SIZE as usize;
        let g_time = start_timer!(|| "Generating powers of G");
        let g_table = FixedBase::get_window_table(scalar_bits, window_size, g);
        let powers_of_g =
            FixedBase::msm::<E::G1>(scalar_bits, window_size, &g_table, &powers_of_beta);
        end_timer!(g_time);
        let gamma_g_time = start_timer!(|| "Generating powers of gamma * G");
        let gamma_g_table = FixedBase::get_window_table(scalar_bits, window_size, gamma_g);
        let mut powers_of_gamma_g =
            FixedBase::msm::<E::G1>(scalar_bits, window_size, &gamma_g_table, &powers_of_beta);
        // Add an additional power of gamma_g, because we want to be able to support
        // up to D queries.
        powers_of_gamma_g.push(powers_of_gamma_g.last().unwrap().mul(&beta));
        end_timer!(gamma_g_time);

        let powers_of_g = E::G1::normalize_batch(&powers_of_g);
        let powers_of_gamma_g = E::G1::normalize_batch(&powers_of_gamma_g)
            .into_iter()
            .enumerate()
            .collect();

        let neg_powers_of_h_time = start_timer!(|| "Generating negative powers of h in G2");
        let neg_powers_of_h = if produce_g2_powers {
            let mut neg_powers_of_beta = vec![E::ScalarField::one()];
            let mut cur = E::ScalarField::one() / &beta;
            for _ in 0..max_degree {
                neg_powers_of_beta.push(cur);
                cur /= &beta;
            }

            let neg_h_table = FixedBase::get_window_table(scalar_bits, window_size, h);
            let neg_powers_of_h = FixedBase::msm::<E::G2>(
                scalar_bits,
                window_size,
                &neg_h_table,
                &neg_powers_of_beta,
            );

            let affines = E::G2::normalize_batch(&neg_powers_of_h);
            let mut affines_map = BTreeMap::new();
            affines.into_iter().enumerate().for_each(|(i, a)| {
                affines_map.insert(i, a);
            });
            affines_map
        } else {
            BTreeMap::new()
        };

        end_timer!(neg_powers_of_h_time);

        let h = h.into_affine();
        let beta_h = h.mul(beta).into_affine();
        let prepared_h = h.into();
        let prepared_beta_h = beta_h.into();

        let pp = UniversalParams {
            powers_of_g,
            powers_of_gamma_g,
            h,
            beta_h,
            neg_powers_of_h,
            prepared_h,
            prepared_beta_h,
        };
        end_timer!(setup_time);
        Ok(pp)
    }

    /// Outputs a commitment to `polynomial`.
    ///
    /// # Examples
    ///
    /// ```
    /// use ark_poly_commit::kzg10::{KZG10, Powers};
    /// use ark_bls12_381::Bls12_381;
    /// use ark_bls12_381::Fr;
    /// use ark_poly::DenseUVPolynomial;
    /// use ark_poly::univariate::DensePolynomial;
    /// use ark_ec::pairing::Pairing;
    /// use ark_ec::AffineRepr;
    /// use ark_std::test_rng;
    /// use ark_std::Zero;
    /// type UniPoly_381 = DensePolynomial<<Bls12_381 as Pairing>::ScalarField>;
    ///
    /// let rng = &mut test_rng();
    /// let params = KZG10::<Bls12_381, UniPoly_381>::setup(10, false, rng).expect("Setup failed");
    /// let powers_of_g = params.powers_of_g[..=10].to_vec();
    /// let powers_of_gamma_g = (0..=10)
    ///     .map(|i| params.powers_of_gamma_g[&i])
    ///     .collect();
    /// let powers = Powers {
    ///     powers_of_g: ark_std::borrow::Cow::Owned(powers_of_g),
    ///     powers_of_gamma_g: ark_std::borrow::Cow::Owned(powers_of_gamma_g),
    /// };
    /// let secret_poly = UniPoly_381::rand(10, rng);
    /// let (comm, r) = KZG10::<Bls12_381, UniPoly_381>::commit(&powers, &secret_poly, None, None).expect("Commitment failed");
    /// assert!(!comm.0.is_zero(), "Commitment should not be zero");
    /// assert!(!r.is_hiding(), "Commitment should not be hiding");
    /// ```
    pub fn commit(
        powers: &Powers<E>,
        polynomial: &P,
        hiding_bound: Option<usize>,
        rng: Option<&mut dyn RngCore>,
    ) -> Result<(Commitment<E>, Randomness<E::ScalarField, P>), Error> {
        Self::check_degree_is_too_large(polynomial.degree(), powers.size())?;

        let commit_time = start_timer!(|| format!(
            "Committing to polynomial of degree {} with hiding_bound: {:?}",
            polynomial.degree(),
            hiding_bound,
        ));

        let (num_leading_zeros, plain_coeffs) =
            skip_leading_zeros_and_convert_to_bigints(polynomial);

        let msm_time = start_timer!(|| "MSM to compute commitment to plaintext poly");
        let mut commitment = <E::G1 as VariableBaseMSM>::msm_bigint(
            &powers.powers_of_g[num_leading_zeros..],
            &plain_coeffs,
        );
        end_timer!(msm_time);

        let mut randomness = Randomness::<E::ScalarField, P>::empty();
        if let Some(hiding_degree) = hiding_bound {
            let mut rng = rng.ok_or(Error::MissingRng)?;
            let sample_random_poly_time = start_timer!(|| format!(
                "Sampling a random polynomial of degree {}",
                hiding_degree
            ));

            randomness = Randomness::rand(hiding_degree, false, None, &mut rng);
            Self::check_hiding_bound(
                randomness.blinding_polynomial.degree(),
                powers.powers_of_gamma_g.len(),
            )?;
            end_timer!(sample_random_poly_time);
        }

        let random_ints = convert_to_bigints(&randomness.blinding_polynomial.coeffs());
        let msm_time = start_timer!(|| "MSM to compute commitment to random poly");
        let random_commitment = <E::G1 as VariableBaseMSM>::msm_bigint(
            &powers.powers_of_gamma_g,
            random_ints.as_slice(),
        )
        .into_affine();
        end_timer!(msm_time);

        commitment += &random_commitment;

        end_timer!(commit_time);
        Ok((Commitment(commitment.into()), randomness))
    }

    /// Compute witness polynomial.
    ///
    /// The witness polynomial w(x) the quotient of the division (p(x) - p(z)) / (x - z)
    /// Observe that this quotient does not change with z because
    /// p(z) is the remainder term. We can therefore omit p(z) when computing the quotient.
    pub fn compute_witness_polynomial(
        p: &P,
        point: P::Point,
        randomness: &Randomness<E::ScalarField, P>,
    ) -> Result<(P, Option<P>), Error> {
        let divisor = P::from_coefficients_vec(vec![-point, E::ScalarField::one()]);

        let witness_time = start_timer!(|| "Computing witness polynomial");
        let witness_polynomial = p / &divisor;
        end_timer!(witness_time);

        let random_witness_polynomial = if randomness.is_hiding() {
            let random_p = &randomness.blinding_polynomial;

            let witness_time = start_timer!(|| "Computing random witness polynomial");
            let random_witness_polynomial = random_p / &divisor;
            end_timer!(witness_time);
            Some(random_witness_polynomial)
        } else {
            None
        };

        Ok((witness_polynomial, random_witness_polynomial))
    }

    pub(crate) fn open_with_witness_polynomial<'a>(
        powers: &Powers<E>,
        point: P::Point,
        randomness: &Randomness<E::ScalarField, P>,
        witness_polynomial: &P,
        hiding_witness_polynomial: Option<&P>,
    ) -> Result<Proof<E>, Error> {
        Self::check_degree_is_too_large(witness_polynomial.degree(), powers.size())?;
        let (num_leading_zeros, witness_coeffs) =
            skip_leading_zeros_and_convert_to_bigints(witness_polynomial);

        let witness_comm_time = start_timer!(|| "Computing commitment to witness polynomial");
        let mut w = <E::G1 as VariableBaseMSM>::msm_bigint(
            &powers.powers_of_g[num_leading_zeros..],
            &witness_coeffs,
        );
        end_timer!(witness_comm_time);

        let random_v = if let Some(hiding_witness_polynomial) = hiding_witness_polynomial {
            let blinding_p = &randomness.blinding_polynomial;
            let blinding_eval_time = start_timer!(|| "Evaluating random polynomial");
            let blinding_evaluation = blinding_p.evaluate(&point);
            end_timer!(blinding_eval_time);

            let random_witness_coeffs = convert_to_bigints(&hiding_witness_polynomial.coeffs());
            let witness_comm_time =
                start_timer!(|| "Computing commitment to random witness polynomial");
            w += &<E::G1 as VariableBaseMSM>::msm_bigint(
                &powers.powers_of_gamma_g,
                &random_witness_coeffs,
            );
            end_timer!(witness_comm_time);
            Some(blinding_evaluation)
        } else {
            None
        };

        Ok(Proof {
            w: w.into_affine(),
            random_v,
        })
    }

    /// On input a polynomial `p` and a point `point`, outputs a proof for the same.
    pub(crate) fn open<'a>(
        powers: &Powers<E>,
        p: &P,
        point: P::Point,
        rand: &Randomness<E::ScalarField, P>,
    ) -> Result<Proof<E>, Error> {
        Self::check_degree_is_too_large(p.degree(), powers.size())?;
        let open_time = start_timer!(|| format!("Opening polynomial of degree {}", p.degree()));

        let witness_time = start_timer!(|| "Computing witness polynomials");
        let (witness_poly, hiding_witness_poly) = Self::compute_witness_polynomial(p, point, rand)?;
        end_timer!(witness_time);

        let proof = Self::open_with_witness_polynomial(
            powers,
            point,
            rand,
            &witness_poly,
            hiding_witness_poly.as_ref(),
        );

        end_timer!(open_time);
        proof
    }

    /// Verifies that `value` is the evaluation at `point` of the polynomial
    /// committed inside `comm`.
    pub fn check(
        vk: &VerifierKey<E>,
        comm: &Commitment<E>,
        point: E::ScalarField,
        value: E::ScalarField,
        proof: &Proof<E>,
    ) -> Result<bool, Error> {
        let check_time = start_timer!(|| "Checking evaluation");
        let mut inner = comm.0.into_group() - &vk.g.mul(value);
        if let Some(random_v) = proof.random_v {
            inner -= &vk.gamma_g.mul(random_v);
        }
        let lhs = E::pairing(inner, vk.h);

        let inner = vk.beta_h.into_group() - &vk.h.mul(point);
        let rhs = E::pairing(proof.w, inner);

        end_timer!(check_time, || format!("Result: {}", lhs == rhs));
        Ok(lhs == rhs)
    }

    /// Check that each `proof_i` in `proofs` is a valid proof of evaluation for
    /// `commitment_i` at `point_i`.
    pub fn batch_check<R: RngCore>(
        vk: &VerifierKey<E>,
        commitments: &[Commitment<E>],
        points: &[E::ScalarField],
        values: &[E::ScalarField],
        proofs: &[Proof<E>],
        rng: &mut R,
    ) -> Result<bool, Error> {
        let check_time =
            start_timer!(|| format!("Checking {} evaluation proofs", commitments.len()));

        let mut total_c = <E::G1>::zero();
        let mut total_w = <E::G1>::zero();

        let combination_time = start_timer!(|| "Combining commitments and proofs");
        let mut randomizer = E::ScalarField::one();
        // Instead of multiplying g and gamma_g in each turn, we simply accumulate
        // their coefficients and perform a final multiplication at the end.
        let mut g_multiplier = E::ScalarField::zero();
        let mut gamma_g_multiplier = E::ScalarField::zero();
        for (((c, z), v), proof) in commitments.iter().zip(points).zip(values).zip(proofs) {
            let w = proof.w;
            let mut temp = w.mul(*z);
            temp += &c.0;
            let c = temp;
            g_multiplier += &(randomizer * v);
            if let Some(random_v) = proof.random_v {
                gamma_g_multiplier += &(randomizer * &random_v);
            }
            total_c += &c.mul(randomizer);
            total_w += &w.mul(randomizer);
            // We don't need to sample randomizers from the full field,
            // only from 128-bit strings.
            randomizer = u128::rand(rng).into();
        }
        total_c -= &vk.g.mul(g_multiplier);
        total_c -= &vk.gamma_g.mul(gamma_g_multiplier);
        end_timer!(combination_time);

        let to_affine_time = start_timer!(|| "Converting results to affine for pairing");
        let affine_points = E::G1::normalize_batch(&[-total_w, total_c]);
        let (total_w, total_c) = (affine_points[0], affine_points[1]);
        end_timer!(to_affine_time);

        let pairing_time = start_timer!(|| "Performing product of pairings");
        let result = E::multi_pairing(
            [total_w, total_c],
            [vk.prepared_beta_h.clone(), vk.prepared_h.clone()],
        )
        .0
        .is_one();
        end_timer!(pairing_time);
        end_timer!(check_time, || format!("Result: {}", result));
        Ok(result)
    }

    pub(crate) fn check_degree_is_too_large(degree: usize, num_powers: usize) -> Result<(), Error> {
        let num_coefficients = degree + 1;
        if num_coefficients > num_powers {
            Err(Error::TooManyCoefficients {
                num_coefficients,
                num_powers,
            })
        } else {
            Ok(())
        }
    }

    pub(crate) fn check_hiding_bound(
        hiding_poly_degree: usize,
        num_powers: usize,
    ) -> Result<(), Error> {
        if hiding_poly_degree == 0 {
            Err(Error::HidingBoundIsZero)
        } else if hiding_poly_degree >= num_powers {
            // The above check uses `>=` because committing to a hiding poly with
            // degree `hiding_poly_degree` requires `hiding_poly_degree + 1`
            // powers.
            Err(Error::HidingBoundToolarge {
                hiding_poly_degree,
                num_powers,
            })
        } else {
            Ok(())
        }
    }

    pub(crate) fn check_degrees_and_bounds<'a>(
        supported_degree: usize,
        max_degree: usize,
        enforced_degree_bounds: Option<&[usize]>,
        p: &'a LabeledPolynomial<E::ScalarField, P>,
    ) -> Result<(), Error> {
        if let Some(bound) = p.degree_bound() {
            let enforced_degree_bounds =
                enforced_degree_bounds.ok_or(Error::UnsupportedDegreeBound(bound))?;

            if enforced_degree_bounds.binary_search(&bound).is_err() {
                Err(Error::UnsupportedDegreeBound(bound))
            } else if bound < p.degree() || bound > max_degree {
                return Err(Error::IncorrectDegreeBound {
                    poly_degree: p.degree(),
                    degree_bound: p.degree_bound().unwrap(),
                    supported_degree,
                    label: p.label().to_string(),
                });
            } else {
                Ok(())
            }
        } else {
            Ok(())
        }
    }
}

fn skip_leading_zeros_and_convert_to_bigints<F: PrimeField, P: DenseUVPolynomial<F>>(
    p: &P,
) -> (usize, Vec<F::BigInt>) {
    let mut num_leading_zeros = 0;
    while num_leading_zeros < p.coeffs().len() && p.coeffs()[num_leading_zeros].is_zero() {
        num_leading_zeros += 1;
    }
    let coeffs = convert_to_bigints(&p.coeffs()[num_leading_zeros..]);
    (num_leading_zeros, coeffs)
}

fn convert_to_bigints<F: PrimeField>(p: &[F]) -> Vec<F::BigInt> {
    let to_bigint_time = start_timer!(|| "Converting polynomial coeffs to bigints");
    let coeffs = ark_std::cfg_iter!(p)
        .map(|s| s.into_bigint())
        .collect::<Vec<_>>();
    end_timer!(to_bigint_time);
    coeffs
}

#[cfg(test)]
mod tests {
    #![allow(non_camel_case_types)]
    use crate::kzg10::*;
    use crate::*;

    use ark_bls12_377::Bls12_377;
    use ark_bls12_381::Bls12_381;
    use ark_bls12_381::Fr;
    use ark_ec::pairing::Pairing;
    use ark_poly::univariate::DensePolynomial as DensePoly;
    use ark_std::test_rng;

    type UniPoly_381 = DensePoly<<Bls12_381 as Pairing>::ScalarField>;
    type UniPoly_377 = DensePoly<<Bls12_377 as Pairing>::ScalarField>;
    type KZG_Bls12_381 = KZG10<Bls12_381, UniPoly_381>;

    impl<E: Pairing, P: DenseUVPolynomial<E::ScalarField>> KZG10<E, P> {
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
            let powers_of_gamma_g = (0..=supported_degree)
                .map(|i| pp.powers_of_gamma_g[&i])
                .collect();

            let powers = Powers {
                powers_of_g: ark_std::borrow::Cow::Owned(powers_of_g),
                powers_of_gamma_g: ark_std::borrow::Cow::Owned(powers_of_gamma_g),
            };
            let vk = VerifierKey {
                g: pp.powers_of_g[0],
                gamma_g: pp.powers_of_gamma_g[&0],
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
        let rng = &mut test_rng();
        let p = DensePoly::from_coefficients_slice(&[
            Fr::rand(rng),
            Fr::rand(rng),
            Fr::rand(rng),
            Fr::rand(rng),
            Fr::rand(rng),
        ]);
        let f = Fr::rand(rng);
        let mut f_p = DensePoly::zero();
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

    fn end_to_end_test_template<E, P>() -> Result<(), Error>
    where
        E: Pairing,
        P: DenseUVPolynomial<E::ScalarField, Point = E::ScalarField>,
        for<'a, 'b> &'a P: Div<&'b P, Output = P>,
    {
        let rng = &mut test_rng();
        for _ in 0..100 {
            let mut degree = 0;
            while degree <= 1 {
                degree = usize::rand(rng) % 20;
            }
            let pp = KZG10::<E, P>::setup(degree, false, rng)?;
            let (ck, vk) = KZG10::<E, P>::trim(&pp, degree)?;
            let p = P::rand(degree, rng);
            let hiding_bound = Some(1);
            let (comm, rand) = KZG10::<E, P>::commit(&ck, &p, hiding_bound, Some(rng))?;
            let point = E::ScalarField::rand(rng);
            let value = p.evaluate(&point);
            let proof = KZG10::<E, P>::open(&ck, &p, point, &rand)?;
            assert!(
                KZG10::<E, P>::check(&vk, &comm, point, value, &proof)?,
                "proof was incorrect for max_degree = {}, polynomial_degree = {}, hiding_bound = {:?}",
                degree,
                p.degree(),
                hiding_bound,
            );
        }
        Ok(())
    }

    fn linear_polynomial_test_template<E, P>() -> Result<(), Error>
    where
        E: Pairing,
        P: DenseUVPolynomial<E::ScalarField, Point = E::ScalarField>,
        for<'a, 'b> &'a P: Div<&'b P, Output = P>,
    {
        let rng = &mut test_rng();
        for _ in 0..100 {
            let degree = 50;
            let pp = KZG10::<E, P>::setup(degree, false, rng)?;
            let (ck, vk) = KZG10::<E, P>::trim(&pp, 2)?;
            let p = P::rand(1, rng);
            let hiding_bound = Some(1);
            let (comm, rand) = KZG10::<E, P>::commit(&ck, &p, hiding_bound, Some(rng))?;
            let point = E::ScalarField::rand(rng);
            let value = p.evaluate(&point);
            let proof = KZG10::<E, P>::open(&ck, &p, point, &rand)?;
            assert!(
                KZG10::<E, P>::check(&vk, &comm, point, value, &proof)?,
                "proof was incorrect for max_degree = {}, polynomial_degree = {}, hiding_bound = {:?}",
                degree,
                p.degree(),
                hiding_bound,
            );
        }
        Ok(())
    }

    fn batch_check_test_template<E, P>() -> Result<(), Error>
    where
        E: Pairing,
        P: DenseUVPolynomial<E::ScalarField, Point = E::ScalarField>,
        for<'a, 'b> &'a P: Div<&'b P, Output = P>,
    {
        let rng = &mut test_rng();
        for _ in 0..10 {
            let mut degree = 0;
            while degree <= 1 {
                degree = usize::rand(rng) % 20;
            }
            let pp = KZG10::<E, P>::setup(degree, false, rng)?;
            let (ck, vk) = KZG10::<E, P>::trim(&pp, degree)?;
            let mut comms = Vec::new();
            let mut values = Vec::new();
            let mut points = Vec::new();
            let mut proofs = Vec::new();
            for _ in 0..10 {
                let p = P::rand(degree, rng);
                let hiding_bound = Some(1);
                let (comm, rand) = KZG10::<E, P>::commit(&ck, &p, hiding_bound, Some(rng))?;
                let point = E::ScalarField::rand(rng);
                let value = p.evaluate(&point);
                let proof = KZG10::<E, P>::open(&ck, &p, point, &rand)?;

                assert!(KZG10::<E, P>::check(&vk, &comm, point, value, &proof)?);
                comms.push(comm);
                values.push(value);
                points.push(point);
                proofs.push(proof);
            }
            assert!(KZG10::<E, P>::batch_check(
                &vk, &comms, &points, &values, &proofs, rng
            )?);
        }
        Ok(())
    }

    #[test]
    fn end_to_end_test() {
        end_to_end_test_template::<Bls12_377, UniPoly_377>().expect("test failed for bls12-377");
        end_to_end_test_template::<Bls12_381, UniPoly_381>().expect("test failed for bls12-381");
    }

    #[test]
    fn linear_polynomial_test() {
        linear_polynomial_test_template::<Bls12_377, UniPoly_377>()
            .expect("test failed for bls12-377");
        linear_polynomial_test_template::<Bls12_381, UniPoly_381>()
            .expect("test failed for bls12-381");
    }
    #[test]
    fn batch_check_test() {
        batch_check_test_template::<Bls12_377, UniPoly_377>().expect("test failed for bls12-377");
        batch_check_test_template::<Bls12_381, UniPoly_381>().expect("test failed for bls12-381");
    }

    #[test]
    fn test_degree_is_too_large() {
        let rng = &mut test_rng();

        let max_degree = 123;
        let pp = KZG_Bls12_381::setup(max_degree, false, rng).unwrap();
        let (powers, _) = KZG_Bls12_381::trim(&pp, max_degree).unwrap();

        let p = DensePoly::<Fr>::rand(max_degree + 1, rng);
        assert!(p.degree() > max_degree);
        assert!(KZG_Bls12_381::check_degree_is_too_large(p.degree(), powers.size()).is_err());
    }
}
