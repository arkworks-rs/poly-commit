use crate::{PCUniversalParams, PCRandomness, Polynomial, PolynomialCommitment};
use crate::{QuerySetError, EquationError, QuerySet, Evaluations};
use crate::{LabeledPolynomial, LabeledCommitment, Equation};
use crate::kzg10;

use algebra::{AffineCurve, Field, PairingEngine, ProjectiveCurve, One, Zero};
use rand_core::RngCore;
use std::marker::PhantomData;
use std::collections::{BTreeMap, BTreeSet};

mod data_structures;
pub use data_structures::*;

mod error;
pub use error::*;

/// `MarlinKZG10` is an implementation of the polynomial commitment scheme of
/// [Kate, Zaverucha and Goldbgerg][kzg10], with degree bound enforcement as
/// described in the [Marlin paper][marlin]
///
/// [kzg10]: http://cacr.uwaterloo.ca/techreports/2010/cacr2010-10.pdf
/// [marlin]: https://eprint.iacr.org/2019/1047
pub struct MarlinKZG10<E: PairingEngine> {
    _engine: PhantomData<E>,
}

pub(crate) fn shift_polynomial<E: PairingEngine>(
    ck: &CommitterKey<E>,
    p: &Polynomial<E::Fr>,
    degree_bound: usize,
) -> Polynomial<E::Fr> {
    if p.is_zero() {
        Polynomial::zero()
    } else {
        let enforced_degree_bounds = ck.enforced_degree_bounds.as_ref().expect("Polynomial requires degree bounds, but `ck` does not support any");
        let largest_enforced_degree_bound = enforced_degree_bounds.last().unwrap();

        let mut shifted_polynomial_coeffs = vec![E::Fr::zero(); largest_enforced_degree_bound - degree_bound];
        shifted_polynomial_coeffs.extend_from_slice(&p.coeffs);
        Polynomial::from_coefficients_vec(shifted_polynomial_coeffs)
    }
}

impl<E: PairingEngine> MarlinKZG10<E> {
    /// MSM for `commitments` and `coeffs`
    fn combine_commitments<'a>(
        coeffs_and_comms: impl IntoIterator<Item = (E::Fr, &'a Commitment<E>)>,
    ) -> (E::G1Projective, E::G1Projective) {
        let mut combined_comm = E::G1Projective::zero();
        let mut combined_shifted_comm = E::G1Projective::zero();
        for (coeff, comm) in coeffs_and_comms {
            combined_comm += &comm.comm.0.mul(coeff);
            if let Some(shifted_comm) = &comm.shifted_comm {
                combined_shifted_comm += &shifted_comm.0.mul(coeff);
            }
        }
        (combined_comm, combined_shifted_comm)
    }

    /// Accumulate `commitments` and `values` according to `opening_challenge`.
    fn accumulate_commitments_and_values<'a>(
        vk: &VerifierKey<E>,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Commitment<E>>>,
        values: impl IntoIterator<Item = E::Fr>,
        opening_challenge: E::Fr,
    ) -> Result<(E::G1Projective, E::Fr), Error> {
        let acc_time = start_timer!(|| "Accumulating commitments and values");
        let mut combined_comm = E::G1Projective::zero();
        let mut combined_value = E::Fr::zero();
        let mut challenge_i = E::Fr::one();
        for (labeled_commitment, value) in commitments.into_iter().zip(values) {
            let degree_bound = labeled_commitment.degree_bound();
            let commitment = labeled_commitment.commitment();
            assert_eq!(degree_bound.is_some(), commitment.shifted_comm.is_some());

            combined_comm += &commitment.comm.0.mul(challenge_i);
            combined_value += &(value * &challenge_i);

            if let Some(degree_bound) = degree_bound {
                let challenge_i_1 = challenge_i * &opening_challenge;
                let shifted_comm = commitment.shifted_comm.as_ref().unwrap().0.into_projective();

                let shift_power = vk.get_shift_power(degree_bound).ok_or(Error::UnsupportedDegreeBound(degree_bound))?;
                let mut adjusted_comm = shifted_comm - &shift_power.mul(value);
                adjusted_comm.mul_assign(challenge_i_1);
                combined_comm += &adjusted_comm;
            }
            challenge_i *= &opening_challenge.square();
        }

        end_timer!(acc_time);
        Ok((combined_comm, combined_value))
    }

}

impl<E: PairingEngine> PolynomialCommitment<E::Fr> for MarlinKZG10<E> {
    type UniversalParams = UniversalParams<E>;
    type CommitterKey = CommitterKey<E>;
    type VerifierKey = VerifierKey<E>;
    type Commitment = Commitment<E>;
    type Randomness = Randomness<E>;
    type Proof = kzg10::Proof<E>;
    type BatchProof = Vec<Self::Proof>;
    type Error = Error;

    /// Constructs public parameters when given as input the maximum degree `max_degree`
    /// for the polynomial commitment scheme.
    fn setup<R: RngCore>(
        max_degree: usize,
        rng: &mut R,
    ) -> Result<Self::UniversalParams, Self::Error> {
        kzg10::KZG10::setup(max_degree, false, rng).map_err(Into::into)
    }

    // TODO: should trim also take in the hiding_bounds? That way we don't
    // have to store many powers of gamma_g.
    // TODO: add an optional hiding_bound.
    fn trim(
        pp: &Self::UniversalParams,
        supported_degree: usize,
        enforced_degree_bounds: Option<&[usize]>,
    ) -> Result<(Self::CommitterKey, Self::VerifierKey), Self::Error> {
        let max_degree = pp.max_degree();
        if supported_degree > max_degree {
            return Err(Error::TrimmingDegreeTooLarge)
        }

        // Construct the KZG10 committer key for committing to unshifted polynomials.
        let ck_time = start_timer!(|| format!("Constructing `powers` of size {} for unshifted polys", supported_degree));
        let powers = pp.powers_of_g[..=supported_degree].to_vec();
        // We want to support making up to supported_degree queries to committed
        // polynomials.
        let powers_of_gamma_g = pp.powers_of_gamma_g[..=(supported_degree + 1)].to_vec();
        end_timer!(ck_time);

        // Construct the core KZG10 verifier key.
        let vk = kzg10::VerifierKey {
            g: pp.powers_of_g[0],
            gamma_g: pp.powers_of_gamma_g[0],
            h: pp.h,
            beta_h: pp.beta_h,
            prepared_h: pp.prepared_h.clone(),
            prepared_beta_h: pp.prepared_beta_h.clone(),
        };

        let enforced_degree_bounds = enforced_degree_bounds.map(|v| {
            let mut v = v.to_vec();
            v.sort();
            v.dedup();
            v
        });

        // Check whether we have some degree bounds to enforce
        let (shifted_powers, degree_bounds_and_shift_powers) =
            if let Some(enforced_degree_bounds) = enforced_degree_bounds.as_ref() {
                if enforced_degree_bounds.is_empty() {
                    (None, None)
                } else {
                    let lowest_shifted_power = max_degree - enforced_degree_bounds
                        .last()
                        .ok_or(Error::EmptyDegreeBounds)?;

                    let shifted_ck_time = start_timer!(|| format!(
                            "Constructing `shifted_powers` of size {}",
                            max_degree - lowest_shifted_power + 1
                    ));

                    let shifted_powers = pp.powers_of_g[lowest_shifted_power..].to_vec();
                    end_timer!(shifted_ck_time);

                    let degree_bounds_and_shift_powers = enforced_degree_bounds
                        .iter()
                        .map(|d| (*d, pp.powers_of_g[max_degree - d]))
                        .collect();
                    (Some(shifted_powers), Some(degree_bounds_and_shift_powers))
                }
            } else {
                (None, None)
            };


        let ck = CommitterKey {
            powers,
            shifted_powers,
            powers_of_gamma_g,
            enforced_degree_bounds: enforced_degree_bounds,
            max_degree,
        };

        let vk = VerifierKey {
            vk,
            degree_bounds_and_shift_powers,
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
        let commit_time = start_timer!(|| "Committing to polynomials");

        let mut commitments = Vec::new();
        let mut randomness = Vec::new();

        let rng = &mut kzg10::optional_rng::OptionalRng(rng);
        for p in polynomials {
            let label = p.label();
            let degree_bound = p.degree_bound();
            let hiding_bound = p.hiding_bound();
            let polynomial = p.polynomial();

            Self::Error::check_degrees_and_bounds(&ck, &p)?;

            let commit_time = start_timer!(|| format!(
                "Polynomial {} of degree {}, degree bound {:?}, and hiding bound {:?}",
                label,
                polynomial.degree(),
                degree_bound,
                hiding_bound,
            ));

            let (comm, rand) = kzg10::KZG10::commit(&ck.powers(), polynomial, hiding_bound, Some(rng))?;
            let (shifted_comm, shifted_rand) = if let Some(degree_bound) = degree_bound {
                let shifted_powers = ck.shifted_powers(degree_bound).ok_or(Error::UnsupportedDegreeBound(degree_bound))?;
                let (shifted_comm, shifted_rand) =
                    kzg10::KZG10::commit(&shifted_powers, &polynomial, hiding_bound, Some(rng))?;
                (Some(shifted_comm), Some(shifted_rand))
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
        let mut p = Polynomial::zero();
        let mut r = kzg10::Randomness::empty();
        let mut shifted_w = Polynomial::zero();
        let mut shifted_r = kzg10::Randomness::empty();
        let mut shifted_r_witness = Polynomial::zero();

        let mut enforce_degree_bound = false;
        for (j, (polynomial, rand)) in labeled_polynomials.into_iter().zip(rands).enumerate() {
            let degree_bound = polynomial.degree_bound();

            Self::Error::check_degrees_and_bounds(&ck, &polynomial)?;

            // compute challenge^j and challenge^{j+1}.
            let challenge_j = opening_challenge.pow([2 * j as u64]);

            assert_eq!(degree_bound.is_some(), rand.shifted_rand.is_some());

            p += (challenge_j, polynomial.polynomial());
            r += (challenge_j, &rand.rand);

            if let Some(degree_bound) = degree_bound {
                enforce_degree_bound = true;
                let shifted_rand = rand.shifted_rand.as_ref().unwrap();
                let (witness, shifted_rand_witness) = kzg10::KZG10::compute_witness_polynomial(
                    polynomial.polynomial(),
                    point,
                    &shifted_rand,
                )?;
                let challenge_j_1 = challenge_j * &opening_challenge;

                let shifted_witness = shift_polynomial(ck, &witness, degree_bound);

                shifted_w += (challenge_j_1, &shifted_witness);
                shifted_r += (challenge_j_1, shifted_rand);
                if let Some(shifted_rand_witness) = shifted_rand_witness {
                    shifted_r_witness += (challenge_j_1, &shifted_rand_witness);
                }
            }
        }
        let proof_time = start_timer!(|| "Creating proof for unshifted polynomials");
        let proof = kzg10::KZG10::open(&ck.powers(), &p, point, &r)?;
        let mut w = proof.w.into_projective();
        let mut random_v = proof.random_v;
        end_timer!(proof_time);

        if enforce_degree_bound {
            let proof_time = start_timer!(|| "Creating proof for shifted polynomials");
            let shifted_proof = kzg10::KZG10::open_with_witness_polynomial(
                &ck.shifted_powers(None).unwrap(),
                point,
                &shifted_r,
                &shifted_w,
                Some(&shifted_r_witness),
            )?;
            end_timer!(proof_time);

            w += &shifted_proof.w.into_projective();
            random_v += &shifted_proof.random_v;
        }

        Ok(kzg10::Proof {
            w: w.into_affine(),
            random_v,
        })
    }


    /// Verifies that `value` is the evaluation at `x` of the polynomial
    /// committed inside `comm`.
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
        let check_time = start_timer!(|| "Checking evaluations");
        let (combined_comm, combined_value) =
            Self::accumulate_commitments_and_values(
                vk,
                commitments,
                values,
                opening_challenge,
            )?;
        let combined_comm = kzg10::Commitment(combined_comm.into());
        let result = kzg10::KZG10::check(&vk.vk, &combined_comm, point, combined_value, proof)?;
        end_timer!(check_time);
        Ok(result)
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

        let mut combined_comms = Vec::new();
        let mut combined_queries = Vec::new();
        let mut combined_evals = Vec::new();
        for (query, labels) in query_to_labels_map.into_iter() {
            let lc_time =
                start_timer!(|| format!("Randomly combining {} commitments", labels.len()));
            let mut comms_to_combine: Vec<&'_ LabeledCommitment<_>> = Vec::new();
            let mut values_to_combine = Vec::new();
            for label in labels.into_iter() {
                let commitment = commitments
                    .get(label)
                    .ok_or(QuerySetError::MissingPolynomial { label: label.to_string() })?;
                let degree_bound = commitment.degree_bound();
                assert_eq!(degree_bound.is_some(), commitment.commitment().shifted_comm.is_some());

                let v_i = values
                    .get(&(label.clone(), *query))
                    .ok_or(QuerySetError::MissingEvaluation { label: label.to_string() })?;

                comms_to_combine.push(commitment);
                values_to_combine.push(*v_i);
            }
            let (c, v) = Self::accumulate_commitments_and_values(
                vk,
                comms_to_combine,
                values_to_combine,
                opening_challenge,
            )?;
            end_timer!(lc_time);
            combined_comms.push(c);
            combined_queries.push(*query);
            combined_evals.push(v);
        }
        E::G1Projective::batch_normalization(&mut combined_comms);
        let combined_comms = combined_comms.into_iter().map(|c| kzg10::Commitment(c.into())).collect::<Vec<_>>();
        let proof_time = start_timer!(|| "Checking SinglePC::Proof");
        let result = kzg10::KZG10::batch_check(
            &vk.vk,
            &combined_comms,
            &combined_queries,
            &combined_evals,
            &proof,
            rng,
        )?;
        end_timer!(proof_time);
        Ok(result)
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
        let label_poly_rand_map =
            labeled_polynomials
            .into_iter()
            .zip(rands)
            .map(|(p, r)| (p.label(), (p, r)))
            .collect::<BTreeMap<_, _>>();

        let mut equation_polynomials = Vec::new();
        let mut equation_randomness = Vec::new();
        let mut query_set = QuerySet::new();
        let mut evaluations = Evaluations::new();
        for equation in equations {
            let equation_label = equation.label.clone();
            let mut poly = Polynomial::zero();
            let mut degree_bound = None;
            let mut hiding_bound = None;
            let mut randomness = Self::Randomness::empty();
            assert!(randomness.shifted_rand.is_none());

            let num_polys = equation.lhs.len();
            if equation.lhs_is_empty() {
                Err(EquationError::MissingLHS { label: equation_label.clone() })?;
            }

            for (coeff, label) in &equation.lhs {
                let &(cur_poly, cur_rand) = label_poly_rand_map
                    .get(label.as_str())
                    .ok_or(QuerySetError::MissingPolynomial { label: label.to_string() })?;

                if num_polys == 1 && cur_poly.degree_bound().is_some() {
                    assert!(coeff.is_one(), "Coefficient must be one for degree-bounded equations");
                    degree_bound = cur_poly.degree_bound();
                } else if cur_poly.degree_bound().is_some() {
                    eprintln!("Degree bound when number of equations is non-zero");
                    return Err(Self::Error::EquationHasDegreeBounds(equation_label));
                }

                // Some(_) > None, always.
                hiding_bound = std::cmp::max(hiding_bound, cur_poly.hiding_bound());
                poly += (*coeff, cur_poly.polynomial());
                randomness += (*coeff, cur_rand);

                if degree_bound.is_none() {
                    assert!(randomness.shifted_rand.is_none());
                }
            }
            let equation_poly = LabeledPolynomial::new_owned(equation_label, poly, degree_bound, hiding_bound);
            assert_eq!(equation_poly.evaluate(equation.evaluation_point), equation.rhs);

            equation_polynomials.push(equation_poly);
            equation_randomness.push(randomness);
            query_set.insert((&equation.label, equation.evaluation_point));
            evaluations.insert((&equation.label, equation.evaluation_point), equation.rhs);
        }

        Self::batch_open(ck, equation_polynomials.iter(), &query_set, opening_challenge, equation_randomness.iter())
    }

    /// Checks that `values` are the true evaluations at `query_set` of the polynomials
    /// committed in `labeled_commitments`.
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
        let label_comm_map =
            commitments
            .into_iter()
            .map(|c| (c.label(), c))
            .collect::<BTreeMap<_, _>>();

        let mut commitments = Vec::new();
        let mut equation_info = Vec::new();
        let mut query_set = QuerySet::new();
        let mut evaluations = Evaluations::new();
        for equation in equations {
            let equation_label = equation.label.clone();
            let num_polys = equation.lhs.len();

            let mut degree_bound = None;
            let mut coeffs_and_comms = Vec::new();

            for (coeff, label) in &equation.lhs {
                let &cur_comm = label_comm_map
                    .get(label.as_str())
                    .ok_or(QuerySetError::MissingPolynomial { label: label.to_string() })?;

                if num_polys == 1 && cur_comm.degree_bound().is_some() {
                    assert!(coeff.is_one(), "Coefficient must be one for degree-bounded equations");
                    degree_bound = cur_comm.degree_bound();
                } else if cur_comm.degree_bound().is_some() {
                    return Err(Self::Error::EquationHasDegreeBounds(equation_label));
                }
                coeffs_and_comms.push((*coeff, cur_comm.commitment()));
            }
            commitments.push(Self::combine_commitments(coeffs_and_comms));
            equation_info.push((equation_label, degree_bound));

            query_set.insert((&equation.label, equation.evaluation_point));
            evaluations.insert((&equation.label, equation.evaluation_point), equation.rhs);
        }
        let (mut comms, mut shifted_comms): (Vec<_>, Vec<_>) = commitments.into_iter().unzip();
        E::G1Projective::batch_normalization(&mut comms);
        E::G1Projective::batch_normalization(&mut shifted_comms);
        let equation_commitments = equation_info
            .into_iter()
            .zip(comms)
            .zip(shifted_comms)
            .map(|(((label, degree_bound), c), s_c)| {
                let shifted_comm = if !s_c.is_zero() {
                    Some(kzg10::Commitment(s_c.into()))
                } else {
                    None
                };
                let commitment = Commitment {
                    comm: kzg10::Commitment(c.into()),
                    shifted_comm: shifted_comm,
                };
                LabeledCommitment::new(label, commitment, degree_bound)
            })
        .collect::<Vec<_>>();

        Self::batch_check(vk, equation_commitments.iter(), &query_set, &evaluations, proof, opening_challenge, rng)
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
        single_poly_test::<_, PC_Bls12_381>().expect("test failed for bls12-381");
        single_poly_test::<_, PC_MNT6>().expect("test failed for MNT6");
        single_poly_test::<_, PC_SW6>().expect("test failed for SW6");
    }


    #[test]
    fn quadratic_poly_degree_bound_multiple_queries_test() {
        use crate::tests::*;
        quadratic_poly_degree_bound_multiple_queries_test::<_, PC_Bls12_377>().expect("test failed for bls12-377");
        quadratic_poly_degree_bound_multiple_queries_test::<_, PC_Bls12_381>().expect("test failed for bls12-381");
        quadratic_poly_degree_bound_multiple_queries_test::<_, PC_MNT6>().expect("test failed for MNT6");
        quadratic_poly_degree_bound_multiple_queries_test::<_, PC_SW6>().expect("test failed for SW6");
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
        single_poly_degree_bound_multiple_queries_test::<_, PC_SW6>()
            .expect("test failed for SW6");
    }

    #[test]
    fn two_polys_degree_bound_single_query_test() {
        use crate::tests::*;
        two_polys_degree_bound_single_query_test::<_, PC_Bls12_377>()
            .expect("test failed for bls12-377");
        two_polys_degree_bound_single_query_test::<_, PC_Bls12_381>()
            .expect("test failed for bls12-381");
        two_polys_degree_bound_single_query_test::<_, PC_MNT6>()
            .expect("test failed for MNT6");
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
