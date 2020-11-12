use crate::{kzg10, PCCommitterKey};
use crate::{BTreeMap, BTreeSet, String, ToString, Vec};
use crate::{BatchLCProof, Error, Evaluations, QuerySet};
use crate::{LabeledCommitment, LabeledPolynomial, LinearCombination};
use crate::{PCRandomness, PCUniversalParams, PolynomialCommitment, UVPolynomial};

use ark_ec::{AffineCurve, PairingEngine, ProjectiveCurve};
use ark_ff::{Field, One, Zero};
use ark_std::{convert::TryInto, marker::PhantomData, ops::Div, vec};
use rand_core::RngCore;

mod data_structures;
pub use data_structures::*;

mod constraints;
pub use constraints::*;

/// Polynomial commitment based on [[KZG10]][kzg], with degree enforcement, batching,
/// and (optional) hiding property taken from [[CHMMVW20, “Marlin”]][marlin].
///
/// Degree bound enforcement requires that (at least one of) the points at
/// which a committed polynomial is evaluated are from a distribution that is
/// random conditioned on the polynomial. This is because degree bound
/// enforcement relies on checking a polynomial identity at this point.
/// More formally, the points must be sampled from an admissible query sampler,
/// as detailed in [[CHMMVW20]][marlin].
///
/// [kzg]: http://cacr.uwaterloo.ca/techreports/2010/cacr2010-10.pdf
/// [marlin]: https://eprint.iacr.org/2019/104
pub struct MarlinKZG10<E: PairingEngine, P: UVPolynomial<E::Fr>> {
    _engine: PhantomData<E>,
    _poly: PhantomData<P>,
}

pub(crate) fn shift_polynomial<E: PairingEngine, P: UVPolynomial<E::Fr>>(
    ck: &CommitterKey<E>,
    p: &P,
    degree_bound: usize,
) -> P {
    if p.is_zero() {
        P::zero()
    } else {
        let enforced_degree_bounds = ck
            .enforced_degree_bounds
            .as_ref()
            .expect("Polynomial requires degree bounds, but `ck` does not support any");
        let largest_enforced_degree_bound = enforced_degree_bounds.last().unwrap();

        let mut shifted_polynomial_coeffs =
            vec![E::Fr::zero(); largest_enforced_degree_bound - degree_bound];
        shifted_polynomial_coeffs.extend_from_slice(&p.coeffs());
        P::from_coefficients_vec(shifted_polynomial_coeffs)
    }
}

impl<E: PairingEngine, P: UVPolynomial<E::Fr>> MarlinKZG10<E, P> {
    /// MSM for `commitments` and `coeffs`
    fn combine_commitments<'a>(
        coeffs_and_comms: impl IntoIterator<Item = (E::Fr, &'a Commitment<E>)>,
    ) -> (E::G1Projective, Option<E::G1Projective>) {
        let mut combined_comm = E::G1Projective::zero();
        let mut combined_shifted_comm = None;
        for (coeff, comm) in coeffs_and_comms {
            if coeff.is_one() {
                combined_comm.add_assign_mixed(&comm.comm.0);
            } else {
                combined_comm += &comm.comm.0.mul(coeff);
            }

            if let Some(shifted_comm) = &comm.shifted_comm {
                let cur = shifted_comm.0.mul(coeff);
                combined_shifted_comm = Some(combined_shifted_comm.map_or(cur, |c| c + cur));
            }
        }
        (combined_comm, combined_shifted_comm)
    }

    fn normalize_commitments(
        commitments: Vec<(E::G1Projective, Option<E::G1Projective>)>,
    ) -> Vec<Commitment<E>> {
        let mut comms = Vec::with_capacity(commitments.len());
        let mut s_comms = Vec::with_capacity(commitments.len());
        let mut s_flags = Vec::with_capacity(commitments.len());
        for (comm, s_comm) in commitments {
            comms.push(comm);
            if let Some(c) = s_comm {
                s_comms.push(c);
                s_flags.push(true);
            } else {
                s_comms.push(E::G1Projective::zero());
                s_flags.push(false);
            }
        }
        let comms = E::G1Projective::batch_normalization_into_affine(&comms);
        let s_comms = E::G1Projective::batch_normalization_into_affine(&s_comms);
        comms
            .into_iter()
            .zip(s_comms)
            .zip(s_flags)
            .map(|((c, s_c), flag)| {
                let shifted_comm = if flag {
                    Some(kzg10::Commitment(s_c))
                } else {
                    None
                };
                Commitment {
                    comm: kzg10::Commitment(c),
                    shifted_comm,
                }
            })
            .collect()
    }

    /// Accumulate `commitments` and `values` according to `opening_challenge`.
    fn accumulate_commitments_and_values_individual_opening_challenges<'a>(
        vk: &VerifierKey<E>,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Commitment<E>>>,
        values: impl IntoIterator<Item = E::Fr>,
        opening_challenges: &dyn Fn(u64) -> E::Fr,
    ) -> Result<(E::G1Projective, E::Fr), Error> {
        let acc_time = start_timer!(|| "Accumulating commitments and values");
        let mut combined_comm = E::G1Projective::zero();
        let mut combined_value = E::Fr::zero();
        let mut opening_challenge_counter = 0;
        for (labeled_commitment, value) in commitments.into_iter().zip(values) {
            let degree_bound = labeled_commitment.degree_bound();
            let commitment = labeled_commitment.commitment();
            assert_eq!(degree_bound.is_some(), commitment.shifted_comm.is_some());

            let challenge_i = opening_challenges(opening_challenge_counter);
            opening_challenge_counter += 1;

            combined_comm += &commitment.comm.0.mul(challenge_i);
            combined_value += &(value * &challenge_i);

            if let Some(degree_bound) = degree_bound {
                let challenge_i_1 = opening_challenges(opening_challenge_counter);
                opening_challenge_counter += 1;

                let shifted_comm = commitment
                    .shifted_comm
                    .as_ref()
                    .unwrap()
                    .0
                    .into_projective();

                let shift_power = vk
                    .get_shift_power(degree_bound)
                    .ok_or(Error::UnsupportedDegreeBound(degree_bound))?;

                let mut adjusted_comm = shifted_comm - &shift_power.mul(value);

                adjusted_comm *= challenge_i_1;
                combined_comm += &adjusted_comm;
            }
        }

        end_timer!(acc_time);
        Ok((combined_comm, combined_value))
    }
}

impl<E, P> PolynomialCommitment<E::Fr, P> for MarlinKZG10<E, P>
where
    E: PairingEngine,
    P: UVPolynomial<E::Fr, Point = E::Fr>,
    for<'a, 'b> &'a P: Div<&'b P, Output = P>,
{
    type UniversalParams = UniversalParams<E>;
    type CommitterKey = CommitterKey<E>;
    type VerifierKey = VerifierKey<E>;
    type PreparedVerifierKey = PreparedVerifierKey<E>;
    type Commitment = Commitment<E>;
    type PreparedCommitment = PreparedCommitment<E>;
    type Randomness = Randomness<E::Fr, P>;
    type Proof = kzg10::Proof<E>;
    type BatchProof = Vec<Self::Proof>;
    type Error = Error;

    /// Constructs public parameters when given as input the maximum degree `max_degree`
    /// for the polynomial commitment scheme.
    fn setup<R: RngCore>(
        max_degree: usize,
        _num_vars: Option<usize>,
        rng: &mut R,
    ) -> Result<Self::UniversalParams, Self::Error> {
        kzg10::KZG10::setup(max_degree, false, rng).map_err(Into::into)
    }

    fn trim(
        pp: &Self::UniversalParams,
        supported_degree: usize,
        supported_hiding_bound: usize,
        enforced_degree_bounds: Option<&[usize]>,
    ) -> Result<(Self::CommitterKey, Self::VerifierKey), Self::Error> {
        let max_degree = pp.max_degree();
        if supported_degree > max_degree {
            return Err(Error::TrimmingDegreeTooLarge);
        }

        // Construct the KZG10 committer key for committing to unshifted polynomials.
        let ck_time = start_timer!(|| format!(
            "Constructing `powers` of size {} for unshifted polys",
            supported_degree
        ));
        let powers = pp.powers_of_g[..=supported_degree].to_vec();
        // We want to support making up to `supported_hiding_bound` queries to committed
        // polynomials.
        let powers_of_gamma_g = (0..=supported_hiding_bound + 1)
            .map(|i| pp.powers_of_gamma_g[&i])
            .collect::<Vec<_>>();

        end_timer!(ck_time);

        // Construct the core KZG10 verifier key.
        let vk = kzg10::VerifierKey {
            g: pp.powers_of_g[0],
            gamma_g: pp.powers_of_gamma_g[&0],
            h: pp.h,
            beta_h: pp.beta_h,
            prepared_h: pp.prepared_h.clone(),
            prepared_beta_h: pp.prepared_beta_h.clone(),
        };

        let enforced_degree_bounds = enforced_degree_bounds.map(|v| {
            let mut v = v.to_vec();
            v.sort_unstable();
            v.dedup();
            v
        });

        // Check whether we have some degree bounds to enforce
        let (shifted_powers, degree_bounds_and_shift_powers) =
            if let Some(enforced_degree_bounds) = enforced_degree_bounds.as_ref() {
                if enforced_degree_bounds.is_empty() {
                    (None, None)
                } else {
                    let mut sorted_enforced_degree_bounds = enforced_degree_bounds.clone();
                    sorted_enforced_degree_bounds.sort_unstable();

                    let lowest_shifted_power = max_degree
                        - sorted_enforced_degree_bounds
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
                        .map(|d| (*d, pp.powers_of_g[max_degree - *d]))
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
            enforced_degree_bounds,
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
    #[allow(clippy::type_complexity)]
    fn commit<'a>(
        ck: &Self::CommitterKey,
        polynomials: impl IntoIterator<Item = &'a LabeledPolynomial<E::Fr, P>>,
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
        let rng = &mut crate::optional_rng::OptionalRng(rng);
        let commit_time = start_timer!(|| "Committing to polynomials");

        let mut commitments = Vec::new();
        let mut randomness = Vec::new();

        for p in polynomials {
            let label = p.label();
            let degree_bound = p.degree_bound();
            let hiding_bound = p.hiding_bound();
            let polynomial: &P = p.polynomial();

            let enforced_degree_bounds: Option<&[usize]> = ck.enforced_degree_bounds.as_deref();
            kzg10::KZG10::<E, P>::check_degrees_and_bounds(
                ck.supported_degree(),
                ck.max_degree,
                enforced_degree_bounds,
                &p,
            )?;

            let commit_time = start_timer!(|| format!(
                "Polynomial {} of degree {}, degree bound {:?}, and hiding bound {:?}",
                label,
                polynomial.degree(),
                degree_bound,
                hiding_bound,
            ));

            let (comm, rand) =
                kzg10::KZG10::commit(&ck.powers(), polynomial, hiding_bound, Some(rng))?;
            let (shifted_comm, shifted_rand) = if let Some(degree_bound) = degree_bound {
                let shifted_powers = ck
                    .shifted_powers(degree_bound)
                    .ok_or(Error::UnsupportedDegreeBound(degree_bound))?;
                let (shifted_comm, shifted_rand) =
                    kzg10::KZG10::commit(&shifted_powers, &polynomial, hiding_bound, Some(rng))?;
                (Some(shifted_comm), Some(shifted_rand))
            } else {
                (None, None)
            };

            let comm = Commitment { comm, shifted_comm };
            let rand = Randomness { rand, shifted_rand };
            commitments.push(LabeledCommitment::new(
                label.to_string(),
                comm,
                degree_bound,
            ));
            randomness.push(rand);
            end_timer!(commit_time);
        }
        end_timer!(commit_time);
        Ok((commitments, randomness))
    }

    /// On input a polynomial `p` and a point `point`, outputs a proof for the same.
    fn open<'a>(
        ck: &Self::CommitterKey,
        labeled_polynomials: impl IntoIterator<Item = &'a LabeledPolynomial<E::Fr, P>>,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Self::Commitment>>,
        point: &'a P::Point,
        opening_challenge: E::Fr,
        rands: impl IntoIterator<Item = &'a Self::Randomness>,
        rng: Option<&mut dyn RngCore>,
    ) -> Result<Self::Proof, Self::Error>
    where
        Randomness<E::Fr, P>: 'a,
        Commitment<E>: 'a,
    {
        let opening_challenges = |j| opening_challenge.pow([j]);
        Self::open_individual_opening_challenges(
            ck,
            labeled_polynomials,
            commitments,
            &point,
            &opening_challenges,
            rands,
            rng,
        )
    }

    /// Verifies that `value` is the evaluation at `x` of the polynomial
    /// committed inside `comm`.
    fn check<'a>(
        vk: &Self::VerifierKey,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Self::Commitment>>,
        point: &'a P::Point,
        values: impl IntoIterator<Item = E::Fr>,
        proof: &Self::Proof,
        opening_challenge: E::Fr,
        rng: Option<&mut dyn RngCore>,
    ) -> Result<bool, Self::Error>
    where
        Self::Commitment: 'a,
    {
        let opening_challenges = |j| opening_challenge.pow([j]);
        Self::check_individual_opening_challenges(
            vk,
            commitments,
            &point,
            values,
            proof,
            &opening_challenges,
            rng,
        )
    }

    fn batch_check<'a, R: RngCore>(
        vk: &Self::VerifierKey,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Self::Commitment>>,
        query_set: &QuerySet<E::Fr>,
        evaluations: &Evaluations<P::Point, E::Fr>,
        proof: &Self::BatchProof,
        opening_challenge: E::Fr,
        rng: &mut R,
    ) -> Result<bool, Self::Error>
    where
        Self::Commitment: 'a,
    {
        let opening_challenges = |j| opening_challenge.pow([j]);
        Self::batch_check_individual_opening_challenges(
            vk,
            commitments,
            query_set,
            evaluations,
            proof,
            &opening_challenges,
            rng,
        )
    }

    fn open_combinations<'a>(
        ck: &Self::CommitterKey,
        lc_s: impl IntoIterator<Item = &'a LinearCombination<E::Fr>>,
        polynomials: impl IntoIterator<Item = &'a LabeledPolynomial<E::Fr, P>>,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Self::Commitment>>,
        query_set: &QuerySet<E::Fr>,
        opening_challenge: E::Fr,
        rands: impl IntoIterator<Item = &'a Self::Randomness>,
        rng: Option<&mut dyn RngCore>,
    ) -> Result<BatchLCProof<E::Fr, P, Self>, Self::Error>
    where
        P: 'a,
        Self::Randomness: 'a,
        Self::Commitment: 'a,
    {
        let opening_challenges = |j| opening_challenge.pow([j]);
        Self::open_combinations_individual_opening_challenges(
            ck,
            lc_s,
            polynomials,
            commitments,
            query_set,
            &opening_challenges,
            rands,
            rng,
        )
    }

    /// Checks that `values` are the true evaluations at `query_set` of the polynomials
    /// committed in `labeled_commitments`.
    fn check_combinations<'a, R: RngCore>(
        vk: &Self::VerifierKey,
        lc_s: impl IntoIterator<Item = &'a LinearCombination<E::Fr>>,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Self::Commitment>>,
        query_set: &QuerySet<E::Fr>,
        evaluations: &Evaluations<P::Point, E::Fr>,
        proof: &BatchLCProof<E::Fr, P, Self>,
        opening_challenge: E::Fr,
        rng: &mut R,
    ) -> Result<bool, Self::Error>
    where
        Self::Commitment: 'a,
    {
        let opening_challenges = |j| opening_challenge.pow([j]);
        Self::check_combinations_individual_opening_challenges(
            vk,
            lc_s,
            commitments,
            query_set,
            evaluations,
            proof,
            &opening_challenges,
            rng,
        )
    }

    // On input a list of labeled polynomials and a query set, `open` outputs a proof of evaluation
    /// of the polynomials at the points in the query set.
    fn batch_open<'a>(
        ck: &Self::CommitterKey,
        labeled_polynomials: impl IntoIterator<Item = &'a LabeledPolynomial<E::Fr, P>>,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Self::Commitment>>,
        query_set: &QuerySet<E::Fr>,
        opening_challenge: E::Fr,
        rands: impl IntoIterator<Item = &'a Self::Randomness>,
        rng: Option<&mut dyn RngCore>,
    ) -> Result<Self::BatchProof, Self::Error>
    where
        P: 'a,
        Self::Randomness: 'a,
        Self::Commitment: 'a,
    {
        let opening_challenges = |j| opening_challenge.pow([j]);
        Self::batch_open_individual_opening_challenges(
            ck,
            labeled_polynomials,
            commitments,
            query_set,
            &opening_challenges,
            rands,
            rng,
        )
    }

    /// On input a polynomial `p` and a point `point`, outputs a proof for the same.
    fn open_individual_opening_challenges<'a>(
        ck: &CommitterKey<E>,
        labeled_polynomials: impl IntoIterator<Item = &'a LabeledPolynomial<E::Fr, P>>,
        _commitments: impl IntoIterator<Item = &'a LabeledCommitment<Commitment<E>>>,
        point: &'a P::Point,
        opening_challenges: &dyn Fn(u64) -> E::Fr,
        rands: impl IntoIterator<Item = &'a Randomness<E::Fr, P>>,
        _rng: Option<&mut dyn RngCore>,
    ) -> Result<kzg10::Proof<E>, Error>
    where
        P: 'a,
        Randomness<E::Fr, P>: 'a,
        Commitment<E>: 'a,
    {
        let mut p = P::zero();
        let mut r = kzg10::Randomness::empty();
        let mut shifted_w = P::zero();
        let mut shifted_r = kzg10::Randomness::empty();
        let mut shifted_r_witness = P::zero();

        let mut enforce_degree_bound = false;
        let mut opening_challenge_counter = 0;
        for (polynomial, rand) in labeled_polynomials.into_iter().zip(rands) {
            let degree_bound = polynomial.degree_bound();
            assert_eq!(degree_bound.is_some(), rand.shifted_rand.is_some());

            let enforced_degree_bounds: Option<&[usize]> = ck.enforced_degree_bounds.as_deref();
            kzg10::KZG10::<E, P>::check_degrees_and_bounds(
                ck.supported_degree(),
                ck.max_degree,
                enforced_degree_bounds,
                &polynomial,
            )?;

            // compute challenge^j and challenge^{j+1}.
            let challenge_j = opening_challenges(opening_challenge_counter);
            opening_challenge_counter += 1;

            assert_eq!(degree_bound.is_some(), rand.shifted_rand.is_some());

            p += (challenge_j, polynomial.polynomial());
            r += (challenge_j, &rand.rand);

            if let Some(degree_bound) = degree_bound {
                enforce_degree_bound = true;
                let shifted_rand = rand.shifted_rand.as_ref().unwrap();
                let (witness, shifted_rand_witness) =
                    kzg10::KZG10::<E, P>::compute_witness_polynomial(
                        polynomial.polynomial(),
                        *point,
                        &shifted_rand,
                    )?;
                let challenge_j_1 = opening_challenges(opening_challenge_counter);
                opening_challenge_counter += 1;

                let shifted_witness = shift_polynomial(ck, &witness, degree_bound);

                shifted_w += (challenge_j_1, &shifted_witness);
                shifted_r += (challenge_j_1, shifted_rand);
                if let Some(shifted_rand_witness) = shifted_rand_witness {
                    shifted_r_witness += (challenge_j_1, &shifted_rand_witness);
                }
            }
        }
        let proof_time = start_timer!(|| "Creating proof for unshifted polynomials");
        let proof = kzg10::KZG10::open(&ck.powers(), &p, *point, &r)?;
        let mut w = proof.w.into_projective();
        let mut random_v = proof.random_v;
        end_timer!(proof_time);

        if enforce_degree_bound {
            let proof_time = start_timer!(|| "Creating proof for shifted polynomials");
            let shifted_proof = kzg10::KZG10::open_with_witness_polynomial(
                &ck.shifted_powers(None).unwrap(),
                *point,
                &shifted_r,
                &shifted_w,
                Some(&shifted_r_witness),
            )?;
            end_timer!(proof_time);

            w += &shifted_proof.w.into_projective();
            if let Some(shifted_random_v) = shifted_proof.random_v {
                random_v = random_v.map(|v| v + &shifted_random_v);
            }
        }

        Ok(kzg10::Proof {
            w: w.into_affine(),
            random_v,
        })
    }

    /// Verifies that `value` is the evaluation at `x` of the polynomial
    /// committed inside `comm`.
    fn check_individual_opening_challenges<'a>(
        vk: &VerifierKey<E>,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Commitment<E>>>,
        point: &'a P::Point,
        values: impl IntoIterator<Item = E::Fr>,
        proof: &kzg10::Proof<E>,
        opening_challenges: &dyn Fn(u64) -> E::Fr,
        _rng: Option<&mut dyn RngCore>,
    ) -> Result<bool, Error>
    where
        Commitment<E>: 'a,
    {
        let check_time = start_timer!(|| "Checking evaluations");
        let (combined_comm, combined_value) =
            Self::accumulate_commitments_and_values_individual_opening_challenges(
                vk,
                commitments,
                values,
                opening_challenges,
            )?;
        let combined_comm = kzg10::Commitment(combined_comm.into());
        let result = kzg10::KZG10::check(&vk.vk, &combined_comm, *point, combined_value, proof)?;
        end_timer!(check_time);
        Ok(result)
    }

    fn batch_check_individual_opening_challenges<'a, R: RngCore>(
        vk: &VerifierKey<E>,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Commitment<E>>>,
        query_set: &QuerySet<P::Point>,
        evaluations: &Evaluations<E::Fr, P::Point>,
        proof: &Vec<kzg10::Proof<E>>,
        opening_challenges: &dyn Fn(u64) -> E::Fr,
        rng: &mut R,
    ) -> Result<bool, Error>
    where
        Commitment<E>: 'a,
    {
        let commitments: BTreeMap<_, _> = commitments.into_iter().map(|c| (c.label(), c)).collect();
        let mut query_to_labels_map = BTreeMap::new();

        for (label, (point_label, point)) in query_set.iter() {
            let labels = query_to_labels_map
                .entry(point_label)
                .or_insert((point, BTreeSet::new()));
            labels.1.insert(label);
        }
        assert_eq!(proof.len(), query_to_labels_map.len());

        let mut combined_comms = Vec::new();
        let mut combined_queries = Vec::new();
        let mut combined_evals = Vec::new();
        for (_, (point, labels)) in query_to_labels_map.into_iter() {
            let lc_time =
                start_timer!(|| format!("Randomly combining {} commitments", labels.len()));
            let mut comms_to_combine: Vec<&'_ LabeledCommitment<_>> = Vec::new();
            let mut values_to_combine = Vec::new();
            for label in labels.into_iter() {
                let commitment = commitments.get(label).ok_or(Error::MissingPolynomial {
                    label: label.to_string(),
                })?;
                let degree_bound = commitment.degree_bound();
                assert_eq!(
                    degree_bound.is_some(),
                    commitment.commitment().shifted_comm.is_some()
                );

                let v_i =
                    evaluations
                        .get(&(label.clone(), *point))
                        .ok_or(Error::MissingEvaluation {
                            label: label.to_string(),
                        })?;

                comms_to_combine.push(commitment);
                values_to_combine.push(*v_i);
            }

            let (c, v) = Self::accumulate_commitments_and_values_individual_opening_challenges(
                vk,
                comms_to_combine,
                values_to_combine,
                opening_challenges,
            )?;
            end_timer!(lc_time);

            combined_comms.push(c);
            combined_queries.push(*point);
            combined_evals.push(v);
        }
        let norm_time = start_timer!(|| "Normalizaing combined commitments");
        E::G1Projective::batch_normalization(&mut combined_comms);
        let combined_comms = combined_comms
            .into_iter()
            .map(|c| kzg10::Commitment(c.into()))
            .collect::<Vec<_>>();
        end_timer!(norm_time);
        let proof_time = start_timer!(|| "Checking KZG10::Proof");
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

    fn open_combinations_individual_opening_challenges<'a>(
        ck: &CommitterKey<E>,
        lc_s: impl IntoIterator<Item = &'a LinearCombination<E::Fr>>,
        polynomials: impl IntoIterator<Item = &'a LabeledPolynomial<E::Fr, P>>,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Commitment<E>>>,
        query_set: &QuerySet<P::Point>,
        opening_challenges: &dyn Fn(u64) -> E::Fr,
        rands: impl IntoIterator<Item = &'a Randomness<E::Fr, P>>,
        rng: Option<&mut dyn RngCore>,
    ) -> Result<BatchLCProof<E::Fr, P, Self>, Error>
    where
        P: 'a,
        Randomness<E::Fr, P>: 'a,
        Commitment<E>: 'a,
    {
        let label_map = polynomials
            .into_iter()
            .zip(rands)
            .zip(commitments)
            .map(|((p, r), c)| (p.label(), (p, r, c)))
            .collect::<BTreeMap<_, _>>();

        let mut lc_polynomials = Vec::new();
        let mut lc_randomness = Vec::new();
        let mut lc_commitments = Vec::new();
        let mut lc_info = Vec::new();

        for lc in lc_s {
            let lc_label = lc.label().clone();
            let mut poly = P::zero();
            let mut degree_bound = None;
            let mut hiding_bound = None;

            let mut randomness = Randomness::<E::Fr, P>::empty();
            assert!(randomness.shifted_rand.is_none());

            let mut coeffs_and_comms = Vec::new();

            let num_polys = lc.len();
            for (coeff, label) in lc.iter().filter(|(_, l)| !l.is_one()) {
                let label: &String = label.try_into().expect("cannot be one!");
                let &(cur_poly, cur_rand, cur_comm) =
                    label_map.get(label).ok_or(Error::MissingPolynomial {
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
                    return Err(Error::EquationHasDegreeBounds(lc_label));
                }

                // Some(_) > None, always.
                hiding_bound = core::cmp::max(hiding_bound, cur_poly.hiding_bound());
                poly += (*coeff, cur_poly.polynomial());
                randomness += (*coeff, cur_rand);
                coeffs_and_comms.push((*coeff, cur_comm.commitment()));

                if degree_bound.is_none() {
                    assert!(randomness.shifted_rand.is_none());
                }
            }

            let lc_poly =
                LabeledPolynomial::new(lc_label.clone(), poly, degree_bound, hiding_bound);
            lc_polynomials.push(lc_poly);
            lc_randomness.push(randomness);
            lc_commitments.push(Self::combine_commitments(coeffs_and_comms));
            lc_info.push((lc_label, degree_bound));
        }

        let comms = Self::normalize_commitments(lc_commitments);
        let lc_commitments = lc_info
            .into_iter()
            .zip(comms)
            .map(|((label, d), c)| LabeledCommitment::new(label, c, d))
            .collect::<Vec<_>>();

        let proof = Self::batch_open_individual_opening_challenges(
            ck,
            lc_polynomials.iter(),
            lc_commitments.iter(),
            &query_set,
            opening_challenges,
            lc_randomness.iter(),
            rng,
        )?;

        Ok(BatchLCProof { proof, evals: None })
    }

    /// Checks that `values` are the true evaluations at `query_set` of the polynomials
    /// committed in `labeled_commitments`.
    fn check_combinations_individual_opening_challenges<'a, R: RngCore>(
        vk: &VerifierKey<E>,
        lc_s: impl IntoIterator<Item = &'a LinearCombination<E::Fr>>,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Commitment<E>>>,
        query_set: &QuerySet<P::Point>,
        evaluations: &Evaluations<E::Fr, P::Point>,
        proof: &BatchLCProof<E::Fr, P, Self>,
        opening_challenges: &dyn Fn(u64) -> E::Fr,
        rng: &mut R,
    ) -> Result<bool, Error>
    where
        Commitment<E>: 'a,
    {
        let BatchLCProof { proof, .. } = proof;
        let label_comm_map = commitments
            .into_iter()
            .map(|c| (c.label(), c))
            .collect::<BTreeMap<_, _>>();

        let mut lc_commitments = Vec::new();
        let mut lc_info = Vec::new();
        let mut evaluations = evaluations.clone();

        let lc_processing_time = start_timer!(|| "Combining commitments");
        for lc in lc_s {
            let lc_label = lc.label().clone();
            let num_polys = lc.len();

            let mut degree_bound = None;
            let mut coeffs_and_comms = Vec::new();

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
                        return Err(Error::EquationHasDegreeBounds(lc_label));
                    }
                    coeffs_and_comms.push((*coeff, cur_comm.commitment()));
                }
            }
            let lc_time =
                start_timer!(|| format!("Combining {} commitments for {}", num_polys, lc_label));
            lc_commitments.push(Self::combine_commitments(coeffs_and_comms));
            end_timer!(lc_time);
            lc_info.push((lc_label, degree_bound));
        }
        end_timer!(lc_processing_time);
        let combined_comms_norm_time = start_timer!(|| "Normalizing commitments");
        let comms = Self::normalize_commitments(lc_commitments);
        let lc_commitments = lc_info
            .into_iter()
            .zip(comms)
            .map(|((label, d), c)| LabeledCommitment::new(label, c, d))
            .collect::<Vec<_>>();
        end_timer!(combined_comms_norm_time);

        Self::batch_check_individual_opening_challenges(
            vk,
            &lc_commitments,
            &query_set,
            &evaluations,
            proof,
            opening_challenges,
            rng,
        )
    }

    /// On input a list of labeled polynomials and a query set, `open` outputs a proof of evaluation
    /// of the polynomials at the points in the query set.
    fn batch_open_individual_opening_challenges<'a>(
        ck: &CommitterKey<E>,
        labeled_polynomials: impl IntoIterator<Item = &'a LabeledPolynomial<E::Fr, P>>,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Commitment<E>>>,
        query_set: &QuerySet<E::Fr>,
        opening_challenges: &dyn Fn(u64) -> E::Fr,
        rands: impl IntoIterator<Item = &'a Randomness<E::Fr, P>>,
        rng: Option<&mut dyn RngCore>,
    ) -> Result<Vec<kzg10::Proof<E>>, Error>
    where
        P: 'a,
        Randomness<E::Fr, P>: 'a,
        Commitment<E>: 'a,
    {
        let rng = &mut crate::optional_rng::OptionalRng(rng);
        let poly_rand_comm: BTreeMap<_, _> = labeled_polynomials
            .into_iter()
            .zip(rands)
            .zip(commitments.into_iter())
            .map(|((poly, r), comm)| (poly.label(), (poly, r, comm)))
            .collect();

        let open_time = start_timer!(|| format!(
            "Opening {} polynomials at query set of size {}",
            poly_rand_comm.len(),
            query_set.len(),
        ));

        let mut query_to_labels_map = BTreeMap::new();

        for (label, (point_label, point)) in query_set.iter() {
            let labels = query_to_labels_map
                .entry(point_label)
                .or_insert((point, BTreeSet::new()));
            labels.1.insert(label);
        }

        let mut proofs = Vec::new();
        for (_point_label, (point, labels)) in query_to_labels_map.into_iter() {
            let mut query_polys: Vec<&'a LabeledPolynomial<_, _>> = Vec::new();
            let mut query_rands: Vec<&'a Randomness<E::Fr, P>> = Vec::new();
            let mut query_comms: Vec<&'a LabeledCommitment<Commitment<E>>> = Vec::new();

            for label in labels {
                let (polynomial, rand, comm) =
                    poly_rand_comm.get(&label).ok_or(Error::MissingPolynomial {
                        label: label.to_string(),
                    })?;

                query_polys.push(polynomial);
                query_rands.push(rand);
                query_comms.push(comm);
            }

            let proof_time = start_timer!(|| "Creating proof");
            let proof = Self::open_individual_opening_challenges(
                ck,
                query_polys,
                query_comms,
                point,
                opening_challenges,
                query_rands,
                Some(rng),
            )?;

            end_timer!(proof_time);

            proofs.push(proof);
        }
        end_timer!(open_time);

        Ok(proofs)
    }
}

#[cfg(test)]
mod tests {
    #![allow(non_camel_case_types)]
    use super::MarlinKZG10;
    use ark_bls12_377::Bls12_377;
    use ark_bls12_381::Bls12_381;
    use ark_ec::PairingEngine;
    use ark_ff::UniformRand;
    use ark_poly::{univariate::DensePolynomial as DensePoly, UVPolynomial};

    type UniPoly_381 = DensePoly<<Bls12_381 as PairingEngine>::Fr>;
    type UniPoly_377 = DensePoly<<Bls12_377 as PairingEngine>::Fr>;

    type PC<E, P> = MarlinKZG10<E, P>;
    type PC_Bls12_381 = PC<Bls12_381, UniPoly_381>;
    type PC_Bls12_377 = PC<Bls12_377, UniPoly_377>;

    fn rand_poly<E: PairingEngine>(
        degree: usize,
        _: Option<usize>,
        rng: &mut rand::prelude::StdRng,
    ) -> DensePoly<E::Fr> {
        DensePoly::<E::Fr>::rand(degree, rng)
    }

    fn rand_point<E: PairingEngine>(_: Option<usize>, rng: &mut rand::prelude::StdRng) -> E::Fr {
        E::Fr::rand(rng)
    }

    #[test]
    fn single_poly_test() {
        use crate::tests::*;
        single_poly_test::<_, _, PC_Bls12_377>(
            None,
            rand_poly::<Bls12_377>,
            rand_point::<Bls12_377>,
        )
        .expect("test failed for bls12-377");
        single_poly_test::<_, _, PC_Bls12_381>(
            None,
            rand_poly::<Bls12_381>,
            rand_point::<Bls12_381>,
        )
        .expect("test failed for bls12-381");
    }

    #[test]
    fn quadratic_poly_degree_bound_multiple_queries_test() {
        use crate::tests::*;
        quadratic_poly_degree_bound_multiple_queries_test::<_, _, PC_Bls12_377>(
            rand_poly::<Bls12_377>,
            rand_point::<Bls12_377>,
        )
        .expect("test failed for bls12-377");
        quadratic_poly_degree_bound_multiple_queries_test::<_, _, PC_Bls12_381>(
            rand_poly::<Bls12_381>,
            rand_point::<Bls12_381>,
        )
        .expect("test failed for bls12-381");
    }

    #[test]
    fn linear_poly_degree_bound_test() {
        use crate::tests::*;
        linear_poly_degree_bound_test::<_, _, PC_Bls12_377>(
            rand_poly::<Bls12_377>,
            rand_point::<Bls12_377>,
        )
        .expect("test failed for bls12-377");
        linear_poly_degree_bound_test::<_, _, PC_Bls12_381>(
            rand_poly::<Bls12_381>,
            rand_point::<Bls12_381>,
        )
        .expect("test failed for bls12-381");
    }

    #[test]
    fn single_poly_degree_bound_test() {
        use crate::tests::*;
        single_poly_degree_bound_test::<_, _, PC_Bls12_377>(
            rand_poly::<Bls12_377>,
            rand_point::<Bls12_377>,
        )
        .expect("test failed for bls12-377");
        single_poly_degree_bound_test::<_, _, PC_Bls12_381>(
            rand_poly::<Bls12_381>,
            rand_point::<Bls12_381>,
        )
        .expect("test failed for bls12-381");
    }

    #[test]
    fn single_poly_degree_bound_multiple_queries_test() {
        use crate::tests::*;
        single_poly_degree_bound_multiple_queries_test::<_, _, PC_Bls12_377>(
            rand_poly::<Bls12_377>,
            rand_point::<Bls12_377>,
        )
        .expect("test failed for bls12-377");
        single_poly_degree_bound_multiple_queries_test::<_, _, PC_Bls12_381>(
            rand_poly::<Bls12_381>,
            rand_point::<Bls12_381>,
        )
        .expect("test failed for bls12-381");
    }

    #[test]
    fn two_polys_degree_bound_single_query_test() {
        use crate::tests::*;
        two_polys_degree_bound_single_query_test::<_, _, PC_Bls12_377>(
            rand_poly::<Bls12_377>,
            rand_point::<Bls12_377>,
        )
        .expect("test failed for bls12-377");
        two_polys_degree_bound_single_query_test::<_, _, PC_Bls12_381>(
            rand_poly::<Bls12_381>,
            rand_point::<Bls12_381>,
        )
        .expect("test failed for bls12-381");
    }

    #[test]
    fn full_end_to_end_test() {
        use crate::tests::*;
        full_end_to_end_test::<_, _, PC_Bls12_377>(
            None,
            rand_poly::<Bls12_377>,
            rand_point::<Bls12_377>,
        )
        .expect("test failed for bls12-377");
        println!("Finished bls12-377");
        full_end_to_end_test::<_, _, PC_Bls12_381>(
            None,
            rand_poly::<Bls12_381>,
            rand_point::<Bls12_381>,
        )
        .expect("test failed for bls12-381");
        println!("Finished bls12-381");
    }

    #[test]
    fn single_equation_test() {
        use crate::tests::*;
        single_equation_test::<_, _, PC_Bls12_377>(
            None,
            rand_poly::<Bls12_377>,
            rand_point::<Bls12_377>,
        )
        .expect("test failed for bls12-377");
        println!("Finished bls12-377");
        single_equation_test::<_, _, PC_Bls12_381>(
            None,
            rand_poly::<Bls12_381>,
            rand_point::<Bls12_381>,
        )
        .expect("test failed for bls12-381");
        println!("Finished bls12-381");
    }

    #[test]
    fn two_equation_test() {
        use crate::tests::*;
        two_equation_test::<_, _, PC_Bls12_377>(
            None,
            rand_poly::<Bls12_377>,
            rand_point::<Bls12_377>,
        )
        .expect("test failed for bls12-377");
        println!("Finished bls12-377");
        two_equation_test::<_, _, PC_Bls12_381>(
            None,
            rand_poly::<Bls12_381>,
            rand_point::<Bls12_381>,
        )
        .expect("test failed for bls12-381");
        println!("Finished bls12-381");
    }

    #[test]
    fn two_equation_degree_bound_test() {
        use crate::tests::*;
        two_equation_degree_bound_test::<_, _, PC_Bls12_377>(
            rand_poly::<Bls12_377>,
            rand_point::<Bls12_377>,
        )
        .expect("test failed for bls12-377");
        println!("Finished bls12-377");
        two_equation_degree_bound_test::<_, _, PC_Bls12_381>(
            rand_poly::<Bls12_381>,
            rand_point::<Bls12_381>,
        )
        .expect("test failed for bls12-381");
        println!("Finished bls12-381");
    }

    #[test]
    fn full_end_to_end_equation_test() {
        use crate::tests::*;
        full_end_to_end_equation_test::<_, _, PC_Bls12_377>(
            None,
            rand_poly::<Bls12_377>,
            rand_point::<Bls12_377>,
        )
        .expect("test failed for bls12-377");
        println!("Finished bls12-377");
        full_end_to_end_equation_test::<_, _, PC_Bls12_381>(
            None,
            rand_poly::<Bls12_381>,
            rand_point::<Bls12_381>,
        )
        .expect("test failed for bls12-381");
        println!("Finished bls12-381");
    }

    #[test]
    #[should_panic]
    fn bad_degree_bound_test() {
        use crate::tests::*;
        bad_degree_bound_test::<_, _, PC_Bls12_377>(
            rand_poly::<Bls12_377>,
            rand_point::<Bls12_377>,
        )
        .expect("test failed for bls12-377");
        println!("Finished bls12-377");
        bad_degree_bound_test::<_, _, PC_Bls12_381>(
            rand_poly::<Bls12_381>,
            rand_point::<Bls12_381>,
        )
        .expect("test failed for bls12-381");
        println!("Finished bls12-381");
    }
}
