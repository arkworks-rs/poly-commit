use crate::{kzg10, Error};
use crate::{BTreeMap, BTreeSet, Debug, RngCore, String, ToString, Vec};
use crate::{BatchLCProof, LabeledPolynomial, LinearCombination};
use crate::{Evaluations, LabeledCommitment, QuerySet};
use crate::{PCRandomness, Polynomial, PolynomialCommitment};
use ark_ec::{AffineCurve, PairingEngine, ProjectiveCurve};
use ark_ff::{One, Zero};
use ark_std::{convert::TryInto, hash::Hash, ops::AddAssign};

/// Polynomial commitment scheme from [[KZG10]][kzg] that enforces
/// strict degree bounds and (optionally) enables hiding commitments by
/// following the approach outlined in [[CHMMVW20, "Marlin"]][marlin].
///
/// [kzg]: http://cacr.uwaterloo.ca/techreports/2010/cacr2010-10.pdf
/// [marlin]: https://eprint.iacr.org/2019/1047
pub mod marlin_pc;

/// Multivariate polynomial commitment based on the construction in
/// [[PST13]][pst] with batching and (optional) hiding property inspired
/// by the univariate scheme in [[CHMMVW20, "Marlin"]][marlin]
///
/// [pst]: https://eprint.iacr.org/2011/587.pdf
/// [marlin]: https://eprint.iacr.org/2019/104
pub mod marlin_pst13_pc;

/// Common functionalities between `marlin_pc` and `marlin_pst13_pc`
struct Marlin<E: PairingEngine> {
    _engine: core::marker::PhantomData<E>,
}

impl<E: PairingEngine> Marlin<E> {
    /// MSM for `commitments` and `coeffs`
    fn combine_commitments<'a>(
        coeffs_and_comms: impl IntoIterator<Item = (E::Fr, &'a marlin_pc::Commitment<E>)>,
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

    /// Normalize a list of commitments
    fn normalize_commitments<'a>(
        commitments: Vec<(E::G1Projective, Option<E::G1Projective>)>,
    ) -> Vec<marlin_pc::Commitment<E>> {
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
        let s_comms = E::G1Projective::batch_normalization_into_affine(&mut s_comms);
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
                marlin_pc::Commitment {
                    comm: kzg10::Commitment(c),
                    shifted_comm,
                }
            })
            .collect()
    }

    /// Accumulate `commitments` and `values` according to `opening_challenge`.
    fn accumulate_commitments_and_values_individual_opening_challenges<'a>(
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<marlin_pc::Commitment<E>>>,
        values: impl IntoIterator<Item = E::Fr>,
        opening_challenges: &dyn Fn(u64) -> E::Fr,
        vk: Option<&marlin_pc::VerifierKey<E>>,
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
                    .unwrap()
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

    /// Combine and normalize a set of commitments
    fn combine_and_normalize<'a, D: Clone + Ord + Sync>(
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<marlin_pc::Commitment<E>>>,
        query_set: &QuerySet<D>,
        evaluations: &Evaluations<D, E::Fr>,
        opening_challenges: &dyn Fn(u64) -> E::Fr,
        vk: Option<&marlin_pc::VerifierKey<E>>,
    ) -> Result<(Vec<kzg10::Commitment<E>>, Vec<D>, Vec<E::Fr>), Error>
    where
        marlin_pc::Commitment<E>: 'a,
    {
        let commitments: BTreeMap<_, _> = commitments.into_iter().map(|c| (c.label(), c)).collect();
        let mut query_to_labels_map = BTreeMap::new();

        for (label, (point_label, point)) in query_set.iter() {
            let labels = query_to_labels_map
                .entry(point_label)
                .or_insert((point, BTreeSet::new()));
            labels.1.insert(label);
        }

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

                let v_i = evaluations.get(&(label.clone(), point.clone())).ok_or(
                    Error::MissingEvaluation {
                        label: label.to_string(),
                    },
                )?;

                comms_to_combine.push(commitment);
                values_to_combine.push(*v_i);
            }

            let (c, v) = Marlin::accumulate_commitments_and_values_individual_opening_challenges(
                comms_to_combine,
                values_to_combine,
                opening_challenges,
                vk,
            )?;
            end_timer!(lc_time);

            combined_comms.push(c);
            combined_queries.push(point.clone());
            combined_evals.push(v);
        }
        let norm_time = start_timer!(|| "Normalizing combined commitments");
        E::G1Projective::batch_normalization(&mut combined_comms);
        let combined_comms = combined_comms
            .into_iter()
            .map(|c| kzg10::Commitment(c.into()))
            .collect::<Vec<_>>();
        end_timer!(norm_time);
        Ok((combined_comms, combined_queries, combined_evals))
    }

    /// On input a list of polynomials, linear combinations of those polynomials,
    /// and a query set, `open_combination` outputs a proof of evaluation of
    /// the combinations at the points in the query set.
    fn open_combinations_individual_opening_challenges<'a, P, D, PC>(
        ck: &PC::CommitterKey,
        lc_s: impl IntoIterator<Item = &'a LinearCombination<E::Fr>>,
        polynomials: impl IntoIterator<Item = &'a LabeledPolynomial<E::Fr, P>>,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<PC::Commitment>>,
        query_set: &QuerySet<D>,
        opening_challenges: &dyn Fn(u64) -> E::Fr,
        rands: impl IntoIterator<Item = &'a PC::Randomness>,
        rng: Option<&mut dyn RngCore>,
    ) -> Result<BatchLCProof<E::Fr, P, PC>, Error>
    where
        P: 'a + Polynomial<E::Fr, Point = D>,
        D: Debug + Clone + Hash + Ord + Sync,
        PC: PolynomialCommitment<
            E::Fr,
            P,
            Commitment = marlin_pc::Commitment<E>,
            PreparedCommitment = marlin_pc::PreparedCommitment<E>,
            Error = Error,
        >,
        PC::Randomness: 'a + AddAssign<(E::Fr, &'a PC::Randomness)>,
        PC::Commitment: 'a,
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

            let mut randomness = PC::Randomness::empty();
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
                    return Err(Error::EquationHasDegreeBounds(lc_label));
                }
                // Some(_) > None, always.
                hiding_bound = core::cmp::max(hiding_bound, cur_poly.hiding_bound());
                poly += (*coeff, cur_poly.polynomial());
                randomness += (*coeff, cur_rand);
                coeffs_and_comms.push((*coeff, cur_comm.commitment()));
            }

            let lc_poly =
                LabeledPolynomial::new(lc_label.clone(), poly, degree_bound, hiding_bound);
            lc_polynomials.push(lc_poly);
            lc_randomness.push(randomness);
            lc_commitments.push(Marlin::combine_commitments(coeffs_and_comms));
            lc_info.push((lc_label, degree_bound));
        }

        let comms = Marlin::normalize_commitments(lc_commitments);
        let lc_commitments = lc_info
            .into_iter()
            .zip(comms)
            .map(|((label, d), c)| LabeledCommitment::new(label, c, d))
            .collect::<Vec<_>>();

        let proof = PC::batch_open_individual_opening_challenges(
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

    fn check_combinations_individual_opening_challenges<'a, R, P, D, PC>(
        vk: &PC::VerifierKey,
        lc_s: impl IntoIterator<Item = &'a LinearCombination<E::Fr>>,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<PC::Commitment>>,
        query_set: &QuerySet<P::Point>,
        evaluations: &Evaluations<P::Point, E::Fr>,
        proof: &BatchLCProof<E::Fr, P, PC>,
        opening_challenges: &dyn Fn(u64) -> E::Fr,
        rng: &mut R,
    ) -> Result<bool, Error>
    where
        R: RngCore,
        P: Polynomial<E::Fr, Point = D>,
        D: Debug + Clone + Hash + Ord + Sync,
        PC: PolynomialCommitment<
            E::Fr,
            P,
            Commitment = marlin_pc::Commitment<E>,
            PreparedCommitment = marlin_pc::PreparedCommitment<E>,
            Error = Error,
        >,
        PC::Commitment: 'a,
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
            lc_commitments.push(Marlin::combine_commitments(coeffs_and_comms));
            end_timer!(lc_time);
            lc_info.push((lc_label, degree_bound));
        }
        end_timer!(lc_processing_time);
        let combined_comms_norm_time = start_timer!(|| "Normalizing commitments");
        let comms = Marlin::normalize_commitments(lc_commitments);
        let lc_commitments = lc_info
            .into_iter()
            .zip(comms)
            .map(|((label, d), c)| LabeledCommitment::new(label, c, d))
            .collect::<Vec<_>>();
        end_timer!(combined_comms_norm_time);

        PC::batch_check_individual_opening_challenges(
            vk,
            &lc_commitments,
            &query_set,
            &evaluations,
            proof,
            opening_challenges,
            rng,
        )
    }
}
