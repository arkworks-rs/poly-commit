use crate::{kzg10, PCCommitterKey, CHALLENGE_SIZE};
use crate::{BTreeMap, BTreeSet, String, ToString, Vec};
use crate::{BatchLCProof, DenseUVPolynomial, Error, Evaluations, QuerySet};
use crate::{LabeledCommitment, LabeledPolynomial, LinearCombination};
use crate::{PCCommitmentState, PCUniversalParams, PolynomialCommitment};
use ark_ec::AffineRepr;
use ark_ec::CurveGroup;

use ark_ec::pairing::Pairing;
use ark_ff::{One, UniformRand, Zero};
use ark_std::rand::RngCore;
use ark_std::{convert::TryInto, marker::PhantomData, ops::Div, ops::Mul, vec};

mod data_structures;
use crate::challenge::ChallengeGenerator;
use ark_crypto_primitives::sponge::CryptographicSponge;
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
pub struct SonicKZG10<E: Pairing, P: DenseUVPolynomial<E::ScalarField>, S: CryptographicSponge> {
    _engine: PhantomData<E>,
    _poly: PhantomData<P>,
    _sponge: PhantomData<S>,
}

impl<E, P, S> SonicKZG10<E, P, S>
where
    E: Pairing,
    P: DenseUVPolynomial<E::ScalarField>,
    S: CryptographicSponge,
{
    fn accumulate_elems<'a>(
        combined_comms: &mut BTreeMap<Option<usize>, E::G1>,
        combined_witness: &mut E::G1,
        combined_adjusted_witness: &mut E::G1,
        vk: &VerifierKey<E>,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Commitment<E>>>,
        point: P::Point,
        values: impl IntoIterator<Item = E::ScalarField>,
        proof: &kzg10::Proof<E>,
        opening_challenges: &mut ChallengeGenerator<E::ScalarField, S>,
        randomizer: Option<E::ScalarField>,
    ) {
        let acc_time = start_timer!(|| "Accumulating elements");

        let mut curr_challenge = opening_challenges.try_next_challenge_of_size(CHALLENGE_SIZE);

        // Keeps track of running combination of values
        let mut combined_values = E::ScalarField::zero();

        // Iterates through all of the commitments and accumulates common degree_bound elements in a BTreeMap
        for (labeled_comm, value) in commitments.into_iter().zip(values) {
            combined_values += &(value * &curr_challenge);

            let comm = labeled_comm.commitment();
            let degree_bound = labeled_comm.degree_bound();

            // Applying opening challenge and randomness (used in batch_checking)
            let mut comm_with_challenge: E::G1 = comm.0.mul(curr_challenge);

            if let Some(randomizer) = randomizer {
                comm_with_challenge = comm_with_challenge.mul(&randomizer);
            }

            // Accumulate values in the BTreeMap
            *combined_comms.entry(degree_bound).or_insert(E::G1::zero()) += &comm_with_challenge;
            curr_challenge = opening_challenges.try_next_challenge_of_size(CHALLENGE_SIZE);
        }

        // Push expected results into list of elems. Power will be the negative of the expected power
        let mut witness: E::G1 = proof.w.into_group();
        let mut adjusted_witness = vk.g.mul(combined_values) - &proof.w.mul(point);
        if let Some(random_v) = proof.random_v {
            adjusted_witness += &vk.gamma_g.mul(random_v);
        }

        if let Some(randomizer) = randomizer {
            witness = proof.w.mul(randomizer);
            adjusted_witness = adjusted_witness.mul(&randomizer);
        }

        *combined_witness += &witness;
        *combined_adjusted_witness += &adjusted_witness;
        end_timer!(acc_time);
    }

    fn check_elems(
        combined_comms: BTreeMap<Option<usize>, E::G1>,
        combined_witness: E::G1,
        combined_adjusted_witness: E::G1,
        vk: &VerifierKey<E>,
    ) -> Result<bool, Error> {
        let check_time = start_timer!(|| "Checking elems");
        let mut g1_projective_elems: Vec<E::G1> = Vec::new();
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

        let g1_prepared_elems_iter: Vec<E::G1Prepared> =
            E::G1::normalize_batch(g1_projective_elems.as_slice())
                .into_iter()
                .map(|a| a.into())
                .collect::<Vec<_>>();

        let is_one: bool = E::multi_pairing(g1_prepared_elems_iter, g2_prepared_elems)
            .0
            .is_one();
        end_timer!(check_time);
        Ok(is_one)
    }
}

impl<E, P, S> PolynomialCommitment<E::ScalarField, P, S> for SonicKZG10<E, P, S>
where
    E: Pairing,
    P: DenseUVPolynomial<E::ScalarField, Point = E::ScalarField>,
    S: CryptographicSponge,
    for<'a, 'b> &'a P: Div<&'b P, Output = P>,
{
    type UniversalParams = UniversalParams<E>;
    type CommitterKey = CommitterKey<E>;
    type VerifierKey = VerifierKey<E>;
    type Commitment = Commitment<E>;
    type CommitmentState = Randomness<E::ScalarField, P>;
    type Proof = kzg10::Proof<E>;
    type BatchProof = Vec<Self::Proof>;
    type Error = Error;

    fn setup<R: RngCore>(
        max_degree: usize,
        _: Option<usize>,
        rng: &mut R,
    ) -> Result<Self::UniversalParams, Self::Error> {
        kzg10::KZG10::<E, P>::setup(max_degree, true, rng).map_err(Into::into)
    }

    fn trim(
        pp: &Self::UniversalParams,
        supported_degree: usize,
        supported_hiding_bound: usize,
        enforced_degree_bounds: Option<&[usize]>,
    ) -> Result<(Self::CommitterKey, Self::VerifierKey), Self::Error> {
        let trim_time = start_timer!(|| "Trimming public parameters");
        let neg_powers_of_h = &pp.neg_powers_of_h;
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

        let (shifted_powers_of_g, shifted_powers_of_gamma_g, degree_bounds_and_neg_powers_of_h) =
            if let Some(enforced_degree_bounds) = enforced_degree_bounds.as_ref() {
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
                    let mut shifted_powers_of_gamma_g = BTreeMap::new();
                    // Also add degree 0.
                    for degree_bound in enforced_degree_bounds {
                        let shift_degree = max_degree - degree_bound;
                        let mut powers_for_degree_bound = vec![];
                        for i in 0..=(supported_hiding_bound + 1) {
                            // We have an additional degree in `powers_of_gamma_g` beyond `powers_of_g`.
                            if shift_degree + i < max_degree + 2 {
                                powers_for_degree_bound
                                    .push(pp.powers_of_gamma_g[&(shift_degree + i)]);
                            }
                        }
                        shifted_powers_of_gamma_g.insert(*degree_bound, powers_for_degree_bound);
                    }

                    end_timer!(shifted_ck_time);

                    let neg_powers_of_h_time = start_timer!(|| format!(
                        "Constructing `neg_powers_of_h` of size {}",
                        enforced_degree_bounds.len()
                    ));

                    let degree_bounds_and_neg_powers_of_h = enforced_degree_bounds
                        .iter()
                        .map(|bound| (*bound, neg_powers_of_h[&(max_degree - *bound)].clone()))
                        .collect();

                    end_timer!(neg_powers_of_h_time);

                    (
                        Some(shifted_powers_of_g),
                        Some(shifted_powers_of_gamma_g),
                        Some(degree_bounds_and_neg_powers_of_h),
                    )
                }
            } else {
                (None, None, None)
            };

        let powers_of_g = pp.powers_of_g[..=supported_degree].to_vec();
        let powers_of_gamma_g = (0..=(supported_hiding_bound + 1))
            .map(|i| pp.powers_of_gamma_g[&i])
            .collect();

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
        let gamma_g = pp.powers_of_gamma_g[&0];
        let prepared_h = (&pp.prepared_h).clone();
        let prepared_beta_h = (&pp.prepared_beta_h).clone();

        let vk = VerifierKey {
            g,
            gamma_g,
            h,
            beta_h,
            prepared_h,
            prepared_beta_h,
            degree_bounds_and_neg_powers_of_h,
            supported_degree,
            max_degree,
        };

        end_timer!(trim_time);
        Ok((ck, vk))
    }

    /// Outputs a commitment to `polynomial`.
    fn commit<'a>(
        ck: &Self::CommitterKey,
        polynomials: impl IntoIterator<Item = &'a LabeledPolynomial<E::ScalarField, P>>,
        rng: Option<&mut dyn RngCore>,
    ) -> Result<
        (
            Vec<LabeledCommitment<Self::Commitment>>,
            Vec<Self::CommitmentState>,
        ),
        Self::Error,
    >
    where
        P: 'a,
    {
        let rng = &mut crate::optional_rng::OptionalRng(rng);
        let commit_time = start_timer!(|| "Committing to polynomials");
        let mut labeled_comms: Vec<LabeledCommitment<Self::Commitment>> = Vec::new();
        let mut randomness: Vec<Self::CommitmentState> = Vec::new();

        for labeled_polynomial in polynomials {
            let enforced_degree_bounds: Option<&[usize]> = ck
                .enforced_degree_bounds
                .as_ref()
                .map(|bounds| bounds.as_slice());

            kzg10::KZG10::<E, P>::check_degrees_and_bounds(
                ck.supported_degree(),
                ck.max_degree,
                enforced_degree_bounds,
                &labeled_polynomial,
            )?;

            let polynomial: &P = labeled_polynomial.polynomial();
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

            let (comm, rand) = kzg10::KZG10::commit(&powers, polynomial, hiding_bound, Some(rng))?;

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
        labeled_polynomials: impl IntoIterator<Item = &'a LabeledPolynomial<E::ScalarField, P>>,
        _commitments: impl IntoIterator<Item = &'a LabeledCommitment<Self::Commitment>>,
        point: &'a P::Point,
        opening_challenges: &mut ChallengeGenerator<E::ScalarField, S>,
        rands: impl IntoIterator<Item = &'a Self::CommitmentState>,
        _rng: Option<&mut dyn RngCore>,
    ) -> Result<Self::Proof, Self::Error>
    where
        Self::CommitmentState: 'a,
        Self::Commitment: 'a,
        P: 'a,
    {
        let mut combined_polynomial = P::zero();
        let mut combined_rand = kzg10::Randomness::empty();

        let mut curr_challenge = opening_challenges.try_next_challenge_of_size(CHALLENGE_SIZE);

        for (polynomial, rand) in labeled_polynomials.into_iter().zip(rands) {
            let enforced_degree_bounds: Option<&[usize]> = ck
                .enforced_degree_bounds
                .as_ref()
                .map(|bounds| bounds.as_slice());

            kzg10::KZG10::<E, P>::check_degrees_and_bounds(
                ck.supported_degree(),
                ck.max_degree,
                enforced_degree_bounds,
                &polynomial,
            )?;

            combined_polynomial += (curr_challenge, polynomial.polynomial());
            combined_rand += (curr_challenge, rand);
            curr_challenge = opening_challenges.try_next_challenge_of_size(CHALLENGE_SIZE);
        }

        let proof_time = start_timer!(|| "Creating proof for polynomials");
        let proof = kzg10::KZG10::open(&ck.powers(), &combined_polynomial, *point, &combined_rand)?;
        end_timer!(proof_time);

        Ok(proof)
    }

    fn check<'a>(
        vk: &Self::VerifierKey,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Self::Commitment>>,
        point: &'a P::Point,
        values: impl IntoIterator<Item = E::ScalarField>,
        proof: &Self::Proof,
        opening_challenges: &mut ChallengeGenerator<E::ScalarField, S>,
        _rng: Option<&mut dyn RngCore>,
    ) -> Result<bool, Self::Error>
    where
        Self::Commitment: 'a,
    {
        let check_time = start_timer!(|| "Checking evaluations");
        let mut combined_comms: BTreeMap<Option<usize>, E::G1> = BTreeMap::new();
        let mut combined_witness: E::G1 = E::G1::zero();
        let mut combined_adjusted_witness: E::G1 = E::G1::zero();

        Self::accumulate_elems(
            &mut combined_comms,
            &mut combined_witness,
            &mut combined_adjusted_witness,
            vk,
            commitments,
            *point,
            values,
            proof,
            opening_challenges,
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
        query_set: &QuerySet<P::Point>,
        values: &Evaluations<E::ScalarField, P::Point>,
        proof: &Self::BatchProof,
        opening_challenges: &mut ChallengeGenerator<E::ScalarField, S>,
        rng: &mut R,
    ) -> Result<bool, Self::Error>
    where
        Self::Commitment: 'a,
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

        let mut randomizer = E::ScalarField::one();

        let mut combined_comms: BTreeMap<Option<usize>, E::G1> = BTreeMap::new();
        let mut combined_witness: E::G1 = E::G1::zero();
        let mut combined_adjusted_witness: E::G1 = E::G1::zero();

        for ((_point_label, (point, labels)), p) in query_to_labels_map.into_iter().zip(proof) {
            let mut comms_to_combine: Vec<&'_ LabeledCommitment<_>> = Vec::new();
            let mut values_to_combine = Vec::new();
            for label in labels.into_iter() {
                let commitment = commitments.get(label).ok_or(Error::MissingPolynomial {
                    label: label.to_string(),
                })?;

                let v_i = values
                    .get(&(label.clone(), *point))
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
                *point,
                values_to_combine.into_iter(),
                p,
                opening_challenges,
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

    fn open_combinations<'a>(
        ck: &Self::CommitterKey,
        linear_combinations: impl IntoIterator<Item = &'a LinearCombination<E::ScalarField>>,
        polynomials: impl IntoIterator<Item = &'a LabeledPolynomial<E::ScalarField, P>>,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Self::Commitment>>,
        query_set: &QuerySet<P::Point>,
        opening_challenges: &mut ChallengeGenerator<E::ScalarField, S>,
        rands: impl IntoIterator<Item = &'a Self::CommitmentState>,
        rng: Option<&mut dyn RngCore>,
    ) -> Result<BatchLCProof<E::ScalarField, Self::BatchProof>, Self::Error>
    where
        Self::CommitmentState: 'a,
        Self::Commitment: 'a,
        P: 'a,
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

        for lc in linear_combinations {
            let lc_label = lc.label().clone();
            let mut poly = P::zero();
            let mut degree_bound = None;
            let mut hiding_bound = None;
            let mut randomness = Self::CommitmentState::empty();
            let mut comm = E::G1::zero();

            let num_polys = lc.len();
            for (coeff, label) in lc.iter().filter(|(_, l)| !l.is_one()) {
                let label: &String = label.try_into().expect("cannot be one!");
                let &(cur_poly, cur_rand, curr_comm) =
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
                    return Err(Self::Error::EquationHasDegreeBounds(lc_label));
                }

                // Some(_) > None, always.
                hiding_bound = core::cmp::max(hiding_bound, cur_poly.hiding_bound());
                poly += (*coeff, cur_poly.polynomial());
                randomness += (*coeff, cur_rand);
                comm += &curr_comm.commitment().0.mul(*coeff);
            }

            let lc_poly =
                LabeledPolynomial::new(lc_label.clone(), poly, degree_bound, hiding_bound);
            lc_polynomials.push(lc_poly);
            lc_randomness.push(randomness);
            lc_commitments.push(comm);
            lc_info.push((lc_label, degree_bound));
        }

        let comms: Vec<Self::Commitment> = E::G1::normalize_batch(&lc_commitments)
            .into_iter()
            .map(|c| kzg10::Commitment::<E>(c))
            .collect();

        let lc_commitments = lc_info
            .into_iter()
            .zip(comms)
            .map(|((label, d), c)| LabeledCommitment::new(label, c, d))
            .collect::<Vec<_>>();

        let proof = Self::batch_open(
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
    fn check_combinations<'a, R: RngCore>(
        vk: &Self::VerifierKey,
        linear_combinations: impl IntoIterator<Item = &'a LinearCombination<E::ScalarField>>,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Self::Commitment>>,
        eqn_query_set: &QuerySet<P::Point>,
        eqn_evaluations: &Evaluations<P::Point, E::ScalarField>,
        proof: &BatchLCProof<E::ScalarField, Self::BatchProof>,
        opening_challenges: &mut ChallengeGenerator<E::ScalarField, S>,
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
        let mut evaluations = eqn_evaluations.clone();
        for lc in linear_combinations {
            let lc_label = lc.label().clone();
            let num_polys = lc.len();

            let mut degree_bound = None;
            let mut combined_comm = E::G1::zero();

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

        let comms: Vec<Self::Commitment> = E::G1::normalize_batch(&lc_commitments)
            .into_iter()
            .map(|c| kzg10::Commitment(c))
            .collect();

        let lc_commitments = lc_info
            .into_iter()
            .zip(comms)
            .map(|((label, d), c)| LabeledCommitment::new(label, c, d))
            .collect::<Vec<_>>();

        Self::batch_check(
            vk,
            &lc_commitments,
            &eqn_query_set,
            &evaluations,
            proof,
            opening_challenges,
            rng,
        )
    }
}

#[cfg(test)]
mod tests {
    #![allow(non_camel_case_types)]
    use super::SonicKZG10;
    use ark_bls12_377::Bls12_377;
    use ark_bls12_381::Bls12_381;
    use ark_crypto_primitives::sponge::poseidon::PoseidonSponge;
    use ark_ec::pairing::Pairing;
    use ark_ff::UniformRand;
    use ark_poly::{univariate::DensePolynomial as DensePoly, DenseUVPolynomial};
    use rand_chacha::ChaCha20Rng;

    type UniPoly_381 = DensePoly<<Bls12_381 as Pairing>::ScalarField>;
    type UniPoly_377 = DensePoly<<Bls12_377 as Pairing>::ScalarField>;

    type PC<E, P, S> = SonicKZG10<E, P, S>;
    type Sponge_Bls12_377 = PoseidonSponge<<Bls12_377 as Pairing>::ScalarField>;
    type Sponge_Bls12_381 = PoseidonSponge<<Bls12_381 as Pairing>::ScalarField>;
    type PC_Bls12_377 = PC<Bls12_377, UniPoly_377, Sponge_Bls12_377>;
    type PC_Bls12_381 = PC<Bls12_381, UniPoly_381, Sponge_Bls12_381>;

    fn rand_poly<E: Pairing>(
        degree: usize,
        _: Option<usize>,
        rng: &mut ChaCha20Rng,
    ) -> DensePoly<E::ScalarField> {
        DensePoly::<E::ScalarField>::rand(degree, rng)
    }

    fn rand_point<E: Pairing>(_: Option<usize>, rng: &mut ChaCha20Rng) -> E::ScalarField {
        E::ScalarField::rand(rng)
    }

    #[test]
    fn single_poly_test() {
        use crate::tests::*;
        single_poly_test::<_, _, PC_Bls12_377, _>(
            None,
            rand_poly::<Bls12_377>,
            rand_point::<Bls12_377>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bls12-377");
        single_poly_test::<_, _, PC_Bls12_381, _>(
            None,
            rand_poly::<Bls12_381>,
            rand_point::<Bls12_381>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bls12-381");
    }

    #[test]
    fn quadratic_poly_degree_bound_multiple_queries_test() {
        use crate::tests::*;
        quadratic_poly_degree_bound_multiple_queries_test::<_, _, PC_Bls12_377, _>(
            rand_poly::<Bls12_377>,
            rand_point::<Bls12_377>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bls12-377");
        quadratic_poly_degree_bound_multiple_queries_test::<_, _, PC_Bls12_381, _>(
            rand_poly::<Bls12_381>,
            rand_point::<Bls12_381>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bls12-381");
    }

    #[test]
    fn linear_poly_degree_bound_test() {
        use crate::tests::*;
        linear_poly_degree_bound_test::<_, _, PC_Bls12_377, _>(
            rand_poly::<Bls12_377>,
            rand_point::<Bls12_377>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bls12-377");
        linear_poly_degree_bound_test::<_, _, PC_Bls12_381, _>(
            rand_poly::<Bls12_381>,
            rand_point::<Bls12_381>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bls12-381");
    }

    #[test]
    fn single_poly_degree_bound_test() {
        use crate::tests::*;
        single_poly_degree_bound_test::<_, _, PC_Bls12_377, _>(
            rand_poly::<Bls12_377>,
            rand_point::<Bls12_377>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bls12-377");
        single_poly_degree_bound_test::<_, _, PC_Bls12_381, _>(
            rand_poly::<Bls12_381>,
            rand_point::<Bls12_381>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bls12-381");
    }

    #[test]
    fn single_poly_degree_bound_multiple_queries_test() {
        use crate::tests::*;
        single_poly_degree_bound_multiple_queries_test::<_, _, PC_Bls12_377, _>(
            rand_poly::<Bls12_377>,
            rand_point::<Bls12_377>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bls12-377");
        single_poly_degree_bound_multiple_queries_test::<_, _, PC_Bls12_381, _>(
            rand_poly::<Bls12_381>,
            rand_point::<Bls12_381>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bls12-381");
    }

    #[test]
    fn two_polys_degree_bound_single_query_test() {
        use crate::tests::*;
        two_polys_degree_bound_single_query_test::<_, _, PC_Bls12_377, _>(
            rand_poly::<Bls12_377>,
            rand_point::<Bls12_377>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bls12-377");
        two_polys_degree_bound_single_query_test::<_, _, PC_Bls12_381, _>(
            rand_poly::<Bls12_381>,
            rand_point::<Bls12_381>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bls12-381");
    }

    #[test]
    fn full_end_to_end_test() {
        use crate::tests::*;
        full_end_to_end_test::<_, _, PC_Bls12_377, _>(
            None,
            rand_poly::<Bls12_377>,
            rand_point::<Bls12_377>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bls12-377");
        println!("Finished bls12-377");
        full_end_to_end_test::<_, _, PC_Bls12_381, _>(
            None,
            rand_poly::<Bls12_381>,
            rand_point::<Bls12_381>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bls12-381");
        println!("Finished bls12-381");
    }

    #[test]
    fn single_equation_test() {
        use crate::tests::*;
        single_equation_test::<_, _, PC_Bls12_377, _>(
            None,
            rand_poly::<Bls12_377>,
            rand_point::<Bls12_377>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bls12-377");
        println!("Finished bls12-377");
        single_equation_test::<_, _, PC_Bls12_381, _>(
            None,
            rand_poly::<Bls12_381>,
            rand_point::<Bls12_381>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bls12-381");
        println!("Finished bls12-381");
    }

    #[test]
    fn two_equation_test() {
        use crate::tests::*;
        two_equation_test::<_, _, PC_Bls12_377, _>(
            None,
            rand_poly::<Bls12_377>,
            rand_point::<Bls12_377>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bls12-377");
        println!("Finished bls12-377");
        two_equation_test::<_, _, PC_Bls12_381, _>(
            None,
            rand_poly::<Bls12_381>,
            rand_point::<Bls12_381>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bls12-381");
        println!("Finished bls12-381");
    }

    #[test]
    fn two_equation_degree_bound_test() {
        use crate::tests::*;
        two_equation_degree_bound_test::<_, _, PC_Bls12_377, _>(
            rand_poly::<Bls12_377>,
            rand_point::<Bls12_377>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bls12-377");
        println!("Finished bls12-377");
        two_equation_degree_bound_test::<_, _, PC_Bls12_381, _>(
            rand_poly::<Bls12_381>,
            rand_point::<Bls12_381>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bls12-381");
        println!("Finished bls12-381");
    }

    #[test]
    fn full_end_to_end_equation_test() {
        use crate::tests::*;
        full_end_to_end_equation_test::<_, _, PC_Bls12_377, _>(
            None,
            rand_poly::<Bls12_377>,
            rand_point::<Bls12_377>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bls12-377");
        println!("Finished bls12-377");
        full_end_to_end_equation_test::<_, _, PC_Bls12_381, _>(
            None,
            rand_poly::<Bls12_381>,
            rand_point::<Bls12_381>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bls12-381");
        println!("Finished bls12-381");
    }

    #[test]
    #[should_panic]
    fn bad_degree_bound_test() {
        use crate::tests::*;
        bad_degree_bound_test::<_, _, PC_Bls12_377, _>(
            rand_poly::<Bls12_377>,
            rand_point::<Bls12_377>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bls12-377");
        println!("Finished bls12-377");
        bad_degree_bound_test::<_, _, PC_Bls12_381, _>(
            rand_poly::<Bls12_381>,
            rand_point::<Bls12_381>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bls12-381");
        println!("Finished bls12-381");
    }
}
