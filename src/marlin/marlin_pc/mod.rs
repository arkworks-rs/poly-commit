use crate::{kzg10, marlin::Marlin, PCCommitterKey};
use crate::{BTreeMap, BTreeSet, ToString, Vec};
use crate::{BatchLCProof, Error, Evaluations, QuerySet};
use crate::{LabeledCommitment, LabeledPolynomial, LinearCombination};
use crate::{PCRandomness, PCUniversalParams, PolynomialCommitment};
use ark_ec::{AffineCurve, PairingEngine, ProjectiveCurve};
use ark_ff::Zero;
use ark_poly::UVPolynomial;
use ark_std::rand::RngCore;
use ark_std::{marker::PhantomData, ops::Div, vec};

mod data_structures;
pub use data_structures::*;

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
            g: pp.powers_of_g[0].clone(),
            gamma_g: pp.powers_of_gamma_g[&0],
            h: pp.h.clone(),
            beta_h: pp.beta_h.clone(),
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
                    let mut sorted_enforced_degree_bounds = enforced_degree_bounds.clone();
                    sorted_enforced_degree_bounds.sort();

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

            let enforced_degree_bounds: Option<&[usize]> = ck
                .enforced_degree_bounds
                .as_ref()
                .map(|bounds| bounds.as_slice());
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
    fn open_individual_opening_challenges<'a>(
        ck: &Self::CommitterKey,
        labeled_polynomials: impl IntoIterator<Item = &'a LabeledPolynomial<E::Fr, P>>,
        _commitments: impl IntoIterator<Item = &'a LabeledCommitment<Self::Commitment>>,
        point: &'a P::Point,
        opening_challenges: &dyn Fn(u64) -> E::Fr,
        rands: impl IntoIterator<Item = &'a Self::Randomness>,
        _rng: Option<&mut dyn RngCore>,
    ) -> Result<Self::Proof, Self::Error>
    where
        P: 'a,
        Self::Randomness: 'a,
        Self::Commitment: 'a,
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
        vk: &Self::VerifierKey,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Self::Commitment>>,
        point: &'a P::Point,
        values: impl IntoIterator<Item = E::Fr>,
        proof: &Self::Proof,
        opening_challenges: &dyn Fn(u64) -> E::Fr,
        _rng: Option<&mut dyn RngCore>,
    ) -> Result<bool, Self::Error>
    where
        Self::Commitment: 'a,
    {
        let check_time = start_timer!(|| "Checking evaluations");
        let (combined_comm, combined_value) =
            Marlin::accumulate_commitments_and_values_individual_opening_challenges(
                commitments,
                values,
                opening_challenges,
                Some(vk),
            )?;
        let combined_comm = kzg10::Commitment(combined_comm.into());
        let result = kzg10::KZG10::check(&vk.vk, &combined_comm, *point, combined_value, proof)?;
        end_timer!(check_time);
        Ok(result)
    }

    fn batch_check_individual_opening_challenges<'a, R: RngCore>(
        vk: &Self::VerifierKey,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Self::Commitment>>,
        query_set: &QuerySet<P::Point>,
        values: &Evaluations<E::Fr, P::Point>,
        proof: &Self::BatchProof,
        opening_challenges: &dyn Fn(u64) -> E::Fr,
        rng: &mut R,
    ) -> Result<bool, Self::Error>
    where
        Self::Commitment: 'a,
    {
        let (combined_comms, combined_queries, combined_evals) = Marlin::combine_and_normalize(
            commitments,
            query_set,
            values,
            opening_challenges,
            Some(vk),
        )?;
        assert_eq!(proof.len(), combined_queries.len());
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
        ck: &Self::CommitterKey,
        lc_s: impl IntoIterator<Item = &'a LinearCombination<E::Fr>>,
        polynomials: impl IntoIterator<Item = &'a LabeledPolynomial<E::Fr, P>>,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Self::Commitment>>,
        query_set: &QuerySet<P::Point>,
        opening_challenges: &dyn Fn(u64) -> E::Fr,
        rands: impl IntoIterator<Item = &'a Self::Randomness>,
        rng: Option<&mut dyn RngCore>,
    ) -> Result<BatchLCProof<E::Fr, P, Self>, Self::Error>
    where
        P: 'a,
        Self::Randomness: 'a,
        Self::Commitment: 'a,
    {
        Marlin::open_combinations_individual_opening_challenges(
            ck,
            lc_s,
            polynomials,
            commitments,
            query_set,
            opening_challenges,
            rands,
            rng,
        )
    }

    /// Checks that `values` are the true evaluations at `query_set` of the polynomials
    /// committed in `labeled_commitments`.
    fn check_combinations_individual_opening_challenges<'a, R: RngCore>(
        vk: &Self::VerifierKey,
        lc_s: impl IntoIterator<Item = &'a LinearCombination<E::Fr>>,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Self::Commitment>>,
        query_set: &QuerySet<P::Point>,
        evaluations: &Evaluations<E::Fr, P::Point>,
        proof: &BatchLCProof<E::Fr, P, Self>,
        opening_challenges: &dyn Fn(u64) -> E::Fr,
        rng: &mut R,
    ) -> Result<bool, Self::Error>
    where
        Self::Commitment: 'a,
    {
        Marlin::check_combinations_individual_opening_challenges(
            vk,
            lc_s,
            commitments,
            query_set,
            evaluations,
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
        query_set: &QuerySet<P::Point>,
        opening_challenges: &dyn Fn(u64) -> E::Fr,
        rands: impl IntoIterator<Item = &'a Self::Randomness>,
        rng: Option<&mut dyn RngCore>,
    ) -> Result<Vec<kzg10::Proof<E>>, Error>
    where
        P: 'a,
        Self::Randomness: 'a,
        Self::Commitment: 'a,
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
            let mut query_rands: Vec<&'a Self::Randomness> = Vec::new();
            let mut query_comms: Vec<&'a LabeledCommitment<Self::Commitment>> = Vec::new();

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

        Ok(proofs.into())
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
    use ark_std::rand::rngs::StdRng;

    type UniPoly_381 = DensePoly<<Bls12_381 as PairingEngine>::Fr>;
    type UniPoly_377 = DensePoly<<Bls12_377 as PairingEngine>::Fr>;

    type PC<E, P> = MarlinKZG10<E, P>;
    type PC_Bls12_381 = PC<Bls12_381, UniPoly_381>;
    type PC_Bls12_377 = PC<Bls12_377, UniPoly_377>;

    fn rand_poly<E: PairingEngine>(
        degree: usize,
        _: Option<usize>,
        rng: &mut StdRng,
    ) -> DensePoly<E::Fr> {
        DensePoly::<E::Fr>::rand(degree, rng)
    }

    fn constant_poly<E: PairingEngine>(
        _: usize,
        _: Option<usize>,
        rng: &mut StdRng,
    ) -> DensePoly<E::Fr> {
        DensePoly::<E::Fr>::from_coefficients_slice(&[E::Fr::rand(rng)])
    }

    fn rand_point<E: PairingEngine>(_: Option<usize>, rng: &mut StdRng) -> E::Fr {
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
    fn constant_poly_test() {
        use crate::tests::*;
        single_poly_test::<_, _, PC_Bls12_377>(
            None,
            constant_poly::<Bls12_377>,
            rand_point::<Bls12_377>,
        )
        .expect("test failed for bls12-377");
        single_poly_test::<_, _, PC_Bls12_381>(
            None,
            constant_poly::<Bls12_381>,
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
