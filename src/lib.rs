#![cfg_attr(not(feature = "std"), no_std)]
//! A crate for polynomial commitment schemes.
#![deny(unused_import_braces, unused_qualifications, trivial_casts)]
#![deny(trivial_numeric_casts, private_in_public, variant_size_differences)]
#![deny(stable_features, unreachable_pub, non_shorthand_field_patterns)]
#![deny(unused_attributes, unused_mut)]
#![deny(missing_docs)]
#![deny(unused_imports)]
#![deny(renamed_and_removed_lints, stable_features, unused_allocation)]
#![deny(unused_comparisons, bare_trait_objects, unused_must_use, const_err)]
#![forbid(unsafe_code)]

#[macro_use]
extern crate derivative;
#[macro_use]
extern crate ark_std;

use ark_ff::{Field, PrimeField};
pub use ark_poly::{Polynomial, UVPolynomial};
use ark_std::rand::RngCore;

use ark_std::{
    collections::{BTreeMap, BTreeSet},
    fmt::Debug,
    hash::Hash,
    iter::FromIterator,
    rc::Rc,
    string::{String, ToString},
    vec::Vec,
};

/// Data structures used by a polynomial commitment scheme.
pub mod data_structures;
pub use data_structures::*;

/// R1CS constraints for polynomial constraints.
#[cfg(feature = "r1cs")]
mod constraints;
#[cfg(feature = "r1cs")]
pub use constraints::*;

/// Errors pertaining to query sets.
pub mod error;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
pub use error::*;

/// Univariate and multivariate polynomial commitment schemes
/// which (optionally) enable hiding commitments by following
/// the approach outlined in [[CHMMVW20, "Marlin"]][marlin].
///
/// [marlin]: https://eprint.iacr.org/2019/1047
pub mod marlin;

/// A random number generator that bypasses some limitations of the Rust borrow
/// checker.
pub mod optional_rng;

#[cfg(not(feature = "std"))]
macro_rules! eprintln {
    () => {};
    ($($arg: tt)*) => {};
}
#[cfg(not(feature = "std"))]
macro_rules! println {
    () => {};
    ($($arg: tt)*) => {};
}
/// The core [[KZG10]][kzg] construction.
///
/// [kzg]: http://cacr.uwaterloo.ca/techreports/2010/cacr2010-10.pdf
// pub mod kzg10; // TODO: refactor me!

/// Polynomial commitment scheme from [[KZG10]][kzg] that enforces
/// strict degree bounds and (optionally) enables hiding commitments by
/// following the approach outlined in [[CHMMVW20, "Marlin"]][marlin].
///
/// [kzg]: http://cacr.uwaterloo.ca/techreports/2010/cacr2010-10.pdf
/// [marlin]: https://eprint.iacr.org/2019/1047
// pub use marlin::marlin_pc; // TODO: refactor me!

/// Polynomial commitment scheme based on the construction in [[KZG10]][kzg],
/// modified to obtain batching and to enforce strict
/// degree bounds by following the approach outlined in [[MBKM19,
/// “Sonic”]][sonic] (more precisely, via the variant in
/// [[Gabizon19, “AuroraLight”]][al] that avoids negative G1 powers).
///
/// [kzg]: http://cacr.uwaterloo.ca/techreports/2010/cacr2010-10.pdf
/// [sonic]: https://eprint.iacr.org/2019/099
/// [al]: https://eprint.iacr.org/2019/601
/// [marlin]: https://eprint.iacr.org/2019/1047
// pub mod sonic_pc; // TODO: refactor me!

/// A polynomial commitment scheme based on the hardness of the
/// discrete logarithm problem in prime-order groups.
/// The construction is detailed in [[BCMS20]][pcdas].
///
/// [pcdas]: https://eprint.iacr.org/2020/499
// pub mod ipa_pc; // TODO: refactor me!

/// Defines the challenge strategies and challenge generator.
pub mod challenge;
/// A multilinear polynomial commitment scheme that converts n-variate multilinear polynomial into
/// n quotient UV polynomial. This scheme is based on hardness of the discrete logarithm
/// in prime-order groups. Construction is detailed in [[XZZPD19]][xzzpd19] and [[ZGKPP18]][zgkpp18]
///
/// [xzzpd19]: https://eprint.iacr.org/2019/317
/// [zgkpp]: https://ieeexplore.ieee.org/document/8418645
pub mod multilinear_pc; // TODO: no need to refactor, but still give a double check.

/// Multivariate polynomial commitment based on the construction in
/// [[PST13]][pst] with batching and (optional) hiding property inspired
/// by the univariate scheme in [[CHMMVW20, "Marlin"]][marlin]
///
/// [pst]: https://eprint.iacr.org/2011/587.pdf
/// [marlin]: https://eprint.iacr.org/2019/104
pub use marlin::marlin_pst13_pc;
use crate::challenge::ChallengeGenerator;
use ark_sponge::FieldBasedCryptographicSponge;

/// `QuerySet` is the set of queries that are to be made to a set of labeled polynomials/equations
/// `p` that have previously been committed to. Each element of a `QuerySet` is a pair of
/// `(label, (point_label, point))`, where `label` is the label of a polynomial in `p`,
/// `point_label` is the label for the point (e.g., "beta"), and  and `point` is the location
/// that `p[label]` is to be queried at.
pub type QuerySet<T> = BTreeSet<(String, (String, T))>;

/// `Evaluations` is the result of querying a set of labeled polynomials or equations
/// `p` at a `QuerySet` `Q`. It maps each element of `Q` to the resulting evaluation.
/// That is, if `(label, query)` is an element of `Q`, then `evaluation.get((label, query))`
/// should equal `p[label].evaluate(query)`.
pub type Evaluations<T, F> = BTreeMap<(String, T), F>;

/// Describes the interface for a polynomial commitment scheme that allows
/// a sender to commit to multiple polynomials and later provide a succinct proof
/// of evaluation for the corresponding commitments at a query set `Q`, while
/// enforcing per-polynomial degree bounds.
pub trait PolynomialCommitment<F: PrimeField, P: Polynomial<F>, S: FieldBasedCryptographicSponge<F>>: Sized {
    /// The universal parameters for the commitment scheme. These are "trimmed"
    /// down to `Self::CommitterKey` and `Self::VerifierKey` by `Self::trim`.
    type UniversalParams: PCUniversalParams;
    /// The committer key for the scheme; used to commit to a polynomial and then
    /// open the commitment to produce an evaluation proof.
    type CommitterKey: PCCommitterKey;
    /// The verifier key for the scheme; used to check an evaluation proof.
    type VerifierKey: PCVerifierKey;
    /// The prepared verifier key for the scheme; used to check an evaluation proof.
    type PreparedVerifierKey: PCPreparedVerifierKey<Self::VerifierKey> + Clone;
    /// The commitment to a polynomial.
    type Commitment: PCCommitment + Default;
    /// The prepared commitment to a polynomial.
    type PreparedCommitment: PCPreparedCommitment<Self::Commitment>;
    /// The commitment randomness.
    type Randomness: PCRandomness;
    /// The evaluation proof for a single point.
    type Proof: PCProof + Clone;
    /// The evaluation proof for a query set.
    type BatchProof: Clone
        + From<Vec<Self::Proof>>
        + Into<Vec<Self::Proof>>
        + CanonicalSerialize
        + CanonicalDeserialize;
    /// The error type for the scheme.
    type Error: ark_std::error::Error + From<Error>;

    /// Constructs public parameters when given as input the maximum degree `degree`
    /// for the polynomial commitment scheme. `num_vars` specifies the number of
    /// variables for multivariate setup
    fn setup<R: RngCore>(
        max_degree: usize,
        num_vars: Option<usize>,
        rng: &mut R,
    ) -> Result<Self::UniversalParams, Self::Error>;

    /// Specializes the public parameters for polynomials up to the given `supported_degree`
    /// and for enforcing degree bounds in the range `1..=supported_degree`.
    fn trim(
        pp: &Self::UniversalParams,
        supported_degree: usize,
        supported_hiding_bound: usize,
        enforced_degree_bounds: Option<&[usize]>,
    ) -> Result<(Self::CommitterKey, Self::VerifierKey), Self::Error>;

    /// Outputs a commitments to `polynomials`. If `polynomials[i].is_hiding()`,
    /// then the `i`-th commitment is hiding up to `polynomials.hiding_bound()` queries.
    /// `rng` should not be `None` if `polynomials[i].is_hiding() == true` for any `i`.
    ///
    /// If for some `i`, `polynomials[i].is_hiding() == false`, then the
    /// corresponding randomness is `Self::Randomness::empty()`.
    ///
    /// If for some `i`, `polynomials[i].degree_bound().is_some()`, then that
    /// polynomial will have the corresponding degree bound enforced.
    fn commit<'a>(
        ck: &Self::CommitterKey,
        polynomials: impl IntoIterator<Item = &'a LabeledPolynomial<F, P>>,
        rng: Option<&mut dyn RngCore>,
    ) -> Result<
        (
            Vec<LabeledCommitment<Self::Commitment>>,
            Vec<Self::Randomness>,
        ),
        Self::Error,
    >
    where
        P: 'a;

    /// open but with individual challenges
    fn open<'a>(
        ck: &Self::CommitterKey,
        labeled_polynomials: impl IntoIterator<Item = &'a LabeledPolynomial<F, P>>,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Self::Commitment>>,
        point: &'a P::Point,
        opening_challenges: &mut ChallengeGenerator<F, S>,
        rands: impl IntoIterator<Item = &'a Self::Randomness>,
        rng: Option<&mut dyn RngCore>,
    ) -> Result<Self::Proof, Self::Error>
    where
        P: 'a,
        Self::Randomness: 'a,
        Self::Commitment: 'a;

    /// check but with individual challenges
    fn check<'a>(
        vk: &Self::VerifierKey,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Self::Commitment>>,
        point: &'a P::Point,
        values: impl IntoIterator<Item = F>,
        proof: &Self::Proof,
        opening_challenges: &mut ChallengeGenerator<F, S>,
        rng: Option<&mut dyn RngCore>,
    ) -> Result<bool, Self::Error>
    where
        Self::Commitment: 'a;

    /// batch_check but with individual challenges
    fn batch_check<'a, R: RngCore>(
        vk: &Self::VerifierKey,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Self::Commitment>>,
        query_set: &QuerySet<P::Point>,
        evaluations: &Evaluations<P::Point, F>,
        proof: &Self::BatchProof,
        opening_challenges: &mut ChallengeGenerator<F, S>,
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

        // Implicit assumption: proofs are order in same manner as queries in
        // `query_to_labels_map`.
        let proofs: Vec<_> = proof.clone().into();
        assert_eq!(proofs.len(), query_to_labels_map.len());

        let mut result = true;
        for ((_point_label, (point, labels)), proof) in query_to_labels_map.into_iter().zip(proofs)
        {
            let mut comms: Vec<&'_ LabeledCommitment<_>> = Vec::new();
            let mut values = Vec::new();
            for label in labels.into_iter() {
                let commitment = commitments.get(label).ok_or(Error::MissingPolynomial {
                    label: label.to_string(),
                })?;

                let v_i = evaluations.get(&(label.clone(), point.clone())).ok_or(
                    Error::MissingEvaluation {
                        label: label.to_string(),
                    },
                )?;

                comms.push(commitment);
                values.push(*v_i);
            }

            let proof_time = start_timer!(|| "Checking per-query proof");
            result &= Self::check(
                vk,
                comms,
                &point,
                values,
                &proof,
                opening_challenges,
                Some(rng),
            )?;
            end_timer!(proof_time);
        }
        Ok(result)
    }

    /// open_combinations but with individual challenges
    fn open_combinations<'a>(
        ck: &Self::CommitterKey,
        linear_combinations: impl IntoIterator<Item = &'a LinearCombination<F>>,
        polynomials: impl IntoIterator<Item = &'a LabeledPolynomial<F, P>>,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Self::Commitment>>,
        query_set: &QuerySet<P::Point>,
        opening_challenges: &mut ChallengeGenerator<F, S>,
        rands: impl IntoIterator<Item = &'a Self::Randomness>,
        rng: Option<&mut dyn RngCore>,
    ) -> Result<BatchLCProof<F, P, Self>, Self::Error>
    where
        Self::Randomness: 'a,
        Self::Commitment: 'a,
        P: 'a,
    {
        let linear_combinations: Vec<_> = linear_combinations.into_iter().collect();
        let polynomials: Vec<_> = polynomials.into_iter().collect();
        let poly_query_set =
            lc_query_set_to_poly_query_set(linear_combinations.iter().copied(), query_set);
        let poly_evals = evaluate_query_set(polynomials.iter().copied(), &poly_query_set);
        let proof = Self::batch_open(
            ck,
            polynomials,
            commitments,
            &poly_query_set,
            opening_challenges,
            rands,
            rng,
        )?;
        Ok(BatchLCProof {
            proof,
            evals: Some(poly_evals.values().copied().collect()),
        })
    }

    /// check_combinations with individual challenges
    fn check_combinations<'a, R: RngCore>(
        vk: &Self::VerifierKey,
        linear_combinations: impl IntoIterator<Item = &'a LinearCombination<F>>,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Self::Commitment>>,
        eqn_query_set: &QuerySet<P::Point>,
        eqn_evaluations: &Evaluations<P::Point, F>,
        proof: &BatchLCProof<F, P, Self>,
        opening_challenges: &mut ChallengeGenerator<F, S>,
        rng: &mut R,
    ) -> Result<bool, Self::Error>
    where
        Self::Commitment: 'a,
    {
        let BatchLCProof { proof, evals } = proof;

        let lc_s = BTreeMap::from_iter(linear_combinations.into_iter().map(|lc| (lc.label(), lc)));

        let poly_query_set = lc_query_set_to_poly_query_set(lc_s.values().copied(), eqn_query_set);
        let poly_evals = Evaluations::from_iter(
            poly_query_set
                .iter()
                .map(|(_, point)| point)
                .cloned()
                .zip(evals.clone().unwrap()),
        );

        for &(ref lc_label, (_, ref point)) in eqn_query_set {
            if let Some(lc) = lc_s.get(lc_label) {
                let claimed_rhs = *eqn_evaluations
                    .get(&(lc_label.clone(), point.clone()))
                    .ok_or(Error::MissingEvaluation {
                        label: lc_label.to_string(),
                    })?;

                let mut actual_rhs = F::zero();

                for (coeff, label) in lc.iter() {
                    let eval = match label {
                        LCTerm::One => F::one(),
                        LCTerm::PolyLabel(l) => *poly_evals
                            .get(&(l.clone().into(), point.clone()))
                            .ok_or(Error::MissingEvaluation { label: l.clone() })?,
                    };

                    actual_rhs += &(*coeff * eval);
                }
                if claimed_rhs != actual_rhs {
                    eprintln!("Claimed evaluation of {} is incorrect", lc.label());
                    return Ok(false);
                }
            }
        }

        let pc_result = Self::batch_check(
            vk,
            commitments,
            &poly_query_set,
            &poly_evals,
            proof,
            opening_challenges,
            rng,
        )?;
        if !pc_result {
            eprintln!("Evaluation proofs failed to verify");
            return Ok(false);
        }

        Ok(true)
    }

    /// batch_open with individual challenges
    fn batch_open<'a>(
        ck: &Self::CommitterKey,
        labeled_polynomials: impl IntoIterator<Item = &'a LabeledPolynomial<F, P>>,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Self::Commitment>>,
        query_set: &QuerySet<P::Point>,
        opening_challenges: &mut ChallengeGenerator<F, S>,
        rands: impl IntoIterator<Item = &'a Self::Randomness>,
        rng: Option<&mut dyn RngCore>,
    ) -> Result<Self::BatchProof, Self::Error>
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
                    poly_rand_comm.get(label).ok_or(Error::MissingPolynomial {
                        label: label.to_string(),
                    })?;

                query_polys.push(polynomial);
                query_rands.push(rand);
                query_comms.push(comm);
            }

            let proof_time = start_timer!(|| "Creating proof");
            let proof = Self::open(
                ck,
                query_polys,
                query_comms,
                &point,
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

/// Evaluate the given polynomials at `query_set`.
pub fn evaluate_query_set<'a, F, P, T>(
    polys: impl IntoIterator<Item = &'a LabeledPolynomial<F, P>>,
    query_set: &QuerySet<T>,
) -> Evaluations<T, F>
where
    F: Field,
    P: 'a + Polynomial<F, Point = T>,
    T: Clone + Debug + Hash + Ord + Sync,
{
    let polys = BTreeMap::from_iter(polys.into_iter().map(|p| (p.label(), p)));
    let mut evaluations = Evaluations::new();
    for (label, (_, point)) in query_set {
        let poly = polys
            .get(label)
            .expect("polynomial in evaluated lc is not found");
        let eval = poly.evaluate(&point);
        evaluations.insert((label.clone(), point.clone()), eval);
    }
    evaluations
}

fn lc_query_set_to_poly_query_set<'a, F: Field, T: Clone + Ord>(
    linear_combinations: impl IntoIterator<Item = &'a LinearCombination<F>>,
    query_set: &QuerySet<T>,
) -> QuerySet<T> {
    let mut poly_query_set = QuerySet::<T>::new();
    let lc_s = linear_combinations.into_iter().map(|lc| (lc.label(), lc));
    let linear_combinations = BTreeMap::from_iter(lc_s);
    for (lc_label, (point_label, point)) in query_set {
        if let Some(lc) = linear_combinations.get(lc_label) {
            for (_, poly_label) in lc.iter().filter(|(_, l)| !l.is_one()) {
                if let LCTerm::PolyLabel(l) = poly_label {
                    poly_query_set.insert((l.into(), (point_label.clone(), point.clone())));
                }
            }
        }
    }
    poly_query_set
}

#[cfg(test)]
pub mod tests {
    use crate::*;
    use ark_ff::Field;
    use ark_poly::Polynomial;
    use ark_std::rand::{
        distributions::{Distribution, Uniform},
        rngs::StdRng,
        Rng,
    };
    use ark_std::test_rng;
    use ark_sponge::poseidon::PoseidonParameters;

    struct TestInfo<F: Field, P: Polynomial<F>> {
        num_iters: usize,
        max_degree: Option<usize>,
        supported_degree: Option<usize>,
        num_vars: Option<usize>,
        num_polynomials: usize,
        enforce_degree_bounds: bool,
        max_num_queries: usize,
        num_equations: Option<usize>,
        rand_poly: fn(usize, Option<usize>, &mut StdRng) -> P,
        rand_point: fn(Option<usize>, &mut StdRng) -> P::Point,
        opening_challenge: fn() -> ChallengeGenerator<F, P>
    }

    pub fn bad_degree_bound_test<F, P, PC, S>(
        rand_poly: fn(usize, Option<usize>, &mut StdRng) -> P,
        rand_point: fn(Option<usize>, &mut StdRng) -> P::Point,
        opening_challenge: fn() -> ChallengeGenerator<F, S>
    ) -> Result<(), PC::Error>
    where
        F: PrimeField,
        P: Polynomial<F>,
        PC: PolynomialCommitment<F, P, S>,
        S: FieldBasedCryptographicSponge<F>
    {
        let rng = &mut test_rng();
        let max_degree = 100;
        let pp = PC::setup(max_degree, None, rng)?;
        for _ in 0..10 {
            let supported_degree = Uniform::from(1..=max_degree).sample(rng);
            assert!(
                max_degree >= supported_degree,
                "max_degree < supported_degree"
            );

            let mut labels = Vec::new();
            let mut polynomials = Vec::new();
            let mut degree_bounds = Vec::new();

            for i in 0..10 {
                let label = format!("Test{}", i);
                labels.push(label.clone());
                let degree_bound = 1usize;
                let hiding_bound = Some(1);
                degree_bounds.push(degree_bound);

                polynomials.push(LabeledPolynomial::new(
                    label,
                    rand_poly(supported_degree, None, rng),
                    Some(degree_bound),
                    hiding_bound,
                ));
            }

            let supported_hiding_bound = polynomials
                .iter()
                .map(|p| p.hiding_bound().unwrap_or(0))
                .max()
                .unwrap_or(0);
            println!("supported degree: {:?}", supported_degree);
            println!("supported hiding bound: {:?}", supported_hiding_bound);
            let (ck, vk) = PC::trim(
                &pp,
                supported_degree,
                supported_hiding_bound,
                Some(degree_bounds.as_slice()),
            )?;
            println!("Trimmed");

            let (comms, rands) = PC::commit(&ck, &polynomials, Some(rng))?;

            let mut query_set = QuerySet::new();
            let mut values = Evaluations::new();
            let point = rand_point(None, rng);
            for (i, label) in labels.iter().enumerate() {
                query_set.insert((label.clone(), (format!("{}", i), point.clone())));
                let value = polynomials[i].evaluate(&point);
                values.insert((label.clone(), point.clone()), value);
            }
            println!("Generated query set");

            let proof = PC::batch_open(
                &ck,
                &polynomials,
                &comms,
                &query_set,
                &mut opening_challenge(),
                &rands,
                Some(rng),
            )?;
            let result = PC::batch_check(
                &vk,
                &comms,
                &query_set,
                &values,
                &proof,
                &mut opening_challenge(),
                rng,
            )?;
            assert!(result, "proof was incorrect, Query set: {:#?}", query_set);
        }
        Ok(())
    }

    fn test_template<F, P, PC, S>(info: TestInfo<F, P>) -> Result<(), PC::Error>
    where
        F: Field,
        P: Polynomial<F>,
        PC: PolynomialCommitment<F, P, S>,
        S: FieldBasedCryptographicSponge<F>
    {
        let TestInfo {
            num_iters,
            max_degree,
            supported_degree,
            num_vars,
            num_polynomials,
            enforce_degree_bounds,
            max_num_queries,
            num_equations: _,
            rand_poly,
            rand_point,
            opening_challenge
        } = info;

        let rng = &mut test_rng();
        // If testing multivariate polynomials, make the max degree lower
        let max_degree = match num_vars {
            Some(_) => max_degree.unwrap_or(Uniform::from(2..=10).sample(rng)),
            None => max_degree.unwrap_or(Uniform::from(2..=64).sample(rng)),
        };
        let pp = PC::setup(max_degree, num_vars, rng)?;

        for _ in 0..num_iters {
            let supported_degree =
                supported_degree.unwrap_or(Uniform::from(1..=max_degree).sample(rng));
            assert!(
                max_degree >= supported_degree,
                "max_degree < supported_degree"
            );
            let mut polynomials: Vec<LabeledPolynomial<F, P>> = Vec::new();
            let mut degree_bounds = if enforce_degree_bounds {
                Some(Vec::new())
            } else {
                None
            };

            let mut labels = Vec::new();
            println!("Sampled supported degree");

            // Generate polynomials
            let num_points_in_query_set = Uniform::from(1..=max_num_queries).sample(rng);
            for i in 0..num_polynomials {
                let label = format!("Test{}", i);
                labels.push(label.clone());
                let degree = Uniform::from(1..=supported_degree).sample(rng);
                let degree_bound = if let Some(degree_bounds) = &mut degree_bounds {
                    let range = Uniform::from(degree..=supported_degree);
                    let degree_bound = range.sample(rng);
                    degree_bounds.push(degree_bound);
                    Some(degree_bound)
                } else {
                    None
                };

                let hiding_bound = if num_points_in_query_set >= degree {
                    Some(degree)
                } else {
                    Some(num_points_in_query_set)
                };

                polynomials.push(LabeledPolynomial::new(
                    label,
                    rand_poly(degree, num_vars, rng).into(),
                    degree_bound,
                    hiding_bound,
                ))
            }
            let supported_hiding_bound = polynomials
                .iter()
                .map(|p| p.hiding_bound().unwrap_or(0))
                .max()
                .unwrap_or(0);
            println!("supported degree: {:?}", supported_degree);
            println!("supported hiding bound: {:?}", supported_hiding_bound);
            println!("num_points_in_query_set: {:?}", num_points_in_query_set);
            let (ck, vk) = PC::trim(
                &pp,
                supported_degree,
                supported_hiding_bound,
                degree_bounds.as_ref().map(|s| s.as_slice()),
            )?;
            println!("Trimmed");

            let (comms, rands) = PC::commit(&ck, &polynomials, Some(rng))?;

            // Construct query set
            let mut query_set = QuerySet::new();
            let mut values = Evaluations::new();
            for _ in 0..num_points_in_query_set {
                let point = rand_point(num_vars, rng);
                for (i, label) in labels.iter().enumerate() {
                    query_set.insert((label.clone(), (format!("{}", i), point.clone())));
                    let value = polynomials[i].evaluate(&point);
                    values.insert((label.clone(), point.clone()), value);
                }
            }
            println!("Generated query set");

            let proof = PC::batch_open(
                &ck,
                &polynomials,
                &comms,
                &query_set,
                &mut opening_challenge(),
                &rands,
                Some(rng),
            )?;
            let result = PC::batch_check(
                &vk,
                &comms,
                &query_set,
                &values,
                &proof,
                &mut opening_challenge(),
                rng,
            )?;
            if !result {
                println!(
                    "Failed with {} polynomials, num_points_in_query_set: {:?}",
                    num_polynomials, num_points_in_query_set
                );
                println!("Degree of polynomials:",);
                for poly in polynomials {
                    println!("Degree: {:?}", poly.degree());
                }
            }
            assert!(result, "proof was incorrect, Query set: {:#?}", query_set);
        }
        Ok(())
    }

    fn equation_test_template<F, P, PC, S>(info: TestInfo<F, P>) -> Result<(), PC::Error>
    where
        F: Field,
        P: Polynomial<F>,
        PC: PolynomialCommitment<F, P, S>,
        S: FieldBasedCryptographicSponge<F>
    {
        let TestInfo {
            num_iters,
            max_degree,
            supported_degree,
            num_vars,
            num_polynomials,
            enforce_degree_bounds,
            max_num_queries,
            num_equations,
            rand_poly,
            rand_point,
            opening_challenge
        } = info;

        let rng = &mut test_rng();
        // If testing multivariate polynomials, make the max degree lower
        let max_degree = match num_vars {
            Some(_) => max_degree.unwrap_or(Uniform::from(2..=10).sample(rng)),
            None => max_degree.unwrap_or(Uniform::from(2..=64).sample(rng)),
        };
        let pp = PC::setup(max_degree, num_vars, rng)?;

        for _ in 0..num_iters {
            let supported_degree =
                supported_degree.unwrap_or(Uniform::from(1..=max_degree).sample(rng));
            assert!(
                max_degree >= supported_degree,
                "max_degree < supported_degree"
            );
            let mut polynomials = Vec::new();
            let mut degree_bounds = if enforce_degree_bounds {
                Some(Vec::new())
            } else {
                None
            };

            let mut labels = Vec::new();
            println!("Sampled supported degree");

            // Generate polynomials
            let num_points_in_query_set = Uniform::from(1..=max_num_queries).sample(rng);
            for i in 0..num_polynomials {
                let label = format!("Test{}", i);
                labels.push(label.clone());
                let degree = Uniform::from(1..=supported_degree).sample(rng);
                let degree_bound = if let Some(degree_bounds) = &mut degree_bounds {
                    if rng.gen() {
                        let range = Uniform::from(degree..=supported_degree);
                        let degree_bound = range.sample(rng);
                        degree_bounds.push(degree_bound);
                        Some(degree_bound)
                    } else {
                        None
                    }
                } else {
                    None
                };

                let hiding_bound = if num_points_in_query_set >= degree {
                    Some(degree)
                } else {
                    Some(num_points_in_query_set)
                };
                println!("Hiding bound: {:?}", hiding_bound);

                polynomials.push(LabeledPolynomial::new(
                    label,
                    rand_poly(degree, num_vars, rng),
                    degree_bound,
                    hiding_bound,
                ))
            }
            println!("supported degree: {:?}", supported_degree);
            println!("num_points_in_query_set: {:?}", num_points_in_query_set);
            println!("{:?}", degree_bounds);
            println!("{}", num_polynomials);
            println!("{}", enforce_degree_bounds);

            let (ck, vk) = PC::trim(
                &pp,
                supported_degree,
                supported_degree,
                degree_bounds.as_ref().map(|s| s.as_slice()),
            )?;
            println!("Trimmed");

            let (comms, rands) = PC::commit(&ck, &polynomials, Some(rng))?;

            // Let's construct our equations
            let mut linear_combinations = Vec::new();
            let mut query_set = QuerySet::new();
            let mut values = Evaluations::new();
            for i in 0..num_points_in_query_set {
                let point = rand_point(num_vars, rng);
                for j in 0..num_equations.unwrap() {
                    let label = format!("query {} eqn {}", i, j);
                    let mut lc = LinearCombination::empty(label.clone());

                    let mut value = F::zero();
                    let should_have_degree_bounds: bool = rng.gen();
                    for (k, label) in labels.iter().enumerate() {
                        if should_have_degree_bounds {
                            value += &polynomials[k].evaluate(&point);
                            lc.push((F::one(), label.to_string().into()));
                            break;
                        } else {
                            let poly = &polynomials[k];
                            if poly.degree_bound().is_some() {
                                continue;
                            } else {
                                assert!(poly.degree_bound().is_none());
                                let coeff = F::rand(rng);
                                value += &(coeff * poly.evaluate(&point));
                                lc.push((coeff, label.to_string().into()));
                            }
                        }
                    }
                    values.insert((label.clone(), point.clone()), value);
                    if !lc.is_empty() {
                        linear_combinations.push(lc);
                        // Insert query
                        query_set.insert((label.clone(), (format!("{}", i), point.clone())));
                    }
                }
            }
            if linear_combinations.is_empty() {
                continue;
            }
            println!("Generated query set");
            println!("Linear combinations: {:?}", linear_combinations);

            let proof = PC::open_combinations(
                &ck,
                &linear_combinations,
                &polynomials,
                &comms,
                &query_set,
                &mut opening_challenge(),
                &rands,
                Some(rng),
            )?;
            println!("Generated proof");
            let result = PC::check_combinations(
                &vk,
                &linear_combinations,
                &comms,
                &query_set,
                &values,
                &proof,
                &mut opening_challenge(),
                rng,
            )?;
            if !result {
                println!(
                    "Failed with {} polynomials, num_points_in_query_set: {:?}",
                    num_polynomials, num_points_in_query_set
                );
                println!("Degree of polynomials:",);
                for poly in polynomials {
                    println!("Degree: {:?}", poly.degree());
                }
            }
            assert!(
                result,
                "proof was incorrect, equations: {:#?}",
                linear_combinations
            );
        }
        Ok(())
    }

    pub fn single_poly_test<F, P, PC, S>(
        num_vars: Option<usize>,
        rand_poly: fn(usize, Option<usize>, &mut StdRng) -> P,
        rand_point: fn(Option<usize>, &mut StdRng) -> P::Point,
        opening_challenge: fn() -> ChallengeGenerator<F, P>
    ) -> Result<(), PC::Error>
    where
        F: PrimeField,
        P: Polynomial<F>,
        PC: PolynomialCommitment<F, P, S>,
        S: FieldBasedCryptographicSponge<F>
    {
        let info = TestInfo {
            num_iters: 100,
            max_degree: None,
            supported_degree: None,
            num_vars,
            num_polynomials: 1,
            enforce_degree_bounds: false,
            max_num_queries: 1,
            num_equations: None,
            rand_poly,
            rand_point,
            opening_challenge
        };
        test_template::<F, P, PC, S>(info)
    }

    pub fn linear_poly_degree_bound_test<F, P, PC, S>(
        rand_poly: fn(usize, Option<usize>, &mut StdRng) -> P,
        rand_point: fn(Option<usize>, &mut StdRng) -> P::Point,
        opening_challenge: fn() -> ChallengeGenerator<F,P>
    ) -> Result<(), PC::Error>
    where
        F: Field,
        P: Polynomial<F>,
        PC: PolynomialCommitment<F, P, S>,
        S: FieldBasedCryptographicSponge<F>
    {
        let info = TestInfo {
            num_iters: 100,
            max_degree: Some(2),
            supported_degree: Some(1),
            num_vars: None,
            num_polynomials: 1,
            enforce_degree_bounds: true,
            max_num_queries: 1,
            num_equations: None,
            rand_poly,
            rand_point,
            opening_challenge
        };
        test_template::<F, P, PC, S>(info)
    }

    pub fn single_poly_degree_bound_test<F, P, PC, S>(
        rand_poly: fn(usize, Option<usize>, &mut StdRng) -> P,
        rand_point: fn(Option<usize>, &mut StdRng) -> P::Point,
        opening_challenge: fn() -> ChallengeGenerator<F,P>
    ) -> Result<(), PC::Error>
        where
            F: Field,
            P: Polynomial<F>,
            PC: PolynomialCommitment<F, P, S>,
            S: FieldBasedCryptographicSponge<F>
    {
        let info = TestInfo {
            num_iters: 100,
            max_degree: None,
            supported_degree: None,
            num_vars: None,
            num_polynomials: 1,
            enforce_degree_bounds: true,
            max_num_queries: 1,
            num_equations: None,
            rand_poly,
            rand_point,
            opening_challenge
        };
        test_template::<F, P, PC, S>(info)
    }

    pub fn quadratic_poly_degree_bound_multiple_queries_test<F, P, PC, S>(
        rand_poly: fn(usize, Option<usize>, &mut StdRng) -> P,
        rand_point: fn(Option<usize>, &mut StdRng) -> P::Point,
        opening_challenge: fn() -> ChallengeGenerator<F,P>
    ) -> Result<(), PC::Error>
        where
            F: Field,
            P: Polynomial<F>,
            PC: PolynomialCommitment<F, P, S>,
            S: FieldBasedCryptographicSponge<F>
    {
        let info = TestInfo {
            num_iters: 100,
            max_degree: Some(3),
            supported_degree: Some(2),
            num_vars: None,
            num_polynomials: 1,
            enforce_degree_bounds: true,
            max_num_queries: 2,
            num_equations: None,
            rand_poly,
            rand_point,
            opening_challenge
        };
        test_template::<F, P, PC, S>(info)
    }

    pub fn single_poly_degree_bound_multiple_queries_test<F, P, PC, S>(
        rand_poly: fn(usize, Option<usize>, &mut StdRng) -> P,
        rand_point: fn(Option<usize>, &mut StdRng) -> P::Point,
        opening_challenge: fn() -> ChallengeGenerator<F,P>
    ) -> Result<(), PC::Error>
        where
            F: Field,
            P: Polynomial<F>,
            PC: PolynomialCommitment<F, P, S>,
            S: FieldBasedCryptographicSponge<F>
    {
        let info = TestInfo {
            num_iters: 100,
            max_degree: None,
            supported_degree: None,
            num_vars: None,
            num_polynomials: 1,
            enforce_degree_bounds: true,
            max_num_queries: 2,
            num_equations: None,
            rand_poly,
            rand_point,
            opening_challenge
        };
        test_template::<F, P, PC, S>(info)
    }

    pub fn two_polys_degree_bound_single_query_test<F, P, PC, S>(
        rand_poly: fn(usize, Option<usize>, &mut StdRng) -> P,
        rand_point: fn(Option<usize>, &mut StdRng) -> P::Point,
        opening_challenge: fn() -> ChallengeGenerator<F,P>
    ) -> Result<(), PC::Error>
        where
            F: Field,
            P: Polynomial<F>,
            PC: PolynomialCommitment<F, P, S>,
            S: FieldBasedCryptographicSponge<F>
    {
        let info = TestInfo {
            num_iters: 100,
            max_degree: None,
            supported_degree: None,
            num_vars: None,
            num_polynomials: 2,
            enforce_degree_bounds: true,
            max_num_queries: 1,
            num_equations: None,
            rand_poly,
            rand_point,
            opening_challenge
        };
        test_template::<F, P, PC, S>(info)
    }

    pub fn full_end_to_end_test<F, P, PC, S>(
        num_vars: Option<usize>,
        rand_poly: fn(usize, Option<usize>, &mut StdRng) -> P,
        rand_point: fn(Option<usize>, &mut StdRng) -> P::Point,
        opening_challenge: fn() -> ChallengeGenerator<F,P>
    ) -> Result<(), PC::Error>
    where
        F: Field,
        P: Polynomial<F>,
        PC: PolynomialCommitment<F, P, S>,
    {
        let info = TestInfo {
            num_iters: 100,
            max_degree: None,
            supported_degree: None,
            num_vars,
            num_polynomials: 10,
            enforce_degree_bounds: true,
            max_num_queries: 5,
            num_equations: None,
            rand_poly,
            rand_point,
            opening_challenge,
        };
        test_template::<F, P, PC, S>(info)
    }

    pub fn full_end_to_end_equation_test<F, P, PC, S>(
        num_vars: Option<usize>,
        rand_poly: fn(usize, Option<usize>, &mut StdRng) -> P,
        rand_point: fn(Option<usize>, &mut StdRng) -> P::Point,
        opening_challenge: fn() -> ChallengeGenerator<F,P>
    ) -> Result<(), PC::Error>
    where
        F: Field,
        P: Polynomial<F>,
        PC: PolynomialCommitment<F, P, S>,
    {
        let info = TestInfo {
            num_iters: 100,
            max_degree: None,
            supported_degree: None,
            num_vars,
            num_polynomials: 10,
            enforce_degree_bounds: true,
            max_num_queries: 5,
            num_equations: Some(10),
            rand_poly,
            rand_point,
            opening_challenge
        };
        equation_test_template::<F, P, PC, S>(info)
    }

    pub fn single_equation_test<F, P, PC, S>(
        num_vars: Option<usize>,
        rand_poly: fn(usize, Option<usize>, &mut StdRng) -> P,
        rand_point: fn(Option<usize>, &mut StdRng) -> P::Point,
        opening_challenge: fn() -> ChallengeGenerator<F,P>
    ) -> Result<(), PC::Error>
    where
        F: Field,
        P: Polynomial<F>,
        PC: PolynomialCommitment<F, P, S>,
    {
        let info = TestInfo {
            num_iters: 100,
            max_degree: None,
            supported_degree: None,
            num_vars,
            num_polynomials: 1,
            enforce_degree_bounds: false,
            max_num_queries: 1,
            num_equations: Some(1),
            rand_poly,
            rand_point,
            opening_challenge
        };
        equation_test_template::<F, P, PC, S>(info)
    }

    pub fn two_equation_test<F, P, PC, S>(
        num_vars: Option<usize>,
        rand_poly: fn(usize, Option<usize>, &mut StdRng) -> P,
        rand_point: fn(Option<usize>, &mut StdRng) -> P::Point,
        opening_challenge: fn() -> ChallengeGenerator<F,P>
    ) -> Result<(), PC::Error>
    where
        F: Field,
        P: Polynomial<F>,
        PC: PolynomialCommitment<F, P, S>,
    {
        let info = TestInfo {
            num_iters: 100,
            max_degree: None,
            supported_degree: None,
            num_vars,
            num_polynomials: 2,
            enforce_degree_bounds: false,
            max_num_queries: 1,
            num_equations: Some(2),
            rand_poly,
            rand_point,
            opening_challenge
        };
        equation_test_template::<F, P, PC, S>(info)
    }

    pub fn two_equation_degree_bound_test<F, P, PC, S>(
        rand_poly: fn(usize, Option<usize>, &mut StdRng) -> P,
        rand_point: fn(Option<usize>, &mut StdRng) -> P::Point,
        opening_challenge: fn() -> ChallengeGenerator<F,P>
    ) -> Result<(), PC::Error>
    where
        F: Field,
        P: Polynomial<F>,
        PC: PolynomialCommitment<F, P, S>,
    {
        let info = TestInfo {
            num_iters: 100,
            max_degree: None,
            supported_degree: None,
            num_vars: None,
            num_polynomials: 2,
            enforce_degree_bounds: true,
            max_num_queries: 1,
            num_equations: Some(2),
            rand_poly,
            rand_point,
            opening_challenge
        };
        equation_test_template::<F, P, PC, S>(info)
    }

    /// Generate default parameters (bls381-fr-only) for alpha = 17, state-size = 8
    pub(crate) fn poseidon_parameters_for_test<F: PrimeField>() -> PoseidonParameters<F> {
        let alpha = 17;
        let mds = vec![
            vec![
                F::from_str(
                    "43228725308391137369947362226390319299014033584574058394339561338097152657858",
                )
                    .map_err(|_| ())
                    .unwrap(),
                F::from_str(
                    "20729134655727743386784826341366384914431326428651109729494295849276339718592",
                )
                    .map_err(|_| ())
                    .unwrap(),
                F::from_str(
                    "14275792724825301816674509766636153429127896752891673527373812580216824074377",
                )
                    .map_err(|_| ())
                    .unwrap(),
            ],
            vec![
                F::from_str(
                    "3039440043015681380498693766234886011876841428799441709991632635031851609481",
                )
                    .map_err(|_| ())
                    .unwrap(),
                F::from_str(
                    "6678863357926068615342013496680930722082156498064457711885464611323928471101",
                )
                    .map_err(|_| ())
                    .unwrap(),
                F::from_str(
                    "37355038393562575053091209735467454314247378274125943833499651442997254948957",
                )
                    .map_err(|_| ())
                    .unwrap(),
            ],
            vec![
                F::from_str(
                    "26481612700543967643159862864328231943993263806649000633819754663276818191580",
                )
                    .map_err(|_| ())
                    .unwrap(),
                F::from_str(
                    "30103264397473155564098369644643015994024192377175707604277831692111219371047",
                )
                    .map_err(|_| ())
                    .unwrap(),
                F::from_str(
                    "5712721806190262694719203887224391960978962995663881615739647362444059585747",
                )
                    .map_err(|_| ())
                    .unwrap(),
            ],
        ];
        let ark = vec![
            vec![
                F::from_str(
                    "44595993092652566245296379427906271087754779418564084732265552598173323099784",
                )
                    .map_err(|_| ())
                    .unwrap(),
                F::from_str(
                    "23298463296221002559050231199021122673158929708101049474262017406235785365706",
                )
                    .map_err(|_| ())
                    .unwrap(),
                F::from_str(
                    "34212491019164671611180318500074499609633402631511849759183986060951187784466",
                )
                    .map_err(|_| ())
                    .unwrap(),
            ],
            vec![
                F::from_str(
                    "19098051134080182375553680073525644187968170656591203562523489333616681350367",
                )
                    .map_err(|_| ())
                    .unwrap(),
                F::from_str(
                    "7027675418691353855077049716619550622043312043660992344940177187528247727783",
                )
                    .map_err(|_| ())
                    .unwrap(),
                F::from_str(
                    "47642753235356257928619065424282314733361764347085604019867862722762702755609",
                )
                    .map_err(|_| ())
                    .unwrap(),
            ],
            vec![
                F::from_str(
                    "24281836129477728386327945482863886685457469794572168729834072693507088619997",
                )
                    .map_err(|_| ())
                    .unwrap(),
                F::from_str(
                    "12624893078331920791384400430193929292743809612452779381349824703573823883410",
                )
                    .map_err(|_| ())
                    .unwrap(),
                F::from_str(
                    "22654862987689323504199204643771547606936339944127455903448909090318619188561",
                )
                    .map_err(|_| ())
                    .unwrap(),
            ],
            vec![
                F::from_str(
                    "27229172992560143399715985732065737093562061782414043625359531774550940662372",
                )
                    .map_err(|_| ())
                    .unwrap(),
                F::from_str(
                    "13224952063922250960936823741448973692264041750100990569445192064567307041002",
                )
                    .map_err(|_| ())
                    .unwrap(),
                F::from_str(
                    "40380869235216625717296601204704413215735530626882135230693823362552484855508",
                )
                    .map_err(|_| ())
                    .unwrap(),
            ],
            vec![
                F::from_str(
                    "4245751157938905689397184705633683893932492370323323780371834663438472308145",
                )
                    .map_err(|_| ())
                    .unwrap(),
                F::from_str(
                    "8252156875535418429533049587170755750275631534314711502253775796882240991261",
                )
                    .map_err(|_| ())
                    .unwrap(),
                F::from_str(
                    "32910829712934971129644416249914075073083903821282503505466324428991624789936",
                )
                    .map_err(|_| ())
                    .unwrap(),
            ],
            vec![
                F::from_str(
                    "49412601297460128335642438246716127241669915737656789613664349252868389975962",
                )
                    .map_err(|_| ())
                    .unwrap(),
                F::from_str(
                    "841661305510340459373323516098909074520942972558284146843779636353111592117",
                )
                    .map_err(|_| ())
                    .unwrap(),
                F::from_str(
                    "37926489020263024391336570420006226544461516787280929232555625742588667303947",
                )
                    .map_err(|_| ())
                    .unwrap(),
            ],
            vec![
                F::from_str(
                    "18433043696013996573551852847056868761017170818820490351056924728720017242180",
                )
                    .map_err(|_| ())
                    .unwrap(),
                F::from_str(
                    "45376910275288438312773930242803223482318753992595269901397542214841496212310",
                )
                    .map_err(|_| ())
                    .unwrap(),
                F::from_str(
                    "47854349410014339708332226068958253098964727682486278458389508597930796651514",
                )
                    .map_err(|_| ())
                    .unwrap(),
            ],
            vec![
                F::from_str(
                    "32638426693771251366613055506166587312642876874690861030672730491779486904360",
                )
                    .map_err(|_| ())
                    .unwrap(),
                F::from_str(
                    "19105439281696418043426755774110765432959446684037017837894045255490581318047",
                )
                    .map_err(|_| ())
                    .unwrap(),
                F::from_str(
                    "13484299981373196201166722380389594773562113262309564134825386266765751213853",
                )
                    .map_err(|_| ())
                    .unwrap(),
            ],
            vec![
                F::from_str(
                    "63360321133852659797114062808297090090814531427710842859827725871241144161",
                )
                    .map_err(|_| ())
                    .unwrap(),
                F::from_str(
                    "42427543035537409467993338717379268954936885184662765745740070438835506287271",
                )
                    .map_err(|_| ())
                    .unwrap(),
                F::from_str(
                    "149101987103211771991327927827692640556911620408176100290586418839323044234",
                )
                    .map_err(|_| ())
                    .unwrap(),
            ],
            vec![
                F::from_str(
                    "8341764062226826803887898710015561861526081583071950015446833446251359696930",
                )
                    .map_err(|_| ())
                    .unwrap(),
                F::from_str(
                    "45635980415044299013530304465786867101223925975971912073759959440335364441441",
                )
                    .map_err(|_| ())
                    .unwrap(),
                F::from_str(
                    "49833261156201520743834327917353893365097424877680239796845398698940689734850",
                )
                    .map_err(|_| ())
                    .unwrap(),
            ],
            vec![
                F::from_str(
                    "26764715016591436228000634284249890185894507497739511725029482580508707525029",
                )
                    .map_err(|_| ())
                    .unwrap(),
                F::from_str(
                    "25054530812095491217523557726611612265064441619646263299990388543372685322499",
                )
                    .map_err(|_| ())
                    .unwrap(),
                F::from_str(
                    "47654590955096246997622155031169641628093104787883934397920286718814889326452",
                )
                    .map_err(|_| ())
                    .unwrap(),
            ],
            vec![
                F::from_str(
                    "16463825890556752307085325855351334996898686633642574805918056141310194135796",
                )
                    .map_err(|_| ())
                    .unwrap(),
                F::from_str(
                    "17473961341633494489168064889016732306117097771640351649096482400214968053040",
                )
                    .map_err(|_| ())
                    .unwrap(),
                F::from_str(
                    "49914603434867854893558366922996753035832008639512305549839666311012232077468",
                )
                    .map_err(|_| ())
                    .unwrap(),
            ],
            vec![
                F::from_str(
                    "17122578514152308432111470949473865420090463026624297565504381163777697818362",
                )
                    .map_err(|_| ())
                    .unwrap(),
                F::from_str(
                    "34870689836420861427379101859113225049736283485335674111421609473028315711541",
                )
                    .map_err(|_| ())
                    .unwrap(),
                F::from_str(
                    "4622082908476410083286670201138165773322781640914243047922441301693321472984",
                )
                    .map_err(|_| ())
                    .unwrap(),
            ],
            vec![
                F::from_str(
                    "6079244375752010013798561155333454682564824861645642293573415833483620500976",
                )
                    .map_err(|_| ())
                    .unwrap(),
                F::from_str(
                    "2635090520059500019661864086615522409798872905401305311748231832709078452746",
                )
                    .map_err(|_| ())
                    .unwrap(),
                F::from_str(
                    "19070766579582338321241892986615538320421651429118757507174186491084617237586",
                )
                    .map_err(|_| ())
                    .unwrap(),
            ],
            vec![
                F::from_str(
                    "12622420533971517050761060317049369208980632120901481436392835424625664738526",
                )
                    .map_err(|_| ())
                    .unwrap(),
                F::from_str(
                    "8965101225657199137904506150282256568170501907667138404080397024857524386266",
                )
                    .map_err(|_| ())
                    .unwrap(),
                F::from_str(
                    "27085091008069524593196374148553176565775450537072498305327481366756159319838",
                )
                    .map_err(|_| ())
                    .unwrap(),
            ],
            vec![
                F::from_str(
                    "45929056591150668409624595495643698205830429971690813312608217341940499221218",
                )
                    .map_err(|_| ())
                    .unwrap(),
                F::from_str(
                    "50361689160518167880500080025023064746137161030119436080957023803101861300846",
                )
                    .map_err(|_| ())
                    .unwrap(),
                F::from_str(
                    "6722586346537620732668048024627882970582133613352245923413730968378696371065",
                )
                    .map_err(|_| ())
                    .unwrap(),
            ],
            vec![
                F::from_str(
                    "7340485916200743279276570085958556798507770452421357119145466906520506506342",
                )
                    .map_err(|_| ())
                    .unwrap(),
                F::from_str(
                    "25946733168219652706630789514519162148860502996914241011500280690204368174083",
                )
                    .map_err(|_| ())
                    .unwrap(),
                F::from_str(
                    "9962367658743163006517635070396368828381757404628822422306438427554934645464",
                )
                    .map_err(|_| ())
                    .unwrap(),
            ],
            vec![
                F::from_str(
                    "7221669722700687417346373353960536661883467014204005276831020252277657076044",
                )
                    .map_err(|_| ())
                    .unwrap(),
                F::from_str(
                    "21487980358388383563030903293359140836304488103090321183948009095669344637431",
                )
                    .map_err(|_| ())
                    .unwrap(),
                F::from_str(
                    "44389482047246878765773958430749333249729101516826571588063797358040130313157",
                )
                    .map_err(|_| ())
                    .unwrap(),
            ],
            vec![
                F::from_str(
                    "32887270862917330820874162842519225370447850172085449103568878409533683733185",
                )
                    .map_err(|_| ())
                    .unwrap(),
                F::from_str(
                    "15453393396765207016379045014101989306173462885430532298601655955681532648226",
                )
                    .map_err(|_| ())
                    .unwrap(),
                F::from_str(
                    "5478929644476681096437469958231489102974161353940993351588559414552523375472",
                )
                    .map_err(|_| ())
                    .unwrap(),
            ],
            vec![
                F::from_str(
                    "41981370411247590312677561209178363054744730805951096631186178388981705304138",
                )
                    .map_err(|_| ())
                    .unwrap(),
                F::from_str(
                    "3474136981645476955784428843999869229067282976757744542648188369810577298585",
                )
                    .map_err(|_| ())
                    .unwrap(),
                F::from_str(
                    "26251477770740399889956219915654371915771248171098220204692699710414817081869",
                )
                    .map_err(|_| ())
                    .unwrap(),
            ],
            vec![
                F::from_str(
                    "51916561889718854106125837319509539220778634838409949714061033196765117231752",
                )
                    .map_err(|_| ())
                    .unwrap(),
                F::from_str(
                    "25355145802812435959748831835587713214179184608408449220418373832038339021974",
                )
                    .map_err(|_| ())
                    .unwrap(),
                F::from_str(
                    "31950684570730625275416731570246297947385359051792335826965013637877068017530",
                )
                    .map_err(|_| ())
                    .unwrap(),
            ],
            vec![
                F::from_str(
                    "40966378914980473680181850710703295982197782082391794594149984057481543436879",
                )
                    .map_err(|_| ())
                    .unwrap(),
                F::from_str(
                    "1141315130963422417761731263662398620858625339733452795772225916965481730059",
                )
                    .map_err(|_| ())
                    .unwrap(),
                F::from_str(
                    "9812100862165422922235757591915383485338044715409891361026651619010947646011",
                )
                    .map_err(|_| ())
                    .unwrap(),
            ],
            vec![
                F::from_str(
                    "25276091996614379065765602410190790163396484122487585763380676888280427744737",
                )
                    .map_err(|_| ())
                    .unwrap(),
                F::from_str(
                    "18512694312063606403196469408971540495273694846641903978723927656359350642619",
                )
                    .map_err(|_| ())
                    .unwrap(),
                F::from_str(
                    "5791584766415439694303685437881192048262049244830616851865505314899699012588",
                )
                    .map_err(|_| ())
                    .unwrap(),
            ],
            vec![
                F::from_str(
                    "34501536331706470927069149344450300773777486993504673779438188495686129846168",
                )
                    .map_err(|_| ())
                    .unwrap(),
                F::from_str(
                    "10797737565565774079718466476236831116206064650762676383469703413649447678207",
                )
                    .map_err(|_| ())
                    .unwrap(),
                F::from_str(
                    "42599392747310354323136214835734307933597896695637215127297036595538235868368",
                )
                    .map_err(|_| ())
                    .unwrap(),
            ],
            vec![
                F::from_str(
                    "1336670998775417133322626564820911986969949054454812685145275612519924150700",
                )
                    .map_err(|_| ())
                    .unwrap(),
                F::from_str(
                    "2630141283339761901081411552890260088516693208402906795133548756078952896770",
                )
                    .map_err(|_| ())
                    .unwrap(),
                F::from_str(
                    "5206688943117414740600380377278238268309952400341418217132724749372435975215",
                )
                    .map_err(|_| ())
                    .unwrap(),
            ],
            vec![
                F::from_str(
                    "10739264253827005683370721104077252560524362323422172665530191908848354339715",
                )
                    .map_err(|_| ())
                    .unwrap(),
                F::from_str(
                    "48010640624945719826344492755710886355389194986527731603685956726907395779674",
                )
                    .map_err(|_| ())
                    .unwrap(),
                F::from_str(
                    "47880724693177306044229143357252697148359033158394459365791331000715957339701",
                )
                    .map_err(|_| ())
                    .unwrap(),
            ],
            vec![
                F::from_str(
                    "51658938856669444737833983076793759752280196674149218924101718974926964118996",
                )
                    .map_err(|_| ())
                    .unwrap(),
                F::from_str(
                    "27558055650076329657496888512074319504342606463881203707330358472954748913263",
                )
                    .map_err(|_| ())
                    .unwrap(),
                F::from_str(
                    "38886981777859313701520424626728402175860609948757992393598285291689196608037",
                )
                    .map_err(|_| ())
                    .unwrap(),
            ],
            vec![
                F::from_str(
                    "17152756165118461969542990684402410297675979513690903033350206658079448802479",
                )
                    .map_err(|_| ())
                    .unwrap(),
                F::from_str(
                    "43766946932033687220387514221943418338304186408056458476301583041390483707207",
                )
                    .map_err(|_| ())
                    .unwrap(),
                F::from_str(
                    "24324495647041812436929170644873622904287038078113808264580396461953421400343",
                )
                    .map_err(|_| ())
                    .unwrap(),
            ],
            vec![
                F::from_str(
                    "6935839211798937659784055008131602708847374430164859822530563797964932598700",
                )
                    .map_err(|_| ())
                    .unwrap(),
                F::from_str(
                    "42126767398190942911395299419182514513368023621144776598842282267908712110039",
                )
                    .map_err(|_| ())
                    .unwrap(),
                F::from_str(
                    "5702364486091252903915715761606014714345316580946072019346660327857498603375",
                )
                    .map_err(|_| ())
                    .unwrap(),
            ],
            vec![
                F::from_str(
                    "28184981699552917714085740963279595942132561155181044254318202220270242523053",
                )
                    .map_err(|_| ())
                    .unwrap(),
                F::from_str(
                    "27078204494010940048327822707224393686245007379331357330801926151074766130790",
                )
                    .map_err(|_| ())
                    .unwrap(),
                F::from_str(
                    "5004172841233947987988267535285080365124079140142987718231874743202918551203",
                )
                    .map_err(|_| ())
                    .unwrap(),
            ],
            vec![
                F::from_str(
                    "7974360962120296064882769128577382489451060235999590492215336103105134345602",
                )
                    .map_err(|_| ())
                    .unwrap(),
                F::from_str(
                    "48062035869818179910046292951628308709251170031813126950740044942870578526376",
                )
                    .map_err(|_| ())
                    .unwrap(),
                F::from_str(
                    "26361151154829600651603985995297072258262605598910254660032612019129606811983",
                )
                    .map_err(|_| ())
                    .unwrap(),
            ],
            vec![
                F::from_str(
                    "46973867849986280770641828877435510444176572688208439836496241838832695841519",
                )
                    .map_err(|_| ())
                    .unwrap(),
                F::from_str(
                    "1219439673853113792340300173186247996249367102884530407862469123523013083971",
                )
                    .map_err(|_| ())
                    .unwrap(),
                F::from_str(
                    "8063356002935671186275773257019749639571745240775941450161086349727882957042",
                )
                    .map_err(|_| ())
                    .unwrap(),
            ],
            vec![
                F::from_str(
                    "8815571992701260640209942886673939234666734294275300852283020522390608544536",
                )
                    .map_err(|_| ())
                    .unwrap(),
                F::from_str(
                    "36384568984671043678320545346945893232044626942887414733675890845013312931948",
                )
                    .map_err(|_| ())
                    .unwrap(),
                F::from_str(
                    "7493936589040764830842760521372106574503511314427857201860148571929278344956",
                )
                    .map_err(|_| ())
                    .unwrap(),
            ],
            vec![
                F::from_str(
                    "26516538878265871822073279450474977673130300973488209984756372331392531193948",
                )
                    .map_err(|_| ())
                    .unwrap(),
                F::from_str(
                    "3872858659373466814413243601289105962248870842202907364656526273784217311104",
                )
                    .map_err(|_| ())
                    .unwrap(),
                F::from_str(
                    "8291822807524000248589997648893671538524566700364221355689839490238724479848",
                )
                    .map_err(|_| ())
                    .unwrap(),
            ],
            vec![
                F::from_str(
                    "32842548776827046388198955038089826231531188946525483251252938248379132381248",
                )
                    .map_err(|_| ())
                    .unwrap(),
                F::from_str(
                    "10749428410907700061565796335489079278748501945557710351216806276547834974736",
                )
                    .map_err(|_| ())
                    .unwrap(),
                F::from_str(
                    "43342287917341177925402357903832370099402579088513884654598017447701677948416",
                )
                    .map_err(|_| ())
                    .unwrap(),
            ],
            vec![
                F::from_str(
                    "29658571352070370791360499299098360881857072189358092237807807261478461425147",
                )
                    .map_err(|_| ())
                    .unwrap(),
                F::from_str(
                    "7805182565862454238315452208989152534554369855020544477885853141626690738363",
                )
                    .map_err(|_| ())
                    .unwrap(),
                F::from_str(
                    "30699555847500141715826240743138908521140760599479365867708690318477369178275",
                )
                    .map_err(|_| ())
                    .unwrap(),
            ],
            vec![
                F::from_str(
                    "1231951350103545216624376889222508148537733140742167414518514908719103925687",
                )
                    .map_err(|_| ())
                    .unwrap(),
                F::from_str(
                    "24784260089125933876714702247471508077514206350883487938806451152907502751770",
                )
                    .map_err(|_| ())
                    .unwrap(),
                F::from_str(
                    "36563542611079418454711392295126742705798573252480028863133394504154697924536",
                )
                    .map_err(|_| ())
                    .unwrap(),
            ],
        ];
        let full_rounds = 8;
        let total_rounds = 37;
        let partial_rounds = total_rounds - full_rounds;
        PoseidonParameters {
            full_rounds,
            partial_rounds,
            alpha,
            ark,
            mds,
        }
    }
}
