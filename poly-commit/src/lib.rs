#![cfg_attr(not(feature = "std"), no_std)]
//! A crate for polynomial commitment schemes.
#![deny(unused_import_braces, unused_qualifications, trivial_casts)]
#![deny(trivial_numeric_casts, variant_size_differences)]
#![deny(stable_features, unreachable_pub, non_shorthand_field_patterns)]
#![deny(unused_attributes, unused_mut)]
#![deny(missing_docs)]
#![deny(unused_imports)]
#![deny(renamed_and_removed_lints, stable_features, unused_allocation)]
#![deny(unused_comparisons, bare_trait_objects, unused_must_use)]
#![forbid(unsafe_code)]
#![doc = include_str!("../../README.md")]

#[allow(unused)]
#[macro_use]
extern crate derivative;
#[macro_use]
extern crate ark_std;

use ark_ff::{Field, PrimeField};
pub use ark_poly::{DenseUVPolynomial, Polynomial};
use ark_std::rand::RngCore;

use ark_std::{
    collections::{BTreeMap, BTreeSet},
    fmt::Debug,
    hash::Hash,
    iter::FromIterator,
    string::{String, ToString},
    vec::Vec,
};

/// Data structures used by a polynomial commitment scheme.
pub mod data_structures;
pub use data_structures::*;

/// Useful functions
pub(crate) mod utils;

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
#[cfg(all(test, not(feature = "std")))]
macro_rules! println {
    () => {};
    ($($arg: tt)*) => {};
}
/// The core [[KZG10]][kzg] construction.
///
/// [kzg]: http://cacr.uwaterloo.ca/techreports/2010/cacr2010-10.pdf
pub mod kzg10;

/// Polynomial commitment scheme from [[KZG10]][kzg] that enforces
/// strict degree bounds and (optionally) enables hiding commitments by
/// following the approach outlined in [[CHMMVW20, "Marlin"]][marlin].
///
/// [kzg]: http://cacr.uwaterloo.ca/techreports/2010/cacr2010-10.pdf
/// [marlin]: https://eprint.iacr.org/2019/1047
pub use marlin::marlin_pc;

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
pub mod sonic_pc;

/// A polynomial commitment scheme based on the hardness of the
/// discrete logarithm problem in prime-order groups.
/// The construction is detailed in [[BCMS20]][pcdas].
///
/// [pcdas]: https://eprint.iacr.org/2020/499
pub mod ipa_pc;

/// A multilinear polynomial commitment scheme that converts n-variate multilinear polynomial into
/// n quotient UV polynomial. This scheme is based on hardness of the discrete logarithm
/// in prime-order groups. Construction is detailed in [[XZZPD19]][xzzpd19] and [[ZGKPP18]][zgkpp18]
///
/// [xzzpd19]: https://eprint.iacr.org/2019/317
/// [zgkpp]: https://ieeexplore.ieee.org/document/8418645
pub mod multilinear_pc;

use ark_crypto_primitives::sponge::{CryptographicSponge, FieldElementSize};
/// Multivariate polynomial commitment based on the construction in
/// [[PST13]][pst] with batching and (optional) hiding property inspired
/// by the univariate scheme in [[CHMMVW20, "Marlin"]][marlin]
///
/// [pst]: https://eprint.iacr.org/2011/587.pdf
/// [marlin]: https://eprint.iacr.org/2019/1047
pub use marlin::marlin_pst13_pc;

/// Streaming polynomial commitment based on the construction in
/// [[BCHO22, "Gemini"]][gemini] with batching techniques inspired
/// by [[BDFG20]][bdfg].
///
/// [gemini]:
/// [bdfg]: https://eprint.iacr.org/2020/081.pdf
pub mod streaming_kzg;

/// A polynomial commitment scheme based on the hardness of the
/// discrete logarithm problem in prime-order groups. This is a
/// Fiat-Shamired version of the PCS described in the Hyrax paper
/// [[WTsTW17]][hyrax], with the difference that, unlike in the
/// cited reference, the evaluation of the polynomial at the point
/// of interest is indeed revealed to the verifier at the end.
///
/// [hyrax]: https://eprint.iacr.org/2017/1132.pdf
pub mod hyrax;

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
pub trait PolynomialCommitment<F: PrimeField, P: Polynomial<F>, S: CryptographicSponge>:
    Sized
{
    /// The universal parameters for the commitment scheme. These are "trimmed"
    /// down to `Self::CommitterKey` and `Self::VerifierKey` by `Self::trim`.
    type UniversalParams: PCUniversalParams;
    /// The committer key for the scheme; used to commit to a polynomial and then
    /// open the commitment to produce an evaluation proof.
    type CommitterKey: PCCommitterKey;
    /// The verifier key for the scheme; used to check an evaluation proof.
    type VerifierKey: PCVerifierKey;
    /// The commitment to a polynomial.
    type Commitment: PCCommitment + Default;
    /// Auxiliary state of the commitment, output by the `commit` phase.
    /// It contains information that can be reused by the committer
    /// during the `open` phase, such as the commitment randomness.
    /// Not to be shared with the verifier.
    type CommitmentState: PCCommitmentState;
    /// The evaluation proof for a single point.
    type Proof: Clone;
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

    /// Outputs a list of commitments to `polynomials`. If `polynomials[i].is_hiding()`,
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
            Vec<Self::CommitmentState>,
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
        sponge: &mut S,
        states: impl IntoIterator<Item = &'a Self::CommitmentState>,
        rng: Option<&mut dyn RngCore>,
    ) -> Result<Self::Proof, Self::Error>
    where
        P: 'a,
        Self::CommitmentState: 'a,
        Self::Commitment: 'a;

    /// check but with individual challenges
    fn check<'a>(
        vk: &Self::VerifierKey,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Self::Commitment>>,
        point: &'a P::Point,
        values: impl IntoIterator<Item = F>,
        proof: &Self::Proof,
        sponge: &mut S,
        rng: Option<&mut dyn RngCore>,
    ) -> Result<bool, Self::Error>
    where
        Self::Commitment: 'a;

    /// Open several polynomials at one or more points each (possibly different
    /// for each polynomial). Each entry in the in the query set of points
    /// contains the label of the polynomial which should be queried at that
    /// point.
    ///
    /// Behaviour is undefined if `query_set` contains the entries with the
    /// same point label but different actual points.
    ///
    /// The opening challenges are independent for each batch of polynomials.
    fn batch_open<'a>(
        ck: &Self::CommitterKey,
        labeled_polynomials: impl IntoIterator<Item = &'a LabeledPolynomial<F, P>>,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Self::Commitment>>,
        query_set: &QuerySet<P::Point>,
        sponge: &mut S,
        states: impl IntoIterator<Item = &'a Self::CommitmentState>,
        rng: Option<&mut dyn RngCore>,
    ) -> Result<Self::BatchProof, Self::Error>
    where
        P: 'a,
        Self::CommitmentState: 'a,
        Self::Commitment: 'a,
    {
        // The default implementation achieves proceeds by rearranging the queries in
        // order to gather (i.e. batch) all polynomials that should be queried at
        // the same point, then opening their commitments simultaneously with a
        // single call to `open` (per point)
        let rng = &mut crate::optional_rng::OptionalRng(rng);
        let poly_st_comm: BTreeMap<_, _> = labeled_polynomials
            .into_iter()
            .zip(states)
            .zip(commitments.into_iter())
            .map(|((poly, st), comm)| (poly.label(), (poly, st, comm)))
            .collect();

        let open_time = start_timer!(|| format!(
            "Opening {} polynomials at query set of size {}",
            poly_st_comm.len(),
            query_set.len(),
        ));

        let mut query_to_labels_map = BTreeMap::new();

        // `label` is the label of the polynomial the query refers to
        // `point_label` is the label of the point being queried
        // `point` is the actual point
        for (label, (point_label, point)) in query_set.iter() {
            // For each point label in `query_set`, we define an entry in
            // `query_to_labels_map` containing a pair whose first element is
            // the actual point and the second one is the set of labels of the
            // polynomials being queried at that point
            let labels = query_to_labels_map
                .entry(point_label)
                .or_insert((point, BTreeSet::new()));
            labels.1.insert(label);
        }

        let mut proofs = Vec::new();
        for (_point_label, (point, labels)) in query_to_labels_map.into_iter() {
            let mut query_polys: Vec<&'a LabeledPolynomial<_, _>> = Vec::new();
            let mut query_states: Vec<&'a Self::CommitmentState> = Vec::new();
            let mut query_comms: Vec<&'a LabeledCommitment<Self::Commitment>> = Vec::new();

            // Constructing matching vectors with the polynomial, commitment
            // randomness and actual commitment for each polynomial being
            // queried at `point`
            for label in labels {
                let (polynomial, state, comm) =
                    poly_st_comm.get(label).ok_or(Error::MissingPolynomial {
                        label: label.to_string(),
                    })?;

                query_polys.push(polynomial);
                query_states.push(state);
                query_comms.push(comm);
            }

            let proof_time = start_timer!(|| "Creating proof");

            // Simultaneously opening the commitments of all polynomials that
            // refer to the the current point using the plain `open` function
            let proof = Self::open(
                ck,
                query_polys,
                query_comms,
                &point,
                sponge,
                query_states,
                Some(rng),
            )?;

            end_timer!(proof_time);

            proofs.push(proof);
        }
        end_timer!(open_time);

        Ok(proofs.into())
    }

    /// Verify opening proofs for several polynomials at one or more points
    /// each (possibly different for each polynomial). Each entry in
    /// the query set of points contains the label of the polynomial which
    /// was queried at that point.
    ///
    /// Behaviour is undefined if `query_set` contains the entries with the
    /// same point label but different points.
    ///
    /// Behaviour is also undefined if proofs are not ordered the same way as
    /// queries in `query_to_labels_map` (this is the outcome of calling
    /// `batch_open` for the same commitment list and query set).H
    ///
    /// The opening challenges are independent for each batch of polynomials.
    fn batch_check<'a, R: RngCore>(
        vk: &Self::VerifierKey,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Self::Commitment>>,
        query_set: &QuerySet<P::Point>,
        evaluations: &Evaluations<P::Point, F>,
        proof: &Self::BatchProof,
        sponge: &mut S,
        rng: &mut R,
    ) -> Result<bool, Self::Error>
    where
        Self::Commitment: 'a,
    {
        // The default implementation proceeds by rearranging the queries in
        // order to gather (i.e. batch) the proofs of all polynomials that should
        // have been opened at the same point, then verifying those proofs
        // simultaneously with a single call to `check` (per point).
        let commitments: BTreeMap<_, _> = commitments.into_iter().map(|c| (c.label(), c)).collect();
        let mut query_to_labels_map = BTreeMap::new();

        // `label` is the label of the polynomial the query refers to
        // `point_label` is the label of the point being queried
        // `point` is the actual point
        for (label, (point_label, point)) in query_set.iter() {
            // For each point label in `query_set`, we define an entry in
            // `query_to_labels_map` containing a pair whose first element is
            // the actual point and the second one is the set of labels of the
            // polynomials being queried at that point
            let labels = query_to_labels_map
                .entry(point_label)
                .or_insert((point, BTreeSet::new()));
            labels.1.insert(label);
        }

        // Implicit assumption: proofs are ordered in same manner as queries in
        // `query_to_labels_map`.
        let proofs: Vec<_> = proof.clone().into();
        assert_eq!(proofs.len(), query_to_labels_map.len());

        let mut result = true;
        for ((_point_label, (point, labels)), proof) in query_to_labels_map.into_iter().zip(proofs)
        {
            // Constructing matching vectors with the commitment and claimed
            // value of each polynomial being queried at `point`
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

            // Verify all proofs referring to the current point simultaneously
            // with a single call to `check`
            result &= Self::check(vk, comms, &point, values, &proof, sponge, Some(rng))?;
            end_timer!(proof_time);
        }
        Ok(result)
    }

    /// Open commitments to all polynomials involved in a number of linear
    /// combinations (LC) simultaneously.
    fn open_combinations<'a>(
        ck: &Self::CommitterKey,
        linear_combinations: impl IntoIterator<Item = &'a LinearCombination<F>>,
        polynomials: impl IntoIterator<Item = &'a LabeledPolynomial<F, P>>,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Self::Commitment>>,
        query_set: &QuerySet<P::Point>,
        sponge: &mut S,
        states: impl IntoIterator<Item = &'a Self::CommitmentState>,
        rng: Option<&mut dyn RngCore>,
    ) -> Result<BatchLCProof<F, Self::BatchProof>, Self::Error>
    where
        Self::CommitmentState: 'a,
        Self::Commitment: 'a,
        P: 'a,
    {
        // The default implementation proceeds by batch-opening all polynomials
        // appearing in those LC that are queried at the same point.
        let linear_combinations: Vec<_> = linear_combinations.into_iter().collect();
        let polynomials: Vec<_> = polynomials.into_iter().collect();

        // Rearrange the information about queries on linear combinations into
        // information about queries on individual polynomials.
        let poly_query_set =
            lc_query_set_to_poly_query_set(linear_combinations.iter().copied(), query_set);
        let poly_evals = evaluate_query_set(polynomials.iter().copied(), &poly_query_set);

        // Batch-open all polynomials that refer to each individual point in `query_set`
        let proof = Self::batch_open(
            ck,
            polynomials,
            commitments,
            &poly_query_set,
            sponge,
            states,
            rng,
        )?;
        Ok(BatchLCProof {
            proof,
            evals: Some(poly_evals.values().copied().collect()),
        })
    }

    /// Verify opening proofs for all polynomials involved in a number of
    /// linear combinations (LC) simultaneously.
    fn check_combinations<'a, R: RngCore>(
        vk: &Self::VerifierKey,
        linear_combinations: impl IntoIterator<Item = &'a LinearCombination<F>>,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Self::Commitment>>,
        eqn_query_set: &QuerySet<P::Point>,
        eqn_evaluations: &Evaluations<P::Point, F>,
        proof: &BatchLCProof<F, Self::BatchProof>,
        sponge: &mut S,
        rng: &mut R,
    ) -> Result<bool, Self::Error>
    where
        Self::Commitment: 'a,
    {
        // The default implementation does this by batch-checking each
        // batch-opening proof of polynomials appearing in those LC that were
        // queried at the same point, then computing the evaluations of each LC
        // using the proved polynomial evaluations.
        let BatchLCProof { proof, evals } = proof;
        let lc_s = BTreeMap::from_iter(linear_combinations.into_iter().map(|lc| (lc.label(), lc)));

        // Rearrange the information about queries on linear combinations into
        // information about queries on individual polynomials.
        let poly_query_set = lc_query_set_to_poly_query_set(lc_s.values().copied(), eqn_query_set);
        let sorted_by_poly_and_query_label: BTreeSet<_> = poly_query_set
            .clone()
            .into_iter()
            .map(|(poly_label, v)| ((poly_label.clone(), v.1), v.0))
            .collect();

        let poly_evals = Evaluations::from_iter(
            sorted_by_poly_and_query_label
                .into_iter()
                .zip(evals.clone().unwrap())
                .map(|(((poly_label, point), _query_label), eval)| ((poly_label, point), eval)),
        );

        for &(ref lc_label, (_, ref point)) in eqn_query_set {
            if let Some(lc) = lc_s.get(lc_label) {
                let claimed_rhs = *eqn_evaluations
                    .get(&(lc_label.clone(), point.clone()))
                    .ok_or(Error::MissingEvaluation {
                        label: lc_label.to_string(),
                    })?;

                let mut actual_rhs = F::zero();

                // Compute the value of the linear combination by adding the
                // claimed value for each polynomial in it (to be proved later)
                // scaled by the corresponding coefficient.
                for (coeff, label) in lc.iter() {
                    let eval = match label {
                        LCTerm::One => F::one(),
                        LCTerm::PolyLabel(l) => *poly_evals
                            .get(&(l.clone().into(), point.clone()))
                            .ok_or(Error::MissingEvaluation {
                                label: format!("{}-{:?}", l.clone(), point.clone()),
                            })?,
                    };

                    actual_rhs += &(*coeff * eval);
                }

                // Checking the computed evaluation matches the claimed one
                if claimed_rhs != actual_rhs {
                    eprintln!("Claimed evaluation of {} is incorrect", lc.label());
                    return Ok(false);
                }
            }
        }

        // Verify the claimed evaluation for each polynomial appearing in the
        // linear combinations, batched by point
        let pc_result = Self::batch_check(
            vk,
            commitments,
            &poly_query_set,
            &poly_evals,
            proof,
            sponge,
            rng,
        )?;
        if !pc_result {
            eprintln!("Evaluation proofs failed to verify");
            return Ok(false);
        }

        Ok(true)
    }
}

/// The size of opening challenges in bits.
pub const CHALLENGE_SIZE: FieldElementSize = FieldElementSize::Truncated(128);

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

// Separate the information about queries on linear combinations into
// information about queries on individual polynomials.
//
// For instance, if `linear_combinations` is
// [
//  ("average", 1/2 * pol_1 + 1/2 * pol_2),
//  ("weighted", 1/2 * pol_1 + 1/2 * pol_3)
// ]
// and `query_set` is
// [
//  ("average", ("three", 3))
//  ("weighted", ("three", 3))
// ],
// then the output is
// {
//  ("pol_1", ("three", 3)),
//  ("pol_2", ("three", 3)),
//  ("pol_3", ("three", 3)),
// }
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
    use ark_crypto_primitives::sponge::poseidon::{PoseidonConfig, PoseidonSponge};
    use ark_poly::Polynomial;
    use ark_std::rand::{
        distributions::{Distribution, Uniform},
        Rng, SeedableRng,
    };
    use ark_std::test_rng;
    use rand_chacha::ChaCha20Rng;

    struct TestInfo<F: PrimeField, P: Polynomial<F>, S: CryptographicSponge> {
        num_iters: usize,
        max_degree: Option<usize>,
        supported_degree: Option<usize>,
        num_vars: Option<usize>,
        num_polynomials: usize,
        enforce_degree_bounds: bool,
        max_num_queries: usize,
        num_equations: Option<usize>,
        rand_poly: fn(usize, Option<usize>, &mut ChaCha20Rng) -> P,
        rand_point: fn(Option<usize>, &mut ChaCha20Rng) -> P::Point,
        sponge: fn() -> S,
    }

    pub fn bad_degree_bound_test<F, P, PC, S>(
        rand_poly: fn(usize, Option<usize>, &mut ChaCha20Rng) -> P,
        rand_point: fn(Option<usize>, &mut ChaCha20Rng) -> P::Point,
        sponge: fn() -> S,
    ) -> Result<(), PC::Error>
    where
        F: PrimeField,
        P: Polynomial<F>,
        PC: PolynomialCommitment<F, P, S>,
        S: CryptographicSponge,
    {
        let sponge = sponge();

        let rng = &mut ChaCha20Rng::from_rng(test_rng()).unwrap();
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
                &mut (sponge.clone()),
                &rands,
                Some(rng),
            )?;
            let result = PC::batch_check(
                &vk,
                &comms,
                &query_set,
                &values,
                &proof,
                &mut (sponge.clone()),
                rng,
            )?;
            assert!(result, "proof was incorrect, Query set: {:#?}", query_set);
        }

        Ok(())
    }

    fn test_template<F, P, PC, S>(info: TestInfo<F, P, S>) -> Result<(), PC::Error>
    where
        F: PrimeField,
        P: Polynomial<F>,
        PC: PolynomialCommitment<F, P, S>,
        S: CryptographicSponge,
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
            sponge,
        } = info;

        let sponge = sponge();

        let rng = &mut ChaCha20Rng::from_rng(test_rng()).unwrap();
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
                &mut (sponge.clone()),
                &rands,
                Some(rng),
            )?;
            let result = PC::batch_check(
                &vk,
                &comms,
                &query_set,
                &values,
                &proof,
                &mut (sponge.clone()),
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

    fn equation_test_template<F, P, PC, S>(info: TestInfo<F, P, S>) -> Result<(), PC::Error>
    where
        F: PrimeField,
        P: Polynomial<F>,
        PC: PolynomialCommitment<F, P, S>,
        S: CryptographicSponge,
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
            sponge,
        } = info;

        let sponge = sponge();

        let rng = &mut ChaCha20Rng::from_rng(test_rng()).unwrap();
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
                &mut (sponge.clone()),
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
                &mut (sponge.clone()),
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
        rand_poly: fn(usize, Option<usize>, &mut ChaCha20Rng) -> P,
        rand_point: fn(Option<usize>, &mut ChaCha20Rng) -> P::Point,
        sponge: fn() -> S,
    ) -> Result<(), PC::Error>
    where
        F: PrimeField,
        P: Polynomial<F>,
        PC: PolynomialCommitment<F, P, S>,
        S: CryptographicSponge,
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
            sponge,
        };
        test_template::<F, P, PC, S>(info)
    }

    pub fn linear_poly_degree_bound_test<F, P, PC, S>(
        rand_poly: fn(usize, Option<usize>, &mut ChaCha20Rng) -> P,
        rand_point: fn(Option<usize>, &mut ChaCha20Rng) -> P::Point,
        sponge: fn() -> S,
    ) -> Result<(), PC::Error>
    where
        F: PrimeField,
        P: Polynomial<F>,
        PC: PolynomialCommitment<F, P, S>,
        S: CryptographicSponge,
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
            sponge,
        };
        test_template::<F, P, PC, S>(info)
    }

    pub fn single_poly_degree_bound_test<F, P, PC, S>(
        rand_poly: fn(usize, Option<usize>, &mut ChaCha20Rng) -> P,
        rand_point: fn(Option<usize>, &mut ChaCha20Rng) -> P::Point,
        sponge: fn() -> S,
    ) -> Result<(), PC::Error>
    where
        F: PrimeField,
        P: Polynomial<F>,
        PC: PolynomialCommitment<F, P, S>,
        S: CryptographicSponge,
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
            sponge,
        };
        test_template::<F, P, PC, S>(info)
    }

    pub fn quadratic_poly_degree_bound_multiple_queries_test<F, P, PC, S>(
        rand_poly: fn(usize, Option<usize>, &mut ChaCha20Rng) -> P,
        rand_point: fn(Option<usize>, &mut ChaCha20Rng) -> P::Point,
        sponge: fn() -> S,
    ) -> Result<(), PC::Error>
    where
        F: PrimeField,
        P: Polynomial<F>,
        PC: PolynomialCommitment<F, P, S>,
        S: CryptographicSponge,
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
            sponge,
        };
        test_template::<F, P, PC, S>(info)
    }

    pub fn single_poly_degree_bound_multiple_queries_test<F, P, PC, S>(
        rand_poly: fn(usize, Option<usize>, &mut ChaCha20Rng) -> P,
        rand_point: fn(Option<usize>, &mut ChaCha20Rng) -> P::Point,
        sponge: fn() -> S,
    ) -> Result<(), PC::Error>
    where
        F: PrimeField,
        P: Polynomial<F>,
        PC: PolynomialCommitment<F, P, S>,
        S: CryptographicSponge,
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
            sponge,
        };
        test_template::<F, P, PC, S>(info)
    }

    pub fn two_polys_degree_bound_single_query_test<F, P, PC, S>(
        rand_poly: fn(usize, Option<usize>, &mut ChaCha20Rng) -> P,
        rand_point: fn(Option<usize>, &mut ChaCha20Rng) -> P::Point,
        sponge: fn() -> S,
    ) -> Result<(), PC::Error>
    where
        F: PrimeField,
        P: Polynomial<F>,
        PC: PolynomialCommitment<F, P, S>,
        S: CryptographicSponge,
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
            sponge,
        };
        test_template::<F, P, PC, S>(info)
    }

    pub fn full_end_to_end_test<F, P, PC, S>(
        num_vars: Option<usize>,
        rand_poly: fn(usize, Option<usize>, &mut ChaCha20Rng) -> P,
        rand_point: fn(Option<usize>, &mut ChaCha20Rng) -> P::Point,
        sponge: fn() -> S,
    ) -> Result<(), PC::Error>
    where
        F: PrimeField,
        P: Polynomial<F>,
        PC: PolynomialCommitment<F, P, S>,
        S: CryptographicSponge,
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
            sponge,
        };
        test_template::<F, P, PC, S>(info)
    }

    pub fn full_end_to_end_equation_test<F, P, PC, S>(
        num_vars: Option<usize>,
        rand_poly: fn(usize, Option<usize>, &mut ChaCha20Rng) -> P,
        rand_point: fn(Option<usize>, &mut ChaCha20Rng) -> P::Point,
        sponge: fn() -> S,
    ) -> Result<(), PC::Error>
    where
        F: PrimeField,
        P: Polynomial<F>,
        PC: PolynomialCommitment<F, P, S>,
        S: CryptographicSponge,
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
            sponge,
        };
        equation_test_template::<F, P, PC, S>(info)
    }

    pub fn single_equation_test<F, P, PC, S>(
        num_vars: Option<usize>,
        rand_poly: fn(usize, Option<usize>, &mut ChaCha20Rng) -> P,
        rand_point: fn(Option<usize>, &mut ChaCha20Rng) -> P::Point,
        sponge: fn() -> S,
    ) -> Result<(), PC::Error>
    where
        F: PrimeField,
        P: Polynomial<F>,
        PC: PolynomialCommitment<F, P, S>,
        S: CryptographicSponge,
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
            sponge,
        };
        equation_test_template::<F, P, PC, S>(info)
    }

    pub fn two_equation_test<F, P, PC, S>(
        num_vars: Option<usize>,
        rand_poly: fn(usize, Option<usize>, &mut ChaCha20Rng) -> P,
        rand_point: fn(Option<usize>, &mut ChaCha20Rng) -> P::Point,
        sponge: fn() -> S,
    ) -> Result<(), PC::Error>
    where
        F: PrimeField,
        P: Polynomial<F>,
        PC: PolynomialCommitment<F, P, S>,
        S: CryptographicSponge,
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
            sponge,
        };
        equation_test_template::<F, P, PC, S>(info)
    }

    pub fn two_equation_degree_bound_test<F, P, PC, S>(
        rand_poly: fn(usize, Option<usize>, &mut ChaCha20Rng) -> P,
        rand_point: fn(Option<usize>, &mut ChaCha20Rng) -> P::Point,
        sponge: fn() -> S,
    ) -> Result<(), PC::Error>
    where
        F: PrimeField,
        P: Polynomial<F>,
        PC: PolynomialCommitment<F, P, S>,
        S: CryptographicSponge,
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
            sponge,
        };
        equation_test_template::<F, P, PC, S>(info)
    }

    pub(crate) fn poseidon_sponge_for_test<F: PrimeField>() -> PoseidonSponge<F> {
        PoseidonSponge::new(&poseidon_parameters_for_test())
    }

    /// Generate default parameters for alpha = 17, state-size = 8
    ///
    /// WARNING: This poseidon parameter is not secure. Please generate
    /// your own parameters according the field you use.
    pub(crate) fn poseidon_parameters_for_test<F: PrimeField>() -> PoseidonConfig<F> {
        let full_rounds = 8;
        let partial_rounds = 31;
        let alpha = 17;

        let mds = vec![
            vec![F::one(), F::zero(), F::one()],
            vec![F::one(), F::one(), F::zero()],
            vec![F::zero(), F::one(), F::one()],
        ];

        let mut ark = Vec::new();
        let mut ark_rng = test_rng();

        for _ in 0..(full_rounds + partial_rounds) {
            let mut res = Vec::new();

            for _ in 0..3 {
                res.push(F::rand(&mut ark_rng));
            }
            ark.push(res);
        }
        PoseidonConfig::new(full_rounds, partial_rounds, alpha, mds, ark, 2, 1)
    }
}
