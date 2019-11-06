//! A crate for polynomial commitment schemes.
#![deny(unused_import_braces, unused_qualifications, trivial_casts)]
#![deny(trivial_numeric_casts, private_in_public, variant_size_differences)]
#![deny(stable_features, unreachable_pub, non_shorthand_field_patterns)]
#![deny(unused_attributes, unused_imports, unused_mut, missing_docs)]
#![deny(renamed_and_removed_lints, stable_features, unused_allocation)]
#![deny(unused_comparisons, bare_trait_objects, unused_must_use, const_err)]
#![forbid(unsafe_code)]

#[macro_use]
extern crate algebra;
#[macro_use]
extern crate derivative;
#[macro_use]
extern crate bench_utils;

use algebra::Field;
pub use ff_fft::DensePolynomial as Polynomial;
use rand_core::RngCore;
use std::collections::{BTreeMap, BTreeSet};


/// Data structures used by a polynomial commitment scheme.
pub mod data_structures;
pub use data_structures::*;

/// The core KZG10 construction.
pub mod kzg10;

/// An adaptation of the KZG10 construction that uses the method outlined in 
/// the [Marlin paper][marlin] to enforce degree bounds.
///
/// [marlin]: https://eprint.iacr.org/2019/1047
pub mod marlin_kzg10;

/// `QuerySet` is the set of queries that are to be made to a set of labeled polynomials
/// `p` that have previously been committed. Each element of a `QuerySet` is a `(label, query)`
/// pair, where `label` is the label of a polynomial in `p`, and `query` is the field element
/// that `p[label]` is to be queried at.
pub type QuerySet<'a, F> = BTreeSet<(&'a str, F)>;

/// `Evaluations` is the result of querying a set of labeled polynomials `p` at a `QuerySet`
/// `Q`. It maps each element of `Q` to the resulting evaluation. That is,
/// if `(label, query)` is an element of `Q`, then `evaluation.get((label, query))`
/// should equal `p[label].evaluate(query)`.
pub type Evaluations<'a, F> = BTreeMap<(&'a str, F), F>;

/// Describes the interface for a polynomial commitment scheme that allows
/// a sender to commit to multiple polynomials and later provide a succinct proof
/// of evaluation for the corresponding commitments at a query set `Q`, while
/// enforcing per-polynomial degree bounds.
pub trait PolynomialCommitment<F: Field> {
    /// The universal parameters for the commitment scheme. These are "trimmed"
    /// down to `Self::CommitterKey` and `Self::VerifierKey` by `Self::trim`.
    type UniversalParams: PCUniversalParams;
    /// The committer key for the scheme; used to commit to a polynomial and then
    /// open the commitment to produce an evaluation proof.
    type CommitterKey: PCCommitterKey;
    /// The verifier key for the scheme; used to check an evaluation proof.
    type VerifierKey: PCVerifierKey;
    /// The commitment to a polynomial.
    type Commitment: PCCommitment;
    /// The commitment randomness.
    type Randomness: PCRandomness;
    /// The evaluation proof for a single point.
    type Proof: Clone;
    /// The evaluation proof for a query set.
    type BatchProof: Clone + From<Vec<Self::Proof>>;
    /// The error type for the scheme.
    type Error: std::error::Error;

    /// Constructs public parameters when given as input the maximum degree `degree`
    /// for the polynomial commitment scheme.
    fn setup<R: RngCore>(
        degree: usize,
        rng: &mut R,
    ) -> Result<Self::UniversalParams, Self::Error>;

    /// Specializes the public parameters for a given support for coefficients
    /// for polynomials and degree bounds.
    fn trim(
        pp: &Self::UniversalParams,
        max_degree: usize,
        degree_bounds_to_support: Option<&[usize]>,
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
        polynomials: impl IntoIterator<Item = &'a LabeledPolynomial<'a, F>>,
        rng: Option<&mut dyn RngCore>,
    ) -> Result<
        (
            Vec<LabeledCommitment<Self::Commitment>>,
            Vec<Self::Randomness>,
        ),
        Self::Error,
    >;

    /// On input a list of labeled polynomials and a query pont, `open` outputs a proof of evaluation
    /// of the polynomials at the query point.
    fn open<'a>(
        ck: &Self::CommitterKey,
        labeled_polynomials: impl IntoIterator<Item = &'a LabeledPolynomial<'a, F>>,
        point: F,
        opening_challenge: F,
        rands: impl IntoIterator<Item = &'a Self::Randomness>,
    ) -> Result<Self::Proof, Self::Error>
        where 
            Self::Randomness: 'a;


    /// On input a list of labeled polynomials and a query set, `open` outputs a proof of evaluation
    /// of the polynomials at the points in the query set.
    fn batch_open<'a>(
        ck: &Self::CommitterKey,
        labeled_polynomials: impl IntoIterator<Item = &'a LabeledPolynomial<'a, F>>,
        query_set: &QuerySet<F>,
        opening_challenge: F,
        rands: impl IntoIterator<Item = &'a Self::Randomness>,
    ) -> Result<Self::BatchProof, Self::Error> 
    where 
        Self::Randomness: 'a
    {
        let polynomials_with_rands: BTreeMap<_, _> = labeled_polynomials
            .into_iter()
            .zip(rands)
            .map(|(poly, r)| (poly.label(), (poly, r)))
            .collect();

        let open_time = start_timer!(|| format!(
            "Opening {} polynomials at query set of size {}",
            polynomials_with_rands.len(),
            query_set.len(),
        ));

        let mut query_to_labels_map = BTreeMap::new();

        for (label, point) in query_set.iter() {
            let labels = query_to_labels_map.entry(point).or_insert(BTreeSet::new());
            labels.insert(label);
        }

        let mut proofs = Vec::new();
        for (query, labels) in query_to_labels_map.into_iter() {
            let mut query_polys: Vec<&'a LabeledPolynomial<'a, _>> = Vec::new();
            let mut query_rands: Vec<&'a Self::Randomness> = Vec::new();
            for label in labels {
                let (polynomial, rand) = polynomials_with_rands.get(label).expect(
                    "query set references polynomial with incorrect label"
                );
                query_polys.push(polynomial);
                query_rands.push(rand);
            }
            let proof_time = start_timer!(|| "Creating proof");
            let proof = Self::open(ck, query_polys, *query, opening_challenge, query_rands)?;
            end_timer!(proof_time);

            proofs.push(proof);
        }
        end_timer!(open_time);

        Ok(proofs.into())
    }

    /// Verifies that `values` are the evaluations at `point` of the polynomials
    /// committed inside `commitments`.
    fn check<'a>(
        vk: &Self::VerifierKey,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Self::Commitment>>,
        point: F,
        values: impl IntoIterator<Item = F>,
        proof: &Self::Proof,
        opening_challenge: F,
    ) -> Result<bool, Self::Error>
    where
        Self::Commitment: 'a;

    /// Checks that `values` are the true evaluations at `query_set` of the polynomials
    /// committed in `labeled_commitments`.
    fn batch_check<'a, R: RngCore>(
        vk: &Self::VerifierKey,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Self::Commitment>>,
        query_set: &QuerySet<F>,
        values: &Evaluations<F>,
        proof: &Self::BatchProof,
        opening_challenge: F,
        rng: &mut R,
    ) -> Result<bool, Self::Error>
    where
        Self::Commitment: 'a;
}

#[cfg(test)]
pub mod tests {
    use crate::*;
    use algebra::Field;
    use algebra::UniformRand;
    use rand::{distributions::Distribution, thread_rng};

    pub fn single_poly_test<F, PC>() -> Result<(), PC::Error>
    where
        F: Field,
        PC: PolynomialCommitment<F>,
    {
        let rng = &mut thread_rng();
        let max_degree = rand::distributions::Uniform::from(2..64).sample(rng);
        let (ck, vk) = PC::setup(max_degree, rng)?;
        for _ in 0..100 {
            let mut polynomials = Vec::new();
            // Generate polynomials
            let num_points_in_query_set = 1;
            for _ in 0..1 {
                let degree = max_degree;
                let polynomial = Polynomial::rand(degree, rng);
                let degree_bound = None;
                let hiding_bound = Some(num_points_in_query_set);
                let l_poly = LabeledPolynomial::new_owned(
                    String::from("Test"),
                    polynomial,
                    degree_bound,
                    hiding_bound,
                );
                polynomials.push(l_poly);
            }

            let (comms, rands) = PC::commit(&ck, &polynomials, Some(rng))?;

            // Construct query set
            let mut query_set = QuerySet::new();
            let mut values = Evaluations::new();
            for _ in 0..num_points_in_query_set {
                let point = F::rand(rng);
                for i in 0..1 {
                    query_set.insert(("Test", point));
                    let value = polynomials[i].evaluate(point);
                    values.insert(("Test", point), value);
                }
            }

            let opening_challenge = F::rand(rng);
            let proof = PC::open(&ck, &polynomials, &query_set, opening_challenge, &rands)?;
            assert!(
                PC::check(
                    &vk,
                    &comms,
                    &query_set,
                    &values,
                    &proof,
                    opening_challenge,
                    rng
                )?,
                "proof was incorrect"
            );
        }
        Ok(())
    }

    pub fn single_poly_degree_bound_test<F, PC>() -> Result<(), PC::Error>
    where
        F: Field,
        PC: PolynomialCommitment<F>,
    {
        let rng = &mut thread_rng();
        let max_degree = rand::distributions::Uniform::from(2..64).sample(rng);
        let (ck, vk) = PC::setup(max_degree, rng)?;
        for _ in 0..100 {
            let mut polynomials = Vec::new();
            // Generate polynomials
            let num_points_in_query_set = 1;
            for _ in 0..1 {
                let degree = rand::distributions::Uniform::from(1..max_degree).sample(rng);
                let poly = Polynomial::rand(degree, rng);

                let range = rand::distributions::Uniform::from(degree..max_degree);
                let degree_bound = Some(range.sample(rng));
                let hiding_bound = Some(num_points_in_query_set);
                let l_poly = LabeledPolynomial::new_owned(
                    String::from("Test"),
                    poly,
                    degree_bound,
                    hiding_bound,
                );
                polynomials.push(l_poly);
            }

            let (comms, rands) = PC::commit(&ck, &polynomials, Some(rng))?;
            println!("Committed");

            // Construct query set
            let mut query_set = QuerySet::new();
            let mut values = Evaluations::new();
            for _ in 0..num_points_in_query_set {
                let point = F::rand(rng);
                for i in 0..1 {
                    query_set.insert(("Test", point));
                    let value = polynomials[i].evaluate(point);
                    values.insert(("Test", point), value);
                }
            }

            println!("Generated query set");
            let opening_challenge = F::rand(rng);
            let proof = PC::open(&ck, &polynomials, &query_set, opening_challenge, &rands)?;
            println!("Opened");
            assert!(
                PC::check(
                    &vk,
                    &comms,
                    &query_set,
                    &values,
                    &proof,
                    opening_challenge,
                    rng
                )?,
                "proof was incorrect"
            );
        }
        Ok(())
    }

    pub fn single_poly_degree_bound_multiple_queries_test<F, PC>() -> Result<(), PC::Error>
    where
        F: Field,
        PC: PolynomialCommitment<F>,
    {
        let rng = &mut thread_rng();
        let max_degree = rand::distributions::Uniform::from(2..64).sample(rng);
        let (ck, vk) = PC::setup(max_degree, rng)?;
        for _ in 0..100 {
            let mut polynomials = Vec::new();
            // Generate polynomials
            let num_points_in_query_set = 2;
            for _ in 0..1 {
                let degree = rand::distributions::Uniform::from(1..max_degree).sample(rng);
                let poly = Polynomial::rand(degree, rng);

                let range = rand::distributions::Uniform::from(degree..max_degree);
                let degree_bound = Some(range.sample(rng));
                let hiding_bound = Some(num_points_in_query_set);
                let l_poly = LabeledPolynomial::new_owned(
                    String::from("Test"),
                    poly,
                    degree_bound,
                    hiding_bound,
                );
                polynomials.push(l_poly);
            }

            let (comms, rands) = PC::commit(&ck, &polynomials, Some(rng))?;

            // Construct query set
            let mut query_set = QuerySet::new();
            let mut values = Evaluations::new();
            for _ in 0..num_points_in_query_set {
                let point = F::rand(rng);
                for i in 0..1 {
                    query_set.insert(("Test", point));
                    let value = polynomials[i].evaluate(point);
                    values.insert(("Test", point), value);
                }
            }

            let opening_challenge = F::rand(rng);
            let proof = PC::open(&ck, &polynomials, &query_set, opening_challenge, &rands)?;
            assert!(
                PC::check(
                    &vk,
                    &comms,
                    &query_set,
                    &values,
                    &proof,
                    opening_challenge,
                    rng
                )?,
                "proof was incorrect"
            );
        }
        Ok(())
    }

    pub fn two_polys_degree_bound_single_query_test<F, PC>() -> Result<(), PC::Error>
    where
        F: Field,
        PC: PolynomialCommitment<F>,
    {
        let rng = &mut thread_rng();
        let max_degree = 10;
        let (ck, vk) = PC::setup(max_degree, rng)?;
        for _ in 0..100 {
            let mut polynomials = Vec::new();
            let mut labels = Vec::new();
            let num_points_in_query_set = 1;
            // Generate polynomials
            for i in 0..2 {

                let label = format!("Test{}", i);
                labels.push(label.clone());

                let degree = rand::distributions::Uniform::from(1..max_degree).sample(rng);
                let range = rand::distributions::Uniform::from(degree..max_degree);
                let poly = Polynomial::rand(degree, rng);
                let degree_bound = Some(range.sample(rng));
                let hiding_bound = Some(num_points_in_query_set);
                polynomials.push(LabeledPolynomial::new_owned(
                    label,
                    poly,
                    degree_bound,
                    hiding_bound,
                ))
            }

            let (comms, rands) = PC::commit(&ck, &polynomials, Some(rng))?;

            // Construct query set
            let mut query_set = QuerySet::new();
            let mut values = Evaluations::new();
            for _ in 0..num_points_in_query_set {
                let point = F::rand(rng);
                for (i, label) in labels.iter().enumerate() {
                    query_set.insert((label, point));
                    let value = polynomials[i].evaluate(point);
                    values.insert((&label, point), value);
                }
            }

            let opening_challenge = F::rand(rng);
            let proof = PC::open(&ck, &polynomials, &query_set, opening_challenge, &rands)?;
            assert!(
                PC::check(
                    &vk,
                    &comms,
                    &query_set,
                    &values,
                    &proof,
                    opening_challenge,
                    rng
                )?,
                "proof was incorrect"
            );
        }
        Ok(())
    }

    pub fn full_end_to_end_test<F, PC>() -> Result<(), PC::Error>
    where
        F: Field,
        PC: PolynomialCommitment<F>,
    {
        let rng = &mut thread_rng();
        let max_degree = rand::distributions::Uniform::from(2..64).sample(rng);
        let (ck, vk) = PC::setup(max_degree, rng)?;
        for _ in 0..50 {
            let mut polynomials = Vec::new();
            let mut labels = Vec::new();

            // Generate polynomials
            let num_points_in_query_set = std::cmp::min(max_degree - 1, rand::distributions::Uniform::from(1..6).sample(rng));
            for i in 0..10 {
                let label = format!("Test{}", i);
                labels.push(label.clone());
                let degree = rand::distributions::Uniform::from(1..max_degree).sample(rng);
                let poly = Polynomial::rand(degree, rng);
                let range = rand::distributions::Uniform::from(degree..max_degree);
                let degree_bound = Some(range.sample(rng));
                let hiding_bound = Some(num_points_in_query_set);
                polynomials.push(LabeledPolynomial::new_owned(
                    label,
                    poly,
                    degree_bound,
                    hiding_bound,
                ))
            }

            let (comms, rands) = PC::commit(&ck, &polynomials, Some(rng))?;

            // Construct query set
            let mut query_set = QuerySet::new();
            let mut values = Evaluations::new();
            // let mut point = F::one();
            for _ in 0..num_points_in_query_set {
                let point = F::rand(rng);
                for (i, label) in labels.iter().enumerate() {
                    let should_be_queried = bool::rand(rng);
                    if should_be_queried {
                        query_set.insert((label, point));
                        let value = polynomials[i].evaluate(point);
                        values.insert((label, point), value);
                    }
                }
            }

            println!();

            let opening_challenge = F::rand(rng);
            let proof = PC::open(&ck, &polynomials, &query_set, opening_challenge, &rands)?;
            let result = PC::check(
                &vk,
                &comms,
                &query_set,
                &values,
                &proof,
                opening_challenge,
                rng,
            )?;
            if !result {
                println!(
                    "Failed with 10 polynomials, num_points_in_query_set: {:?}",
                    num_points_in_query_set
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
}
