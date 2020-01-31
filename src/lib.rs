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


/// Errors pertaining to query sets.
pub mod error;
pub use error::*;

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
    type Proof: PCProof + Clone;
    /// The evaluation proof for a query set.
    type BatchProof: Clone + From<Vec<Self::Proof>> + Into<Vec<Self::Proof>>;
    /// The error type for the scheme.
    type Error: std::error::Error + From<QuerySetError> + From<EquationError>;

    /// Constructs public parameters when given as input the maximum degree `degree`
    /// for the polynomial commitment scheme.
    fn setup<R: RngCore>(
        max_degree: usize,
        rng: &mut R,
    ) -> Result<Self::UniversalParams, Self::Error>;

    /// Specializes the public parameters for polynomials up to the given `supported_degree`
    /// and for enforcing degree bounds in the range `1..=supported_degree`.
    fn trim(
        pp: &Self::UniversalParams,
        supported_degree: usize,
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
        polynomials: impl IntoIterator<Item = &'a LabeledPolynomial<'a, F>>,
        rng: Option<&mut dyn RngCore>,
    ) -> Result<
        (
            Vec<LabeledCommitment<Self::Commitment>>,
            Vec<Self::Randomness>,
        ),
        Self::Error,
    >;

    /// On input a list of labeled polynomials and a query point, `open` outputs a proof of evaluation
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
                let (polynomial, rand) = polynomials_with_rands
                    .get(label)
                    .ok_or(QuerySetError::MissingPolynomial { label: label.to_string() })?;
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
        evaluations: &Evaluations<F>,
        proof: &Self::BatchProof,
        opening_challenge: F,
        _rng: &mut R,
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

        // Implicit assumption: proofs are order in same manner as queries in
        // `query_to_labels_map`.
        let proofs: Vec<_> = proof.clone().into();
        assert_eq!(proofs.len(), query_to_labels_map.len());

        let mut result = true;
        for ((query, labels), proof) in query_to_labels_map.into_iter().zip(proofs) {
            let mut comms: Vec<&'_ LabeledCommitment<_>> = Vec::new();
            let mut values = Vec::new();
            for label in labels.into_iter() {
                let commitment = commitments
                    .get(label)
                    .ok_or(QuerySetError::MissingPolynomial { label: label.to_string() })?;

                let v_i = evaluations
                    .get(&(label, *query))
                    .ok_or(QuerySetError::MissingEvaluation { label: label.to_string() })?;

                comms.push(commitment);
                values.push(*v_i);
            }

            let proof_time = start_timer!(|| "Checking per-query proof");
            result &= Self::check(vk, comms, *query, values, &proof, opening_challenge)?;
            end_timer!(proof_time);
        }
        Ok(result)
    }

    /// On input a list of labeled polynomials and a query set, `open` outputs a proof of evaluation
    /// of the polynomials at the points in the query set.
    fn open_equations<'a>(
        ck: &Self::CommitterKey,
        equations: impl IntoIterator<Item = &'a Equation<F>>,
        labeled_polynomials: impl IntoIterator<Item = &'a LabeledPolynomial<'a, F>>,
        opening_challenge: F,
        rands: impl IntoIterator<Item = &'a Self::Randomness>,
    ) -> Result<Self::BatchProof, Self::Error>
    where
        Self::Randomness: 'a
    {
        let query_set = Equation::query_set(equations);
        Self::batch_open(ck, labeled_polynomials, &query_set, opening_challenge, rands)
    }

    /// Checks that `values` are the true evaluations at `query_set` of the polynomials
    /// committed in `labeled_commitments`.
    fn check_equations<'a, R: RngCore>(
        vk: &Self::VerifierKey,
        equations: impl IntoIterator<Item = &'a Equation<F>>,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Self::Commitment>>,
        evaluations: Option<&Evaluations<F>>,
        proof: &Self::BatchProof,
        opening_challenge: F,
        rng: &mut R,
    ) -> Result<bool, Self::Error>
    where
        Self::Commitment: 'a
    {
        let evaluations = evaluations.unwrap();
        let equations = equations.into_iter().collect::<Vec<_>>();
        for eqn in &equations {
            let claimed_rhs = eqn.rhs;
            let mut actual_rhs = F::zero();
            let eval_point = eqn.evaluation_point;

            for (coeff, label) in &eqn.lhs {
                let eval = evaluations
                    .get(&(label.as_str(), eval_point))
                    .ok_or(QuerySetError::MissingEvaluation { label: label.clone() })?;
                actual_rhs += &(*coeff * eval);
            }
            if claimed_rhs != actual_rhs {
                eprintln!("Equation {} failed to verify", eqn.label);
                return Ok(false);
            }
        }
        let query_set = Equation::query_set(equations);

        let pc_result = Self::batch_check(vk, commitments, &query_set, evaluations, proof, opening_challenge, rng)?;
        if !pc_result {
            eprintln!("Evaluation proofs failed to verify");
            return Ok(false);
        }

        Ok(true)
    }
}

#[cfg(test)]
pub mod tests {
    use crate::*;
    use algebra::Field;
    use rand::{distributions::Distribution, thread_rng, Rng};

    #[derive(Default)]
    struct TestInfo {
        num_iters: usize,
        max_degree: Option<usize>,
        supported_degree: Option<usize>,
        num_polynomials: usize,
        enforce_degree_bounds: bool,
        max_num_queries: usize,
        num_equations: Option<usize>
    }

    fn test_template<F, PC>(
        info: TestInfo,
    ) -> Result<(), PC::Error>
        where
            F: Field,
            PC: PolynomialCommitment<F>
    {
        let TestInfo {
            num_iters,
            max_degree,
            supported_degree,
            num_polynomials,
            enforce_degree_bounds,
            max_num_queries,
            ..
        } = info;

        let rng = &mut thread_rng();
        let max_degree = max_degree.unwrap_or(rand::distributions::Uniform::from(2..=64).sample(rng));
        let pp = PC::setup(max_degree, rng)?;

        for _ in 0..num_iters {
            let supported_degree =
                supported_degree.unwrap_or(rand::distributions::Uniform::from(1..=max_degree).sample(rng));
            assert!(max_degree >= supported_degree, "max_degree < supported_degree");
            let mut polynomials = Vec::new();
            let mut degree_bounds = if enforce_degree_bounds {
                Some(Vec::new())
            } else {
                None
            };

            let mut labels = Vec::new();
            println!("Sampled supported degree");

            // Generate polynomials
            let num_points_in_query_set = rand::distributions::Uniform::from(1..=max_num_queries).sample(rng);
            for i in 0..num_polynomials {
                let label = format!("Test{}", i);
                labels.push(label.clone());
                let degree = rand::distributions::Uniform::from(1..=supported_degree).sample(rng);
                let poly = Polynomial::rand(degree, rng);

                let degree_bound = if let Some(degree_bounds) = &mut degree_bounds {
                    let range = rand::distributions::Uniform::from(degree..=max_degree);
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
                println!("Hiding bound: {:?}", hiding_bound);

                polynomials.push(LabeledPolynomial::new_owned(
                    label,
                    poly,
                    degree_bound,
                    hiding_bound,
                ))
            }
            println!("supported degree: {:?}", supported_degree);
            println!("num_points_in_query_set: {:?}", num_points_in_query_set);
            let (ck, vk) = PC::trim(&pp, supported_degree, degree_bounds.as_ref().map(|s| s.as_slice()))?;
            println!("Trimmed");

            let (comms, rands) = PC::commit(&ck, &polynomials, Some(rng))?;

            // Construct query set
            let mut query_set = QuerySet::new();
            let mut values = Evaluations::new();
            // let mut point = F::one();
            for _ in 0..num_points_in_query_set {
                let point = F::rand(rng);
                for (i, label) in labels.iter().enumerate() {
                    query_set.insert((label, point));
                    let value = polynomials[i].evaluate(point);
                    values.insert((label, point), value);
                }
            }
            println!("Generated query set");

            let opening_challenge = F::rand(rng);
            let proof = PC::batch_open(&ck, &polynomials, &query_set, opening_challenge, &rands)?;
            let result = PC::batch_check(
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
                    "Failed with {} polynomials, num_points_in_query_set: {:?}",
                    num_polynomials,
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

    fn equation_test_template<F, PC>(
        info: TestInfo,
    ) -> Result<(), PC::Error>
        where
            F: Field,
            PC: PolynomialCommitment<F>
    {
        let TestInfo {
            num_iters,
            max_degree,
            supported_degree,
            num_polynomials,
            enforce_degree_bounds,
            max_num_queries,
            num_equations,
        } = info;

        let rng = &mut thread_rng();
        let max_degree = max_degree.unwrap_or(rand::distributions::Uniform::from(2..=64).sample(rng));
        let pp = PC::setup(max_degree, rng)?;

        for _ in 0..num_iters {
            let supported_degree =
                supported_degree.unwrap_or(rand::distributions::Uniform::from(1..=max_degree).sample(rng));
            assert!(max_degree >= supported_degree, "max_degree < supported_degree");
            let mut polynomials = Vec::new();
            let mut degree_bounds = if enforce_degree_bounds {
                Some(Vec::new())
            } else {
                None
            };

            let mut labels = Vec::new();
            println!("Sampled supported degree");

            // Generate polynomials
            let num_points_in_query_set = rand::distributions::Uniform::from(1..=max_num_queries).sample(rng);
            for i in 0..num_polynomials {
                let label = format!("Test{}", i);
                labels.push(label.clone());
                let degree = rand::distributions::Uniform::from(1..=supported_degree).sample(rng);
                let poly = Polynomial::rand(degree, rng);

                let degree_bound = if let Some(degree_bounds) = &mut degree_bounds {
                    if rng.gen() {
                        let range = rand::distributions::Uniform::from(degree..=max_degree);
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

                polynomials.push(LabeledPolynomial::new_owned(
                    label,
                    poly,
                    degree_bound,
                    hiding_bound,
                ))
            }
            println!("supported degree: {:?}", supported_degree);
            println!("num_points_in_query_set: {:?}", num_points_in_query_set);
            let (ck, vk) = PC::trim(&pp, supported_degree, degree_bounds.as_ref().map(|s| s.as_slice()))?;
            println!("Trimmed");

            let (comms, rands) = PC::commit(&ck, &polynomials, Some(rng))?;

            // Let's construct our equations
            let mut equations = Vec::new();
            for i in 0..num_points_in_query_set {
                let point = F::rand(rng);
                for j in 0..num_equations.unwrap() {
                    let mut equation = Equation::empty(format!("query {} eqn {}", i, j), point);
                    let should_have_degree_bounds: bool = rng.gen();
                    for (k, label) in labels.iter().enumerate() {
                        if should_have_degree_bounds {
                            let eval = polynomials[k].evaluate(point);
                            equation.push((F::one(), label.to_string()), eval);
                            break;
                        } else {
                            let poly = &polynomials[k];
                            if poly.degree_bound().is_some() {
                                continue;
                            } else {
                                assert!(poly.degree_bound().is_none());
                                let coeff = F::rand(rng);
                                let eval = coeff * poly.evaluate(point);
                                equation.push((coeff, label.to_string()), eval);
                            }
                        }
                    }
                    if !equation.lhs.is_empty() {
                        equations.push(equation);
                    }
                }
            }
            if equations.is_empty() {
                continue
            }
            println!("Generated query set");

            let opening_challenge = F::rand(rng);
            let proof = PC::open_equations(&ck, &equations, &polynomials, opening_challenge, &rands)?;
            let result = PC::check_equations(
                &vk,
                &equations,
                &comms,
                None,
                &proof,
                opening_challenge,
                rng,
            )?;
            if !result {
                println!(
                    "Failed with {} polynomials, num_points_in_query_set: {:?}",
                    num_polynomials,
                    num_points_in_query_set
                );
                println!("Degree of polynomials:",);
                for poly in polynomials {
                    println!("Degree: {:?}", poly.degree());

                }

            }
            assert!(result, "proof was incorrect, equations: {:#?}", equations);
        }
        Ok(())
    }



    pub fn single_poly_test<F, PC>() -> Result<(), PC::Error>
    where
        F: Field,
        PC: PolynomialCommitment<F>,
    {
        let info = TestInfo {
            num_iters: 100,
            max_degree: None,
            supported_degree: None,
            num_polynomials: 1,
            enforce_degree_bounds: false,
            max_num_queries: 1,
            ..Default::default()
        };
        test_template::<F, PC>(info)
    }

    pub fn linear_poly_degree_bound_test<F, PC>() -> Result<(), PC::Error>
    where
        F: Field,
        PC: PolynomialCommitment<F>,
    {
        let info = TestInfo {
            num_iters: 100,
            max_degree: Some(2),
            supported_degree: Some(1),
            num_polynomials: 1,
            enforce_degree_bounds: true,
            max_num_queries: 1,
            ..Default::default()
        };
        test_template::<F, PC>(info)
    }


    pub fn single_poly_degree_bound_test<F, PC>() -> Result<(), PC::Error>
    where
        F: Field,
        PC: PolynomialCommitment<F>,
    {
        let info = TestInfo {
            num_iters: 100,
            max_degree: None,
            supported_degree: None,
            num_polynomials: 1,
            enforce_degree_bounds: true,
            max_num_queries: 1,
            ..Default::default()
        };
        test_template::<F, PC>(info)
    }

    pub fn quadratic_poly_degree_bound_multiple_queries_test<F, PC>() -> Result<(), PC::Error>
    where
        F: Field,
        PC: PolynomialCommitment<F>,
    {
        let info = TestInfo {
            num_iters: 100,
            max_degree: Some(3),
            supported_degree: Some(2),
            num_polynomials: 1,
            enforce_degree_bounds: true,
            max_num_queries: 2,
            ..Default::default()
        };
        test_template::<F, PC>(info)
    }

    pub fn single_poly_degree_bound_multiple_queries_test<F, PC>() -> Result<(), PC::Error>
    where
        F: Field,
        PC: PolynomialCommitment<F>,
    {
        let info = TestInfo {
            num_iters: 100,
            max_degree: None,
            supported_degree: None,
            num_polynomials: 1,
            enforce_degree_bounds: true,
            max_num_queries: 2,
            ..Default::default()
        };
        test_template::<F, PC>(info)
    }

    pub fn two_polys_degree_bound_single_query_test<F, PC>() -> Result<(), PC::Error>
    where
        F: Field,
        PC: PolynomialCommitment<F>,
    {
        let info = TestInfo {
            num_iters: 100,
            max_degree: None,
            supported_degree: None,
            num_polynomials: 2,
            enforce_degree_bounds: true,
            max_num_queries: 1,
            ..Default::default()
        };
        test_template::<F, PC>(info)
    }

    pub fn full_end_to_end_test<F, PC>() -> Result<(), PC::Error>
    where
        F: Field,
        PC: PolynomialCommitment<F>,
    {
        let info = TestInfo {
            num_iters: 100,
            max_degree: None,
            supported_degree: None,
            num_polynomials: 10,
            enforce_degree_bounds: true,
            max_num_queries: 5,
            ..Default::default()
        };
        test_template::<F, PC>(info)
    }

    pub fn full_end_to_end_equation_test<F, PC>() -> Result<(), PC::Error>
    where
        F: Field,
        PC: PolynomialCommitment<F>,
    {
        let info = TestInfo {
            num_iters: 100,
            max_degree: None,
            supported_degree: None,
            num_polynomials: 10,
            enforce_degree_bounds: true,
            max_num_queries: 5,
            num_equations: Some(10),
        };
        equation_test_template::<F, PC>(info)
    }

    pub fn single_equation_test<F, PC>() -> Result<(), PC::Error>
    where
        F: Field,
        PC: PolynomialCommitment<F>,
    {
        let info = TestInfo {
            num_iters: 100,
            max_degree: None,
            supported_degree: None,
            num_polynomials: 1,
            enforce_degree_bounds: false,
            max_num_queries: 1,
            num_equations: Some(1),
        };
        equation_test_template::<F, PC>(info)
    }

    pub fn two_equation_test<F, PC>() -> Result<(), PC::Error>
    where
        F: Field,
        PC: PolynomialCommitment<F>,
    {
        let info = TestInfo {
            num_iters: 100,
            max_degree: None,
            supported_degree: None,
            num_polynomials: 2,
            enforce_degree_bounds: false,
            max_num_queries: 1,
            num_equations: Some(2),
        };
        equation_test_template::<F, PC>(info)
    }

    pub fn two_equation_degree_bound_test<F, PC>() -> Result<(), PC::Error>
    where
        F: Field,
        PC: PolynomialCommitment<F>,
    {
        let info = TestInfo {
            num_iters: 100,
            max_degree: None,
            supported_degree: None,
            num_polynomials: 2,
            enforce_degree_bounds: true,
            max_num_queries: 1,
            num_equations: Some(2),
        };
        equation_test_template::<F, PC>(info)
    }
}
