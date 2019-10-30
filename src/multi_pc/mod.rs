use algebra::Field;
use rand::RngCore;
use std::collections::{BTreeMap, BTreeSet};

use crate::*;

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
pub trait MultiPolynomialCommitment<F: Field> {
    /// The committer key for the scheme; used to commit to a polynomial and then
    /// open the commitment to produce an evaluation proof.
    type CommitterKey: PCCommitterKey;
    /// The verifier key for the scheme; used to check an evaluation proof.
    type VerifierKey: PCVerifierKey;
    /// The commitment to a polynomial.
    type Commitment: PCCommitment;
    /// The commitment randomness.
    type Randomness: PCRandomness;
    /// The evaluation proof.
    type Proof: Clone;
    /// The error type for the scheme.
    type Error: std::error::Error;

    /// Constructs public parameters when given as input the maximum degree `degree`
    /// for the polynomial commitment scheme.
    fn setup<R: RngCore>(
        degree: usize,
        rng: &mut R,
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

    /// On input a list of labeled polynomials and a query set, `open` outputs a proof of evaluation
    /// of the polynomials at the points in the query set.
    fn open<'a>(
        ck: &Self::CommitterKey,
        polynomials: impl IntoIterator<Item = &'a LabeledPolynomial<'a, F>>,
        query_set: &QuerySet<F>,
        opening_challenge: F,
        rands: &[Self::Randomness],
    ) -> Result<Self::Proof, Self::Error>;

    /// Checks that `values` are the true evaluations at `query_set` of the polynomials
    /// committed in `comm`.
    fn check<R: RngCore>(
        vk: &Self::VerifierKey,
        labeled_commitments: &[LabeledCommitment<Self::Commitment>],
        query_set: &QuerySet<F>,
        values: &Evaluations<F>,
        proof: &Self::Proof,
        opening_challenge: F,
        rng: &mut R,
    ) -> Result<bool, Self::Error>;
}

/// Generic construction of a `MultiPolynomialCommitment` scheme from a
/// `SinglePolynomialCommitment` scheme whenever the commitment and randomness of the
/// `SinglePolynomialCommitment` scheme are additively homomorphic.
/// Specifically, we require `C = MultiPolynomialCommitment::Commitment`
/// to satisfy `for<'a> C: AddAssign<(F, &'a C)>.
///
/// The construction follows the blueprint laid out in [CHMMVW19](insert eprint link).
pub mod mpc_from_spc;

#[cfg(test)]
pub mod tests {
    use crate::multi_pc::*;
    use algebra::Field;
    use algebra::UniformRand;
    use rand::{distributions::Distribution, thread_rng};

    pub fn single_poly_test<F, MultiPC>() -> Result<(), MultiPC::Error>
    where
        F: Field,
        MultiPC: MultiPolynomialCommitment<F>,
    {
        let rng = &mut thread_rng();
        let max_degree = rand::distributions::Uniform::from(2..64).sample(rng);
        let (ck, vk) = MultiPC::setup(max_degree, rng)?;
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

            let (comms, rands) = MultiPC::commit(&ck, &polynomials, Some(rng))?;

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
            let proof = MultiPC::open(&ck, &polynomials, &query_set, opening_challenge, &rands)?;
            assert!(
                MultiPC::check(
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

    pub fn single_poly_degree_bound_test<F, MultiPC>() -> Result<(), MultiPC::Error>
    where
        F: Field,
        MultiPC: MultiPolynomialCommitment<F>,
    {
        let rng = &mut thread_rng();
        let max_degree = rand::distributions::Uniform::from(2..64).sample(rng);
        let (ck, vk) = MultiPC::setup(max_degree, rng)?;
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

            let (comms, rands) = MultiPC::commit(&ck, &polynomials, Some(rng))?;
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
            let proof = MultiPC::open(&ck, &polynomials, &query_set, opening_challenge, &rands)?;
            println!("Opened");
            assert!(
                MultiPC::check(
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

    pub fn single_poly_degree_bound_multiple_queries_test<F, MultiPC>() -> Result<(), MultiPC::Error>
    where
        F: Field,
        MultiPC: MultiPolynomialCommitment<F>,
    {
        let rng = &mut thread_rng();
        let max_degree = rand::distributions::Uniform::from(2..64).sample(rng);
        let (ck, vk) = MultiPC::setup(max_degree, rng)?;
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

            let (comms, rands) = MultiPC::commit(&ck, &polynomials, Some(rng))?;

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
            let proof = MultiPC::open(&ck, &polynomials, &query_set, opening_challenge, &rands)?;
            assert!(
                MultiPC::check(
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

    pub fn two_polys_degree_bound_single_query_test<F, MultiPC>() -> Result<(), MultiPC::Error>
    where
        F: Field,
        MultiPC: MultiPolynomialCommitment<F>,
    {
        let rng = &mut thread_rng();
        let max_degree = 10;
        let (ck, vk) = MultiPC::setup(max_degree, rng)?;
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

            let (comms, rands) = MultiPC::commit(&ck, &polynomials, Some(rng))?;

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
            let proof = MultiPC::open(&ck, &polynomials, &query_set, opening_challenge, &rands)?;
            assert!(
                MultiPC::check(
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

    pub fn full_end_to_end_test<F, MultiPC>() -> Result<(), MultiPC::Error>
    where
        F: Field,
        MultiPC: MultiPolynomialCommitment<F>,
    {
        let rng = &mut thread_rng();
        let max_degree = rand::distributions::Uniform::from(2..64).sample(rng);
        let (ck, vk) = MultiPC::setup(max_degree, rng)?;
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

            let (comms, rands) = MultiPC::commit(&ck, &polynomials, Some(rng))?;

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
            let proof = MultiPC::open(&ck, &polynomials, &query_set, opening_challenge, &rands)?;
            let result = MultiPC::check(
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
