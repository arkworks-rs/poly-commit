use algebra::Field;
use rand::Rng;
use std::borrow::Borrow;
use std::collections::{BTreeMap, BTreeSet};

use crate::*;

/// `QuerySet` is the set of queries that are to be made to a list of polynomials
/// `p` that have previously been committed. Each element of a `QuerySet` is a `(index, query)`
/// pair, where `index` is the index into `p`, and `query` is the field element
/// that `p[index]` is to be queried at.
pub type QuerySet<F> = BTreeSet<(usize, F)>;

/// `Evaluations` is the result of querying a list of polynomials `p` at a `QuerySet`
/// `Q`. It maps each element of `Q` to the resulting evaluation. That is,
/// if `(index, query)` is an element of `Q`, then `evaluation.get((index, query))`
/// should equal `p[index].evaluate(query)`.
pub type Evaluations<F> = BTreeMap<(usize, F), F>;

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
    fn setup<R: Rng>(degree: usize, rng: &mut R) -> Result<(Self::CommitterKey, Self::VerifierKey), Self::Error>;

    /// Outputs a commitments to `polynomials`. If `hiding_bounds[i].is_some()`,
    /// then the `i`-th commitment is hiding up to `hiding_bounds[i]` number of queries.
    /// `rng` should not be `None` if `hiding_bounds[i].is_some()`
    /// is true.
    ///
    /// If `hiding_bounds[i].is_none()`, then the `i`-th randomness is
    /// `Self::Randomness::empty()`.
    // TODO: technically each iterator below can have different lifetimes, but
    // that becomes quite noisy. Also, lifetime elision in argument position
    // should fix this (https://github.com/rust-lang/rust/issues/49287)
    fn commit<'a>(
        ck: &Self::CommitterKey,
        polynomials: impl IntoIterator<Item = &'a Polynomial<F>>,
        degree_bounds: impl IntoIterator<Item = &'a Option<usize>>,
        hiding_bounds: impl IntoIterator<Item = &'a Option<usize>>,
        rng: Option<&mut dyn Rng>,
    ) -> Result<(Vec<Self::Commitment>, Vec<Self::Randomness>), Self::Error>;

    /// On input a list of polynomials and a query set, `open` outputs a proof of evaluation
    /// of the polynomials at the points in the query set.
    fn open(
        ck: &Self::CommitterKey,
        polynomials: &[impl Borrow<Polynomial<F>>],
        degree_bounds: &[Option<usize>],
        query_set: &QuerySet<F>,
        opening_challenge: F,
        r: &[impl Borrow<Self::Randomness>],
    ) -> Result<Self::Proof, Self::Error>;

    /// Checks that `values` are the true evaluations at `query_set` of the polynomials
    /// committed in `comm`.
    fn check<R: Rng>(
        vk: &Self::VerifierKey,
        comm: &[Self::Commitment],
        degree_bounds: &[Option<usize>],
        query_set: &QuerySet<F>,
        values: &Evaluations<F>,
        proof: &Self::Proof,
        opening_challenge: F,
        rng: &mut R,
    ) -> Result<bool, Self::Error>;


    /// Commit to labeled polynomials.
    fn commit_labeled<'a>(
        ck: &Self::CommitterKey,
        labeled_polynomials: impl Iterator<Item = &'a LabeledPolynomial<'a, F>>,
        rng: Option<&mut dyn Rng>,
    ) -> Result<(Vec<Self::Commitment>, Vec<Self::Randomness>), Self::Error> {
        let mut polynomials = Vec::new();
        let mut degree_bounds = Vec::new();
        let mut hiding_bounds = Vec::new();
        for labeled_poly in labeled_polynomials {
            polynomials.push(labeled_poly.polynomial());
            degree_bounds.push(labeled_poly.degree_bound());
            hiding_bounds.push(labeled_poly.hiding_bound());
        }
        Self::commit(
            ck,
            polynomials,
            &degree_bounds,
            &hiding_bounds,
            rng,
        )
    }

    /// Open labeled polynomials.
    fn open_labeled<'a>(
        ck: &Self::CommitterKey,
        labeled_polynomials: impl Iterator<Item = &'a LabeledPolynomial<'a, F>>,
        query_set: &QuerySet<F>,
        opening_challenge: F,
        rands: &[Self::Randomness],
    ) -> Result<Self::Proof, Self::Error> {
        let mut polynomials = Vec::new();
        let mut degree_bounds = Vec::new();
        let mut hiding_bounds = Vec::new();
        for labeled_poly in labeled_polynomials {
            polynomials.push(labeled_poly.polynomial());
            degree_bounds.push(labeled_poly.degree_bound());
            hiding_bounds.push(labeled_poly.hiding_bound());
        }

        Self::open(
            &ck,
            &polynomials,
            &degree_bounds,
            &query_set,
            opening_challenge,
            &rands,
        )
    }
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
    use rand::{distributions::Sample, thread_rng, Rand};

    pub fn single_poly_test<F, MultiPC>() -> Result<(), MultiPC::Error>
    where
        F: Field,
        MultiPC: MultiPolynomialCommitment<F>,
    {
        let rng = &mut thread_rng();
        let max_degree = rand::distributions::Range::new(1, 64).sample(rng);;
        let (ck, vk) = MultiPC::setup(max_degree, rng)?;
        for _ in 0..100 {
            let mut polynomials = Vec::new();
            let mut degree_bounds = Vec::new();
            // Generate polynomials
            for _ in 0..1 {
                let degree = max_degree;
                polynomials.push(Polynomial::rand(degree, rng));
                degree_bounds.push(None);
            }

            let num_points_in_query_set = 1;
            let hiding_bounds = vec![Some(num_points_in_query_set); 1];

            let (comms, rands) = MultiPC::commit(
                &ck,
                &polynomials,
                &degree_bounds,
                &hiding_bounds,
                Some(rng),
            )?;

            // Construct query set
            let mut query_set = QuerySet::new();
            let mut values = Evaluations::new();
            for _ in 0..num_points_in_query_set {
                let point = F::rand(rng);
                for i in 0..1 {
                    query_set.insert((i, point));
                    let value = polynomials[i].evaluate(point);
                    values.insert((i, point), value);
                }
            }

            let opening_challenge = F::rand(rng);
            let proof = MultiPC::open(
                &ck,
                &polynomials,
                &degree_bounds,
                &query_set,
                opening_challenge,
                &rands,
            )?;
            assert!(
                MultiPC::check(
                    &vk,
                    &comms,
                    &degree_bounds,
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
        let max_degree = rand::distributions::Range::new(1, 64).sample(rng);;
        let (ck, vk) = MultiPC::setup(max_degree, rng)?;
        for _ in 0..100 {
            let mut polynomials = Vec::new();
            let mut degree_bounds = Vec::new();
            // Generate polynomials
            for _ in 0..1 {
                let degree = rand::distributions::Range::new(1, max_degree).sample(rng);
                polynomials.push(Polynomial::rand(degree, rng));

                let mut range = rand::distributions::Range::new(degree, max_degree);
                let degree_bound = Some(range.sample(rng));
                degree_bounds.push(degree_bound);
            }

            let num_points_in_query_set = 1;
            let hiding_bounds = vec![Some(num_points_in_query_set); 1];

            let (comms, rands) = MultiPC::commit(
                &ck,
                &polynomials,
                &degree_bounds,
                &hiding_bounds,
                Some(rng),
            )?;
            println!("Committed");

            // Construct query set
            let mut query_set = QuerySet::new();
            let mut values = Evaluations::new();
            for _ in 0..num_points_in_query_set {
                let point = F::rand(rng);
                for i in 0..1 {
                    query_set.insert((i, point));
                    let value = polynomials[i].evaluate(point);
                    values.insert((i, point), value);
                }
            }

            println!("Generated query set");
            let opening_challenge = F::rand(rng);
            let proof = MultiPC::open(
                &ck,
                &polynomials,
                &degree_bounds,
                &query_set,
                opening_challenge,
                &rands,
            )?;
            println!("Opened");
            assert!(
                MultiPC::check(
                    &vk,
                    &comms,
                    &degree_bounds,
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
        let max_degree = rand::distributions::Range::new(1, 64).sample(rng);;
        let (ck, vk) = MultiPC::setup(max_degree, rng)?;
        for _ in 0..100 {
            let mut polynomials = Vec::new();
            let mut degree_bounds = Vec::new();
            // Generate polynomials
            for _ in 0..1 {
                let degree = rand::distributions::Range::new(1, max_degree).sample(rng);
                polynomials.push(Polynomial::rand(degree, rng));

                let mut range = rand::distributions::Range::new(degree, max_degree);
                let degree_bound = Some(range.sample(rng));
                degree_bounds.push(degree_bound);
            }

            let num_points_in_query_set = 2;
            let hiding_bounds = vec![Some(num_points_in_query_set); 1];

            let (comms, rands) = MultiPC::commit(
                &ck,
                &polynomials,
                &degree_bounds,
                &hiding_bounds,
                Some(rng),
            )?;

            // Construct query set
            let mut query_set = QuerySet::new();
            let mut values = Evaluations::new();
            for _ in 0..num_points_in_query_set {
                let point = F::rand(rng);
                for i in 0..1 {
                    query_set.insert((i, point));
                    let value = polynomials[i].evaluate(point);
                    values.insert((i, point), value);
                }
            }

            let opening_challenge = F::rand(rng);
            let proof = MultiPC::open(
                &ck,
                &polynomials,
                &degree_bounds,
                &query_set,
                opening_challenge,
                &rands,
            )?;
            assert!(
                MultiPC::check(
                    &vk,
                    &comms,
                    &degree_bounds,
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
            let mut degree_bounds = Vec::new();
            // Generate polynomials
            for _ in 0..2 {
                let degree = rand::distributions::Range::new(1, max_degree).sample(rng);
                polynomials.push(Polynomial::rand(degree, rng));

                let mut range = rand::distributions::Range::new(degree, max_degree);
                let degree_bound = Some(range.sample(rng));
                degree_bounds.push(degree_bound);
            }

            let num_points_in_query_set = 1;
            let hiding_bounds = vec![Some(num_points_in_query_set); 2];

            let (comms, rands) = MultiPC::commit(
                &ck,
                &polynomials,
                &degree_bounds,
                &hiding_bounds,
                Some(rng),
            )?;

            // Construct query set
            let mut query_set = QuerySet::new();
            let mut values = Evaluations::new();
            for _ in 0..num_points_in_query_set {
                let point = F::rand(rng);
                for i in 0..2 {
                    query_set.insert((i, point));
                    let value = polynomials[i].evaluate(point);
                    values.insert((i, point), value);
                }
            }

            let opening_challenge = F::rand(rng);
            let proof = MultiPC::open(
                &ck,
                &polynomials,
                &degree_bounds,
                &query_set,
                opening_challenge,
                &rands,
            )?;
            assert!(
                MultiPC::check(
                    &vk,
                    &comms,
                    &degree_bounds,
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
        let max_degree = rand::distributions::Range::new(2, 64).sample(rng);;
        let (ck, vk) = MultiPC::setup(max_degree, rng)?;
        for _ in 0..10 {
            let mut polynomials = Vec::new();
            let mut degree_bounds = Vec::new();
            // Generate polynomials
            for _ in 0..10 {
                let degree = rand::distributions::Range::new(1, max_degree).sample(rng);
                polynomials.push(Polynomial::rand(degree, rng));

                let mut range = rand::distributions::Range::new(degree, max_degree);
                let degree_bound = Some(range.sample(rng));
                degree_bounds.push(degree_bound);
            }

            let num_points_in_query_set = rand::distributions::Range::new(1, 5).sample(rng);
            let hiding_bounds = vec![Some(num_points_in_query_set); 10];

            let (comms, rands) = MultiPC::commit(
                &ck,
                &polynomials,
                &degree_bounds,
                &hiding_bounds,
                Some(rng),
            )?;

            // Construct query set
            let mut query_set = QuerySet::new();
            let mut values = Evaluations::new();
            for _ in 0..num_points_in_query_set {
                let point = F::rand(rng);
                for i in 0..10 {
                    let should_be_queried = bool::rand(rng);
                    if should_be_queried {
                        query_set.insert((i, point));
                        let value = polynomials[i].evaluate(point);
                        values.insert((i, point), value);
                    }
                }
            }

            let opening_challenge = F::rand(rng);
            let proof = MultiPC::open(
                &ck,
                &polynomials,
                &degree_bounds,
                &query_set,
                opening_challenge,
                &rands,
            )?;
            let result = MultiPC::check(
                &vk,
                &comms,
                &degree_bounds,
                &query_set,
                &values,
                &proof,
                opening_challenge,
                rng
            )?;
            if !result {
                println!("Failed with 10 polynomials, num_points_in_query_set: {:?}", num_points_in_query_set);
            }
            assert!(result, "proof was incorrect");
        }
        Ok(())
    }
}
