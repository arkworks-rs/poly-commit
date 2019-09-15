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
use rand::Rng;
use std::borrow::Cow;

/// Defines `SinglePolynomialCommitment` schemes that allow one to commit to
/// a single polynomial, and then provide an evaluation proof for that polynomial
/// at a single point.
pub mod single_pc;

/// Defines `MultiPolynomialCommitment` schemes that allow one to commit to
/// multiple polynomials, and then provide evaluation proofs for these polynomials
/// at many points.
pub mod multi_pc;

pub use multi_pc::MultiPolynomialCommitment;
pub use single_pc::SinglePolynomialCommitment;

/// Defines the minimal interface of committer keys for any polynomial
/// commitment scheme.
pub trait PCCommitterKey: Clone {
    /// Outputs the maximum degree supported by the committer key.
    fn max_degree(&self) -> usize;
}

/// Defines the minimal interface of verifier keys for any polynomial
/// commitment scheme.
pub trait PCVerifierKey: Clone {
    /// Outputs the maximum degree supported by the verifier key.
    fn max_degree(&self) -> usize;
}

/// Defines the minimal interface of commitments for any polynomial
/// commitment scheme.
pub trait PCCommitment: Clone + algebra::ToBytes {
    /// Outputs a non-hiding commitment to the zero polynomial.
    fn empty() -> Self;

    /// Does this commitment have a degree bound?
    fn has_degree_bound(&self) -> bool;

    /// Size in bytes
    fn size_in_bytes(&self) -> usize;
}

/// Defines the minimal interface of commitment randomness for any polynomial
/// commitment scheme.
pub trait PCRandomness: Clone {
    /// Outputs empty randomness that does not hide the commitment.
    fn empty() -> Self;

    /// Samples randomness for commitments;
    /// `num_queries` specifies the number of queries that the commitment will be opened at.
    fn rand<R: Rng>(num_queries: usize, rng: &mut R) -> Self;
}

/// A polynomial along with other information necessary for the HIOP protocol
/// and for the polynomial commitment scheme.
#[derive(Clone)]
pub struct LabeledPolynomial<'a, F: Field> {
    polynomial: Cow<'a, Polynomial<F>>,
    degree_bound: Option<usize>,
    hiding_bound: Option<usize>,
}

impl<'a, F: Field> std::ops::Deref for LabeledPolynomial<'a, F> {
    type Target = Polynomial<F>;

    fn deref(&self) -> &Self::Target {
        &self.polynomial
    }
}

impl<'a, F: Field> LabeledPolynomial<'a, F> {
    /// Instantiate a new polynomial_context.
    pub fn new_owned(
        polynomial: Polynomial<F>,
        degree_bound: Option<usize>,
        hiding_bound: Option<usize>,
    ) -> Self {
        Self {
            polynomial: Cow::Owned(polynomial),
            degree_bound,
            hiding_bound,
        }
    }

    /// Instantiate a new polynomial_context.
    pub fn new(
        polynomial: &'a Polynomial<F>,
        degree_bound: Option<usize>,
        hiding_bound: Option<usize>,
    ) -> Self {
        Self {
            polynomial: Cow::Borrowed(polynomial),
            degree_bound,
            hiding_bound,
        }
    }

    /// Retrieve the polynomial from `self`.
    pub fn polynomial(&self) -> &Polynomial<F> {
        &self.polynomial
    }

    /// Evaluate the polynomial in `self`.
    pub fn evaluate(&self, point: F) -> F {
        self.polynomial.evaluate(point)
    }

    /// Retrieve the degree bound in `self`.
    pub fn degree_bound(&self) -> Option<usize> {
        self.degree_bound
    }

    /// Retrieve whether the polynomial in `self` should be hidden.
    pub fn is_hiding(&self) -> bool {
        self.hiding_bound.is_some()
    }

    /// Retrieve the hiding bound for the polynomial in `self`.
    pub fn hiding_bound(&self) -> Option<usize> {
        self.hiding_bound
    }
}
