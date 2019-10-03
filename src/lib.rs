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

/// A polynomial along with information about its degree bound (if any), and the
/// maximum number of queries that will be made to it. This latter number determines
/// the amount of protection that will be provided to a commitment for this polynomial.
#[derive(Clone)]
pub struct LabeledPolynomial<'a, F: Field> {
    label: String,
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
    /// Construct a new labeled polynomial by consuming `polynomial`.
    pub fn new_owned(
        label: String,
        polynomial: Polynomial<F>,
        degree_bound: Option<usize>,
        hiding_bound: Option<usize>,
    ) -> Self {
        Self {
            label,
            polynomial: Cow::Owned(polynomial),
            degree_bound,
            hiding_bound,
        }
    }

    /// Construct a new labeled polynomial.
    pub fn new(
        label: String,
        polynomial: &'a Polynomial<F>,
        degree_bound: Option<usize>,
        hiding_bound: Option<usize>,
    ) -> Self {
        Self {
            label,
            polynomial: Cow::Borrowed(polynomial),
            degree_bound,
            hiding_bound,
        }
    }

    /// Return the label for `self`.
    pub fn label(&self) -> &str {
        &self.label
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

/// A commitment along with information about its degree bound (if any).
#[derive(Clone)]
pub struct LabeledCommitment<C: PCCommitment> {
    label: String,
    commitment: C,
    degree_bound: Option<usize>,
}

impl<C: PCCommitment> LabeledCommitment<C> {
    /// Instantiate a new polynomial_context.
    pub fn new(label: String, commitment: C, degree_bound: Option<usize>) -> Self {
        Self {
            label,
            commitment,
            degree_bound,
        }
    }

    /// Return the label for `self`.
    pub fn label(&self) -> &str {
        &self.label
    }

    /// Retrieve the polynomial from `self`.
    pub fn commitment(&self) -> &C {
        &self.commitment
    }

    /// Retrieve the degree bound in `self`.
    pub fn degree_bound(&self) -> Option<usize> {
        self.degree_bound
    }
}

impl<C: PCCommitment> algebra::ToBytes for LabeledCommitment<C> {
    #[inline]
    fn write<W: std::io::Write>(&self, writer: W) -> std::io::Result<()> {
        self.commitment.write(writer)
    }
}
