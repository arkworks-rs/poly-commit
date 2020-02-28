use crate::{Cow, String, Vec};
use algebra_core::Field;
pub use ff_fft::DensePolynomial as Polynomial;
use rand_core::RngCore;

/// Labels a `LabeledPolynomial` or a `LabeledCommitment`.
pub type PolynomialLabel = String;

/// Defines the minimal interface for public params for any polynomial
/// commitment scheme.
pub trait PCUniversalParams: Clone + core::fmt::Debug {
    /// Outputs the maximum degree supported by the committer key.
    fn max_degree(&self) -> usize;
}

/// Defines the minimal interface of committer keys for any polynomial
/// commitment scheme.
pub trait PCCommitterKey: Clone + core::fmt::Debug {
    /// Outputs the maximum degree supported by the universal parameters
    /// `Self` was derived from.
    fn max_degree(&self) -> usize;

    /// Outputs the maximum degree supported by the committer key.
    fn supported_degree(&self) -> usize;
}

/// Defines the minimal interface of verifier keys for any polynomial
/// commitment scheme.
pub trait PCVerifierKey: Clone + core::fmt::Debug {
    /// Outputs the maximum degree supported by the universal parameters
    /// `Self` was derived from.
    fn max_degree(&self) -> usize;

    /// Outputs the maximum degree supported by the verifier key.
    fn supported_degree(&self) -> usize;
}

/// Defines the minimal interface of commitments for any polynomial
/// commitment scheme.
pub trait PCCommitment: Clone + algebra_core::ToBytes {
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
    /// `has_degree_bound` indicates that the corresponding commitment has an enforced
    /// strict degree bound.
    fn rand<R: RngCore>(num_queries: usize, has_degree_bound: bool, rng: &mut R) -> Self;
}

/// Defines the minimal interface of evaluation proofs for any polynomial
/// commitment scheme.
pub trait PCProof: Clone + algebra_core::ToBytes {
    /// Size in bytes
    fn size_in_bytes(&self) -> usize;
}

/// A polynomial along with information about its degree bound (if any), and the
/// maximum number of queries that will be made to it. This latter number determines
/// the amount of protection that will be provided to a commitment for this polynomial.
#[derive(Debug, Clone)]
pub struct LabeledPolynomial<'a, F: Field> {
    label: PolynomialLabel,
    polynomial: Cow<'a, Polynomial<F>>,
    degree_bound: Option<usize>,
    hiding_bound: Option<usize>,
}

impl<'a, F: Field> core::ops::Deref for LabeledPolynomial<'a, F> {
    type Target = Polynomial<F>;

    fn deref(&self) -> &Self::Target {
        &self.polynomial
    }
}

impl<'a, F: Field> LabeledPolynomial<'a, F> {
    /// Construct a new labeled polynomial by consuming `polynomial`.
    pub fn new_owned(
        label: PolynomialLabel,
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
        label: PolynomialLabel,
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
    label: PolynomialLabel,
    commitment: C,
    degree_bound: Option<usize>,
}

impl<C: PCCommitment> LabeledCommitment<C> {
    /// Instantiate a new polynomial_context.
    pub fn new(label: PolynomialLabel, commitment: C, degree_bound: Option<usize>) -> Self {
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

impl<C: PCCommitment> algebra_core::ToBytes for LabeledCommitment<C> {
    #[inline]
    fn write<W: algebra_core::io::Write>(&self, writer: W) -> algebra_core::io::Result<()> {
        self.commitment.write(writer)
    }
}

/// A linear equation where the LHS consists of linear combinations of polynomials,
/// while the RHS contains a claimed evaluation of the LHS at a challenge
/// point.
#[derive(Clone, Debug)]
pub struct Equation<F> {
    /// The label for the equation.
    pub label: String,
    /// The RHS of the equation, consisting of `(coeff, poly_label)` pairs.
    pub lhs: Vec<(F, PolynomialLabel)>,
    /// The LHS of the equation, consisting of the evaluation of `self.lhs` at
    /// `self.evaluation_point`.
    pub rhs: F,
    /// The point that satisfies the equation.
    pub evaluation_point: F,
}

impl<F: Field> Equation<F> {
    /// Construct a new labeled equation.
    pub fn empty(label: String, evaluation_point: F) -> Self {
        Self {
            label,
            lhs: Vec::new(),
            rhs: F::zero(),
            evaluation_point,
        }
    }

    /// Returns `true` if the LHS of `self` is empty.
    pub fn lhs_is_empty(&self) -> bool {
        self.lhs.is_empty()
    }

    /// Add a term to the equation, updating the LHS and RHS in the process.
    /// It *must* be the case that `eval` is the evaluation of the provided term
    /// at `self.evaluation_point`.
    pub fn push(&mut self, term: (F, PolynomialLabel), eval: F) {
        self.rhs += &eval;
        self.lhs.push(term);
    }

    /// Obtain a query set from the given equations.
    /// This method simply maps each polynomial in an equation into its own
    /// entry in the query set.
    pub fn query_set<'a>(equations: impl IntoIterator<Item = &'a Self>) -> crate::QuerySet<'a, F> {
        equations
            .into_iter()
            .flat_map(|eqn| {
                let point = eqn.evaluation_point;
                eqn.lhs.iter().map(move |poly| (poly.1.as_str(), point))
            })
            .collect()
    }
}
