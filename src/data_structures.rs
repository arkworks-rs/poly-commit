use crate::{Cow, MultiPolynomial, String, UniPolynomial, Vec};
use algebra_core::Field;
use core::borrow::Borrow;
use core::ops::{AddAssign, MulAssign, SubAssign};
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
    /// `num_vars` specifies the number of variables for multivariate commitment.
    /// strict degree bound.
    fn rand<R: RngCore>(
        num_queries: usize,
        has_degree_bound: bool,
        num_vars: Option<usize>,
        rng: &mut R,
    ) -> Self;
}

/// Defines the minimal interface of evaluation proofs for any polynomial
/// commitment scheme.
pub trait PCProof: Clone + algebra_core::ToBytes {
    /// Size in bytes
    fn size_in_bytes(&self) -> usize;
}

/// Wrapper for univariate or multivariate polynomial
#[derive(Debug, Clone)]
pub enum Polynomial<'a, F: Field> {
    /// Univariate polynomial type
    Uni(Cow<'a, UniPolynomial<F>>),
    /// Multivariate polynomial type
    Mul(Cow<'a, MultiPolynomial<F>>),
}

impl<'a, F: 'a + Field> From<UniPolynomial<F>> for Polynomial<'a, F> {
    fn from(other: UniPolynomial<F>) -> Self {
        Self::Uni(Cow::Owned(other))
    }
}

impl<'a, F: 'a + Field> From<&'a UniPolynomial<F>> for Polynomial<'a, F> {
    fn from(other: &'a UniPolynomial<F>) -> Self {
        Self::Uni(Cow::Borrowed(other))
    }
}

impl<'a, F: 'a + Field> From<MultiPolynomial<F>> for Polynomial<'a, F> {
    fn from(other: MultiPolynomial<F>) -> Self {
        Self::Mul(Cow::Owned(other))
    }
}

impl<'a, F: 'a + Field> From<&'a MultiPolynomial<F>> for Polynomial<'a, F> {
    fn from(other: &'a MultiPolynomial<F>) -> Self {
        Self::Mul(Cow::Borrowed(other))
    }
}

impl<'a, F: Field> core::convert::TryInto<&'a UniPolynomial<F>> for &'a Polynomial<'a, F> {
    type Error = crate::Error;

    fn try_into(self) -> Result<&'a UniPolynomial<F>, Self::Error> {
        match self {
            Polynomial::Uni(p) => Ok(&p),
            Polynomial::Mul(_) => Err(Self::Error::InvalidPolynomialType),
        }
    }
}

impl<'a, F: Field> core::convert::TryInto<&'a MultiPolynomial<F>> for &'a Polynomial<'a, F> {
    type Error = crate::Error;

    fn try_into(self) -> Result<&'a MultiPolynomial<F>, Self::Error> {
        match self {
            Polynomial::Uni(_) => Err(Self::Error::InvalidPolynomialType),
            Polynomial::Mul(p) => Ok(&p),
        }
    }
}

/// A polynomial along with information about its degree bound (if any), and the
/// maximum number of queries that will be made to it. This latter number determines
/// the amount of protection that will be provided to a commitment for this polynomial.
#[derive(Debug, Clone)]
pub struct LabeledPolynomial<'a, F: Field> {
    label: PolynomialLabel,
    polynomial: Polynomial<'a, F>,
    degree_bound: Option<usize>,
    hiding_bound: Option<usize>,
}

impl<'a, F: Field> core::ops::Deref for LabeledPolynomial<'a, F> {
    type Target = Polynomial<'a, F>;

    fn deref(&self) -> &Self::Target {
        &self.polynomial
    }
}

impl<'a, F: Field> LabeledPolynomial<'a, F> {
    /// Construct a new labeled polynomial.
    pub fn new(
        label: PolynomialLabel,
        polynomial: Polynomial<'a, F>,
        degree_bound: Option<usize>,
        hiding_bound: Option<usize>,
    ) -> Self {
        Self {
            label,
            polynomial: polynomial,
            degree_bound,
            hiding_bound,
        }
    }

    /// Return the label for `self`.
    pub fn label(&self) -> &String {
        &self.label
    }

    /// Retrieve the polynomial from `self`
    pub fn polynomial(&self) -> &Polynomial<'a, F> {
        &self.polynomial
    }

    /// Evaluate the polynomial in `self`.
    pub fn evaluate(&self, point: &[F]) -> F {
        match &self.polynomial {
            Polynomial::Uni(p) => p.evaluate(point[0]),
            Polynomial::Mul(p) => p.evaluate(point),
        }
    }

    /// Retrieve the degree of the polynomial in `self`.
    pub fn degree(&self) -> usize {
        match &self.polynomial {
            Polynomial::Uni(p) => p.degree(),
            Polynomial::Mul(p) => p.degree(),
        }
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
    pub fn label(&self) -> &String {
        &self.label
    }

    /// Retrieve the commitment from `self`.
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

/// A term in a linear combination.
#[derive(Hash, Ord, PartialOrd, Clone, Eq, PartialEq, Debug)]
pub enum LCTerm {
    /// The constant term representing `one`.
    One,
    /// Label for a polynomial.
    PolyLabel(String),
}

impl LCTerm {
    /// Returns `true` if `self == LCTerm::One`
    #[inline]
    pub fn is_one(&self) -> bool {
        if let LCTerm::One = self {
            true
        } else {
            false
        }
    }
}

impl From<PolynomialLabel> for LCTerm {
    fn from(other: PolynomialLabel) -> Self {
        Self::PolyLabel(other)
    }
}

impl<'a> From<&'a str> for LCTerm {
    fn from(other: &str) -> Self {
        Self::PolyLabel(other.into())
    }
}

impl core::convert::TryInto<PolynomialLabel> for LCTerm {
    type Error = ();
    fn try_into(self) -> Result<PolynomialLabel, ()> {
        match self {
            Self::One => Err(()),
            Self::PolyLabel(l) => Ok(l),
        }
    }
}

impl<'a> core::convert::TryInto<&'a PolynomialLabel> for &'a LCTerm {
    type Error = ();

    fn try_into(self) -> Result<&'a PolynomialLabel, ()> {
        match self {
            LCTerm::One => Err(()),
            LCTerm::PolyLabel(l) => Ok(l),
        }
    }
}

impl<B: Borrow<String>> PartialEq<B> for LCTerm {
    fn eq(&self, other: &B) -> bool {
        match self {
            Self::One => false,
            Self::PolyLabel(l) => l == other.borrow(),
        }
    }
}

/// A labeled linear combinations of polynomials.
#[derive(Clone, Debug)]
pub struct LinearCombination<F> {
    /// The label.
    pub label: String,
    /// The linear combination of `(coeff, poly_label)` pairs.
    terms: Vec<(F, LCTerm)>,
}

impl<F: Field> LinearCombination<F> {
    /// Construct an empty labeled linear combination.
    pub fn empty(label: impl Into<String>) -> Self {
        Self {
            label: label.into(),
            terms: Vec::new(),
        }
    }

    /// Construct a new labeled linear combination.
    /// with the terms specified in `term`.
    pub fn new(label: impl Into<String>, terms: Vec<(F, impl Into<LCTerm>)>) -> Self {
        let terms = terms.into_iter().map(|(c, t)| (c, t.into())).collect();
        Self {
            label: label.into(),
            terms: terms,
        }
    }

    /// Returns the label of the linear combination.
    pub fn label(&self) -> &String {
        &self.label
    }

    /// Returns `true` if the linear combination has no terms.
    pub fn is_empty(&self) -> bool {
        self.terms.is_empty()
    }

    /// Add a term to the linear combination.
    pub fn push(&mut self, term: (F, LCTerm)) -> &mut Self {
        self.terms.push(term);
        self
    }
}

impl<'a, F: Field> AddAssign<(F, &'a LinearCombination<F>)> for LinearCombination<F> {
    fn add_assign(&mut self, (coeff, other): (F, &'a LinearCombination<F>)) {
        self.terms
            .extend(other.terms.iter().map(|(c, t)| (coeff * c, t.clone())));
    }
}

impl<'a, F: Field> SubAssign<(F, &'a LinearCombination<F>)> for LinearCombination<F> {
    fn sub_assign(&mut self, (coeff, other): (F, &'a LinearCombination<F>)) {
        self.terms
            .extend(other.terms.iter().map(|(c, t)| (-coeff * c, t.clone())));
    }
}

impl<'a, F: Field> AddAssign<&'a LinearCombination<F>> for LinearCombination<F> {
    fn add_assign(&mut self, other: &'a LinearCombination<F>) {
        self.terms.extend(other.terms.iter().cloned());
    }
}

impl<'a, F: Field> SubAssign<&'a LinearCombination<F>> for LinearCombination<F> {
    fn sub_assign(&mut self, other: &'a LinearCombination<F>) {
        self.terms
            .extend(other.terms.iter().map(|(c, t)| (-*c, t.clone())));
    }
}

impl<F: Field> AddAssign<F> for LinearCombination<F> {
    fn add_assign(&mut self, coeff: F) {
        self.terms.push((coeff, LCTerm::One));
    }
}

impl<F: Field> SubAssign<F> for LinearCombination<F> {
    fn sub_assign(&mut self, coeff: F) {
        self.terms.push((-coeff, LCTerm::One));
    }
}

impl<F: Field> MulAssign<F> for LinearCombination<F> {
    fn mul_assign(&mut self, coeff: F) {
        self.terms.iter_mut().for_each(|(c, _)| *c *= coeff);
    }
}

impl<F: Field> core::ops::Deref for LinearCombination<F> {
    type Target = [(F, LCTerm)];

    fn deref(&self) -> &Self::Target {
        &self.terms
    }
}
