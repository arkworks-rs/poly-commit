use crate::String;

/// The error type for `PolynomialCommitment`.
#[derive(Debug)]
pub enum Error {
    /// The query set contains a label for a polynomial that was not provided as
    /// input to the `PC::open`.
    MissingPolynomial {
        /// The label of the missing polynomial.
        label: String,
    },

    /// `Evaluations` does not contain an evaluation for the polynomial labelled
    /// `label` at a particular query.
    MissingEvaluation {
        /// The label of the missing polynomial.
        label: String,
    },

    /// The LHS of the equation is empty.
    MissingLHS {
        /// The label of the equation.
        label: String,
    },

    /// The provided polynomial was meant to be hiding, but `rng` was `None`.
    MissingRng,

    /// The degree provided in setup was too small; degree 0 polynomials
    /// are not supported.
    DegreeIsZero,

    /// The degree of the polynomial passed to `commit` or `open`
    /// was too large.
    TooManyCoefficients {
        /// The number of coefficients in the polynomial.
        num_coefficients: usize,
        /// The maximum number of powers provided in `Powers`.
        num_powers: usize,
    },

    /// The hiding bound was not `None`, but the hiding bound was zero.
    HidingBoundIsZero,

    /// The hiding bound was too large for the given `Powers`.
    HidingBoundToolarge {
        /// The hiding bound
        hiding_poly_degree: usize,
        /// The number of powers.
        num_powers: usize,
    },

    /// The degree provided to `trim` was too large.
    TrimmingDegreeTooLarge,

    /// The provided `enforced_degree_bounds` was `Some<&[]>`.
    EmptyDegreeBounds,

    /// The provided equation contained multiple polynomials, of which least one
    /// had a strict degree bound.
    EquationHasDegreeBounds(String),

    /// The required degree bound is not supported by ck/vk
    UnsupportedDegreeBound(usize),

    /// The degree bound for the `index`-th polynomial passed to `commit`, `open`
    /// or `check` was incorrect, that is, `degree_bound >= poly_degree` or
    /// `degree_bound <= max_degree`.
    IncorrectDegreeBound {
        /// Degree of the polynomial.
        poly_degree: usize,
        /// Degree bound.
        degree_bound: usize,
        /// Maximum supported degree.
        supported_degree: usize,
        /// Index of the offending polynomial.
        label: String,
    },

    /// The inputs to `commit`, `open` or `verify` had incorrect lengths.
    IncorrectInputLength(String),

    /// An invalid number of variables was provided to `setup`
    InvalidNumberOfVariables,

    /// The degree of the `index`-th polynomial passed to `commit`, `open`
    /// or `check` was incorrect, that is, `supported_degree <= poly_degree`
    PolynomialDegreeTooLarge {
        /// Degree of the polynomial.
        poly_degree: usize,
        /// Maximum supported degree.
        supported_degree: usize,
        /// Index of the offending polynomial.
        label: String,
    },

    /// This means a failure in verifying the commitment or the opening.
    InvalidCommitment,

    /// For PCS which rely on Fiat-Shamir to be rendered non-interactive,
    /// these are errors that result from incorrect transcript manipulation.
    TranscriptError,

    /// This means the required soundness error bound is inherently impossible.
    /// E.g., the field is not big enough.
    InvalidParameters(String),

    /// Error resulting from hashing in linear code - based PCS.
    HashingError,

    /// Shows that encoding is not feasible
    EncodingError,
}

impl core::fmt::Display for Error {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Error::MissingPolynomial { label } => write!(
                f,
                "`QuerySet` refers to polynomial \"{}\", but it was not provided.",
                label
            ),
            Error::MissingEvaluation { label } => write!(
                f,
                "`QuerySet` refers to polynomial \"{}\", but `Evaluations` does not contain an evaluation for it.",
                label
            ),
            Error::MissingLHS { label } => {
                write!(f, "Equation \"{}\" does not have a LHS.", label)
            },
            Error::MissingRng => write!(f, "hiding commitments require `Some(rng)`"),
            Error::DegreeIsZero => write!(
                f,
                "this scheme does not support committing to degree 0 polynomials"
            ),
            Error::TooManyCoefficients {
                num_coefficients,
                num_powers,
            } => write!(
                f,
                "the number of coefficients in the polynomial ({:?}) is greater than\
                 the maximum number of powers in `Powers` ({:?})",
                num_coefficients, num_powers
            ),
            Error::HidingBoundIsZero => write!(
                f,
                "this scheme does not support non-`None` hiding bounds that are 0"
            ),
            Error::HidingBoundToolarge {
                hiding_poly_degree,
                num_powers,
            } => write!(
                f,
                "the degree of the hiding poly ({:?}) is not less than the maximum number of powers in `Powers` ({:?})",
                hiding_poly_degree, num_powers
            ),
            Error::TrimmingDegreeTooLarge => {
                write!(f, "the degree provided to `trim` was too large")
            }
            Error::EmptyDegreeBounds => {
                write!(f, "provided `enforced_degree_bounds` was `Some<&[]>`")
            }
            Error::EquationHasDegreeBounds(e) => write!(
                f,
                "the eqaution \"{}\" contained degree-bounded polynomials",
                e
            ),
            Error::UnsupportedDegreeBound(bound) => write!(
                f,
                "the degree bound ({:?}) is not supported by the parameters",
                bound,
            ),
            Error::IncorrectDegreeBound {
                poly_degree,
                degree_bound,
                supported_degree,
                label,
            } => write!(
                f,
                "the degree bound ({:?}) for the polynomial {} \
                 (having degree {:?}) is greater than the maximum \
                 supported degree ({:?})",
                degree_bound, label, poly_degree, supported_degree
            ),
            Error::InvalidNumberOfVariables => write!(
                f,
                "An invalid number of variables was provided to `setup`"
            ),
            Error::PolynomialDegreeTooLarge {
                poly_degree,
                supported_degree,
                label,
            } => write!(
                f,
                "the polynomial {} has degree {:?}, but parameters only
                support up to degree ({:?})", label, poly_degree, supported_degree
            ),
            Error::IncorrectInputLength(err) => write!(f, "{}", err),
            Error::InvalidCommitment => write!(f, "Failed to verify the commitment"),
            Error::TranscriptError => write!(f, "Incorrect transcript manipulation"),
            Error::InvalidParameters(err) => write!(f, "{}", err),
            Error::HashingError => write!(f, "Error resulting from hashing"),
            Error::EncodingError => write!(f, "Encoding failed"),
        }
    }
}

impl ark_std::error::Error for Error {}
