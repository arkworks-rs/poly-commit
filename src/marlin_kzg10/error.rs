use crate::kzg10;
use crate::{
    EquationError as EqError, LabeledPolynomial, PCCommitterKey, QuerySetError as QSError,
};
use crate::{String, ToString};

/// Error type for `MultiPCFromSinglePC`.
#[derive(Debug)]
pub enum Error {
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

    /// An error related to the `QuerySet`.
    QuerySetError(QSError),

    /// An error related to the `QuerySet`.
    EquationError(EqError),

    /// An error from the underlying `KZG10`.
    KZG10Error(kzg10::Error),
}

impl From<kzg10::Error> for Error {
    fn from(other: kzg10::Error) -> Self {
        Error::KZG10Error(other)
    }
}

impl From<QSError> for Error {
    fn from(other: QSError) -> Self {
        Error::QuerySetError(other)
    }
}

impl From<EqError> for Error {
    fn from(other: EqError) -> Self {
        Error::EquationError(other)
    }
}

impl Error {
    pub(crate) fn check_degrees_and_bounds<'a, E: algebra_core::PairingEngine>(
        ck: &super::CommitterKey<E>,
        p: &'a LabeledPolynomial<'a, E::Fr>,
    ) -> Result<(), Self> {
        if let Some(bound) = p.degree_bound() {
            let enforced_degree_bounds = ck
                .enforced_degree_bounds
                .as_ref()
                .ok_or(Self::UnsupportedDegreeBound(bound))?;

            if enforced_degree_bounds.binary_search(&bound).is_err() {
                Err(Self::UnsupportedDegreeBound(bound))
            } else if bound < p.degree() || bound > ck.max_degree() {
                return Err(Error::IncorrectDegreeBound {
                    poly_degree: p.degree(),
                    degree_bound: p.degree_bound().unwrap(),
                    supported_degree: ck.supported_degree(),
                    label: p.label().to_string(),
                });
            } else {
                Ok(())
            }
        } else {
            Ok(())
        }
    }
}

impl core::fmt::Display for Error {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
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
            Error::IncorrectInputLength(err) => write!(f, "{}", err),
            Error::QuerySetError(err) => write!(f, "{}", err),
            Error::EquationError(err) => write!(f, "{}", err),
            Error::KZG10Error(err) => write!(f, "KZG10 error: {}", err),
        }
    }
}

impl algebra_core::Error for Error {}
