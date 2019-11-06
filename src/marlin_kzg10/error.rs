use crate::kzg10;

/// Error type for `MultiPCFromSinglePC`.
#[derive(Debug)]
pub enum Error {
    /// The degree of the `index`-th polynomial passed to `commit` or `open`
    /// was too large.
    PolynomialDegreeTooLarge {
        /// Degree of the polynomial.
        poly_degree: usize,
        /// Maximum supported degree.
        max_degree: usize,
        /// Index of the offending polynomial.
        label: String,
    },
    /// The degree bound for the `index`-th polynomial passed to `commit`, `open`
    /// or `check` was incorrect, that is, `degree_bound >= poly_degree` or
    /// `degree_bound <= max_degree`.
    IncorrectDegreeBound {
        /// Degree of the polynomial.
        poly_degree: usize,
        /// Degree bound.
        degree_bound: usize,
        /// Maximum supported degree.
        max_degree: usize,
        /// Index of the offending polynomial.
        label: String,
    },
    /// The inputs to `commit`, `open` or `verify` had incorrect lengths.
    IncorrectInputLength(String),
    /// The query set referenced a non-existent polynomial.
    IncorrectQuerySet(&'static str),
    /// The evaluations referenced a non-existent member of the query set.
    IncorrectEvaluation(&'static str),
    /// An error from the underlying `SinglePC`.
    KZG10Error(kzg10::Error),
}

impl From<kzg10::Error> for Error {
    fn from(other: kzg10::Error) -> Self {
        Error::KZG10Error(other)
    }
}

impl Error {
    pub(crate) fn poly_degree_too_large(poly_degree: usize, max_degree: usize, label: String) -> Self {
        Error::PolynomialDegreeTooLarge {
            poly_degree,
            max_degree,
            label,
        }
    }

    pub(crate) fn incorrect_bound(
        poly_degree: usize,
        degree_bound: usize,
        max_degree: usize,
        label: String,
    ) -> Self {
        Error::IncorrectDegreeBound {
            poly_degree,
            degree_bound,
            max_degree,
            label,
        }
    }

    pub(crate) fn check_degrees(
        d: usize,
        bound: Option<usize>,
        max_degree: usize,
        label: String,
    ) -> Result<(), Self> {
        if let Some(bound) = bound {
            if d > max_degree {
                Err(Error::poly_degree_too_large(d, max_degree, label))
            } else if bound < d || bound > max_degree {
                Err(Error::incorrect_bound(d, bound, max_degree, label))
            } else {
                Ok(())
            }
        } else {
            Ok(())
        }
    }
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Error::PolynomialDegreeTooLarge {
                poly_degree,
                max_degree,
                label,
            } => write!(
                f,
                "the degree of the polynomial {} ({:?}) is greater than\
                 the maximum supported degree ({:?})",
                label, poly_degree, max_degree
            ),
            Error::IncorrectDegreeBound {
                poly_degree,
                degree_bound,
                max_degree,
                label,
            } => write!(
                f,
                "the degree bound ({:?}) for the polynomial {} \
                 (having degree {:?}) is greater than the maximum \
                 supported degree ({:?})",
                degree_bound, label, poly_degree, max_degree
            ),
            Error::IncorrectInputLength(err) => write!(f, "{}", err),
            Error::IncorrectQuerySet(err) => write!(f, "{}", err),
            Error::IncorrectEvaluation(err) => write!(f, "{}", err),
            Error::KZG10Error(err) => write!(f, "KZG10 error: {}", err),
        }
    }
}

impl std::error::Error for Error {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        None
    }
}
