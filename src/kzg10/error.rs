/// The error type for `KZG10`.
#[derive(Debug)]
pub enum Error {
    /// The provided polynomial was meant to be hiding, but `rng` was `None`.
    MissingRng,
    /// The degree provided in setup was too small; degree 0 polynomials
    /// are not supported.
    DegreeIsZero,
    /// The degree of the polynomial passed to `KZG10::commit` or `KZG10::open`
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
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
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
            )
        }
    }
}

impl std::error::Error for Error {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        None
    }
}

impl Error {
    pub(crate) fn check_degree_is_within_bounds(num_coefficients: usize, num_powers: usize) -> Result<(), Self> {
        if num_coefficients < 1 {
            Err(Error::DegreeIsZero)
        } else {
            Self::check_degree_is_too_large(num_coefficients, num_powers)
        }
    }


    pub(crate) fn check_degree_is_too_large(num_coefficients: usize, num_powers: usize) -> Result<(), Self> {
        if num_coefficients > num_powers {
            Err(Error::TooManyCoefficients {
                num_coefficients,
                num_powers,
            })
        } else {
            Ok(())
        }
    }

    pub(crate) fn check_hiding_bound(hiding_poly_degree: usize, num_powers: usize) -> Result<(), Self> {
        if hiding_poly_degree == 0 {
            Err(Error::HidingBoundIsZero)
        } else if hiding_poly_degree >= num_powers {
            // The above check uses `>=` because committing to a hiding poly with
            // degree `hiding_poly_degree` requires `hiding_poly_degree + 1` 
            // powers.
            Err(Error::HidingBoundToolarge {
                hiding_poly_degree,
                num_powers,
            })
        } else {
            Ok(())
        }
    }
}
