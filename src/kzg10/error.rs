/// The error type for `KZG10`.
#[derive(Debug)]
pub enum Error {
    /// The provided polynomial was meant to be hiding, but `rng` was `None`.
    MissingRng,
    /// The degree provided in setup was too small; degree 0 polynomials
    /// are not supported.
    UnsupportedDegree,
    /// The degree of the polynomial passed to `KZG10::commit` or `KZG10::open`
    /// was too large.
    PolynomialDegreeTooLarge {
        /// The degree of the polynomial.
        poly_degree: usize,
        /// The maximum supported degree.
        max_degree: usize,
    },
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Error::MissingRng => write!(f, "hiding commitments require `Some(rng)`"),
            Error::UnsupportedDegree => write!(
                f,
                "this scheme does not support committing to degree 0 polynomials"
            ),
            Error::PolynomialDegreeTooLarge {
                poly_degree,
                max_degree,
            } => write!(
                f,
                "the degree of the polynomial ({:?}) is greater than\
                 the maximum supported degree ({:?})",
                poly_degree, max_degree
            ),
        }
    }
}

impl std::error::Error for Error {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        None
    }
}

impl Error {
    fn check_degree(d: usize, max_degree: usize) -> Result<(), Self> {
        if d < 1 {
            Err(Error::UnsupportedDegree)
        } else if d > max_degree {
            Err(Error::PolynomialDegreeTooLarge {
                poly_degree: d,
                max_degree,
            })
        } else {
            Ok(())
        }
    }
}


