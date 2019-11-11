/// The error type for `KZG10`.
#[derive(Debug)]
pub enum QuerySetError {
    /// The query set contains a label for a polynomial that was not provided as 
    /// input to the `PC::open`.
    MissingPolynomial {
        /// The label of the missing polynomial.
        label: String
    },
    /// `Evaluations` does not contain an evaluation for the polynomial labelled 
    /// `label` at a particular query.
    MissingEvaluation {
        /// The label of the missing polynomial.
        label: String
    },
}

impl std::fmt::Display for QuerySetError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            QuerySetError::MissingPolynomial { label } => write!(
                f,
                "`QuerySet` refers to polynomial \"{}\", but it was not provided.",
                label
            ),
            QuerySetError::MissingEvaluation { label } => write!(
                f,
                "`QuerySet` refers to polynomial \"{}\", but `Evaluations` does not contain an evaluation for it.",
                label
            ),
        }
    }
}

impl std::error::Error for QuerySetError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        None
    }
}
