use crate::String;
/// Errors that arise when dealing with query sets.
#[derive(Debug)]
pub enum QuerySetError {
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
}

impl core::fmt::Display for QuerySetError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
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

impl algebra_core::Error for QuerySetError {}

/// Equation errors that arise when dealing with query sets.
#[derive(Debug)]
pub enum EquationError {
    /// The LHS of the equation is empty.
    MissingLHS {
        /// The label of the equation.
        label: String,
    },
}

impl core::fmt::Display for EquationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            EquationError::MissingLHS { label } => {
                write!(f, "Equation \"{}\" does not have a LHS.", label)
            }
        }
    }
}

impl algebra_core::Error for EquationError {}
