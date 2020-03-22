use crate::Error;
use crate::{
    EquationError as EqError, LabeledPolynomial, PCCommitterKey, QuerySetError as QSError,
};
use crate::{String, ToString};

pub(crate) fn check_degree_is_within_bounds(
    num_coefficients: usize,
    num_powers: usize,
) -> Result<(), Error> {
    if num_coefficients < 1 {
        Err(Error::DegreeIsZero)
    } else {
        check_degree_is_too_large(num_coefficients, num_powers)
    }
}

pub(crate) fn check_degree_is_too_large(
    num_coefficients: usize,
    num_powers: usize,
) -> Result<(), Error> {
    if num_coefficients > num_powers {
        Err(Error::TooManyCoefficients {
            num_coefficients,
            num_powers,
        })
    } else {
        Ok(())
    }
}

pub(crate) fn check_hiding_bound(
    hiding_poly_degree: usize,
    num_powers: usize,
) -> Result<(), Error> {
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

pub(crate) fn check_degrees_and_bounds<'a, E: algebra_core::PairingEngine>(
    supported_degree: usize,
    max_degree: usize,
    enforced_degree_bounds: Option<&[usize]>,
    p: &'a LabeledPolynomial<'a, E::Fr>,
) -> Result<(), Error> {
    if let Some(bound) = p.degree_bound() {
        let enforced_degree_bounds =
            enforced_degree_bounds.ok_or(Error::UnsupportedDegreeBound(bound))?;

        if enforced_degree_bounds.binary_search(&bound).is_err() {
            Err(Error::UnsupportedDegreeBound(bound))
        } else if bound < p.degree() || bound > max_degree {
            return Err(Error::IncorrectDegreeBound {
                poly_degree: p.degree(),
                degree_bound: p.degree_bound().unwrap(),
                supported_degree,
                label: p.label().to_string(),
            });
        } else {
            Ok(())
        }
    } else {
        Ok(())
    }
}
