
#[cfg(test)]
mod tests {

    use ark_bls12_381::Fq as F;
    use ark_poly::{domain::general::GeneralEvaluationDomain, EvaluationDomain, DenseUVPolynomial, univariate::DensePolynomial, Polynomial};

    use crate::ligero::utils::*;

    #[test]
    fn test_matrix_constructor_flat() {
        let entries: Vec<F> = to_field(vec![10, 100, 4, 67, 44, 50]);
        let mat = Matrix::new_from_flat(2, 3, &entries);
        assert_eq!(mat.entry(1, 2), F::from(50));
    }

    #[test]
    #[should_panic]
    fn test_matrix_constructor_flat_panic() {
        let entries: Vec<F> = to_field(vec![10, 100, 4, 67, 44]);
        Matrix::new_from_flat(2, 3, &entries);
    }

    #[test]
    fn test_matrix_constructor_rows() {
        let rows: Vec<Vec<F>> = vec!(
            to_field(vec![10, 100, 4]),
            to_field(vec![23, 1, 0]),
            to_field(vec![55, 58, 9]),
        );
        let mat = Matrix::new_from_rows(rows);
        assert_eq!(mat.entry(2, 0), F::from(55));
    }

    #[test]
    #[should_panic]
    fn test_matrix_constructor_rows_panic() {
        let rows: Vec<Vec<F>> = vec!(
            to_field(vec![10, 100, 4]),
            to_field(vec![23, 1, 0]),
            to_field(vec![55, 58]),
        );
        Matrix::new_from_rows(rows);
    }

    #[test]
    fn test_cols() {
        let rows: Vec<Vec<F>> = vec!(
            to_field(vec![4, 76]),
            to_field(vec![14, 92,]),
            to_field(vec![17, 89]),
        );

        let mat = Matrix::new_from_rows(rows);

        assert_eq!(mat.cols()[1], to_field(vec![76, 92, 89]));
    }

    #[test]
    fn test_row_mul() {
        let rows: Vec<Vec<F>> = vec!(
            to_field(vec![10, 100, 4]),
            to_field(vec![23, 1, 0]),
            to_field(vec![55, 58, 9]),
        );

        let mat = Matrix::new_from_rows(rows);
        let v: Vec<F> = to_field(vec![12, 41, 55]);
        // by giving the result in the integers and then converting to F
        // we ensure the test will still pass even if F changes
        assert_eq!(mat.row_mul(&v), to_field::<F>(vec![4088, 4431, 543]));
    }

    #[test]
    fn test_fft_interface() {

        // we use this polynomial to generate the the values we will ask the fft to interpolate
        let pol_coeffs: Vec<F> = to_field(vec![30, 2, 91]);
        let pol: DensePolynomial<F> = DensePolynomial::from_coefficients_slice(&pol_coeffs);
        
        let fft_domain = GeneralEvaluationDomain::<F>::new(pol_coeffs.len()).unwrap();

        // generating the values
        let mut vals = Vec::new();

        for i in 0..4 {
            vals.push(pol.evaluate(&fft_domain.element(i)));
        }

        // the fft should recover the original polynomial
        assert_eq!(fft_domain.ifft(&vals), pol_coeffs);
    }

    #[test]
    fn test_get_num_bytes() {
        assert_eq!(get_num_bytes(0), 0);
        assert_eq!(get_num_bytes(1), 1);
        assert_eq!(get_num_bytes(9), 1);
        assert_eq!(get_num_bytes(1 << 11), 2);
        assert_eq!(get_num_bytes(1 << 32 - 1), 4);
        assert_eq!(get_num_bytes(1 << 32), 5);
        assert_eq!(get_num_bytes(1 << 32 + 1), 5);
    }
}