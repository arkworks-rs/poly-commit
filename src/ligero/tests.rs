
#[cfg(test)]
mod tests {

    use ark_bls12_381::Fq as F;

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

}