use ark_ff::Field;
#[cfg(not(feature = "std"))]
use ark_std::vec::Vec;

#[cfg(feature = "parallel")]
use rayon::prelude::*;

/// Transforms a flat vector into a n*m matrix in column-major order. The
/// latter is given as a list of rows.
///
/// For example, if flat = [1, 2, 3, 4, 5, 6] and n = 3, m = 2, then
/// the output is [[1, 3, 5], [2, 4, 6]].
pub(crate) fn flat_to_matrix_column_major<T: Copy>(flat: &[T], n: usize, m: usize) -> Vec<Vec<T>> {
    assert_eq!(flat.len(), n * m, "n * m should coincide with flat.len()");
    let mut res = Vec::new();

    for row in 0..n {
        res.push((0..m).map(|col| flat[col * n + row]).collect())
    }
    res
}

// This function computes all evaluations of the MLE EQ(i, values) for i
// between 0...0 and 1...1 (n-bit strings). This results in essentially
// the same as the tensor_vec function in the `linear_codes/utils.rs`,
// the difference being the endianness of the order of the output.
pub(crate) fn tensor_prime<F: Field>(values: &[F]) -> Vec<F> {
    if values.is_empty() {
        return vec![F::one()];
    }

    let tail = tensor_prime(&values[1..]);
    let val = values[0];

    cfg_iter!(tail)
        .map(|v| *v * (F::one() - val))
        .chain(cfg_iter!(tail).map(|v| *v * val))
        .collect()
}
