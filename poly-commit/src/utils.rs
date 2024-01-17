#[cfg(not(feature = "std"))]
use num_traits::Float;

#[cfg(feature = "parallel")]
use rayon::{
    iter::{IntoParallelRefIterator, ParallelIterator},
    prelude::IndexedParallelIterator,
};

use ark_ff::Field;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::vec::Vec;

/// Takes as input a struct, and converts them to a series of bytes. All traits
/// that implement `CanonicalSerialize` can be automatically converted to bytes
/// in this manner.
/// From jellyfish lib
#[macro_export]
macro_rules! to_bytes {
    ($x:expr) => {{
        let mut buf = ark_std::vec![];
        ark_serialize::CanonicalSerialize::serialize_compressed($x, &mut buf).map(|_| buf)
    }};
}

/// Entropy function
pub(crate) fn ent(x: f64) -> f64 {
    assert!(0f64 <= x && x <= 1f64);
    if x == 0f64 || x == 1f64 {
        0f64
    } else {
        -x * x.log2() - (1.0 - x) * (1.0 - x).log2()
    }
}

/// ceil of a * b, where a is integer and b is a rational number
#[inline]
pub(crate) fn ceil_mul(a: usize, b: (usize, usize)) -> usize {
    (a * b.0 + b.1 - 1) / b.1
}

/// Return ceil(x / y).
pub(crate) fn ceil_div(x: usize, y: usize) -> usize {
    // XXX. warning: this expression can overflow.
    (x + y - 1) / y
}

#[derive(Derivative, CanonicalSerialize, CanonicalDeserialize)]
#[derivative(Default(bound = ""), Clone(bound = ""), Debug(bound = ""))]
pub struct Matrix<F: Field> {
    pub(crate) n: usize,
    pub(crate) m: usize,
    entries: Vec<Vec<F>>,
}

impl<F: Field> Matrix<F> {
    /// Returns a Matrix of dimensions n x m given a list of n * m field elements.
    /// The list should be ordered row-first, i.e. [a11, ..., a1m, a21, ..., a2m, ...].
    ///
    /// # Panics
    /// Panics if the dimensions do not match the length of the list
    pub(crate) fn new_from_flat(n: usize, m: usize, entry_list: &[F]) -> Self {
        assert_eq!(
            entry_list.len(),
            n * m,
            "Invalid matrix construction: dimensions are {} x {} but entry vector has {} entries",
            n,
            m,
            entry_list.len()
        );

        // TODO more efficient to run linearly?
        let entries: Vec<Vec<F>> = (0..n)
            .map(|row| (0..m).map(|col| entry_list[m * row + col]).collect())
            .collect();

        Self { n, m, entries }
    }

    /// Returns a Matrix given a list of its rows, each in turn represented as a list of field elements.
    ///
    /// # Panics
    /// Panics if the sub-lists do not all have the same length.
    pub(crate) fn new_from_rows(row_list: Vec<Vec<F>>) -> Self {
        let m = row_list[0].len();

        for row in row_list.iter().skip(1) {
            assert_eq!(
                row.len(),
                m,
                "Invalid matrix construction: not all rows have the same length"
            );
        }

        Self {
            n: row_list.len(),
            m,
            entries: row_list,
        }
    }

    /// Returns the entry in position (i, j). **Indexing starts at 0 in both coordinates**,
    /// i.e. the first element is in position (0, 0) and the last one in (n - 1, j - 1),
    /// where n and m are the number of rows and columns, respectively.
    ///
    /// Index bound checks are waived for efficiency and behaviour under invalid indexing is undefined
    #[cfg(test)]
    pub(crate) fn entry(&self, i: usize, j: usize) -> F {
        self.entries[i][j]
    }

    /// Returns self as a list of rows
    pub(crate) fn rows(&self) -> Vec<Vec<F>> {
        self.entries.clone()
    }

    /// Returns self as a list of columns
    pub(crate) fn cols(&self) -> Vec<Vec<F>> {
        (0..self.m)
            .map(|col| (0..self.n).map(|row| self.entries[row][col]).collect())
            .collect()
    }

    /// Returns the product v * self, where v is interpreted as a row vector. In other words,
    /// it returns a linear combination of the rows of self with coefficients given by v.
    ///
    /// Panics if the length of v is different from the number of rows of self.
    pub(crate) fn row_mul(&self, v: &[F]) -> Vec<F> {
        assert_eq!(
            v.len(),
            self.n,
            "Invalid row multiplication: vector has {} elements whereas each matrix column has {}",
            v.len(),
            self.n
        );

        (0..self.m)
            .map(|col| {
                inner_product(
                    v,
                    &(0..self.n)
                        .map(|row| self.entries[row][col])
                        .collect::<Vec<F>>(),
                )
            })
            .collect()
    }
}

#[inline]
pub(crate) fn inner_product<F: Field>(v1: &[F], v2: &[F]) -> F {
    ark_std::cfg_iter!(v1)
        .zip(v2)
        .map(|(li, ri)| *li * ri)
        .sum()
}

#[inline]
#[cfg(test)]
pub(crate) fn to_field<F: Field>(v: Vec<u64>) -> Vec<F> {
    v.iter().map(|x| F::from(*x)).collect::<Vec<F>>()
}

// TODO: replace by https://github.com/arkworks-rs/crypto-primitives/issues/112.
#[cfg(test)]
use ark_crypto_primitives::sponge::poseidon::PoseidonSponge;
#[cfg(test)]
use ark_ff::PrimeField;

#[cfg(test)]
pub(crate) fn test_sponge<F: PrimeField>() -> PoseidonSponge<F> {
    use ark_crypto_primitives::sponge::{poseidon::PoseidonConfig, CryptographicSponge};
    use ark_std::test_rng;

    let full_rounds = 8;
    let partial_rounds = 31;
    let alpha = 17;

    let mds = vec![
        vec![F::one(), F::zero(), F::one()],
        vec![F::one(), F::one(), F::zero()],
        vec![F::zero(), F::one(), F::one()],
    ];

    let mut v = Vec::new();
    let mut ark_rng = test_rng();

    for _ in 0..(full_rounds + partial_rounds) {
        let mut res = Vec::new();

        for _ in 0..3 {
            res.push(F::rand(&mut ark_rng));
        }
        v.push(res);
    }
    let config = PoseidonConfig::new(full_rounds, partial_rounds, alpha, mds, v, 2, 1);
    PoseidonSponge::new(&config)
}

#[cfg(test)]
pub(crate) mod tests {

    use super::*;

    use ark_bls12_377::Fr;

    #[test]
    fn test_matrix_constructor_flat() {
        let entries: Vec<Fr> = to_field(vec![10, 100, 4, 67, 44, 50]);
        let mat = Matrix::new_from_flat(2, 3, &entries);
        assert_eq!(mat.entry(1, 2), Fr::from(50));
    }

    #[test]
    fn test_matrix_constructor_flat_square() {
        let entries: Vec<Fr> = to_field(vec![10, 100, 4, 67]);
        let mat = Matrix::new_from_flat(2, 2, &entries);
        assert_eq!(mat.entry(1, 1), Fr::from(67));
    }

    #[test]
    #[should_panic(expected = "dimensions are 2 x 3 but entry vector has 5 entries")]
    fn test_matrix_constructor_flat_panic() {
        let entries: Vec<Fr> = to_field(vec![10, 100, 4, 67, 44]);
        Matrix::new_from_flat(2, 3, &entries);
    }

    #[test]
    fn test_matrix_constructor_rows() {
        let rows: Vec<Vec<Fr>> = vec![
            to_field(vec![10, 100, 4]),
            to_field(vec![23, 1, 0]),
            to_field(vec![55, 58, 9]),
        ];
        let mat = Matrix::new_from_rows(rows);
        assert_eq!(mat.entry(2, 0), Fr::from(55));
    }

    #[test]
    #[should_panic(expected = "not all rows have the same length")]
    fn test_matrix_constructor_rows_panic() {
        let rows: Vec<Vec<Fr>> = vec![
            to_field(vec![10, 100, 4]),
            to_field(vec![23, 1, 0]),
            to_field(vec![55, 58]),
        ];
        Matrix::new_from_rows(rows);
    }

    #[test]
    fn test_cols() {
        let rows: Vec<Vec<Fr>> = vec![
            to_field(vec![4, 76]),
            to_field(vec![14, 92]),
            to_field(vec![17, 89]),
        ];

        let mat = Matrix::new_from_rows(rows);

        assert_eq!(mat.cols()[1], to_field(vec![76, 92, 89]));
    }

    #[test]
    fn test_row_mul() {
        let rows: Vec<Vec<Fr>> = vec![
            to_field(vec![10, 100, 4]),
            to_field(vec![23, 1, 0]),
            to_field(vec![55, 58, 9]),
        ];

        let mat = Matrix::new_from_rows(rows);
        let v: Vec<Fr> = to_field(vec![12, 41, 55]);
        // by giving the result in the integers and then converting to Fr
        // we ensure the test will still pass even if Fr changes
        assert_eq!(mat.row_mul(&v), to_field::<Fr>(vec![4088, 4431, 543]));
    }
}
