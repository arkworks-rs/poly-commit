use ark_ff::PrimeField;

use rayon::{iter::{ParallelIterator, IntoParallelRefIterator}, prelude::IndexedParallelIterator};

pub(crate) struct Matrix<F: PrimeField> {
    n: usize,
    m: usize,
    entries: Vec<Vec<F>>
}

impl<F: PrimeField> Matrix<F> {
    
    /// Returns a Matrix of dimensions n x m given a list of n * m field elements.
    /// The list should be ordered row-first, i.e. [a11, ..., a1m, a21, ..., a2m, ...].
    ///
    /// # Panics
    /// Panics if the dimensions do not match the length of the list
    pub(crate) fn new_from_flat(n: usize, m: usize, entry_list: &[F]) -> Self {
        assert_eq!(entry_list.len(), n * m, "Invalid matrix construction: dimensions are {} x {} but entry vector has {} entries", m, n, entry_list.len());

        // TODO more efficient to run linearly?
        let entries: Vec<Vec<F>> = (0..n).map(|row| (0..m).map(|col| entry_list[n * row + m]).collect()).collect();

        Self {n, m, entries}
    }

    /// Returns a Matrix given a list of its rows, each in turn represented as a list of field elements.
    ///
    /// # Panics
    /// Panics if the sub-lists do not all have the same length.
    pub(crate) fn new_from_row_list(n: usize, m: usize, row_list: Vec<Vec<F>>) -> Self {
        let mut m = row_list[0].len();

        for row in row_list {
            assert_eq!(row.len(), m, "Invalid matrix construction: not all rows have the same length");
            m = row.len()
        }

        Self {n: row_list.len(), m, entries: row_list}
    }

    // Returns the number of rows of self.
    pub(crate) fn num_rows(&self) -> usize {
        self.n
    }

    // Returns the number of columns of self.
    pub(crate) fn num_cols(&self) -> usize {
        self.m
    }

    /// Returns the entry in position (i, j). **Indexing starts at 0 in both coordinates**,
    /// i.e. the first element is in position (0, 0) and the last one in (n - 1, j - 1),
    /// where n and m are the number of rows and columns, respectively.
    /// 
    /// Index bound checks are waived for efficiency and behaviour under invalid indexing is undefined
    pub(crate) fn entry(&self, i: usize, j: usize) -> F {
        return self.entries[i][j]
    }

    /// Returns the product x * self, where x is interpreted as a row vector. In other words,
    /// it returns a linear combination of the rows of self with coefficients given by f.
    ///
    /// Panics if the length of x is different from the number of rows of self.
    pub(crate) fn row_mul(&self, x: &[F]) -> Vec<F> {
        assert_eq!(x.len(), self.n, "Invalid row multiplication x has {} elements whereas the matrix has {}", x.len(), self.n);
        
        (0..self.m).map(|col| inner_product(
            x, &(0..self.n).map(|row| self.entries[row][col]).collect::<Vec<F>>()
        )).collect()
    }

}

#[inline]
fn inner_product<F: PrimeField>(v1: &[F], v2: &[F]) -> F {
    ark_std::cfg_iter!(v1).zip(v2).map(|(li, ri)| *li * ri).sum()
}