use ark_ff::{FftField, Field, PrimeField};

use ark_poly::{
    domain::general::GeneralElements, univariate::DensePolynomial, DenseUVPolynomial,
    EvaluationDomain, GeneralEvaluationDomain, Polynomial,
};
use ark_serialize::CanonicalSerialize;
use digest::Digest;
use rayon::{
    iter::{IntoParallelRefIterator, ParallelIterator},
    prelude::IndexedParallelIterator,
};

use crate::streaming_kzg::ceil_div;

pub(crate) struct Matrix<F: Field> {
    n: usize,
    m: usize,
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
            m,
            n,
            entry_list.len()
        );

        // TODO more efficient to run linearly?
        let entries: Vec<Vec<F>> = (0..n)
            .map(|row| (0..m).map(|col| entry_list[n * row + m]).collect())
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
    pub(crate) fn entry(&self, i: usize, j: usize) -> F {
        return self.entries[i][j];
    }

    /// Returns self as a list of rows
    pub(crate) fn rows(&self) -> Vec<Vec<F>> {
        return self.entries.clone();
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
            "Invalid row multiplication x has {} elements whereas the matrix has {}",
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

// TODO batch for all rows? possibly minimal savings due to not having to create fft_domain every time
// here we restrict F to PrimeField
pub(crate) fn reed_solomon<F: FftField>(
    msg: &[F],
    rho_inverse: usize,
    fft_domain: GeneralEvaluationDomain<F>,
    domain_iter: &mut GeneralElements<F>,
) -> Vec<F> {
    // TODO is this check worth it?
    // rho_inverse = 0 should never happen; rho_inverse = 1 means no expansion
    if rho_inverse <= 1 {
        return msg.to_vec();
    }

    let m = msg.len();

    let pol: DensePolynomial<F> = DensePolynomial::from_coefficients_slice(&fft_domain.ifft(msg));

    let mut encoding = msg.to_vec();

    domain_iter.nth(m - 1);

    for _ in 0..(rho_inverse - 1) * m {
        // TODO make sure the domain has enough elements in the caller or check here
        let zeta = domain_iter.next().unwrap();
        encoding.push(pol.evaluate(&zeta));
    }

    encoding
}

/* DummyCK<F> {
    t: int;
}

impl<F: PrimeField> DummyCK<F> {
    fn new() -> Self {
        println!("WARNING: You are using dummy parameters"),

    }
}
 */
#[inline]
pub(crate) fn inner_product<F: Field>(v1: &[F], v2: &[F]) -> F {
    ark_std::cfg_iter!(v1)
        .zip(v2)
        .map(|(li, ri)| *li * ri)
        .sum()
}

#[inline]
pub(crate) fn to_field<F: Field>(v: Vec<u64>) -> Vec<F> {
    v.iter().map(|x| F::from(*x)).collect::<Vec<F>>()
}

#[inline]
pub(crate) fn get_num_bytes(n: usize) -> usize {
    ceil_div((usize::BITS - n.leading_zeros()) as usize, 8)
}

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

#[inline]
pub(crate) fn hash_array<D: Digest, F: PrimeField + CanonicalSerialize>(array: &[F]) -> Vec<u8> {

    let mut dig = D::new();
    for elem in array {
        dig.update(to_bytes!(elem).unwrap());
    }    
    dig.finalize().to_vec()
}
