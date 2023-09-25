use ark_ff::{FftField, Field, PrimeField};

use ark_poly::{EvaluationDomain, GeneralEvaluationDomain};
use ark_serialize::CanonicalSerialize;
use digest::Digest;
use jf_primitives::pcs::transcript::IOPTranscript;
use rayon::{
    iter::{IntoParallelRefIterator, ParallelIterator},
    prelude::IndexedParallelIterator,
};

use crate::streaming_kzg::ceil_div;
use crate::Error;

#[derive(Debug)]
pub(crate) struct Matrix<F: Field> {
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
            "Invalid row multiplication: vectir has {} elements whereas each matrix column has {}",
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

/// Compute the dimensions of an FFT-friendly (over F) matrix with at least n entries.
/// The return pair (n, m) corresponds to the dimensions n x m.
pub(crate) fn compute_dimensions<F: FftField>(n: usize) -> (usize, usize) {
    assert_eq!(
        (n as f64) as usize,
        n,
        "n cannot be converted to f64: aborting"
    );

    let aux = (n as f64).sqrt().ceil() as usize;
    let n_cols = GeneralEvaluationDomain::<F>::new(aux)
        .expect("Field F does not admit FFT with m elements")
        .size();

    (ceil_div(n, n_cols), n_cols)
}

/// Apply reed-solomon encoding to msg.
/// Assumes msg.len() is equal to the order of an FFT domain in F.
/// Returns a vector of length equal to the smallest FFT domain of size at least msg.len() * RHO_INV.
pub(crate) fn reed_solomon<F: FftField>(
    // msg, of length m, is interpreted as a vector of coefficients of a polynomial of degree m - 1
    msg: &[F],
    rho_inv: usize,
) -> Vec<F> {
    let m = msg.len();

    let domain = GeneralEvaluationDomain::<F>::new(m).unwrap();
    assert_eq!(
        m,
        domain.size(),
        "The evaluation vector has length {} elements \\
        but the smallest FFT domain admitting that many elements has order {}",
        m,
        domain.size()
    );
    let poly_coeffs = domain.ifft(msg).to_vec();

    let extended_domain = GeneralEvaluationDomain::<F>::new(m * rho_inv).unwrap_or_else(|| {
        panic!(
            "The field F cannot accomodate FFT for msg.len() * RHO_INV = {} elements (too many)",
            m * rho_inv
        )
    });

    extended_domain.fft(&poly_coeffs)
}

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
pub(crate) fn hash_column<D: Digest, F: PrimeField + CanonicalSerialize>(array: &[F]) -> Vec<u8> {
    let mut dig = D::new();
    for elem in array {
        dig.update(to_bytes!(elem).unwrap());
    }
    dig.finalize().to_vec()
}

/// Generate `t` (not necessarily distinct) random points in `[0, n)` using the current state of `transcript`
pub(crate) fn get_indices_from_transcript<F: PrimeField>(
    n: usize,
    t: usize,
    transcript: &mut IOPTranscript<F>,
) -> Result<Vec<usize>, Error> {
    let bytes_to_squeeze = get_num_bytes(n);
    let mut indices = Vec::with_capacity(t);
    for _ in 0..t {
        let mut bytes: Vec<u8> = vec![0; bytes_to_squeeze];
        transcript
            .get_and_append_byte_challenge(b"i", &mut bytes)
            .map_err(|_| Error::TranscriptError)?;

        // get the usize from Vec<u8>:
        let ind = bytes.iter().fold(0, |acc, &x| (acc << 8) + x as usize);
        // modulo the number of columns in the encoded matrix
        indices.push(ind % n);
    }
    Ok(indices)
}

#[inline]
pub(crate) fn calculate_t(rho_inv: usize, sec_param: usize) -> usize {
    // Double-check with BCI+20 if you can simply replace
    // $\delta = \frac{1-\rho}{3}$ with $\frac{1-\rho}{2}$. In that case, we
    // will find the smallest t such that
    // $(1-\delta)^t + (\rho+\delta)^t + n/F < 2^(-\lambda)$. Since we do not
    // have $n/F$ here and and it is negligible for security less than 230 bits,
    // we eliminate it here. With \delta = \frac{1-\rho}{2}, the expreesion is
    // 2 * (1-\delta)^t < 2^(-\lambda). Then
    // $t * log2 (1-\delta) < - \lambda - 1$.

    // TODO: Maybe we should not eliminate $n/F$. In original Ligero, this was
    // $d/F$ for $\delta = \frac{1-\rho}{3}$. But from expression 1.1 from
    // BCI+20, I wrote $n/F$ for $\delta = \frac{1-\rho}{3}$. It is negligible
    // anyways.
    let sec_param = sec_param as f64;
    let nom = -sec_param - 1.0;
    let denom = (0.5 + 0.5 / rho_inv as f64).log2();
    (nom / denom).ceil() as usize // This is the `t`
}
