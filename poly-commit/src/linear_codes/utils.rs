use crate::{utils::ceil_div, Error};
use ark_crypto_primitives::sponge::CryptographicSponge;
use ark_ff::{FftField, Field, PrimeField};
use ark_poly::{EvaluationDomain, GeneralEvaluationDomain};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
#[cfg(not(feature = "std"))]
use ark_std::{string::ToString, vec::Vec};

#[cfg(all(not(feature = "std"), target_arch = "aarch64"))]
use num_traits::Float;

#[cfg(test)]
use {
    crate::to_bytes,
    ark_crypto_primitives::crh::CRHScheme,
    ark_std::{borrow::Borrow, marker::PhantomData, rand::RngCore},
    digest::Digest,
};

/// This is CSC format
/// https://en.wikipedia.org/wiki/Sparse_matrix#Compressed_sparse_column_(CSC_or_CCS)
#[derive(Derivative, CanonicalSerialize, CanonicalDeserialize)]
#[derivative(Clone(bound = ""), Debug(bound = ""))]
pub struct SprsMat<F: Field> {
    /// Number of rows.
    pub(crate) n: usize,
    /// Number of columns.
    pub(crate) m: usize,
    /// Number of non-zero entries in each row.
    pub(crate) d: usize,
    /// Numbers of non-zero elements in each columns.
    ind_ptr: Vec<usize>,
    /// The indices in each columns where exists a non-zero element.
    col_ind: Vec<usize>,
    // The values of non-zero entries.
    val: Vec<F>,
}

impl<F: Field> SprsMat<F> {
    /// Calulates v.M
    pub(crate) fn row_mul(&self, v: &[F]) -> Vec<F> {
        (0..self.m)
            .map(|j| {
                let ij = self.ind_ptr[j]..self.ind_ptr[j + 1];
                self.col_ind[ij.clone()]
                    .iter()
                    .zip(&self.val[ij])
                    .map(|(&idx, x)| v[idx] * x)
                    .sum::<F>()
            })
            .collect::<Vec<_>>()
    }
    /// Create a new `SprsMat` from list of elements that represents the
    /// matrix in column major order. `n` is the number of rows, `m` is
    /// the number of columns, and `d` is NNZ in each row.
    pub fn new_from_flat(n: usize, m: usize, d: usize, list: &[F]) -> Self {
        let nnz = d * n;
        let mut ind_ptr = vec![0; m + 1];
        let mut col_ind = Vec::<usize>::with_capacity(nnz);
        let mut val = Vec::<F>::with_capacity(nnz);
        assert!(list.len() == m * n, "The dimension is incorrect.");
        for i in 0..m {
            for (c, &v) in list[i * n..(i + 1) * n].iter().enumerate() {
                if v != F::zero() {
                    ind_ptr[i + 1] += 1;
                    col_ind.push(c);
                    val.push(v);
                }
            }
            ind_ptr[i + 1] += ind_ptr[i];
        }
        assert!(ind_ptr[m] <= nnz, "The dimension or NNZ is incorrect.");
        Self {
            n,
            m,
            d,
            ind_ptr,
            col_ind,
            val,
        }
    }
    pub fn new_from_columns(n: usize, m: usize, d: usize, list: &[Vec<(usize, F)>]) -> Self {
        let nnz = d * n;
        let mut ind_ptr = vec![0; m + 1];
        let mut col_ind = Vec::<usize>::with_capacity(nnz);
        let mut val = Vec::<F>::with_capacity(nnz);
        assert!(list.len() == m, "The dimension is incorrect.");
        for j in 0..m {
            for (i, v) in list[j].iter() {
                ind_ptr[j + 1] += 1;
                col_ind.push(*i);
                val.push(*v);
            }
            assert!(list[j].len() <= n, "The dimension is incorrect.");
            ind_ptr[j + 1] += ind_ptr[j];
        }
        assert!(ind_ptr[m] <= nnz, "The dimension or NNZ is incorrect.");
        Self {
            n,
            m,
            d,
            ind_ptr,
            col_ind,
            val,
        }
    }
}

/// Apply reed-solomon encoding to msg.
/// Assumes msg.len() is equal to the order of some FFT domain in F.
/// Returns a vector of length equal to the smallest FFT domain of size at least msg.len() * RHO_INV.
pub(crate) fn reed_solomon<F: FftField>(
    // msg, of length m, is interpreted as a vector of coefficients of a polynomial of degree m - 1
    msg: &[F],
    rho_inv: usize,
) -> Vec<F> {
    let m = msg.len();

    let extended_domain = GeneralEvaluationDomain::<F>::new(m * rho_inv).unwrap_or_else(|| {
        panic!(
            "The field F cannot accomodate FFT for msg.len() * RHO_INV = {} elements (too many)",
            m * rho_inv
        )
    });

    extended_domain.fft(msg)
}

#[inline]
pub(crate) fn get_num_bytes(n: usize) -> usize {
    ceil_div((usize::BITS - n.leading_zeros()) as usize, 8)
}

/// Generate `t` (not necessarily distinct) random points in `[0, n)`
/// using the current state of the `transcript`.
pub(crate) fn get_indices_from_sponge<S: CryptographicSponge>(
    n: usize,
    t: usize,
    sponge: &mut S,
) -> Result<Vec<usize>, Error> {
    let bytes_to_squeeze = get_num_bytes(n);
    let mut indices = Vec::with_capacity(t);
    for _ in 0..t {
        let bytes = sponge.squeeze_bytes(bytes_to_squeeze);
        sponge.absorb(&bytes);

        // get the usize from Vec<u8>:
        let ind = bytes.iter().fold(0, |acc, &x| (acc << 8) + x as usize);
        // modulo the number of columns in the encoded matrix
        indices.push(ind % n);
    }
    Ok(indices)
}

#[inline]
pub(crate) fn calculate_t<F: PrimeField>(
    sec_param: usize,
    distance: (usize, usize),
    codeword_len: usize,
) -> Result<usize, Error> {
    // Took from the analysis by BCI+20 and Ligero
    // We will find the smallest $t$ such that
    // $(1-\delta)^t + (\rho+\delta)^t + \frac{n}{F} < 2^{-\lambda}$.
    // With $\delta = \frac{1-\rho}{2}$, the expreesion is
    // $2 * (\frac{1+\rho}{2})^t + \frac{n}{F} < 2^(-\lambda)$.

    let field_bits = F::MODULUS_BIT_SIZE as i32;
    let sec_param = sec_param as i32;

    let residual = codeword_len as f64 / 2.0_f64.powi(field_bits);
    let rhs = (2.0_f64.powi(-sec_param) - residual).log2();
    if !(rhs.is_normal()) {
        return Err(Error::InvalidParameters("For the given codeword length and the required security guarantee, the field is not big enough.".to_string()));
    }
    let nom = rhs - 1.0;
    let denom = (1.0 - 0.5 * distance.0 as f64 / distance.1 as f64).log2();
    if !(denom.is_normal()) {
        return Err(Error::InvalidParameters(
            "The distance is wrong".to_string(),
        ));
    }
    let t = (nom / denom).ceil() as usize;
    Ok(if t < codeword_len { t } else { codeword_len })
}

#[cfg(test)]
pub(crate) struct LeafIdentityHasher;

#[cfg(test)]
impl CRHScheme for LeafIdentityHasher {
    type Input = Vec<u8>;
    type Output = Vec<u8>;
    type Parameters = ();

    fn setup<R: RngCore>(_: &mut R) -> Result<Self::Parameters, ark_crypto_primitives::Error> {
        Ok(())
    }

    fn evaluate<T: Borrow<Self::Input>>(
        _: &Self::Parameters,
        input: T,
    ) -> Result<Self::Output, ark_crypto_primitives::Error> {
        Ok(input.borrow().to_vec().into())
    }
}

#[cfg(test)]
pub(crate) struct FieldToBytesColHasher<F, D>
where
    F: PrimeField + CanonicalSerialize,
    D: Digest,
{
    _phantom: PhantomData<(F, D)>,
}

#[cfg(test)]
impl<F, D> CRHScheme for FieldToBytesColHasher<F, D>
where
    F: PrimeField + CanonicalSerialize,
    D: Digest,
{
    type Input = Vec<F>;
    type Output = Vec<u8>;
    type Parameters = ();

    fn setup<R: RngCore>(_rng: &mut R) -> Result<Self::Parameters, ark_crypto_primitives::Error> {
        Ok(())
    }

    fn evaluate<T: Borrow<Self::Input>>(
        _parameters: &Self::Parameters,
        input: T,
    ) -> Result<Self::Output, ark_crypto_primitives::Error> {
        let mut dig = D::new();
        dig.update(to_bytes!(input.borrow()).unwrap());
        Ok(dig.finalize().to_vec())
    }
}

pub(crate) fn tensor_vec<F: PrimeField>(values: &[F]) -> Vec<F> {
    let one = F::one();
    let anti_values: Vec<F> = values.iter().map(|v| one - *v).collect();

    let mut layer: Vec<F> = vec![one];

    for i in 0..values.len() {
        let mut new_layer = Vec::new();
        for v in &layer {
            new_layer.push(*v * anti_values[i]);
        }
        for v in &layer {
            new_layer.push(*v * values[i]);
        }
        layer = new_layer;
    }

    layer
}

#[cfg(test)]
pub(crate) mod tests {
    use crate::linear_codes::utils::{calculate_t, get_num_bytes, reed_solomon, SprsMat};
    use crate::utils::to_field;
    use ark_bls12_377::Fq;
    use ark_bls12_377::Fr;
    use ark_ff::PrimeField;
    use ark_poly::{
        domain::general::GeneralEvaluationDomain, univariate::DensePolynomial, DenseUVPolynomial,
        EvaluationDomain, Polynomial,
    };
    use ark_std::test_rng;
    use rand_chacha::{rand_core::SeedableRng, ChaCha20Rng};

    #[test]
    fn test_sprs_row_mul() {
        // The columns major representation of a matrix.
        let mat: Vec<Fr> = to_field(vec![10, 23, 55, 100, 1, 58, 4, 0, 9]);

        let mat = SprsMat::new_from_flat(3, 3, 3, &mat);
        let v: Vec<Fr> = to_field(vec![12, 41, 55]);
        // by giving the result in the integers and then converting to Fr
        // we ensure the test will still pass even if Fr changes
        assert_eq!(mat.row_mul(&v), to_field::<Fr>(vec![4088, 4431, 543]));
    }

    #[test]
    fn test_sprs_row_mul_sparse_mat() {
        // The columns major representation of a matrix.
        let mat: Vec<Fr> = to_field(vec![10, 23, 55, 100, 1, 58, 4, 0, 9]);
        let mat = vec![
            vec![(0usize, mat[0]), (1usize, mat[1]), (2usize, mat[2])],
            vec![(0usize, mat[3]), (1usize, mat[4]), (2usize, mat[5])],
            vec![(0usize, mat[6]), (1usize, mat[7]), (2usize, mat[8])],
        ];

        let mat = SprsMat::new_from_columns(3, 3, 3, &mat);
        let v: Vec<Fr> = to_field(vec![12, 41, 55]);
        // by giving the result in the integers and then converting to Fr
        // we ensure the test will still pass even if Fr changes
        assert_eq!(mat.row_mul(&v), to_field::<Fr>(vec![4088, 4431, 543]));
    }

    #[test]
    fn test_reed_solomon() {
        let rho_inv = 3;
        // `i` is the min number of evaluations we need to interpolate a poly of degree `i - 1`
        for i in 1..10 {
            let deg = (1 << i) - 1;

            let rand_chacha = &mut ChaCha20Rng::from_rng(test_rng()).unwrap();
            let mut pol = DensePolynomial::rand(deg, rand_chacha);

            while pol.degree() != deg {
                pol = DensePolynomial::rand(deg, rand_chacha);
            }

            let coeffs = &pol.coeffs;

            // size of evals might be larger than deg + 1 (the min. number of evals needed to interpolate): we could still do R-S encoding on smaller evals, but the resulting polynomial will differ, so for this test to work we should pass it in full
            let m = deg + 1;

            let encoded = reed_solomon(&coeffs, rho_inv);

            let large_domain = GeneralEvaluationDomain::<Fr>::new(m * rho_inv).unwrap();

            // the encoded elements should agree with the evaluations of the polynomial in the larger domain
            for j in 0..(rho_inv * m) {
                assert_eq!(pol.evaluate(&large_domain.element(j)), encoded[j]);
            }
        }
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

    #[test]
    fn test_calculate_t_with_good_parameters() {
        assert!(calculate_t::<Fq>(128, (3, 4), 2_usize.pow(32)).unwrap() < 200);
        assert!(calculate_t::<Fq>(256, (3, 4), 2_usize.pow(32)).unwrap() < 400);
    }

    #[test]
    fn test_calculate_t_with_bad_parameters() {
        calculate_t::<Fq>(
            (Fq::MODULUS_BIT_SIZE - 60) as usize,
            (3, 4),
            2_usize.pow(60),
        )
        .unwrap_err();
        calculate_t::<Fq>(400, (3, 4), 2_usize.pow(32)).unwrap_err();
    }
}
