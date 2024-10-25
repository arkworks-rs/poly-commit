use crate::{
    linear_codes::{utils::reed_solomon, LigeroPCParams, LinearEncode},
    Error,
};
use ark_crypto_primitives::{
    crh::{CRHScheme, TwoToOneCRHScheme},
    merkle_tree::Config,
};
use ark_ff::PrimeField;
use ark_poly::DenseUVPolynomial;
use ark_std::marker::PhantomData;
#[cfg(not(feature = "std"))]
use ark_std::vec::Vec;

mod tests;

/// The univariate Ligero polynomial commitment scheme based on [[Ligero]][ligero].
/// The scheme defaults to the naive batching strategy.
///
/// Note: The scheme currently does not support hiding.
///
/// [ligero]: https://eprint.iacr.org/2022/1608.pdf
pub struct UnivariateLigero<F: PrimeField, C: Config, P: DenseUVPolynomial<F>, H: CRHScheme> {
    _phantom: PhantomData<(F, C, P, H)>,
}

impl<F, C, P, H> LinearEncode<F, C, P, H> for UnivariateLigero<F, C, P, H>
where
    F: PrimeField,
    C: Config,
    P: DenseUVPolynomial<F>,
    P::Point: Into<F>,
    H: CRHScheme,
{
    type LinCodePCParams = LigeroPCParams<F, C, H>;

    fn setup<R>(
        _max_degree: usize,
        _num_vars: Option<usize>,
        _rng: &mut R,
        leaf_hash_param: <<C as Config>::LeafHash as CRHScheme>::Parameters,
        two_to_one_hash_param: <<C as Config>::TwoToOneHash as TwoToOneCRHScheme>::Parameters,
        col_hash_params: H::Parameters,
    ) -> Self::LinCodePCParams {
        Self::LinCodePCParams::new(
            128,
            4,
            true,
            leaf_hash_param,
            two_to_one_hash_param,
            col_hash_params,
        )
    }

    fn encode(msg: &[F], param: &Self::LinCodePCParams) -> Result<Vec<F>, Error> {
        Ok(reed_solomon(msg, param.rho_inv))
    }

    /// For a univariate polynomial, we simply return the list of coefficients.
    fn poly_to_vec(polynomial: &P) -> Vec<F> {
        polynomial.coeffs().to_vec()
    }

    fn point_to_vec(point: P::Point) -> Vec<F> {
        vec![point]
    }

    /// For a univariate polynomial it returns a tuple:
    /// ((1, z, z^2, ..., z^n), (1, z^n, z^(2n), ..., z^((m-1)n)))
    fn tensor(z: &F, left: usize, right: usize) -> (Vec<F>, Vec<F>) {
        let mut left_out = Vec::with_capacity(left);
        let mut pow_a = F::one();
        for _ in 0..left {
            left_out.push(pow_a);
            pow_a *= z;
        }

        let mut right_out = Vec::with_capacity(right);
        let mut pow_b = F::one();
        for _ in 0..right {
            right_out.push(pow_b);
            pow_b *= pow_a;
        }

        (left_out, right_out)
    }
}
