use super::{
    utils::{reed_solomon, tensor_vec},
    LigeroPCParams, LinearEncode,
};

use ark_crypto_primitives::{
    crh::{CRHScheme, TwoToOneCRHScheme},
    merkle_tree::Config,
    sponge::CryptographicSponge,
};
use ark_ff::{FftField, PrimeField};
use ark_poly::{MultilinearExtension, Polynomial};
use ark_std::log2;
use ark_std::marker::PhantomData;
use ark_std::vec::Vec;

mod tests;

/// The multilinear Ligero polynomial commitment scheme based on [[Ligero]][ligero].
/// The scheme defaults to the naive batching strategy.
///
/// Note: The scheme currently does not support hiding.
///
/// [ligero]: https://eprint.iacr.org/2022/1608.pdf
pub struct MultilinearLigero<
    F: PrimeField,
    C: Config,
    S: CryptographicSponge,
    P: MultilinearExtension<F>,
    H: CRHScheme,
> {
    _phantom: PhantomData<(F, C, S, P, H)>,
}

impl<F, C, S, P, H> LinearEncode<F, C, P, H> for MultilinearLigero<F, C, S, P, H>
where
    F: PrimeField + FftField,
    C: Config,
    S: CryptographicSponge,
    P: MultilinearExtension<F>,
    <P as Polynomial<F>>::Point: Into<Vec<F>>,
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
            2,
            true,
            leaf_hash_param,
            two_to_one_hash_param,
            col_hash_params,
        )
    }

    fn encode(msg: &[F], param: &Self::LinCodePCParams) -> Vec<F> {
        reed_solomon(msg, param.rho_inv)
    }

    fn poly_to_vec(polynomial: &P) -> Vec<F> {
        polynomial.to_evaluations()
    }

    fn point_to_vec(point: <P as Polynomial<F>>::Point) -> Vec<F> {
        point
    }

    /// For a multilinear polynomial in n+m variables it returns a tuple for k={n,m}:
    /// ((1-z_1)*(1-z_2)*...*(1_z_k), z_1*(1-z_2)*...*(1-z_k), ..., z_1*z_2*...*z_k)
    fn tensor(
        point: &<P as Polynomial<F>>::Point,
        left_len: usize,
        _right_len: usize,
    ) -> (Vec<F>, Vec<F>) {
        let point: Vec<F> = Self::point_to_vec(point.clone());

        let split = log2(left_len) as usize;
        let left = &point[..split];
        let right = &point[split..];
        (tensor_vec(left), tensor_vec(right))
    }
}
