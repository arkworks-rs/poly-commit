use super::LigeroPCParams;
use super::LinCodeParametersInfo;
use crate::linear_codes::utils::calculate_t;
use crate::utils::ceil_div;
use crate::{PCCommitterKey, PCUniversalParams, PCVerifierKey};

use ark_crypto_primitives::crh::{CRHScheme, TwoToOneCRHScheme};
use ark_crypto_primitives::merkle_tree::{Config, LeafParam, TwoToOneParam};
use ark_ff::PrimeField;
use ark_std::log2;
use ark_std::marker::PhantomData;
#[cfg(not(feature = "std"))]
use num_traits::Float;

impl<F, C, H> LigeroPCParams<F, C, H>
where
    F: PrimeField,
    C: Config,
    H: CRHScheme,
{
    /// Create new UniversalParams
    pub fn new(
        sec_param: usize,
        rho_inv: usize,
        check_well_formedness: bool,
        leaf_hash_param: LeafParam<C>,
        two_to_one_hash_param: TwoToOneParam<C>,
        col_hash_params: H::Parameters,
    ) -> Self {
        Self {
            _field: PhantomData,
            sec_param,
            rho_inv,
            check_well_formedness,
            leaf_hash_param,
            two_to_one_hash_param,
            col_hash_params,
        }
    }
}

impl<F, C, H> PCUniversalParams for LigeroPCParams<F, C, H>
where
    F: PrimeField,
    C: Config,
    H: CRHScheme,
{
    fn max_degree(&self) -> usize {
        if F::TWO_ADICITY < self.rho_inv as u32 {
            0
        } else if (F::TWO_ADICITY - self.rho_inv as u32) * 2 < 64 {
            2_usize.pow((F::TWO_ADICITY - self.rho_inv as u32) * 2)
        } else {
            usize::MAX
        }
    }
}

impl<F, C, H> PCCommitterKey for LigeroPCParams<F, C, H>
where
    F: PrimeField,
    C: Config,
    H: CRHScheme,
{
    fn max_degree(&self) -> usize {
        if (F::TWO_ADICITY - self.rho_inv as u32) * 2 < 64 {
            2_usize.pow((F::TWO_ADICITY - self.rho_inv as u32) * 2)
        } else {
            usize::MAX
        }
    }

    fn supported_degree(&self) -> usize {
        <LigeroPCParams<F, C, H> as PCCommitterKey>::max_degree(self)
    }
}

impl<F, C, H> PCVerifierKey for LigeroPCParams<F, C, H>
where
    F: PrimeField,
    C: Config,
    H: CRHScheme,
{
    fn max_degree(&self) -> usize {
        if (F::TWO_ADICITY - self.rho_inv as u32) * 2 < 64 {
            2_usize.pow((F::TWO_ADICITY - self.rho_inv as u32) * 2)
        } else {
            usize::MAX
        }
    }

    fn supported_degree(&self) -> usize {
        <LigeroPCParams<F, C, H> as PCVerifierKey>::max_degree(self)
    }
}

impl<F, C, H> LinCodeParametersInfo<C, H> for LigeroPCParams<F, C, H>
where
    F: PrimeField,
    C: Config,
    H: CRHScheme,
{
    fn check_well_formedness(&self) -> bool {
        self.check_well_formedness
    }

    fn distance(&self) -> (usize, usize) {
        (self.rho_inv - 1, self.rho_inv)
    }

    fn sec_param(&self) -> usize {
        self.sec_param
    }

    /// Compute the a suitable (for instance, FFT-friendly over F) matrix with at least poly_len entries.
    /// The return pair (n, m) corresponds to the dimensions n x m.
    /// FIXME: Maybe, there should be some checks for making sure the extended row can have an FFT.
    fn compute_dimensions(&self, poly_len: usize) -> (usize, usize) {
        assert_eq!(
            (poly_len as f64) as usize,
            poly_len,
            "n cannot be converted to f64: aborting"
        );
        let t = calculate_t::<F>(self.sec_param(), self.distance(), poly_len).unwrap();
        let n = 1 << log2((ceil_div(2 * poly_len, t) as f64).sqrt().ceil() as usize);
        let m = ceil_div(poly_len, n);
        (n, m)
    }

    fn leaf_hash_param(&self) -> &<<C as Config>::LeafHash as CRHScheme>::Parameters {
        &self.leaf_hash_param
    }

    fn two_to_one_hash_param(
        &self,
    ) -> &<<C as Config>::TwoToOneHash as TwoToOneCRHScheme>::Parameters {
        &self.two_to_one_hash_param
    }

    fn col_hash_params(&self) -> &<H as CRHScheme>::Parameters {
        &self.col_hash_params
    }
}
