use crate::lh_pc::linear_hash::data_structures::{
    LHCommitment, LHCommitterKey, LHUniversalParameters,
};
use ark_ff::Field;
use rand_core::RngCore;

pub mod data_structures;
pub mod pedersen;

pub trait LinearHashFunction<F: Field> {
    type UniversalParams: LHUniversalParameters;
    type CommitterKey: LHCommitterKey;
    type Commitment: LHCommitment<F>;

    type Error: ark_std::error::Error;

    fn setup<R: RngCore>(
        max_num_elems: usize,
        rng: &mut R,
    ) -> Result<Self::UniversalParams, Self::Error>;

    fn trim(
        pp: &Self::UniversalParams,
        supported_num_elems: usize,
    ) -> Result<Self::CommitterKey, Self::Error>;

    fn commit(ck: &Self::CommitterKey, elems: &[F]) -> Result<Self::Commitment, Self::Error>;
}
