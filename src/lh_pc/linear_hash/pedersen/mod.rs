use crate::lh_pc::linear_hash::LinearHashFunction;
use crate::Error;
use ark_ec::msm::VariableBaseMSM;
use ark_ec::{AffineCurve, ProjectiveCurve};
use ark_ff::{to_bytes, PrimeField};
use ark_std::vec::Vec;
use core::marker::PhantomData;
use digest::Digest;
use rand_core::RngCore;

#[cfg(feature = "parallel")]
use rayon::prelude::*;

pub mod data_structures;
pub use data_structures::*;

pub struct PedersenCommitment<G: AffineCurve, D: Digest> {
    _field: PhantomData<G>,
    _digest: PhantomData<D>,
}

impl<G: AffineCurve, D: Digest> PedersenCommitment<G, D> {
    pub const PROTOCOL_NAME: &'static [u8] = b"pedersen-linear-hash";
}

impl<G: AffineCurve, D: Digest> LinearHashFunction<G::ScalarField> for PedersenCommitment<G, D> {
    type UniversalParams = UniversalParams<G>;
    type CommitterKey = CommitterKey<G>;
    type Commitment = Commitment<G>;
    type Error = Error;

    fn setup<R: RngCore>(
        max_num_elems: usize,
        _rng: &mut R,
    ) -> Result<Self::UniversalParams, Self::Error> {
        let generators: Vec<_> = ark_std::cfg_into_iter!(0..max_num_elems)
            .map(|i| {
                let i = i as u64;
                let mut hash = D::digest(&to_bytes![&Self::PROTOCOL_NAME, i].unwrap());
                let mut g = G::from_random_bytes(&hash);
                let mut j = 0u64;
                while g.is_none() {
                    hash = D::digest(&to_bytes![&Self::PROTOCOL_NAME, i, j].unwrap());
                    g = G::from_random_bytes(&hash);
                    j += 1;
                }
                let generator = g.unwrap();
                generator.mul_by_cofactor_to_projective()
            })
            .collect();

        let generators = G::Projective::batch_normalization_into_affine(&generators);
        let pp = UniversalParams(generators);
        Ok(pp)
    }

    fn trim(
        pp: &Self::UniversalParams,
        supported_num_elems: usize,
    ) -> Result<Self::CommitterKey, Self::Error> {
        let trimmed = pp.0[0..supported_num_elems].to_vec();
        let ck = CommitterKey {
            generators: trimmed,
            max_elems_len: pp.0.len(),
        };

        Ok(ck)
    }

    fn commit(
        ck: &Self::CommitterKey,
        elems: &[G::ScalarField],
    ) -> Result<Self::Commitment, Self::Error> {
        let scalars_bigint = ark_std::cfg_iter!(elems)
            .map(|s| s.into_repr())
            .collect::<Vec<_>>();

        let comm = VariableBaseMSM::multi_scalar_mul(&ck.generators, &scalars_bigint);
        let conversion = G::Projective::batch_normalization_into_affine(&[comm])
            .pop()
            .unwrap();
        Ok(Commitment(conversion))
    }
}
