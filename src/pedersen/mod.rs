use crate::Error;
use ark_ec::msm::VariableBaseMSM;
use ark_ec::{AffineCurve, ProjectiveCurve};
use ark_ff::{to_bytes, PrimeField};
use ark_std::vec::Vec;
use core::marker::PhantomData;
use digest::Digest;
use rand_core::RngCore;
use blake2::Blake2s;

#[cfg(feature = "parallel")]
use rayon::prelude::*;

pub mod data_structures;
pub use data_structures::*;

pub struct PedersenCommitment<G: AffineCurve> {
    _group: PhantomData<G>,
}

impl<G: AffineCurve> PedersenCommitment<G> {
    pub const PROTOCOL_NAME: &'static [u8] = b"pedersen-linear-hash";
}

impl<G: AffineCurve> PedersenCommitment<G> {
    pub fn setup<R: RngCore>(
        max_num_elems: usize,
        _rng: &mut R,
    ) -> Result<UniversalParams<G>, Error> {
        let generators: Vec<_> = ark_std::cfg_into_iter!(0..max_num_elems)
            .map(|i| {
                let i = i as u64;
                let mut hash = Blake2s::digest(&to_bytes![&Self::PROTOCOL_NAME, i].unwrap());
                let mut g = G::from_random_bytes(&hash);
                let mut j = 0u64;
                while g.is_none() {
                    hash = Blake2s::digest(&to_bytes![&Self::PROTOCOL_NAME, i, j].unwrap());
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

    pub fn trim(
        pp: &UniversalParams<G>,
        supported_num_elems: usize,
    ) -> Result<CommitterKey<G>, Error> {
        let trimmed = pp.0[0..supported_num_elems].to_vec();
        let ck = CommitterKey {
            generators: trimmed,
            max_elems_len: pp.0.len(),
        };

        Ok(ck)
    }

    pub fn commit(
        ck: &CommitterKey<G>,
        elems: &[G::ScalarField],
    ) -> Result<Commitment<G>, Error> {
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