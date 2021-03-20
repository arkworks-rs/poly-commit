use crate::trivial_pc::{CommitterKey, UniversalParams};
use ark_ec::msm::VariableBaseMSM;
use ark_ec::{AffineCurve, ProjectiveCurve};
use ark_ff::{to_bytes, PrimeField};
use ark_std::vec::Vec;
use blake2::Blake2s;
use core::marker::PhantomData;
use digest::Digest;

#[cfg(feature = "parallel")]
use rayon::prelude::*;

/// The protocol name for the commitment scheme.
pub(crate) const PROTOCOL_NAME: &[u8] = b"Pedersen-Commitment";

/// A simple polynomial commitment scheme from [[BCLMS20]][bclms] built on Pedersen commitments.
///
/// [bclms]: https://eprint.iacr.org/2020/1618
pub struct PedersenCommitment<G: AffineCurve> {
    _group: PhantomData<G>,
}

impl<G: AffineCurve> PedersenCommitment<G> {
    /// Outputs the universal parameters of a Pedersen commitment.
    pub fn setup(max_num_elems: usize) -> UniversalParams<G> {
        let generators: Vec<_> = ark_std::cfg_into_iter!(0..(max_num_elems + 1))
            .map(|i| {
                let i = i as u64;
                let mut hash = Blake2s::digest(&to_bytes![&PROTOCOL_NAME, i].unwrap());
                let mut g = G::from_random_bytes(&hash);
                let mut j = 0u64;
                while g.is_none() {
                    hash = Blake2s::digest(&to_bytes![&PROTOCOL_NAME, i, j].unwrap());
                    g = G::from_random_bytes(&hash);
                    j += 1;
                }
                let generator = g.unwrap();
                generator.mul_by_cofactor_to_projective()
            })
            .collect();

        let mut generators = G::Projective::batch_normalization_into_affine(&generators);
        let hiding_generator = generators.pop().unwrap();

        let pp = UniversalParams {
            generators,
            hiding_generator,
        };

        pp
    }

    /// Outputs the committer key from the universal parameters.
    /// Assumes that the input `supported_num_elems` is within the universal parameter's
    /// number of supported elements.
    pub fn trim(pp: &UniversalParams<G>, supported_num_elems: usize) -> CommitterKey<G> {
        assert!(supported_num_elems <= pp.generators.len());

        let ck = CommitterKey {
            generators: pp.generators[0..supported_num_elems].to_vec(),
            hiding_generator: pp.hiding_generator,
            max_elems: pp.generators.len(),
        };

        ck
    }

    /// Computes a Pedersen commitment.
    /// Assumes that `ck` has enough elements to commit to `elems`.
    pub fn commit(
        ck: &CommitterKey<G>,
        elems: &[G::ScalarField],
        randomizer: Option<G::ScalarField>,
    ) -> G {
        assert!(elems.len() <= ck.generators.len());

        let scalars_bigint = ark_std::cfg_iter!(elems)
            .map(|s| s.into_repr())
            .collect::<Vec<_>>();

        let mut comm = VariableBaseMSM::multi_scalar_mul(&ck.generators, &scalars_bigint);
        if let Some(randomizer) = randomizer {
            comm += &ck.hiding_generator.mul(randomizer);
        }

        comm.into()
    }
}
