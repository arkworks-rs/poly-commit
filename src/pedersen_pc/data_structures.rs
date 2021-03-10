use crate::pedersen_pc::pedersen::{CommitterKey, UniversalParams};
use crate::{
    LabeledPolynomial, PCCommitment, PCCommitterKey, PCPreparedCommitment, PCPreparedVerifierKey,
    PCProof, PCRandomness, PCUniversalParams, PCVerifierKey,
};
use ark_ec::msm::VariableBaseMSM;
use ark_ec::{AffineCurve, ProjectiveCurve};
use ark_ff::PrimeField;
use ark_ff::{to_bytes, Field, ToBytes, Zero};
use ark_poly::UVPolynomial;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, SerializationError};
use ark_std::io::{Read, Write};
use ark_std::iter::Sum;
use rand_core::RngCore;

#[cfg(feature = "parallel")]
use rayon::prelude::*;

impl<G: AffineCurve> PCUniversalParams for UniversalParams<G> {
    fn max_degree(&self) -> usize {
        self.generators.len() - 1
    }
}

impl<G: AffineCurve> PCCommitterKey for CommitterKey<G> {
    fn max_degree(&self) -> usize {
        self.max_elems - 1
    }

    fn supported_degree(&self) -> usize {
        self.generators.len() - 1
    }
}

impl<G: AffineCurve> PCVerifierKey for CommitterKey<G> {
    fn max_degree(&self) -> usize {
        self.max_elems - 1
    }

    fn supported_degree(&self) -> usize {
        self.generators.len() - 1
    }
}

impl<G: AffineCurve> PCPreparedVerifierKey<CommitterKey<G>> for CommitterKey<G> {
    fn prepare(vk: &CommitterKey<G>) -> Self {
        vk.clone()
    }
}

/// Commitment to a polynomial in PC_LH. This is equal to a Pedersen commitment.
#[derive(Clone, Eq, PartialEq, Debug, CanonicalSerialize, CanonicalDeserialize)]
pub struct Commitment<G: AffineCurve> {
    /// The Pedersen commitment output.
    pub elem: G,
}

impl<G: AffineCurve> PCCommitment for Commitment<G> {
    fn empty() -> Self {
        Self { elem: G::zero() }
    }

    fn has_degree_bound(&self) -> bool {
        false
    }

    fn size_in_bytes(&self) -> usize {
        unimplemented!()
    }
}

impl<G: AffineCurve> PCPreparedCommitment<Commitment<G>> for Commitment<G> {
    fn prepare(comm: &Commitment<G>) -> Self {
        comm.clone()
    }
}

impl<G: AffineCurve> Default for Commitment<G> {
    fn default() -> Self {
        Self::empty()
    }
}

impl<G: AffineCurve> ToBytes for Commitment<G> {
    fn write<W: Write>(&self, writer: W) -> ark_std::io::Result<()> {
        self.write(writer)
    }
}

impl<'a, G: AffineCurve> Sum<(G::ScalarField, Self)> for Commitment<G> {
    fn sum<I: Iterator<Item = (G::ScalarField, Self)>>(iter: I) -> Self {
        let mut scalars = Vec::new();
        let mut comms = Vec::new();

        for (f, comm) in iter {
            scalars.push(f);
            comms.push(comm.elem);
        }

        let scalars_bigint = ark_std::cfg_iter!(scalars.as_slice())
            .map(|s| s.into_repr())
            .collect::<Vec<_>>();

        let sum = VariableBaseMSM::multi_scalar_mul(comms.as_slice(), scalars_bigint.as_slice());
        let mut conversion = G::Projective::batch_normalization_into_affine(&[sum]);

        Self {
            elem: conversion.pop().unwrap(),
        }
    }
}

/// Opening randomness for PC_LH.
#[derive(Clone, CanonicalSerialize, CanonicalDeserialize)]
pub struct Randomness;

impl PCRandomness for Randomness {
    fn empty() -> Self {
        Self
    }

    fn rand<R: RngCore>(
        _num_queries: usize,
        _has_degree_bound: bool,
        _: Option<usize>,
        _rng: &mut R,
    ) -> Self {
        Self
    }
}

/// An evaluation proof in PC_LH. This is just equal to the polynomial itself.
#[derive(Clone, CanonicalSerialize, CanonicalDeserialize)]
pub struct Proof<F: Field, P: UVPolynomial<F>> {
    /// The opened polynomial.
    pub polynomial: LabeledPolynomial<F, P>,
}

impl<F: Field, P: UVPolynomial<F>> ToBytes for Proof<F, P> {
    fn write<W: Write>(&self, mut writer: W) -> ark_std::io::Result<()> {
        writer.write_all((to_bytes!(self.polynomial.coeffs())?).as_slice())
    }
}

impl<F: Field, P: UVPolynomial<F>> PCProof for Proof<F, P> {
    fn size_in_bytes(&self) -> usize {
        to_bytes![self].unwrap().len()
    }
}
