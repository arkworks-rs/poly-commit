use crate::{PCCommitment, PCCommitterKey, PCPreparedCommitment, PCPreparedVerifierKey, PCProof, PCRandomness, PCUniversalParams, PCVerifierKey, LabeledPolynomial};
use ark_ff::{to_bytes, Field, ToBytes, Zero};
use ark_ec::AffineCurve;
use ark_poly::UVPolynomial;
use ark_serialize::{CanonicalSerialize, CanonicalDeserialize, SerializationError};
use ark_std::io::{Read, Write};
use rand_core::RngCore;
use crate::pedersen;

#[derive(Clone, Debug, CanonicalSerialize, CanonicalDeserialize)]
pub struct UniversalParameters<G: AffineCurve>(pub(crate) pedersen::UniversalParams<G>);

impl<G: AffineCurve> PCUniversalParams for UniversalParameters<G> {
    fn max_degree(&self) -> usize {
        (self.0).0.len() - 1
    }
}

#[derive(Clone, Debug, CanonicalSerialize, CanonicalDeserialize)]
pub struct CommitterKey<G: AffineCurve>(pub(crate) pedersen::CommitterKey<G>);

impl<G: AffineCurve> PCCommitterKey for CommitterKey<G> {
    fn max_degree(&self) -> usize {
        self.0.max_elems_len - 1
    }

    fn supported_degree(&self) -> usize {
        self.0.generators.len() - 1
    }
}

pub type VerifierKey<G> = CommitterKey<G>;

impl<G: AffineCurve> PCVerifierKey for VerifierKey<G> {
    fn max_degree(&self) -> usize {
        self.0.max_elems_len - 1
    }

    fn supported_degree(&self) -> usize {
        self.0.generators.len() - 1
    }
}

impl<G: AffineCurve> PCPreparedVerifierKey<VerifierKey<G>>
    for VerifierKey<G>
{
    fn prepare(vk: &VerifierKey<G>) -> Self {
        vk.clone()
    }
}

#[derive(Clone, Eq, PartialEq, Debug, CanonicalSerialize, CanonicalDeserialize)]
pub struct Commitment<G: AffineCurve>(pub pedersen::Commitment<G>);

impl<G: AffineCurve> PCCommitment for Commitment<G> {
    fn empty() -> Self {
        Self(pedersen::Commitment::zero())
    }

    fn has_degree_bound(&self) -> bool {
        false
    }

    fn size_in_bytes(&self) -> usize {
        unimplemented!()
    }
}

impl<G: AffineCurve> PCPreparedCommitment<Commitment<G>>
    for Commitment<G>
{
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
        self.0.write(writer)
    }
}

#[derive(Clone)]
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

#[derive(Clone, CanonicalSerialize, CanonicalDeserialize)]
pub struct Proof<F: Field, P: UVPolynomial<F>>(pub LabeledPolynomial<F, P>);

impl<F: Field, P: UVPolynomial<F>> ToBytes for Proof<F, P> {
    fn write<W: Write>(&self, mut writer: W) -> ark_std::io::Result<()> {
        writer.write_all((to_bytes!(self.0.coeffs())?).as_slice())
    }
}

impl<F: Field, P: UVPolynomial<F>> PCProof for Proof<F, P> {
    fn size_in_bytes(&self) -> usize {
        to_bytes![self].unwrap().len()
    }
}