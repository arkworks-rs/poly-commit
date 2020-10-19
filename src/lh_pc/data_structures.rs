use crate::lh_pc::linear_hash::data_structures::{
    LHCommitment, LHCommitterKey, LHUniversalParameters,
};
use crate::lh_pc::LinearHashFunction;
use crate::{PCCommitment, PCCommitterKey, PCPreparedCommitment, PCPreparedVerifierKey, PCProof, PCRandomness, PCUniversalParams, PCVerifierKey, Polynomial, LabeledPolynomial};
use ark_ff::{to_bytes, Field, ToBytes, Zero};
use ark_poly::UVPolynomial;
use core::marker::PhantomData;
use rand_core::RngCore;

#[derive(Derivative)]
#[derivative(Clone(bound = ""), Debug(bound = ""))]
pub struct UniversalParameters<F: Field, LH: LinearHashFunction<F>>(pub(crate) LH::UniversalParams);

impl<F: Field, LH: LinearHashFunction<F>> PCUniversalParams for UniversalParameters<F, LH> {
    fn max_degree(&self) -> usize {
        self.0.max_elems_len() - 1
    }
}

#[derive(Derivative)]
#[derivative(Clone(bound = ""), Debug(bound = ""))]
pub struct CommitterKey<F: Field, LH: LinearHashFunction<F>>(pub(crate) LH::CommitterKey);

impl<F: Field, LH: LinearHashFunction<F>> PCCommitterKey for CommitterKey<F, LH> {
    fn max_degree(&self) -> usize {
        self.0.max_elems_len() - 1
    }

    fn supported_degree(&self) -> usize {
        self.0.supported_elems_len() - 1
    }
}

#[derive(Derivative)]
#[derivative(Clone(bound = ""), Debug(bound = ""))]
pub struct VerifierKey<F: Field, LH: LinearHashFunction<F>>(pub(crate) LH::CommitterKey);

impl<F: Field, LH: LinearHashFunction<F>> PCVerifierKey for VerifierKey<F, LH> {
    fn max_degree(&self) -> usize {
        self.0.max_elems_len() - 1
    }

    fn supported_degree(&self) -> usize {
        self.0.supported_elems_len() - 1
    }
}

impl<F: Field, LH: LinearHashFunction<F>> PCPreparedVerifierKey<VerifierKey<F, LH>>
    for VerifierKey<F, LH>
{
    fn prepare(vk: &VerifierKey<F, LH>) -> Self {
        vk.clone()
    }
}

#[derive(Derivative)]
#[derivative(Clone(bound = ""), Eq(bound = ""), PartialEq(bound = ""))]
pub struct Commitment<F: Field, LH: LinearHashFunction<F>>(pub LH::Commitment);

impl<F: Field, LH: LinearHashFunction<F>> PCCommitment for Commitment<F, LH> {
    fn empty() -> Self {
        Self(LH::Commitment::zero())
    }

    fn has_degree_bound(&self) -> bool {
        false
    }

    fn size_in_bytes(&self) -> usize {
        LHCommitment::<F>::size_in_bytes(&self.0)
    }
}

impl<F: Field, LH: LinearHashFunction<F>> PCPreparedCommitment<Commitment<F, LH>>
    for Commitment<F, LH>
{
    fn prepare(comm: &Commitment<F, LH>) -> Self {
        comm.clone()
    }
}

impl<F: Field, LH: LinearHashFunction<F>> Default for Commitment<F, LH> {
    fn default() -> Self {
        Self(LH::Commitment::zero())
    }
}

impl<F: Field, LH: LinearHashFunction<F>> ToBytes for Commitment<F, LH> {
    fn write<W: ark_std::io::Write>(&self, writer: W) -> ark_std::io::Result<()> {
        self.0.write(writer)
    }
}

#[derive(Derivative)]
#[derivative(Clone(bound = ""))]
pub struct Randomness(pub(crate) ());

impl PCRandomness for Randomness {
    fn empty() -> Self {
        Self(())
    }

    fn rand<R: RngCore>(
        _num_queries: usize,
        _has_degree_bound: bool,
        _: Option<usize>,
        _rng: &mut R,
    ) -> Self {
        Self(())
    }
}

#[derive(Derivative)]
#[derivative(Clone(bound = ""))]
pub struct Proof<F: Field, P: UVPolynomial<F>>(pub LabeledPolynomial<F, P>);

impl<F: Field, P: UVPolynomial<F>> ToBytes for Proof<F, P> {
    fn write<W: ark_std::io::Write>(&self, mut writer: W) -> ark_std::io::Result<()> {
        writer.write_all((to_bytes!(self.0.coeffs())?).as_slice())
    }
}

impl<F: Field, P: UVPolynomial<F>> PCProof for Proof<F, P> {
    fn size_in_bytes(&self) -> usize {
        to_bytes![self].unwrap().len()
    }
}
