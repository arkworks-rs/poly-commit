use ark_ec::msm::VariableBaseMSM;
use ark_ec::{AffineCurve, ProjectiveCurve};
use ark_ff::{PrimeField, ToBytes, Zero};
use ark_std::{io::{Read, Write}, ops::{Add, Mul}};
use ark_std::vec::Vec;
use ark_serialize::{CanonicalSerialize, CanonicalDeserialize, SerializationError};

use ark_std::iter::Sum;
#[cfg(feature = "parallel")]
use rayon::prelude::*;

#[derive(Clone, Debug, PartialEq, Eq, CanonicalSerialize, CanonicalDeserialize)]
pub struct UniversalParams<G: AffineCurve>(pub(crate) Vec<G>);

#[derive(Clone, Debug, PartialEq, Eq, CanonicalSerialize, CanonicalDeserialize)]
pub struct CommitterKey<G: AffineCurve> {
    pub(crate) generators: Vec<G>,
    pub(crate) max_elems_len: usize,
}

#[derive(Clone, Debug, PartialEq, Eq, CanonicalSerialize, CanonicalDeserialize)]
pub struct Commitment<G: AffineCurve>(pub G);

impl<G: AffineCurve> Zero for Commitment<G> {
    fn zero() -> Self {
        Self(G::zero())
    }

    fn is_zero(&self) -> bool {
        self.0.is_zero()
    }
}

impl<G: AffineCurve> ToBytes for Commitment<G> {
    fn write<W: Write>(&self, writer: W) -> ark_std::io::Result<()> {
        self.0.write(writer)
    }
}

impl<G: AffineCurve> Add<Self> for Commitment<G> {
    type Output = Self;
    fn add(self, rhs: Commitment<G>) -> Self::Output {
        Self(self.0 + rhs.0)
    }
}

impl<G: AffineCurve> Add<&Self> for Commitment<G> {
    type Output = Self;
    fn add(self, rhs: &Commitment<G>) -> Self::Output {
        Self(self.0 + rhs.0.clone())
    }
}

impl<G: AffineCurve> Mul<G::ScalarField> for Commitment<G> {
    type Output = Self;
    fn mul(self, rhs: <G as AffineCurve>::ScalarField) -> Self::Output {
        let product = self.0.mul(rhs);
        let mut conversion = G::Projective::batch_normalization_into_affine(&[product]);
        Self(conversion.pop().unwrap())
    }
}

impl<'a, G: AffineCurve> Sum<(G::ScalarField, &'a Self)> for Commitment<G> {
    fn sum<I: Iterator<Item = (G::ScalarField, &'a Self)>>(iter: I) -> Self {
        let mut scalars = Vec::new();
        let mut comms = Vec::new();

        for (f, comm) in iter {
            scalars.push(f);
            comms.push(comm.0.clone());
        }

        let scalars_bigint = ark_std::cfg_iter!(scalars.as_slice())
            .map(|s| s.into_repr())
            .collect::<Vec<_>>();

        let sum = VariableBaseMSM::multi_scalar_mul(comms.as_slice(), scalars_bigint.as_slice());
        let mut conversion = G::Projective::batch_normalization_into_affine(&[sum]);
        Self(conversion.pop().unwrap())
    }
}