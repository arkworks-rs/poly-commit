use ark_ec::msm::VariableBaseMSM;
use ark_ec::{AffineCurve, ProjectiveCurve};
use ark_ff::{PrimeField, ToBytes, Zero};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, SerializationError};
use ark_std::vec::Vec;
use ark_std::{
    io::{Read, Write},
    ops::{Add, Mul},
};

use ark_std::iter::Sum;

#[cfg(feature = "parallel")]
use rayon::prelude::*;

#[derive(Clone, Debug, PartialEq, Eq, CanonicalSerialize, CanonicalDeserialize)]
pub struct UniversalParams<G: AffineCurve> {
    pub(crate) generators: Vec<G>,
    pub(crate) hiding_generator: G,
}

#[derive(Clone, Debug, PartialEq, Eq, CanonicalSerialize, CanonicalDeserialize)]
pub struct CommitterKey<G: AffineCurve> {
    pub(crate) generators: Vec<G>,
    pub(crate) hiding_generator: G,
    pub(crate) max_elems: usize,
}

impl<G: AffineCurve> CommitterKey<G> {
    pub fn supported_elems_len(&self) -> usize {
        self.generators.len()
    }

    pub fn max_elems_len(&self) -> usize {
        self.max_elems
    }
}
