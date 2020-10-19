use ark_ff::{Field, ToBytes, Zero};
use ark_std::iter::Sum;
use ark_std::ops::{Add, Mul};
use core::fmt::Debug;

pub trait LHUniversalParameters: Clone + Debug {
    fn max_elems_len(&self) -> usize;
}

pub trait LHCommitterKey: Clone + Debug {
    fn max_elems_len(&self) -> usize;
    fn supported_elems_len(&self) -> usize;
}

pub trait LHCommitment<F: Field>:
    Clone
    + ToBytes
    + Zero
    + Eq
    + Add<Self, Output = Self>
    + for<'a> Add<&'a Self, Output = Self>
    + Mul<F, Output = Self>
    + for<'a> Sum<(F, &'a Self)>
{
    fn size_in_bytes(&self) -> usize;
}
