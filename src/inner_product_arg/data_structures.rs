use crate::*;
use crate::{PCCommitterKey, PCVerifierKey, Vec};
use algebra_core::{ToBytes, PrimeField, AffineCurve, ProjectiveCurve};
use std::ops::{AddAssign, Add, Sub};

#[derive(Derivative)]
#[derivative(Default(bound = ""), Clone(bound = ""), Debug(bound = ""))]
pub struct UniversalParams<G: AffineCurve> {
    pub comm_key: Vec<G>,
    pub h: G
}

impl<G: AffineCurve> PCUniversalParams for UniversalParams<G> {
    fn max_degree(&self) -> usize {
        self.comm_key.len() - 1
    }
}

#[derive(Derivative)]
#[derivative(
    Default(bound = ""),
    Hash(bound = ""),
    Clone(bound = ""),
    Debug(bound = "")
)]
pub struct CommitterKey<G: AffineCurve> {
    pub comm_key: Vec<G>,
    pub h: G,
    pub max_degree: usize,
}

impl<G: AffineCurve> PCCommitterKey for CommitterKey<G> {
    fn max_degree(&self) -> usize {
        self.max_degree
    }
    fn supported_degree(&self) -> usize {
        self.comm_key.len() - 1
    }
}

#[derive(Derivative)]
#[derivative(
Default(bound = ""),
Hash(bound = ""),
Clone(bound = ""),
Debug(bound = "")
)]
pub struct VerifierKey<G: AffineCurve> {
    pub comm_key: Vec<G>,
    pub h: G,
    pub max_degree: usize,
}

impl<G: AffineCurve> PCVerifierKey for VerifierKey<G> {
    fn max_degree(&self) -> usize {
        self.max_degree
    }

    fn supported_degree(&self) -> usize {
        self.comm_key.len() - 1
    }
}

#[derive(Derivative)]
#[derivative(
Default(bound = ""),
Hash(bound = ""),
Clone(bound = ""),
Copy(bound = ""),
Debug(bound = ""),
PartialEq(bound = ""),
Eq(bound = "")
)]
pub struct Commitment<G: AffineCurve> (pub G);

impl<G: AffineCurve> PCCommitment for Commitment<G> {
    #[inline]
    fn empty() -> Self {
        Commitment(G::zero())
    }

    fn has_degree_bound(&self) -> bool {
        false
    }

    fn size_in_bytes(&self) -> usize {
        algebra_core::to_bytes![G::zero()].unwrap().len() / 2
    }
}

impl<G: AffineCurve> ToBytes for Commitment<G> {
    #[inline]
    fn write<W: algebra_core::io::Write>(&self, writer: W) -> algebra_core::io::Result<()> {
        self.0.write(writer)
    }
}

/*

impl<G: AffineCurve> Add<G> for Commitment<G> {
    type Output = Self;

    #[inline]
    fn add(mut self, other: G) -> Self {
        self.0 += other;
        self
    }
}

impl<G: AffineCurve> Sub<G> for Commitment<G> {
    type Output = Self;

    #[inline]
    fn sub(mut self, other: G) -> Self {
        self.0 -= other;
        self
    }
}

impl<G: AffineCurve> AddAssign<G> for Commitment<G> {
    #[inline]
    fn add_assign(&mut self, g: G) {
        self.0 = self.0.into_projective().add_mixed(&g).into();
    }
}

impl<G: AffineCurve> AddAssign<(G::ScalarField, &Commitment<G>)> for Commitment<G> {
    #[inline]
    fn add_assign(&mut self, (f, other): (G::ScalarField, &Commitment<G>)) {
        let other_mul = other.0.mul(f.into_repr());
        self.0 = other_mul.add_mixed(&self.0).into();
    }
}
*/

#[derive(Derivative)]
#[derivative(
Default(bound = ""),
Hash(bound = ""),
Clone(bound = ""),
Debug(bound = ""),
PartialEq(bound = ""),
Eq(bound = "")
)]
// TODO: Placeholder randomness
pub struct Randomness();
impl PCRandomness for Randomness {
    fn empty() -> Self {
        Randomness ()
    }
    fn rand <R: RngCore> (_: usize, _: bool, _: &mut R) -> Self {
        Randomness()
    }
}

#[derive(Derivative)]
#[derivative(
Default(bound = ""),
Hash(bound = ""),
Clone(bound = ""),
Debug(bound = "")
)]
pub struct Proof<G: AffineCurve> {
    pub l_vec: Vec<G>,
    pub r_vec: Vec<G>,
    pub final_comm_key: G,
    pub c: G::ScalarField
}

impl<G: AffineCurve> PCProof for Proof<G> {
    fn size_in_bytes(&self) -> usize {
        algebra_core::to_bytes![self].unwrap().len()
    }
}

impl<G: AffineCurve> ToBytes for Proof<G> {
    #[inline]
    fn write<W: algebra_core::io::Write>(&self, mut writer: W) -> algebra_core::io::Result<()> {
        self.l_vec.write(&mut writer)?;
        self.r_vec.write(&mut writer)?;
        self.final_comm_key.write(&mut writer)?;
        self.c.write(&mut writer)
    }
}
