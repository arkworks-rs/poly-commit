use crate::*;
use crate::{PCCommitterKey, PCVerifierKey, Vec};
use algebra_core::{ProjectiveCurve, ToBytes, PrimeField};
use std::ops::AddAssign;

pub struct UniversalParams<G: ProjectiveCurve> {
    pub comm_key: Vec<G::Affine>,
    pub h: G
}

impl<G: ProjectiveCurve> PCUniversalParams for UniversalParams<G> {
    fn max_degree(&self) -> usize {
        self.elems.len() - 1
    }
}

#[derive(Derivative)]
#[derivative(
    Default(bound = ""),
    Hash(bound = ""),
    Clone(bound = ""),
    Debug(bound = "")
)]
pub struct Key<G: ProjectiveCurve> {
    pub comm_key: Vec<G::Affine>,
    pub h: G,
    pub max_degree: usize,
}

pub type VerifierKey<G> = Key<G>;
pub type CommitterKey<G> = Key<G>;

impl<G: ProjectiveCurve> PCCommitterKey for Key<G> {
    fn max_degree(&self) -> usize {
        self.max_degree
    }
    fn supported_degree(&self) -> usize {
        self.comm_key.len() - 1
    }
}

impl<G: ProjectiveCurve> PCVerifierKey for Key<G> {
    fn max_degree(&self) -> usize {
        self.max_degree
    }

    fn supported_degree(&self) -> usize {
        self.comm_key.len() - 1
    }
}

pub struct Commitment<G: ProjectiveCurve> (pub G);

impl<G: ProjectiveCurve> PCCommitment for Commitment<G> {
    #[inline]
    fn empty() -> Self {
        Commitment(G::zero());
    }

    fn has_degree_bound(&self) -> bool {
        false
    }

    fn size_in_bytes(&self) -> usize {
        algebra_core::to_bytes![G::zero()].unwrap().len() / 2
    }
}

impl<G: ProjectiveCurve> AddAssign<(G::ScalarField, &Commitment<G>)> for Commitment<G> {
    #[inline]
    fn add_assign(&mut self, (f, other): (G::ScalarField, &Commitment<G>)) {
        self.0 += other.0.mul(f.into_repr());
    }
}

/// `Randomness` hides the polynomial inside a commitment. It is output by `KZG10::commit`.
#[derive(Derivative)]
#[derivative(
Default(bound = ""),
Hash(bound = ""),
Clone(bound = ""),
Debug(bound = ""),
PartialEq(bound = ""),
Eq(bound = "")
)]
pub struct Randomness();
impl PCRandomness for Randomness {
    fn empty() -> Self {
        Randomness ()
    }

    fn rand<R: RngCore>(hiding_bound: usize, _: bool, rng: &mut R) -> Self {
        Randomness()
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
pub struct Proof<G: ProjectiveCurve> {
    pub L: Vec<G>,
    pub R: Vec<G>,
    pub comm_key: G,
    pub p: G::ScalarField
}

impl<G: ProjectiveCurve> PCProof for Proof<G> {
    fn size_in_bytes(&self) -> usize {
        algebra_core::to_bytes![G::zero()].unwrap().len() / 2
            + algebra_core::to_bytes![G::zero()].unwrap().len()
    }
}

impl<G: ProjectiveCurve> ToBytes for Proof<G> {
    #[inline]
    fn write<W: algebra_core::io::Write>(&self, mut writer: W) -> algebra_core::io::Result<()> {
        self.L.write(&mut writer);
        self.R.write(&mut writer);
        self.comm_key.write(&mut writer);
        self.p.write(&mut writer);
    }
}
