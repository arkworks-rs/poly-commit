use rand_core::RngCore;
use algebra::{ToBytes, PairingEngine};
use crate::{PCCommitment, PCCommitterKey, PCRandomness, PCVerifierKey};
use std::ops::{Add, AddAssign};

use crate::kzg10;

pub struct UniversalParams<E: PairingEngine> {
    // Elements in the form {\beta^i G1} for 0 <= i <= max_degree
    pub powers_of_g: Vec<E::G1Affine>,
    pub powers_of_gamma_g: Vec<E::G1Affine>,

    // Elements in the form {\beta^i G2} for -max_degree <= i <= 0
    pub neg_powers_of_h: Vec<E::G2Affine>,
    pub beta_h: E::G2Affine
}

impl<E: PairingEngine> PCUniversalParams for UniversalParams<E> {
    fn max_degree(&self) -> usize {
        self.powers_of_g.len() - 1
    }
}

pub struct CommitterKey<E: PairingEngine> {
    pub powers_of_g: Vec<E::G1Affine>,
    pub powers_of_gamma_g: Vec<E::G1Affine>,

    pub shifted_powers_of_g: Option<Vec<E::G1Affine>>,
    pub shifted_powers_of_gamma_g: Option<Vec<E::G1Affine>>,
    pub enforced_degree_bounds: Option<Vec<usize>>,

    pub max_degree: usize,
}

impl<E: PairingEngine> CommitterKey<E> {
    pub fn powers<'a>(&'a self) -> kzg10::Powers<'a, E> {
        kzg10::Powers {
            powers_of_g: self.powers_of_g.as_slice().into(),
            powers_of_gamma_g: self.powers_of_gamma_g.as_slice().into(),
        }
    }

    pub fn shifted_powers<'a>(
        &'a self,
        degree_bound: impl Into<Option<usize>>
    ) -> Option<kzg10::Powers<'a, E>> {
        self.shifted_powers_of_g.as_ref().map(|shifted_powers_of_g| {
            self.shifted_powers_of_gamma_g.as_ref().map(|shifted_powers_of_gamma_g| {
                let powers_range = if let Some(degree_bound) = degree_bound.into() {
                    assert!(self.enforced_degree_bounds.as_ref().unwrap().contains(&degree_bound));
                    let max_bound = self.enforced_degree_bounds.as_ref().unwrap().last().unwrap();
                    (max_bound - degree_bound)..
                } else {
                    0..
                };

                let ck = kzg10::Powers {
                    powers_of_g: (&shifted_powers_of_g[powers_range]).into(),
                    powers_of_gamma_g: (&shifted_powers_of_gamma_g[powers_range]).into(),
                };

                ck
            })
        })
    }
}

impl<E: PairingEngine> PCCommitterKey for CommitterKey<E> {
    fn max_degree(&self) -> usize {
        self.max_degree
    }

    fn supported_degree(&self) -> usize {
        self.powers_of_g.len() - 1
    }
}

#[derive(Derivative)]
#[derivative(Default(bound = ""), Clone(bound = ""), Debug(bound = ""))]
pub struct VerifierKey<E: PairingEngine> {
    pub g: E::G1Affine,
    pub gamma_g: E::G1Affine,
    pub h: E::G2Affine,
    pub degree_bounds_and_neg_powers_of_h: Option<Vec<(usize, E::G2Affine)>>,

    pub supported_degree: usize,
    pub max_degree: usize,
}

impl<E: PairingEngine> VerifierKey<E> {

}

impl<E: PairingEngine> PCVerifierKey for VerifierKey<E> {
    fn max_degree(&self) -> usize {
        self.max_degree
    }

    fn supported_degree(&self) -> usize {
        self.supported_degree
    }
}

/// Commitment to a polynomial that optionally enforces a degree bound.
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
pub struct Commitment<E: PairingEngine> {
    pub(crate) comm: kzg10::Commitment<E>,
    pub(crate) has_degree_bound: bool
}

impl<E: PairingEngine> ToBytes for Commitment<E> {

}

impl<E: PairingEngine> PCCommitment for Commitment<E> {
    #[inline]
    fn empty() -> Self {
        Self {
            comm: kzg10::Commitment::empty(),
            has_degree_bound: false,
        }
    }

    fn has_degree_bound(&self) -> bool {
        self.has_degree_bound
    }

    fn size_in_bytes(&self) -> usize {

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
pub struct Randomness<E: PairingEngine> {
    pub(crate) rand: kzg10::Randomness<E>,
}

impl<'a, E: PairingEngine> Add<&'a Self> for Randomness<E> {
    type Output = Self;

    fn add(mut self, other: &'a Self) -> Self {
        self += other;
        self
    }
}

impl<'a, E: PairingEngine> AddAssign<&'a Self> for Randomness<E> {
    #[inline]
    fn add_assign(&mut self, other: &'a Self) {
        self.rand += &other.rand;
    }
}

impl<'a, E: PairingEngine> Add<(E::Fr, &'a Randomness<E>)> for Randomness<E> {
    type Output = Self;

    #[inline]
    fn add(mut self, other: (E::Fr, &'a Randomness<E>)) -> Self {
        self += other;
        self
    }
}

impl<'a, E: PairingEngine> AddAssign<(E::Fr, &'a Randomness<E>)> for Randomness<E> {
    #[inline]
    fn add_assign(&mut self, (f, other): (E::Fr, &'a Randomness<E>)) {
        self.rand += (f, &other.rand);
    }
}

impl<E: PairingEngine> PCRandomness for Randomness<E> {
    fn empty() -> Self {
        Self {
            rand: kzg10::Randomness::empty(),
        }
    }

    fn rand<R: RngCore>(hiding_bound: usize, has_degree_bound: bool, rng: &mut R) -> Self {
        Self {
            rand: kzg10::Randomness::rand(hiding_bound, false, rng),
        }
    }
}

/// Evaluation proof output by `MultiPCFromSinglePC::batch_open`.
#[derive(Derivative)]
#[derivative(
    Default(bound = ""),
    Hash(bound = ""),
    Clone(bound = ""),
    Debug(bound = ""),
    PartialEq(bound = ""),
    Eq(bound = "")
)]
pub struct BatchProof<E: PairingEngine>(pub(crate) Vec<kzg10::Proof<E>>);
