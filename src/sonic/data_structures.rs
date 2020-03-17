use crate::{PCCommitment, PCCommitterKey, PCRandomness, PCVerifierKey, Vec};
use algebra_core::{PairingEngine, ToBytes};
use core::ops::{Add, AddAssign};
use rand_core::RngCore;

use crate::kzg10;
use core::mem;

pub type UniversalParams<E> = kzg10::UniversalParams<E>;
pub type Randomness<E> = kzg10::Randomness<E>;

/// `ComitterKey` is used to commit to and create evaluation proofs for a given
/// polynomial.
#[derive(Derivative)]
#[derivative(
Default(bound = ""),
Hash(bound = ""),
Clone(bound = ""),
Debug(bound = "")
)]
pub struct CommitterKey<E: PairingEngine> {
    /// The key used to commit to polynomials.
    pub powers_of_g: Vec<E::G1Affine>,

    /// The key used to commit to hiding polynomials.
    pub powers_of_gamma_g: Vec<E::G1Affine>,

    /// The powers used to commit to shifted polynomials.
    /// This is `None` if `self` does not support enforcing any degree bounds.
    pub shifted_powers_of_g: Option<Vec<E::G1Affine>>,

    /// The powers used to commit to shifted hiding polynomials.
    /// This is `None` if `self` does not support enforcing any degree bounds.
    pub shifted_powers_of_gamma_g: Option<Vec<E::G1Affine>>,

    /// The degree bounds that are supported by `self`.
    /// In ascending order from smallest to largest.
    /// This is `None` if `self` does not support enforcing any degree bounds.
    pub enforced_degree_bounds: Option<Vec<usize>>,

    /// The maximum degree supported by the `UniversalParams` from which `self` was derived
    pub max_degree: usize,
}

impl<E: PairingEngine> CommitterKey<E> {
    /// Obtain powers for the underlying KZG10 construction
    pub fn powers(&self) -> kzg10::Powers<E> {
        kzg10::Powers {
            powers_of_g: self.powers_of_g.as_slice().into(),
            powers_of_gamma_g: self.powers_of_gamma_g.as_slice().into(),
        }
    }

    /// Obtain powers for committing to shifted polynomials.
    pub fn shifted_powers(
        &self,
        degree_bound: impl Into<Option<usize>>
    ) -> Option<kzg10::Powers<E>> {
        match (&self.shifted_powers_of_g, &self.shifted_powers_of_gamma_g) {
            (Some(shifted_powers_of_g), Some(shifted_powers_of_gamma_g)) => {
                let powers_range =
                    if let Some(degree_bound) = degree_bound.into() {
                        assert!(self.enforced_degree_bounds.as_ref().unwrap().contains(&degree_bound));
                        let max_bound = self.enforced_degree_bounds.as_ref().unwrap().last().unwrap();
                        (max_bound - degree_bound)..
                    } else {
                        0..
                    };

                let ck = kzg10::Powers {
                    powers_of_g: shifted_powers_of_g[powers_range.clone()].into(),
                    powers_of_gamma_g: shifted_powers_of_gamma_g[powers_range].into(),
                };

                Some(ck)
            }

            (_, _) => None,
        }
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

/// `VerifierKey` is used to check evaluation proofs for a given commitment.
#[derive(Derivative)]
#[derivative(Default(bound = ""), Clone(bound = ""), Debug(bound = ""))]
pub struct VerifierKey<E: PairingEngine> {

    /// The generator of G1.
    pub g: E::G1Affine,

    /// The generator of G1 that is used for making a commitment hiding.
    pub gamma_g: E::G1Affine,

    /// The generator of G2, prepared for use in pairings.
    pub prepared_h: E::G2Prepared,

    /// The \beta times the generator of G2, prepared for use in pairings.
    pub prepared_beta_h: E::G2Prepared,

    /// Pairs a degree_bound with its corresponding G2 element, which has been prepared for use in pairings.
    /// Each pair is in the form `(degree_bound, \beta^{degree_bound - max_degree} h),` where `h` is the generator of G2 above
    pub degree_bounds_and_prepared_neg_powers_of_h: Option<Vec<(usize, E::G2Prepared)>>,

    pub supported_degree: usize,
    pub max_degree: usize,
}

impl<E: PairingEngine> VerifierKey<E> {
    pub fn get_shift_power(&self, degree_bound: usize) -> Option<E::G2Prepared> {
        self.degree_bounds_and_prepared_neg_powers_of_h.as_ref().and_then(|v| {
            v.binary_search_by(|(d, _)| d.cmp(&degree_bound))
                .ok()
                .map(|i| v[i].1.clone())
        })
    }
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
    fn write<W: algebra_core::io::Write>(&self, mut writer: W) -> algebra_core::io::Result<()> {
        self.comm.write(&mut writer)?;
        self.has_degree_bound.write(&mut writer)
    }
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
        self.comm.size_in_bytes() + mem::size_of::<bool>()

    }
}

/// Evaluation proof at a query set.
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
