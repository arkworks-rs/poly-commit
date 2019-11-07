use rand_core::RngCore;
use algebra::{ToBytes, PairingEngine};
use crate::{PCCommitment, PCCommitterKey, PCRandomness, PCVerifierKey};

use crate::kzg10;
/// `UniversalParams` are the universal parameters for the KZG10 scheme.
pub type UniversalParams<E> = kzg10::UniversalParams<E>;

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
    pub powers: kzg10::Powers<E>,
    /// The key used to commit to shifted polynomials.
    /// This is `None` if `self` does not support enforcing any degree bounds.
    pub shifted_powers: Option<kzg10::Powers<E>>,
    /// The degree bounds that are supported by `self`.
    /// In ascending order from smallest to largest.
    /// This is `None` if `self` does not support enforcing any degree bounds.
    pub enforced_degree_bounds: Option<Vec<usize>>,
    /// The maximum degree supported by the `UniversalParams` `self` was derived
    /// from.
    pub max_degree: usize,
}

impl<E: PairingEngine> PCCommitterKey for CommitterKey<E> {
    fn max_degree(&self) -> usize {
        self.max_degree
    }

    fn supported_degree(&self) -> usize {
        self.powers.size()
    }
}

/// `VerifierKey` is used to check evaluation proofs for a given commitment.
#[derive(Derivative)]
#[derivative(Default(bound = ""), Clone(bound = ""), Debug(bound = ""))]
pub struct VerifierKey<E: PairingEngine> {
    /// The verification key for the underlying KZG10 scheme.
    pub vk: kzg10::VerifierKey<E>,
    /// Information required to enforce degree bounds. Each pair
    /// is of the form `(degree_bound, shifting_advice)`.
    /// The vector is sorted in ascending order of `degree_bound`.
    /// This is `None` if `self` does not support enforcing any degree bounds.
    pub degree_bounds_and_shift_powers: Option<Vec<(usize, E::G1Affine)>>,
    /// The maximum degree supported by the `UniversalParams` `self` was derived
    /// from.
    pub max_degree: usize,
    /// The maximum degree supported by the trimmed parameters that `self` is
    /// a part of.
    pub supported_degree: usize,
}

impl<E: PairingEngine> VerifierKey<E> {
    /// Find the appropriate shift for the degree bound.
    pub fn get_shift_power(&self, bound: usize) -> Option<E::G1Affine> {
        self.degree_bounds_and_shift_powers.as_ref().and_then(|v| {
            v.binary_search_by(|(d, _)| d.cmp(&bound)).ok().map(|i| v[i].1)
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
    pub(crate) shifted_comm: Option<kzg10::Commitment<E>>,
}

impl<E: PairingEngine> ToBytes for Commitment<E> {
    #[inline]
    fn write<W: std::io::Write>(&self, mut writer: W) -> std::io::Result<()> {
        self.comm.write(&mut writer)?;
        let shifted_exists = self.shifted_comm.is_some();
        shifted_exists.write(&mut writer)?;
        self.shifted_comm
            .as_ref()
            .unwrap_or(&kzg10::Commitment::empty())
            .write(&mut writer)
    }
}

impl<E: PairingEngine> PCCommitment for Commitment<E> {
    #[inline]
    fn empty() -> Self {
        Self {
            comm: kzg10::Commitment::empty(),
            shifted_comm: Some(kzg10::Commitment::empty()),
        }
    }

    fn has_degree_bound(&self) -> bool {
        self.shifted_comm.is_some()
    }

    fn size_in_bytes(&self) -> usize {
        self.comm.size_in_bytes() + self.shifted_comm.as_ref().map_or(0, |c| c.size_in_bytes())
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
    pub(crate) shifted_rand: Option<kzg10::Randomness<E>>,
}

impl<E: PairingEngine> PCRandomness for Randomness<E> {
    fn empty() -> Self {
        Self {
            rand: kzg10::Randomness::empty(),
            shifted_rand: Some(kzg10::Randomness::empty()),
        }
    }

    fn rand<R: RngCore>(hiding_bound: usize, rng: &mut R) -> Self {
        Self {
            rand: kzg10::Randomness::rand(hiding_bound, rng),
            shifted_rand: Some(kzg10::Randomness::rand(hiding_bound, rng)),
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
