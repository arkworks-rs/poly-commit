use crate::kzg10;
use crate::{
    BTreeMap, PCCommitterKey, PCPreparedCommitment, PCPreparedVerifierKey, PCVerifierKey, Vec,
};
use ark_ec::{PairingEngine, ProjectiveCurve};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, SerializationError};
use ark_std::io::{Read, Write};

/// `UniversalParams` are the universal parameters for the KZG10 scheme.
pub type UniversalParams<E> = kzg10::UniversalParams<E>;

/// `Randomness` is the randomness for the KZG10 scheme.
pub type Randomness<F, P> = kzg10::Randomness<F, P>;

/// `Commitment` is the commitment for the KZG10 scheme.
pub type Commitment<E> = kzg10::Commitment<E>;

/// `PreparedCommitment` is the prepared commitment for the KZG10 scheme.
pub type PreparedCommitment<E> = kzg10::PreparedCommitment<E>;

impl<E: PairingEngine> PCPreparedCommitment<Commitment<E>> for PreparedCommitment<E> {
    /// prepare `PreparedCommitment` from `Commitment`
    fn prepare(comm: &Commitment<E>) -> Self {
        let mut prepared_comm = Vec::<E::G1Affine>::new();
        let mut cur = E::G1Projective::from(comm.0.clone());
        for _ in 0..128 {
            prepared_comm.push(cur.clone().into());
            cur.double_in_place();
        }

        Self { 0: prepared_comm }
    }
}

/// `ComitterKey` is used to commit to, and create evaluation proofs for, a given
/// polynomial.
#[derive(Derivative, CanonicalSerialize, CanonicalDeserialize)]
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
    pub shifted_powers_of_gamma_g: Option<BTreeMap<usize, Vec<E::G1Affine>>>,

    /// The degree bounds that are supported by `self`.
    /// Sorted in ascending order from smallest bound to largest bound.
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
        degree_bound: impl Into<Option<usize>>,
    ) -> Option<kzg10::Powers<E>> {
        match (&self.shifted_powers_of_g, &self.shifted_powers_of_gamma_g) {
            (Some(shifted_powers_of_g), Some(shifted_powers_of_gamma_g)) => {
                let max_bound = self
                    .enforced_degree_bounds
                    .as_ref()
                    .unwrap()
                    .last()
                    .unwrap();
                let (bound, powers_range) = if let Some(degree_bound) = degree_bound.into() {
                    assert!(self
                        .enforced_degree_bounds
                        .as_ref()
                        .unwrap()
                        .contains(&degree_bound));
                    (degree_bound, (max_bound - degree_bound)..)
                } else {
                    (*max_bound, 0..)
                };

                let ck = kzg10::Powers {
                    powers_of_g: shifted_powers_of_g[powers_range.clone()].into(),
                    powers_of_gamma_g: shifted_powers_of_gamma_g[&bound].clone().into(),
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

    /// The generator of G2.
    pub h: E::G2Affine,

    /// \beta times the generator of G2.
    pub beta_h: E::G2Affine,

    /// The generator of G2, prepared for use in pairings.
    pub prepared_h: E::G2Prepared,

    /// The \beta times the generator of G2, prepared for use in pairings.
    pub prepared_beta_h: E::G2Prepared,

    /// Pairs a degree_bound with its corresponding G2 element, which has been prepared for use in pairings.
    /// Each pair is in the form `(degree_bound, \beta^{degree_bound - max_degree} h),` where `h` is the generator of G2 above
    pub degree_bounds_and_neg_powers_of_h: Option<Vec<(usize, E::G2Affine)>>,

    /// The maximum degree supported by the trimmed parameters that `self` is
    /// a part of.
    pub supported_degree: usize,

    /// The maximum degree supported by the `UniversalParams` `self` was derived
    /// from.
    pub max_degree: usize,
}

impl<E: PairingEngine> VerifierKey<E> {
    /// Find the appropriate shift for the degree bound.
    pub fn get_shift_power(&self, degree_bound: usize) -> Option<E::G2Prepared> {
        self.degree_bounds_and_neg_powers_of_h
            .as_ref()
            .and_then(|v| {
                v.binary_search_by(|(d, _)| d.cmp(&degree_bound))
                    .ok()
                    .map(|i| v[i].1.clone().into())
            })
    }
}

impl<E: PairingEngine> CanonicalSerialize for VerifierKey<E> {
    fn serialize<W: Write>(&self, mut writer: W) -> Result<(), SerializationError> {
        self.g.serialize(&mut writer)?;
        self.gamma_g.serialize(&mut writer)?;
        self.h.serialize(&mut writer)?;
        self.beta_h.serialize(&mut writer)?;
        self.degree_bounds_and_neg_powers_of_h
            .serialize(&mut writer)?;
        self.supported_degree.serialize(&mut writer)?;
        self.max_degree.serialize(&mut writer)
    }

    fn serialized_size(&self) -> usize {
        self.g.serialized_size()
            + self.gamma_g.serialized_size()
            + self.h.serialized_size()
            + self.beta_h.serialized_size()
            + self.degree_bounds_and_neg_powers_of_h.serialized_size()
            + self.supported_degree.serialized_size()
            + self.max_degree.serialized_size()
    }

    fn serialize_uncompressed<W: Write>(&self, mut writer: W) -> Result<(), SerializationError> {
        self.g.serialize_uncompressed(&mut writer)?;
        self.gamma_g.serialize_uncompressed(&mut writer)?;
        self.h.serialize_uncompressed(&mut writer)?;
        self.beta_h.serialize_uncompressed(&mut writer)?;
        self.degree_bounds_and_neg_powers_of_h
            .serialize_uncompressed(&mut writer)?;
        self.supported_degree.serialize_uncompressed(&mut writer)?;
        self.max_degree.serialize_uncompressed(&mut writer)
    }

    fn serialize_unchecked<W: Write>(&self, mut writer: W) -> Result<(), SerializationError> {
        self.g.serialize_unchecked(&mut writer)?;
        self.gamma_g.serialize_unchecked(&mut writer)?;
        self.h.serialize_unchecked(&mut writer)?;
        self.beta_h.serialize_unchecked(&mut writer)?;
        self.degree_bounds_and_neg_powers_of_h
            .serialize_unchecked(&mut writer)?;
        self.supported_degree.serialize_unchecked(&mut writer)?;
        self.max_degree.serialize_unchecked(&mut writer)
    }

    fn uncompressed_size(&self) -> usize {
        self.g.uncompressed_size()
            + self.gamma_g.uncompressed_size()
            + self.h.uncompressed_size()
            + self.beta_h.uncompressed_size()
            + self.degree_bounds_and_neg_powers_of_h.uncompressed_size()
            + self.supported_degree.uncompressed_size()
            + self.max_degree.uncompressed_size()
    }
}

impl<E: PairingEngine> CanonicalDeserialize for VerifierKey<E> {
    fn deserialize<R: Read>(mut reader: R) -> Result<Self, SerializationError> {
        let g = E::G1Affine::deserialize(&mut reader)?;
        let gamma_g = E::G1Affine::deserialize(&mut reader)?;
        let h = E::G2Affine::deserialize(&mut reader)?;
        let beta_h = E::G2Affine::deserialize(&mut reader)?;
        let degree_bounds_and_neg_powers_of_h =
            Option::<Vec<(usize, E::G2Affine)>>::deserialize(&mut reader)?;
        let supported_degree = usize::deserialize(&mut reader)?;
        let max_degree = usize::deserialize(&mut reader)?;

        let prepared_h = E::G2Prepared::from(h.clone());
        let prepared_beta_h = E::G2Prepared::from(beta_h.clone());

        Ok(Self {
            g,
            gamma_g,
            h,
            beta_h,
            prepared_h,
            prepared_beta_h,
            degree_bounds_and_neg_powers_of_h,
            supported_degree,
            max_degree,
        })
    }

    fn deserialize_uncompressed<R: Read>(mut reader: R) -> Result<Self, SerializationError> {
        let g = E::G1Affine::deserialize_uncompressed(&mut reader)?;
        let gamma_g = E::G1Affine::deserialize_uncompressed(&mut reader)?;
        let h = E::G2Affine::deserialize_uncompressed(&mut reader)?;
        let beta_h = E::G2Affine::deserialize_uncompressed(&mut reader)?;
        let degree_bounds_and_neg_powers_of_h =
            Option::<Vec<(usize, E::G2Affine)>>::deserialize_uncompressed(&mut reader)?;
        let supported_degree = usize::deserialize_uncompressed(&mut reader)?;
        let max_degree = usize::deserialize_uncompressed(&mut reader)?;

        let prepared_h = E::G2Prepared::from(h.clone());
        let prepared_beta_h = E::G2Prepared::from(beta_h.clone());

        Ok(Self {
            g,
            gamma_g,
            h,
            beta_h,
            prepared_h,
            prepared_beta_h,
            degree_bounds_and_neg_powers_of_h,
            supported_degree,
            max_degree,
        })
    }

    fn deserialize_unchecked<R: Read>(mut reader: R) -> Result<Self, SerializationError> {
        let g = E::G1Affine::deserialize_unchecked(&mut reader)?;
        let gamma_g = E::G1Affine::deserialize_unchecked(&mut reader)?;
        let h = E::G2Affine::deserialize_unchecked(&mut reader)?;
        let beta_h = E::G2Affine::deserialize_unchecked(&mut reader)?;
        let degree_bounds_and_neg_powers_of_h =
            Option::<Vec<(usize, E::G2Affine)>>::deserialize_unchecked(&mut reader)?;
        let supported_degree = usize::deserialize_unchecked(&mut reader)?;
        let max_degree = usize::deserialize_unchecked(&mut reader)?;

        let prepared_h = E::G2Prepared::from(h.clone());
        let prepared_beta_h = E::G2Prepared::from(beta_h.clone());

        Ok(Self {
            g,
            gamma_g,
            h,
            beta_h,
            prepared_h,
            prepared_beta_h,
            degree_bounds_and_neg_powers_of_h,
            supported_degree,
            max_degree,
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

/// Nothing to do to prepare this verifier key (for now).
pub type PreparedVerifierKey<E> = VerifierKey<E>;

impl<E: PairingEngine> PCPreparedVerifierKey<VerifierKey<E>> for PreparedVerifierKey<E> {
    /// prepare `PreparedVerifierKey` from `VerifierKey`
    fn prepare(vk: &VerifierKey<E>) -> Self {
        vk.clone()
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
