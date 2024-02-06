use crate::kzg10;
use crate::{
    BTreeMap, PCCommitterKey, PCPreparedCommitment, PCPreparedVerifierKey, PCVerifierKey, Vec,
};
use ark_crypto_primitives::sponge::Absorb;
use ark_ec::pairing::Pairing;
use ark_ec::AdditiveGroup;
use ark_serialize::{
    CanonicalDeserialize, CanonicalSerialize, Compress, SerializationError, Valid, Validate,
};
use ark_std::io::{Read, Write};

/// `UniversalParams` are the universal parameters for the KZG10 scheme.
pub type UniversalParams<E> = kzg10::UniversalParams<E>;

/// `Randomness` is the randomness for the KZG10 scheme.
pub type Randomness<F, P> = kzg10::Randomness<F, P>;

/// `Commitment` is the commitment for the KZG10 scheme.
pub type Commitment<E> = kzg10::Commitment<E>;

/// `PreparedCommitment` is the prepared commitment for the KZG10 scheme.
pub type PreparedCommitment<E> = kzg10::PreparedCommitment<E>;

impl<E> PCPreparedCommitment<Commitment<E>> for PreparedCommitment<E>
where
    E: Pairing,
    E::G1Affine: Absorb,
{
    /// prepare `PreparedCommitment` from `Commitment`
    fn prepare(comm: &Commitment<E>) -> Self {
        let mut prepared_comm = Vec::<E::G1Affine>::new();
        let mut cur = E::G1::from(comm.0.clone());
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
pub struct CommitterKey<E: Pairing> {
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

impl<E: Pairing> CommitterKey<E> {
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

impl<E: Pairing> PCCommitterKey for CommitterKey<E> {
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
pub struct VerifierKey<E: Pairing> {
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

impl<E: Pairing> VerifierKey<E> {
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

impl<E: Pairing> Valid for VerifierKey<E> {
    fn check(&self) -> Result<(), SerializationError> {
        self.g.check()?;
        self.gamma_g.check()?;
        self.h.check()?;
        self.beta_h.check()?;
        self.degree_bounds_and_neg_powers_of_h.check()?;
        if self.supported_degree > self.max_degree {
            return Err(SerializationError::InvalidData);
        }
        Ok(())
    }
}

impl<E: Pairing> CanonicalSerialize for VerifierKey<E> {
    fn serialize_with_mode<W: Write>(
        &self,
        mut writer: W,
        compress: Compress,
    ) -> Result<(), SerializationError> {
        self.g.serialize_with_mode(&mut writer, compress)?;
        self.gamma_g.serialize_with_mode(&mut writer, compress)?;
        self.h.serialize_with_mode(&mut writer, compress)?;
        self.beta_h.serialize_with_mode(&mut writer, compress)?;
        self.degree_bounds_and_neg_powers_of_h
            .serialize_with_mode(&mut writer, compress)?;
        self.supported_degree
            .serialize_with_mode(&mut writer, compress)?;
        self.max_degree.serialize_with_mode(&mut writer, compress)
    }

    fn serialized_size(&self, compress: Compress) -> usize {
        self.g.serialized_size(compress)
            + self.gamma_g.serialized_size(compress)
            + self.h.serialized_size(compress)
            + self.beta_h.serialized_size(compress)
            + self
                .degree_bounds_and_neg_powers_of_h
                .serialized_size(compress)
            + self.supported_degree.serialized_size(compress)
            + self.max_degree.serialized_size(compress)
    }
}

impl<E: Pairing> CanonicalDeserialize for VerifierKey<E> {
    fn deserialize_with_mode<R: Read>(
        mut reader: R,
        compress: Compress,
        validate: Validate,
    ) -> Result<Self, SerializationError> {
        let g = E::G1Affine::deserialize_with_mode(&mut reader, compress, Validate::No)?;
        let gamma_g = E::G1Affine::deserialize_with_mode(&mut reader, compress, Validate::No)?;
        let h = E::G2Affine::deserialize_with_mode(&mut reader, compress, Validate::No)?;
        let beta_h = E::G2Affine::deserialize_with_mode(&mut reader, compress, Validate::No)?;
        let degree_bounds_and_neg_powers_of_h =
            Option::<Vec<(usize, E::G2Affine)>>::deserialize_with_mode(
                &mut reader,
                compress,
                Validate::No,
            )?;
        let supported_degree = usize::deserialize_with_mode(&mut reader, compress, Validate::No)?;
        let max_degree = usize::deserialize_with_mode(&mut reader, compress, Validate::No)?;

        let prepared_h = E::G2Prepared::from(h.clone());
        let prepared_beta_h = E::G2Prepared::from(beta_h.clone());

        let result = Self {
            g,
            gamma_g,
            h,
            beta_h,
            prepared_h,
            prepared_beta_h,
            degree_bounds_and_neg_powers_of_h,
            supported_degree,
            max_degree,
        };

        if let Validate::Yes = validate {
            result.check()?;
        }

        Ok(result)
    }
}

impl<E: Pairing> PCVerifierKey for VerifierKey<E> {
    fn max_degree(&self) -> usize {
        self.max_degree
    }

    fn supported_degree(&self) -> usize {
        self.supported_degree
    }
}

/// Nothing to do to prepare this verifier key (for now).
pub type PreparedVerifierKey<E> = VerifierKey<E>;

impl<E: Pairing> PCPreparedVerifierKey<VerifierKey<E>> for PreparedVerifierKey<E> {
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
pub struct BatchProof<E: Pairing>(pub(crate) Vec<kzg10::Proof<E>>);
