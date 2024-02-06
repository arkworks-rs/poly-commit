use crate::*;
use ark_crypto_primitives::sponge::Absorb;
use ark_ec::pairing::Pairing;
use ark_ec::AdditiveGroup;
use ark_ec::AffineRepr;
use ark_ff::{PrimeField, ToConstraintField};
use ark_serialize::{
    CanonicalDeserialize, CanonicalSerialize, Compress, SerializationError, Valid, Validate,
};
use ark_std::{
    borrow::Cow,
    io::{Read, Write},
    marker::PhantomData,
    ops::{Add, AddAssign},
};

/// `UniversalParams` are the universal parameters for the KZG10 scheme.
#[derive(Derivative)]
#[derivative(
    Clone(bound = ""),
    Debug(bound = ""),
    PartialEq(bound = ""),
    Eq(bound = "")
)]
pub struct UniversalParams<E: Pairing> {
    /// Group elements of the form `{ \beta^i G }`, where `i` ranges from 0 to `degree`.
    pub powers_of_g: Vec<E::G1Affine>,
    /// Group elements of the form `{ \beta^i \gamma G }`, where `i` ranges from 0 to `degree`.
    pub powers_of_gamma_g: BTreeMap<usize, E::G1Affine>,
    /// The generator of G2.
    pub h: E::G2Affine,
    /// \beta times the above generator of G2.
    pub beta_h: E::G2Affine,
    /// Group elements of the form `{ \beta^i G2 }`, where `i` ranges from `0` to `-degree`.
    pub neg_powers_of_h: BTreeMap<usize, E::G2Affine>,
    /// The generator of G2, prepared for use in pairings.
    #[derivative(Debug = "ignore", PartialEq = "ignore")]
    pub prepared_h: E::G2Prepared,
    /// \beta times the above generator of G2, prepared for use in pairings.
    #[derivative(Debug = "ignore", PartialEq = "ignore")]
    pub prepared_beta_h: E::G2Prepared,
}

impl<E: Pairing> Valid for UniversalParams<E> {
    fn check(&self) -> Result<(), SerializationError> {
        self.powers_of_g.check()?;
        self.powers_of_gamma_g.check()?;
        self.h.check()?;
        self.beta_h.check()?;
        self.neg_powers_of_h.check()?;
        Ok(())
    }
}
impl<E: Pairing> PCUniversalParams for UniversalParams<E> {
    fn max_degree(&self) -> usize {
        self.powers_of_g.len() - 1
    }
}

impl<E: Pairing> CanonicalSerialize for UniversalParams<E> {
    fn serialize_with_mode<W: Write>(
        &self,
        mut writer: W,
        compress: Compress,
    ) -> Result<(), SerializationError> {
        self.powers_of_g
            .serialize_with_mode(&mut writer, compress)?;
        self.powers_of_gamma_g
            .serialize_with_mode(&mut writer, compress)?;
        self.h.serialize_with_mode(&mut writer, compress)?;
        self.beta_h.serialize_with_mode(&mut writer, compress)?;
        self.neg_powers_of_h
            .serialize_with_mode(&mut writer, compress)
    }

    fn serialized_size(&self, compress: Compress) -> usize {
        self.powers_of_g.serialized_size(compress)
            + self.powers_of_gamma_g.serialized_size(compress)
            + self.h.serialized_size(compress)
            + self.beta_h.serialized_size(compress)
            + self.neg_powers_of_h.serialized_size(compress)
    }
}

impl<E: Pairing> CanonicalDeserialize for UniversalParams<E> {
    fn deserialize_with_mode<R: Read>(
        mut reader: R,
        compress: Compress,
        validate: Validate,
    ) -> Result<Self, SerializationError> {
        let powers_of_g = Vec::deserialize_with_mode(&mut reader, compress, Validate::No)?;
        let powers_of_gamma_g =
            BTreeMap::deserialize_with_mode(&mut reader, compress, Validate::No)?;
        let h = E::G2Affine::deserialize_with_mode(&mut reader, compress, Validate::No)?;
        let beta_h = E::G2Affine::deserialize_with_mode(&mut reader, compress, Validate::No)?;
        let neg_powers_of_h = BTreeMap::deserialize_with_mode(&mut reader, compress, Validate::No)?;

        let prepared_h = E::G2Prepared::from(h.clone());
        let prepared_beta_h = E::G2Prepared::from(beta_h.clone());
        let result = Self {
            powers_of_g,
            powers_of_gamma_g,
            h,
            beta_h,
            neg_powers_of_h,
            prepared_h,
            prepared_beta_h,
        };
        if let Validate::Yes = validate {
            result.check()?;
        }

        Ok(result)
    }
}

/// `Powers` is used to commit to and create evaluation proofs for a given
/// polynomial.
#[derive(Derivative)]
#[derivative(
    Default(bound = ""),
    Hash(bound = ""),
    Clone(bound = ""),
    Debug(bound = ""),
    PartialEq
)]
pub struct Powers<'a, E: Pairing> {
    /// Group elements of the form `β^i G`, for different values of `i`.
    pub powers_of_g: Cow<'a, [E::G1Affine]>,
    /// Group elements of the form `β^i γG`, for different values of `i`.
    pub powers_of_gamma_g: Cow<'a, [E::G1Affine]>,
}

impl<E: Pairing> Powers<'_, E> {
    /// The number of powers in `self`.
    pub fn size(&self) -> usize {
        self.powers_of_g.len()
    }
}
impl<'a, E: Pairing> Valid for Powers<'a, E> {
    fn check(&self) -> Result<(), SerializationError> {
        Ok(())
    }
}
impl<'a, E: Pairing> CanonicalSerialize for Powers<'a, E> {
    fn serialize_with_mode<W: Write>(
        &self,
        mut writer: W,
        compress: Compress,
    ) -> Result<(), SerializationError> {
        self.powers_of_g
            .serialize_with_mode(&mut writer, compress)?;
        self.powers_of_gamma_g
            .serialize_with_mode(&mut writer, compress)
    }

    fn serialized_size(&self, compress: Compress) -> usize {
        self.powers_of_g.serialized_size(compress)
            + self.powers_of_gamma_g.serialized_size(compress)
    }
}

impl<'a, E: Pairing> CanonicalDeserialize for Powers<'a, E> {
    fn deserialize_with_mode<R: Read>(
        mut reader: R,
        compress: Compress,
        validate: Validate,
    ) -> Result<Self, SerializationError> {
        let powers_of_g = Vec::deserialize_with_mode(&mut reader, compress, validate)?;
        let powers_of_gamma_g = Vec::deserialize_with_mode(&mut reader, compress, validate)?;
        let result = Self {
            powers_of_g: Cow::Owned(powers_of_g),
            powers_of_gamma_g: Cow::Owned(powers_of_gamma_g),
        };
        if let Validate::Yes = validate {
            result.check()?;
        }
        Ok(result)
    }
}
/// `VerifierKey` is used to check evaluation proofs for a given commitment.
#[derive(Derivative)]
#[derivative(
    Default(bound = ""),
    Clone(bound = ""),
    Debug(bound = ""),
    PartialEq(bound = ""),
    Eq(bound = "")
)]
pub struct VerifierKey<E: Pairing> {
    /// The generator of G1.
    pub g: E::G1Affine,
    /// The generator of G1 that is used for making a commitment hiding.
    pub gamma_g: E::G1Affine,
    /// The generator of G2.
    pub h: E::G2Affine,
    /// \beta times the above generator of G2.
    pub beta_h: E::G2Affine,
    /// The generator of G2, prepared for use in pairings.
    #[derivative(Debug = "ignore", PartialEq = "ignore")]
    pub prepared_h: E::G2Prepared,
    /// \beta times the above generator of G2, prepared for use in pairings.
    #[derivative(Debug = "ignore", PartialEq = "ignore")]
    pub prepared_beta_h: E::G2Prepared,
}

impl<E: Pairing> Valid for VerifierKey<E> {
    fn check(&self) -> Result<(), SerializationError> {
        self.g.check()?;
        self.gamma_g.check()?;
        self.h.check()?;
        self.beta_h.check()?;

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
        self.beta_h.serialize_with_mode(&mut writer, compress)
    }

    fn serialized_size(&self, compress: Compress) -> usize {
        self.g.serialized_size(compress)
            + self.gamma_g.serialized_size(compress)
            + self.h.serialized_size(compress)
            + self.beta_h.serialized_size(compress)
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

        let prepared_h = E::G2Prepared::from(h.clone());
        let prepared_beta_h = E::G2Prepared::from(beta_h.clone());
        let result = Self {
            g,
            gamma_g,
            h,
            beta_h,
            prepared_h,
            prepared_beta_h,
        };
        if let Validate::Yes = validate {
            result.check()?;
        }

        Ok(result)
    }
}

impl<E: Pairing> ToConstraintField<<E::TargetField as Field>::BasePrimeField> for VerifierKey<E>
where
    E::G1Affine: ToConstraintField<<E::TargetField as Field>::BasePrimeField>,
    E::G2Affine: ToConstraintField<<E::TargetField as Field>::BasePrimeField>,
{
    fn to_field_elements(&self) -> Option<Vec<<E::TargetField as Field>::BasePrimeField>> {
        let mut res = Vec::new();

        res.extend_from_slice(&self.g.to_field_elements().unwrap());
        res.extend_from_slice(&self.gamma_g.to_field_elements().unwrap());
        res.extend_from_slice(&self.h.to_field_elements().unwrap());
        res.extend_from_slice(&self.beta_h.to_field_elements().unwrap());

        Some(res)
    }
}

/// `PreparedVerifierKey` is the fully prepared version for checking evaluation proofs for a given commitment.
/// We omit gamma here for simplicity.
#[derive(Derivative)]
#[derivative(Default(bound = ""), Clone(bound = ""), Debug(bound = ""))]
pub struct PreparedVerifierKey<E: Pairing> {
    /// The generator of G1, prepared for power series.
    pub prepared_g: Vec<E::G1Affine>,
    /// The generator of G2, prepared for use in pairings.
    pub prepared_h: E::G2Prepared,
    /// \beta times the above generator of G2, prepared for use in pairings.
    pub prepared_beta_h: E::G2Prepared,
}

impl<E: Pairing> PreparedVerifierKey<E> {
    /// prepare `PreparedVerifierKey` from `VerifierKey`
    pub fn prepare(vk: &VerifierKey<E>) -> Self {
        let supported_bits = E::ScalarField::MODULUS_BIT_SIZE as usize;

        let mut prepared_g = Vec::<E::G1Affine>::new();
        let mut g = E::G1::from(vk.g.clone());
        for _ in 0..supported_bits {
            prepared_g.push(g.clone().into());
            g.double_in_place();
        }

        Self {
            prepared_g,
            prepared_h: vk.prepared_h.clone(),
            prepared_beta_h: vk.prepared_beta_h.clone(),
        }
    }
}

/// `Commitment` commits to a polynomial. It is output by `KZG10::commit`.
#[derive(Derivative, CanonicalSerialize, CanonicalDeserialize, Absorb)]
#[derivative(
    Default(bound = ""),
    Hash(bound = ""),
    Clone(bound = ""),
    Copy(bound = ""),
    Debug(bound = ""),
    PartialEq(bound = ""),
    Eq(bound = "")
)]
pub struct Commitment<E>(
    /// The commitment is a group element.
    pub E::G1Affine,
)
where
    E: Pairing,
    E::G1Affine: Absorb;

impl<E> PCCommitment for Commitment<E>
where
    E: Pairing,
    E::G1Affine: Absorb,
{
    #[inline]
    fn empty() -> Self {
        Commitment(E::G1Affine::zero())
    }

    fn has_degree_bound(&self) -> bool {
        false
    }
}

impl<E> ToConstraintField<<E::TargetField as Field>::BasePrimeField> for Commitment<E>
where
    E::G1Affine: ToConstraintField<<E::TargetField as Field>::BasePrimeField> + Absorb,
    E: Pairing,
{
    fn to_field_elements(&self) -> Option<Vec<<E::TargetField as Field>::BasePrimeField>> {
        self.0.to_field_elements()
    }
}

impl<'a, E> AddAssign<(E::ScalarField, &'a Commitment<E>)> for Commitment<E>
where
    E: Pairing,
    E::G1Affine: Absorb,
{
    #[inline]
    fn add_assign(&mut self, (f, other): (E::ScalarField, &'a Commitment<E>)) {
        let mut other = other.0 * f;
        other.add_assign(&self.0);
        self.0 = other.into();
    }
}

/// `PreparedCommitment` commits to a polynomial and prepares for mul_bits.
#[derive(Derivative)]
#[derivative(
    Default(bound = ""),
    Hash(bound = ""),
    Clone(bound = ""),
    Debug(bound = ""),
    PartialEq(bound = ""),
    Eq(bound = "")
)]
pub struct PreparedCommitment<E: Pairing>(
    /// The commitment is a group element.
    pub Vec<E::G1Affine>,
);

impl<E> PreparedCommitment<E>
where
    E: Pairing,
    E::G1Affine: Absorb,
{
    /// prepare `PreparedCommitment` from `Commitment`
    pub fn prepare(comm: &Commitment<E>) -> Self {
        let mut prepared_comm = Vec::<E::G1Affine>::new();
        let mut cur = E::G1::from(comm.0.clone());

        let supported_bits = E::ScalarField::MODULUS_BIT_SIZE as usize;

        for _ in 0..supported_bits {
            prepared_comm.push(cur.clone().into());
            cur.double_in_place();
        }

        Self { 0: prepared_comm }
    }
}

/// `Randomness` hides the polynomial inside a commitment. It is output by `KZG10::commit`.
#[derive(Derivative, CanonicalSerialize, CanonicalDeserialize)]
#[derivative(
    Hash(bound = ""),
    Clone(bound = ""),
    Debug(bound = ""),
    PartialEq(bound = ""),
    Eq(bound = "")
)]
pub struct Randomness<F: PrimeField, P: DenseUVPolynomial<F>> {
    /// For KZG10, the commitment randomness is a random polynomial.
    pub blinding_polynomial: P,
    _field: PhantomData<F>,
}

impl<F: PrimeField, P: DenseUVPolynomial<F>> Randomness<F, P> {
    /// Does `self` provide any hiding properties to the corresponding commitment?
    /// `self.is_hiding() == true` only if the underlying polynomial is non-zero.
    #[inline]
    pub fn is_hiding(&self) -> bool {
        !self.blinding_polynomial.is_zero()
    }

    /// What is the degree of the hiding polynomial for a given hiding bound?
    #[inline]
    pub fn calculate_hiding_polynomial_degree(hiding_bound: usize) -> usize {
        hiding_bound + 1
    }
}

impl<F: PrimeField, P: DenseUVPolynomial<F>> PCCommitmentState for Randomness<F, P> {
    type Randomness = Self;
    fn empty() -> Self {
        Self {
            blinding_polynomial: P::zero(),
            _field: PhantomData,
        }
    }

    fn rand<R: RngCore>(hiding_bound: usize, _: bool, _: Option<usize>, rng: &mut R) -> Self {
        let mut randomness = Randomness::empty();
        let hiding_poly_degree = Self::calculate_hiding_polynomial_degree(hiding_bound);
        randomness.blinding_polynomial = P::rand(hiding_poly_degree, rng);
        randomness
    }
}

impl<'a, F: PrimeField, P: DenseUVPolynomial<F>> Add<&'a Randomness<F, P>> for Randomness<F, P> {
    type Output = Self;

    #[inline]
    fn add(mut self, other: &'a Self) -> Self {
        self.blinding_polynomial += &other.blinding_polynomial;
        self
    }
}

impl<'a, F: PrimeField, P: DenseUVPolynomial<F>> Add<(F, &'a Randomness<F, P>)>
    for Randomness<F, P>
{
    type Output = Self;

    #[inline]
    fn add(mut self, other: (F, &'a Randomness<F, P>)) -> Self {
        self += other;
        self
    }
}

impl<'a, F: PrimeField, P: DenseUVPolynomial<F>> AddAssign<&'a Randomness<F, P>>
    for Randomness<F, P>
{
    #[inline]
    fn add_assign(&mut self, other: &'a Self) {
        self.blinding_polynomial += &other.blinding_polynomial;
    }
}

impl<'a, F: PrimeField, P: DenseUVPolynomial<F>> AddAssign<(F, &'a Randomness<F, P>)>
    for Randomness<F, P>
{
    #[inline]
    fn add_assign(&mut self, (f, other): (F, &'a Randomness<F, P>)) {
        self.blinding_polynomial += (f, &other.blinding_polynomial);
    }
}

/// `Proof` is an evaluation proof that is output by `KZG10::open`.
#[derive(Derivative, CanonicalSerialize, CanonicalDeserialize)]
#[derivative(
    Default(bound = ""),
    Hash(bound = ""),
    Clone(bound = ""),
    Copy(bound = ""),
    Debug(bound = ""),
    PartialEq(bound = ""),
    Eq(bound = "")
)]
pub struct Proof<E: Pairing> {
    /// This is a commitment to the witness polynomial; see [KZG10] for more details.
    pub w: E::G1Affine,
    /// This is the evaluation of the random polynomial at the point for which
    /// the evaluation proof was produced.
    pub random_v: Option<E::ScalarField>,
}
