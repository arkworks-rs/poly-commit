use crate::*;
use ark_ec::{AffineCurve, PairingEngine, ProjectiveCurve};
use ark_ff::{PrimeField, ToBytes, ToConstraintField, Zero};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, SerializationError};
use ark_std::{
    borrow::Cow,
    io::{Read, Write},
    marker::PhantomData,
    ops::{Add, AddAssign},
};

/// `UniversalParams` are the universal parameters for the KZG10 scheme.
#[derive(Derivative)]
#[derivative(Clone(bound = ""), Debug(bound = ""))]
pub struct UniversalParams<E: PairingEngine> {
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
    #[derivative(Debug = "ignore")]
    pub prepared_h: E::G2Prepared,
    /// \beta times the above generator of G2, prepared for use in pairings.
    #[derivative(Debug = "ignore")]
    pub prepared_beta_h: E::G2Prepared,
}

impl<E: PairingEngine> PCUniversalParams for UniversalParams<E> {
    fn max_degree(&self) -> usize {
        self.powers_of_g.len() - 1
    }
}

impl<E: PairingEngine> CanonicalSerialize for UniversalParams<E> {
    fn serialize<W: Write>(&self, mut writer: W) -> Result<(), SerializationError> {
        self.powers_of_g.serialize(&mut writer)?;
        self.powers_of_gamma_g.serialize(&mut writer)?;
        self.h.serialize(&mut writer)?;
        self.beta_h.serialize(&mut writer)?;
        self.neg_powers_of_h.serialize(&mut writer)
    }

    fn serialized_size(&self) -> usize {
        self.powers_of_g.serialized_size()
            + self.powers_of_gamma_g.serialized_size()
            + self.h.serialized_size()
            + self.beta_h.serialized_size()
            + self.neg_powers_of_h.serialized_size()
    }

    fn serialize_unchecked<W: Write>(&self, mut writer: W) -> Result<(), SerializationError> {
        self.powers_of_g.serialize_unchecked(&mut writer)?;
        self.powers_of_gamma_g.serialize_unchecked(&mut writer)?;
        self.h.serialize_unchecked(&mut writer)?;
        self.beta_h.serialize_unchecked(&mut writer)?;
        self.neg_powers_of_h.serialize_unchecked(&mut writer)
    }

    fn serialize_uncompressed<W: Write>(&self, mut writer: W) -> Result<(), SerializationError> {
        self.powers_of_g.serialize_uncompressed(&mut writer)?;
        self.powers_of_gamma_g.serialize_uncompressed(&mut writer)?;
        self.h.serialize_uncompressed(&mut writer)?;
        self.beta_h.serialize_uncompressed(&mut writer)?;
        self.neg_powers_of_h.serialize_uncompressed(&mut writer)
    }

    fn uncompressed_size(&self) -> usize {
        self.powers_of_g.uncompressed_size()
            + self.powers_of_gamma_g.uncompressed_size()
            + self.h.uncompressed_size()
            + self.beta_h.uncompressed_size()
            + self.neg_powers_of_h.uncompressed_size()
    }
}

impl<E: PairingEngine> CanonicalDeserialize for UniversalParams<E> {
    fn deserialize<R: Read>(mut reader: R) -> Result<Self, SerializationError> {
        let powers_of_g = Vec::<E::G1Affine>::deserialize(&mut reader)?;
        let powers_of_gamma_g = BTreeMap::<usize, E::G1Affine>::deserialize(&mut reader)?;
        let h = E::G2Affine::deserialize(&mut reader)?;
        let beta_h = E::G2Affine::deserialize(&mut reader)?;
        let neg_powers_of_h = BTreeMap::<usize, E::G2Affine>::deserialize(&mut reader)?;

        let prepared_h = E::G2Prepared::from(h.clone());
        let prepared_beta_h = E::G2Prepared::from(beta_h.clone());

        Ok(Self {
            powers_of_g,
            powers_of_gamma_g,
            h,
            beta_h,
            neg_powers_of_h,
            prepared_h,
            prepared_beta_h,
        })
    }

    fn deserialize_uncompressed<R: Read>(mut reader: R) -> Result<Self, SerializationError> {
        let powers_of_g = Vec::<E::G1Affine>::deserialize_uncompressed(&mut reader)?;
        let powers_of_gamma_g =
            BTreeMap::<usize, E::G1Affine>::deserialize_uncompressed(&mut reader)?;
        let h = E::G2Affine::deserialize_uncompressed(&mut reader)?;
        let beta_h = E::G2Affine::deserialize_uncompressed(&mut reader)?;
        let neg_powers_of_h =
            BTreeMap::<usize, E::G2Affine>::deserialize_uncompressed(&mut reader)?;

        let prepared_h = E::G2Prepared::from(h.clone());
        let prepared_beta_h = E::G2Prepared::from(beta_h.clone());

        Ok(Self {
            powers_of_g,
            powers_of_gamma_g,
            h,
            beta_h,
            neg_powers_of_h,
            prepared_h,
            prepared_beta_h,
        })
    }

    fn deserialize_unchecked<R: Read>(mut reader: R) -> Result<Self, SerializationError> {
        let powers_of_g = Vec::<E::G1Affine>::deserialize_unchecked(&mut reader)?;
        let powers_of_gamma_g = BTreeMap::<usize, E::G1Affine>::deserialize_unchecked(&mut reader)?;
        let h = E::G2Affine::deserialize_unchecked(&mut reader)?;
        let beta_h = E::G2Affine::deserialize_unchecked(&mut reader)?;
        let neg_powers_of_h = BTreeMap::<usize, E::G2Affine>::deserialize_unchecked(&mut reader)?;

        let prepared_h = E::G2Prepared::from(h.clone());
        let prepared_beta_h = E::G2Prepared::from(beta_h.clone());

        Ok(Self {
            powers_of_g,
            powers_of_gamma_g,
            h,
            beta_h,
            neg_powers_of_h,
            prepared_h,
            prepared_beta_h,
        })
    }
}

/// `Powers` is used to commit to and create evaluation proofs for a given
/// polynomial.
#[derive(Derivative)]
#[derivative(
    Default(bound = ""),
    Hash(bound = ""),
    Clone(bound = ""),
    Debug(bound = "")
)]
pub struct Powers<'a, E: PairingEngine> {
    /// Group elements of the form `β^i G`, for different values of `i`.
    pub powers_of_g: Cow<'a, [E::G1Affine]>,
    /// Group elements of the form `β^i γG`, for different values of `i`.
    pub powers_of_gamma_g: Cow<'a, [E::G1Affine]>,
}

impl<E: PairingEngine> Powers<'_, E> {
    /// The number of powers in `self`.
    pub fn size(&self) -> usize {
        self.powers_of_g.len()
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
    /// \beta times the above generator of G2.
    pub beta_h: E::G2Affine,
    /// The generator of G2, prepared for use in pairings.
    #[derivative(Debug = "ignore")]
    pub prepared_h: E::G2Prepared,
    /// \beta times the above generator of G2, prepared for use in pairings.
    #[derivative(Debug = "ignore")]
    pub prepared_beta_h: E::G2Prepared,
}

impl<E: PairingEngine> CanonicalSerialize for VerifierKey<E> {
    fn serialize<W: Write>(&self, mut writer: W) -> Result<(), SerializationError> {
        self.g.serialize(&mut writer)?;
        self.gamma_g.serialize(&mut writer)?;
        self.h.serialize(&mut writer)?;
        self.beta_h.serialize(&mut writer)
    }

    fn serialized_size(&self) -> usize {
        self.g.serialized_size()
            + self.gamma_g.serialized_size()
            + self.h.serialized_size()
            + self.beta_h.serialized_size()
    }

    fn serialize_uncompressed<W: Write>(&self, mut writer: W) -> Result<(), SerializationError> {
        self.g.serialize_uncompressed(&mut writer)?;
        self.gamma_g.serialize_uncompressed(&mut writer)?;
        self.h.serialize_uncompressed(&mut writer)?;
        self.beta_h.serialize_uncompressed(&mut writer)
    }

    fn serialize_unchecked<W: Write>(&self, mut writer: W) -> Result<(), SerializationError> {
        self.g.serialize_unchecked(&mut writer)?;
        self.gamma_g.serialize_unchecked(&mut writer)?;
        self.h.serialize_unchecked(&mut writer)?;
        self.beta_h.serialize_unchecked(&mut writer)
    }

    fn uncompressed_size(&self) -> usize {
        self.g.uncompressed_size()
            + self.gamma_g.uncompressed_size()
            + self.h.uncompressed_size()
            + self.beta_h.uncompressed_size()
    }
}

impl<E: PairingEngine> CanonicalDeserialize for VerifierKey<E> {
    fn deserialize<R: Read>(mut reader: R) -> Result<Self, SerializationError> {
        let g = E::G1Affine::deserialize(&mut reader)?;
        let gamma_g = E::G1Affine::deserialize(&mut reader)?;
        let h = E::G2Affine::deserialize(&mut reader)?;
        let beta_h = E::G2Affine::deserialize(&mut reader)?;

        let prepared_h = E::G2Prepared::from(h.clone());
        let prepared_beta_h = E::G2Prepared::from(beta_h.clone());

        Ok(Self {
            g,
            gamma_g,
            h,
            beta_h,
            prepared_h,
            prepared_beta_h,
        })
    }

    fn deserialize_uncompressed<R: Read>(mut reader: R) -> Result<Self, SerializationError> {
        let g = E::G1Affine::deserialize_uncompressed(&mut reader)?;
        let gamma_g = E::G1Affine::deserialize_uncompressed(&mut reader)?;
        let h = E::G2Affine::deserialize_uncompressed(&mut reader)?;
        let beta_h = E::G2Affine::deserialize_uncompressed(&mut reader)?;

        let prepared_h = E::G2Prepared::from(h.clone());
        let prepared_beta_h = E::G2Prepared::from(beta_h.clone());

        Ok(Self {
            g,
            gamma_g,
            h,
            beta_h,
            prepared_h,
            prepared_beta_h,
        })
    }

    fn deserialize_unchecked<R: Read>(mut reader: R) -> Result<Self, SerializationError> {
        let g = E::G1Affine::deserialize_unchecked(&mut reader)?;
        let gamma_g = E::G1Affine::deserialize_unchecked(&mut reader)?;
        let h = E::G2Affine::deserialize_unchecked(&mut reader)?;
        let beta_h = E::G2Affine::deserialize_unchecked(&mut reader)?;

        let prepared_h = E::G2Prepared::from(h.clone());
        let prepared_beta_h = E::G2Prepared::from(beta_h.clone());

        Ok(Self {
            g,
            gamma_g,
            h,
            beta_h,
            prepared_h,
            prepared_beta_h,
        })
    }
}

impl<E: PairingEngine> ToBytes for VerifierKey<E> {
    #[inline]
    fn write<W: Write>(&self, mut writer: W) -> ark_std::io::Result<()> {
        self.g.write(&mut writer)?;
        self.gamma_g.write(&mut writer)?;
        self.h.write(&mut writer)?;
        self.beta_h.write(&mut writer)?;
        self.prepared_h.write(&mut writer)?;
        self.prepared_beta_h.write(&mut writer)
    }
}

impl<E: PairingEngine> ToConstraintField<<E::Fq as Field>::BasePrimeField> for VerifierKey<E>
where
    E::G1Affine: ToConstraintField<<E::Fq as Field>::BasePrimeField>,
    E::G2Affine: ToConstraintField<<E::Fq as Field>::BasePrimeField>,
{
    fn to_field_elements(&self) -> Option<Vec<<E::Fq as Field>::BasePrimeField>> {
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
pub struct PreparedVerifierKey<E: PairingEngine> {
    /// The generator of G1, prepared for power series.
    pub prepared_g: Vec<E::G1Affine>,
    /// The generator of G2, prepared for use in pairings.
    pub prepared_h: E::G2Prepared,
    /// \beta times the above generator of G2, prepared for use in pairings.
    pub prepared_beta_h: E::G2Prepared,
}

impl<E: PairingEngine> PreparedVerifierKey<E> {
    /// prepare `PreparedVerifierKey` from `VerifierKey`
    pub fn prepare(vk: &VerifierKey<E>) -> Self {
        let supported_bits = E::Fr::size_in_bits();

        let mut prepared_g = Vec::<E::G1Affine>::new();
        let mut g = E::G1Projective::from(vk.g.clone());
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
pub struct Commitment<E: PairingEngine>(
    /// The commitment is a group element.
    pub E::G1Affine,
);

impl<E: PairingEngine> ToBytes for Commitment<E> {
    #[inline]
    fn write<W: Write>(&self, writer: W) -> ark_std::io::Result<()> {
        self.0.write(writer)
    }
}

impl<E: PairingEngine> PCCommitment for Commitment<E> {
    #[inline]
    fn empty() -> Self {
        Commitment(E::G1Affine::zero())
    }

    fn has_degree_bound(&self) -> bool {
        false
    }

    fn size_in_bytes(&self) -> usize {
        ark_ff::to_bytes![E::G1Affine::zero()].unwrap().len() / 2
    }
}

impl<E: PairingEngine> ToConstraintField<<E::Fq as Field>::BasePrimeField> for Commitment<E>
where
    E::G1Affine: ToConstraintField<<E::Fq as Field>::BasePrimeField>,
{
    fn to_field_elements(&self) -> Option<Vec<<E::Fq as Field>::BasePrimeField>> {
        self.0.to_field_elements()
    }
}

impl<'a, E: PairingEngine> AddAssign<(E::Fr, &'a Commitment<E>)> for Commitment<E> {
    #[inline]
    fn add_assign(&mut self, (f, other): (E::Fr, &'a Commitment<E>)) {
        let mut other = other.0.mul(f.into_repr());
        other.add_assign_mixed(&self.0);
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
pub struct PreparedCommitment<E: PairingEngine>(
    /// The commitment is a group element.
    pub Vec<E::G1Affine>,
);

impl<E: PairingEngine> PreparedCommitment<E> {
    /// prepare `PreparedCommitment` from `Commitment`
    pub fn prepare(comm: &Commitment<E>) -> Self {
        let mut prepared_comm = Vec::<E::G1Affine>::new();
        let mut cur = E::G1Projective::from(comm.0.clone());

        let supported_bits = E::Fr::size_in_bits();

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
pub struct Randomness<F: PrimeField, P: UVPolynomial<F>> {
    /// For KZG10, the commitment randomness is a random polynomial.
    pub blinding_polynomial: P,
    _field: PhantomData<F>,
}

impl<F: PrimeField, P: UVPolynomial<F>> Randomness<F, P> {
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

impl<F: PrimeField, P: UVPolynomial<F>> PCRandomness for Randomness<F, P> {
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

impl<'a, F: PrimeField, P: UVPolynomial<F>> Add<&'a Randomness<F, P>> for Randomness<F, P> {
    type Output = Self;

    #[inline]
    fn add(mut self, other: &'a Self) -> Self {
        self.blinding_polynomial += &other.blinding_polynomial;
        self
    }
}

impl<'a, F: PrimeField, P: UVPolynomial<F>> Add<(F, &'a Randomness<F, P>)> for Randomness<F, P> {
    type Output = Self;

    #[inline]
    fn add(mut self, other: (F, &'a Randomness<F, P>)) -> Self {
        self += other;
        self
    }
}

impl<'a, F: PrimeField, P: UVPolynomial<F>> AddAssign<&'a Randomness<F, P>> for Randomness<F, P> {
    #[inline]
    fn add_assign(&mut self, other: &'a Self) {
        self.blinding_polynomial += &other.blinding_polynomial;
    }
}

impl<'a, F: PrimeField, P: UVPolynomial<F>> AddAssign<(F, &'a Randomness<F, P>)>
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
pub struct Proof<E: PairingEngine> {
    /// This is a commitment to the witness polynomial; see [KZG10] for more details.
    pub w: E::G1Affine,
    /// This is the evaluation of the random polynomial at the point for which
    /// the evaluation proof was produced.
    pub random_v: Option<E::Fr>,
}

impl<E: PairingEngine> PCProof for Proof<E> {
    fn size_in_bytes(&self) -> usize {
        let hiding_size = if self.random_v.is_some() {
            ark_ff::to_bytes![E::Fr::zero()].unwrap().len()
        } else {
            0
        };
        ark_ff::to_bytes![E::G1Affine::zero()].unwrap().len() / 2 + hiding_size
    }
}

impl<E: PairingEngine> ToBytes for Proof<E> {
    #[inline]
    fn write<W: Write>(&self, mut writer: W) -> ark_std::io::Result<()> {
        self.w.write(&mut writer)?;
        self.random_v
            .as_ref()
            .unwrap_or(&E::Fr::zero())
            .write(&mut writer)
    }
}
