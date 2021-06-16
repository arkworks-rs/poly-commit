use crate::{BTreeMap, Vec};
use crate::{
    PCCommitterKey, PCPreparedVerifierKey, PCProof, PCRandomness, PCUniversalParams, PCVerifierKey,
};
use ark_ec::PairingEngine;
use ark_ff::{ToBytes, Zero};
use ark_poly::MVPolynomial;
use ark_std::{
    io::{Read, Write},
    marker::PhantomData,
    ops::{Add, AddAssign, Index},
};

use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, SerializationError};
use ark_std::rand::RngCore;

/// `UniversalParams` are the universal parameters for the MarlinPST13 scheme.
#[derive(Derivative)]
#[derivative(Default(bound = ""), Clone(bound = ""), Debug(bound = ""))]
pub struct UniversalParams<E, P>
where
    E: PairingEngine,
    P: MVPolynomial<E::Fr>,
    P::Point: Index<usize, Output = E::Fr>,
{
    /// Contains group elements corresponding to all possible monomials with
    /// `num_vars` and maximum degree `max_degree` evaluated at `\beta`
    pub powers_of_g: BTreeMap<P::Term, E::G1Affine>,
    /// `\gamma` times the generater of G1
    pub gamma_g: E::G1Affine,
    /// Group elements of the form `{ \beta_i^j \gamma G }`, where `i` ranges
    /// from 0 to `num_vars-1` and `j` ranges from `1` to `max_degree+1`.
    pub powers_of_gamma_g: Vec<Vec<E::G1Affine>>,
    /// The generator of G2.
    pub h: E::G2Affine,
    /// Group elements of the form `{ \beta_i H }`, where `i` ranges from 0 to `num_vars-1`
    pub beta_h: Vec<E::G2Affine>,
    /// The generator of G2, prepared for use in pairings.
    #[derivative(Debug = "ignore")]
    pub prepared_h: E::G2Prepared,
    /// Group elements of the form `{ \beta_i H }`, where `i` ranges from 0 to `num_vars-1`,
    /// prepared for use in pairings
    #[derivative(Debug = "ignore")]
    pub prepared_beta_h: Vec<E::G2Prepared>,
    /// The number of variables `self` is initialized for
    pub num_vars: usize,
    /// The maximum degree supported by `self`
    pub max_degree: usize,
}

impl<E, P> CanonicalSerialize for UniversalParams<E, P>
where
    E: PairingEngine,
    P: MVPolynomial<E::Fr>,
    P::Point: Index<usize, Output = E::Fr>,
{
    fn serialize<W: Write>(&self, mut writer: W) -> Result<(), SerializationError> {
        self.powers_of_g.serialize(&mut writer)?;
        self.gamma_g.serialize(&mut writer)?;
        self.powers_of_gamma_g.serialize(&mut writer)?;
        self.h.serialize(&mut writer)?;
        self.beta_h.serialize(&mut writer)?;
        self.num_vars.serialize(&mut writer)?;
        self.max_degree.serialize(&mut writer)
    }

    fn serialized_size(&self) -> usize {
        self.powers_of_g.serialized_size()
            + self.gamma_g.serialized_size()
            + self.powers_of_gamma_g.serialized_size()
            + self.h.serialized_size()
            + self.beta_h.serialized_size()
            + self.num_vars.serialized_size()
            + self.max_degree.serialized_size()
    }

    fn serialize_uncompressed<W: Write>(&self, mut writer: W) -> Result<(), SerializationError> {
        self.powers_of_g.serialize_uncompressed(&mut writer)?;
        self.gamma_g.serialize_uncompressed(&mut writer)?;
        self.powers_of_gamma_g.serialize_uncompressed(&mut writer)?;
        self.h.serialize_uncompressed(&mut writer)?;
        self.beta_h.serialize_uncompressed(&mut writer)?;
        self.num_vars.serialize_uncompressed(&mut writer)?;
        self.max_degree.serialize_uncompressed(&mut writer)
    }

    fn serialize_unchecked<W: Write>(&self, mut writer: W) -> Result<(), SerializationError> {
        self.powers_of_g.serialize_unchecked(&mut writer)?;
        self.gamma_g.serialize_unchecked(&mut writer)?;
        self.powers_of_gamma_g.serialize_unchecked(&mut writer)?;
        self.h.serialize_unchecked(&mut writer)?;
        self.beta_h.serialize_unchecked(&mut writer)?;
        self.num_vars.serialize_unchecked(&mut writer)?;
        self.max_degree.serialize_unchecked(&mut writer)
    }

    fn uncompressed_size(&self) -> usize {
        self.powers_of_g.uncompressed_size()
            + self.gamma_g.uncompressed_size()
            + self.powers_of_gamma_g.uncompressed_size()
            + self.h.uncompressed_size()
            + self.beta_h.uncompressed_size()
            + self.num_vars.uncompressed_size()
            + self.max_degree.uncompressed_size()
    }
}

impl<E, P> CanonicalDeserialize for UniversalParams<E, P>
where
    E: PairingEngine,
    P: MVPolynomial<E::Fr>,
    P::Point: Index<usize, Output = E::Fr>,
{
    fn deserialize<R: Read>(mut reader: R) -> Result<Self, SerializationError> {
        let powers_of_g = BTreeMap::<P::Term, E::G1Affine>::deserialize(&mut reader)?;
        let gamma_g = E::G1Affine::deserialize(&mut reader)?;
        let powers_of_gamma_g = Vec::<Vec<E::G1Affine>>::deserialize(&mut reader)?;
        let h = E::G2Affine::deserialize(&mut reader)?;
        let beta_h = Vec::<E::G2Affine>::deserialize(&mut reader)?;
        let num_vars = usize::deserialize(&mut reader)?;
        let max_degree = usize::deserialize(&mut reader)?;

        let prepared_beta_h = beta_h.iter().map(|x| x.clone().into()).collect();
        Ok(Self {
            powers_of_g,
            gamma_g,
            powers_of_gamma_g,
            h,
            beta_h,
            prepared_h: h.into(),
            prepared_beta_h,
            num_vars,
            max_degree,
        })
    }

    fn deserialize_uncompressed<R: Read>(mut reader: R) -> Result<Self, SerializationError> {
        let powers_of_g = BTreeMap::<P::Term, E::G1Affine>::deserialize_uncompressed(&mut reader)?;
        let gamma_g = E::G1Affine::deserialize_uncompressed(&mut reader)?;
        let powers_of_gamma_g = Vec::<Vec<E::G1Affine>>::deserialize_uncompressed(&mut reader)?;
        let h = E::G2Affine::deserialize_uncompressed(&mut reader)?;
        let beta_h = Vec::<E::G2Affine>::deserialize_uncompressed(&mut reader)?;
        let num_vars = usize::deserialize_uncompressed(&mut reader)?;
        let max_degree = usize::deserialize_uncompressed(&mut reader)?;

        let prepared_beta_h = beta_h.iter().map(|x| x.clone().into()).collect();
        Ok(Self {
            powers_of_g,
            gamma_g,
            powers_of_gamma_g,
            h,
            beta_h,
            prepared_h: h.into(),
            prepared_beta_h,
            num_vars,
            max_degree,
        })
    }

    fn deserialize_unchecked<R: Read>(mut reader: R) -> Result<Self, SerializationError> {
        let powers_of_g = BTreeMap::<P::Term, E::G1Affine>::deserialize_unchecked(&mut reader)?;
        let gamma_g = E::G1Affine::deserialize_unchecked(&mut reader)?;
        let powers_of_gamma_g = Vec::<Vec<E::G1Affine>>::deserialize_unchecked(&mut reader)?;
        let h = E::G2Affine::deserialize_unchecked(&mut reader)?;
        let beta_h = Vec::<E::G2Affine>::deserialize_unchecked(&mut reader)?;
        let num_vars = usize::deserialize_unchecked(&mut reader)?;
        let max_degree = usize::deserialize_unchecked(&mut reader)?;

        let prepared_beta_h = beta_h.iter().map(|x| x.clone().into()).collect();
        Ok(Self {
            powers_of_g,
            gamma_g,
            powers_of_gamma_g,
            h,
            beta_h,
            prepared_h: h.into(),
            prepared_beta_h,
            num_vars,
            max_degree,
        })
    }
}

impl<E, P> PCUniversalParams for UniversalParams<E, P>
where
    E: PairingEngine,
    P: MVPolynomial<E::Fr>,
    P::Point: Index<usize, Output = E::Fr>,
{
    fn max_degree(&self) -> usize {
        self.max_degree
    }
}

/// `CommitterKey` is used to commit to and create evaluation proofs for a given
/// polynomial.
#[derive(Derivative, CanonicalSerialize, CanonicalDeserialize)]
#[derivative(Hash(bound = ""), Clone(bound = ""), Debug(bound = ""))]
pub struct CommitterKey<E, P>
where
    E: PairingEngine,
    P: MVPolynomial<E::Fr>,
    P::Point: Index<usize, Output = E::Fr>,
{
    /// Contains group elements corresponding to all possible monomials with
    /// `num_vars` and maximum degree `supported_degree` evaluated at `\beta`
    pub powers_of_g: BTreeMap<P::Term, E::G1Affine>,
    /// `\gamma` times the generater of G1
    pub gamma_g: E::G1Affine,
    /// Group elements of the form `{ \beta_i^j \gamma G }`, where `i` ranges
    /// from 0 to `num_vars-1` and `j` ranges from `1` to `supported_degree+1`.
    pub powers_of_gamma_g: Vec<Vec<E::G1Affine>>,
    /// The number of variables `self` is initialized for
    pub num_vars: usize,
    /// The maximum degree supported by the trimmed parameters that `self` is
    /// a part of
    pub supported_degree: usize,
    /// The maximum degree supported by the `UniversalParams` `self` was derived
    /// from.
    pub max_degree: usize,
}

impl<E, P> PCCommitterKey for CommitterKey<E, P>
where
    E: PairingEngine,
    P: MVPolynomial<E::Fr>,
    P::Point: Index<usize, Output = E::Fr>,
{
    fn max_degree(&self) -> usize {
        self.max_degree
    }

    fn supported_degree(&self) -> usize {
        self.supported_degree
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
    /// `\beta_i` times the above generator of G2 where `i` ranges from 0 to `num_vars-1`
    pub beta_h: Vec<E::G2Affine>,
    /// The generator of G2, prepared for use in pairings.
    #[derivative(Debug = "ignore")]
    pub prepared_h: E::G2Prepared,
    /// `\beta_i` times the above generator of G2 where `i` ranges from 0 to `num_vars-1`,
    /// prepared for use in pairings
    #[derivative(Debug = "ignore")]
    pub prepared_beta_h: Vec<E::G2Prepared>,
    /// The number of variables `self` is initialized for
    pub num_vars: usize,
    /// The maximum degree supported by the trimmed parameters that `self` is
    /// a part of.
    pub supported_degree: usize,
    /// The maximum degree supported by the `UniversalParams` `self` was derived
    /// from.
    pub max_degree: usize,
}

impl<E: PairingEngine> CanonicalSerialize for VerifierKey<E> {
    fn serialize<W: Write>(&self, mut writer: W) -> Result<(), SerializationError> {
        self.g.serialize(&mut writer)?;
        self.gamma_g.serialize(&mut writer)?;
        self.h.serialize(&mut writer)?;
        self.beta_h.serialize(&mut writer)?;
        self.num_vars.serialize(&mut writer)?;
        self.supported_degree.serialize(&mut writer)?;
        self.max_degree.serialize(&mut writer)
    }

    fn serialized_size(&self) -> usize {
        self.g.serialized_size()
            + self.gamma_g.serialized_size()
            + self.h.serialized_size()
            + self.beta_h.serialized_size()
            + self.num_vars.serialized_size()
            + self.supported_degree.serialized_size()
            + self.max_degree.serialized_size()
    }

    fn serialize_uncompressed<W: Write>(&self, mut writer: W) -> Result<(), SerializationError> {
        self.g.serialize_uncompressed(&mut writer)?;
        self.gamma_g.serialize_uncompressed(&mut writer)?;
        self.h.serialize_uncompressed(&mut writer)?;
        self.beta_h.serialize_uncompressed(&mut writer)?;
        self.num_vars.serialize_uncompressed(&mut writer)?;
        self.supported_degree.serialize_uncompressed(&mut writer)?;
        self.max_degree.serialize_uncompressed(&mut writer)
    }

    fn serialize_unchecked<W: Write>(&self, mut writer: W) -> Result<(), SerializationError> {
        self.g.serialize_unchecked(&mut writer)?;
        self.gamma_g.serialize_unchecked(&mut writer)?;
        self.h.serialize_unchecked(&mut writer)?;
        self.beta_h.serialize_unchecked(&mut writer)?;
        self.num_vars.serialize_unchecked(&mut writer)?;
        self.supported_degree.serialize_unchecked(&mut writer)?;
        self.max_degree.serialize_unchecked(&mut writer)
    }

    fn uncompressed_size(&self) -> usize {
        self.g.uncompressed_size()
            + self.gamma_g.uncompressed_size()
            + self.h.uncompressed_size()
            + self.beta_h.uncompressed_size()
            + self.num_vars.uncompressed_size()
            + self.supported_degree.uncompressed_size()
            + self.max_degree.uncompressed_size()
    }
}

impl<E: PairingEngine> CanonicalDeserialize for VerifierKey<E> {
    fn deserialize<R: Read>(mut reader: R) -> Result<Self, SerializationError> {
        let g = E::G1Affine::deserialize(&mut reader)?;
        let gamma_g = E::G1Affine::deserialize(&mut reader)?;
        let h = E::G2Affine::deserialize(&mut reader)?;
        let beta_h = Vec::<E::G2Affine>::deserialize(&mut reader)?;
        let num_vars = usize::deserialize(&mut reader)?;
        let supported_degree = usize::deserialize(&mut reader)?;
        let max_degree = usize::deserialize(&mut reader)?;

        let prepared_beta_h = beta_h.iter().map(|x| x.clone().into()).collect();
        Ok(Self {
            g,
            gamma_g,
            h,
            beta_h,
            prepared_h: h.into(),
            prepared_beta_h,
            num_vars,
            supported_degree,
            max_degree,
        })
    }

    fn deserialize_uncompressed<R: Read>(mut reader: R) -> Result<Self, SerializationError> {
        let g = E::G1Affine::deserialize_uncompressed(&mut reader)?;
        let gamma_g = E::G1Affine::deserialize_uncompressed(&mut reader)?;
        let h = E::G2Affine::deserialize_uncompressed(&mut reader)?;
        let beta_h = Vec::<E::G2Affine>::deserialize_uncompressed(&mut reader)?;
        let num_vars = usize::deserialize_uncompressed(&mut reader)?;
        let supported_degree = usize::deserialize_uncompressed(&mut reader)?;
        let max_degree = usize::deserialize_uncompressed(&mut reader)?;

        let prepared_beta_h = beta_h.iter().map(|x| x.clone().into()).collect();
        Ok(Self {
            g,
            gamma_g,
            h,
            beta_h,
            prepared_h: h.into(),
            prepared_beta_h,
            num_vars,
            supported_degree,
            max_degree,
        })
    }

    fn deserialize_unchecked<R: Read>(mut reader: R) -> Result<Self, SerializationError> {
        let g = E::G1Affine::deserialize_unchecked(&mut reader)?;
        let gamma_g = E::G1Affine::deserialize_unchecked(&mut reader)?;
        let h = E::G2Affine::deserialize_unchecked(&mut reader)?;
        let beta_h = Vec::<E::G2Affine>::deserialize_unchecked(&mut reader)?;
        let num_vars = usize::deserialize_unchecked(&mut reader)?;
        let supported_degree = usize::deserialize_unchecked(&mut reader)?;
        let max_degree = usize::deserialize_unchecked(&mut reader)?;

        let prepared_beta_h = beta_h.iter().map(|x| x.clone().into()).collect();
        Ok(Self {
            g,
            gamma_g,
            h,
            beta_h,
            prepared_h: h.into(),
            prepared_beta_h,
            num_vars,
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

/// `Randomness` hides the polynomial inside a commitment`.
#[derive(Derivative, CanonicalSerialize, CanonicalDeserialize)]
#[derivative(
    Hash(bound = ""),
    Clone(bound = ""),
    Debug(bound = ""),
    PartialEq(bound = ""),
    Eq(bound = "")
)]
pub struct Randomness<E, P>
where
    E: PairingEngine,
    P: MVPolynomial<E::Fr>,
    P::Point: Index<usize, Output = E::Fr>,
{
    /// A multivariate polynomial where each monomial is univariate with random coefficient
    pub blinding_polynomial: P,
    _engine: PhantomData<E>,
}

impl<E, P> Randomness<E, P>
where
    E: PairingEngine,
    P: MVPolynomial<E::Fr>,
    P::Point: Index<usize, Output = E::Fr>,
{
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

impl<E, P> PCRandomness for Randomness<E, P>
where
    E: PairingEngine,
    P: MVPolynomial<E::Fr>,
    P::Point: Index<usize, Output = E::Fr>,
{
    fn empty() -> Self {
        Self {
            blinding_polynomial: P::zero(),
            _engine: PhantomData,
        }
    }

    fn rand<R: RngCore>(
        hiding_bound: usize,
        _: bool,
        num_vars: Option<usize>,
        rng: &mut R,
    ) -> Self {
        let hiding_poly_degree = Self::calculate_hiding_polynomial_degree(hiding_bound);
        Randomness {
            blinding_polynomial: P::rand(hiding_poly_degree, num_vars.unwrap(), rng),
            _engine: PhantomData,
        }
    }
}

impl<'a, E: PairingEngine, P: MVPolynomial<E::Fr>> Add<&'a Randomness<E, P>> for Randomness<E, P>
where
    E: PairingEngine,
    P: MVPolynomial<E::Fr>,
    P::Point: Index<usize, Output = E::Fr>,
{
    type Output = Self;

    #[inline]
    fn add(mut self, other: &'a Self) -> Self {
        self.blinding_polynomial += &other.blinding_polynomial;
        self
    }
}

impl<'a, E, P> Add<(E::Fr, &'a Randomness<E, P>)> for Randomness<E, P>
where
    E: PairingEngine,
    P: MVPolynomial<E::Fr>,
    P::Point: Index<usize, Output = E::Fr>,
{
    type Output = Self;

    #[inline]
    fn add(mut self, other: (E::Fr, &'a Randomness<E, P>)) -> Self {
        self += other;
        self
    }
}

impl<'a, E, P> AddAssign<&'a Randomness<E, P>> for Randomness<E, P>
where
    E: PairingEngine,
    P: MVPolynomial<E::Fr>,
    P::Point: Index<usize, Output = E::Fr>,
{
    #[inline]
    fn add_assign(&mut self, other: &'a Self) {
        self.blinding_polynomial += &other.blinding_polynomial;
    }
}

impl<'a, E, P> AddAssign<(E::Fr, &'a Randomness<E, P>)> for Randomness<E, P>
where
    E: PairingEngine,
    P: MVPolynomial<E::Fr>,
    P::Point: Index<usize, Output = E::Fr>,
{
    #[inline]
    fn add_assign(&mut self, (f, other): (E::Fr, &'a Randomness<E, P>)) {
        self.blinding_polynomial += (f, &other.blinding_polynomial);
    }
}

/// `Proof` is an evaluation proof that is output by `KZG10::open`.
#[derive(Derivative, CanonicalSerialize, CanonicalDeserialize)]
#[derivative(
    Default(bound = ""),
    Hash(bound = ""),
    Clone(bound = ""),
    Debug(bound = ""),
    PartialEq(bound = ""),
    Eq(bound = "")
)]
pub struct Proof<E: PairingEngine> {
    /// Commitments to the witness polynomials
    pub w: Vec<E::G1Affine>,
    /// Evaluation of the random polynomial at the point for which
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
        (self.w.len() * ark_ff::to_bytes![E::G1Affine::zero()].unwrap().len()) / 2 + hiding_size
    }
}

impl<E: PairingEngine> ToBytes for Proof<E> {
    #[inline]
    fn write<W: Write>(&self, mut writer: W) -> ark_std::io::Result<()> {
        self.w
            .iter()
            .map(|e| e.write(&mut writer))
            .collect::<Result<_, _>>()?;
        self.random_v
            .as_ref()
            .unwrap_or(&E::Fr::zero())
            .write(&mut writer)
    }
}
