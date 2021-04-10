use crate::{
    PCCommitment, PCCommitterKey, PCPreparedCommitment, PCPreparedVerifierKey, PCRandomness,
    PCVerifierKey, UVPolynomial, Vec,
};
use ark_ec::{PairingEngine, ProjectiveCurve};
use ark_ff::{Field, PrimeField, ToBytes, ToConstraintField};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, SerializationError};
use ark_std::io::{Read, Write};
use ark_std::ops::{Add, AddAssign};
use ark_std::rand::RngCore;

use crate::kzg10;
/// `UniversalParams` are the universal parameters for the KZG10 scheme.
pub type UniversalParams<E> = kzg10::UniversalParams<E>;

/// `CommitterKey` is used to commit to and create evaluation proofs for a given
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
    pub powers: Vec<E::G1Affine>,

    /// The key used to commit to shifted polynomials.
    /// This is `None` if `self` does not support enforcing any degree bounds.
    pub shifted_powers: Option<Vec<E::G1Affine>>,

    /// The key used to commit to hiding polynomials.
    pub powers_of_gamma_g: Vec<E::G1Affine>,

    /// The degree bounds that are supported by `self`.
    /// In ascending order from smallest to largest.
    /// This is `None` if `self` does not support enforcing any degree bounds.
    pub enforced_degree_bounds: Option<Vec<usize>>,
    /// The maximum degree supported by the `UniversalParams` `self` was derived
    /// from.
    pub max_degree: usize,
}

impl<E: PairingEngine> CommitterKey<E> {
    /// Obtain powers for the underlying KZG10 construction
    pub fn powers<'a>(&'a self) -> kzg10::Powers<'a, E> {
        kzg10::Powers {
            powers_of_g: self.powers.as_slice().into(),
            powers_of_gamma_g: self.powers_of_gamma_g.as_slice().into(),
        }
    }

    /// Obtain powers for committing to shifted polynomials.
    pub fn shifted_powers<'a>(
        &'a self,
        degree_bound: impl Into<Option<usize>>,
    ) -> Option<kzg10::Powers<'a, E>> {
        self.shifted_powers.as_ref().map(|shifted_powers| {
            let powers_range = if let Some(degree_bound) = degree_bound.into() {
                assert!(self
                    .enforced_degree_bounds
                    .as_ref()
                    .unwrap()
                    .contains(&degree_bound));
                let max_bound = self
                    .enforced_degree_bounds
                    .as_ref()
                    .unwrap()
                    .last()
                    .unwrap();
                (max_bound - degree_bound)..
            } else {
                0..
            };
            let ck = kzg10::Powers {
                powers_of_g: (&shifted_powers[powers_range]).into(),
                powers_of_gamma_g: self.powers_of_gamma_g.as_slice().into(),
            };
            ck
        })
    }
}

impl<E: PairingEngine> PCCommitterKey for CommitterKey<E> {
    fn max_degree(&self) -> usize {
        self.max_degree
    }

    fn supported_degree(&self) -> usize {
        self.powers.len() - 1
    }
}

/// `VerifierKey` is used to check evaluation proofs for a given commitment.
#[derive(Derivative, CanonicalSerialize, CanonicalDeserialize)]
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
            v.binary_search_by(|(d, _)| d.cmp(&bound))
                .ok()
                .map(|i| v[i].1)
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

impl<E: PairingEngine> ToBytes for VerifierKey<E> {
    #[inline]
    fn write<W: Write>(&self, mut writer: W) -> ark_std::io::Result<()> {
        self.vk.write(&mut writer)?;
        if let Some(degree_bounds_and_shift_powers) = &self.degree_bounds_and_shift_powers {
            writer.write_all(&degree_bounds_and_shift_powers.len().to_le_bytes())?;
            for (degree_bound, shift_power) in degree_bounds_and_shift_powers {
                writer.write_all(&degree_bound.to_le_bytes())?;
                shift_power.write(&mut writer)?;
            }
        }
        writer.write_all(&self.supported_degree.to_le_bytes())?;
        writer.write_all(&self.max_degree.to_le_bytes())
    }
}

impl<E: PairingEngine> ToConstraintField<<E::Fq as Field>::BasePrimeField> for VerifierKey<E>
where
    E::G1Affine: ToConstraintField<<E::Fq as Field>::BasePrimeField>,
    E::G2Affine: ToConstraintField<<E::Fq as Field>::BasePrimeField>,
{
    fn to_field_elements(&self) -> Option<Vec<<E::Fq as Field>::BasePrimeField>> {
        let mut res = Vec::new();
        res.extend_from_slice(&self.vk.to_field_elements().unwrap());

        if let Some(degree_bounds_and_shift_powers) = &self.degree_bounds_and_shift_powers {
            for (d, shift_power) in degree_bounds_and_shift_powers.iter() {
                let d_elem: <E::Fq as Field>::BasePrimeField = (*d as u64).into();

                res.push(d_elem);
                res.extend_from_slice(&shift_power.to_field_elements().unwrap());
            }
        }

        Some(res)
    }
}

/// `PreparedVerifierKey` is used to check evaluation proofs for a given commitment.
#[derive(Derivative)]
#[derivative(Clone(bound = ""), Debug(bound = ""))]
pub struct PreparedVerifierKey<E: PairingEngine> {
    /// The verification key for the underlying KZG10 scheme.
    pub prepared_vk: kzg10::PreparedVerifierKey<E>,
    /// Information required to enforce degree bounds. Each pair
    /// is of the form `(degree_bound, shifting_advice)`.
    /// This is `None` if `self` does not support enforcing any degree bounds.
    pub prepared_degree_bounds_and_shift_powers: Option<Vec<(usize, Vec<E::G1Affine>)>>,
    /// The maximum degree supported by the `UniversalParams` `self` was derived
    /// from.
    pub max_degree: usize,
    /// The maximum degree supported by the trimmed parameters that `self` is
    /// a part of.
    pub supported_degree: usize,
}

impl<E: PairingEngine> PCPreparedVerifierKey<VerifierKey<E>> for PreparedVerifierKey<E> {
    /// prepare `PreparedVerifierKey` from `VerifierKey`
    fn prepare(vk: &VerifierKey<E>) -> Self {
        let prepared_vk = kzg10::PreparedVerifierKey::<E>::prepare(&vk.vk);

        let supported_bits = E::Fr::size_in_bits();

        let prepared_degree_bounds_and_shift_powers: Option<Vec<(usize, Vec<E::G1Affine>)>> =
            if vk.degree_bounds_and_shift_powers.is_some() {
                let mut res = Vec::<(usize, Vec<E::G1Affine>)>::new();

                let degree_bounds_and_shift_powers =
                    vk.degree_bounds_and_shift_powers.as_ref().unwrap();

                for (d, shift_power) in degree_bounds_and_shift_powers {
                    let mut prepared_shift_power = Vec::<E::G1Affine>::new();

                    let mut cur = E::G1Projective::from(shift_power.clone());
                    for _ in 0..supported_bits {
                        prepared_shift_power.push(cur.clone().into());
                        cur.double_in_place();
                    }

                    res.push((d.clone(), prepared_shift_power));
                }

                Some(res)
            } else {
                None
            };

        Self {
            prepared_vk,
            prepared_degree_bounds_and_shift_powers,
            max_degree: vk.max_degree,
            supported_degree: vk.supported_degree,
        }
    }
}

/// Commitment to a polynomial that optionally enforces a degree bound.
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
pub struct Commitment<E: PairingEngine> {
    /// A KZG10 commitment to the polynomial.
    pub comm: kzg10::Commitment<E>,

    /// A KZG10 commitment to the shifted polynomial.
    /// This is `none` if the committed polynomial does not
    /// enforce a strict degree bound.
    pub shifted_comm: Option<kzg10::Commitment<E>>,
}

impl<E: PairingEngine> ToBytes for Commitment<E> {
    #[inline]
    fn write<W: Write>(&self, mut writer: W) -> ark_std::io::Result<()> {
        self.comm.write(&mut writer)?;
        let shifted_exists = self.shifted_comm.is_some();
        shifted_exists.write(&mut writer)?;
        self.shifted_comm
            .as_ref()
            .unwrap_or(&kzg10::Commitment::empty())
            .write(&mut writer)
    }
}

impl<E: PairingEngine> ToConstraintField<<E::Fq as Field>::BasePrimeField> for Commitment<E>
where
    E::G1Affine: ToConstraintField<<E::Fq as Field>::BasePrimeField>,
{
    fn to_field_elements(&self) -> Option<Vec<<E::Fq as Field>::BasePrimeField>> {
        let mut res = Vec::new();
        res.extend_from_slice(&self.comm.to_field_elements().unwrap());

        if let Some(shifted_comm) = &self.shifted_comm {
            res.extend_from_slice(&shifted_comm.to_field_elements().unwrap());
        }

        Some(res)
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

/// Prepared commitment to a polynomial that optionally enforces a degree bound.
#[derive(Derivative)]
#[derivative(
    Hash(bound = ""),
    Clone(bound = ""),
    Debug(bound = ""),
    PartialEq(bound = ""),
    Eq(bound = "")
)]
pub struct PreparedCommitment<E: PairingEngine> {
    pub(crate) prepared_comm: kzg10::PreparedCommitment<E>,
    pub(crate) shifted_comm: Option<kzg10::Commitment<E>>,
}

impl<E: PairingEngine> PCPreparedCommitment<Commitment<E>> for PreparedCommitment<E> {
    /// Prepare commitment to a polynomial that optionally enforces a degree bound.
    fn prepare(comm: &Commitment<E>) -> Self {
        let prepared_comm = kzg10::PreparedCommitment::<E>::prepare(&comm.comm);

        let shifted_comm = comm.shifted_comm.clone();

        Self {
            prepared_comm,
            shifted_comm,
        }
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
    /// Commitment randomness for a KZG10 commitment.
    pub rand: kzg10::Randomness<F, P>,
    /// Commitment randomness for a KZG10 commitment to the shifted polynomial.
    /// This is `None` if the committed polynomial does not enforce a strict
    /// degree bound.
    pub shifted_rand: Option<kzg10::Randomness<F, P>>,
}

impl<'a, F: PrimeField, P: UVPolynomial<F>> Add<&'a Self> for Randomness<F, P> {
    type Output = Self;

    fn add(mut self, other: &'a Self) -> Self {
        self += other;
        self
    }
}

impl<'a, F: PrimeField, P: UVPolynomial<F>> AddAssign<&'a Self> for Randomness<F, P> {
    #[inline]
    fn add_assign(&mut self, other: &'a Self) {
        self.rand += &other.rand;
        if let Some(r1) = &mut self.shifted_rand {
            *r1 += other
                .shifted_rand
                .as_ref()
                .unwrap_or(&kzg10::Randomness::empty());
        } else {
            self.shifted_rand = other.shifted_rand.as_ref().map(|r| r.clone());
        }
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

impl<'a, F: PrimeField, P: UVPolynomial<F>> AddAssign<(F, &'a Randomness<F, P>)>
    for Randomness<F, P>
{
    #[inline]
    fn add_assign(&mut self, (f, other): (F, &'a Randomness<F, P>)) {
        self.rand += (f, &other.rand);
        let empty = kzg10::Randomness::empty();
        if let Some(r1) = &mut self.shifted_rand {
            *r1 += (f, other.shifted_rand.as_ref().unwrap_or(&empty));
        } else {
            self.shifted_rand = other.shifted_rand.as_ref().map(|r| empty + (f, r));
        }
    }
}

impl<F: PrimeField, P: UVPolynomial<F>> PCRandomness for Randomness<F, P> {
    fn empty() -> Self {
        Self {
            rand: kzg10::Randomness::empty(),
            shifted_rand: None,
        }
    }

    fn rand<R: RngCore>(
        hiding_bound: usize,
        has_degree_bound: bool,
        _: Option<usize>,
        rng: &mut R,
    ) -> Self {
        let shifted_rand = if has_degree_bound {
            Some(kzg10::Randomness::rand(hiding_bound, false, None, rng))
        } else {
            None
        };
        Self {
            rand: kzg10::Randomness::rand(hiding_bound, false, None, rng),
            shifted_rand,
        }
    }
}
