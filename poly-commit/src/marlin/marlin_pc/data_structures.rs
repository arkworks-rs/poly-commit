use crate::{
    DenseUVPolynomial, PCCommitment, PCCommitmentState, PCCommitterKey, PCPreparedCommitment,
    PCPreparedVerifierKey, PCVerifierKey, Vec,
};
use ark_crypto_primitives::sponge::Absorb;
use ark_ec::pairing::Pairing;
use ark_ec::AdditiveGroup;
use ark_ff::{Field, PrimeField, ToConstraintField};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
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
pub struct CommitterKey<E: Pairing> {
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

impl<E: Pairing> CommitterKey<E> {
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

impl<E: Pairing> PCCommitterKey for CommitterKey<E> {
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
pub struct VerifierKey<E: Pairing> {
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

impl<E: Pairing> VerifierKey<E> {
    /// Find the appropriate shift for the degree bound.
    pub fn get_shift_power(&self, bound: usize) -> Option<E::G1Affine> {
        self.degree_bounds_and_shift_powers.as_ref().and_then(|v| {
            v.binary_search_by(|(d, _)| d.cmp(&bound))
                .ok()
                .map(|i| v[i].1)
        })
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

impl<E: Pairing> ToConstraintField<<E::TargetField as Field>::BasePrimeField> for VerifierKey<E>
where
    E::G1Affine: ToConstraintField<<E::TargetField as Field>::BasePrimeField>,
    E::G2Affine: ToConstraintField<<E::TargetField as Field>::BasePrimeField>,
{
    fn to_field_elements(&self) -> Option<Vec<<E::TargetField as Field>::BasePrimeField>> {
        let mut res = Vec::new();
        res.extend_from_slice(&self.vk.to_field_elements().unwrap());

        if let Some(degree_bounds_and_shift_powers) = &self.degree_bounds_and_shift_powers {
            for (d, shift_power) in degree_bounds_and_shift_powers.iter() {
                let d_elem: <E::TargetField as Field>::BasePrimeField = (*d as u64).into();

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
pub struct PreparedVerifierKey<E: Pairing> {
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

impl<E: Pairing> PCPreparedVerifierKey<VerifierKey<E>> for PreparedVerifierKey<E> {
    /// prepare `PreparedVerifierKey` from `VerifierKey`
    fn prepare(vk: &VerifierKey<E>) -> Self {
        let prepared_vk = kzg10::PreparedVerifierKey::<E>::prepare(&vk.vk);

        let supported_bits = E::ScalarField::MODULUS_BIT_SIZE as usize;

        let prepared_degree_bounds_and_shift_powers: Option<Vec<(usize, Vec<E::G1Affine>)>> =
            if vk.degree_bounds_and_shift_powers.is_some() {
                let mut res = Vec::<(usize, Vec<E::G1Affine>)>::new();

                let degree_bounds_and_shift_powers =
                    vk.degree_bounds_and_shift_powers.as_ref().unwrap();

                for (d, shift_power) in degree_bounds_and_shift_powers {
                    let mut prepared_shift_power = Vec::<E::G1Affine>::new();

                    let mut cur = E::G1::from(shift_power.clone());
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
pub struct Commitment<E>
where
    E: Pairing,
    E::G1Affine: Absorb,
{
    /// A KZG10 commitment to the polynomial.
    pub comm: kzg10::Commitment<E>,

    /// A KZG10 commitment to the shifted polynomial.
    /// This is `none` if the committed polynomial does not
    /// enforce a strict degree bound.
    pub shifted_comm: Option<kzg10::Commitment<E>>,
}

impl<E: Pairing> ToConstraintField<<E::TargetField as Field>::BasePrimeField> for Commitment<E>
where
    E::G1Affine: ToConstraintField<<E::TargetField as Field>::BasePrimeField> + Absorb,
{
    fn to_field_elements(&self) -> Option<Vec<<E::TargetField as Field>::BasePrimeField>> {
        let mut res = Vec::new();
        res.extend_from_slice(&self.comm.to_field_elements().unwrap());

        if let Some(shifted_comm) = &self.shifted_comm {
            res.extend_from_slice(&shifted_comm.to_field_elements().unwrap());
        }

        Some(res)
    }
}

impl<E> PCCommitment for Commitment<E>
where
    E: Pairing,
    E::G1Affine: Absorb,
{
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
pub struct PreparedCommitment<E>
where
    E: Pairing,
    E::G1Affine: Absorb,
{
    pub(crate) prepared_comm: kzg10::PreparedCommitment<E>,
    pub(crate) shifted_comm: Option<kzg10::Commitment<E>>,
}

impl<E> PCPreparedCommitment<Commitment<E>> for PreparedCommitment<E>
where
    E: Pairing,
    E::G1Affine: Absorb,
{
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
pub struct Randomness<F: PrimeField, P: DenseUVPolynomial<F>> {
    /// Commitment randomness for a KZG10 commitment.
    pub rand: kzg10::Randomness<F, P>,
    /// Commitment randomness for a KZG10 commitment to the shifted polynomial.
    /// This is `None` if the committed polynomial does not enforce a strict
    /// degree bound.
    pub shifted_rand: Option<kzg10::Randomness<F, P>>,
}

impl<'a, F: PrimeField, P: DenseUVPolynomial<F>> Add<&'a Self> for Randomness<F, P> {
    type Output = Self;

    fn add(mut self, other: &'a Self) -> Self {
        self += other;
        self
    }
}

impl<'a, F: PrimeField, P: DenseUVPolynomial<F>> AddAssign<&'a Self> for Randomness<F, P> {
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

impl<'a, F: PrimeField, P: DenseUVPolynomial<F>> AddAssign<(F, &'a Randomness<F, P>)>
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

impl<F: PrimeField, P: DenseUVPolynomial<F>> PCCommitmentState for Randomness<F, P> {
    type Randomness = Self;
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
