use crate::*;
use algebra_core::{AffineCurve, PairingEngine, PrimeField, ProjectiveCurve, ToBytes, Zero};
use core::marker::PhantomData;
use core::ops::{Add, AddAssign};

/// `UniversalParams` are the universal parameters for the KZG10 scheme.
#[derive(Derivative)]
#[derivative(Default(bound = ""), Clone(bound = ""), Debug(bound = ""))]
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
    pub prepared_neg_powers_of_h: BTreeMap<usize, E::G2Prepared>,
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

/// `Commitment` commits to a polynomial. It is output by `KZG10::commit`.
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
pub struct Commitment<E: PairingEngine>(
    /// The commitment is a group element.
    pub E::G1Affine,
);

impl<E: PairingEngine> ToBytes for Commitment<E> {
    #[inline]
    fn write<W: algebra_core::io::Write>(&self, writer: W) -> algebra_core::io::Result<()> {
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
        algebra_core::to_bytes![E::G1Affine::zero()].unwrap().len() / 2
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
pub struct Randomness<E: PairingEngine, P: UVPolynomial<E::Fr>> {
    /// For KZG10, the commitment randomness is a random polynomial.
    pub blinding_polynomial: P,
    _engine: PhantomData<E>,
}

impl<E: PairingEngine, P: UVPolynomial<E::Fr>> Randomness<E, P> {
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

impl<E: PairingEngine, P: UVPolynomial<E::Fr>> PCRandomness for Randomness<E, P> {
    fn empty() -> Self {
        Self {
            blinding_polynomial: P::zero(),
            _engine: PhantomData,
        }
    }

    fn rand<R: RngCore>(hiding_bound: usize, _: bool, _: Option<usize>, rng: &mut R) -> Self {
        let mut randomness = Randomness::empty();
        let hiding_poly_degree = Self::calculate_hiding_polynomial_degree(hiding_bound);
        randomness.blinding_polynomial = P::rand(hiding_poly_degree, None, rng);
        randomness
    }
}

impl<'a, E: PairingEngine, P: UVPolynomial<E::Fr>> Add<&'a Randomness<E, P>> for Randomness<E, P> {
    type Output = Self;

    #[inline]
    fn add(mut self, other: &'a Self) -> Self {
        self.blinding_polynomial += &other.blinding_polynomial;
        self
    }
}

impl<'a, E: PairingEngine, P: UVPolynomial<E::Fr>> Add<(E::Fr, &'a Randomness<E, P>)>
    for Randomness<E, P>
{
    type Output = Self;

    #[inline]
    fn add(mut self, other: (E::Fr, &'a Randomness<E, P>)) -> Self {
        self += other;
        self
    }
}

impl<'a, E: PairingEngine, P: UVPolynomial<E::Fr>> AddAssign<&'a Randomness<E, P>>
    for Randomness<E, P>
{
    #[inline]
    fn add_assign(&mut self, other: &'a Self) {
        self.blinding_polynomial += &other.blinding_polynomial;
    }
}

impl<'a, E: PairingEngine, P: UVPolynomial<E::Fr>> AddAssign<(E::Fr, &'a Randomness<E, P>)>
    for Randomness<E, P>
{
    #[inline]
    fn add_assign(&mut self, (f, other): (E::Fr, &'a Randomness<E, P>)) {
        self.blinding_polynomial += (f, &other.blinding_polynomial);
    }
}

/// `Proof` is an evaluation proof that is output by `KZG10::open`.
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
            algebra_core::to_bytes![E::Fr::zero()].unwrap().len()
        } else {
            0
        };
        algebra_core::to_bytes![E::G1Affine::zero()].unwrap().len() / 2 + hiding_size
    }
}

impl<E: PairingEngine> ToBytes for Proof<E> {
    #[inline]
    fn write<W: algebra_core::io::Write>(&self, mut writer: W) -> algebra_core::io::Result<()> {
        self.w.write(&mut writer)?;
        self.random_v
            .as_ref()
            .unwrap_or(&E::Fr::zero())
            .write(&mut writer)
    }
}
