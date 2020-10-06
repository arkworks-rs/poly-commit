use crate::{BTreeMap, MVPolynomial};
use crate::{PCCommitterKey, PCProof, PCRandomness, PCUniversalParams, PCVerifierKey};
use algebra_core::{PairingEngine, ToBytes, Zero};
use core::{
    marker::PhantomData,
    ops::{Add, AddAssign, Index},
};
use rand_core::RngCore;

/// `UniversalParams` are the universal parameters for the MarlinPST13 scheme.
#[derive(Derivative)]
#[derivative(Default(bound = ""), Clone(bound = ""), Debug(bound = ""))]
pub struct UniversalParams<E, P>
where
    E: PairingEngine,
    P: MVPolynomial<E::Fr>,
    P::Domain: Index<usize, Output = E::Fr>,
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

impl<E, P> PCUniversalParams for UniversalParams<E, P>
where
    E: PairingEngine,
    P: MVPolynomial<E::Fr>,
    P::Domain: Index<usize, Output = E::Fr>,
{
    fn max_degree(&self) -> usize {
        self.max_degree
    }
}

/// `CommitterKey` is used to commit to and create evaluation proofs for a given
/// polynomial.
#[derive(Derivative)]
#[derivative(
    Default(bound = ""),
    Hash(bound = ""),
    Clone(bound = ""),
    Debug(bound = "")
)]
pub struct CommitterKey<E, P>
where
    E: PairingEngine,
    P: MVPolynomial<E::Fr>,
    P::Domain: Index<usize, Output = E::Fr>,
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
    P::Domain: Index<usize, Output = E::Fr>,
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

impl<E: PairingEngine> PCVerifierKey for VerifierKey<E> {
    fn max_degree(&self) -> usize {
        self.max_degree
    }

    fn supported_degree(&self) -> usize {
        self.supported_degree
    }
}

/// `Randomness` hides the polynomial inside a commitment`.
#[derive(Derivative)]
#[derivative(
    Default(bound = ""),
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
    P::Domain: Index<usize, Output = E::Fr>,
{
    /// A multivariate polynomial where each monomial is univariate with random coefficient
    pub blinding_polynomial: P,
    _engine: PhantomData<E>,
}

impl<E, P> Randomness<E, P>
where
    E: PairingEngine,
    P: MVPolynomial<E::Fr>,
    P::Domain: Index<usize, Output = E::Fr>,
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
    P::Domain: Index<usize, Output = E::Fr>,
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
            blinding_polynomial: P::rand(hiding_poly_degree, num_vars, rng),
            _engine: PhantomData,
        }
    }
}

impl<'a, E: PairingEngine, P: MVPolynomial<E::Fr>> Add<&'a Randomness<E, P>> for Randomness<E, P>
where
    E: PairingEngine,
    P: MVPolynomial<E::Fr>,
    P::Domain: Index<usize, Output = E::Fr>,
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
    P::Domain: Index<usize, Output = E::Fr>,
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
    P::Domain: Index<usize, Output = E::Fr>,
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
    P::Domain: Index<usize, Output = E::Fr>,
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
            algebra_core::to_bytes![E::Fr::zero()].unwrap().len()
        } else {
            0
        };
        (self.w.len() * algebra_core::to_bytes![E::G1Affine::zero()].unwrap().len()) / 2
            + hiding_size
    }
}

impl<E: PairingEngine> ToBytes for Proof<E> {
    #[inline]
    fn write<W: algebra_core::io::Write>(&self, mut writer: W) -> algebra_core::io::Result<()> {
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
