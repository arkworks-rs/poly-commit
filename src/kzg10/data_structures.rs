use algebra::PairingEngine;

/// `UniversalParams` are the universal parameters for the KZG10 scheme.
#[derive(Derivative)]
#[derivative(
    Default(bound = ""),
    Clone(bound = ""),
    Debug(bound = "")
)]
pub struct UniversalParams<E: PairingEngine> {
    /// Group elements of the form `{ \beta^i G }`, where `i` ranges from 0 to `degree`.
    pub powers_of_g: Vec<E::G1Affine>,
    /// Group elements of the form `{ \beta^i \gamma G }`, where `i` ranges from 0 to `degree`.
    pub powers_of_gamma_g: Vec<E::G1Affine>,
    /// The generator of G2.
    pub h: E::G2Affine,
    /// \beta times the above generator of G2.
    pub beta_h: E::G2Affine,
    /// The generator of G2, prepared for use in pairings.
    pub prepared_h: <E::G2Affine as PairingCurve>::Prepared,
    /// \beta times the above generator of G2, prepared for use in pairings.
    pub prepared_beta_h: <E::G2Affine as PairingCurve>::Prepared,
    /// Maximum degree of polynomials supported by these parameters
    pub max_degree: usize,
}

impl<E: PairingEngine> PCUniversalParams for UniversalParams<E> {
    fn max_degree(&self) -> usize {
        self.max_degree
    }
}

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
    /// Group elements of the form `{ \beta^i G }`, where `i` ranges from 0 to `degree`.
    pub powers_of_g: Vec<E::G1Affine>,
    /// Group elements of the form `{ \beta^i \gamma G }`, where `i` ranges from 0 to `degree`.
    pub powers_of_gamma_g: Vec<E::G1Affine>,
    /// Maximum degree of polynomials supported by these parameters
    pub max_degree: usize,
}

impl<E: PairingEngine> PCCommitterKey for CommitterKey<E> {
    fn max_degree(&self) -> usize {
        self.max_degree
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
    pub prepared_h: <E::G2Affine as PairingCurve>::Prepared,
    /// \beta times the above generator of G2, prepared for use in pairings.
    pub prepared_beta_h: <E::G2Affine as PairingCurve>::Prepared,
    /// Maximum degree of polynomials supported by these parameters
    pub max_degree: usize,
}

impl<E: PairingEngine> PCVerifierKey for VerifierKey<E> {
    fn max_degree(&self) -> usize {
        self.max_degree
    }
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
    fn write<W: std::io::Write>(&self, writer: W) -> std::io::Result<()> {
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
        to_bytes![E::G1Affine::zero()].unwrap().len() / 2
    }
}

impl<'a, E: PairingEngine> AddAssign<(E::Fr, &'a Commitment<E>)> for Commitment<E> {
    #[inline]
    fn add_assign(&mut self, (f, other): (E::Fr, &'a Commitment<E>)) {
        let other = other.0.into_projective().mul(&f);
        self.0 = (self.0.into_projective() + &other).into();
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
    /// For KZG10, the commitment randomness is a random polynomial.
    blinding_polynomial: Polynomial<E::Fr>,
}

impl<E: PairingEngine> Randomness<E> {
    #[inline]
    fn is_hiding(&self) -> bool {
        !self.blinding_polynomial.is_zero()
    }
}

impl<E: PairingEngine> PCRandomness for Randomness<E> {
    fn empty() -> Self {
        Self {
            blinding_polynomial: Polynomial::zero(),
        }
    }

    fn rand<R: RngCore>(d: usize, rng: &mut R) -> Self {
        let mut randomness = Randomness::empty();
        randomness.blinding_polynomial = Polynomial::rand(d + 1, rng);
        randomness
    }
}

impl<'a, E: PairingEngine> AddAssign<(E::Fr, &'a Randomness<E>)> for Randomness<E> {
    #[inline]
    fn add_assign(&mut self, (f, other): (E::Fr, &'a Randomness<E>)) {
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
    pub random_v: E::Fr,
}
