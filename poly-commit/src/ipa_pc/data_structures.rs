use crate::*;
use crate::{PCCommitterKey, PCVerifierKey, Vec};
use ark_crypto_primitives::sponge::Absorb;
use ark_ec::AffineRepr;
use ark_ff::{Field, UniformRand, Zero};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::rand::RngCore;
use ark_std::vec;

/// `UniversalParams` are the universal parameters for the inner product arg scheme.
#[derive(Derivative, CanonicalSerialize, CanonicalDeserialize)]
#[derivative(Default(bound = ""), Clone(bound = ""), Debug(bound = ""))]
pub struct UniversalParams<G: AffineRepr> {
    /// The key used to commit to polynomials.
    pub comm_key: Vec<G>,

    /// Some group generator.
    pub h: G,

    /// Some group generator specifically used for hiding.
    pub s: G,
}

impl<G: AffineRepr> PCUniversalParams for UniversalParams<G> {
    fn max_degree(&self) -> usize {
        self.comm_key.len() - 1
    }
}

/// `CommitterKey` is used to commit to, and create evaluation proofs for, a given
/// polynomial.
#[derive(Derivative, CanonicalSerialize, CanonicalDeserialize)]
#[derivative(
    Default(bound = ""),
    Hash(bound = ""),
    Clone(bound = ""),
    Debug(bound = "")
)]
pub struct CommitterKey<G: AffineRepr> {
    /// The key used to commit to polynomials.
    pub comm_key: Vec<G>,

    /// A random group generator.
    pub h: G,

    /// A random group generator that is to be used to make
    /// a commitment hiding.
    pub s: G,

    /// The maximum degree supported by the parameters
    /// this key was derived from.
    pub max_degree: usize,
}

impl<G: AffineRepr> PCCommitterKey for CommitterKey<G> {
    fn max_degree(&self) -> usize {
        self.max_degree
    }
    fn supported_degree(&self) -> usize {
        self.comm_key.len() - 1
    }
}

/// `VerifierKey` is used to check evaluation proofs for a given commitment.
pub type VerifierKey<G> = CommitterKey<G>;

impl<G: AffineRepr> PCVerifierKey for VerifierKey<G> {
    fn max_degree(&self) -> usize {
        self.max_degree
    }

    fn supported_degree(&self) -> usize {
        self.comm_key.len() - 1
    }
}

/// Nothing to do to prepare this verifier key (for now).
pub type PreparedVerifierKey<G> = VerifierKey<G>;

impl<G: AffineRepr> PCPreparedVerifierKey<VerifierKey<G>> for PreparedVerifierKey<G> {
    /// prepare `PreparedVerifierKey` from `VerifierKey`
    fn prepare(vk: &VerifierKey<G>) -> Self {
        vk.clone()
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
pub struct Commitment<G: AffineRepr + Absorb> {
    /// A Pedersen commitment to the polynomial.
    pub comm: G,

    /// A Pedersen commitment to the shifted polynomial.
    /// This is `none` if the committed polynomial does not
    /// enforce a strict degree bound.
    pub shifted_comm: Option<G>,
}

impl<G: AffineRepr + Absorb> PCCommitment for Commitment<G> {
    #[inline]
    fn empty() -> Self {
        Commitment {
            comm: G::zero(),
            shifted_comm: None,
        }
    }

    fn has_degree_bound(&self) -> bool {
        false
    }
}

/// Nothing to do to prepare this commitment (for now).
pub type PreparedCommitment<E> = Commitment<E>;

impl<G: AffineRepr + Absorb> PCPreparedCommitment<Commitment<G>> for PreparedCommitment<G> {
    /// prepare `PreparedCommitment` from `Commitment`
    fn prepare(vk: &Commitment<G>) -> Self {
        vk.clone()
    }
}

/// `Randomness` hides the polynomial inside a commitment and is outputted by `InnerProductArg::commit`.
#[derive(Derivative, CanonicalSerialize, CanonicalDeserialize)]
#[derivative(
    Default(bound = ""),
    Hash(bound = ""),
    Clone(bound = ""),
    Debug(bound = ""),
    PartialEq(bound = ""),
    Eq(bound = "")
)]
pub struct Randomness<G: AffineRepr> {
    /// Randomness is some scalar field element.
    pub rand: G::ScalarField,

    /// Randomness applied to the shifted commitment is some scalar field element.
    pub shifted_rand: Option<G::ScalarField>,
}

impl<G: AffineRepr> PCCommitmentState for Randomness<G> {
    type Randomness = Self;
    fn empty() -> Self {
        Self {
            rand: G::ScalarField::zero(),
            shifted_rand: None,
        }
    }

    fn rand<R: RngCore>(_: usize, has_degree_bound: bool, _: Option<usize>, rng: &mut R) -> Self {
        let rand = G::ScalarField::rand(rng);
        let shifted_rand = if has_degree_bound {
            Some(G::ScalarField::rand(rng))
        } else {
            None
        };

        Self { rand, shifted_rand }
    }
}

/// `Proof` is an evaluation proof that is output by `InnerProductArg::open`.
#[derive(Derivative, CanonicalSerialize, CanonicalDeserialize)]
#[derivative(
    Default(bound = ""),
    Hash(bound = ""),
    Clone(bound = ""),
    Debug(bound = "")
)]
pub struct Proof<G: AffineRepr> {
    /// Vector of left elements for each of the log_d iterations in `open`
    pub l_vec: Vec<G>,

    /// Vector of right elements for each of the log_d iterations within `open`
    pub r_vec: Vec<G>,

    /// Committer key from the last iteration within `open`
    pub final_comm_key: G,

    /// Coefficient from the last iteration within withinopen`
    pub c: G::ScalarField,

    /// Commitment to the blinding polynomial.
    pub hiding_comm: Option<G>,

    /// Linear combination of all the randomness used for commitments
    /// to the opened polynomials, along with the randomness used for the
    /// commitment to the hiding polynomial.
    pub rand: Option<G::ScalarField>,
}

/// `SuccinctCheckPolynomial` is a succinctly-representated polynomial
/// generated from the `log_d` random oracle challenges generated in `open`.
/// It has the special property that can be evaluated in `O(log_d)` time.
pub struct SuccinctCheckPolynomial<F: Field>(pub Vec<F>);

impl<F: Field> SuccinctCheckPolynomial<F> {
    /// Computes the coefficients of the underlying degree `d` polynomial.
    pub fn compute_coeffs(&self) -> Vec<F> {
        let challenges = &self.0;
        let log_d = challenges.len();

        let mut coeffs = vec![F::one(); 1 << log_d];
        for (i, challenge) in challenges.iter().enumerate() {
            let i = i + 1;
            let elem_degree = 1 << (log_d - i);
            for start in (elem_degree..coeffs.len()).step_by(elem_degree * 2) {
                for offset in 0..elem_degree {
                    coeffs[start + offset] *= challenge;
                }
            }
        }

        coeffs
    }

    /// Evaluate `self` at `point` in time `O(log_d)`.
    pub fn evaluate(&self, point: F) -> F {
        let challenges = &self.0;
        let log_d = challenges.len();

        let mut product = F::one();
        for (i, challenge) in challenges.iter().enumerate() {
            let i = i + 1;
            let elem_degree: u64 = (1 << (log_d - i)) as u64;
            let elem = point.pow([elem_degree]);
            product *= &(F::one() + &(elem * challenge));
        }

        product
    }
}
