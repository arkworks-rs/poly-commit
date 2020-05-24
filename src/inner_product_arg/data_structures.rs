use crate::*;
use crate::{PCCommitterKey, PCVerifierKey, Vec};
use algebra_core::{AffineCurve, Field, ToBytes, UniformRand, Zero};
use rand_core::RngCore;

/// `UniversalParams` are the universal parameters for the inner product arg scheme.
#[derive(Derivative)]
#[derivative(Default(bound = ""), Clone(bound = ""), Debug(bound = ""))]
pub struct UniversalParams<G: AffineCurve> {
    /// The key used to commit to polynomials.
    pub comm_key: Vec<G>,

    /// Some group generator.
    pub h: G,

    /// Some group generator specifically used for hiding.
    pub s: G,
}

impl<G: AffineCurve> PCUniversalParams for UniversalParams<G> {
    fn max_degree(&self) -> usize {
        self.comm_key.len() - 1
    }
}

/// `CommitterKey` is used to commit to, and create evaluation proofs for, a given
/// polynomial.
#[derive(Derivative)]
#[derivative(
    Default(bound = ""),
    Hash(bound = ""),
    Clone(bound = ""),
    Debug(bound = "")
)]
pub struct CommitterKey<G: AffineCurve> {
    /// The key used to commit to polynomials.
    pub comm_key: Vec<G>,

    /// Some group generator.
    pub h: G,

    /// Some group generator specifically used for hiding.
    pub s: G,

    /// Max degree supported by the universal parameters.
    pub max_degree: usize,
}

impl<G: AffineCurve> PCCommitterKey for CommitterKey<G> {
    fn max_degree(&self) -> usize {
        self.max_degree
    }
    fn supported_degree(&self) -> usize {
        self.comm_key.len() - 1
    }
}

/// `VerifierKey` is used to check evaluation proofs for a given commitment.
#[derive(Derivative)]
#[derivative(
    Default(bound = ""),
    Hash(bound = ""),
    Clone(bound = ""),
    Debug(bound = "")
)]
pub struct VerifierKey<G: AffineCurve> {
    /// The key used to commit to polynomials.
    pub comm_key: Vec<G>,

    /// Some group generator.
    pub h: G,

    /// Some group generator specifically used for hiding.
    pub s: G,

    /// Max degree supported by the universal parameters.
    pub max_degree: usize,
}

impl<G: AffineCurve> PCVerifierKey for VerifierKey<G> {
    fn max_degree(&self) -> usize {
        self.max_degree
    }

    fn supported_degree(&self) -> usize {
        self.comm_key.len() - 1
    }
}

/// Commitment to a polynomial that optionally enforces a degree bound.
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
pub struct Commitment<G: AffineCurve> {
    /// Commitment is some group element.
    pub comm: G,

    /// Commitment to the shifted polynomial is some group element.
    /// Is `none` if polynomial has no degree bound.
    pub shifted_comm: Option<G>,
}

impl<G: AffineCurve> PCCommitment for Commitment<G> {
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

    fn size_in_bytes(&self) -> usize {
        algebra_core::to_bytes![G::zero()].unwrap().len() / 2
    }
}

impl<G: AffineCurve> ToBytes for Commitment<G> {
    #[inline]
    fn write<W: algebra_core::io::Write>(&self, mut writer: W) -> algebra_core::io::Result<()> {
        self.comm.write(&mut writer)?;
        let shifted_exists = self.shifted_comm.is_some();
        shifted_exists.write(&mut writer)?;
        self.shifted_comm
            .as_ref()
            .unwrap_or(&G::zero())
            .write(&mut writer)
    }
}

/// `Randomness` hides the polynomial inside a commitment and is outputted by `InnerProductArg::commit`.
#[derive(Derivative)]
#[derivative(
    Default(bound = ""),
    Hash(bound = ""),
    Clone(bound = ""),
    Debug(bound = ""),
    PartialEq(bound = ""),
    Eq(bound = "")
)]
pub struct Randomness<G: AffineCurve> {
    /// Randomness is some scalar field element.
    pub rand: G::ScalarField,

    /// Randomness applied to the shifted commitment is some scalar field element.
    pub shifted_rand: Option<G::ScalarField>,
}

impl<G: AffineCurve> PCRandomness for Randomness<G> {
    fn empty() -> Self {
        Self {
            rand: G::ScalarField::zero(),
            shifted_rand: None,
        }
    }

    fn rand<R: RngCore>(_hiding_bound: usize, has_degree_bound: bool, rng: &mut R) -> Self {
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
#[derive(Derivative)]
#[derivative(
    Default(bound = ""),
    Hash(bound = ""),
    Clone(bound = ""),
    Debug(bound = "")
)]
pub struct Proof<G: AffineCurve> {
    /// Vector of left elements from each of the log_d sub-procedures in `open`
    pub l_vec: Vec<G>,

    /// Vector of right elements from each of the log_d sub-procedures of `open`
    pub r_vec: Vec<G>,

    /// Last committer key from the last of the log_d sub-procedures of `open`
    pub final_comm_key: G,

    /// Last coefficient from the last of the log_d sub-procedures of `open`
    pub c: G::ScalarField,

    /// Commitment to the hiding polynomial
    pub hiding_comm: Option<G>,

    /// Sum of all randomness in commitments and hiding polynomial
    pub rand: Option<G::ScalarField>,
}

impl<G: AffineCurve> PCProof for Proof<G> {
    fn size_in_bytes(&self) -> usize {
        algebra_core::to_bytes![self].unwrap().len()
    }
}

impl<G: AffineCurve> ToBytes for Proof<G> {
    #[inline]
    fn write<W: algebra_core::io::Write>(&self, mut writer: W) -> algebra_core::io::Result<()> {
        self.l_vec.write(&mut writer)?;
        self.r_vec.write(&mut writer)?;
        self.final_comm_key.write(&mut writer)?;
        self.c.write(&mut writer)?;
        self.hiding_comm
            .as_ref()
            .unwrap_or(&G::zero())
            .write(&mut writer)?;
        self.rand
            .as_ref()
            .unwrap_or(&G::ScalarField::zero())
            .write(&mut writer)
    }
}

/// `SuccinctCheckPolynomial` is a special check polynomial generated
/// from a size `log_d` vector of scalar field elements.
/// Can be evaluated in `O(log_d)` but can optionally be represented as a degree `d` polynomial.
pub struct SuccinctCheckPolynomial<F: Field>(pub Vec<F>);

impl<F: Field> SuccinctCheckPolynomial<F> {
    /// Computes the underlying degree `d` polynomial coefficients
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

    /// Evaluates the underlying degree `d` polynomial coefficients in `O(log_d)`
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
