use crate::*;
use crate::{PCCommitterKey, PCVerifierKey, Vec};
use algebra_core::{One, ToBytes, PrimeField, AffineCurve, ProjectiveCurve};
use std::ops::{AddAssign, Add, Sub};

#[derive(Derivative)]
#[derivative(Default(bound = ""), Clone(bound = ""), Debug(bound = ""))]
pub struct UniversalParams<G: AffineCurve> {
    pub comm_key: Vec<G>,
    pub h: G
}

impl<G: AffineCurve> PCUniversalParams for UniversalParams<G> {
    fn max_degree(&self) -> usize {
        self.comm_key.len() - 1
    }
}

#[derive(Derivative)]
#[derivative(
    Default(bound = ""),
    Hash(bound = ""),
    Clone(bound = ""),
    Debug(bound = "")
)]
pub struct CommitterKey<G: AffineCurve> {
    pub comm_key: Vec<G>,
    pub h: G,
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

#[derive(Derivative)]
#[derivative(
Default(bound = ""),
Hash(bound = ""),
Clone(bound = ""),
Debug(bound = "")
)]
pub struct VerifierKey<G: AffineCurve> {
    pub comm_key: Vec<G>,
    pub h: G,
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
    pub comm: G,
    pub shifted_comm: Option<G>
}

impl<G: AffineCurve> PCCommitment for Commitment<G> {
    #[inline]
    fn empty() -> Self {
        Commitment{
            comm: G::zero(),
            shifted_comm: None
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

#[derive(Derivative)]
#[derivative(
Default(bound = ""),
Hash(bound = ""),
Clone(bound = ""),
Debug(bound = ""),
PartialEq(bound = ""),
Eq(bound = "")
)]
// TODO: Placeholder randomness
pub struct Randomness();
impl PCRandomness for Randomness {
    fn empty() -> Self {
        Randomness ()
    }
    fn rand <R: RngCore> (_: usize, _: bool, _: &mut R) -> Self {
        Randomness()
    }
}

#[derive(Derivative)]
#[derivative(
Default(bound = ""),
Hash(bound = ""),
Clone(bound = ""),
Debug(bound = "")
)]
pub struct Proof<G: AffineCurve> {
    pub l_vec: Vec<G>,
    pub r_vec: Vec<G>,
    pub final_comm_key: G,
    pub c: G::ScalarField
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
        self.c.write(&mut writer)
    }
}

pub struct SuccinctCheckPolynomial<F: Field> (pub Vec<F>);

impl<F: Field> SuccinctCheckPolynomial<F> {
    pub fn compute_coeffs (&self) -> Vec<F> {
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

    pub fn evaluate (&self, point: F) -> F {
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
