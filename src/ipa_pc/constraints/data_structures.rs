use crate::ipa_pc::{Commitment, CommitterKey, Proof};
use crate::LabeledCommitment;
use ark_ec::AffineCurve;
use ark_ff::vec::Vec;
use ark_ff::PrimeField;
use ark_ff::{BitIteratorLE, Field};
use ark_nonnative_field::NonNativeFieldVar;
use ark_r1cs_std::alloc::{AllocVar, AllocationMode};
use ark_r1cs_std::bits::boolean::Boolean;
use ark_r1cs_std::fields::FieldVar;
use ark_r1cs_std::groups::CurveVar;
use ark_r1cs_std::ToBitsGadget;
use ark_relations::r1cs::{Namespace, SynthesisError};
use ark_std::borrow::Borrow;
use ark_std::marker::PhantomData;

pub type ConstraintF<G> = <<G as AffineCurve>::BaseField as Field>::BasePrimeField;
pub type NNFieldVar<G> = NonNativeFieldVar<<G as AffineCurve>::ScalarField, ConstraintF<G>>;

pub struct CommitterKeyVar<G, C>
where
    G: AffineCurve,
    C: CurveVar<G::Projective, ConstraintF<G>>,
{
    /// The key used to commit to polynomials.
    pub comm_key_var: Vec<C>,

    /// A random group generator.
    pub h_var: C,

    /// A random group generator that is to be used to make
    /// a commitment hiding.
    pub s_var: C,

    /// Phantom data
    pub _affine: PhantomData<G>,
}

impl<G, C> CommitterKeyVar<G, C>
where
    G: AffineCurve,
    C: CurveVar<G::Projective, ConstraintF<G>>,
{
    pub fn supported_degree(&self) -> usize {
        self.comm_key_var.len() - 1
    }
}

impl<G, C> AllocVar<CommitterKey<G>, ConstraintF<G>> for CommitterKeyVar<G, C>
where
    G: AffineCurve,
    C: CurveVar<G::Projective, ConstraintF<G>>,
{
    fn new_variable<T: Borrow<CommitterKey<G>>>(
        cs: impl Into<Namespace<ConstraintF<G>>>,
        f: impl FnOnce() -> Result<T, SynthesisError>,
        mode: AllocationMode,
    ) -> Result<Self, SynthesisError> {
        let ns = cs.into();
        let key = f()?;

        let comm_key_var = key
            .borrow()
            .comm_key
            .iter()
            .map(|elem| C::new_variable(ns.clone(), || Ok(elem.clone()), mode))
            .collect::<Result<Vec<C>, SynthesisError>>()?;

        let h_var = C::new_variable(ns.clone(), || Ok(key.borrow().h.clone()), mode)?;
        let s_var = C::new_variable(ns.clone(), || Ok(key.borrow().s.clone()), mode)?;

        Ok(Self {
            comm_key_var,
            h_var,
            s_var,
            _affine: PhantomData,
        })
    }
}

pub type VerifierKeyVar<G, C> = CommitterKeyVar<G, C>;

#[derive(Derivative)]
#[derivative(Clone(
    bound = "G: AffineCurve, C: CurveVar<G::Projective, <G::BaseField as Field>::BasePrimeField>"
))]
pub struct CommitmentVar<G, C>
where
    G: AffineCurve,
    C: CurveVar<G::Projective, ConstraintF<G>>,
{
    /// A Pedersen commitment to the polynomial.
    pub comm_var: C,

    /// The degree bound and a Pedersen commitment to the shifted polynomial.
    /// This is `none` if the committed polynomial does not
    /// enforce a strict degree bound.
    pub shifted_comm_var: Option<(usize, C)>,

    pub _affine: PhantomData<G>,
}

impl<G, C> AllocVar<LabeledCommitment<Commitment<G>>, ConstraintF<G>> for CommitmentVar<G, C>
where
    G: AffineCurve,
    C: CurveVar<G::Projective, ConstraintF<G>>,
{
    fn new_variable<T: Borrow<LabeledCommitment<Commitment<G>>>>(
        cs: impl Into<Namespace<ConstraintF<G>>>,
        f: impl FnOnce() -> Result<T, SynthesisError>,
        mode: AllocationMode,
    ) -> Result<Self, SynthesisError> {
        let ns = cs.into();
        let lc = f()?;

        let comm_var = C::new_variable(
            ns.clone(),
            || Ok(lc.borrow().commitment().comm.clone()),
            mode,
        )?;

        assert_eq!(
            lc.borrow().degree_bound().is_some(),
            lc.borrow().commitment().shifted_comm.is_some()
        );

        let shifted_comm_var = if lc.borrow().degree_bound().is_some() {
            Some((
                lc.borrow().degree_bound().unwrap(),
                C::new_variable(
                    ns.clone(),
                    || Ok(lc.borrow().commitment().shifted_comm.clone().unwrap()),
                    mode,
                )?,
            ))
        } else {
            None
        };

        Ok(Self {
            comm_var,
            shifted_comm_var,
            _affine: PhantomData,
        })
    }
}

pub struct ProofVar<G, C>
where
    G: AffineCurve,
    C: CurveVar<G::Projective, ConstraintF<G>>,
{
    /// Vector of left elements for each of the log_d iterations in `open`
    pub l_var_vec: Vec<C>,

    /// Vector of right elements for each of the log_d iterations within `open`
    pub r_var_vec: Vec<C>,

    /// Committer key from the last iteration within `open`
    pub final_comm_key_var: C,

    /// Coefficient from the last iteration within `open`
    pub c_var: NNFieldVar<G>,

    /// The first element is the commitment to the blinding polynomial.
    /// The second element is the linear combination of all the randomness
    /// used for commitments to the opened polynomials, along with the
    /// randomness used for the commitment to the hiding polynomial.
    pub hiding_var: Option<(C, Vec<Boolean<ConstraintF<G>>>)>,
}

impl<G, C> AllocVar<Proof<G>, ConstraintF<G>> for ProofVar<G, C>
where
    G: AffineCurve,
    C: CurveVar<G::Projective, ConstraintF<G>>,
{
    fn new_variable<T: Borrow<Proof<G>>>(
        cs: impl Into<Namespace<ConstraintF<G>>>,
        f: impl FnOnce() -> Result<T, SynthesisError>,
        mode: AllocationMode,
    ) -> Result<Self, SynthesisError> {
        let ns = cs.into();
        let proof = f()?;

        let l_var_vec: Vec<C> = proof
            .borrow()
            .l_vec
            .iter()
            .map(|elem| C::new_variable(ns.clone(), || Ok(elem.clone()), mode))
            .collect::<Result<Vec<C>, SynthesisError>>()?;

        let r_var_vec: Vec<C> = proof
            .borrow()
            .r_vec
            .iter()
            .map(|elem| C::new_variable(ns.clone(), || Ok(elem.clone()), mode))
            .collect::<Result<Vec<C>, SynthesisError>>()?;

        let final_comm_key_var = C::new_variable(
            ns.clone(),
            || Ok(proof.borrow().final_comm_key.clone()),
            mode,
        )?;

        let c_var =
            NNFieldVar::<G>::new_variable(ns.clone(), || Ok(proof.borrow().c.clone()), mode)?;

        assert_eq!(
            proof.borrow().hiding_comm.is_some(),
            proof.borrow().rand.is_some()
        );

        let hiding_var = if proof.borrow().hiding_comm.is_some() {
            let hiding_comm_var = C::new_variable(
                ns.clone(),
                || Ok(proof.borrow().hiding_comm.clone().unwrap()),
                mode,
            )?;

            let rand_bits_var =
                BitIteratorLE::without_trailing_zeros((&proof.borrow().rand.unwrap()).into_repr())
                    .map(|b| Boolean::new_variable(ns.clone(), || Ok(b), mode))
                    .collect::<Result<Vec<_>, SynthesisError>>()?;

            Some((hiding_comm_var, rand_bits_var))
        } else {
            None
        };

        Ok(Self {
            l_var_vec,
            r_var_vec,
            final_comm_key_var,
            c_var,
            hiding_var,
        })
    }
}

/// `SuccinctCheckPolynomial` is a succinctly-represented polynomial
/// generated from the `log_d` random oracle challenges generated in `open`.
/// It has the special property that can be evaluated in `O(log_d)` time.
pub struct SuccinctCheckPolynomialVar<G: AffineCurve>(pub Vec<NNFieldVar<G>>);

impl<G: AffineCurve> SuccinctCheckPolynomialVar<G> {
    /// Computes the coefficients of the underlying degree `d` polynomial.
    pub fn compute_coeff_vars(&self) -> Vec<NNFieldVar<G>> {
        let challenge_vars = &self.0;
        let log_d = challenge_vars.len();

        let mut coeff_vars = vec![NNFieldVar::<G>::one(); 1 << log_d];
        for (i, challenge) in challenge_vars.iter().enumerate() {
            let i = i + 1;
            let elem_degree = 1 << (log_d - i);
            for start in (elem_degree..coeff_vars.len()).step_by(elem_degree * 2) {
                for offset in 0..elem_degree {
                    coeff_vars[start + offset] *= challenge;
                }
            }
        }

        coeff_vars
    }

    /// Evaluate `self` at `point` in time `O(log_d)`.
    pub fn evaluate(&self, point: &NNFieldVar<G>) -> Result<NNFieldVar<G>, SynthesisError> {
        let challenges = &self.0;
        let log_d = challenges.len();

        let mut product = NNFieldVar::<G>::one();
        for (i, challenge) in challenges.iter().enumerate() {
            let i = i + 1;
            let elem_degree: u64 = (1 << (log_d - i)) as u64;
            let elem = point.pow_by_constant(&[elem_degree])?;
            product *= &(NNFieldVar::<G>::one() + &(elem * challenge));
        }

        Ok(product)
    }
}

pub struct CMCommitGadget<G, C>
where
    G: AffineCurve,
    C: CurveVar<G::Projective, ConstraintF<G>>,
{
    pub _affine: PhantomData<G>,
    pub _curve_var: PhantomData<C>,
}

impl<G, C> CMCommitGadget<G, C>
where
    G: AffineCurve,
    C: CurveVar<G::Projective, ConstraintF<G>>,
{
    pub fn commit(
        comm_key: &[C],
        scalars: &[Vec<Boolean<ConstraintF<G>>>],
        hiding: Option<(&C, &Vec<Boolean<ConstraintF<G>>>)>,
    ) -> Result<C, SynthesisError> {
        assert!(scalars.len() <= comm_key.len());

        let mut comm = C::zero();
        for (c, s) in comm_key.iter().zip(scalars) {
            comm += c.scalar_mul_le(s.iter())?;
        }

        if let Some((hiding_generator, randomizer)) = &hiding {
            comm += &hiding_generator.scalar_mul_le(randomizer.iter())?;
        }

        Ok(comm)
    }
}
