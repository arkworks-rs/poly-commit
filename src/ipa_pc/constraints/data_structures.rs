use crate::ipa_pc::VerifierKey;
use crate::ipa_pc::{Commitment, CommitterKey, Proof, SuccinctVerifierKey};
use crate::LabeledCommitment;
use ark_ec::AffineCurve;
use ark_ff::vec::Vec;
use ark_ff::One;
use ark_ff::PrimeField;
use ark_ff::{BitIteratorLE, Field};
use ark_nonnative_field::{NonNativeFieldMulResultVar, NonNativeFieldVar};
use ark_r1cs_std::alloc::{AllocVar, AllocationMode};
use ark_r1cs_std::bits::boolean::Boolean;
use ark_r1cs_std::fields::FieldVar;
use ark_r1cs_std::groups::CurveVar;
use ark_relations::r1cs::{Namespace, SynthesisError};
use ark_std::borrow::Borrow;
use ark_std::marker::PhantomData;

pub(crate) type ConstraintF<G> = <<G as AffineCurve>::BaseField as Field>::BasePrimeField;
pub(crate) type NNFieldVar<G> = NonNativeFieldVar<<G as AffineCurve>::ScalarField, ConstraintF<G>>;

/// The R1CS equivalent of `SuccinctVerifierKey`.
pub struct SuccinctVerifierKeyVar<G, C>
where
    G: AffineCurve,
    C: CurveVar<G::Projective, ConstraintF<G>>,
{
    /// A random group generator.
    pub h: C,

    /// A random group generator that is to be used to make a commitment hiding.
    pub s: C,

    /// The supported degree of the verifier key
    pub supported_degree: usize,

    #[doc(hidden)]
    pub _affine: PhantomData<G>,
}

impl<G, C> AllocVar<SuccinctVerifierKey<G>, ConstraintF<G>> for SuccinctVerifierKeyVar<G, C>
where
    G: AffineCurve,
    C: CurveVar<G::Projective, ConstraintF<G>>,
{
    fn new_variable<T: Borrow<SuccinctVerifierKey<G>>>(
        cs: impl Into<Namespace<ConstraintF<G>>>,
        f: impl FnOnce() -> Result<T, SynthesisError>,
        mode: AllocationMode,
    ) -> Result<Self, SynthesisError> {
        let ns = cs.into();
        let key = f()?;

        let supported_degree = key.borrow().supported_degree;
        let h = C::new_variable(ns.clone(), || Ok(key.borrow().h.clone()), mode)?;
        let s = C::new_variable(ns.clone(), || Ok(key.borrow().s.clone()), mode)?;

        Ok(Self {
            h,
            s,
            supported_degree,
            _affine: PhantomData,
        })
    }
}

/// `VerifierKey` is used to check evaluation proofs for a given commitment.
pub struct VerifierKeyVar<G, C>
where
    G: AffineCurve,
    C: CurveVar<G::Projective, ConstraintF<G>>,
{
    /// The key used to commit to polynomials.
    pub comm_key: Vec<C>,

    /// The verifier key for succinct check.
    pub svk: SuccinctVerifierKeyVar<G, C>,

    #[doc(hidden)]
    pub _affine: PhantomData<G>,
}

impl<G, C> VerifierKeyVar<G, C>
where
    G: AffineCurve,
    C: CurveVar<G::Projective, ConstraintF<G>>,
{
    /// Returns the maximum supported polynomial degree
    pub fn supported_degree(&self) -> usize {
        self.comm_key.len() - 1
    }
}

impl<G, C> AllocVar<VerifierKey<G>, ConstraintF<G>> for VerifierKeyVar<G, C>
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

        let comm_key = key
            .borrow()
            .comm_key
            .iter()
            .map(|elem| C::new_variable(ns.clone(), || Ok(elem.clone()), mode))
            .collect::<Result<Vec<C>, SynthesisError>>()?;

        let svk = SuccinctVerifierKeyVar::new_variable(
            ns.clone(),
            || Ok(key.borrow().svk.clone()),
            mode,
        )?;

        Ok(Self {
            comm_key,
            svk,
            _affine: PhantomData,
        })
    }
}

/// Commitment to a polynomial that optionally enforces a degree bound.
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
    pub comm: C,

    /// The degree bound and a Pedersen commitment to the shifted polynomial.
    /// This is `none` if the committed polynomial does not
    /// enforce a strict degree bound.
    pub shifted_comm: Option<(usize, C)>,

    #[doc(hidden)]
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

        let comm = C::new_variable(
            ns.clone(),
            || Ok(lc.borrow().commitment().comm.clone()),
            mode,
        )?;

        assert_eq!(
            lc.borrow().degree_bound().is_some(),
            lc.borrow().commitment().shifted_comm.is_some()
        );

        let shifted_comm = if lc.borrow().degree_bound().is_some() {
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
            comm,
            shifted_comm,
            _affine: PhantomData,
        })
    }
}

/// `Proof` is an evaluation proof that is output by `InnerProductArg::open`.
pub struct ProofVar<G, C>
where
    G: AffineCurve,
    C: CurveVar<G::Projective, ConstraintF<G>>,
{
    /// Vector of left elements for each of the log_d iterations in `open`
    pub l_vec: Vec<C>,

    /// Vector of right elements for each of the log_d iterations within `open`
    pub r_vec: Vec<C>,

    /// Committer key from the last iteration within `open`
    pub final_comm_key: C,

    /// Coefficient from the last iteration within `open`
    pub c: NNFieldVar<G>,

    /// The first element is the commitment to the blinding polynomial.
    /// The second element is the linear combination of all the randomness
    /// used for commitments to the opened polynomials, along with the
    /// randomness used for the commitment to the hiding polynomial.
    pub hiding: Option<(C, Vec<Boolean<ConstraintF<G>>>)>,
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

        let l_vec: Vec<C> = proof
            .borrow()
            .l_vec
            .iter()
            .map(|elem| C::new_variable(ns.clone(), || Ok(elem.clone()), mode))
            .collect::<Result<Vec<C>, SynthesisError>>()?;

        let r_vec: Vec<C> = proof
            .borrow()
            .r_vec
            .iter()
            .map(|elem| C::new_variable(ns.clone(), || Ok(elem.clone()), mode))
            .collect::<Result<Vec<C>, SynthesisError>>()?;

        let final_comm_key = C::new_variable(
            ns.clone(),
            || Ok(proof.borrow().final_comm_key.clone()),
            mode,
        )?;

        let c = NNFieldVar::<G>::new_variable(ns.clone(), || Ok(proof.borrow().c.clone()), mode)?;

        assert_eq!(
            proof.borrow().hiding_comm.is_some(),
            proof.borrow().rand.is_some()
        );

        let hiding = if proof.borrow().hiding_comm.is_some() {
            let hiding_comm = C::new_variable(
                ns.clone(),
                || Ok(proof.borrow().hiding_comm.clone().unwrap()),
                mode,
            )?;

            let rand_bits =
                BitIteratorLE::without_trailing_zeros((&proof.borrow().rand.unwrap()).into_repr())
                    .map(|b| Boolean::new_variable(ns.clone(), || Ok(b), mode))
                    .collect::<Result<Vec<_>, SynthesisError>>()?;

            Some((hiding_comm, rand_bits))
        } else {
            None
        };

        Ok(Self {
            l_vec,
            r_vec,
            final_comm_key,
            c,
            hiding,
        })
    }
}

/// `SuccinctCheckPolynomial` is a succinctly-represented polynomial
/// generated from the `log_d` random oracle challenges generated in `open`.
/// It has the special property that can be evaluated in `O(log_d)` time.
pub struct SuccinctCheckPolynomialVar<G: AffineCurve>(pub Vec<NNFieldVar<G>>);

impl<G: AffineCurve> SuccinctCheckPolynomialVar<G> {
    /// Computes the coefficients of the underlying degree `d` polynomial.
    pub fn compute_coeffs(&self) -> Vec<NNFieldVar<G>> {
        let challenges = &self.0;
        let log_d = challenges.len();

        let mut coeffs = vec![NNFieldVar::<G>::one(); 1 << log_d];
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
    pub fn evaluate(&self, point: &NNFieldVar<G>) -> Result<NNFieldVar<G>, SynthesisError> {
        let challenges = &self.0;
        let log_d = challenges.len();

        let mut point_2 = point.clone();
        let mut powers: Vec<NNFieldVar<G>> = Vec::with_capacity(log_d);
        powers.push(point_2.clone());
        for _ in 1..log_d {
            powers.push(point_2.square_in_place()?.clone());
        }
        powers.reverse();

        let one = NonNativeFieldMulResultVar::<G::ScalarField, ConstraintF<G>>::constant(
            G::ScalarField::one(),
        );
        let mut product = NNFieldVar::<G>::one();
        for (power, challenge) in challenges.iter().zip(&powers) {
            product *= &(&one + &(power.mul_without_reduce(challenge)?)).reduce()?;
        }

        Ok(product)
    }
}

/// The gadget for the [`IpaPC cm_commit`][cm_commit] function.
///
/// [cm_commit]: crate::ipa_pc::InnerProductArgPC::cm_commit
pub struct CMCommitGadget<G, C>
where
    G: AffineCurve,
    C: CurveVar<G::Projective, ConstraintF<G>>,
{
    _affine: PhantomData<G>,
    _curve: PhantomData<C>,
}

impl<G, C> CMCommitGadget<G, C>
where
    G: AffineCurve,
    C: CurveVar<G::Projective, ConstraintF<G>>,
{
    /// Commit to `scalars` using `comm_key`.
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
