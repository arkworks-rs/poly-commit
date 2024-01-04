use crate::{
    data_structures::LabeledCommitment, BatchLCProof, LCTerm, LinearCombination,
    PCPreparedCommitment, PCPreparedVerifierKey, PolynomialCommitment, String, Vec,
};
use ark_crypto_primitives::sponge::CryptographicSponge;
use ark_ff::PrimeField;
use ark_poly::Polynomial;
use ark_r1cs_std::fields::emulated_fp::EmulatedFpVar;
use ark_r1cs_std::{fields::fp::FpVar, prelude::*};
use ark_relations::r1cs::{ConstraintSystemRef, Namespace, Result as R1CSResult, SynthesisError};
use ark_std::{borrow::Borrow, cmp::Eq, cmp::PartialEq, hash::Hash, marker::Sized};
use hashbrown::{HashMap, HashSet};

/// Define the minimal interface of prepared allocated structures.
pub trait PrepareGadget<Unprepared, ConstraintF: PrimeField>: Sized {
    /// Prepare from an unprepared element.
    fn prepare(unprepared: &Unprepared) -> R1CSResult<Self>;
}

/// A coefficient of `LinearCombination`.
#[derive(Clone)]
pub enum LinearCombinationCoeffVar<TargetField: PrimeField, BaseField: PrimeField> {
    /// Coefficient 1.
    One,
    /// Coefficient -1.
    MinusOne,
    /// Other coefficient, represented as a "emulated" field element.
    Var(EmulatedFpVar<TargetField, BaseField>),
}

/// An allocated version of `LinearCombination`.
#[derive(Clone)]
pub struct LinearCombinationVar<TargetField: PrimeField, BaseField: PrimeField> {
    /// The label.
    pub label: String,
    /// The linear combination of `(coeff, poly_label)` pairs.
    pub terms: Vec<(LinearCombinationCoeffVar<TargetField, BaseField>, LCTerm)>,
}

impl<TargetField: PrimeField, BaseField: PrimeField>
    AllocVar<LinearCombination<TargetField>, BaseField>
    for LinearCombinationVar<TargetField, BaseField>
{
    fn new_variable<T>(
        cs: impl Into<Namespace<BaseField>>,
        val: impl FnOnce() -> Result<T, SynthesisError>,
        mode: AllocationMode,
    ) -> R1CSResult<Self>
    where
        T: Borrow<LinearCombination<TargetField>>,
    {
        let LinearCombination { label, terms } = val()?.borrow().clone();

        let ns = cs.into();
        let cs = ns.cs();

        let new_terms: Vec<(LinearCombinationCoeffVar<TargetField, BaseField>, LCTerm)> = terms
            .iter()
            .map(|term| {
                let (f, lc_term) = term;

                let fg =
                    EmulatedFpVar::new_variable(ark_relations::ns!(cs, "term"), || Ok(f), mode)
                        .unwrap();

                (LinearCombinationCoeffVar::Var(fg), lc_term.clone())
            })
            .collect();

        Ok(Self {
            label,
            terms: new_terms,
        })
    }
}

#[derive(Clone, Debug)]
/// A collection of random data used in the polynomial commitment checking.
pub struct PCCheckRandomDataVar<TargetField: PrimeField, BaseField: PrimeField> {
    /// Opening challenges.
    /// The prover and the verifier MUST use the same opening challenges.
    pub opening_challenges: Vec<EmulatedFpVar<TargetField, BaseField>>,
    /// Bit representations of the opening challenges.
    pub opening_challenges_bits: Vec<Vec<Boolean<BaseField>>>,
    /// Batching random numbers.
    /// The verifier can choose these numbers freely, as long as they are random.
    pub batching_rands: Vec<EmulatedFpVar<TargetField, BaseField>>,
    /// Bit representations of the batching random numbers.
    pub batching_rands_bits: Vec<Vec<Boolean<BaseField>>>,
}

/// Describes the interface for a gadget for a `PolynomialCommitment`
/// verifier.
pub trait PCCheckVar<
    PCF: PrimeField,
    P: Polynomial<PCF>,
    PC: PolynomialCommitment<PCF, P, S>,
    ConstraintF: PrimeField,
    S: CryptographicSponge,
>: Clone
{
    /// The prepared verifier key for the scheme; used to check an evaluation proof.
    type PreparedVerifierKey: PCPreparedVerifierKey<PC::VerifierKey> + Clone;
    /// The prepared commitment to a polynomial.
    type PreparedCommitment: PCPreparedCommitment<PC::Commitment>;
    /// An allocated version of `PC::VerifierKey`.
    type VerifierKeyVar: AllocVar<PC::VerifierKey, ConstraintF> + Clone;
    /// An allocated version of `PC::PreparedVerifierKey`.
    type PreparedVerifierKeyVar: AllocVar<Self::PreparedVerifierKey, ConstraintF>
        + Clone
        + PrepareGadget<Self::VerifierKeyVar, ConstraintF>;
    /// An allocated version of `PC::Commitment`.
    type CommitmentVar: AllocVar<PC::Commitment, ConstraintF> + Clone;
    /// An allocated version of `PC::PreparedCommitment`.
    type PreparedCommitmentVar: AllocVar<Self::PreparedCommitment, ConstraintF>
        + PrepareGadget<Self::CommitmentVar, ConstraintF>
        + Clone;
    /// An allocated version of `LabeledCommitment<PC::Commitment>`.
    type LabeledCommitmentVar: AllocVar<LabeledCommitment<PC::Commitment>, ConstraintF> + Clone;
    /// A prepared, allocated version of `LabeledCommitment<PC::Commitment>`.
    type PreparedLabeledCommitmentVar: Clone;
    /// An allocated version of `PC::Proof`.
    type ProofVar: AllocVar<PC::Proof, ConstraintF> + Clone;

    /// An allocated version of `PC::BatchLCProof`.
    type BatchLCProofVar: AllocVar<BatchLCProof<PCF, PC::BatchProof>, ConstraintF> + Clone;

    /// Add to `ConstraintSystemRef<ConstraintF>` new constraints that check that `proof_i` is a valid evaluation
    /// proof at `point_i` for the polynomial in `commitment_i`.
    fn batch_check_evaluations(
        cs: ConstraintSystemRef<ConstraintF>,
        verification_key: &Self::VerifierKeyVar,
        commitments: &[Self::LabeledCommitmentVar],
        query_set: &QuerySetVar<PCF, ConstraintF>,
        evaluations: &EvaluationsVar<PCF, ConstraintF>,
        proofs: &[Self::ProofVar],
        rand_data: &PCCheckRandomDataVar<PCF, ConstraintF>,
    ) -> R1CSResult<Boolean<ConstraintF>>;

    /// Add to `ConstraintSystemRef<ConstraintF>` new constraints that conditionally check that `proof` is a valid evaluation
    /// proof at the points in `query_set` for the combinations `linear_combinations`.
    fn prepared_check_combinations(
        cs: ConstraintSystemRef<ConstraintF>,
        prepared_verification_key: &Self::PreparedVerifierKeyVar,
        linear_combinations: &[LinearCombinationVar<PCF, ConstraintF>],
        prepared_commitments: &[Self::PreparedLabeledCommitmentVar],
        query_set: &QuerySetVar<PCF, ConstraintF>,
        evaluations: &EvaluationsVar<PCF, ConstraintF>,
        proof: &Self::BatchLCProofVar,
        rand_data: &PCCheckRandomDataVar<PCF, ConstraintF>,
    ) -> R1CSResult<Boolean<ConstraintF>>;

    /// Create the labeled commitment gadget from the commitment gadget
    fn create_labeled_commitment(
        label: String,
        commitment: Self::CommitmentVar,
        degree_bound: Option<FpVar<ConstraintF>>,
    ) -> Self::LabeledCommitmentVar;

    /// Create the prepared labeled commitment gadget from the commitment gadget
    fn create_prepared_labeled_commitment(
        label: String,
        commitment: Self::PreparedCommitmentVar,
        degree_bound: Option<FpVar<ConstraintF>>,
    ) -> Self::PreparedLabeledCommitmentVar;
}

#[derive(Clone, Hash, PartialEq, Eq)]
/// A labeled point variable, for queries to a polynomial commitment.
pub struct LabeledPointVar<TargetField: PrimeField, BaseField: PrimeField> {
    /// The label of the point.
    /// MUST be a unique identifier in a query set.
    pub name: String,
    /// The point value.
    pub value: EmulatedFpVar<TargetField, BaseField>,
}

/// An allocated version of `QuerySet`.
#[derive(Clone)]
pub struct QuerySetVar<TargetField: PrimeField, BaseField: PrimeField>(
    pub HashSet<(String, LabeledPointVar<TargetField, BaseField>)>,
);

/// An allocated version of `Evaluations`.
#[derive(Clone)]
pub struct EvaluationsVar<TargetField: PrimeField, BaseField: PrimeField>(
    pub HashMap<LabeledPointVar<TargetField, BaseField>, EmulatedFpVar<TargetField, BaseField>>,
);

impl<TargetField: PrimeField, BaseField: PrimeField> EvaluationsVar<TargetField, BaseField> {
    /// find the evaluation result
    pub fn get_lc_eval(
        &self,
        lc_string: &str,
        point: &EmulatedFpVar<TargetField, BaseField>,
    ) -> Result<EmulatedFpVar<TargetField, BaseField>, SynthesisError> {
        let key = LabeledPointVar::<TargetField, BaseField> {
            name: String::from(lc_string),
            value: point.clone(),
        };
        Ok(self.0.get(&key).map(|v| (*v).clone()).unwrap())
    }
}
