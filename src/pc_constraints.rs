use algebra::PrimeField;
use core::marker::Sized;
use r1cs_core::{ConstraintSystem, SynthesisError};
use r1cs_std::prelude::*;

use crate::data_structures::{LabeledCommitment, LinearCombination};
use crate::{BTreeMap, BTreeSet, String, Vec};
use crate::{BatchLCProof, PolynomialCommitment};
use nonnative::{NonNativeFieldGadget, NonNativeFieldParams};

/// A generic gadget for the prepared* structures
pub trait PrepareGadget<UNPREPARED, ConstraintF: PrimeField>: Sized {
    /// prepare from an unprepared element
    fn prepare<CS: ConstraintSystem<ConstraintF>>(
        cs: CS,
        unprepared: &UNPREPARED,
    ) -> Result<Self, SynthesisError>;
}

/// Describes the interface for a gadget for a `PolynomialCommitment`
/// verifier.
pub trait PCCheckGadget<
    PCF: PrimeField,
    PC: PolynomialCommitment<PCF>,
    NP: NonNativeFieldParams<BaseField = ConstraintF, TargetField = PCF>,
    ConstraintF: PrimeField,
>
{
    /// An allocated version of `PC::VerifierKey`.
    type VerifierKeyGadget: AllocGadget<PC::VerifierKey, ConstraintF>
        + Clone
        + ToConstraintFieldGadget<ConstraintF>
        + ToBytesGadget<ConstraintF>;
    /// An allocated version of `PC::PreparedVerifierKey`.
    type PreparedVerifierKeyGadget: AllocGadget<PC::PreparedVerifierKey, ConstraintF>
        + PrepareGadget<Self::VerifierKeyGadget, ConstraintF>;
    /// An allocated version of `PC::Commitment`.
    type CommitmentGadget: AllocGadget<PC::Commitment, ConstraintF>
        + Clone
        + ToConstraintFieldGadget<ConstraintF>
        + ToBytesGadget<ConstraintF>;
    /// An allocated version of `PC::PreparedCommitment`.
    type PreparedCommitmentGadget: AllocGadget<PC::PreparedCommitment, ConstraintF>
        + PrepareGadget<Self::CommitmentGadget, ConstraintF>
        + Clone;
    /// An allocated version of `LabeledCommitment<PC::Commitment>`.
    type LabeledCommitmentGadget: AllocGadget<LabeledCommitment<PC::Commitment>, ConstraintF>
        + Clone;
    /// A prepared, allocated version of `LabeledCommitment<PC::Commitment>`.
    type PreparedLabeledCommitmentGadget: Clone;
    /// An allocated version of `PC::Proof`.
    type ProofGadget: AllocGadget<PC::Proof, ConstraintF> + Clone;

    /// An allocated version of `LinearCombination`.
    type LinearCombinationGadget: AllocGadget<LinearCombination<PCF>, ConstraintF> + Clone;
    /// An allocated version of `PC::BatchLCProof`.
    type BatchLCProofGadget: AllocGadget<BatchLCProof<PCF, PC>, ConstraintF> + Clone;

    /// Add to `CS` new constraints that check that `proof_i` is a valid evaluation
    /// proof at `point_i` for the polynomial in `commitment_i`.
    fn batch_check_evaluations<CS: ConstraintSystem<ConstraintF>>(
        cs: CS,
        verification_key: &Self::VerifierKeyGadget,
        commitments: &Vec<Self::LabeledCommitmentGadget>,
        query_set: &QuerySetGadget<NP>,
        evaluations: &EvaluationsGadget<NP>,
        proofs: &Vec<Self::ProofGadget>,
        opening_challenges: &Vec<NonNativeFieldGadget<NP>>,
        opening_challenges_bits: &Vec<Vec<Boolean>>,
        batching_rands: &Vec<NonNativeFieldGadget<NP>>,
        batching_rands_bits: &Vec<Vec<Boolean>>,
    ) -> Result<(), SynthesisError>;

    /// The prepared version of `batch_check_evaluations`.
    fn prepared_batch_check_evaluations<CS: ConstraintSystem<ConstraintF>>(
        cs: CS,
        prepared_verification_key: &Self::PreparedVerifierKeyGadget,
        prepared_commitments: &Vec<Self::PreparedLabeledCommitmentGadget>,
        query_set: &QuerySetGadget<NP>,
        evaluations: &EvaluationsGadget<NP>,
        proofs: &Vec<Self::ProofGadget>,
        opening_challenges: &Vec<NonNativeFieldGadget<NP>>,
        opening_challenges_bits: &Vec<Vec<Boolean>>,
        batching_rands: &Vec<NonNativeFieldGadget<NP>>,
        batching_rands_bits: &Vec<Vec<Boolean>>,
    ) -> Result<(), SynthesisError>;

    /// The prepared, allocated version of `check_combinations`.
    fn prepared_check_combinations<CS: ConstraintSystem<ConstraintF>>(
        cs: CS,
        prepared_verification_key: &Self::PreparedVerifierKeyGadget,
        linear_combinations: &Vec<Self::LinearCombinationGadget>,
        prepared_commitments: &Vec<Self::PreparedLabeledCommitmentGadget>,
        query_set: &QuerySetGadget<NP>,
        evaluations: &EvaluationsGadget<NP>,
        proof: &Self::BatchLCProofGadget,
        opening_challenges: &Vec<NonNativeFieldGadget<NP>>,
        opening_challenges_bits: &Vec<Vec<Boolean>>,
        batching_rands: &Vec<NonNativeFieldGadget<NP>>,
        batching_rands_bits: &Vec<Vec<Boolean>>,
    ) -> Result<(), SynthesisError>;

    /// Create the labeled commitment gadget from the commitment gadget
    fn create_labeled_commitment_gadget(
        label: String,
        commitment: Self::CommitmentGadget,
        degree_bound: Option<UInt32>,
    ) -> Self::LabeledCommitmentGadget;

    /// Create the prepared labeled commitment gadget from the commitment gadget
    fn create_prepared_labeled_commitment_gadget(
        label: String,
        commitment: Self::PreparedCommitmentGadget,
        degree_bound: Option<UInt32>,
    ) -> Self::PreparedLabeledCommitmentGadget;
}

/// An allocated version of `QuerySet`.
pub type QuerySetGadget<NP> = BTreeSet<(String, NonNativeFieldGadget<NP>)>;
/// An allocated version of `Evaluations`.
pub type EvaluationsGadget<NP> =
    BTreeMap<(String, NonNativeFieldGadget<NP>), NonNativeFieldGadget<NP>>;
