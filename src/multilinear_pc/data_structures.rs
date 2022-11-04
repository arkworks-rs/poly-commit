use ark_ec::pairing::Pairing;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::vec::Vec;
#[allow(type_alias_bounds)]
/// Evaluations over {0,1}^n for G1
pub type EvaluationHyperCubeOnG1<E: Pairing> = Vec<E::G1Affine>;
#[allow(type_alias_bounds)]
/// Evaluations over {0,1}^n for G2
pub type EvaluationHyperCubeOnG2<E: Pairing> = Vec<E::G2Affine>;

/// Public Parameter used by prover
#[derive(CanonicalSerialize, CanonicalDeserialize, Clone, Debug)]
pub struct UniversalParams<E: Pairing> {
    /// number of variables
    pub num_vars: usize,
    /// `pp_{num_vars}`, `pp_{num_vars - 1}`, `pp_{num_vars - 2}`, ..., defined by XZZPD19
    pub powers_of_g: Vec<EvaluationHyperCubeOnG1<E>>,
    /// `pp_{num_vars}`, `pp_{num_vars - 1}`, `pp_{num_vars - 2}`, ..., defined by XZZPD19
    pub powers_of_h: Vec<EvaluationHyperCubeOnG2<E>>,
    /// generator for G1
    pub g: E::G1Affine,
    /// generator for G2
    pub h: E::G2Affine,
    /// g^randomness
    pub g_mask: Vec<E::G1Affine>,
}

/// Public Parameter used by prover
#[derive(CanonicalSerialize, CanonicalDeserialize, Clone, Debug)]
pub struct CommitterKey<E: Pairing> {
    /// number of variables
    pub nv: usize,
    /// pp_k defined by libra
    pub powers_of_g: Vec<EvaluationHyperCubeOnG1<E>>,
    /// pp_h defined by libra
    pub powers_of_h: Vec<EvaluationHyperCubeOnG2<E>>,
    /// generator for G1
    pub g: E::G1Affine,
    /// generator for G2
    pub h: E::G2Affine,
}

/// Public Parameter used by prover
#[derive(CanonicalSerialize, CanonicalDeserialize, Clone, Debug)]
pub struct VerifierKey<E: Pairing> {
    /// number of variables
    pub nv: usize,
    /// generator of G1
    pub g: E::G1Affine,
    /// generator of G2
    pub h: E::G2Affine,
    /// g^t1, g^t2, ...
    pub g_mask_random: Vec<E::G1Affine>,
}

#[derive(CanonicalSerialize, CanonicalDeserialize, Clone, Debug)]
/// commitment
pub struct Commitment<E: Pairing> {
    /// number of variables
    pub nv: usize,
    /// product of g as described by the vRAM paper
    pub g_product: E::G1Affine,
}

#[derive(CanonicalSerialize, CanonicalDeserialize, Clone, Debug)]
/// proof of opening
pub struct Proof<E: Pairing> {
    /// Evaluation of quotients
    pub proofs: Vec<E::G2Affine>,
}
