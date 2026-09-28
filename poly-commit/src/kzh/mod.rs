//! The KZH-`k` multilinear polynomial commitment family.
//!
//! This module implements the non-hiding construction in Appendix C.1 of
//! [KZH-Fold][kzh], with the generic-opening auxiliary tables described in
//! Appendix E of [IronDict][irondict]. The const generic `K` is the tensor
//! arity. This documentation records only where the implementation differs
//! from the papers or from the crate's default PCS APIs.
//!
//! Tensor blocks follow arkworks' little-endian MLE order from `x_0` upward,
//! and partial evaluation uses [`MultilinearExtension::fix_variables`]. Same-
//! point batching is native. Multi-point and linear-combination queries use the
//! crate's default APIs, which issue one same-point opening per distinct point.
//! Hiding commitments, strict degree bounds, a Boolean-specialized opening API,
//! and an R1CS gadget are not implemented and are rejected where they would
//! otherwise apply.
//!
//! Commitment state stores only the auxiliary tables used by generic openings.
//! Tables that exist solely for free Boolean openings in the papers are omitted.
//!
//! # Security and setup
//!
//! `setup` produces an honestly generated, trusted SRS for one exact number of
//! variables. Every operation boundary validates dimensions and family
//! metadata, but this does not prove the algebraic consistency of an externally
//! supplied SRS or implement an updatable ceremony. The caller must use a
//! cryptographically secure RNG and destroy its recoverable state; materialized
//! trapdoor buffers are zeroized on drop. The security argument is the one in
//! the paper's algebraic-group model, under its `(q1, q2)` discrete-log and
//! Setup-find-representation assumptions. Batched openings additionally use
//! the sponge as a Fiat--Shamir random oracle with 128-bit challenges.
//!
//! Figure 12 presents equal tensor sides. This implementation uses the
//! corresponding balanced rectangular decomposition when the number of
//! variables is not divisible by `K`. In particular, `KZH2` is the `K = 2`
//! member of Appendix C.1's uniform family, not a claim of SRS compatibility
//! with the separately parameterized presentation in the paper's main text.
//!
//! [kzh]: https://eprint.iacr.org/2025/144
//! [irondict]: https://eprint.iacr.org/2025/1580
//! [`MultilinearExtension::fix_variables`]: ark_poly::MultilinearExtension::fix_variables

use crate::{Error, LabeledCommitment, LabeledPolynomial, PolynomialCommitment, CHALLENGE_SIZE};
use ark_crypto_primitives::sponge::{Absorb, CryptographicSponge};
use ark_ec::{
    pairing::Pairing, scalar_mul::BatchMulPreprocessing, AffineRepr, CurveGroup, VariableBaseMSM,
};
use ark_ff::{One, PrimeField, Zero};
use ark_poly::{DenseMultilinearExtension, MultilinearExtension};
use ark_serialize::serialize_to_vec;
use ark_std::{marker::PhantomData, rand::RngCore, string::ToString, vec::Vec, UniformRand};
#[cfg(feature = "parallel")]
use rayon::prelude::*;
use zeroize::Zeroizing;

mod data_structures;
pub use data_structures::*;
mod utils;
use utils::{
    arkworks_lagrange_evaluations, auxiliary_prefix_lengths, balanced_block_sizes,
    block_dimensions, opening_work_plan, product, split_point_low_to_high,
    validate_commitment_shape, validate_committer_key_shape, validate_params_shape,
    validate_proof_shape, validate_state_shape, validate_verifier_key_shape, OpeningLayerSource,
};

#[cfg(test)]
mod tests;

const BATCH_OPENING_DOMAIN_SEPARATOR: &[u8] = b"ark-poly-commit/KZH-k/batch-opening/v2";

struct PairingInputs<'a, Q> {
    v: &'a Q,
    v_tau: &'a [Vec<Q>],
}

/// KZH-`K`, a pairing-based multilinear polynomial commitment scheme.
///
/// `K` must satisfy `2 <= K <= num_vars`. It is part of the Rust type so that
/// different members of the KZH family cannot accidentally share keys. A
/// runtime copy of `K` is also serialized in every public object and checked
/// at each operation boundary. Setup, opening, and verification all use exactly
/// these `K` balanced tensor blocks.
pub struct KZH<E: Pairing, P: MultilinearExtension<E::ScalarField>, const K: usize> {
    _phantom: PhantomData<(E, P)>,
}

/// The two-dimensional member of the KZH family.
pub type KZH2<E, P> = KZH<E, P, 2>;

/// The three-dimensional member of the KZH family.
pub type KZH3<E, P> = KZH<E, P, 3>;

/// The four-dimensional member of the KZH family.
pub type KZH4<E, P> = KZH<E, P, 4>;

impl<E, P, const K: usize> KZH<E, P, K>
where
    E: Pairing,
    E::ScalarField: Absorb,
    P: MultilinearExtension<E::ScalarField>,
{
    fn invalid_input_length(message: &str) -> Error {
        Error::IncorrectInputLength(message.to_string())
    }

    fn msm(bases: &[E::G1Affine], scalars: &[E::ScalarField]) -> Result<E::G1, Error> {
        <E::G1 as VariableBaseMSM>::msm(bases, scalars).map_err(|_| {
            Self::invalid_input_length("KZH MSM bases and scalars have different lengths")
        })
    }

    fn msm_bigint(
        bases: &[E::G1Affine],
        scalars: &[<E::ScalarField as PrimeField>::BigInt],
    ) -> Result<E::G1, Error> {
        if bases.len() != scalars.len() {
            return Err(Self::invalid_input_length(
                "KZH MSM bases and scalars have different lengths",
            ));
        }
        Ok(<E::G1 as VariableBaseMSM>::msm_bigint(bases, scalars))
    }

    fn compute_auxiliary_tables(
        ck: &CommitterKey<E, K>,
        evaluations: &[<E::ScalarField as PrimeField>::BigInt],
    ) -> Result<Vec<Vec<E::G1Affine>>, Error> {
        let dimensions = block_dimensions(&ck.num_vars_per_block)?;
        let auxiliary_lengths = auxiliary_prefix_lengths(&dimensions)?;
        let mut auxiliary_tables = Vec::with_capacity(auxiliary_lengths.len());
        for (level, &expected_table_len) in auxiliary_lengths.iter().enumerate() {
            let bases = &ck.h[level + 1];
            let suffix_len = bases.len();
            if suffix_len == 0 || evaluations.len() % suffix_len != 0 {
                return Err(Self::invalid_input_length(
                    "KZH evaluation table does not match the commitment key",
                ));
            }

            let table_len = evaluations.len() / suffix_len;
            if table_len != expected_table_len {
                return Err(Self::invalid_input_length(
                    "KZH auxiliary table dimensions do not match the commitment key",
                ));
            }

            let compute_row =
                |row: usize, coefficients: &mut Vec<<E::ScalarField as PrimeField>::BigInt>| {
                    // Arkworks' low variables are the least-significant tensor
                    // axes. Hence a fixed low-prefix assignment is strided
                    // across the higher-variable suffix.
                    coefficients.clear();
                    coefficients.extend(
                        (0..suffix_len)
                            .map(|suffix_index| evaluations[row + table_len * suffix_index]),
                    );
                    Self::msm_bigint(bases, coefficients)
                };

            #[cfg(feature = "parallel")]
            let row_commitments = {
                // Arkworks' MSM already runs its bucket windows in parallel.
                // Parallelize across rows only for tables made of many smaller
                // MSMs, and bound the number of scratch buffers to
                // approximately one per worker.
                const MAX_PARALLEL_SUFFIX_LEN: usize = 1 << 12;
                let num_threads = rayon::current_num_threads();
                if num_threads > 1
                    && table_len >= num_threads
                    && suffix_len <= MAX_PARALLEL_SUFFIX_LEN
                {
                    let min_rows_per_task = table_len.div_ceil(num_threads);
                    (0..table_len)
                        .into_par_iter()
                        .with_min_len(min_rows_per_task)
                        .map_init(
                            || Vec::with_capacity(suffix_len),
                            |coefficients, row| compute_row(row, coefficients),
                        )
                        .collect::<Result<Vec<_>, _>>()?
                } else {
                    let mut rows = Vec::with_capacity(table_len);
                    let mut coefficients = Vec::with_capacity(suffix_len);
                    for row in 0..table_len {
                        rows.push(compute_row(row, &mut coefficients)?);
                    }
                    rows
                }
            };

            #[cfg(not(feature = "parallel"))]
            let row_commitments = {
                let mut rows = Vec::with_capacity(table_len);
                let mut coefficients = Vec::with_capacity(suffix_len);
                for row in 0..table_len {
                    rows.push(compute_row(row, &mut coefficients)?);
                }
                rows
            };
            auxiliary_tables.push(E::G1::normalize_batch(&row_commitments));
        }
        Ok(auxiliary_tables)
    }

    fn commitment_matches_num_vars(commitment: &Commitment<E, K>, num_vars: usize) -> bool {
        commitment.num_vars == num_vars || (commitment.num_vars == 0 && commitment.comm.is_zero())
    }

    fn cached_opening_layer(
        states: &[&CommitmentState<E, K>],
        challenges: &[E::ScalarField],
        prefix_weights: &[E::ScalarField],
        level: usize,
        dimension: usize,
        expected_scalar_terms: usize,
    ) -> Result<Vec<E::G1Affine>, Error> {
        // The cached table stores one contiguous column per current block
        // index. Fuse the polynomial-batching challenges with
        // eq(point_prefix), so no aggregated auxiliary table is allocated or
        // contracted a second time.
        let expected_len = prefix_weights.len().checked_mul(dimension).ok_or_else(|| {
            Error::InvalidParameters("KZH prefix contraction size overflow".to_string())
        })?;
        for state in states {
            if state.auxiliary_tables[level].len() != expected_len {
                return Err(Self::invalid_input_length(
                    "KZH auxiliary table does not match its prefix weights",
                ));
            }
        }

        if states.len() == 1 && level == 0 {
            // There is no prior point prefix to contract at the first level,
            // so the cached table is already the first proof layer.
            return Ok(states[0].auxiliary_tables[0].clone());
        }

        if states.len() == 1 {
            if expected_len != expected_scalar_terms {
                return Err(Error::InvalidParameters(
                    "KZH cached opening plan has an incorrect size".to_string(),
                ));
            }
            let projective = states[0].auxiliary_tables[level]
                .chunks_exact(prefix_weights.len())
                .map(|column| Self::msm(column, prefix_weights))
                .collect::<Result<Vec<_>, _>>()?;
            return Ok(E::G1::normalize_batch(&projective));
        }

        let combined_len = states
            .len()
            .checked_mul(prefix_weights.len())
            .ok_or_else(|| {
                Error::InvalidParameters("KZH batched prefix contraction size overflow".to_string())
            })?;
        if combined_len.checked_mul(dimension) != Some(expected_scalar_terms) {
            return Err(Error::InvalidParameters(
                "KZH batched opening plan has an incorrect size".to_string(),
            ));
        }

        let mut combined_weights = Vec::with_capacity(combined_len);
        for challenge in challenges {
            combined_weights.extend(
                prefix_weights
                    .iter()
                    .map(|prefix_weight| *challenge * prefix_weight),
            );
        }

        let mut projective = Vec::with_capacity(dimension);
        let mut combined_bases = Vec::with_capacity(combined_len);
        for current_index in 0..dimension {
            let start = current_index * prefix_weights.len();
            let end = start + prefix_weights.len();
            combined_bases.clear();
            for state in states {
                combined_bases.extend_from_slice(&state.auxiliary_tables[level][start..end]);
            }
            projective.push(Self::msm(&combined_bases, &combined_weights)?);
        }
        Ok(E::G1::normalize_batch(&projective))
    }

    fn direct_opening_layer(
        ck: &CommitterKey<E, K>,
        partially_evaluated: &DenseMultilinearExtension<E::ScalarField>,
        level: usize,
        dimension: usize,
        expected_scalar_terms: usize,
    ) -> Result<Vec<E::G1Affine>, Error> {
        if partially_evaluated.evaluations.len() != expected_scalar_terms {
            return Err(Error::InvalidParameters(
                "KZH direct opening plan has an incorrect size".to_string(),
            ));
        }

        let bases = &ck.h[level + 1];
        let suffix_len = bases.len();
        if suffix_len == 0
            || partially_evaluated.evaluations.len() % suffix_len != 0
            || partially_evaluated.evaluations.len() / suffix_len != dimension
        {
            return Err(Self::invalid_input_length(
                "KZH partial evaluation does not match the commitment key",
            ));
        }

        let mut projective = Vec::with_capacity(dimension);
        let mut coefficients = Vec::with_capacity(suffix_len);
        for current_index in 0..dimension {
            coefficients.clear();
            coefficients.extend((0..suffix_len).map(|suffix_index| {
                partially_evaluated.evaluations[current_index + dimension * suffix_index]
            }));
            projective.push(Self::msm(bases, &coefficients)?);
        }
        Ok(E::G1::normalize_batch(&projective))
    }

    fn open_evaluations(
        ck: &CommitterKey<E, K>,
        evaluations: Vec<E::ScalarField>,
        point: &[E::ScalarField],
        states: &[&CommitmentState<E, K>],
        challenges: &[E::ScalarField],
    ) -> Result<Proof<E, K>, Error> {
        validate_committer_key_shape(ck)?;
        if states.is_empty() || states.len() != challenges.len() {
            return Err(Self::invalid_input_length(
                "KZH states and batching challenges have different lengths or are empty",
            ));
        }
        for state in states {
            validate_state_shape(state)?;
            if state.num_vars != ck.num_vars {
                return Err(Self::invalid_input_length(
                    "KZH state and committer key support different numbers of variables",
                ));
            }
        }

        let dimensions = block_dimensions(&ck.num_vars_per_block)?;
        let expected_evaluations = product(&dimensions)?;
        if evaluations.len() != expected_evaluations {
            return Err(Self::invalid_input_length(
                "KZH polynomial has an incorrect evaluation-table length",
            ));
        }
        let point_blocks = split_point_low_to_high(point, &ck.num_vars_per_block)?;
        let work_plan = opening_work_plan(&dimensions, states.len())?;

        let mut proof_layers = Vec::with_capacity(work_plan.len());
        let mut partially_evaluated =
            DenseMultilinearExtension::from_evaluations_vec(ck.num_vars, evaluations);
        let mut prefix_weights = vec![E::ScalarField::one()];
        for (level, (point_block, layer_plan)) in point_blocks.iter().zip(&work_plan).enumerate() {
            let dimension = dimensions[level];
            let proof_layer = match layer_plan.source {
                OpeningLayerSource::Cached => Self::cached_opening_layer(
                    states,
                    challenges,
                    &prefix_weights,
                    level,
                    dimension,
                    layer_plan.scalar_terms,
                )?,
                OpeningLayerSource::Direct => Self::direct_opening_layer(
                    ck,
                    &partially_evaluated,
                    level,
                    dimension,
                    layer_plan.scalar_terms,
                )?,
            };
            proof_layers.push(proof_layer);
            partially_evaluated = partially_evaluated.fix_variables(point_block);

            if matches!(
                work_plan.get(level + 1),
                Some(next) if next.source == OpeningLayerSource::Cached
            ) {
                let block_weights = arkworks_lagrange_evaluations(point_block)?;
                let capacity = prefix_weights
                    .len()
                    .checked_mul(block_weights.len())
                    .ok_or_else(|| {
                        Error::InvalidParameters(
                            "KZH prefix weight table size overflow".to_string(),
                        )
                    })?;
                let mut next_weights = Vec::with_capacity(capacity);
                // In arkworks' little-endian order the already-fixed prefix is
                // the inner (faster) tensor axis and this new block is outer.
                for block_weight in &block_weights {
                    for prefix_weight in &prefix_weights {
                        next_weights.push(*prefix_weight * block_weight);
                    }
                }
                prefix_weights = next_weights;
            }
        }

        let proof = Proof {
            k: K,
            num_vars: ck.num_vars,
            layer_commitments: proof_layers,
            final_evaluations: partially_evaluated.evaluations,
        };
        validate_proof_shape(&proof)?;
        Ok(proof)
    }

    fn verify_opening_with_g2<Q>(
        vk: &VerifierKey<E, K>,
        pairing_v: &Q,
        pairing_v_tau: &[Vec<Q>],
        commitment: E::G1Affine,
        point: &[E::ScalarField],
        value: E::ScalarField,
        proof: &Proof<E, K>,
    ) -> Result<bool, Error>
    where
        Q: Clone + Into<E::G2Prepared>,
    {
        validate_verifier_key_shape(vk)?;
        if pairing_v_tau.len() != vk.v_tau.len()
            || pairing_v_tau
                .iter()
                .zip(&vk.v_tau)
                .any(|(prepared, affine)| prepared.len() != affine.len())
        {
            return Err(Error::InvalidParameters(
                "incorrect KZH prepared verifier-key dimensions".to_string(),
            ));
        }
        validate_proof_shape(proof)?;
        if proof.num_vars != vk.num_vars {
            return Err(Self::invalid_input_length(
                "KZH proof and verifier key support different numbers of variables",
            ));
        }

        let point_blocks = split_point_low_to_high(point, &vk.num_vars_per_block)?;
        let mut current_commitment = commitment;
        for (level, point_block) in point_blocks
            .iter()
            .take(point_blocks.len().saturating_sub(1))
            .enumerate()
        {
            let proof_layer = &proof.layer_commitments[level];

            let mut g1_terms = Vec::with_capacity(proof_layer.len() + 1);
            let mut g2_terms = Vec::with_capacity(proof_layer.len() + 1);
            g1_terms.push(current_commitment);
            g2_terms.push(pairing_v.clone());
            for (d_i, v_i) in proof_layer.iter().zip(&pairing_v_tau[level]) {
                g1_terms.push((-d_i.into_group()).into_affine());
                g2_terms.push(v_i.clone());
            }

            if !E::multi_pairing(g1_terms, g2_terms).is_zero() {
                return Ok(false);
            }

            let equality_vector = arkworks_lagrange_evaluations(point_block)?;
            current_commitment = Self::msm(proof_layer, &equality_vector)?.into_affine();
        }

        let alleged_last_commitment =
            Self::msm(&vk.h_last, &proof.final_evaluations)?.into_affine();
        if current_commitment != alleged_last_commitment {
            return Ok(false);
        }

        let last_point_block = point_blocks.last().ok_or(Error::InvalidNumberOfVariables)?;
        let final_polynomial = DenseMultilinearExtension::from_evaluations_vec(
            last_point_block.len(),
            proof.final_evaluations.clone(),
        );
        let alleged_value = final_polynomial.fix_variables(last_point_block)[0];
        Ok(alleged_value == value)
    }

    fn absorb_length_prefixed_bytes(sponge: &mut impl CryptographicSponge, bytes: &[u8]) {
        sponge.absorb(&(bytes.len() as u64).to_le_bytes().to_vec());
        sponge.absorb(&bytes.to_vec());
    }

    fn batch_challenges(
        sponge: &mut impl CryptographicSponge,
        commitments: &[&LabeledCommitment<Commitment<E, K>>],
        point: &[E::ScalarField],
        values: &[E::ScalarField],
        num_vars: usize,
    ) -> Result<Vec<E::ScalarField>, Error> {
        if commitments.is_empty() || commitments.len() != values.len() {
            return Err(Self::invalid_input_length(
                "KZH commitments and evaluations have different lengths or are empty",
            ));
        }

        Self::absorb_length_prefixed_bytes(sponge, BATCH_OPENING_DOMAIN_SEPARATOR);
        Self::absorb_length_prefixed_bytes(sponge, &(K as u64).to_le_bytes());
        Self::absorb_length_prefixed_bytes(sponge, &(num_vars as u64).to_le_bytes());
        Self::absorb_length_prefixed_bytes(sponge, &(commitments.len() as u64).to_le_bytes());
        for (commitment, value) in commitments.iter().zip(values) {
            Self::absorb_length_prefixed_bytes(sponge, commitment.label().as_bytes());
            let encoded = serialize_to_vec!(commitment.commitment().comm)
                .map_err(|_| Error::TranscriptError)?;
            Self::absorb_length_prefixed_bytes(sponge, &encoded);
            let encoded_value = serialize_to_vec!(*value).map_err(|_| Error::TranscriptError)?;
            Self::absorb_length_prefixed_bytes(sponge, &encoded_value);
        }
        sponge.absorb(&point.to_vec());

        let mut challenges = Vec::with_capacity(commitments.len());
        challenges.push(E::ScalarField::one());
        if commitments.len() > 1 {
            challenges.extend(sponge.squeeze_field_elements_with_sizes::<E::ScalarField>(
                &vec![CHALLENGE_SIZE; commitments.len() - 1],
            ));
        }
        Ok(challenges)
    }

    fn check_with_g2<'a, Q>(
        vk: &VerifierKey<E, K>,
        pairing_inputs: PairingInputs<'_, Q>,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Commitment<E, K>>>,
        point: &'a P::Point,
        values: impl IntoIterator<Item = E::ScalarField>,
        proof: &Proof<E, K>,
        sponge: &mut impl CryptographicSponge,
    ) -> Result<bool, Error>
    where
        Q: Clone + Into<E::G2Prepared>,
        Commitment<E, K>: 'a,
    {
        validate_verifier_key_shape(vk)?;
        if pairing_inputs.v_tau.len() != vk.v_tau.len()
            || pairing_inputs
                .v_tau
                .iter()
                .zip(&vk.v_tau)
                .any(|(prepared, affine)| prepared.len() != affine.len())
        {
            return Err(Error::InvalidParameters(
                "incorrect KZH prepared verifier-key dimensions".to_string(),
            ));
        }

        let commitments: Vec<_> = commitments.into_iter().collect();
        let values: Vec<_> = values.into_iter().collect();
        if commitments.is_empty() || commitments.len() != values.len() {
            return Err(Self::invalid_input_length(
                "KZH commitments and claimed values have different lengths or are empty",
            ));
        }
        for commitment in &commitments {
            if let Some(bound) = commitment.degree_bound() {
                return Err(Error::UnsupportedDegreeBound(bound));
            }
            validate_commitment_shape(commitment.commitment())?;
            if !Self::commitment_matches_num_vars(commitment.commitment(), vk.num_vars) {
                return Err(Self::invalid_input_length(
                    "KZH commitment and verifier key support different numbers of variables",
                ));
            }
        }

        split_point_low_to_high(point, &vk.num_vars_per_block)?;
        validate_proof_shape(proof)?;
        if proof.num_vars != vk.num_vars {
            return Err(Self::invalid_input_length(
                "KZH proof and verifier key support different numbers of variables",
            ));
        }

        // All validation is complete before the caller's sponge is mutated.
        // This path is shared by ordinary and prepared verification so their
        // Fiat--Shamir transcript and commitment aggregation cannot diverge.
        let challenges = Self::batch_challenges(sponge, &commitments, point, &values, vk.num_vars)?;
        let (aggregate_commitment, aggregate_value) = if commitments.len() == 1 {
            (commitments[0].commitment().comm.into_group(), values[0])
        } else {
            let mut aggregate_commitment = E::G1::zero();
            let mut aggregate_value = E::ScalarField::zero();
            for ((commitment, value), challenge) in commitments.iter().zip(values).zip(challenges) {
                aggregate_commitment += commitment.commitment().comm * challenge;
                aggregate_value += value * challenge;
            }
            (aggregate_commitment, aggregate_value)
        };

        Self::verify_opening_with_g2(
            vk,
            pairing_inputs.v,
            pairing_inputs.v_tau,
            aggregate_commitment.into_affine(),
            point,
            aggregate_value,
            proof,
        )
    }

    /// Checks a KZH opening using a verifier key whose `G2` inputs have
    /// already been prepared for pairing.
    ///
    /// Uses the same transcript, aggregation, pairing equations, and result as
    /// [`PolynomialCommitment::check`]. Prepare the key once when it will be
    /// reused.
    pub fn check_prepared<'a>(
        prepared_vk: &PreparedVerifierKey<E, K>,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Commitment<E, K>>>,
        point: &'a P::Point,
        values: impl IntoIterator<Item = E::ScalarField>,
        proof: &Proof<E, K>,
        sponge: &mut impl CryptographicSponge,
        _rng: Option<&mut dyn RngCore>,
    ) -> Result<bool, Error>
    where
        Commitment<E, K>: 'a,
    {
        Self::check_with_g2(
            prepared_vk.verifier_key(),
            PairingInputs {
                v: prepared_vk.prepared_v(),
                v_tau: prepared_vk.prepared_v_tau(),
            },
            commitments,
            point,
            values,
            proof,
            sponge,
        )
    }
}

impl<E, P, const K: usize> PolynomialCommitment<E::ScalarField, P> for KZH<E, P, K>
where
    E: Pairing,
    E::ScalarField: Absorb,
    P: MultilinearExtension<E::ScalarField>,
{
    type UniversalParams = UniversalParams<E, K>;
    type CommitterKey = CommitterKey<E, K>;
    type VerifierKey = VerifierKey<E, K>;
    type Commitment = Commitment<E, K>;
    type CommitmentState = CommitmentState<E, K>;
    type Proof = Proof<E, K>;
    type BatchProof = Vec<Self::Proof>;
    type Error = Error;

    fn setup<R: RngCore>(
        max_degree: usize,
        num_vars: Option<usize>,
        rng: &mut R,
    ) -> Result<Self::UniversalParams, Self::Error> {
        if max_degree == 0 {
            return Err(Error::DegreeIsZero);
        }
        if max_degree != 1 {
            return Err(Error::InvalidParameters(
                "KZH supports only multilinear polynomials of individual degree one".to_string(),
            ));
        }

        let num_vars = num_vars.ok_or(Error::InvalidNumberOfVariables)?;
        let num_vars_per_block = balanced_block_sizes(num_vars, K)?;
        let dimensions = block_dimensions(&num_vars_per_block)?;
        let num_blocks = dimensions.len();
        let total_evaluations = product(&dimensions)?;

        let g = E::G1::rand(rng);
        let v = E::G2::rand(rng);
        let tau = Zeroizing::new(
            dimensions
                .iter()
                .map(|&dimension| (0..dimension).map(|_| E::ScalarField::rand(rng)).collect())
                .collect::<Vec<Vec<_>>>(),
        );

        // Every H_j is a fixed-base table at g. Build each suffix in arkworks'
        // little-endian order, with its first/lowest tensor axis varying
        // fastest.
        let g_table = BatchMulPreprocessing::new(g, total_evaluations);
        let mut h = Vec::with_capacity(num_blocks);
        for suffix_start in 0..num_blocks {
            let mut suffix_exponents = Zeroizing::new(vec![E::ScalarField::one()]);
            for block_trapdoors in &tau[suffix_start..] {
                let mut next = Zeroizing::new(Vec::with_capacity(
                    block_trapdoors.len() * suffix_exponents.len(),
                ));
                for trapdoor in block_trapdoors {
                    for suffix_exponent in suffix_exponents.iter() {
                        next.push(*trapdoor * suffix_exponent);
                    }
                }
                suffix_exponents = next;
            }
            h.push(g_table.batch_mul(&suffix_exponents));
        }

        let transition_count = num_blocks.saturating_sub(1);
        let max_dimension = dimensions[..transition_count]
            .iter()
            .copied()
            .max()
            .unwrap_or(1);
        let v_affine = v.into_affine();
        let v_table = BatchMulPreprocessing::new(v, max_dimension);
        let v_tau = tau
            .iter()
            .take(transition_count)
            .map(|trapdoors| v_table.batch_mul(trapdoors))
            .collect();

        let params = UniversalParams {
            k: K,
            num_vars,
            num_vars_per_block,
            h,
            v: v_affine,
            v_tau,
        };
        validate_params_shape(&params)?;
        Ok(params)
    }

    fn trim(
        pp: &Self::UniversalParams,
        supported_degree: usize,
        supported_hiding_bound: usize,
        enforced_degree_bounds: Option<&[usize]>,
    ) -> Result<(Self::CommitterKey, Self::VerifierKey), Self::Error> {
        validate_params_shape(pp)?;
        if supported_degree == 0 {
            return Err(Error::InvalidParameters(
                "KZH supported degree must be one".to_string(),
            ));
        }
        if supported_degree > 1 {
            return Err(Error::TrimmingDegreeTooLarge);
        }
        if supported_hiding_bound != 0 {
            return Err(Error::InvalidParameters(
                "KZH does not support hiding commitments".to_string(),
            ));
        }
        if let Some(bounds) = enforced_degree_bounds {
            if bounds.is_empty() {
                return Err(Error::EmptyDegreeBounds);
            }
            return Err(Error::UnsupportedDegreeBound(bounds[0]));
        }

        let ck = CommitterKey {
            k: pp.k,
            num_vars: pp.num_vars,
            num_vars_per_block: pp.num_vars_per_block.clone(),
            h: pp.h.clone(),
        };
        let vk = VerifierKey {
            k: pp.k,
            num_vars: pp.num_vars,
            num_vars_per_block: pp.num_vars_per_block.clone(),
            h_last: pp.h.last().ok_or(Error::InvalidNumberOfVariables)?.clone(),
            v: pp.v,
            v_tau: pp.v_tau.clone(),
        };
        validate_committer_key_shape(&ck)?;
        validate_verifier_key_shape(&vk)?;
        Ok((ck, vk))
    }

    fn commit<'a>(
        ck: &Self::CommitterKey,
        polynomials: impl IntoIterator<Item = &'a LabeledPolynomial<E::ScalarField, P>>,
        _rng: Option<&mut dyn RngCore>,
    ) -> Result<
        (
            Vec<LabeledCommitment<Self::Commitment>>,
            Vec<Self::CommitmentState>,
        ),
        Self::Error,
    >
    where
        P: 'a,
    {
        validate_committer_key_shape(ck)?;
        let mut commitments = Vec::new();
        let mut states = Vec::new();
        for polynomial in polynomials {
            if let Some(bound) = polynomial.degree_bound() {
                return Err(Error::UnsupportedDegreeBound(bound));
            }
            if polynomial.is_hiding() {
                return Err(Error::InvalidParameters(
                    "KZH does not support hiding commitments".to_string(),
                ));
            }
            if polynomial.num_vars() != ck.num_vars {
                return Err(Error::MismatchedNumVars {
                    poly_nv: polynomial.num_vars(),
                    point_nv: ck.num_vars,
                });
            }

            let evaluations = polynomial.to_evaluations();
            let evaluation_bigints = ark_std::cfg_into_iter!(evaluations)
                .map(|evaluation| evaluation.into_bigint())
                .collect::<Vec<_>>();
            let commitment = Commitment {
                k: K,
                num_vars: ck.num_vars,
                comm: Self::msm_bigint(&ck.h[0], &evaluation_bigints)?.into_affine(),
            };
            let state = CommitmentState {
                k: K,
                num_vars: ck.num_vars,
                auxiliary_tables: Self::compute_auxiliary_tables(ck, &evaluation_bigints)?,
            };
            validate_commitment_shape(&commitment)?;
            validate_state_shape(&state)?;
            commitments.push(LabeledCommitment::new(
                polynomial.label().clone(),
                commitment,
                None,
            ));
            states.push(state);
        }
        Ok((commitments, states))
    }

    fn open<'a>(
        ck: &Self::CommitterKey,
        labeled_polynomials: impl IntoIterator<Item = &'a LabeledPolynomial<E::ScalarField, P>>,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Self::Commitment>>,
        point: &'a P::Point,
        sponge: &mut impl CryptographicSponge,
        states: impl IntoIterator<Item = &'a Self::CommitmentState>,
        _rng: Option<&mut dyn RngCore>,
    ) -> Result<Self::Proof, Self::Error>
    where
        P: 'a,
        Self::CommitmentState: 'a,
        Self::Commitment: 'a,
    {
        validate_committer_key_shape(ck)?;
        let polynomials: Vec<_> = labeled_polynomials.into_iter().collect();
        let commitments: Vec<_> = commitments.into_iter().collect();
        let states: Vec<_> = states.into_iter().collect();
        if polynomials.is_empty()
            || polynomials.len() != commitments.len()
            || polynomials.len() != states.len()
        {
            return Err(Self::invalid_input_length(
                "KZH opening inputs have different lengths or are empty",
            ));
        }
        for ((polynomial, commitment), state) in polynomials.iter().zip(&commitments).zip(&states) {
            if polynomial.label() != commitment.label() {
                return Err(Error::MismatchedLabels {
                    commitment_label: commitment.label().clone(),
                    polynomial_label: polynomial.label().clone(),
                });
            }
            if let Some(bound) = polynomial.degree_bound() {
                return Err(Error::UnsupportedDegreeBound(bound));
            }
            if polynomial.is_hiding() {
                return Err(Error::InvalidParameters(
                    "KZH does not support hiding commitments".to_string(),
                ));
            }
            if let Some(bound) = commitment.degree_bound() {
                return Err(Error::UnsupportedDegreeBound(bound));
            }
            if polynomial.num_vars() != ck.num_vars {
                return Err(Error::MismatchedNumVars {
                    poly_nv: polynomial.num_vars(),
                    point_nv: ck.num_vars,
                });
            }
            validate_commitment_shape(commitment.commitment())?;
            validate_state_shape(state)?;
            if !Self::commitment_matches_num_vars(commitment.commitment(), ck.num_vars)
                || state.num_vars != ck.num_vars
            {
                return Err(Self::invalid_input_length(
                    "KZH opening inputs support different numbers of variables",
                ));
            }
        }

        split_point_low_to_high(point, &ck.num_vars_per_block)?;
        let dimensions = block_dimensions(&ck.num_vars_per_block)?;
        let evaluation_count = product(&dimensions)?;
        let mut evaluation_tables = Vec::with_capacity(polynomials.len());
        for polynomial in &polynomials {
            let evaluations = polynomial.to_evaluations();
            if evaluations.len() != evaluation_count {
                return Err(Self::invalid_input_length(
                    "KZH polynomial has an incorrect evaluation-table length",
                ));
            }
            evaluation_tables.push(evaluations);
        }

        // Bind each claimed evaluation into the Fiat--Shamir challenges. If
        // values were omitted, a prover could offset false claims within a
        // batch while preserving the aggregate claim. All input validation is
        // complete before the caller's sponge is mutated.
        let claimed_values: Vec<_> = polynomials
            .iter()
            .map(|polynomial| polynomial.evaluate(point))
            .collect();
        let challenges =
            Self::batch_challenges(sponge, &commitments, point, &claimed_values, ck.num_vars)?;

        // The first Fiat--Shamir coefficient is one, so use its evaluation
        // table as the accumulator instead of multiplying it or allocating a
        // separate zero-filled table.
        let mut evaluation_tables = evaluation_tables.into_iter();
        let mut aggregate_evaluations = evaluation_tables.next().ok_or_else(|| {
            Self::invalid_input_length("KZH opening requires at least one polynomial")
        })?;
        for (evaluations, challenge) in evaluation_tables.zip(challenges.iter().skip(1)) {
            for (target, evaluation) in aggregate_evaluations.iter_mut().zip(evaluations) {
                *target += *challenge * evaluation;
            }
        }
        Self::open_evaluations(ck, aggregate_evaluations, point, &states, &challenges)
    }

    fn check<'a>(
        vk: &Self::VerifierKey,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Self::Commitment>>,
        point: &'a P::Point,
        values: impl IntoIterator<Item = E::ScalarField>,
        proof: &Self::Proof,
        sponge: &mut impl CryptographicSponge,
        _rng: Option<&mut dyn RngCore>,
    ) -> Result<bool, Self::Error>
    where
        Self::Commitment: 'a,
    {
        Self::check_with_g2(
            vk,
            PairingInputs {
                v: &vk.v,
                v_tau: &vk.v_tau,
            },
            commitments,
            point,
            values,
            proof,
            sponge,
        )
    }
}
