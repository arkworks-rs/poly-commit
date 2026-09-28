use crate::{
    kzh::data_structures::{
        Commitment, CommitmentState, CommitterKey, Proof, UniversalParams, VerifierKey,
    },
    Error,
};
use ark_ec::{pairing::Pairing, AffineRepr};
use ark_ff::Field;
use ark_poly::{DenseMultilinearExtension, MultilinearExtension};
use ark_std::vec::Vec;

/// Splits `num_vars` into `k` nonempty, balanced tensor blocks ordered
/// low-to-high, following ark-poly's variable order.
///
/// When the division is uneven, the earliest (lowest-variable) blocks receive
/// one additional variable.
pub(crate) fn balanced_block_sizes(num_vars: usize, k: usize) -> Result<Vec<usize>, Error> {
    if k < 2 || num_vars < k {
        return Err(Error::InvalidNumberOfVariables);
    }

    let quotient = num_vars / k;
    let remainder = num_vars % k;
    let mut sizes = vec![quotient; k];
    for size in sizes.iter_mut().take(remainder) {
        *size += 1;
    }
    Ok(sizes)
}

/// Converts variable counts into the corresponding Boolean-hypercube dimensions.
pub(crate) fn block_dimensions(num_vars_per_block: &[usize]) -> Result<Vec<usize>, Error> {
    num_vars_per_block
        .iter()
        .map(|&num_vars| {
            if num_vars == 0 || num_vars >= usize::BITS as usize {
                Err(Error::InvalidParameters(
                    "KZH block dimension does not fit in usize".into(),
                ))
            } else {
                Ok(1usize << num_vars)
            }
        })
        .collect()
}

/// Computes a checked product of dimensions.
pub(crate) fn product(values: &[usize]) -> Result<usize, Error> {
    values.iter().try_fold(1usize, |accumulator, &value| {
        accumulator.checked_mul(value).ok_or_else(|| {
            Error::InvalidParameters("KZH tensor dimension does not fit in usize".into())
        })
    })
}

/// Returns the lengths of the auxiliary tables that reduce single-opening work.
///
/// At level `j`, contracting a cached prefix table costs one group-scalar term
/// per entry in `prod(dimensions[..=j])`. Recommitting the current partially
/// evaluated tensor costs `prod(dimensions[j..])` terms. The first quantity is
/// increasing and the second decreasing, so useful auxiliary layers form a
/// prefix. A tie is left to the direct path to avoid equal-cost preprocessing
/// and storage.
pub(crate) fn auxiliary_prefix_lengths(dimensions: &[usize]) -> Result<Vec<usize>, Error> {
    opening_auxiliary_prefix_lengths(dimensions, 1)
}

/// Returns the auxiliary-table prefix worth using for a batched opening.
///
/// The commitment state is independent of the eventual batch size, so it
/// stores every table useful to a single opening. When opening `batch_size`
/// polynomials together, directly fusing the batching challenge and point
/// contraction costs `batch_size * prod(dimensions[..=j])` group-scalar terms
/// at level `j`. This function selects only the stored prefix for which that
/// remains strictly cheaper than recommitting the aggregated field tensor.
pub(crate) fn opening_auxiliary_prefix_lengths(
    dimensions: &[usize],
    batch_size: usize,
) -> Result<Vec<usize>, Error> {
    if batch_size == 0 {
        return Err(Error::InvalidParameters(
            "KZH opening batch must not be empty".into(),
        ));
    }

    let mut current_len = product(dimensions)?;
    let mut prefix_len = 1usize;
    let mut lengths = Vec::new();

    for &dimension in dimensions.iter().take(dimensions.len().saturating_sub(1)) {
        prefix_len = prefix_len.checked_mul(dimension).ok_or_else(|| {
            Error::InvalidParameters("KZH auxiliary table size does not fit in usize".into())
        })?;
        // `batch_size * prefix_len < current_len`, written with division to
        // avoid overflowing usize for large batches.
        if prefix_len > (current_len - 1) / batch_size {
            break;
        }
        lengths.push(prefix_len);
        current_len /= dimension;
    }
    Ok(lengths)
}

/// How an opening proof layer is constructed.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum OpeningLayerSource {
    /// Contract a commitment-time auxiliary table.
    Cached,
    /// Commit the current partially evaluated field tensor directly.
    Direct,
}

/// The exact online group-scalar work selected for one opening layer.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct OpeningLayerPlan {
    pub(crate) source: OpeningLayerSource,
    pub(crate) scalar_terms: usize,
}

/// Plans the cached/direct opening crossover used by the prover.
pub(crate) fn opening_work_plan(
    dimensions: &[usize],
    batch_size: usize,
) -> Result<Vec<OpeningLayerPlan>, Error> {
    let cached_levels = opening_auxiliary_prefix_lengths(dimensions, batch_size)?.len();
    let mut current_len = product(dimensions)?;
    let mut prefix_len = 1usize;
    let mut plan = Vec::with_capacity(dimensions.len().saturating_sub(1));

    for (level, &dimension) in dimensions
        .iter()
        .take(dimensions.len().saturating_sub(1))
        .enumerate()
    {
        prefix_len = prefix_len
            .checked_mul(dimension)
            .ok_or_else(|| Error::InvalidParameters("KZH opening prefix size overflow".into()))?;
        if level < cached_levels {
            // For a single polynomial, the first cached table already is the
            // first proof layer, so opening performs no online MSM here.
            let scalar_terms = if batch_size == 1 && level == 0 {
                0
            } else {
                batch_size.checked_mul(prefix_len).ok_or_else(|| {
                    Error::InvalidParameters("KZH cached opening work overflow".into())
                })?
            };
            plan.push(OpeningLayerPlan {
                source: OpeningLayerSource::Cached,
                scalar_terms,
            });
        } else {
            plan.push(OpeningLayerPlan {
                source: OpeningLayerSource::Direct,
                scalar_terms: current_len,
            });
        }
        current_len /= dimension;
    }
    Ok(plan)
}

/// Validates serialized family metadata and returns the block dimensions.
pub(crate) fn validate_family_metadata<const K: usize>(
    encoded_k: usize,
    num_vars: usize,
    num_vars_per_block: &[usize],
) -> Result<Vec<usize>, Error> {
    if encoded_k != K {
        return Err(Error::InvalidParameters(
            "serialized KZH family parameter does not match its Rust type".into(),
        ));
    }

    let expected = balanced_block_sizes(num_vars, K)?;
    if num_vars_per_block != expected.as_slice() {
        return Err(Error::InvalidParameters(
            "invalid KZH variable-block decomposition".into(),
        ));
    }
    block_dimensions(num_vars_per_block)
}

fn object_dimensions<const K: usize>(
    encoded_k: usize,
    num_vars: usize,
) -> Result<Vec<usize>, Error> {
    if encoded_k != K {
        return Err(Error::InvalidParameters(
            "serialized KZH family parameter does not match its Rust type".into(),
        ));
    }
    let sizes = balanced_block_sizes(num_vars, K)?;
    block_dimensions(&sizes)
}

fn validate_h_layers<G>(h: &[Vec<G>], dimensions: &[usize]) -> Result<(), Error> {
    if h.len() != dimensions.len() {
        return Err(Error::InvalidParameters(
            "incorrect number of KZH commitment-key layers".into(),
        ));
    }

    for (layer, suffix) in h
        .iter()
        .zip((0..dimensions.len()).map(|j| &dimensions[j..]))
    {
        let expected = product(suffix)?;
        if layer.len() != expected {
            return Err(Error::IncorrectCommitmentSize {
                encountered: layer.len(),
                expected,
            });
        }
    }
    Ok(())
}

fn validate_v_tau<G>(v_tau: &[Vec<G>], dimensions: &[usize]) -> Result<(), Error> {
    if v_tau.len() != dimensions.len() {
        return Err(Error::InvalidParameters(
            "incorrect number of KZH verifier-key layers".into(),
        ));
    }

    for (layer, &expected) in v_tau.iter().zip(dimensions) {
        if layer.len() != expected {
            return Err(Error::IncorrectCommitmentSize {
                encountered: layer.len(),
                expected,
            });
        }
    }
    Ok(())
}

/// Validates the dimensions of KZH universal parameters.
pub(crate) fn validate_params_shape<E: Pairing, const K: usize>(
    params: &UniversalParams<E, K>,
) -> Result<(), Error> {
    let dimensions =
        validate_family_metadata::<K>(params.k, params.num_vars, &params.num_vars_per_block)?;
    validate_h_layers(&params.h, &dimensions)?;
    validate_v_tau(
        &params.v_tau,
        &dimensions[..dimensions.len().saturating_sub(1)],
    )
}

/// Validates the dimensions of a KZH committer key.
pub(crate) fn validate_committer_key_shape<E: Pairing, const K: usize>(
    ck: &CommitterKey<E, K>,
) -> Result<(), Error> {
    let dimensions = validate_family_metadata::<K>(ck.k, ck.num_vars, &ck.num_vars_per_block)?;
    validate_h_layers(&ck.h, &dimensions)
}

/// Validates the dimensions of a KZH verifier key.
pub(crate) fn validate_verifier_key_shape<E: Pairing, const K: usize>(
    vk: &VerifierKey<E, K>,
) -> Result<(), Error> {
    let dimensions = validate_family_metadata::<K>(vk.k, vk.num_vars, &vk.num_vars_per_block)?;
    let expected_last = *dimensions.last().ok_or(Error::InvalidNumberOfVariables)?;
    if vk.h_last.len() != expected_last {
        return Err(Error::IncorrectCommitmentSize {
            encountered: vk.h_last.len(),
            expected: expected_last,
        });
    }
    validate_v_tau(&vk.v_tau, &dimensions[..dimensions.len().saturating_sub(1)])
}

/// Validates the metadata of a KZH commitment.
pub(crate) fn validate_commitment_shape<E: Pairing, const K: usize>(
    commitment: &Commitment<E, K>,
) -> Result<(), Error> {
    // `PCCommitment::empty()` cannot know the key's number of variables. A
    // zero group element with `num_vars == 0` is therefore an explicit empty
    // sentinel that is valid for every key in the same KZH family.
    if commitment.num_vars == 0 {
        if commitment.k != K {
            return Err(Error::InvalidParameters(
                "serialized KZH family parameter does not match its Rust type".into(),
            ));
        }
        if !commitment.comm.is_zero() {
            return Err(Error::InvalidParameters(
                "a KZH commitment with zero variables must be empty".into(),
            ));
        }
        return Ok(());
    }
    object_dimensions::<K>(commitment.k, commitment.num_vars).map(|_| ())
}

/// Validates the dimensions of cached KZH commitment state.
pub(crate) fn validate_state_shape<E: Pairing, const K: usize>(
    state: &CommitmentState<E, K>,
) -> Result<(), Error> {
    let dimensions = object_dimensions::<K>(state.k, state.num_vars)?;
    let expected_lengths = auxiliary_prefix_lengths(&dimensions)?;
    if state.auxiliary_tables.len() != expected_lengths.len() {
        return Err(Error::InvalidParameters(
            "incorrect number of KZH auxiliary tables".into(),
        ));
    }

    for (table, &expected) in state.auxiliary_tables.iter().zip(&expected_lengths) {
        if table.len() != expected {
            return Err(Error::IncorrectCommitmentSize {
                encountered: table.len(),
                expected,
            });
        }
    }
    Ok(())
}

/// Validates the dimensions of a KZH opening proof.
pub(crate) fn validate_proof_shape<E: Pairing, const K: usize>(
    proof: &Proof<E, K>,
) -> Result<(), Error> {
    let dimensions = object_dimensions::<K>(proof.k, proof.num_vars)?;
    let transition_count = dimensions.len().saturating_sub(1);
    if proof.layer_commitments.len() != transition_count {
        return Err(Error::InvalidParameters(
            "incorrect number of KZH proof layers".into(),
        ));
    }

    for (layer, &expected) in proof
        .layer_commitments
        .iter()
        .zip(&dimensions[..transition_count])
    {
        if layer.len() != expected {
            return Err(Error::IncorrectCommitmentSize {
                encountered: layer.len(),
                expected,
            });
        }
    }

    let expected_final = *dimensions.last().ok_or(Error::InvalidNumberOfVariables)?;
    if proof.final_evaluations.len() != expected_final {
        return Err(Error::IncorrectCommitmentSize {
            encountered: proof.final_evaluations.len(),
            expected: expected_final,
        });
    }
    Ok(())
}

/// Splits a little-endian MLE point into low-to-high variable blocks.
///
/// This is the same order used by ark-poly's `fix_variables`. For example,
/// sizes `[3, 2]` split `[x0, x1, x2, x3, x4]` into
/// `[[x0, x1, x2], [x3, x4]]`.
pub(crate) fn split_point_low_to_high<'a, F>(
    point: &'a [F],
    sizes: &[usize],
) -> Result<Vec<&'a [F]>, Error> {
    let total = sizes.iter().try_fold(0usize, |accumulator, &size| {
        if size == 0 {
            return Err(Error::InvalidNumberOfVariables);
        }
        accumulator.checked_add(size).ok_or_else(|| {
            Error::InvalidParameters("KZH point length does not fit in usize".into())
        })
    })?;
    if total != point.len() {
        return Err(Error::MismatchedNumVars {
            poly_nv: total,
            point_nv: point.len(),
        });
    }

    let mut start = 0usize;
    let mut blocks = Vec::with_capacity(sizes.len());
    for &size in sizes {
        let end = start
            .checked_add(size)
            .filter(|&end| end <= point.len())
            .ok_or_else(|| Error::InvalidParameters("invalid KZH point decomposition".into()))?;
        blocks.push(&point[start..end]);
        start = end;
    }
    Ok(blocks)
}

/// Builds the Boolean Lagrange evaluations using ark-poly's dense-MLE layout.
///
/// Ark-poly does not expose its internal equality-table routine. Constructing
/// the table through `DenseMultilinearExtension::concat` makes each point
/// coordinate the next variable in ark-poly's canonical little-endian order.
pub(crate) fn arkworks_lagrange_evaluations<F: Field>(point: &[F]) -> Result<Vec<F>, Error> {
    if point.len() >= usize::BITS as usize {
        return Err(Error::InvalidParameters(
            "KZH equality table dimension does not fit in usize".into(),
        ));
    }

    let mut equality = DenseMultilinearExtension::from_evaluations_vec(0, vec![F::one()]);
    for &coordinate in point {
        let num_vars = equality.num_vars();
        let dimension = 1usize << num_vars;
        let low_scalar = F::one() - coordinate;

        // ark-poly represents a zero polynomial as a zero-variate object when
        // using scalar multiplication by zero. Preserve the current arity so
        // `concat` still introduces exactly one new variable at Boolean points.
        let low = if low_scalar.is_zero() {
            DenseMultilinearExtension::from_evaluations_vec(num_vars, vec![F::zero(); dimension])
        } else {
            &equality * &low_scalar
        };
        let high = if coordinate.is_zero() {
            DenseMultilinearExtension::from_evaluations_vec(num_vars, vec![F::zero(); dimension])
        } else {
            &equality * &coordinate
        };
        equality = DenseMultilinearExtension::concat(&[&low, &high]);
    }
    Ok(equality.evaluations)
}

#[cfg(test)]
mod tests {
    use super::{
        auxiliary_prefix_lengths, balanced_block_sizes, opening_auxiliary_prefix_lengths,
        split_point_low_to_high,
    };
    use ark_bls12_381::Fr;

    #[test]
    fn balances_blocks_and_splits_low_to_high() {
        assert_eq!(balanced_block_sizes(8, 3).unwrap(), [3, 3, 2]);
        assert_eq!(balanced_block_sizes(6, 3).unwrap(), [2, 2, 2]);
        let sizes = balanced_block_sizes(8, 3).unwrap();

        let point: Vec<_> = (0u64..8).map(Fr::from).collect();
        let blocks = split_point_low_to_high(&point, &sizes).unwrap();
        assert_eq!(blocks[0], &point[0..3]);
        assert_eq!(blocks[1], &point[3..6]);
        assert_eq!(blocks[2], &point[6..8]);
    }

    #[test]
    fn caches_only_strictly_beneficial_auxiliary_prefixes() {
        assert_eq!(auxiliary_prefix_lengths(&[4, 4, 4, 4]).unwrap(), [4, 16]);
        assert_eq!(auxiliary_prefix_lengths(&[4, 4, 4]).unwrap(), [4]);
        assert_eq!(auxiliary_prefix_lengths(&[4, 4, 2, 2]).unwrap(), [4]);

        assert_eq!(
            opening_auxiliary_prefix_lengths(&[4, 4, 4, 4], 3).unwrap(),
            [4, 16]
        );
        assert_eq!(
            opening_auxiliary_prefix_lengths(&[4, 4, 4, 4], 4).unwrap(),
            [4]
        );
        assert!(opening_auxiliary_prefix_lengths(&[4, 4], 0).is_err());
    }
}
