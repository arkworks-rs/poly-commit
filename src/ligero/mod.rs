use ark_crypto_primitives::crh::{CRHScheme, TwoToOneCRHScheme};
use ark_crypto_primitives::{
    merkle_tree::Config,
    sponge::{Absorb, CryptographicSponge},
};
use ark_ff::PrimeField;
use ark_poly::DenseUVPolynomial;
use ark_std::borrow::Borrow;
use ark_std::fmt::Debug;
use ark_std::marker::PhantomData;
use ark_std::rand::RngCore;
use digest::Digest;

use crate::data_structures::PCRandomness;
use crate::ligero::utils::{inner_product, reed_solomon, IOPTranscript};
use crate::{Error, LabeledCommitment, LabeledPolynomial, PCUniversalParams, PolynomialCommitment};

mod utils;

mod data_structures;
use data_structures::*;

pub use data_structures::{
    Ligero, LigeroPCCommitterKey, LigeroPCProof, LigeroPCUniversalParams, LigeroPCVerifierKey,
};

use utils::{calculate_t, get_indices_from_transcript, hash_column};

mod tests;

impl<F, P, S, C, D> PolynomialCommitment<F, P, S> for Ligero<F, C, D, S, P>
where
    F: PrimeField,
    P: DenseUVPolynomial<F>,
    S: CryptographicSponge,
    C: Config + 'static,
    Vec<u8>: Borrow<C::Leaf>,
    <<C as Config>::TwoToOneHash as TwoToOneCRHScheme>::Parameters: Debug,
    <<C as Config>::LeafHash as CRHScheme>::Parameters: Debug,
    C::InnerDigest: Absorb,
    D: Digest,
{
    type UniversalParams = LigeroPCUniversalParams<F, C>;

    type CommitterKey = LigeroPCCommitterKey<F, C>;

    type VerifierKey = LigeroPCVerifierKey<F, C>;

    type PreparedVerifierKey = LigeroPCPreparedVerifierKey;

    type Commitment = LigeroPCCommitment<C>;

    type PreparedCommitment = LigeroPCPreparedCommitment;

    type Randomness = LigeroPCRandomness;

    type Proof = LPCPArray<F, C>;

    type BatchProof = Vec<Self::Proof>;

    type Error = Error;

    /// This is only a default setup with reasonable parameters.
    /// To create your own public parameters (from which vk/ck can be derived by `trim`),
    /// see the documentation for `LigeroPCUniversalParams`.
    fn setup<R: RngCore>(
        max_degree: usize,
        _num_vars: Option<usize>,
        rng: &mut R,
    ) -> Result<Self::UniversalParams, Self::Error> {
        let leaf_hash_params = <C::LeafHash as CRHScheme>::setup(rng).unwrap();
        let two_to_one_params = <C::TwoToOneHash as TwoToOneCRHScheme>::setup(rng)
            .unwrap()
            .clone();
        let pp = Self::UniversalParams::new(128, 4, true, leaf_hash_params, two_to_one_params);
        let real_max_degree = pp.max_degree();
        if max_degree > real_max_degree || real_max_degree == 0 {
            return Err(Error::InvalidParameters(
                "This field is not suitable for the proposed parameters".to_string(),
            ));
        }
        Ok(pp)
    }

    fn trim(
        pp: &Self::UniversalParams,
        _supported_degree: usize,
        _supported_hiding_bound: usize,
        _enforced_degree_bounds: Option<&[usize]>,
    ) -> Result<(Self::CommitterKey, Self::VerifierKey), Self::Error> {
        if pp.max_degree() == 0 {
            return Err(Error::InvalidParameters(
                "This field is not suitable for the proposed parameters".to_string(),
            ));
        }
        let ck = LigeroPCCommitterKey::<F, C> {
            _field: PhantomData,
            sec_param: pp.sec_param,
            rho_inv: pp.rho_inv,
            leaf_hash_params: pp.leaf_hash_params.clone(),
            two_to_one_params: pp.two_to_one_params.clone(),
            check_well_formedness: pp.check_well_formedness,
        };
        let vk = LigeroPCVerifierKey::<F, C> {
            _field: PhantomData,
            sec_param: pp.sec_param,
            rho_inv: pp.rho_inv,
            leaf_hash_params: pp.leaf_hash_params.clone(),
            two_to_one_params: pp.two_to_one_params.clone(),
            check_well_formedness: pp.check_well_formedness,
        };
        Ok((ck, vk))
    }

    fn commit<'a>(
        ck: &Self::CommitterKey,
        polynomials: impl IntoIterator<Item = &'a LabeledPolynomial<F, P>>,
        _rng: Option<&mut dyn RngCore>,
    ) -> Result<
        (
            Vec<LabeledCommitment<Self::Commitment>>,
            Vec<Self::Randomness>,
        ),
        Self::Error,
    >
    where
        P: 'a,
    {
        let mut commitments = Vec::new();

        for labeled_polynomial in polynomials.into_iter() {
            let polynomial = labeled_polynomial.polynomial();

            // 1. Arrange the coefficients of the polynomial into a matrix,
            // and apply Reed-Solomon encoding to get `ext_mat`.
            let (mat, ext_mat) = Self::compute_matrices(polynomial, ck.rho_inv);

            // 2. Create the Merkle tree from the hashes of each column.
            let col_tree =
                Self::create_merkle_tree(&ext_mat, &ck.leaf_hash_params, &ck.two_to_one_params);

            // 3. Obtain the MT root and add it to the transcript.
            let root = col_tree.root();

            let mut transcript: IOPTranscript<F> = IOPTranscript::new(b"transcript");

            transcript
                .append_serializable_element(b"root", &root)
                .map_err(|_| Error::TranscriptError)?;

            let n_rows = mat.n;
            let n_cols = mat.m;
            let n_ext_cols = ext_mat.m;

            // 4. The commitment is just the root, but since each commitment could be to a differently-sized polynomial, we also add some metadata.
            let commitment = LigeroPCCommitment {
                metadata: Metadata {
                    n_rows,
                    n_cols,
                    n_ext_cols,
                },
                root,
            };

            commitments.push(LabeledCommitment::new(
                labeled_polynomial.label().clone(),
                commitment,
                None,
            ));
        }
        let com_len = &commitments.len();
        Ok((commitments, vec![Self::Randomness::empty(); *com_len]))
    }

    fn open<'a>(
        ck: &Self::CommitterKey,
        labeled_polynomials: impl IntoIterator<Item = &'a LabeledPolynomial<F, P>>,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Self::Commitment>>,
        point: &'a P::Point,
        _challenge_generator: &mut crate::challenge::ChallengeGenerator<F, S>,
        _rands: impl IntoIterator<Item = &'a Self::Randomness>,
        _rng: Option<&mut dyn RngCore>,
    ) -> Result<Self::Proof, Self::Error>
    where
        P: 'a,
        Self::Randomness: 'a,
        Self::Commitment: 'a,
    {
        let mut proof_array = LPCPArray::default();
        let labeled_commitments: Vec<&'a LabeledCommitment<Self::Commitment>> =
            commitments.into_iter().collect();
        let labeled_polynomials: Vec<&'a LabeledPolynomial<F, P>> =
            labeled_polynomials.into_iter().collect();

        if labeled_commitments.len() != labeled_polynomials.len() {
            return Err(Error::IncorrectInputLength(format!(
                "Mismatched lengths: {} commitments, {} polynomials",
                labeled_commitments.len(),
                labeled_polynomials.len()
            )));
        }

        for i in 0..labeled_polynomials.len() {
            let polynomial = labeled_polynomials[i].polynomial();
            let commitment = labeled_commitments[i].commitment();
            let n_rows = commitment.metadata.n_rows;
            let n_cols = commitment.metadata.n_cols;
            let root = &commitment.root;

            // 1. Arrange the coefficients of the polynomial into a matrix,
            // and apply Reed-Solomon encoding to get `ext_mat`.
            let (mat, ext_mat) = Self::compute_matrices(polynomial, ck.rho_inv);

            // 2. Create the Merkle tree from the hashes of each column.
            let col_tree =
                Self::create_merkle_tree(&ext_mat, &ck.leaf_hash_params, &ck.two_to_one_params);

            // 3. Generate vector `b = [1, z^m, z^(2m), ..., z^((m-1)m)]`
            let mut b = Vec::new();
            // This could potentially fail when n_cols > 1<<64, but `ck` won't allow commiting to such polynomials.
            let point_pow = point.pow([n_cols as u64]);
            let mut pow_b = F::one();
            for _ in 0..n_rows {
                b.push(pow_b);
                pow_b *= point_pow;
            }

            let mut transcript = IOPTranscript::new(b"transcript");
            transcript
                .append_serializable_element(b"root", root)
                .map_err(|_| Error::TranscriptError)?;

            // If we are checking well-formedness, we need to compute the well-formedness proof (which is just r.M) and append it to the transcript.
            let well_formedness = if ck.check_well_formedness {
                let mut r = Vec::new();
                for _ in 0..n_rows {
                    r.push(
                        transcript
                            .get_and_append_challenge(b"r")
                            .map_err(|_| Error::TranscriptError)?,
                    );
                }
                let v = mat.row_mul(&r);

                transcript
                    .append_serializable_element(b"v", &v)
                    .map_err(|_| Error::TranscriptError)?;
                Some(v)
            } else {
                None
            };

            transcript
                .append_serializable_element(b"point", point)
                .map_err(|_| Error::TranscriptError)?;

            proof_array.push(LigeroPCProof {
                // compute the opening proof and append b.M to the transcript
                opening: Self::generate_proof(
                    ck.sec_param,
                    ck.rho_inv,
                    &b,
                    &mat,
                    &ext_mat,
                    &col_tree,
                    &mut transcript,
                )?,
                well_formedness,
            });
        }

        Ok(proof_array)
    }

    fn check<'a>(
        vk: &Self::VerifierKey,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Self::Commitment>>,
        point: &'a P::Point,
        values: impl IntoIterator<Item = F>,
        proof_array: &Self::Proof,
        _challenge_generator: &mut crate::challenge::ChallengeGenerator<F, S>,
        _rng: Option<&mut dyn RngCore>,
    ) -> Result<bool, Self::Error>
    where
        Self::Commitment: 'a,
    {
        let labeled_commitments: Vec<&'a LabeledCommitment<Self::Commitment>> =
            commitments.into_iter().collect();
        let values: Vec<F> = values.into_iter().collect();

        if labeled_commitments.len() != proof_array.len()
            || labeled_commitments.len() != values.len()
        {
            return Err(Error::IncorrectInputLength(
                format!(
                    "Mismatched lengths: {} proofs were provided for {} commitments with {} claimed values",labeled_commitments.len(), proof_array.len(), values.len()
                )
            ));
        }
        let leaf_hash_params: &<<C as Config>::LeafHash as CRHScheme>::Parameters =
            &vk.leaf_hash_params;
        let two_to_one_params: &<<C as Config>::TwoToOneHash as TwoToOneCRHScheme>::Parameters =
            &vk.two_to_one_params;

        for (i, labeled_commitment) in labeled_commitments.iter().enumerate() {
            let commitment = labeled_commitment.commitment();
            let n_rows = commitment.metadata.n_rows;
            let n_cols = commitment.metadata.n_cols;
            let n_ext_cols = commitment.metadata.n_ext_cols;
            let root = &commitment.root;
            let t = calculate_t::<F>(vk.sec_param, vk.rho_inv, n_ext_cols)?;

            let mut transcript = IOPTranscript::new(b"transcript");
            transcript
                .append_serializable_element(b"root", &commitment.root)
                .map_err(|_| Error::TranscriptError)?;

            let out = if vk.check_well_formedness {
                if proof_array[i].well_formedness.is_none() {
                    return Err(Error::InvalidCommitment);
                }
                let tmp = &proof_array[i].well_formedness.as_ref();
                let well_formedness = tmp.unwrap();
                let mut r = Vec::with_capacity(n_rows);
                for _ in 0..n_rows {
                    r.push(
                        transcript
                            .get_and_append_challenge(b"r")
                            .map_err(|_| Error::TranscriptError)?,
                    );
                }
                // Upon sending `v` to the Verifier, add it to the sponge. Claim is that v = r.M
                transcript
                    .append_serializable_element(b"v", well_formedness)
                    .map_err(|_| Error::TranscriptError)?;

                (Some(well_formedness), Some(r))
            } else {
                (None, None)
            };

            // 1. Seed the transcript with the point and the recieved vector
            // TODO Consider removing the evaluation point from the transcript.
            transcript
                .append_serializable_element(b"point", point)
                .map_err(|_| Error::TranscriptError)?;
            transcript
                .append_serializable_element(b"v", &proof_array[i].opening.v)
                .map_err(|_| Error::TranscriptError)?;

            // 2. Ask random oracle for the `t` indices where the checks happen
            let indices = get_indices_from_transcript::<F>(n_ext_cols, t, &mut transcript)?;

            // 3. Hash the received columns into leaf hashes
            let col_hashes: Vec<_> = proof_array[i]
                .opening
                .columns
                .iter()
                .map(|c| hash_column::<D, F>(c))
                .collect();

            // 4. Verify the paths for each of the leaf hashes - this is only run once,
            // even if we have a well-formedness check (i.e., we save sending and checking the columns).
            // See "Concrete optimizations to the commitment scheme", p.12 of [Brakedown](https://eprint.iacr.org/2021/1043.pdf)
            for (j, (leaf, q_j)) in col_hashes.iter().zip(indices.iter()).enumerate() {
                let path = &proof_array[i].opening.paths[j];
                if path.leaf_index != *q_j {
                    return Err(Error::InvalidCommitment);
                }

                path.verify(leaf_hash_params, two_to_one_params, root, leaf.clone())
                    .map_err(|_| Error::InvalidCommitment)?;
            }

            // helper closure: checks if a.b = c
            let check_inner_product = |a, b, c| -> Result<(), Error> {
                if inner_product(a, b) != c {
                    return Err(Error::InvalidCommitment);
                }

                Ok(())
            };

            // 5. Compute the encoding w = E(v)
            let w = reed_solomon(&proof_array[i].opening.v, vk.rho_inv);

            // 6. Compute a = [1, z, z^2, ..., z^(n_cols_1)]
            // where z denotes the query `point`.
            let mut a = Vec::with_capacity(n_cols);
            let mut pow_a = F::one();
            for _ in 0..n_cols {
                a.push(pow_a);
                pow_a *= point;
            }

            // Compute b = [1, z^n_cols, z^(2*n_cols), ..., z^((n_rows-1)*n_cols)]
            let mut b = Vec::with_capacity(n_rows);
            let mut pow_b = F::one();
            for _ in 0..n_rows {
                b.push(pow_b);
                pow_b *= pow_a;
            }
            let coeffs: &[F] = &b;

            // 7. Probabilistic checks that whatever the prover sent,
            // matches with what the verifier computed for himself.
            // Note: we sacrifice some code repetition in order not to repeat execution.
            if let (Some(well_formedness), Some(r)) = out {
                let w_well_formedness = reed_solomon(well_formedness, vk.rho_inv);
                for (transcript_index, matrix_index) in indices.iter().enumerate() {
                    check_inner_product(
                        &r,
                        &proof_array[i].opening.columns[transcript_index],
                        w_well_formedness[*matrix_index],
                    )?;
                    check_inner_product(
                        coeffs,
                        &proof_array[i].opening.columns[transcript_index],
                        w[*matrix_index],
                    )?;
                }
            } else {
                for (transcript_index, matrix_index) in indices.iter().enumerate() {
                    check_inner_product(
                        coeffs,
                        &proof_array[i].opening.columns[transcript_index],
                        w[*matrix_index],
                    )?;
                }
            }

            if inner_product(&proof_array[i].opening.v, &a) != values[i] {
                println!("Function check: claimed value in position {i} does not match the evaluation of the committed polynomial in the same position");
                return Ok(false);
            }
        }

        Ok(true)
    }
}
