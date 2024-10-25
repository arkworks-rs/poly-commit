use crate::hyrax::utils::tensor_prime;
use crate::utils::{inner_product, scalar_by_vector, vector_sum, Matrix};
use crate::{
    hyrax::utils::flat_to_matrix_column_major, Error, LabeledCommitment, LabeledPolynomial,
    PolynomialCommitment,
};
use ark_crypto_primitives::sponge::{Absorb, CryptographicSponge};
use ark_ec::{AffineRepr, CurveGroup, VariableBaseMSM};
use ark_ff::PrimeField;
use ark_poly::MultilinearExtension;
use ark_serialize::serialize_to_vec;
use ark_std::{marker::PhantomData, rand::RngCore, string::ToString, vec::Vec, UniformRand};

use blake2::Blake2s256;
use digest::Digest;

#[cfg(feature = "parallel")]
use rayon::prelude::*;

mod data_structures;
pub use data_structures::*;
#[cfg(test)]
mod tests;
mod utils;
/// String of bytes used to seed the randomness during the setup function.
/// Note that the latter should never be used in production environments.
pub const PROTOCOL_NAME: &'static [u8] = b"Hyrax protocol";

/// Hyrax polynomial committment scheme:
/// A polynomial commitment scheme based on the hardness of the
/// discrete logarithm problem in prime-order groups. This is a
/// Fiat-Shamired version of the PCS described in the Hyrax paper
/// [[WTsTW17]][hyrax].
///
/// [hyrax]: https://eprint.iacr.org/2017/1132.pdf
///
/// ### Future optimisations
///
/// - Add parallelisation. There is at least one natural place where
///   parallelisation could bring performance gains: in essence, the prover
///   commits to the polynomial by expressing it as an evaluation matrix and
///   Pederson-multi-committing to each row. Each of this commitments can be
///   computed independently from the rest, and therefore, in parallel. It is
///   still to be seen how much of an improvement this would entail, since each
///   Pederson multi-commitment boils down to a multi-exponentiation and this
///   operation is itself parallelised.
/// - Due to the homomorphic nature of Pedersen commitments, it is likely
///   some of the following methods can be designed more efficiently than their
///   default implementations: `batch_open`, `batch_check`,
///   `open_combinations`, `check_combinations`. This is not discussed in the
///   reference article, but the IPA and KZG modules might be a good starting
///   point.
/// - On a related note to the previous point, there might be a more
///   efficient way to open several polynomials at a single point (this is the
///   functionality of the `open` method) than the currently implemented
///   technique, where only the computation of the vectors `L` and `R` is
///   shared across polynomials.
/// - The cited article proposes an optimisation in the section _Reducing the
///   cost of proof-of-dot-prod_. It allows for non-square matrices (and hence
///   removes the requirement for the number of variables to be even) and
///   introduces a tradeoff between proof size and verifier time. It is
///   probably worth pursuing.

pub struct HyraxPC<
    // The elliptic curve used for Pedersen commitments (only EC groups are
    // supported as of now).
    G: AffineRepr,
    // A polynomial type representing multilinear polynomials
    P: MultilinearExtension<G::ScalarField>,
    // The sponge used in the protocol as random oracle
    S: CryptographicSponge,
> {
    _phantom: PhantomData<(G, P, S)>,
}

impl<G, P, S> HyraxPC<G, P, S>
where
    G: AffineRepr,
    P: MultilinearExtension<G::ScalarField>,
    S: CryptographicSponge,
{
    /// Pedersen commitment to a vector of scalars as described in appendix A.1
    /// of the reference article.
    /// The function does not add handle hiding term `h * r`.
    /// It is only a wrapper around MSM.
    ///
    /// # Panics
    ///
    /// Panics if `key` and `scalars` do not have the same length
    fn pedersen_commit(key: &[G], scalars: &[G::ScalarField]) -> G::Group {
        assert_eq!(key.len(), scalars.len());
        let scalars_bigint = ark_std::cfg_iter!(scalars)
            .map(|s| s.into_bigint())
            .collect::<Vec<_>>();
        // Multi-exponentiation in the group of points of the EC
        <G::Group as VariableBaseMSM>::msm_bigint(&key, &scalars_bigint)
    }
}

impl<G, P, S> PolynomialCommitment<G::ScalarField, P, S> for HyraxPC<G, P, S>
where
    G: AffineRepr,
    G::ScalarField: Absorb,
    P: MultilinearExtension<G::ScalarField>,
    S: CryptographicSponge,
{
    type UniversalParams = HyraxUniversalParams<G>;
    type CommitterKey = HyraxCommitterKey<G>;
    type VerifierKey = HyraxVerifierKey<G>;
    type Commitment = HyraxCommitment<G>;
    type CommitmentState = HyraxCommitmentState<G::ScalarField>;
    type Proof = Vec<HyraxProof<G>>;
    type BatchProof = Vec<Self::Proof>;
    type Error = Error;

    /// Outputs mock universal parameters for the Hyrax polynomial commitment
    /// scheme. It does *not* return random keys across calls and should never
    /// be used in settings where security is required - it is only useful for
    /// testing.
    ///
    /// # Panics
    ///
    /// Panics if `num_vars` is None or contains an odd value.
    fn setup<R: RngCore>(
        _max_degree: usize,
        num_vars: Option<usize>,
        _rng: &mut R,
    ) -> Result<Self::UniversalParams, Self::Error> {
        if num_vars.is_none() {
            return Err(Error::InvalidNumberOfVariables);
        }

        let n = num_vars.unwrap();

        if n % 2 == 1 {
            // Only polynomials with an even number of variables are
            // supported in this implementation
            return Err(Error::InvalidNumberOfVariables);
        }

        // Number of rows (or, equivalently, colums) of a square matrix
        // containing the coefficients of an n-variate ML polynomial
        let dim = 1 << n / 2;

        // The following block of code is largely taking from the IPA module
        // in this crate. It generates random points (not guaranteed to be
        // generators, since the point at infinity should theoretically occur)
        let points: Vec<_> = ark_std::cfg_into_iter!(0u64..dim + 1)
            .map(|i| {
                let hash =
                    Blake2s256::digest([PROTOCOL_NAME, &i.to_le_bytes()].concat().as_slice());
                let mut p = G::from_random_bytes(&hash);
                let mut j = 0u64;
                while p.is_none() {
                    let mut bytes = PROTOCOL_NAME.to_vec();
                    bytes.extend(i.to_le_bytes());
                    bytes.extend(j.to_le_bytes());
                    let hash = Blake2s256::digest(bytes.as_slice());
                    p = G::from_random_bytes(&hash);
                    j += 1;
                }
                let point = p.unwrap();
                point.mul_by_cofactor_to_group()
            })
            .collect();

        // Converting from projective to affine representation
        let mut points = G::Group::normalize_batch(&points);

        let h: G = points.pop().unwrap();

        Ok(HyraxUniversalParams { com_key: points, h })
    }

    /// Trims a key into a prover key and a verifier key. This should only
    /// amount to discarding some of the points in said key if the prover
    /// and verifier only wish to commit to polynomials with fewer variables
    /// than the key can support. Since the number of variables is not
    /// considered in the prototype, this function currently simply clones the
    /// key.
    fn trim(
        pp: &Self::UniversalParams,
        _supported_degree: usize,
        _supported_hiding_bound: usize,
        _enforced_degree_bounds: Option<&[usize]>,
    ) -> Result<(Self::CommitterKey, Self::VerifierKey), Self::Error> {
        Ok((pp.clone(), pp.clone()))
    }

    /// Produces a list of commitments to the passed polynomials. Cf. the
    /// section "Square-root commitment scheme" from the reference article.
    ///
    /// # Panics
    ///
    /// Panics if `rng` is None, since Hyrax requires randomness in order to
    /// commit to a polynomial
    #[allow(unused_variables)]
    fn commit<'a>(
        ck: &Self::CommitterKey,
        polynomials: impl IntoIterator<Item = &'a LabeledPolynomial<G::ScalarField, P>>,
        rng: Option<&mut dyn RngCore>,
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
        let mut coms = Vec::new();
        let mut states = Vec::new();

        #[cfg(not(feature = "parallel"))]
        let rng_inner = rng.expect("Committing to polynomials requires a random generator");

        for l_poly in polynomials {
            let label = l_poly.label();
            let poly = l_poly.polynomial();

            let n = poly.num_vars();
            let dim = 1 << n / 2;

            if n % 2 == 1 {
                // Only polynomials with an even number of variables are
                // supported in this implementation
                return Err(Error::InvalidNumberOfVariables);
            }

            if n > ck.com_key.len() {
                return Err(Error::InvalidNumberOfVariables);
            }

            let m = flat_to_matrix_column_major(&poly.to_evaluations(), dim, dim);

            // Commiting to the matrix with one multi-commitment per row
            let (row_coms, com_rands): (Vec<_>, Vec<_>) = cfg_iter!(m)
                .map(|row| {
                    #[cfg(not(feature = "parallel"))]
                    let r = G::ScalarField::rand(rng_inner);
                    #[cfg(feature = "parallel")]
                    let r = G::ScalarField::rand(&mut rand::thread_rng());
                    let c = (Self::pedersen_commit(&ck.com_key, row) + ck.h * r).into();
                    (c, r)
                })
                .unzip();

            let com = HyraxCommitment { row_coms };
            let l_comm = LabeledCommitment::new(label.to_string(), com, Some(1));

            coms.push(l_comm);
            states.push(HyraxCommitmentState {
                randomness: com_rands,
                mat: Matrix::new_from_rows(m),
            });
        }

        Ok((coms, states))
    }

    /// Opens a list of polynomial commitments at a desired point. This
    /// requires the list of original polynomials (`labeled_polynomials`) as
    /// well as the random values using by the Pedersen multi-commits during
    /// the commitment phase (`randomness`). Cf. sections "Square-root
    /// commitment scheme" and appendix A.2 from the reference article.
    ///
    /// # Panics
    ///
    /// Panics if
    /// - `rng` is None, since Hyrax requires randomness in order to
    /// open the commitment to a polynomial.
    /// - The point doesn't have an even number of variables.
    /// - The labels of a commitment doesn't match that of the corresponding
    /// polynomial.
    /// - The number of variables of a polynomial doesn't match that of the
    /// point.
    fn open<'a>(
        ck: &Self::CommitterKey,
        labeled_polynomials: impl IntoIterator<Item = &'a LabeledPolynomial<G::ScalarField, P>>,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Self::Commitment>>,
        point: &'a P::Point,
        sponge: &mut S,
        states: impl IntoIterator<Item = &'a Self::CommitmentState>,
        rng: Option<&mut dyn RngCore>,
    ) -> Result<Self::Proof, Self::Error>
    where
        Self::Commitment: 'a,
        Self::CommitmentState: 'a,
        P: 'a,
    {
        let n = point.len();

        if n % 2 == 1 {
            // Only polynomials with an even number of variables are
            // supported in this implementation
            return Err(Error::InvalidNumberOfVariables);
        }

        let dim = 1 << n / 2;

        // Reversing the point is necessary because the MLE interface returns
        // evaluations in little-endian order
        let point_rev: Vec<G::ScalarField> = point.iter().rev().cloned().collect();

        let point_lower = &point_rev[n / 2..];
        let point_upper = &point_rev[..n / 2];

        // Deriving the tensors which result in the evaluation of the polynomial
        // when they are multiplied by the coefficient matrix.
        let l = tensor_prime(point_lower);
        let r = tensor_prime(point_upper);

        let mut proofs = Vec::new();

        let rng_inner = rng.expect("Opening polynomials requires randomness");

        for (l_poly, (l_com, state)) in labeled_polynomials
            .into_iter()
            .zip(commitments.into_iter().zip(states.into_iter()))
        {
            let label = l_poly.label();
            if label != l_com.label() {
                return Err(Error::MismatchedLabels {
                    commitment_label: l_com.label().to_string(),
                    polynomial_label: label.to_string(),
                });
            }

            let poly = l_poly.polynomial();
            let com = l_com.commitment();

            if poly.num_vars() != n {
                return Err(Error::MismatchedNumVars {
                    poly_nv: poly.num_vars(),
                    point_nv: n,
                });
            }

            // Absorbing public parameters
            sponge.absorb(&serialize_to_vec!(*ck).map_err(|_| Error::TranscriptError)?);

            // Absorbing the commitment to the polynomial
            sponge.absorb(&serialize_to_vec!(com.row_coms).map_err(|_| Error::TranscriptError)?);

            // Absorbing the point
            sponge.absorb(point);

            // Commiting to the matrix formed by the polynomial coefficients
            let t = &state.mat;

            let lt = t.row_mul(&l);

            // t_prime coincides witht he Pedersen commitment to lt with the
            // randomnes r_lt computed here
            let r_lt = cfg_iter!(l)
                .zip(&state.randomness)
                .map(|(l, r)| *l * r)
                .sum::<G::ScalarField>();

            let eval = inner_product(&lt, &r);

            // Singleton commit
            let (com_eval, r_eval) = {
                let r = G::ScalarField::rand(rng_inner);
                ((ck.com_key[0] * eval + ck.h * r).into(), r)
            };

            // ******** Dot product argument ********
            // Appendix A.2 in the reference article

            let d: Vec<G::ScalarField> =
                (0..dim).map(|_| G::ScalarField::rand(rng_inner)).collect();

            let b = inner_product(&r, &d);

            // Multi-commit
            let r_d = G::ScalarField::rand(rng_inner);
            let com_d = (Self::pedersen_commit(&ck.com_key, &d) + ck.h * r_d).into();

            // Singleton commit
            let r_b = G::ScalarField::rand(rng_inner);
            let com_b = (ck.com_key[0] * b + ck.h * r_b).into();

            // Absorbing the commitment to the evaluation
            sponge.absorb(&serialize_to_vec!(com_eval).map_err(|_| Error::TranscriptError)?);

            // Absorbing the two auxiliary commitments
            sponge.absorb(&serialize_to_vec!(com_d).map_err(|_| Error::TranscriptError)?);
            sponge.absorb(&serialize_to_vec!(com_b).map_err(|_| Error::TranscriptError)?);

            // Receive the random challenge c from the verifier, i.e. squeeze
            // it from the transcript.
            let c = sponge.squeeze_field_elements(1)[0];

            let z = vector_sum(&d, &scalar_by_vector(c, &lt));
            let z_d = c * r_lt + r_d;
            let z_b = c * r_eval + r_b;

            proofs.push(HyraxProof {
                com_eval,
                com_d,
                com_b,
                z,
                z_d,
                z_b,
            });
        }

        Ok(proofs)
    }

    /// Verifies a list of opening proofs and confirms the evaluation of the
    /// committed polynomials at the desired point.
    ///
    /// # Panics
    /// - If the point doesn't have an even number of variables.
    /// - If the length of a commitment does not correspond to the length of the
    /// point (specifically, commitment length should be 2^(point-length/2)).
    ///
    /// # Disregarded arguments
    /// - `rng`
    fn check<'a>(
        vk: &Self::VerifierKey,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Self::Commitment>>,
        point: &'a P::Point,
        _values: impl IntoIterator<Item = G::ScalarField>,
        proof: &Self::Proof,
        sponge: &mut S,
        _rng: Option<&mut dyn RngCore>,
    ) -> Result<bool, Self::Error>
    where
        Self::Commitment: 'a,
    {
        let n = point.len();

        if n % 2 == 1 {
            // Only polynomials with an even number of variables are
            // supported in this implementation
            return Err(Error::InvalidNumberOfVariables);
        }

        // Reversing the point is necessary because the MLE interface returns
        // evaluations in little-endian order
        let point_rev: Vec<G::ScalarField> = point.iter().rev().cloned().collect();

        let point_lower = &point_rev[n / 2..];
        let point_upper = &point_rev[..n / 2];

        // Deriving the tensors which result in the evaluation of the polynomial
        // when they are multiplied by the coefficient matrix.
        let l = tensor_prime(point_lower);
        let r = tensor_prime(point_upper);

        for (com, h_proof) in commitments.into_iter().zip(proof.iter()) {
            let row_coms = &com.commitment().row_coms;

            // extract each field from h_proof
            let HyraxProof {
                com_eval,
                com_d,
                com_b,
                z,
                z_d,
                z_b,
            } = h_proof;

            if row_coms.len() != 1 << n / 2 {
                return Err(Error::IncorrectCommitmentSize {
                    encountered: row_coms.len(),
                    expected: 1 << n / 2,
                });
            }

            // Absorbing public parameters
            sponge.absorb(&serialize_to_vec!(*vk).map_err(|_| Error::TranscriptError)?);

            // Absorbing the commitment to the polynomial
            sponge.absorb(&serialize_to_vec!(*row_coms).map_err(|_| Error::TranscriptError)?);

            // Absorbing the point
            sponge.absorb(point);

            // Absorbing the commitment to the evaluation
            sponge.absorb(&serialize_to_vec!(*com_eval).map_err(|_| Error::TranscriptError)?);

            // Absorbing the two auxiliary commitments
            sponge.absorb(&serialize_to_vec!(*com_d).map_err(|_| Error::TranscriptError)?);
            sponge.absorb(&serialize_to_vec!(*com_b).map_err(|_| Error::TranscriptError)?);

            // Receive the random challenge c from the verifier, i.e. squeeze
            // it from the transcript.
            let c: G::ScalarField = sponge.squeeze_field_elements(1)[0];

            // Second check from the paper (figure 6, equation (14))
            // Moved here for potential early return
            let com_dp = (vk.com_key[0] * inner_product(&r, z) + vk.h * z_b).into();
            if com_dp != (com_eval.mul(c) + com_b).into() {
                return Ok(false);
            }

            // Computing t_prime with a multi-exponentiation
            let l_bigint = cfg_iter!(l)
                .map(|chi| chi.into_bigint())
                .collect::<Vec<_>>();
            let t_prime: G = <G::Group as VariableBaseMSM>::msm_bigint(&row_coms, &l_bigint).into();

            // First check from the paper (figure 6, equation (13))
            let com_z_zd = (Self::pedersen_commit(&vk.com_key, z) + vk.h * z_d).into();
            if com_z_zd != (t_prime.mul(c) + com_d).into() {
                return Ok(false);
            }
        }

        Ok(true)
    }
}
