use crate::{PCCommitterKey, PCCommitment, PCVerifierKey};
use crate::{BTreeMap, BTreeSet, ToString, Vec};
use crate::{BatchLCProof, Error, Evaluations, QuerySet};
use crate::{LabeledCommitment, LabeledPolynomial, LinearCombination};
use crate::{PCRandomness, PCUniversalParams, Polynomial, PolynomialCommitment};

use algebra_core::{One, PairingEngine, ProjectiveCurve, UniformRand, Zero, VariableBaseMSM, FixedBaseMSM, PrimeField, ToBytes, Field, AffineCurve};
use core::{convert::TryInto, marker::PhantomData};
use rand_core::RngCore;

mod data_structures;
pub use data_structures::*;

#[cfg(feature = "parallel")]
use rayon::prelude::*;

use std::iter::Iterator;
use digest::Digest;

pub struct InnerProductArg<G: ProjectiveCurve, D: Digest> {
    _projective: PhantomData<G>,
    _digest: PhantomData<D>
}

impl<G: ProjectiveCurve, D: Digest> InnerProductArg<G, D> {
    fn dh_commit (
        comm_key: &[G::Affine],
        scalars: &[G::ScalarField]
    ) -> G {
        let scalars_bigint = ff_fft::cfg_iter!(scalars)
            .map(|s| s.into_repr())
            .collect::<Vec<<G::ScalarField as PrimeField>::BigInt>>();

        VariableBaseMSM::multi_scalar_mul(
            comm_key,
            scalars_bigint.as_slice(),
        )
    }

    // Takes three ToByte elements and converts them into a scalar field element
    fn rand_oracle (
        a: impl ToBytes,
        b: impl ToBytes,
        c: impl ToBytes
    ) -> G::ScalarField {
        let mut hash_input = Vec::new();

        hash_input.extend_from_slice (algebra_core::to_bytes!(a).unwrap().as_slice());
        hash_input.extend_from_slice (algebra_core::to_bytes!(b).unwrap().as_slice());
        hash_input.extend_from_slice (algebra_core::to_bytes!(c).unwrap().as_slice());
        hash_input.push(0u8);

        let mut point = None;
        while point.is_none() {
            let n = hash_input.len();
            if hash_input[n-1] == 255u8 {
                hash_input.push(1u8);
            } else {
                hash_input[n-1] += 1u8;
            }

            let hash = D::digest(hash_input.as_slice());
            point = <G::ScalarField as Field>::from_random_bytes (&hash);
        }

        point.unwrap()
    }

    fn inner_product(
        l: &[G::ScalarField],
        r: &[G::ScalarField]
    ) -> G::ScalarField {
        ff_fft::cfg_iter!(l)
            .zip(r)
            .map(|(li, ri)| *li * ri)
            .sum()
    }

    fn rep_poly (
        log_d: usize,
        challenges: &Vec<G::ScalarField>
    ) -> Polynomial<G::ScalarField> {
        // Representative polynomial is degree 2^log_d
        let mut coeffs = vec![G::ScalarField::one(); 1 << log_d];

        for (i, challenge) in challenges.iter().enumerate() {
            let i = i + 1;
            let elem_degree = 1 << (log_d - i);
            for start in (elem_degree..coeffs.len()).step_by(elem_degree * 2) {
                for offset in 0..elem_degree {
                    coeffs[start + offset] *= challenge;
                }
            }
        }

        Polynomial::from_coefficients_vec(coeffs)
    }
}

impl<G: ProjectiveCurve, D: Digest> PolynomialCommitment<G::ScalarField> for InnerProductArg<G, D> {
    type UniversalParams = UniversalParams<G>;
    type CommitterKey = CommitterKey<G>;
    type VerifierKey = VerifierKey<G>;
    type Commitment = Commitment<G>;
    type Randomness = Randomness;
    type Proof = Proof<G>;
    type BatchProof = Vec<Self::Proof>;
    type Error = Error;

    fn setup<R: RngCore>(
        max_degree: usize,
        rng: &mut R,
    ) -> Result<Self::UniversalParams, Self::Error> {
        // TODO: Should we error if max_degree+1 isn't a power of 2 or should we just adjust it be fit the specification like this:
        let max_degree = (1 << ((max_degree + 1) as f32).log2().ceil() as usize) - 1;

        /*
        // TODO: This doesn't work
        let mut affines = Vec::new();
        let mut hash_input = Vec::new();

        // Some random 32-bit seed
        hash_input.extend_from_slice (&[0xb, 0x9, 0xc, 0xd, 0x3, 0x0, 0xf, 0x7]);
        hash_input.push(0u8);

        while affines.len() < max_degree + 2 {
            let n = hash_input.len();
            if hash_input[n-1] == 255u8 {
                hash_input.push(1u8);
            } else {
                hash_input[n-1] += 1u8;
            }

            let hash = D::digest(hash_input.as_slice());
            let affine = G::Affine::from_random_bytes (&hash);
            if affine.is_some() {
                affines.push (affine.unwrap());
            }
        }

        let h: G = affines.pop().unwrap().into();

        let pp = UniversalParams {
            comm_key: affines,
            h
        };

        */

        // TODO: But this does work?
        let mut projectives: Vec<G> = Vec::new();
        while projectives.len() < max_degree + 2 {
            projectives.push(G::rand(rng));
        }

        let h: G = projectives.pop().unwrap();
        let comm_key = G::batch_normalization_into_affine(projectives.as_slice());

        let pp = UniversalParams {
            comm_key,
            h
        };

        Ok(pp)
    }

    fn trim(
        pp: &Self::UniversalParams,
        supported_degree: usize,
        _enforced_degree_bounds: Option<&[usize]>,
    ) -> Result<(Self::CommitterKey, Self::VerifierKey), Self::Error> {
        //TODO: Same question as above with supported degree
        let supported_degree = (1 << ((supported_degree + 1) as f32).log2().ceil() as usize) - 1;
        if supported_degree > pp.max_degree () {
            return Err(Error::TrimmingDegreeTooLarge);
        }

        let ck = CommitterKey {
            comm_key: pp.comm_key[0..(supported_degree + 1)].to_vec(),
            h: pp.h.clone(),
            max_degree: pp.max_degree()
        };

        let vk = VerifierKey {
            comm_key: pp.comm_key[0..(supported_degree + 1)].to_vec(),
            h: pp.h.clone(),
            max_degree: pp.max_degree()
        };

        Ok((ck, vk))
    }

    /// Outputs a commitment to `polynomial`.
    fn commit<'a>(
        ck: &Self::CommitterKey,
        polynomials: impl IntoIterator<Item = &'a LabeledPolynomial<'a, G::ScalarField>>,
        _rng: Option<&mut dyn RngCore>,
    ) -> Result<
        (
            Vec<LabeledCommitment<Self::Commitment>>,
            Vec<Self::Randomness>,
        ),
        Self::Error,
    > {
        let mut comms: Vec<LabeledCommitment<Commitment<G>>> = Vec::new();
        for labeled_polynomial in polynomials {
            let polynomial = labeled_polynomial.polynomial();
            let label = labeled_polynomial.label();
            let degree_bound = labeled_polynomial.degree_bound();

            // TODO: Do I need to check degree bounds? We never use degree bounds in this scheme.
            if polynomial.degree() < 1 {
                return Err(Error::DegreeIsZero);
            } else if polynomial.degree() > ck.supported_degree() {
                // TODO: Should this be supported_degree or number of keys (supported_degree + 1)? In a previous use case, it was number of keys, but that doesn't make too much sense because num_coefficients was in terms of degree.
                return Err(Error::TooManyCoefficients {
                    num_coefficients: polynomial.degree(),
                    num_powers: ck.supported_degree()
                });
            }

            let comm = Self::dh_commit(&ck.comm_key[..(polynomial.degree() + 1)], &polynomial.coeffs);
            let labeled_comm = LabeledCommitment::new(
                label.to_string(),
                Commitment(comm),
                degree_bound
            );

            comms.push(labeled_comm);
        }

        let num_comms = comms.len();
        Ok((comms, vec![Randomness(); num_comms]))
    }

    fn open<'a>(
        ck: &Self::CommitterKey,
        labeled_polynomials: impl IntoIterator<Item = &'a LabeledPolynomial<'a, G::ScalarField>>,
        point: G::ScalarField,
        opening_challenge: G::ScalarField,
        _rands: impl IntoIterator<Item = &'a Self::Randomness>,
        commitments: Option<impl IntoIterator<Item = &'a LabeledCommitment<Self::Commitment>>>,
    ) -> Result<Self::Proof, Self::Error>
    where
        Self::Commitment: 'a,
        Self::Randomness: 'a,
    {
        let commitments = commitments.ok_or(Error::MissingCommitment)?;

        let mut combined_polynomial = Polynomial::zero();
        let mut combined_commitment = Commitment::empty();

        let labeled_polynomials = labeled_polynomials.into_iter();
        let labeled_commitments = commitments.into_iter();

        let mut curr_challenge = opening_challenge;
        for (labeled_polynomial, labeled_commitment) in labeled_polynomials.zip(labeled_commitments) {
            let polynomial = labeled_polynomial.polynomial();

            // TODO: Check degree bounds too?
            if polynomial.degree() < 1 {
                return Err(Error::DegreeIsZero);
            } else if polynomial.degree() > ck.supported_degree() {
                // TODO: Same question as before with num_coefficients and num_powers
                return Err(Error::TooManyCoefficients {
                    num_coefficients: polynomial.degree(),
                    num_powers: ck.supported_degree()
                });
            }

            combined_polynomial += (curr_challenge, labeled_polynomial.polynomial());
            combined_commitment += (curr_challenge, labeled_commitment.commitment());
            curr_challenge *= &opening_challenge;
        }

        let combined_v = combined_polynomial.evaluate(point);

        // ith challenge
        let mut round_challenge = Self::rand_oracle(combined_commitment, point, combined_v);
        let h_prime = ck.h.mul(round_challenge);

        // Pad the coefficients to the appropriate vector size
        let d = ck.supported_degree();
        let mut coeffs = combined_polynomial.coeffs;
        if coeffs.len() < d + 1 {
            for _ in coeffs.len()..(d+1) {
                coeffs.push(G::ScalarField::zero());
            }
        }
        let mut coeffs = coeffs.as_mut_slice();

        // Powers of z
        let mut z: Vec<G::ScalarField> = Vec::with_capacity(d + 1);
        let mut curr_z: G::ScalarField = G::ScalarField::one();
        for _ in 0..(d+1) {
            z.push(curr_z);
            curr_z *= &point;
        }
        let mut z = z.as_mut_slice();

        // Comm key in both affine and projective form
        let mut key = Vec::new();
        let mut key_proj =
            ck.comm_key.iter()
                .map(|x| (*x).into())
                .collect::<Vec<G>>();
        let mut key_proj = key_proj.as_mut_slice();

        let mut l_vec = Vec::new();
        let mut r_vec = Vec::new();

        let mut n = d + 1;
        while n > 1 {
            let (coeffs_l, coeffs_r) = coeffs.split_at_mut(n/2);
            let (z_l, z_r) = z.split_at_mut(n/2);

            let (key_l, key_r) =
                if n != d+1 {
                    key.split_at(n/2)
                } else {
                    ck.comm_key.split_at(n/2)
                };

            let mut l = Self::dh_commit(key_l, coeffs_r);
            l += h_prime.mul(Self::inner_product(coeffs_r, z_l));

            let mut r = Self::dh_commit(key_r, coeffs_l);
            r += h_prime.mul(Self::inner_product(coeffs_l, z_r));

            l_vec.push(l);
            r_vec.push(r);

            round_challenge = Self::rand_oracle (round_challenge, l, r);
            let (key_proj_l, key_proj_r) = key_proj.split_at_mut(n/2);

            // TODO: When can unwrap fail? Is it safe to just unwrap the inverse?
            for i in 0..n/2 {
                coeffs_l[i] += &(round_challenge.inverse().unwrap() * &coeffs_r[i]);
                z_l[i] += &(round_challenge * &z_r[i]);
                key_proj_l[i] += &(key_proj_r[i].mul(round_challenge));
            }

            n /= 2;
            coeffs = coeffs_l;
            z = z_l;
            key_proj = key_proj_l;
            key = G::batch_normalization_into_affine(key_proj);
        }

        Ok(Proof {
            l_vec,
            r_vec,
            final_comm_key: key[0].clone(),
            c: coeffs[0].clone()
        })
    }

    fn check<'a>(
        vk: &Self::VerifierKey,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Self::Commitment>>,
        point: G::ScalarField,
        values: impl IntoIterator<Item = G::ScalarField>,
        proof: &Self::Proof,
        opening_challenge: G::ScalarField,
    ) -> Result<bool, Self::Error>
    where
        Self::Commitment: 'a,
    {
        let d = vk.supported_degree();
        let log_d = ((d + 1) as f32).log2() as usize;

        //TODO: Do we need to check that the number of commitments is equal to the number of values? Similarly, should we check that the proof size is log_(d+1)
        //TODO: In the other implementations, we never explicitly check that the number of commitments or values are the same

        let mut combined_commitment = Commitment::empty();
        let mut combined_v = G::ScalarField::zero();

        let mut curr_challenge = opening_challenge;
        let labeled_commitment = commitments.into_iter();
        let values = values.into_iter();
        for (labeled_commitment, value) in labeled_commitment.zip(values) {
            combined_commitment += (curr_challenge, labeled_commitment.commitment());
            combined_v += &(curr_challenge * &value);
            curr_challenge *= &opening_challenge;
        }

        // Challenge for each round
        let mut round_challenges = Vec::new();
        let mut round_challenge = Self::rand_oracle(combined_commitment, point, combined_v);

        let h_prime = vk.h.mul(round_challenge);

        let mut round_commitment = combined_commitment + h_prime.mul(combined_v);

        let l_iter = proof.l_vec.iter();
        let r_iter = proof.r_vec.iter();

        for (l, r) in l_iter.zip(r_iter) {
            round_challenge = Self::rand_oracle(round_challenge, l, r);
            round_challenges.push(round_challenge);
            // TODO: When can unwrap fail? Is it safe to just unwrap the inverse?
            round_commitment += l.mul(round_challenge.inverse().unwrap()) + r.mul(round_challenge);
        }

        let rep_poly = Self::rep_poly(log_d, &round_challenges);
        let v_prime = rep_poly.evaluate (point) * &proof.c;
        let h_prime = G::batch_normalization_into_affine(&[h_prime]).pop().unwrap();

        let check_commitment_elem: G =
            Self::dh_commit(
                &[proof.final_comm_key.clone(), h_prime],
                &[proof.c.clone(), v_prime]
            );
        if !(round_commitment.0 - check_commitment_elem).is_zero() {
            return Ok(false);
        }

        let final_key = Self::dh_commit(vk.comm_key.as_slice(), rep_poly.coeffs.as_slice());
        if !(final_key - proof.final_comm_key.into()).is_zero() {
            return Ok(false);
        }

        Ok(true)
    }

}

#[cfg(test)]
mod tests {
    #![allow(non_camel_case_types)]

    use crate::inner_product_arg::InnerProductArg;

    use algebra::jubjub;
    use blake2::Blake2s;

    // TODO: Add more curves
    type PC<E, D> = InnerProductArg<E, D>;
    type PC_JJB2S = PC <jubjub::JubJubProjective, Blake2s>;

    #[test]
    fn single_poly_test() {
        use crate::tests::*;
        single_poly_test::<_, PC_JJB2S>().expect("test failed for jubjub-blake2s");
    }

    #[test]
    fn quadratic_poly_degree_bound_multiple_queries_test() {
        use crate::tests::*;
        quadratic_poly_degree_bound_multiple_queries_test::<_, PC_JJB2S>()
            .expect("test failed for jubjub-blake2s");
    }

    #[test]
    fn linear_poly_degree_bound_test() {
        use crate::tests::*;
        linear_poly_degree_bound_test::<_, PC_JJB2S>().expect("test failed for jubjub-blake2s");
    }

    #[test]
    fn single_poly_degree_bound_test() {
        use crate::tests::*;
        single_poly_degree_bound_test::<_, PC_JJB2S>().expect("test failed for jubjub-blake2s");
    }

    #[test]
    fn single_poly_degree_bound_multiple_queries_test() {
        use crate::tests::*;
        single_poly_degree_bound_multiple_queries_test::<_, PC_JJB2S>()
            .expect("test failed for jubjub-blake2s");
    }

    #[test]
    fn two_polys_degree_bound_single_query_test() {
        use crate::tests::*;
        two_polys_degree_bound_single_query_test::<_, PC_JJB2S>()
            .expect("test failed for jubjub-blake2s");
    }

    #[test]
    fn full_end_to_end_test() {
        use crate::tests::*;
        full_end_to_end_test::<_, PC_JJB2S>().expect("test failed for jubjub-blake2s");
        println!("Finished jubjub-blake2s");
    }

    #[test]
    fn single_equation_test() {
        use crate::tests::*;
        single_equation_test::<_, PC_JJB2S>().expect("test failed for jubjub-blake2s");
        println!("Finished jubjub-blake2s");
    }

    #[test]
    fn two_equation_test() {
        use crate::tests::*;
        two_equation_test::<_, PC_JJB2S>().expect("test failed for jubjub-blake2s");
        println!("Finished jubjub-blake2s");
    }

    #[test]
    fn two_equation_degree_bound_test() {
        use crate::tests::*;
        two_equation_degree_bound_test::<_, PC_JJB2S>().expect("test failed for jubjub-blake2s");
        println!("Finished jubjub-blake2s");
    }

    #[test]
    fn full_end_to_end_equation_test() {
        use crate::tests::*;
        full_end_to_end_equation_test::<_, PC_JJB2S>().expect("test failed for jubjub-blake2s");
        println!("Finished jubjub-blake2s");
    }
}
