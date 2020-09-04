use crate::{BTreeMap, BTreeSet, String, ToString, Vec};
use crate::{BatchLCProof, Error, Evaluations, QuerySet};
use crate::{LabeledCommitment, LabeledPolynomial, LinearCombination};
use crate::{PCCommitterKey, PCRandomness, PCUniversalParams, Polynomial, PolynomialCommitment};

use algebra_core::{
    to_bytes, AffineCurve, Field, One, PrimeField, ProjectiveCurve, ToBytes, UniformRand,
    VariableBaseMSM, Zero,
};
use core::{convert::TryInto, marker::PhantomData};
use rand_core::RngCore;

mod data_structures;
pub use data_structures::*;

#[cfg(feature = "parallel")]
use rayon::prelude::*;

use digest::Digest;

/// A polynomial commitment scheme based on the hardness of the
/// discrete logarithm problem in prime-order groups.
/// The construction is described in detail in [[BCMS20]][pcdas].
///
/// Degree bound enforcement requires that (at least one of) the points at
/// which a committed polynomial is evaluated are from a distribution that is
/// random conditioned on the polynomial. This is because degree bound
/// enforcement relies on checking a polynomial identity at this point.
/// More formally, the points must be sampled from an admissible query sampler,
/// as detailed in [[CHMMVW20]][marlin].
///
/// [pcdas]: https://eprint.iacr.org/2020/499
/// [marlin]: https://eprint.iacr.org/2019/1047
pub struct InnerProductArgPC<G: AffineCurve, D: Digest> {
    _projective: PhantomData<G>,
    _digest: PhantomData<D>,
}

impl<G: AffineCurve, D: Digest> InnerProductArgPC<G, D> {
    /// `PROTOCOL_NAME` is used as a seed for the setup function.
    pub const PROTOCOL_NAME: &'static [u8] = b"PC-DL-2020";

    /// Create a Pedersen commitment to `scalars` using the commitment key `comm_key`.
    /// Optionally, randomize the commitment using `hiding_generator` and `randomizer`.
    fn cm_commit(
        comm_key: &[G],
        scalars: &[G::ScalarField],
        hiding_generator: Option<G>,
        randomizer: Option<G::ScalarField>,
    ) -> G::Projective {
        let scalars_bigint = ff_fft::cfg_iter!(scalars)
            .map(|s| s.into_repr())
            .collect::<Vec<_>>();

        let mut comm = VariableBaseMSM::multi_scalar_mul(comm_key, &scalars_bigint);

        if randomizer.is_some() {
            assert!(hiding_generator.is_some());
            comm += &hiding_generator.unwrap().mul(randomizer.unwrap());
        }

        comm
    }

    fn compute_random_oracle_challenge(bytes: &[u8]) -> G::ScalarField {
        let mut i = 0u64;
        let mut challenge = None;
        while challenge.is_none() {
            let hash_input = algebra_core::to_bytes![bytes, i].unwrap();
            let hash = D::digest(&hash_input);
            challenge = <G::ScalarField as Field>::from_random_bytes(&hash);

            i += 1;
        }

        challenge.unwrap()
    }

    #[inline]
    fn inner_product(l: &[G::ScalarField], r: &[G::ScalarField]) -> G::ScalarField {
        ff_fft::cfg_iter!(l).zip(r).map(|(li, ri)| *li * ri).sum()
    }

    /// The succinct portion of `PC::check`. This algorithm runs in time
    /// O(log d), where d is the degree of the committed polynomials.
    fn succinct_check<'a>(
        vk: &VerifierKey<G>,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Commitment<G>>>,
        point: G::ScalarField,
        values: impl IntoIterator<Item = G::ScalarField>,
        proof: &Proof<G>,
        opening_challenge: G::ScalarField,
    ) -> Option<SuccinctCheckPolynomial<G::ScalarField>> {
        let check_time = start_timer!(|| "Succinct checking");

        let d = vk.supported_degree();

        // `log_d` is ceil(log2 (d + 1)), which is the number of steps to compute all of the challenges
        let log_d = algebra_core::log2(d + 1) as usize;

        let mut combined_commitment_proj = G::Projective::zero();
        let mut combined_v = G::ScalarField::zero();

        let mut cur_challenge = opening_challenge;
        let labeled_commitments = commitments.into_iter();
        let values = values.into_iter();

        for (labeled_commitment, value) in labeled_commitments.zip(values) {
            let commitment = labeled_commitment.commitment();
            combined_v += &(cur_challenge * &value);
            combined_commitment_proj += &labeled_commitment.commitment().comm.mul(cur_challenge);
            cur_challenge *= &opening_challenge;

            let degree_bound = labeled_commitment.degree_bound();
            assert_eq!(degree_bound.is_some(), commitment.shifted_comm.is_some());

            if let Some(degree_bound) = degree_bound {
                let shift = point.pow([(vk.supported_degree() - degree_bound) as u64]);
                combined_v += &(cur_challenge * &value * &shift);
                combined_commitment_proj += &commitment.shifted_comm.unwrap().mul(cur_challenge);
            }

            cur_challenge *= &opening_challenge;
        }

        let mut combined_commitment = combined_commitment_proj.into_affine();

        assert_eq!(proof.hiding_comm.is_some(), proof.rand.is_some());
        if proof.hiding_comm.is_some() {
            let hiding_comm = proof.hiding_comm.unwrap();
            let rand = proof.rand.unwrap();

            let hiding_challenge = Self::compute_random_oracle_challenge(
                &algebra_core::to_bytes![combined_commitment, point, combined_v, hiding_comm]
                    .unwrap(),
            );
            combined_commitment_proj += &(hiding_comm.mul(hiding_challenge) - &vk.s.mul(rand));
            combined_commitment = combined_commitment_proj.into_affine();
        }

        // Challenge for each round
        let mut round_challenges = Vec::with_capacity(log_d);
        let mut round_challenge = Self::compute_random_oracle_challenge(
            &algebra_core::to_bytes![combined_commitment, point, combined_v].unwrap(),
        );

        let h_prime = vk.h.mul(round_challenge);

        let mut round_commitment_proj = combined_commitment_proj + &h_prime.mul(combined_v);

        let l_iter = proof.l_vec.iter();
        let r_iter = proof.r_vec.iter();

        for (l, r) in l_iter.zip(r_iter) {
            round_challenge = Self::compute_random_oracle_challenge(
                &algebra_core::to_bytes![round_challenge, l, r].unwrap(),
            );
            round_challenges.push(round_challenge);
            round_commitment_proj +=
                &(l.mul(round_challenge.inverse().unwrap()) + &r.mul(round_challenge));
        }

        let check_poly = SuccinctCheckPolynomial::<G::ScalarField>(round_challenges);
        let v_prime = check_poly.evaluate(point) * &proof.c;
        let h_prime = h_prime.into_affine();

        let check_commitment_elem: G::Projective = Self::cm_commit(
            &[proof.final_comm_key.clone(), h_prime],
            &[proof.c.clone(), v_prime],
            None,
            None,
        );

        if !(round_commitment_proj - &check_commitment_elem).is_zero() {
            end_timer!(check_time);
            return None;
        }

        end_timer!(check_time);
        Some(check_poly)
    }

    fn check_degrees_and_bounds(
        supported_degree: usize,
        p: &LabeledPolynomial<G::ScalarField>,
    ) -> Result<(), Error> {
        if p.degree() < 1 {
            return Err(Error::DegreeIsZero);
        } else if p.degree() > supported_degree {
            return Err(Error::TooManyCoefficients {
                num_coefficients: p.degree() + 1,
                num_powers: supported_degree + 1,
            });
        }

        if let Some(bound) = p.degree_bound() {
            if bound < p.degree() || bound > supported_degree {
                return Err(Error::IncorrectDegreeBound {
                    poly_degree: p.degree(),
                    degree_bound: bound,
                    supported_degree,
                    label: p.label().to_string(),
                });
            }
        }

        Ok(())
    }

    fn shift_polynomial(
        ck: &CommitterKey<G>,
        p: &Polynomial<G::ScalarField>,
        degree_bound: usize,
    ) -> Polynomial<G::ScalarField> {
        if p.is_zero() {
            Polynomial::zero()
        } else {
            let mut shifted_polynomial_coeffs =
                vec![G::ScalarField::zero(); ck.supported_degree() - degree_bound];
            shifted_polynomial_coeffs.extend_from_slice(&p.coeffs);
            Polynomial::from_coefficients_vec(shifted_polynomial_coeffs)
        }
    }

    fn combine_shifted_rand(
        combined_rand: Option<G::ScalarField>,
        new_rand: Option<G::ScalarField>,
        coeff: G::ScalarField,
    ) -> Option<G::ScalarField> {
        if let Some(new_rand) = new_rand {
            let coeff_new_rand = new_rand * &coeff;
            return Some(combined_rand.map_or(coeff_new_rand, |r| r + &coeff_new_rand));
        };

        combined_rand
    }

    fn combine_shifted_comm(
        combined_comm: Option<G::Projective>,
        new_comm: Option<G>,
        coeff: G::ScalarField,
    ) -> Option<G::Projective> {
        if let Some(new_comm) = new_comm {
            let coeff_new_comm = new_comm.mul(coeff);
            return Some(combined_comm.map_or(coeff_new_comm, |c| c + &coeff_new_comm));
        };

        combined_comm
    }

    fn construct_labeled_commitments(
        lc_info: &[(String, Option<usize>)],
        elements: &[G::Projective],
    ) -> Vec<LabeledCommitment<Commitment<G>>> {
        let comms = G::Projective::batch_normalization_into_affine(elements);
        let mut commitments = Vec::new();

        let mut i = 0;
        for info in lc_info.into_iter() {
            let commitment;
            let label = info.0.clone();
            let degree_bound = info.1;

            if degree_bound.is_some() {
                commitment = Commitment {
                    comm: comms[i].clone(),
                    shifted_comm: Some(comms[i + 1].clone()),
                };

                i += 2;
            } else {
                commitment = Commitment {
                    comm: comms[i].clone(),
                    shifted_comm: None,
                };

                i += 1;
            }

            commitments.push(LabeledCommitment::new(label, commitment, degree_bound));
        }

        return commitments;
    }

    fn sample_generators(num_generators: usize) -> Vec<G> {
        let generators: Vec<_> = ff_fft::cfg_into_iter!(0..num_generators)
            .map(|i| {
                let i = i as u64;
                let mut hash = D::digest(&to_bytes![&Self::PROTOCOL_NAME, i].unwrap());
                let mut g = G::from_random_bytes(&hash);
                let mut j = 0u64;
                while g.is_none() {
                    hash = D::digest(&to_bytes![&Self::PROTOCOL_NAME, i, j].unwrap());
                    g = G::from_random_bytes(&hash);
                    j += 1;
                }
                let generator = g.unwrap();
                generator.mul_by_cofactor_to_projective()
            })
            .collect();

        G::Projective::batch_normalization_into_affine(&generators)
    }
}

impl<G: AffineCurve, D: Digest> PolynomialCommitment<G::ScalarField> for InnerProductArgPC<G, D> {
    type UniversalParams = UniversalParams<G>;
    type CommitterKey = CommitterKey<G>;
    type VerifierKey = VerifierKey<G>;
    type Commitment = Commitment<G>;
    type Randomness = Randomness<G>;
    type Proof = Proof<G>;
    type BatchProof = Vec<Self::Proof>;
    type Error = Error;

    fn setup<R: RngCore>(
        max_degree: usize,
        _rng: &mut R,
    ) -> Result<Self::UniversalParams, Self::Error> {
        // Ensure that max_degree + 1 is a power of 2
        let max_degree = (max_degree + 1).next_power_of_two() - 1;

        let setup_time = start_timer!(|| format!("Sampling {} generators", max_degree + 3));
        let mut generators = Self::sample_generators(max_degree + 3);
        end_timer!(setup_time);

        let h = generators.pop().unwrap();
        let s = generators.pop().unwrap();

        let pp = UniversalParams {
            comm_key: generators,
            h,
            s,
        };

        Ok(pp)
    }

    fn trim(
        pp: &Self::UniversalParams,
        supported_degree: usize,
        _enforced_degree_bounds: Option<&[usize]>,
    ) -> Result<(Self::CommitterKey, Self::VerifierKey), Self::Error> {
        // Ensure that supported_degree + 1 is a power of two
        let supported_degree = (supported_degree + 1).next_power_of_two() - 1;
        if supported_degree > pp.max_degree() {
            return Err(Error::TrimmingDegreeTooLarge);
        }

        let trim_time =
            start_timer!(|| format!("Trimming to supported degree of {}", supported_degree));

        let ck = CommitterKey {
            comm_key: pp.comm_key[0..(supported_degree + 1)].to_vec(),
            h: pp.h.clone(),
            s: pp.s.clone(),
            max_degree: pp.max_degree(),
        };

        let vk = VerifierKey {
            comm_key: pp.comm_key[0..(supported_degree + 1)].to_vec(),
            h: pp.h.clone(),
            s: pp.s.clone(),
            max_degree: pp.max_degree(),
        };

        end_timer!(trim_time);

        Ok((ck, vk))
    }

    /// Outputs a commitment to `polynomial`.
    fn commit<'a>(
        ck: &Self::CommitterKey,
        polynomials: impl IntoIterator<Item = &'a LabeledPolynomial<'a, G::ScalarField>>,
        rng: Option<&mut dyn RngCore>,
    ) -> Result<
        (
            Vec<LabeledCommitment<Self::Commitment>>,
            Vec<Self::Randomness>,
        ),
        Self::Error,
    > {
        let rng = &mut crate::optional_rng::OptionalRng(rng);
        let mut comms = Vec::new();
        let mut rands = Vec::new();

        let commit_time = start_timer!(|| "Committing to polynomials");
        for labeled_polynomial in polynomials {
            Self::check_degrees_and_bounds(ck.supported_degree(), labeled_polynomial)?;

            let polynomial = labeled_polynomial.polynomial();
            let label = labeled_polynomial.label();
            let hiding_bound = labeled_polynomial.hiding_bound();
            let degree_bound = labeled_polynomial.degree_bound();

            let commit_time = start_timer!(|| format!(
                "Polynomial {} of degree {}, degree bound {:?}, and hiding bound {:?}",
                label,
                polynomial.degree(),
                degree_bound,
                hiding_bound,
            ));

            let randomness = if let Some(h) = hiding_bound {
                Randomness::rand(h, degree_bound.is_some(), rng)
            } else {
                Randomness::empty()
            };

            let comm = Self::cm_commit(
                &ck.comm_key[..(polynomial.degree() + 1)],
                &polynomial.coeffs,
                Some(ck.s),
                Some(randomness.rand),
            )
            .into();

            let shifted_comm = degree_bound.map(|d| {
                Self::cm_commit(
                    &ck.comm_key[(ck.supported_degree() - d)..],
                    &polynomial.coeffs,
                    Some(ck.s),
                    randomness.shifted_rand,
                )
                .into()
            });

            let commitment = Commitment { comm, shifted_comm };
            let labeled_comm = LabeledCommitment::new(label.to_string(), commitment, degree_bound);

            comms.push(labeled_comm);
            rands.push(randomness);

            end_timer!(commit_time);
        }

        end_timer!(commit_time);
        Ok((comms, rands))
    }

    fn open<'a>(
        ck: &Self::CommitterKey,
        labeled_polynomials: impl IntoIterator<Item = &'a LabeledPolynomial<'a, G::ScalarField>>,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Self::Commitment>>,
        point: G::ScalarField,
        opening_challenge: G::ScalarField,
        rands: impl IntoIterator<Item = &'a Self::Randomness>,
        rng: Option<&mut dyn RngCore>,
    ) -> Result<Self::Proof, Self::Error>
    where
        Self::Commitment: 'a,
        Self::Randomness: 'a,
    {
        let mut combined_polynomial = Polynomial::zero();
        let mut combined_rand = G::ScalarField::zero();
        let mut combined_commitment_proj = G::Projective::zero();

        let mut has_hiding = false;

        let polys_iter = labeled_polynomials.into_iter();
        let rands_iter = rands.into_iter();
        let comms_iter = commitments.into_iter();

        let combine_time = start_timer!(|| "Combining polynomials, randomness, and commitments.");
        let mut cur_challenge = opening_challenge;
        for (labeled_polynomial, (labeled_commitment, randomness)) in
            polys_iter.zip(comms_iter.zip(rands_iter))
        {
            let label = labeled_polynomial.label();
            assert_eq!(labeled_polynomial.label(), labeled_commitment.label());
            Self::check_degrees_and_bounds(ck.supported_degree(), labeled_polynomial)?;

            let polynomial = labeled_polynomial.polynomial();
            let degree_bound = labeled_polynomial.degree_bound();
            let hiding_bound = labeled_polynomial.hiding_bound();
            let commitment = labeled_commitment.commitment();

            combined_polynomial += (cur_challenge, polynomial);
            combined_commitment_proj += &commitment.comm.mul(cur_challenge);

            if hiding_bound.is_some() {
                has_hiding = true;
                combined_rand += &(cur_challenge * &randomness.rand);
            }

            cur_challenge *= &opening_challenge;

            let has_degree_bound = degree_bound.is_some();

            assert_eq!(
                has_degree_bound,
                commitment.shifted_comm.is_some(),
                "shifted_comm mismatch for {}",
                label
            );

            assert_eq!(
                degree_bound,
                labeled_commitment.degree_bound(),
                "labeled_comm degree bound mismatch for {}",
                label
            );
            if let Some(degree_bound) = degree_bound {
                let shifted_polynomial = Self::shift_polynomial(ck, polynomial, degree_bound);
                combined_polynomial += (cur_challenge, &shifted_polynomial);
                combined_commitment_proj += &commitment.shifted_comm.unwrap().mul(cur_challenge);

                if hiding_bound.is_some() {
                    let shifted_rand = randomness.shifted_rand;
                    assert!(
                        shifted_rand.is_some(),
                        "shifted_rand.is_none() for {}",
                        label
                    );
                    combined_rand += &(cur_challenge * &shifted_rand.unwrap());
                }
            }

            cur_challenge *= &opening_challenge;
        }

        end_timer!(combine_time);

        let combined_v = combined_polynomial.evaluate(point);

        // Pad the coefficients to the appropriate vector size
        let d = ck.supported_degree();

        // `log_d` is ceil(log2 (d + 1)), which is the number of steps to compute all of the challenges
        let log_d = algebra_core::log2(d + 1) as usize;

        let mut combined_commitment;
        let mut hiding_commitment = None;

        if has_hiding {
            let mut rng = rng.expect("hiding commitments require randomness");
            let hiding_time = start_timer!(|| "Applying hiding.");
            let mut hiding_polynomial = Polynomial::rand(d, &mut rng);
            hiding_polynomial -=
                &Polynomial::from_coefficients_slice(&[hiding_polynomial.evaluate(point)]);

            let hiding_rand = G::ScalarField::rand(rng);
            let hiding_commitment_proj = Self::cm_commit(
                ck.comm_key.as_slice(),
                hiding_polynomial.coeffs.as_slice(),
                Some(ck.s),
                Some(hiding_rand),
            );

            let mut batch = G::Projective::batch_normalization_into_affine(&[
                combined_commitment_proj,
                hiding_commitment_proj,
            ]);
            hiding_commitment = Some(batch.pop().unwrap());
            combined_commitment = batch.pop().unwrap();

            let hiding_challenge = Self::compute_random_oracle_challenge(
                &algebra_core::to_bytes![
                    combined_commitment,
                    point,
                    combined_v,
                    hiding_commitment.unwrap()
                ]
                .unwrap(),
            );
            combined_polynomial += (hiding_challenge, &hiding_polynomial);
            combined_rand += &(hiding_challenge * &hiding_rand);
            combined_commitment_proj +=
                &(hiding_commitment_proj.mul(hiding_challenge) - &ck.s.mul(combined_rand));

            end_timer!(hiding_time);
        }

        let combined_rand = if has_hiding {
            Some(combined_rand)
        } else {
            None
        };

        let proof_time =
            start_timer!(|| format!("Generating proof for degree {} combined polynomial", d + 1));

        combined_commitment = combined_commitment_proj.into_affine();

        // ith challenge
        let mut round_challenge = Self::compute_random_oracle_challenge(
            &algebra_core::to_bytes![combined_commitment, point, combined_v].unwrap(),
        );

        let h_prime = ck.h.mul(round_challenge).into_affine();

        // Pads the coefficients with zeroes to get the number of coeff to be d+1
        let mut coeffs = combined_polynomial.coeffs;
        if coeffs.len() < d + 1 {
            for _ in coeffs.len()..(d + 1) {
                coeffs.push(G::ScalarField::zero());
            }
        }
        let mut coeffs = coeffs.as_mut_slice();

        // Powers of z
        let mut z: Vec<G::ScalarField> = Vec::with_capacity(d + 1);
        let mut cur_z: G::ScalarField = G::ScalarField::one();
        for _ in 0..(d + 1) {
            z.push(cur_z);
            cur_z *= &point;
        }
        let mut z = z.as_mut_slice();

        // This will be used for transforming the key in each step
        let mut key_proj: Vec<G::Projective> = ck.comm_key.iter().map(|x| (*x).into()).collect();
        let mut key_proj = key_proj.as_mut_slice();

        let mut temp;

        // Key for MSM
        // We initialize this to capacity 0 initially because we want to use the key slice first
        let mut comm_key = &ck.comm_key;

        let mut l_vec = Vec::with_capacity(log_d);
        let mut r_vec = Vec::with_capacity(log_d);

        let mut n = d + 1;
        while n > 1 {
            let (coeffs_l, coeffs_r) = coeffs.split_at_mut(n / 2);
            let (z_l, z_r) = z.split_at_mut(n / 2);
            let (key_l, key_r) = comm_key.split_at(n / 2);
            let (key_proj_l, _) = key_proj.split_at_mut(n / 2);

            let l = Self::cm_commit(key_l, coeffs_r, None, None)
                + &h_prime.mul(Self::inner_product(coeffs_r, z_l));

            let r = Self::cm_commit(key_r, coeffs_l, None, None)
                + &h_prime.mul(Self::inner_product(coeffs_l, z_r));

            let lr = G::Projective::batch_normalization_into_affine(&[l, r]);
            l_vec.push(lr[0]);
            r_vec.push(lr[1]);

            round_challenge = Self::compute_random_oracle_challenge(
                &algebra_core::to_bytes![round_challenge, lr[0], lr[1]].unwrap(),
            );
            let round_challenge_inv = round_challenge.inverse().unwrap();

            ff_fft::cfg_iter_mut!(coeffs_l)
                .zip(coeffs_r)
                .for_each(|(c_l, c_r)| *c_l += &(round_challenge_inv * &c_r));

            ff_fft::cfg_iter_mut!(z_l)
                .zip(z_r)
                .for_each(|(z_l, z_r)| *z_l += &(round_challenge * &z_r));

            ff_fft::cfg_iter_mut!(key_proj_l)
                .zip(key_r)
                .for_each(|(k_l, k_r)| *k_l += &(k_r.mul(round_challenge)));

            coeffs = coeffs_l;
            z = z_l;

            key_proj = key_proj_l;
            temp = G::Projective::batch_normalization_into_affine(key_proj);
            comm_key = &temp;

            n /= 2;
        }

        end_timer!(proof_time);

        Ok(Proof {
            l_vec,
            r_vec,
            final_comm_key: comm_key[0],
            c: coeffs[0],
            hiding_comm: hiding_commitment,
            rand: combined_rand,
        })
    }

    fn check<'a, R: RngCore>(
        vk: &Self::VerifierKey,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Self::Commitment>>,
        point: G::ScalarField,
        values: impl IntoIterator<Item = G::ScalarField>,
        proof: &Self::Proof,
        opening_challenge: G::ScalarField,
        _rng: &mut R,
    ) -> Result<bool, Self::Error>
    where
        Self::Commitment: 'a,
    {
        let check_time = start_timer!(|| "Checking evaluations");
        let d = vk.supported_degree();

        // `log_d` is ceil(log2 (d + 1)), which is the number of steps to compute all of the challenges
        let log_d = algebra_core::log2(d + 1) as usize;

        if proof.l_vec.len() != proof.r_vec.len() || proof.l_vec.len() != log_d {
            return Err(Error::IncorrectInputLength(
                format!(
                    "Expected proof vectors to be {:}. Instead, l_vec size is {:} and r_vec size is {:}",
                    log_d,
                    proof.l_vec.len(),
                    proof.r_vec.len()
                )
            ));
        }

        let check_poly =
            Self::succinct_check(vk, commitments, point, values, proof, opening_challenge);

        if check_poly.is_none() {
            return Ok(false);
        }

        let check_poly_coeffs = check_poly.unwrap().compute_coeffs();
        let final_key = Self::cm_commit(
            vk.comm_key.as_slice(),
            check_poly_coeffs.as_slice(),
            None,
            None,
        );
        if !(final_key - &proof.final_comm_key.into()).is_zero() {
            return Ok(false);
        }

        end_timer!(check_time);
        Ok(true)
    }

    fn batch_check<'a, R: RngCore>(
        vk: &Self::VerifierKey,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Self::Commitment>>,
        query_set: &QuerySet<G::ScalarField>,
        values: &Evaluations<G::ScalarField>,
        proof: &Self::BatchProof,
        opening_challenge: G::ScalarField,
        rng: &mut R,
    ) -> Result<bool, Self::Error>
    where
        Self::Commitment: 'a,
    {
        let commitments: BTreeMap<_, _> = commitments.into_iter().map(|c| (c.label(), c)).collect();
        let mut query_to_labels_map = BTreeMap::new();

        for (label, point) in query_set.iter() {
            let labels = query_to_labels_map.entry(point).or_insert(BTreeSet::new());
            labels.insert(label);
        }

        assert_eq!(proof.len(), query_to_labels_map.len());

        let mut randomizer = G::ScalarField::one();

        let mut combined_check_poly = Polynomial::zero();
        let mut combined_final_key = G::Projective::zero();

        for ((query, labels), p) in query_to_labels_map.into_iter().zip(proof) {
            let lc_time =
                start_timer!(|| format!("Randomly combining {} commitments", labels.len()));
            let mut comms: Vec<&'_ LabeledCommitment<_>> = Vec::new();
            let mut vals = Vec::new();
            for label in labels.into_iter() {
                let commitment = commitments.get(label).ok_or(Error::MissingPolynomial {
                    label: label.to_string(),
                })?;

                let v_i = values
                    .get(&(label.clone(), *query))
                    .ok_or(Error::MissingEvaluation {
                        label: label.to_string(),
                    })?;

                comms.push(commitment);
                vals.push(*v_i);
            }

            let check_poly = Self::succinct_check(
                vk,
                comms.into_iter(),
                *query,
                vals.into_iter(),
                p,
                opening_challenge,
            );

            if check_poly.is_none() {
                return Ok(false);
            }

            let check_poly =
                Polynomial::from_coefficients_vec(check_poly.unwrap().compute_coeffs());
            combined_check_poly += (randomizer, &check_poly);
            combined_final_key += &p.final_comm_key.into_projective().mul(randomizer);

            randomizer = u128::rand(rng).into();
            end_timer!(lc_time);
        }

        let proof_time = start_timer!(|| "Checking batched proof");
        let final_key = Self::cm_commit(
            vk.comm_key.as_slice(),
            combined_check_poly.coeffs.as_slice(),
            None,
            None,
        );
        if !(final_key - &combined_final_key).is_zero() {
            return Ok(false);
        }

        end_timer!(proof_time);

        Ok(true)
    }

    fn open_combinations<'a>(
        ck: &Self::CommitterKey,
        lc_s: impl IntoIterator<Item = &'a LinearCombination<G::ScalarField>>,
        polynomials: impl IntoIterator<Item = &'a LabeledPolynomial<'a, G::ScalarField>>,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Self::Commitment>>,
        query_set: &QuerySet<G::ScalarField>,
        opening_challenge: G::ScalarField,
        rands: impl IntoIterator<Item = &'a Self::Randomness>,
        rng: Option<&mut dyn RngCore>,
    ) -> Result<BatchLCProof<G::ScalarField, Self>, Self::Error>
    where
        Self::Randomness: 'a,
        Self::Commitment: 'a,
    {
        let label_poly_map = polynomials
            .into_iter()
            .zip(rands)
            .zip(commitments)
            .map(|((p, r), c)| (p.label(), (p, r, c)))
            .collect::<BTreeMap<_, _>>();

        let mut lc_polynomials = Vec::new();
        let mut lc_randomness = Vec::new();
        let mut lc_commitments = Vec::new();
        let mut lc_info = Vec::new();

        for lc in lc_s {
            let lc_label = lc.label().clone();
            let mut poly = Polynomial::zero();
            let mut degree_bound = None;
            let mut hiding_bound = None;

            let mut combined_comm = G::Projective::zero();
            let mut combined_shifted_comm: Option<G::Projective> = None;

            let mut combined_rand = G::ScalarField::zero();
            let mut combined_shifted_rand: Option<G::ScalarField> = None;

            let num_polys = lc.len();
            for (coeff, label) in lc.iter().filter(|(_, l)| !l.is_one()) {
                let label: &String = label.try_into().expect("cannot be one!");
                let &(cur_poly, cur_rand, cur_comm) =
                    label_poly_map.get(label).ok_or(Error::MissingPolynomial {
                        label: label.to_string(),
                    })?;

                if num_polys == 1 && cur_poly.degree_bound().is_some() {
                    assert!(
                        coeff.is_one(),
                        "Coefficient must be one for degree-bounded equations"
                    );
                    degree_bound = cur_poly.degree_bound();
                } else if cur_poly.degree_bound().is_some() {
                    eprintln!("Degree bound when number of equations is non-zero");
                    return Err(Self::Error::EquationHasDegreeBounds(lc_label));
                }

                // Some(_) > None, always.
                hiding_bound = core::cmp::max(hiding_bound, cur_poly.hiding_bound());
                poly += (*coeff, cur_poly.polynomial());

                combined_rand += &(cur_rand.rand * coeff);
                combined_shifted_rand = Self::combine_shifted_rand(
                    combined_shifted_rand,
                    cur_rand.shifted_rand,
                    *coeff,
                );

                let commitment = cur_comm.commitment();
                combined_comm += &commitment.comm.mul(*coeff);
                combined_shifted_comm = Self::combine_shifted_comm(
                    combined_shifted_comm,
                    commitment.shifted_comm,
                    *coeff,
                );
            }

            let lc_poly =
                LabeledPolynomial::new_owned(lc_label.clone(), poly, degree_bound, hiding_bound);
            lc_polynomials.push(lc_poly);
            lc_randomness.push(Randomness {
                rand: combined_rand,
                shifted_rand: combined_shifted_rand,
            });

            lc_commitments.push(combined_comm);
            if let Some(combined_shifted_comm) = combined_shifted_comm {
                lc_commitments.push(combined_shifted_comm);
            }

            lc_info.push((lc_label, degree_bound));
        }

        let lc_commitments = Self::construct_labeled_commitments(&lc_info, &lc_commitments);

        let proof = Self::batch_open(
            ck,
            lc_polynomials.iter(),
            lc_commitments.iter(),
            &query_set,
            opening_challenge,
            lc_randomness.iter(),
            rng,
        )?;
        Ok(BatchLCProof { proof, evals: None })
    }

    /// Checks that `values` are the true evaluations at `query_set` of the polynomials
    /// committed in `labeled_commitments`.
    fn check_combinations<'a, R: RngCore>(
        vk: &Self::VerifierKey,
        lc_s: impl IntoIterator<Item = &'a LinearCombination<G::ScalarField>>,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Self::Commitment>>,
        query_set: &QuerySet<G::ScalarField>,
        evaluations: &Evaluations<G::ScalarField>,
        proof: &BatchLCProof<G::ScalarField, Self>,
        opening_challenge: G::ScalarField,
        rng: &mut R,
    ) -> Result<bool, Self::Error>
    where
        Self::Commitment: 'a,
    {
        let BatchLCProof { proof, .. } = proof;
        let label_comm_map = commitments
            .into_iter()
            .map(|c| (c.label(), c))
            .collect::<BTreeMap<_, _>>();

        let mut lc_commitments = Vec::new();
        let mut lc_info = Vec::new();
        let mut evaluations = evaluations.clone();
        for lc in lc_s {
            let lc_label = lc.label().clone();
            let num_polys = lc.len();

            let mut degree_bound = None;
            let mut combined_comm = G::Projective::zero();
            let mut combined_shifted_comm: Option<G::Projective> = None;

            for (coeff, label) in lc.iter() {
                if label.is_one() {
                    for (&(ref label, _), ref mut eval) in evaluations.iter_mut() {
                        if label == &lc_label {
                            **eval -= coeff;
                        }
                    }
                } else {
                    let label: &String = label.try_into().unwrap();
                    let &cur_comm = label_comm_map.get(label).ok_or(Error::MissingPolynomial {
                        label: label.to_string(),
                    })?;

                    if num_polys == 1 && cur_comm.degree_bound().is_some() {
                        assert!(
                            coeff.is_one(),
                            "Coefficient must be one for degree-bounded equations"
                        );
                        degree_bound = cur_comm.degree_bound();
                    } else if cur_comm.degree_bound().is_some() {
                        return Err(Self::Error::EquationHasDegreeBounds(lc_label));
                    }

                    let commitment = cur_comm.commitment();
                    combined_comm += &commitment.comm.mul(*coeff);
                    combined_shifted_comm = Self::combine_shifted_comm(
                        combined_shifted_comm,
                        commitment.shifted_comm,
                        *coeff,
                    );
                }
            }

            lc_commitments.push(combined_comm);

            if let Some(combined_shifted_comm) = combined_shifted_comm {
                lc_commitments.push(combined_shifted_comm);
            }

            lc_info.push((lc_label, degree_bound));
        }

        let lc_commitments = Self::construct_labeled_commitments(&lc_info, &lc_commitments);

        Self::batch_check(
            vk,
            &lc_commitments,
            &query_set,
            &evaluations,
            proof,
            opening_challenge,
            rng,
        )
    }
}

#[cfg(test)]
mod tests {
    #![allow(non_camel_case_types)]

    use super::InnerProductArgPC;

    use algebra::ed_on_bls12_381::EdwardsAffine;
    use blake2::Blake2s;

    type PC<E, D> = InnerProductArgPC<E, D>;
    type PC_JJB2S = PC<EdwardsAffine, Blake2s>;

    #[test]
    fn single_poly_test() {
        use crate::tests::*;
        single_poly_test::<_, PC_JJB2S>().expect("test failed for ed_on_bls12_381-blake2s");
    }

    #[test]
    fn quadratic_poly_degree_bound_multiple_queries_test() {
        use crate::tests::*;
        quadratic_poly_degree_bound_multiple_queries_test::<_, PC_JJB2S>()
            .expect("test failed for ed_on_bls12_381-blake2s");
    }

    #[test]
    fn linear_poly_degree_bound_test() {
        use crate::tests::*;
        linear_poly_degree_bound_test::<_, PC_JJB2S>()
            .expect("test failed for ed_on_bls12_381-blake2s");
    }

    #[test]
    fn single_poly_degree_bound_test() {
        use crate::tests::*;
        single_poly_degree_bound_test::<_, PC_JJB2S>()
            .expect("test failed for ed_on_bls12_381-blake2s");
    }

    #[test]
    fn single_poly_degree_bound_multiple_queries_test() {
        use crate::tests::*;
        single_poly_degree_bound_multiple_queries_test::<_, PC_JJB2S>()
            .expect("test failed for ed_on_bls12_381-blake2s");
    }

    #[test]
    fn two_polys_degree_bound_single_query_test() {
        use crate::tests::*;
        two_polys_degree_bound_single_query_test::<_, PC_JJB2S>()
            .expect("test failed for ed_on_bls12_381-blake2s");
    }

    #[test]
    fn full_end_to_end_test() {
        use crate::tests::*;
        full_end_to_end_test::<_, PC_JJB2S>().expect("test failed for ed_on_bls12_381-blake2s");
        println!("Finished ed_on_bls12_381-blake2s");
    }

    #[test]
    fn single_equation_test() {
        use crate::tests::*;
        single_equation_test::<_, PC_JJB2S>().expect("test failed for ed_on_bls12_381-blake2s");
        println!("Finished ed_on_bls12_381-blake2s");
    }

    #[test]
    fn two_equation_test() {
        use crate::tests::*;
        two_equation_test::<_, PC_JJB2S>().expect("test failed for ed_on_bls12_381-blake2s");
        println!("Finished ed_on_bls12_381-blake2s");
    }

    #[test]
    fn two_equation_degree_bound_test() {
        use crate::tests::*;
        two_equation_degree_bound_test::<_, PC_JJB2S>()
            .expect("test failed for ed_on_bls12_381-blake2s");
        println!("Finished ed_on_bls12_381-blake2s");
    }

    #[test]
    fn full_end_to_end_equation_test() {
        use crate::tests::*;
        full_end_to_end_equation_test::<_, PC_JJB2S>()
            .expect("test failed for ed_on_bls12_381-blake2s");
        println!("Finished ed_on_bls12_381-blake2s");
    }

    #[test]
    #[should_panic]
    fn bad_degree_bound_test() {
        use crate::tests::*;
        bad_degree_bound_test::<_, PC_JJB2S>().expect("test failed for ed_on_bls12_381-blake2s");
        println!("Finished ed_on_bls12_381-blake2s");
    }
}
