//! Here we constuct a polynomial commitment that enables users to commit to a
//! single polynomial `p`, and then later provide an evaluation proof that
//! convinces verifiers that a claimed value `v` is the true evaluation of `p`
//! at a chosen point `x`. Our construction follows the template of the construction
//! proposed by Kate, Zaverucha, and Goldberg ([KZG10](http://cacr.uwaterloo.ca/techreports/2010/cacr2010-10.pdf)).
//! This construction achieves extractability in the algebraic group model (AGM).

use crate::{
    PCUniversalParams, PCCommitment, PCCommitterKey, PCRandomness, PCVerifierKey, 
    Polynomial, PolynomialCommitment,
};
use crate::{LabeledPolynomial, LabeledCommitment};
use algebra::bytes::*;
use algebra::{
    AffineCurve, Field, Group, PairingEngine, ProjectiveCurve, UniformRand};
use rand_core::RngCore;
use rayon::prelude::*;
use std::marker::PhantomData;

use crate::kzg10;

/// `KZG10` is an implementation of the polynomial commitment scheme of
/// [Kate, Zaverucha and Goldbgerg][kzg10]
///
/// [kzg10]: http://cacr.uwaterloo.ca/techreports/2010/cacr2010-10.pdf
pub struct MarlinKZG10<E: PairingEngine> {
    _engine: PhantomData<E>,
}

/// `UniversalParams` are the universal parameters for the KZG10 scheme.
pub type UniversalParams<E> = kzg10::UniversalParams<E>;

/// `ComitterKey` is used to commit to and create evaluation proofs for a given
/// polynomial.
#[derive(Derivative)]
#[derivative(
    Default(bound = ""),
    Hash(bound = ""),
    Clone(bound = ""),
    Debug(bound = "")
)]
pub struct CommitterKey<E: PairingEngine> {
    /// The key used to commit to polynomials.
    pub ck: kzg10::CommitterKey<E>,
    /// The key used to commit to shifted polynomials.
    /// This is `None` if `self` does not support enforcing any degree bounds.
    pub shifted_ck: Option<kzg10::CommitterKey<E>>,
    /// The degree bounds that are supported by `self`.
    /// In ascending order from smallest to largest.
    /// This is `None` if `self` does not support enforcing any degree bounds.
    pub supported_degree_bounds: Option<Vec<usize>>,
}

impl<E: PairingEngine> PCCommitterKey for CommitterKey<E> {
    fn max_degree(&self) -> usize {
        self.ck.max_degree
    }
}

/// `VerifierKey` is used to check evaluation proofs for a given commitment.
#[derive(Derivative)]
#[derivative(Default(bound = ""), Clone(bound = ""), Debug(bound = ""))]
pub struct VerifierKey<E: PairingEngine> {
    /// The verification key for the underlying KZG10 scheme.
    pub vk: kzg10::VerifierKey<E>,
    /// Information required to enforce degree bounds. Each pair
    /// is of the form `(degree_bound, shifting_advice)`.
    /// The vector is sorted in ascending order of `degree_bound`. 
    /// This is `None` if `self` does not support enforcing any degree bounds.
    pub degree_bound_shifts: Option<Vec<(usize, E::G1Affine)>>,
}

impl<E: PairingEngine> PCVerifierKey for VerifierKey<E> {
    fn max_degree(&self) -> usize {
        self.vk.max_degree
    }
}

/// Commitment to a polynomial that optionally enforces a degree bound.
#[derive(Derivative)]
#[derivative(
    Default(bound = ""),
    Hash(bound = ""),
    Clone(bound = ""),
    Copy(bound = ""),
    Debug(bound = ""),
    PartialEq(bound = ""),
    Eq(bound = "")
)]
pub struct Commitment<E: PairingEngine> {
    comm: kzg10::Commitment<E>,
    shifted_comm: Option<kzg10::Commitment<E>>,
}

impl<E: PairingEngine> ToBytes for Commitment<E> {
    #[inline]
    fn write<W: std::io::Write>(&self, mut writer: W) -> std::io::Result<()> {
        self.comm.write(&mut writer)?;
        let shifted_exists = self.shifted_comm.is_some();
        shifted_exists.write(&mut writer)?;
        self.shifted_comm
            .as_ref()
            .unwrap_or(&kzg10::Commitment::empty())
            .write(&mut writer)
    }
}

impl<E: PairingEngine> PCCommitment for Commitment<E> {
    #[inline]
    fn empty() -> Self {
        Self {
            comm: kzg10::Commitment::empty(),
            shifted_comm: Some(kzg10::Commitment::empty()),
        }
    }

    fn has_degree_bound(&self) -> bool {
        self.shifted_comm.is_some()
    }

    fn size_in_bytes(&self) -> usize {
        self.comm.size_in_bytes() + self.shifted_comm.as_ref().map_or(0, |c| c.size_in_bytes())
    }
}

/// `Randomness` hides the polynomial inside a commitment. It is output by `KZG10::commit`.
#[derive(Derivative)]
#[derivative(
    Default(bound = ""),
    Hash(bound = ""),
    Clone(bound = ""),
    Debug(bound = ""),
    PartialEq(bound = ""),
    Eq(bound = "")
)]
pub struct Randomness<E: PairingEngine> {
    rand: kzg10::Randomness<E>,
    shifted_rand: Option<kzg10::Randomness<E>>,
}

impl<E: PairingEngine> PCRandomness for Randomness<E> {
    fn empty() -> Self {
        Self {
            rand: kzg10::Randomness::empty(),
            shifted_rand: Some(kzg10::Randomness::empty()),
        }
    }

    fn rand<R: RngCore>(hiding_bound: usize, rng: &mut R) -> Self {
        Self {
            rand: kzg10::Randomness::rand(hiding_bound, rng),
            shifted_rand: Some(kzg10::Randomness::rand(hiding_bound, rng)),
        }
    }
}

/// Evaluation proof output by `MultiPCFromSinglePC::batch_open`.
#[derive(Derivative)]
#[derivative(
    Default(bound = ""),
    Hash(bound = ""),
    Clone(bound = ""),
    Debug(bound = ""),
    PartialEq(bound = ""),
    Eq(bound = "")
)]
pub struct BatchProof<E: PairingEngine>(Vec<kzg10::Proof<E>>);

pub(crate) fn shift_polynomial<F: Field>(
    p: &Polynomial<F>,
    degree_bound: usize,
    max_degree: usize,
) -> Polynomial<F> {
    if p.is_zero() {
        Polynomial::zero()
    } else {
        let mut shifted_polynomial_coeffs = vec![F::zero(); max_degree - degree_bound];
        shifted_polynomial_coeffs.extend_from_slice(&p.coeffs);
        Polynomial::from_coefficients_vec(shifted_polynomial_coeffs)
    }
}

impl<E: PairingEngine> MarlinKZG10<E> {
    /// Accumulated `commitments` and `values` according to `opening_challenge`.
    fn accumulate_commitments_and_values<'a, R: RngCore>(
        vk: &Self::VerifierKey,
        commitments: impl IntoIterator<Item = &'a Commitment<E>>,
        point: E::Fr,
        values: impl IntoIterator<Item = E::Fr>,
        opening_challenge: E::Fr,
    ) -> Result<(E::G1Projective, E::Fr), kzg10::Error> {
        let acc_time = start_timer!(|| "Accumulating commitments and values");
        let mut combined_comm = E::G1Projective::zero();
        let mut combined_value = E::Fr::zero();
        let mut challenge_i = E::Fr::one();
        for (labeled_commitment, value) in commitments.into_iter().zip(values) {
            let degree_bound = labeled_commitment.degree_bound();
            let commitment = labeled_commitment.commitment();
            assert_eq!(degree_bound.is_some(), commitment.shifted_comm.is_some());

            combined_comm += &commitment.comm.0.mul(challenge_i);
            combined_value += &(value * challenge_i);

            if let Some(degree_bound) = degree_bound {
                let challenge_i_1 = challenge_i * &opening_challenge;
                let shifts = vk.degree_bound_shifts.as_ref().unwrap();
                let shift_index = shifts.binary_search(degree_bound);
                // Is this v_i or point^(D - d_i) * v_i
                let shift = shifts[shift_index].mul(value);
                combined_comm += &(commitment.shifted_comm.as_ref().unwrap().clone() - &shift);
            }
            challenge_i *= &opening_challenge.square();
        }

        end_timer!(acc_time);
        Ok((combined_comm, combined_value))
    }

}

impl<E: PairingEngine> PolynomialCommitment<E::Fr> for MarlinKZG10<E> {
    type UniversalParams = UniversalParams<E>;
    type CommitterKey = CommitterKey<E>;
    type VerifierKey = VerifierKey<E>;
    type Commitment = Commitment<E>;
    type Randomness = Randomness<E>;
    type Proof = kzg10::Proof<E>;
    type BatchProof = Vec<Self::Proof>;
    type Error = kzg10::Error;

    /// Constructs public parameters when given as input the maximum degree `degree`
    /// for the polynomial commitment scheme.
    fn setup<R: RngCore>(
        max_degree: usize,
        rng: &mut R,
    ) -> Result<Self::UniversalParams, Self::Error> {
        kzg10::KZG10::setup(max_degree, rng)
    }

    fn trim(
        pp: &Self::UniversalParams,
        max_degree: usize,
        degree_bounds_to_support: Option<&[usize]>,
    ) -> Result<(Self::CommitterKey, Self::VerifierKey), Self::Error> {
        let pp_max_degree = pp.max_degree();
        // TODO: Make error.
        assert!(max_degree <= pp_max_degree, "Trimming cannot create powers");

        // Construct the KZG10 committer key for committing to unshifted polynomials.
        let ck_time = start_timer!(|| format!("Constructing ck of size {} for unshifted polys", max_degree));
        let mut powers_of_g = Vec::new();
        let mut powers_of_gamma_g = Vec::new();

        let powers_of_g = pp.powers_of_g[..=max_degree].to_vec();
        let powers_of_gamma_g = pp.powers_of_gamma_g[..=max_degree].to_vec();

        let ck = CommitterKey {
            powers_of_g,
            powers_of_gamma_g,
            max_degree,
        };
        end_timer!(ck_time);

        // Construct the core KZG10 verifier key.
        let vk = kzg10::VerifierKey {
            g: pp.powers_of_g[0],
            gamma_g: pp.powers_of_gamma_g[0],
            h: pp.h,
            beta_h: pp.beta_h,
            prepared_h: pp.prepared_h,
            prepared_beta_h: pp.prepared_beta_h,
            max_degree,
        }:

        // Check whether we have some degree bounds to enforce
        let (shifted_ck, degree_bound_shifts) = if let Some(degree_bounds_to_support) = degree_bounds_to_support {
            let mut degree_bounds_to_support = degree_bounds_to_support.to_vec();
            degree_bounds_to_support.sort();
            // TODO: Make error instead of `expect`.
            let lowest_shifted_power = pp_max_degree - degree_bounds_to_support.last().expect("degree_bounds_to_support should not be empty");

            let shifted_ck_time = start_timer!(|| format!(
                    "Constructing shifted_ck of size {} for shifted polys",
                    pp_max_degree - lowest_shifted_power
            ));

            let mut powers_of_g = Vec::new();
            let mut powers_of_gamma_g = Vec::new();

            let powers_of_g = pp.powers_of_g[lowest_shifted_power..].to_vec();
            let powers_of_gamma_g = pp.powers_of_gamma_g[..(pp_max_degree - lowest_shifted_power)].to_vec();

            let shifted_ck = CommitterKey {
                powers_of_g,
                powers_of_gamma_g,
                max_degree: pp_max_degree - lowest_shifted_power,
            };
            
            end_timer!(shifted_ck_time);

            let degree_bound_shifts = degree_bounds_to_support.into_iter().map(|d| (d, pp.powers_of_g[pp_max_degree - d])).collect();
            (Some(shifted_ck), Some(degree_bound_shifts))
        } else {
            (None, None)
        }:

        let ck = {
            ck,
            shifted_ck,
            supported_degree_bounds: degree_bounds_to_support.cloned(),
        };

        let vk = VerifierKey {
            vk,
            degree_bound_shifts,
        }
        Ok(ck, vk)
    }

    /// Outputs a commitment to `polynomial`.
    fn commit<'a>(
        ck: &Self::CommitterKey,
        polynomials: impl IntoIterator<Item = &'a LabeledPolynomial<'a, E::Fr>>,
        rng: Option<&mut dyn RngCore>,
    ) -> Result<
        (
            Vec<LabeledCommitment<Self::Commitment>>,
            Vec<Self::Randomness>,
        ),
        Self::Error,
    > {
        let commit_time = start_timer!(|| "Committing to polynomials");

        let mut commitments = Vec::new();
        let mut randomness = Vec::new();
        let max_degree = ck.max_degree();

        let rng = &mut optional_rng::OptionalRng(rng);
        for polynomial in polynomials {
            let label = polynomial.label();
            let degree_bound = polynomial.degree_bound();
            let hiding_bound = polynomial.hiding_bound();
            let polynomial = polynomial.polynomial();

            Self::Error::check_degrees(
                polynomial.degree(),
                degree_bound,
                max_degree,
                label.to_string(),
            )?;

            let commit_time = start_timer!(|| format!(
                "Polynomial {} of degree {}, degree bound {:?}, and hiding bound {:?}",
                label,
                polynomial.degree(),
                degree_bound,
                hiding_bound,
            ));
            let (comm, rand) = kzg10::KZG10::commit(&ck.ck, polynomial, hiding_bound, Some(rng))?;
            let (shifted_comm, shifted_rand) = if let Some(degree_bound) = degree_bound {
                // TODO: Convert to errors
                let shifted_ck = ck.shifted_ck.as_ref().expect("Polynomial requires degree bounds, but `ck` does not support any");
                let supported_degree_bounds = ck.supported_degree_bounds.as_ref().expect("Polynomial requires degree bounds, but `ck` does not support any");

                let max_degree_bound = supported_degree_bounds.last().unwrap();

                // TODO: Convert to errors
                let _ = supported_degree_bounds.binary_search(&degree_bound).expect("Unsupported degree bound");
                let mut shifted_polynomial_coeffs = vec![E::Fr::zero(); max_degree_bound - degree_bound];
                shifted_polynomial_coeffs.extend_from_slice(&polynomial.coeffs);
                let s_polynomial = Polynomial::from_coefficients_vec(shifted_polynomial_coeffs);
                assert!(
                    s_polynomial.degree() <= shifted_ck.max_degree(), 
                    "polynomial.degree(): {}; s_polynomial.degree(): {}; max_degree: {}.",
                    polynomial.degree(),
                    s_polynomial.degree(),
                    shifted_ck.max_degree()
                );
                let (shifted_comm, shifted_rand) =
                    kzg10::KZG10::commit(&ck.shifted_ck, &s_polynomial, hiding_bound, Some(rng))?;
                (Some(shifted_comm), Some(shifted_rand))
            } else {
                (None, None)
            };

            let comm = Commitment { comm, shifted_comm };
            let rand = Randomness { rand, shifted_rand };
            commitments.push(LabeledCommitment::new(label.to_string(), comm, degree_bound));
            randomness.push(rand);
            end_timer!(commit_time);
        }
        end_timer!(commit_time);
        Ok((commitments, randomness))
    }

    /// On input a polynomial `p` and a point `point`, outputs a proof for the same.
    fn open<'a>(
        ck: &Self::CommitterKey,
        labeled_polynomials: impl IntoIterator<Item = &'a LabeledPolynomial<'a, E::Fr>>,
        point: E::Fr,
        opening_challenge: E::Fr,
        rands: impl IntoIterator<Item = &'a Self::Randomness<E::Fr>>,
    ) -> Result<Self::Proof, Self::Error> {
        let max_degree = ck.max_degree();
        let mut p = Polynomial::zero();
        let mut r = kzg10::Randomness::empty();
        let mut shifted_w = Polynomial::zero();
        let mut shifted_r = kzg10::Randomness::empty();
        let mut shifted_r_witness = kzg10::Randomness::empty();
        for ((polynomial, rand), j) in labeled_polynomials.into_iter().zip(rands).enumerate() {
            let degree_bound = polynomial.degree_bound();
            let label = polynomial.label();

            Self::Error::check_degrees(
                polynomial.degree(),
                degree_bound,
                max_degree,
                label.to_string(),
            )?;

            // compute challenge^j and challenge^{j+1}.
            let challenge_j = opening_challenge.pow([2 * j as u64]);

            assert_eq!(degree_bound.is_some(), rand.shifted_rand.is_some());

            p += (challenge_j, polynomial.polynomial());
            r += (challenge_j, &rand.rand);

            if let Some(degree_bound) = degree_bound {
                let shifted_rand = rand.shifted_rand.as_ref().unwrap();
                let (witness, shifted_rand_witness) = kzg10::KZG10::compute_witness_polynomial(
                    polynomial.polynomial(),
                    point,
                    &rand
                );
                let challenge_j_1 = challenge_j * &opening_challenge;

                let shifted_witness = shift_polynomial(&witness, degree_bound, max_degree);

                shifted_w += (challenge_j_1, &shifted_witness);
                shifted_r += (challenge_j_1, &shifted_rand);
                shifted_r_witness += (challenge_j_1, &shifted_rand_witness);
            }
        }
        let proof_time = start_timer!(|| "Creating proof for unshifted polynomials");
        let proof = kzg10::KZG10::open(ck, &p, point, &r)?;
        end_timer!(proof_time);

        let proof_time = start_timer!(|| "Creating proof for shifted polynomials");
        let shifted_proof = kzg10::KZG10::open_with_witness_polynomial(ck, point, &r, &shifted_w, &shifted_r_witness)?;
        end_timer!(proof_time);

        Ok(kzg10::Proof {
            w: proof.w.into_projective() + &shifted_proof.w.into_projective(),
            random_v: proof.random_v + &shifted_proof.random_v,
        })
    }


    /// Verifies that `value` is the evaluation at `x` of the polynomial
    /// committed inside `comm`.
    fn check<'a, R: RngCore>(
        vk: &Self::VerifierKey,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Self::Commitment>>,
        point: E::Fr,
        values: impl IntoIterator<Item = E::Fr>,
        proof: &Self::Proof,
        opening_challenge: E::Fr,
        rng: &mut R,
    ) -> Result<bool, Self::Error> {
        let check_time = start_timer!(|| "Checking evaluations");
        let (combined_comm, combined_value) = 
            Self::accumulate_commitments_and_values(vk, commitments, point, values, proof, opening_challenge, rng)?;
        let result = kzg10::KZG10::check(&vk.vk, combined_comm.comm, point, combined_value, proof);
        end_timer!(check_time);
        result
    }

    fn batch_check<R: RngCore>(
        vk: &Self::VerifierKey,
        commitments: &[LabeledCommitment<Self::Commitment>],
        query_set: &QuerySet<E::Fr>,
        values: &Evaluations<E::Fr>,
        proof: &Self::BatchProof,
        opening_challenge: E::Fr,
        rng: &mut R,
    ) -> Result<bool, Self::Error> {
        let commitments: BTreeMap<_, _> = commitments.into_iter().map(|c| (c.label(), c)).collect();
        let mut query_to_labels_map = BTreeMap::new();

        for (label, point) in query_set.iter() {
            let labels = query_to_labels_map.entry(point).or_insert(BTreeSet::new());
            labels.insert(label);
        }
        assert_eq!(proof.proofs.len(), query_to_labels_map.len());

        let max_degree = vk.max_degree();
        let mut combined_comms = Vec::new();
        let mut combined_queries = Vec::new();
        let mut combined_evals = Vec::new();
        for (query, labels) in query_to_labels_map.into_iter() {
            let lc_time =
                start_timer!(|| format!("Randomly combining {} commitments", labels.len()));
            let mut comms_to_combine = Vec::new();
            let mut values_to_combine = Vec::new();
            for label in labels.into_iter() {
                let commitment = commitments.get(label).ok_or(Self::Error::IncorrectQuerySet(
                    "query set references commitment with incorrect label",
                ))?;
                let degree_bound = commitment.degree_bound();
                let commitment = commitment.commitment();
                assert_eq!(degree_bound.is_some(), commitment.shifted_comm.is_some());


                let v_i = values
                    .get(&(label.clone(), *query))
                    .ok_or(Self::Error::IncorrectEvaluation("value does not exist"))?;

                comms_to_combine.push(commitment);
                values_to_combine.push(*v_i);
            }
            let (c, v) = Self::accumulate_commitments_and_values(
                vk,
                comms_to_combine,
                query,
                values_to_combine,
                opening_challenge,
            )?;
            end_timer!(lc_time);
            combined_comms.push(c);
            combined_queries.push(*query);
            combined_evals.push(v);
        }
        let proof_time = start_timer!(|| "Checking SinglePC::Proof");
        let result = kzg10::KZG10::batch_check(
            vk,
            &combined_comms,
            &combined_queries,
            &combined_evals,
            &proof.proofs,
            rng,
        )?;
        end_timer!(proof_time);
        Ok(result)
    }
}

#[cfg(test)]
mod tests {
    #![allow(non_camel_case_types)]
    use crate::single_pc::kzg10::*;
    use crate::*;
    use algebra::fields::bls12_381::Fr;

    use algebra::curves::bls12_377::Bls12_377;
    use algebra::curves::bls12_381::Bls12_381;
    use algebra::curves::mnt6::MNT6;
    use algebra::curves::sw6::SW6;

    use rand::thread_rng;

    #[test]
    fn add_commitments_test() {
        let rng = &mut thread_rng();
        let p = Polynomial::from_coefficients_slice(&[
            Fr::rand(rng),
            Fr::rand(rng),
            Fr::rand(rng),
            Fr::rand(rng),
            Fr::rand(rng),
        ]);
        let f = Fr::rand(rng);
        let mut f_p = Polynomial::zero();
        f_p += (f, &p);

        let degree = 4;
        let pp = KZG_Bls12_381::setup(degree, rng).unwrap();
        let support = CoefficientSupport::from_dual_interval(4, 4);
        let (ck, _) = KZG_Bls12_381::trim(degree).unwrap();

        let hiding_bound = None;
        let (comm, _) = KZG10::commit(&ck, &p, hiding_bound, Some(rng)).unwrap();
        let (f_comm, _) = KZG10::commit(&ck, &f_p, hiding_bound, Some(rng)).unwrap();
        let mut f_comm_2 = Commitment::empty();
        f_comm_2 += (f, &comm);

        assert_eq!(f_comm, f_comm_2);
    }
    type KZG_Bls12_381 = KZG10<Bls12_381>;
    type KZG_Bls12_377 = KZG10<Bls12_377>;
    type KZG_MNT6 = KZG10<MNT6>;
    type KZG_SW6 = KZG10<SW6>;

    #[test]
    fn end_to_end_test() {
        use crate::single_pc::tests::*;

        end_to_end_test::<_, KZG_Bls12_377>().expect("test failed for bls12-377");
        end_to_end_test::<_, KZG_Bls12_381>().expect("test failed for bls12-381");
        end_to_end_test::<_, KZG_MNT6>().expect("test failed for MNT6");
        end_to_end_test::<_, KZG_SW6>().expect("test failed for SW6");
    }

    #[test]
    fn linear_polynomial_test() {
        use crate::single_pc::tests::*;

        linear_polynomial_test::<_, KZG_Bls12_377>().expect("test failed for bls12-377");
        linear_polynomial_test::<_, KZG_Bls12_381>().expect("test failed for bls12-381");
        linear_polynomial_test::<_, KZG_MNT6>().expect("test failed for MNT6");
        linear_polynomial_test::<_, KZG_SW6>().expect("test failed for SW6");
    }
    #[test]
    fn batch_check_test() {
        use crate::single_pc::tests::*;

        batch_check_test::<_, KZG_Bls12_377>().expect("test failed for bls12-377");
        batch_check_test::<_, KZG_Bls12_381>().expect("test failed for bls12-381");
        batch_check_test::<_, KZG_MNT6>().expect("test failed for MNT6");
        batch_check_test::<_, KZG_SW6>().expect("test failed for SW6");
    }
}
