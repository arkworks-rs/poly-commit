use crate::{PCCommitterKey, PCCommitment};
use crate::{BTreeMap, BTreeSet, ToString, Vec};
use crate::{BatchLCProof, Error, Evaluations, QuerySet};
use crate::{LabeledCommitment, LabeledPolynomial, LinearCombination};
use crate::{PCRandomness, PCUniversalParams, Polynomial, PolynomialCommitment};

use algebra_core::{One, PairingEngine, ProjectiveCurve, UniformRand, Zero, VariableBaseMSM, FixedBaseMSM, PrimeField, ToBytes};
use core::{convert::TryInto, marker::PhantomData};
use rand_core::RngCore;

mod data_structures;
pub use data_structures::*;

#[cfg(feature = "parallel")]
use rayon::prelude::*;

pub struct InnerProductArg<G: ProjectiveCurve> {
    _group: PhantomData<G>,
}

impl<G: ProjectiveCurve> InnerProductArg<G> {
    fn fr_to_bigints(
        p: &[G::ScalarField]
    ) -> Vec<<G::ScalarField>::BigInt> {
        ff_fft::cfg_iter!(p)
            .map(|s| s.into_repr())
            .collect()
    }

    fn dh_commit (
        comm_key: &[G::Affine],
        scalars: &[G::ScalarField]
    ) -> G {
        VariableBaseMSM::multi_scalar_mul(
            comm_key,
            InnerProductArg::fr_to_bigints(scalars).as_slice(),
        )
    }

    fn rand_oracle (
        challenge: impl ToBytes,
        L: G,
        R: G
    ) -> G::ScalarField {

    }

    fn inner_product(
        l: &[G::ScalarField],
        r: &[G::ScalarField]
    ) -> G::ScalarField {
        let prod = ff_fft::cfg_iter!(l)
            .zip(r)
            .map(|(li, ri)| li * ri)
            .sum();
        prod
    }

    fn split <T> (
        mut v: &[T]
    ) -> (&mut [T], &mut [T]) {
        let len = v.len()/2;
        v.split_at_mut(len)
    }
}

impl<G: ProjectiveCurve> PolynomialCommitment<G::ScalarField> for InnerProductArg<G> {
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
        // TODO: Replace these with deterministic random oracle
        let g = G::rand(rng);
        let h = G::rand(rng);

        let mut powers: Vec<G::ScalarField> = Vec::new();
        let mut curr_power: G::ScalarField = G::ScalarField::one();
        for _ in 0..=max_degree {
            powers.push(curr_power);
            curr_power += G::ScalarField::one();
        };

        let window_size = FixedBaseMSM::get_mul_window_size(max_degree + 1);
        let scalar_bits = G::ScalarField::size_in_bits();
        let g_table = FixedBaseMSM::get_window_table(scalar_bits, window_size, g);
        let comm_key = FixedBaseMSM::multi_scalar_mul::<G::Projective>(
            scalar_bits,
            window_size,
            &g_table,
            &powers,
        );

        let comm_key = G::batch_normalization_into_affine(comm_key.as_slice());

        let pp = UniversalParams {
            comm_key,
            h
        };

        Ok(pp)
    }

    fn trim(
        pp: &Self::UniversalParams,
        supported_degree: usize,
        enforced_degree_bounds: Option<&[usize]>,
    ) -> Result<(Self::CommitterKey, Self::VerifierKey), Self::Error> {
        //TODO: Check that supported degree is smaller than maximum degree
        let ck = CommitterKey {
            comm_key: pp.comm_key[0..=supported_degree].to_vec(),
            h: pp.h.copy(),
            max_degree: pp.max_degree()
        };

        let vk = VerifierKey {
            comm_key: pp.comm_key[0..=supported_degree].to_vec(),
            h: pp.h.copy(),
            max_degree: pp.max_degree()
        };

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
        let mut comms :Vec<LabeledCommitment<G>> = Vec::new();
        for labeled_polynomial in polynomials {
            // TODO: Check that degree of polynomial fits committer key

            let polynomial = labeled_polynomial.polynomial();
            let label = labeled_polynomial.label();
            let degree_bound = labeled_polynomial.degree_bound();

            let comm = Self::dh_commit(ck.comm_key.as_slice(), polynomial.coeffs.as_slice());
            let labeled_comm = LabeledCommitment::new(
                label.to_string(),
                Commitment(comm),
                degree_bound
            );

            comms.push(labeled_comm);
        }

        Ok((comms, Vec::new()))
    }

    fn open<'a>(
        ck: &Self::CommitterKey,
        labeled_polynomials: impl IntoIterator<Item = &'a LabeledPolynomial<'a, G::ScalarField>>,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Self::Commitment>>,
        point: G::ScalarField,
        opening_challenge: G::ScalarField,
        rands: impl IntoIterator<Item = &'a Self::Randomness>,
    ) -> Result<Self::Proof, Self::Error>
    where
        Self::Commitment: 'a,
        Self::Randomness: 'a,
    {
        let mut combined_polynomial = Polynomial::zero();
        let mut combined_commitment = Commitment::empty();

        let mut curr_challenge = opening_challenge;
        for (labeled_polynomial, labeled_commitment) in labeled_polynomials.zip(commitments) {
            combined_polynomial += (curr_challenge, labeled_polynomial.polynomial());
            combined_commitment += (curr_challenge, labeled_commitment.commitment());
            curr_challenge *= opening_challenge;
        }

        let mut combined_v = combined_polynomial.evaluate(point);

        // ith challenge
        let mut x = Self::rand_oracle(combined_commitment, z, );
        let h_prime = ck.h * x;

        // Polynomial coefficients
        let mut p = combined_polynomial.coeffs.as_mut_slice();

        // Comm key in both affine and projective form
        let mut key_aff = ck.comm_key.into_vec().as_mut_slice();
        let mut key_proj = ck.comm_key.iter().map(|x| x.into()).collect().as_mut_slice();

        // Powers of z
        let mut z: Vec<G::ScalarField> = Vec::with_capacity(p.len());
        let mut curr_z: G::ScalarField = G::ScalarField::one();
        for _ in 0..p.len() {
            z.push(curr_z);
            curr_z *= point;
        }
        let mut z = z.as_mut_slice();

        let mut L = Vec::new();
        let mut R = Vec::new();

        while p.len() > 1 {
            let (p_l, p_r) = Self::split(p);
            let (key_aff_l, key_aff_r) = Self::split(key_aff);
            let (z_l, z_r) = Self::split(z);

            let mut L_i = Self::dh_commit(key_aff_l, p_r);
            let mut R_i = Self::dh_commit(key_aff_r, p_l);

            L_i += h_prime.mul(Self::inner_product(p_r, z_l));
            R_i += h_prime.mul(Self::inner_product(p_l, z_r));

            L.push(L_i);
            R.push(R_i);

            x = random_oracle (x, L_i, R_i);
            let x_inverse = x.inverse();

            let (key_proj_l, key_proj_r):(G, G) = Self::split(key_proj);
            for i in 0..p_l.len() {
                p_l[i] += x_inverse * &p_r[i];
                z_l[i] += x * &z_r[i];
                key_proj_l[i] += x * &key_proj_r[i];
            }

            p = p_l;
            z = z_l;

            key_proj = key_proj_l;
            key_aff = G::batch_normalization_into_affine(key_proj).as_mut_slice();
        }

        let p = Proof {
            L,
            R,
            comm_key: key_aff[0].copy(),
            p: p[0].copy()
        };

        Ok(p)
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

    }
}

#[cfg(test)]
mod tests {
    #![allow(non_camel_case_types)]

    use crate::sonic_kzg10::SonicKZG10;
    use algebra::Bls12_377;
    use algebra::Bls12_381;
    use algebra::MNT6;
    use algebra::SW6;

    type PC<E> = SonicKZG10<E>;
    type PC_Bls12_377 = PC<Bls12_377>;
    type PC_Bls12_381 = PC<Bls12_381>;
    type PC_MNT6 = PC<MNT6>;
    type PC_SW6 = PC<SW6>;

    #[test]
    fn single_poly_test() {
        use crate::tests::*;
        single_poly_test::<_, PC_Bls12_377>().expect("test failed for bls12-377");
        single_poly_test::<_, PC_Bls12_381>().expect("test failed for bls12-381");
        single_poly_test::<_, PC_MNT6>().expect("test failed for MNT6");
        single_poly_test::<_, PC_SW6>().expect("test failed for SW6");
    }

    #[test]
    fn quadratic_poly_degree_bound_multiple_queries_test() {
        use crate::tests::*;
        quadratic_poly_degree_bound_multiple_queries_test::<_, PC_Bls12_377>()
            .expect("test failed for bls12-377");
        quadratic_poly_degree_bound_multiple_queries_test::<_, PC_Bls12_381>()
            .expect("test failed for bls12-381");
        quadratic_poly_degree_bound_multiple_queries_test::<_, PC_MNT6>()
            .expect("test failed for MNT6");
        quadratic_poly_degree_bound_multiple_queries_test::<_, PC_SW6>()
            .expect("test failed for SW6");
    }

    #[test]
    fn linear_poly_degree_bound_test() {
        use crate::tests::*;
        linear_poly_degree_bound_test::<_, PC_Bls12_377>().expect("test failed for bls12-377");
        linear_poly_degree_bound_test::<_, PC_Bls12_381>().expect("test failed for bls12-381");
        linear_poly_degree_bound_test::<_, PC_MNT6>().expect("test failed for MNT6");
        linear_poly_degree_bound_test::<_, PC_SW6>().expect("test failed for SW6");
    }

    #[test]
    fn single_poly_degree_bound_test() {
        use crate::tests::*;
        single_poly_degree_bound_test::<_, PC_Bls12_377>().expect("test failed for bls12-377");
        single_poly_degree_bound_test::<_, PC_Bls12_381>().expect("test failed for bls12-381");
        single_poly_degree_bound_test::<_, PC_MNT6>().expect("test failed for MNT6");
        single_poly_degree_bound_test::<_, PC_SW6>().expect("test failed for SW6");
    }

    #[test]
    fn single_poly_degree_bound_multiple_queries_test() {
        use crate::tests::*;
        single_poly_degree_bound_multiple_queries_test::<_, PC_Bls12_377>()
            .expect("test failed for bls12-377");
        single_poly_degree_bound_multiple_queries_test::<_, PC_Bls12_381>()
            .expect("test failed for bls12-381");
        single_poly_degree_bound_multiple_queries_test::<_, PC_MNT6>()
            .expect("test failed for MNT6");
        single_poly_degree_bound_multiple_queries_test::<_, PC_SW6>().expect("test failed for SW6");
    }

    #[test]
    fn two_polys_degree_bound_single_query_test() {
        use crate::tests::*;
        two_polys_degree_bound_single_query_test::<_, PC_Bls12_377>()
            .expect("test failed for bls12-377");
        two_polys_degree_bound_single_query_test::<_, PC_Bls12_381>()
            .expect("test failed for bls12-381");
        two_polys_degree_bound_single_query_test::<_, PC_MNT6>().expect("test failed for MNT6");
        two_polys_degree_bound_single_query_test::<_, PC_SW6>().expect("test failed for SW6");
    }

    #[test]
    fn full_end_to_end_test() {
        use crate::tests::*;
        full_end_to_end_test::<_, PC_Bls12_377>().expect("test failed for bls12-377");
        println!("Finished bls12-377");
        full_end_to_end_test::<_, PC_Bls12_381>().expect("test failed for bls12-381");
        println!("Finished bls12-381");
        full_end_to_end_test::<_, PC_MNT6>().expect("test failed for MNT6");
        println!("Finished mnt6");
        full_end_to_end_test::<_, PC_SW6>().expect("test failed for SW6");
        println!("Finished sw6");
    }

    #[test]
    fn single_equation_test() {
        use crate::tests::*;
        single_equation_test::<_, PC_Bls12_377>().expect("test failed for bls12-377");
        println!("Finished bls12-377");
        single_equation_test::<_, PC_Bls12_381>().expect("test failed for bls12-381");
        println!("Finished bls12-381");
        single_equation_test::<_, PC_MNT6>().expect("test failed for MNT6");
        println!("Finished mnt6");
        single_equation_test::<_, PC_SW6>().expect("test failed for SW6");
        println!("Finished sw6");
    }

    #[test]
    fn two_equation_test() {
        use crate::tests::*;
        two_equation_test::<_, PC_Bls12_377>().expect("test failed for bls12-377");
        println!("Finished bls12-377");
        two_equation_test::<_, PC_Bls12_381>().expect("test failed for bls12-381");
        println!("Finished bls12-381");
        two_equation_test::<_, PC_MNT6>().expect("test failed for MNT6");
        println!("Finished mnt6");
        two_equation_test::<_, PC_SW6>().expect("test failed for SW6");
        println!("Finished sw6");
    }

    #[test]
    fn two_equation_degree_bound_test() {
        use crate::tests::*;
        two_equation_degree_bound_test::<_, PC_Bls12_377>().expect("test failed for bls12-377");
        println!("Finished bls12-377");
        two_equation_degree_bound_test::<_, PC_Bls12_381>().expect("test failed for bls12-381");
        println!("Finished bls12-381");
        two_equation_degree_bound_test::<_, PC_MNT6>().expect("test failed for MNT6");
        println!("Finished mnt6");
        two_equation_degree_bound_test::<_, PC_SW6>().expect("test failed for SW6");
        println!("Finished sw6");
    }

    #[test]
    fn full_end_to_end_equation_test() {
        use crate::tests::*;
        full_end_to_end_equation_test::<_, PC_Bls12_377>().expect("test failed for bls12-377");
        println!("Finished bls12-377");
        full_end_to_end_equation_test::<_, PC_Bls12_381>().expect("test failed for bls12-381");
        println!("Finished bls12-381");
        full_end_to_end_equation_test::<_, PC_MNT6>().expect("test failed for MNT6");
        println!("Finished mnt6");
        full_end_to_end_equation_test::<_, PC_SW6>().expect("test failed for SW6");
        println!("Finished sw6");
    }
}
