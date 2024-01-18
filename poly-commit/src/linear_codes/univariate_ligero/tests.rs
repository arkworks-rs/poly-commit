#[cfg(test)]
mod tests {

    use crate::ark_std::UniformRand;
    use crate::linear_codes::LinearCodePCS;
    use crate::utils::test_sponge;
    use crate::{
        linear_codes::{LigeroPCParams, PolynomialCommitment, UnivariateLigero},
        LabeledPolynomial,
    };
    use ark_bls12_377::Fr;
    use ark_bls12_381::Fr as Fr381;
    use ark_crypto_primitives::{
        crh::{sha256::Sha256, CRHScheme, TwoToOneCRHScheme},
        merkle_tree::{ByteDigestConverter, Config},
        sponge::poseidon::PoseidonSponge,
    };
    use ark_ff::{Field, PrimeField};
    use ark_poly::{univariate::DensePolynomial, DenseUVPolynomial};
    use ark_std::test_rng;
    use blake2::Blake2s256;
    use rand_chacha::{rand_core::SeedableRng, ChaCha20Rng};

    use ark_pcs_bench_templates::{FieldToBytesColHasher, LeafIdentityHasher};

    type LeafH = LeafIdentityHasher;
    type CompressH = Sha256;
    type ColHasher<F, D> = FieldToBytesColHasher<F, D>;

    struct MerkleTreeParams;

    impl Config for MerkleTreeParams {
        type Leaf = Vec<u8>;

        type LeafDigest = <LeafH as CRHScheme>::Output;
        type LeafInnerDigestConverter = ByteDigestConverter<Self::LeafDigest>;
        type InnerDigest = <CompressH as TwoToOneCRHScheme>::Output;

        type LeafHash = LeafH;
        type TwoToOneHash = CompressH;
    }

    type MTConfig = MerkleTreeParams;
    type Sponge<F> = PoseidonSponge<F>;

    type LigeroPCS = LinearCodePCS<
        UnivariateLigero<Fr, MTConfig, Sponge<Fr>, DensePolynomial<Fr>, ColHasher<Fr, Blake2s256>>,
        Fr,
        DensePolynomial<Fr>,
        Sponge<Fr>,
        MTConfig,
        ColHasher<Fr, Blake2s256>,
    >;

    type LigeroPcsF<F> = LinearCodePCS<
        UnivariateLigero<F, MTConfig, Sponge<F>, DensePolynomial<F>, ColHasher<F, Blake2s256>>,
        F,
        DensePolynomial<F>,
        Sponge<F>,
        MTConfig,
        ColHasher<F, Blake2s256>,
    >;

    fn rand_poly<Fr: PrimeField>(
        degree: usize,
        _: Option<usize>,
        rng: &mut ChaCha20Rng,
    ) -> DensePolynomial<Fr> {
        DensePolynomial::rand(degree, rng)
    }

    fn constant_poly<Fr: PrimeField>(
        _: usize,
        _: Option<usize>,
        rng: &mut ChaCha20Rng,
    ) -> DensePolynomial<Fr> {
        DensePolynomial::from_coefficients_slice(&[Fr::rand(rng)])
    }

    #[test]
    fn test_construction() {
        let degree = 4;
        let mut rng = &mut test_rng();
        // just to make sure we have the right degree given the FFT domain for our field
        let leaf_hash_param = <LeafH as CRHScheme>::setup(&mut rng).unwrap();
        let two_to_one_hash_param = <CompressH as TwoToOneCRHScheme>::setup(&mut rng)
            .unwrap()
            .clone();
        let col_hash_params = <ColHasher<Fr, Blake2s256> as CRHScheme>::setup(&mut rng).unwrap();
        let check_well_formedness = true;

        let pp: LigeroPCParams<Fr, MTConfig, ColHasher<_, _>> = LigeroPCParams::new(
            128,
            4,
            check_well_formedness,
            leaf_hash_param,
            two_to_one_hash_param,
            col_hash_params,
        );

        let (ck, vk) = LigeroPCS::trim(&pp, 0, 0, None).unwrap();

        let rand_chacha = &mut ChaCha20Rng::from_rng(test_rng()).unwrap();
        let labeled_poly = LabeledPolynomial::new(
            "test".to_string(),
            rand_poly(degree, None, rand_chacha),
            None,
            None,
        );

        let mut test_sponge = test_sponge::<Fr>();
        let (c, rands) = LigeroPCS::commit(&ck, &[labeled_poly.clone()], None).unwrap();

        let point = Fr::rand(rand_chacha);

        let value = labeled_poly.evaluate(&point);

        let proof = LigeroPCS::open(
            &ck,
            &[labeled_poly],
            &c,
            &point,
            &mut (test_sponge.clone()),
            &rands,
            None,
        )
        .unwrap();
        assert!(
            LigeroPCS::check(&vk, &c, &point, [value], &proof, &mut test_sponge, None).unwrap()
        );
    }

    fn rand_point<F: Field>(_: Option<usize>, rng: &mut ChaCha20Rng) -> F {
        F::rand(rng)
    }

    #[test]
    fn single_poly_test() {
        use crate::tests::*;
        single_poly_test::<_, _, LigeroPCS, _>(
            None,
            rand_poly::<Fr>,
            rand_point::<Fr>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bls12-377");
        single_poly_test::<_, _, LigeroPcsF<Fr381>, _>(
            None,
            rand_poly::<Fr381>,
            rand_point::<Fr381>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bls12-381");
    }

    #[test]
    fn constant_poly_test() {
        use crate::tests::*;
        single_poly_test::<_, _, LigeroPCS, _>(
            None,
            constant_poly::<Fr>,
            rand_point::<Fr>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bls12-377");
        single_poly_test::<_, _, LigeroPcsF<Fr381>, _>(
            None,
            constant_poly::<Fr381>,
            rand_point::<Fr381>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bls12-381");
    }

    #[test]
    fn quadratic_poly_degree_bound_multiple_queries_test() {
        use crate::tests::*;
        quadratic_poly_degree_bound_multiple_queries_test::<_, _, LigeroPCS, _>(
            rand_poly::<Fr>,
            rand_point::<Fr>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bls12-377");
        quadratic_poly_degree_bound_multiple_queries_test::<_, _, LigeroPcsF<Fr381>, _>(
            rand_poly::<Fr381>,
            rand_point::<Fr381>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bls12-381");
    }

    #[test]
    fn linear_poly_degree_bound_test() {
        use crate::tests::*;
        linear_poly_degree_bound_test::<_, _, LigeroPCS, _>(
            rand_poly::<Fr>,
            rand_point::<Fr>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bls12-377");
        linear_poly_degree_bound_test::<_, _, LigeroPcsF<Fr381>, _>(
            rand_poly::<Fr381>,
            rand_point::<Fr381>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bls12-381");
    }

    #[test]
    fn single_poly_degree_bound_test() {
        use crate::tests::*;
        single_poly_degree_bound_test::<_, _, LigeroPCS, _>(
            rand_poly::<Fr>,
            rand_point::<Fr>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bls12-377");
        single_poly_degree_bound_test::<_, _, LigeroPcsF<Fr381>, _>(
            rand_poly::<Fr381>,
            rand_point::<Fr381>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bls12-381");
    }

    #[test]
    fn single_poly_degree_bound_multiple_queries_test() {
        use crate::tests::*;
        single_poly_degree_bound_multiple_queries_test::<_, _, LigeroPCS, _>(
            rand_poly::<Fr>,
            rand_point::<Fr>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bls12-377");
        single_poly_degree_bound_multiple_queries_test::<_, _, LigeroPcsF<Fr381>, _>(
            rand_poly::<Fr381>,
            rand_point::<Fr381>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bls12-381");
    }

    #[test]
    fn two_polys_degree_bound_single_query_test() {
        use crate::tests::*;
        two_polys_degree_bound_single_query_test::<_, _, LigeroPCS, _>(
            rand_poly::<Fr>,
            rand_point::<Fr>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bls12-377");
        two_polys_degree_bound_single_query_test::<_, _, LigeroPcsF<Fr381>, _>(
            rand_poly::<Fr381>,
            rand_point::<Fr381>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bls12-381");
    }

    #[test]
    fn full_end_to_end_test() {
        use crate::tests::*;
        full_end_to_end_test::<_, _, LigeroPCS, _>(
            None,
            rand_poly::<Fr>,
            rand_point::<Fr>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bls12-377");
        println!("Finished bls12-377");
        full_end_to_end_test::<_, _, LigeroPcsF<Fr381>, _>(
            None,
            rand_poly::<Fr381>,
            rand_point::<Fr381>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bls12-381");
        println!("Finished bls12-381");
    }

    #[test]
    fn single_equation_test() {
        use crate::tests::*;
        single_equation_test::<_, _, LigeroPCS, _>(
            None,
            rand_poly::<Fr>,
            rand_point::<Fr>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bls12-377");
        println!("Finished bls12-377");
        single_equation_test::<_, _, LigeroPcsF<Fr381>, _>(
            None,
            rand_poly::<Fr381>,
            rand_point::<Fr381>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bls12-381");
        println!("Finished bls12-381");
    }

    #[test]
    fn two_equation_test() {
        use crate::tests::*;
        two_equation_test::<_, _, LigeroPCS, _>(
            None,
            rand_poly::<Fr>,
            rand_point::<Fr>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bls12-377");
        println!("Finished bls12-377");
        two_equation_test::<_, _, LigeroPcsF<Fr381>, _>(
            None,
            rand_poly::<Fr381>,
            rand_point::<Fr381>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bls12-381");
        println!("Finished bls12-381");
    }

    #[test]
    fn two_equation_degree_bound_test() {
        use crate::tests::*;
        two_equation_degree_bound_test::<_, _, LigeroPCS, _>(
            rand_poly::<Fr>,
            rand_point::<Fr>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bls12-377");
        println!("Finished bls12-377");
        two_equation_degree_bound_test::<_, _, LigeroPcsF<Fr381>, _>(
            rand_poly::<Fr381>,
            rand_point::<Fr381>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bls12-381");
        println!("Finished bls12-381");
    }

    #[test]
    fn full_end_to_end_equation_test() {
        use crate::tests::*;
        full_end_to_end_equation_test::<_, _, LigeroPCS, _>(
            None,
            rand_poly::<Fr>,
            rand_point::<Fr>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bls12-377");
        println!("Finished bls12-377");
        full_end_to_end_equation_test::<_, _, LigeroPcsF<Fr381>, _>(
            None,
            rand_poly::<Fr381>,
            rand_point::<Fr381>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bls12-381");
        println!("Finished bls12-381");
    }

    #[test]
    #[should_panic]
    fn bad_degree_bound_test() {
        use crate::tests::*;
        use ark_bls12_381::Fq as Fq381;
        bad_degree_bound_test::<_, _, LigeroPcsF<Fq381>, _>(
            rand_poly::<Fq381>,
            rand_point::<Fq381>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bls12-377");
        println!("Finished bls12-377");
    }
}
