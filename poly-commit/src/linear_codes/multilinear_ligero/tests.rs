#[cfg(test)]
mod tests {

    use crate::linear_codes::LinearCodePCS;
    use crate::utils::test_sponge;
    use crate::{
        linear_codes::{LigeroPCParams, MultilinearLigero, PolynomialCommitment},
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
    use ark_poly::evaluations::multivariate::{MultilinearExtension, SparseMultilinearExtension};
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

    type LigeroPCS<F> = LinearCodePCS<
        MultilinearLigero<
            F,
            MTConfig,
            Sponge<F>,
            SparseMultilinearExtension<F>,
            ColHasher<F, Blake2s256>,
        >,
        F,
        SparseMultilinearExtension<F>,
        Sponge<F>,
        MTConfig,
        ColHasher<F, Blake2s256>,
    >;

    fn rand_poly<Fr: PrimeField>(
        _: usize,
        num_vars: Option<usize>,
        rng: &mut ChaCha20Rng,
    ) -> SparseMultilinearExtension<Fr> {
        match num_vars {
            Some(n) => SparseMultilinearExtension::rand(n, rng),
            None => unimplemented!(), // should not happen in ML case!
        }
    }

    fn constant_poly<Fr: PrimeField>(
        _: usize,
        num_vars: Option<usize>,
        rng: &mut ChaCha20Rng,
    ) -> SparseMultilinearExtension<Fr> {
        // f1 = (1-x1)(1-x2)(1-x3)(1-x5)[(1-x6)*x4 + 2(1-x4)*x6]
        match num_vars {
            Some(n) => {
                let points = vec![(1, Fr::rand(rng))];
                SparseMultilinearExtension::from_evaluations(n, &points)
            }
            None => unimplemented!(), // should not happen in ML case!
        }
    }

    #[test]
    fn test_construction() {
        let mut rng = &mut test_rng();
        let num_vars = 10;
        // just to make sure we have the right degree given the FFT domain for our field
        let leaf_hash_param = <LeafH as CRHScheme>::setup(&mut rng).unwrap();
        let two_to_one_hash_param = <CompressH as TwoToOneCRHScheme>::setup(&mut rng)
            .unwrap()
            .clone();
        let col_hash_params = <ColHasher<Fr, Blake2s256> as CRHScheme>::setup(&mut rng).unwrap();
        let check_well_formedness = true;

        let pp: LigeroPCParams<Fr, MTConfig, ColHasher<Fr, Blake2s256>> = LigeroPCParams::new(
            128,
            4,
            check_well_formedness,
            leaf_hash_param,
            two_to_one_hash_param,
            col_hash_params,
        );

        let (ck, vk) = LigeroPCS::<Fr>::trim(&pp, 0, 0, None).unwrap();

        let rand_chacha = &mut ChaCha20Rng::from_rng(test_rng()).unwrap();
        let labeled_poly = LabeledPolynomial::new(
            "test".to_string(),
            rand_poly(1, Some(num_vars), rand_chacha),
            Some(num_vars),
            Some(num_vars),
        );

        let mut test_sponge = test_sponge::<Fr>();
        let (c, rands) = LigeroPCS::<Fr>::commit(&ck, &[labeled_poly.clone()], None).unwrap();

        let point = rand_point(Some(num_vars), rand_chacha);

        let value = labeled_poly.evaluate(&point);

        let proof = LigeroPCS::<Fr>::open(
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
            LigeroPCS::<Fr>::check(&vk, &c, &point, [value], &proof, &mut test_sponge, None)
                .unwrap()
        );
    }

    fn rand_point<F: Field>(num_vars: Option<usize>, rng: &mut ChaCha20Rng) -> Vec<F> {
        match num_vars {
            Some(n) => (0..n).map(|_| F::rand(rng)).collect(),
            None => unimplemented!(), // should not happen!
        }
    }

    #[test]
    fn single_poly_test() {
        use crate::tests::*;
        single_poly_test::<_, _, LigeroPCS<Fr>, _>(
            Some(5),
            rand_poly::<Fr>,
            rand_point::<Fr>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bls12-377");
        single_poly_test::<_, _, LigeroPCS<Fr381>, _>(
            Some(10),
            rand_poly::<Fr381>,
            rand_point::<Fr381>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bls12-381");
    }

    #[test]
    fn constant_poly_test() {
        use crate::tests::*;
        single_poly_test::<_, _, LigeroPCS<Fr>, _>(
            Some(10),
            constant_poly::<Fr>,
            rand_point::<Fr>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bls12-377");
        single_poly_test::<_, _, LigeroPCS<Fr381>, _>(
            Some(5),
            constant_poly::<Fr381>,
            rand_point::<Fr381>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bls12-381");
    }

    #[test]
    fn full_end_to_end_test() {
        use crate::tests::*;
        full_end_to_end_test::<_, _, LigeroPCS<Fr>, _>(
            Some(8),
            rand_poly::<Fr>,
            rand_point::<Fr>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bls12-377");
        println!("Finished bls12-377");
        full_end_to_end_test::<_, _, LigeroPCS<Fr381>, _>(
            Some(3),
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
        single_equation_test::<_, _, LigeroPCS<Fr>, _>(
            Some(10),
            rand_poly::<Fr>,
            rand_point::<Fr>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bls12-377");
        println!("Finished bls12-377");
        single_equation_test::<_, _, LigeroPCS<Fr381>, _>(
            Some(5),
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
        two_equation_test::<_, _, LigeroPCS<Fr>, _>(
            Some(5),
            rand_poly::<Fr>,
            rand_point::<Fr>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bls12-377");
        println!("Finished bls12-377");
        two_equation_test::<_, _, LigeroPCS<Fr381>, _>(
            Some(10),
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
        full_end_to_end_equation_test::<_, _, LigeroPCS<Fr>, _>(
            Some(5),
            rand_poly::<Fr>,
            rand_point::<Fr>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bls12-377");
        println!("Finished bls12-377");
        full_end_to_end_equation_test::<_, _, LigeroPCS<Fr381>, _>(
            Some(8),
            rand_poly::<Fr381>,
            rand_point::<Fr381>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bls12-381");
        println!("Finished bls12-381");
    }
}
