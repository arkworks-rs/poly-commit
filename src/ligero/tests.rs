#[cfg(test)]
mod tests {

    use crate::ark_std::UniformRand;
    use crate::{
        challenge::ChallengeGenerator,
        ligero::{
            utils::*, Ligero, LigeroPCCommitterKey, LigeroPCVerifierKey, PolynomialCommitment,
        },
        LabeledPolynomial,
    };
    use ark_bls12_377::Fq;
    use ark_crypto_primitives::sponge::poseidon::PoseidonConfig;
    use ark_crypto_primitives::sponge::CryptographicSponge;
    use ark_crypto_primitives::{
        crh::{pedersen, sha256::Sha256, CRHScheme, TwoToOneCRHScheme},
        merkle_tree::{ByteDigestConverter, Config},
        sponge::poseidon::PoseidonSponge,
    };
    use ark_ff::PrimeField;
    use ark_poly::{
        domain::general::GeneralEvaluationDomain, univariate::DensePolynomial, DenseUVPolynomial,
        EvaluationDomain, Polynomial,
    };
    use ark_std::test_rng;
    use blake2::Blake2s256;
    use rand_chacha::{rand_core::SeedableRng, ChaCha20Rng};

    type UniPoly = DensePolynomial<Fq>;
    #[derive(Clone)]
    pub(super) struct Window4x256;
    impl pedersen::Window for Window4x256 {
        const WINDOW_SIZE: usize = 4;
        const NUM_WINDOWS: usize = 256;
    }

    type LeafH = Sha256; //::CRH<JubJub, Window4x256>;
    type CompressH = Sha256; //pedersen::TwoToOneCRH<JubJub, Window4x256>;

    struct MerkleTreeParams;

    impl Config for MerkleTreeParams {
        type Leaf = [u8];
        // type Leaf = Vec<u8>;

        type LeafDigest = <LeafH as CRHScheme>::Output;
        type LeafInnerDigestConverter = ByteDigestConverter<Self::LeafDigest>;
        type InnerDigest = <CompressH as TwoToOneCRHScheme>::Output;

        type LeafHash = LeafH;
        type TwoToOneHash = CompressH;
    }

    type MTConfig = MerkleTreeParams;
    type Sponge = PoseidonSponge<Fq>;
    type PC<F, C, D, S, P, const rho_inv: usize> = Ligero<F, C, D, S, P, rho_inv, 128>;
    type LigeroPCS<const rho_inv: usize> = PC<Fq, MTConfig, Blake2s256, Sponge, UniPoly, rho_inv>;
    type LigeroPCS_F<const rho_inv: usize, F> =
        PC<F, MTConfig, Blake2s256, Sponge, DensePolynomial<F>, rho_inv>;

    #[test]
    fn test_matrix_constructor_flat() {
        let entries: Vec<Fq> = to_field(vec![10, 100, 4, 67, 44, 50]);
        let mat = Matrix::new_from_flat(2, 3, &entries);
        assert_eq!(mat.entry(1, 2), Fq::from(50));
    }

    #[test]
    fn test_matrix_constructor_flat_square() {
        let entries: Vec<Fq> = to_field(vec![10, 100, 4, 67]);
        let mat = Matrix::new_from_flat(2, 2, &entries);
        assert_eq!(mat.entry(1, 1), Fq::from(67));
    }

    #[test]
    #[should_panic]
    fn test_matrix_constructor_flat_panic() {
        let entries: Vec<Fq> = to_field(vec![10, 100, 4, 67, 44]);
        Matrix::new_from_flat(2, 3, &entries);
    }

    #[test]
    fn test_matrix_constructor_rows() {
        let rows: Vec<Vec<Fq>> = vec![
            to_field(vec![10, 100, 4]),
            to_field(vec![23, 1, 0]),
            to_field(vec![55, 58, 9]),
        ];
        let mat = Matrix::new_from_rows(rows);
        assert_eq!(mat.entry(2, 0), Fq::from(55));
    }

    #[test]
    #[should_panic]
    fn test_matrix_constructor_rows_panic() {
        let rows: Vec<Vec<Fq>> = vec![
            to_field(vec![10, 100, 4]),
            to_field(vec![23, 1, 0]),
            to_field(vec![55, 58]),
        ];
        Matrix::new_from_rows(rows);
    }

    #[test]
    fn test_cols() {
        let rows: Vec<Vec<Fq>> = vec![
            to_field(vec![4, 76]),
            to_field(vec![14, 92]),
            to_field(vec![17, 89]),
        ];

        let mat = Matrix::new_from_rows(rows);

        assert_eq!(mat.cols()[1], to_field(vec![76, 92, 89]));
    }

    #[test]
    fn test_row_mul() {
        let rows: Vec<Vec<Fq>> = vec![
            to_field(vec![10, 100, 4]),
            to_field(vec![23, 1, 0]),
            to_field(vec![55, 58, 9]),
        ];

        let mat = Matrix::new_from_rows(rows);
        let v: Vec<Fq> = to_field(vec![12, 41, 55]);
        // by giving the result in the integers and then converting to Fq
        // we ensure the test will still pass even if Fq changes
        assert_eq!(mat.row_mul(&v), to_field::<Fq>(vec![4088, 4431, 543]));
    }

    #[test]
    fn test_reed_solomon() {
        // we use this polynomial to generate the the values we will ask the fft to interpolate

        let rho_inv = 3;
        // `i` is the min number of evaluations we need to interpolate a poly of degree `i - 1`
        for i in 1..10 {
            let rand_chacha = &mut ChaCha20Rng::from_rng(test_rng()).unwrap();
            let pol = rand_poly::<Fq>(i - 1, None, rand_chacha);

            let coeffs = &pol.coeffs;
            assert_eq!(
                pol.degree(),
                coeffs.len() - 1,
                "degree of poly and coeffs mismatch"
            );
            assert_eq!(coeffs.len(), i, "length of coeffs and m mismatch");

            let small_domain = GeneralEvaluationDomain::<Fq>::new(i).unwrap();

            // size of evals might be larger than i (the min. number of evals needed to interpolate): we could still do R-S encoding on smaller evals, but the resulting polynomial will differ, so for this test to work we should pass it in full
            let evals = small_domain.fft(&coeffs);
            let m = evals.len();

            let coeffs_again = small_domain.ifft(&evals);
            assert_eq!(coeffs_again[..i], *coeffs);

            let encoded = reed_solomon(&evals, rho_inv);
            // Assert that the encoded vector is of the right length
            assert_eq!(encoded.len(), rho_inv * m);

            // first elements of encoded should be itself, since the code is systematic
            assert_eq!(encoded[..m], evals);

            let large_domain = GeneralEvaluationDomain::<Fq>::new(m * (rho_inv - 1)).unwrap();

            // The rest of the elements should agree with the domain
            for j in 0..((rho_inv - 1) * m) {
                assert_eq!(pol.evaluate(&large_domain.element(j)), encoded[j + m]);
            }
        }
    }

    #[test]
    fn test_fft_interface() {
        // This test is probably too verbose, and should be merged with the RS test
        let rho = 2;
        for m in 1..10 {
            // rand poly of degree m
            let domain = GeneralEvaluationDomain::<Fq>::new(m).unwrap();
            let poly = UniPoly::rand(m - 1, &mut test_rng());
            // get its evaluations at the entire domain
            let evals = (0..domain.size())
                .map(|i| poly.evaluate(&domain.element(i)))
                .collect::<Vec<_>>();

            // convert back to the coefficients
            let coeffs = domain.ifft(&evals);
            assert_eq!(coeffs[..m], poly.coeffs);

            let evals2 = domain.fft(&coeffs.to_vec());
            assert_eq!(evals[..m], evals2[..m]);

            // now we try with a larger domain
            let large_domain = GeneralEvaluationDomain::<Fq>::new(m * rho).unwrap();

            let evals3 = large_domain.fft(&coeffs.to_vec());
            let evals4: Vec<_> = (0..large_domain.size())
                .map(|i| poly.evaluate(&large_domain.element(i)))
                .collect::<Vec<_>>();

            assert_eq!(evals3[..m], evals4[..m]);

            let coeffs2 = large_domain.ifft(&evals3);
        }
    }

    #[test]
    fn test_merkle_tree() {
        let mut rng = &mut test_rng();
        let leaf_hash_params = <LeafH as CRHScheme>::setup(&mut rng).unwrap();
        let two_to_one_params = <CompressH as TwoToOneCRHScheme>::setup(&mut rng)
            .unwrap()
            .clone();

        let rows: Vec<Vec<Fq>> = vec![
            to_field(vec![4, 76]),
            to_field(vec![14, 92]),
            to_field(vec![17, 89]),
        ];

        let mat = Matrix::new_from_rows(rows);
        let mt = LigeroPCS::<2>::create_merkle_tree(&mat, &leaf_hash_params, &two_to_one_params);

        let root = mt.root();

        for (i, col) in mat.cols().iter().enumerate() {
            let col_hash = hash_column::<Blake2s256, Fq>(col);

            let proof = mt.generate_proof(i).unwrap();
            assert!(proof
                .verify(
                    &leaf_hash_params,
                    &two_to_one_params,
                    &root,
                    col_hash.clone()
                )
                .unwrap());
        }
    }

    #[test]
    fn test_get_num_bytes() {
        assert_eq!(get_num_bytes(0), 0);
        assert_eq!(get_num_bytes(1), 1);
        assert_eq!(get_num_bytes(9), 1);
        assert_eq!(get_num_bytes(1 << 11), 2);
        assert_eq!(get_num_bytes(1 << 32 - 1), 4);
        assert_eq!(get_num_bytes(1 << 32), 5);
        assert_eq!(get_num_bytes(1 << 32 + 1), 5);
    }

    fn rand_poly<Fq: PrimeField>(
        degree: usize,
        _: Option<usize>,
        rng: &mut ChaCha20Rng,
    ) -> DensePolynomial<Fq> {
        DensePolynomial::rand(degree, rng)
    }

    fn constant_poly<Fq: PrimeField>(
        _: usize,
        _: Option<usize>,
        rng: &mut ChaCha20Rng,
    ) -> DensePolynomial<Fq> {
        DensePolynomial::from_coefficients_slice(&[Fq::rand(rng)])
    }

    // TODO: replace by https://github.com/arkworks-rs/crypto-primitives/issues/112.
    fn test_sponge<F: PrimeField>() -> PoseidonSponge<F> {
        let full_rounds = 8;
        let partial_rounds = 31;
        let alpha = 17;

        let mds = vec![
            vec![F::one(), F::zero(), F::one()],
            vec![F::one(), F::one(), F::zero()],
            vec![F::zero(), F::one(), F::one()],
        ];

        let mut v = Vec::new();
        let mut ark_rng = test_rng();

        for _ in 0..(full_rounds + partial_rounds) {
            let mut res = Vec::new();

            for _ in 0..3 {
                res.push(F::rand(&mut ark_rng));
            }
            v.push(res);
        }
        let config = PoseidonConfig::new(full_rounds, partial_rounds, alpha, mds, v, 2, 1);
        PoseidonSponge::new(&config)
    }

    #[test]
    fn test_setup() {
        let rng = &mut test_rng();
        let _ = LigeroPCS::<2>::setup(1 << 44, None, rng).unwrap();

        assert_eq!(LigeroPCS::<5>::setup(1 << 45, None, rng).is_err(), true);

        // but the base field of bls12_381 doesnt have such large domains
        use ark_bls12_381::Fq as F_381;
        assert_eq!(LigeroPCS_F::<5, F_381>::setup(10, None, rng).is_err(), true);
    }

    #[test]
    fn test_construction() {
        let degree = 4;
        let mut rng = &mut test_rng();
        // just to make sure we have the right degree given the FFT domain for our field
        LigeroPCS::<2>::setup(degree, None, rng).unwrap();
        let leaf_hash_params = <LeafH as CRHScheme>::setup(&mut rng).unwrap();
        let two_to_one_params = <CompressH as TwoToOneCRHScheme>::setup(&mut rng)
            .unwrap()
            .clone();

        let ck: LigeroPCCommitterKey<MTConfig> = LigeroPCCommitterKey {
            leaf_hash_params,
            two_to_one_params,
        };
        let vk: LigeroPCVerifierKey<MTConfig> = LigeroPCVerifierKey {
            leaf_hash_params,
            two_to_one_params,
        };

        let rand_chacha = &mut ChaCha20Rng::from_rng(test_rng()).unwrap();
        let labeled_poly = LabeledPolynomial::new(
            "test".to_string(),
            rand_poly(degree, None, rand_chacha),
            None,
            None,
        );

        let mut test_sponge = test_sponge::<Fq>();
        let (c, rands) = LigeroPCS::<2>::commit(&ck, &[labeled_poly.clone()], None).unwrap();

        let point = Fq::rand(rand_chacha);

        let mut challenge_generator: ChallengeGenerator<Fq, PoseidonSponge<Fq>> =
            ChallengeGenerator::new_univariate(&mut test_sponge);

        let proof = LigeroPCS::<2>::open(
            &ck,
            &[labeled_poly],
            &c,
            &point,
            &mut challenge_generator,
            &rands,
            None,
        )
        .unwrap();
    }
}
