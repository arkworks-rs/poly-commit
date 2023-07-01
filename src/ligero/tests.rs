#[cfg(test)]
mod tests {

    use ark_bls12_381::Fq as F;
    use ark_ff::PrimeField;
    use ark_poly::{
        domain::general::GeneralEvaluationDomain, univariate::DensePolynomial, DenseUVPolynomial,
        EvaluationDomain, Polynomial,
    };
    use ark_std::test_rng;
    use blake2::Blake2s256;
    use rand_chacha::{rand_core::SeedableRng, ChaCha20Rng};

    use crate::{
        ligero::{
            utils::*, Ligero, LigeroPCCommitterKey, LigeroPCVerifierKey, PolynomialCommitment,
        },
        LabeledPolynomial,
    };
    use ark_crypto_primitives::{
        crh::{pedersen, sha256::Sha256, CRHScheme, TwoToOneCRHScheme},
        merkle_tree::{ByteDigestConverter, Config},
        sponge::poseidon::PoseidonSponge,
    };

    type UniPoly = DensePolynomial<F>;
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
    type Sponge = PoseidonSponge<F>;
    type PC<F, C, D, S, P, const rho_inv: usize> = Ligero<F, C, D, S, P, rho_inv, 128>;
    type LigeroPCS<const rho_inv: usize> = PC<F, MTConfig, Blake2s256, Sponge, UniPoly, rho_inv>;

    #[test]
    fn test_matrix_constructor_flat() {
        let entries: Vec<F> = to_field(vec![10, 100, 4, 67, 44, 50]);
        let mat = Matrix::new_from_flat(2, 3, &entries);
        assert_eq!(mat.entry(1, 2), F::from(50));
    }

    #[test]
    fn test_matrix_constructor_flat_square() {
        let entries: Vec<F> = to_field(vec![10, 100, 4, 67]);
        let mat = Matrix::new_from_flat(2, 2, &entries);
        assert_eq!(mat.entry(1, 1), F::from(67));
    }

    #[test]
    #[should_panic]
    fn test_matrix_constructor_flat_panic() {
        let entries: Vec<F> = to_field(vec![10, 100, 4, 67, 44]);
        Matrix::new_from_flat(2, 3, &entries);
    }

    #[test]
    fn test_matrix_constructor_rows() {
        let rows: Vec<Vec<F>> = vec![
            to_field(vec![10, 100, 4]),
            to_field(vec![23, 1, 0]),
            to_field(vec![55, 58, 9]),
        ];
        let mat = Matrix::new_from_rows(rows);
        assert_eq!(mat.entry(2, 0), F::from(55));
    }

    #[test]
    #[should_panic]
    fn test_matrix_constructor_rows_panic() {
        let rows: Vec<Vec<F>> = vec![
            to_field(vec![10, 100, 4]),
            to_field(vec![23, 1, 0]),
            to_field(vec![55, 58]),
        ];
        Matrix::new_from_rows(rows);
    }

    #[test]
    fn test_cols() {
        let rows: Vec<Vec<F>> = vec![
            to_field(vec![4, 76]),
            to_field(vec![14, 92]),
            to_field(vec![17, 89]),
        ];

        let mat = Matrix::new_from_rows(rows);

        assert_eq!(mat.cols()[1], to_field(vec![76, 92, 89]));
    }

    #[test]
    fn test_row_mul() {
        let rows: Vec<Vec<F>> = vec![
            to_field(vec![10, 100, 4]),
            to_field(vec![23, 1, 0]),
            to_field(vec![55, 58, 9]),
        ];

        let mat = Matrix::new_from_rows(rows);
        let v: Vec<F> = to_field(vec![12, 41, 55]);
        // by giving the result in the integers and then converting to F
        // we ensure the test will still pass even if F changes
        assert_eq!(mat.row_mul(&v), to_field::<F>(vec![4088, 4431, 543]));
    }

    #[test]
    fn test_reed_solomon() {
        // we use this polynomial to generate the the values we will ask the fft to interpolate
        // let pol_evals: Vec<F> = to_field(vec![30, 2, 91, 4, 8]);
        // TODO try for different values of m
        // for i in 0..16

        let rho_inv = 3;
        // `i` is the min number of evaluations we need to interpolate a poly of degree `i - 1`
        for i in 1..10 {
            let rand_chacha = &mut ChaCha20Rng::from_rng(test_rng()).unwrap();
            let pol = rand_poly::<F>(i - 1, None, rand_chacha);

            let coeffs = &pol.coeffs;
            assert_eq!(
                pol.degree(),
                coeffs.len() - 1,
                "degree of poly and coeffs mismatch"
            );
            assert_eq!(coeffs.len(), i, "length of coeffs and m mismatch");

            let small_domain = GeneralEvaluationDomain::<F>::new(i).unwrap();

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

            let large_domain = GeneralEvaluationDomain::<F>::new(m * (rho_inv - 1)).unwrap();

            // The rest of the elements should agree with the domain
            for j in 0..((rho_inv - 1) * m) {
                println!("j: {:?}", j);
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
            let domain = GeneralEvaluationDomain::<F>::new(m).unwrap();
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
            let large_domain = GeneralEvaluationDomain::<F>::new(m * rho).unwrap();

            let evals3 = large_domain.fft(&coeffs.to_vec());
            let evals4: Vec<_> = (0..large_domain.size())
                .map(|i| poly.evaluate(&large_domain.element(i)))
                .collect::<Vec<_>>();

            assert_eq!(evals3[..m], evals4[..m]);

            let coeffs2 = large_domain.ifft(&evals3);
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

    fn rand_poly<F: PrimeField>(
        degree: usize,
        _: Option<usize>,
        rng: &mut ChaCha20Rng,
    ) -> DensePolynomial<F> {
        DensePolynomial::rand(degree, rng)
    }

    fn constant_poly<F: PrimeField>(
        _: usize,
        _: Option<usize>,
        rng: &mut ChaCha20Rng,
    ) -> DensePolynomial<F> {
        DensePolynomial::from_coefficients_slice(&[F::rand(rng)])
    }

    #[test]
    fn test_setup() {
        let mut rng = &mut test_rng();
        let _ = LigeroPCS::<2>::setup(10, None, rng).unwrap();

        // the field we use doesnt have such large domains
        assert_eq!(LigeroPCS::<5>::setup(10, None, rng).is_err(), true);
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

        let c = LigeroPCS::<2>::commit(&ck, &[labeled_poly], None).unwrap();
    }
}
