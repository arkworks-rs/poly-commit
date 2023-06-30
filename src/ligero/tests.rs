#[cfg(test)]
mod tests {

    use ark_bls12_381::Fq as F;
    use ark_poly::{
        domain::general::GeneralEvaluationDomain, univariate::DensePolynomial, DenseUVPolynomial,
        EvaluationDomain, Polynomial,
    };
    use ark_std::test_rng;
    use blake2::Blake2s256;

    use crate::ligero::{utils::*, Ligero, PolynomialCommitment};
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
    type PC<F, C, D, S, P> = Ligero<F, C, D, S, P, 2, 128>;
    type LigeroPCS = PC<F, MTConfig, Blake2s256, Sponge, UniPoly>;

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
    fn test_fft_interface() {
        // we use this polynomial to generate the the values we will ask the fft to interpolate
        let pol_coeffs: Vec<F> = to_field(vec![30, 2, 91]);
        let pol: DensePolynomial<F> = DensePolynomial::from_coefficients_slice(&pol_coeffs);

        let fft_domain = GeneralEvaluationDomain::<F>::new(pol_coeffs.len()).unwrap();

        // generating the values
        let mut vals = Vec::new();

        for i in 0..4 {
            vals.push(pol.evaluate(&fft_domain.element(i)));
        }

        // the fft should recover the original polynomial
        assert_eq!(fft_domain.ifft(&vals), pol_coeffs);
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

    #[test]
    fn test_construction() {
        let rng = &mut test_rng();
        let pp = LigeroPCS::setup(2, None, rng).unwrap();
        // This fails since trim is not implemented
        let (ck, vk) = LigeroPCS::trim(&pp, 0, 2, Some(&[0])).unwrap();
    }
}
