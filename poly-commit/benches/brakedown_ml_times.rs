use ark_crypto_primitives::{
    crh::{sha256::Sha256, CRHScheme, TwoToOneCRHScheme},
    merkle_tree::{ByteDigestConverter, Config},
    sponge::poseidon::PoseidonSponge,
};
use ark_pcs_bench_templates::*;
use ark_poly::{DenseMultilinearExtension, MultilinearExtension};

use ark_bn254::Fr;
use ark_ff::PrimeField;

use ark_poly_commit::linear_codes::{LinearCodePCS, MultilinearBrakedown};
use blake2::Blake2s256;
use rand_chacha::ChaCha20Rng;

// Brakedown PCS over BN254
struct MerkleTreeParams;
type LeafH = LeafIdentityHasher;
type CompressH = Sha256;
impl Config for MerkleTreeParams {
    type Leaf = Vec<u8>;

    type LeafDigest = <LeafH as CRHScheme>::Output;
    type LeafInnerDigestConverter = ByteDigestConverter<Self::LeafDigest>;
    type InnerDigest = <CompressH as TwoToOneCRHScheme>::Output;

    type LeafHash = LeafH;
    type TwoToOneHash = CompressH;
}

pub type MLE<F> = DenseMultilinearExtension<F>;
type MTConfig = MerkleTreeParams;
type Sponge<F> = PoseidonSponge<F>;
type ColHasher<F> = FieldToBytesColHasher<F, Blake2s256>;
type Brakedown<F> = LinearCodePCS<
    MultilinearBrakedown<F, MTConfig, Sponge<F>, MLE<F>, ColHasher<F>>,
    F,
    MLE<F>,
    Sponge<F>,
    MTConfig,
    ColHasher<F>,
>;

fn rand_poly_brakedown_ml<F: PrimeField>(
    num_vars: usize,
    rng: &mut ChaCha20Rng,
) -> DenseMultilinearExtension<F> {
    DenseMultilinearExtension::rand(num_vars, rng)
}

fn rand_point_brakedown_ml<F: PrimeField>(num_vars: usize, rng: &mut ChaCha20Rng) -> Vec<F> {
    (0..num_vars).map(|_| F::rand(rng)).collect()
}

const MIN_NUM_VARS: usize = 12;
const MAX_NUM_VARS: usize = 22;

bench!(
    Brakedown<Fr>,
    rand_poly_brakedown_ml,
    rand_point_brakedown_ml
);
