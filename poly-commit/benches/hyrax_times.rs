use ark_crypto_primitives::sponge::poseidon::PoseidonSponge;
use ark_pcs_bench_templates::*;
use ark_poly::{DenseMultilinearExtension, MultilinearExtension};

use ark_bn254::{Fr, G1Affine};
use ark_ff::PrimeField;
use ark_poly_commit::hyrax::HyraxPC;

use rand_chacha::ChaCha20Rng;

// Hyrax PCS over BN254
type Hyrax254 = HyraxPC<G1Affine, DenseMultilinearExtension<Fr>, PoseidonSponge<Fr>>;

fn rand_poly_hyrax<F: PrimeField>(
    num_vars: usize,
    rng: &mut ChaCha20Rng,
) -> DenseMultilinearExtension<F> {
    DenseMultilinearExtension::rand(num_vars, rng)
}

fn rand_point_hyrax<F: PrimeField>(num_vars: usize, rng: &mut ChaCha20Rng) -> Vec<F> {
    (0..num_vars).map(|_| F::rand(rng)).collect()
}

const MIN_NUM_VARS: usize = 12;
const MAX_NUM_VARS: usize = 22;

bench!(Hyrax254, rand_poly_hyrax, rand_point_hyrax);
