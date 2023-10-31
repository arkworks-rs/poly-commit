use ark_ec::AffineRepr;
use ark_pcs_bench_templates::*;
use ark_poly::DenseUVPolynomial;
use blake2::Blake2s256;

use ark_crypto_primitives::sponge::poseidon::PoseidonSponge;
use ark_ed_on_bls12_381::{EdwardsAffine, Fr};
use ark_ff::PrimeField;
use ark_poly::univariate::DensePolynomial as DenseUnivariatePoly;
use ark_poly_commit::ipa_pc::InnerProductArgPC;

use rand_chacha::ChaCha20Rng;

type UniPoly = DenseUnivariatePoly<Fr>;
type Sponge = PoseidonSponge<<EdwardsAffine as AffineRepr>::ScalarField>;

// IPA_PC over the JubJub curve with Blake2s as the hash function
#[allow(non_camel_case_types)]
type IPA_JubJub = InnerProductArgPC<EdwardsAffine, Blake2s256, UniPoly, Sponge>;

fn rand_poly_ipa_pc<F: PrimeField>(degree: usize, rng: &mut ChaCha20Rng) -> DenseUnivariatePoly<F> {
    DenseUnivariatePoly::rand(degree, rng)
}

fn rand_point_ipa_pc<F: PrimeField>(_: usize, rng: &mut ChaCha20Rng) -> F {
    F::rand(rng)
}

const MIN_NUM_VARS: usize = 10;
const MAX_NUM_VARS: usize = 20;

bench!(IPA_JubJub, rand_poly_ipa_pc, rand_point_ipa_pc);
