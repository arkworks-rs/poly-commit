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
type PC<E, D, P, S> = InnerProductArgPC<E, D, P, S>;

// IPA_PC over the JubJub curve with Blake2s as the hash function
#[allow(non_camel_case_types)]
type IPA_JubJub = PC<EdwardsAffine, Blake2s256, UniPoly, Sponge>;

fn rand_poly_ipa_pc<F: PrimeField>(degree: usize, rng: &mut ChaCha20Rng) -> DenseUnivariatePoly<F> {
    DenseUnivariatePoly::rand(degree, rng)
}

const MIN_NUM_VARS: usize = 10;
const MAX_NUM_VARS: usize = 20;

fn main() {
    println!("\nIPA on JubJub: Commitment size");
    for num_vars in (MIN_NUM_VARS..MAX_NUM_VARS).step_by(2) {
        println!(
            "\tnum_vars: {}, size: {} B",
            num_vars,
            commitment_size::<_, _, IPA_JubJub>(num_vars, rand_poly_ipa_pc)
        );
    }

    println!("\nIPA on JubJub: Proof size");
    for num_vars in (MIN_NUM_VARS..MAX_NUM_VARS).step_by(2) {
        println!(
            "\tnum_vars: {}, size: {} B",
            num_vars,
            proof_size::<_, _, IPA_JubJub>(num_vars, rand_poly_ipa_pc)
        );
    }
}
