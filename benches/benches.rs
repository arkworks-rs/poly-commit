// #![cfg(feature = "benches")]
use ark_ec::AffineRepr;
use ark_poly::DenseUVPolynomial;
use blake2::Blake2s256;
use criterion::{criterion_group, criterion_main, Criterion};

use ark_crypto_primitives::sponge::poseidon::PoseidonSponge;
use ark_ed_on_bls12_381::{EdwardsAffine, Fr};
use ark_ff::PrimeField;
use ark_poly::univariate::DensePolynomial as DenseUnivariatePoly;
use ark_poly_commit::bench_templates::*;
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

macro_rules! bench_pcs {
    ($c:expr, $method:ident, $scheme_type:ty, $rand_poly:ident) => {{
        let scheme_type_str = stringify!($scheme_type);
        let bench_name = format!("{} {}", stringify!($method), scheme_type_str);
        bench_pcs_method::<_, _, $scheme_type>(
            $c,
            (MIN_NUM_VARS..MAX_NUM_VARS).step_by(2).collect(),
            &bench_name,
            $method::<_, _, $scheme_type>,
            $rand_poly::<_>,
        );
    }};
}

fn ipa_pc_benches(c: &mut Criterion) {
    bench_pcs!(c, commit, IPA_JubJub, rand_poly_ipa_pc);
    bench_pcs!(c, open, IPA_JubJub, rand_poly_ipa_pc);
    bench_pcs!(c, verify, IPA_JubJub, rand_poly_ipa_pc);
}

criterion_group! {
    name = ipa_pc;
    config = Criterion::default();
    targets = ipa_pc_benches
}

criterion_main!(ipa_pc);
