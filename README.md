<h1 align="center">Polynomial Commitments</h1>

<p align="center">
   <a href="https://github.com/arkworks-rs/poly-commit/blob/master/LICENSE-APACHE"><img src="https://img.shields.io/badge/license-APACHE-blue.svg"></a>
   <a href="https://github.com/arkworks-rs/poly-commit/blob/master/LICENSE-MIT"><img src="https://img.shields.io/badge/license-MIT-blue.svg"></a>
</p>

`poly-commit` is a Rust library that implements *polynomial commitment schemes*. This library was initially developed as part of the [Marlin paper][marlin], and is released under the MIT License and the Apache v2 License (see [License](#license)).

**WARNING:** This is an academic prototype, and in particular has not received careful code review. This implementation is NOT ready for production use.

## Overview

A polynomial commitment scheme is a cryptographic primitive that enables a party to commit to a polynomial over a given finite field, and then, later on, to reveal desired evaluations of the polynomial along with cryptographic proofs attesting to their correctness.

This library provides various constructions of polynomial commitment schemes. These constructions support committing to multiple polynomials at a time with differing degree bounds, batching multiple evaluation proofs for the same evaluation point into a single one, and batch verification of proofs.

The key properties satisfied by the polynomial commitment schemes are **succinctness**, **extractability**, and **hiding**. See [the Marlin paper][marlin] for definitions of these properties.


[kzg10]: http://cacr.uwaterloo.ca/techreports/2010/cacr2010-10.pdf

## Build guide

The library compiles on the `stable` toolchain of the Rust compiler. To install the latest version of Rust, first install `rustup` by following the instructions [here](https://rustup.rs/), or via your platform's package manager. Once `rustup` is installed, install the Rust toolchain by invoking:
```bash
rustup install stable
```

After that, use `cargo` (the standard Rust build tool) to build the library:
```bash
git clone https://github.com/scipr-lab/poly-commit.git
cd poly-commit
cargo build --release
```

This library comes with some unit and integration tests. Run these tests with:
```bash
cargo test
```

Lastly, this library is instrumented with profiling infrastructure that prints detailed traces of execution time. To enable this, compile with `cargo build --features print-trace`.

## Usage

### [`PolynomialCommitment`](https://github.com/arkworks-rs/poly-commit/blob/master/src/lib.rs#L145)

```rust
// In this example, we will commit to a single polynomial, open it at two points, and finally verify the proof.
// We will use the KZG10 polynomial commitment scheme from Marlin: src/marlin/marlin_pc/mod.rs

use ark_poly_commit::{marlin_pc::MarlinKZG10, LabeledPolynomial, PolynomialCommitment, QuerySet, Evaluations, challenge::ChallengeGenerator};
use ark_bls12_377::Bls12_377;
use ark_crypto_primitives::sponge::poseidon::PoseidonSponge;
use ark_ec::pairing::Pairing;
use ark_ff::UniformRand;
use ark_std::test_rng;
use ark_poly::{univariate::DensePolynomial, DenseUVPolynomial};
use rand_chacha::ChaCha20Rng;

type UniPoly_377 = DensePolynomial<<Bls12_377 as Pairing>::ScalarField>;
type Sponge_Bls12_377 = PoseidonSponge<<Bls12_377 as Pairing>::ScalarField>;
type PCS = MarlinKZG10<Bls12_377, UniPoly_377, Sponge_Bls12_377>;

let rng = &mut test_rng();

let max_degree = 12; // max degree supported by the scheme with the given public parameters generated by the setup here.
let pp = PCS::setup(max_degree, None, rng).unwrap();

let degree = 10; //degree of our polynomial
let secret_poly = UniPoly_377::rand(degree, rng);

let point_1 = <Bls12_377 as Pairing>::ScalarField::rand(rng);
let point_2 = <Bls12_377 as Pairing>::ScalarField::rand(rng);

let labeled_poly = LabeledPolynomial::new(
   String::from("secret_poly"),
   secret_poly,
   Some(degree),
   Some(2), // we will open a univariate poly at two points
);

// let sponge_params = <Bls12_377 as Pairing>::ScalarField::get_default_poseidon_parameters(2, false).unwrap();
let (ck, vk) = PCS::trim(&pp, degree, 2, None).unwrap(); // since the setup produced pp with a max degree of 12, and our poly is of degree 10, we can trim the SRS (for efficiency?).

let (comms, rands) = PCS::commit(&ck, &[labeled_poly], Some(rng)).unwrap(); // the prover (aka committer) commits to the polynomial.

// let sponge = poseidon_sponge_for_test::<<Bls12_377 as Pairing>::ScalarField>();
// let challenge_generator = ChallengeGenerator::new_univariate(&mut sponge());
// let's get the opening proof at a single point first
// let proof = PCS::open(&ck, &[labeled_poly], &comms, &point_1, &challenge_generator, None, None).unwrap();

```

## License

This library is licensed under either of the following licenses, at your discretion.

 * [Apache License Version 2.0](LICENSE-APACHE)
 * [MIT License](LICENSE-MIT)

Unless you explicitly state otherwise, any contribution that you submit to this library shall be dual licensed as above (as defined in the Apache v2 License), without any additional terms or conditions.

[marlin]: https://ia.cr/2019/1047
[sonic]: https://ia.cr/2019/099
[aurora-light]: https://ia.cr/2019/601
[pcd-acc]: https://ia.cr/2020/499
[pst]: https://ia.cr/2011/587

## Reference papers

[Polynomial Commitments][kzg10]     
Aniket Kate, Gregory M. Zaverucha, Ian Goldberg     
ASIACRYPT 2010

[Sonic: Zero-Knowledge SNARKs from Linear-Size Universal and Updateable Structured Reference Strings][sonic]     
Mary Maller, Sean Bowe, Markulf Kohlweiss, Sarah Meiklejohn     
CCS 2019

[AuroraLight: Improved Prover Efficiency and SRS Size in a Sonic-Like System][aurora-light]     
Ariel Gabizon     
ePrint, 2019

[Marlin: Preprocessing zkSNARKs with Universal and Updatable SRS][marlin]     
Alessandro Chiesa, Yuncong Hu, Mary Maller, [Pratyush Mishra](https://www.github.com/pratyush), Noah Vesely, [Nicholas Ward](https://www.github.com/npwardberkeley)     
EUROCRYPT 2020

[Proof-Carrying Data from Accumulation Schemes][pcd-acc]     
Benedikt Bünz, Alessandro Chiesa, [Pratyush Mishra](https://www.github.com/pratyush), Nicholas Spooner     
TCC 2020

[Signatures of Correct Computation][pst]    
Charalampos Papamanthou, Elaine Shi, Roberto Tamassia   
TCC 2013

## Acknowledgements

This work was supported by: an Engineering and Physical Sciences Research Council grant; a Google Faculty Award; the RISELab at UC Berkeley; and donations from the Ethereum Foundation and the Interchain Foundation.
