<h1 align="center">`poly-commit`</h1>

<p align="center">
    <a href="https://travis-ci.org/scipr-lab/poly-commit"><img src="https://travis-ci.org/scipr-lab/poly-commit.svg?branch=master"></a>
    <a href="https://github.com/scipr-lab/poly-commit/blob/master/AUTHORS"><img src="https://img.shields.io/badge/authors-SCIPR%20Lab-orange.svg"></a>
    <a href="https://github.com/scipr-lab/poly-commit/blob/master/LICENSE-APACHE"><img src="https://img.shields.io/badge/license-APACHE-blue.svg"></a>
   <a href="https://github.com/scipr-lab/poly-commit/blob/master/LICENSE-MIT"><img src="https://img.shields.io/badge/license-MIT-blue.svg"></a>
</p>

`poly-commit` is a Rust library that implements *polynomial commitment schemes*, a cryptographic primitive that enables parties to "commit" to polynomials and then provide proofs of correct evaluation for these commitments.

This library was initially developed as part of the paper *"[<span style="font-variant:small-caps;">Marlin</span>: Preprocessing zkSNARKs with Universal and Updatable SRS][marlin]"*, and it is released under the MIT License and the Apache v2 License (see [License](#license)).

**WARNING:** This is an academic proof-of-concept prototype, and in particular has not received careful code review. This implementation is NOT ready for production use.

## Overview

This library provides traits for *polynomial commitment schemes*, as well as constructions of these based on the paper *"[Polynomial Commitments][kzg10]"*. Polynomial commitment (PC) schemes allow a user to first commit to polynomials, and to then provide proofs of correct evaluation for the committed polynomials. PC schemes satisfy the following properties:

- **Hiding** - commitments and proofs should reveal no information about the committed polynomial.
- **Succinctness** - commitment and proof size should scale sublinearly with the degree of the committed polynomial. Optionally, verification time should also scale sublinearly.
- **Extractability** - there should exist an extractor *E*, which when given a commitment and a valid evaluation proof, should be able to extract a polynomial that agrees with the evaluation.

[kzg10]: http://cacr.uwaterloo.ca/techreports/2010/cacr2010-10.pdf

## Build guide

The library compiles on the `stable` toolchain of the Rust compiler. To install the latest version of Rust, first install `rustup` by following the instructions [here](https://rustup.rs/), or via your platform's package manager. Once `rustup` is installed, install the Rust toolchain by invoking:
```bash
rustup install stable
```

After that, use `cargo`, the standard Rust build tool, to build the library:
```bash
git clone https://github.com/scipr-lab/poly-commit.git
cd poly-commit
cargo build --release
```

This library comes with some unit and integration tests. Run these tests with:
```bash
cargo test
```

Lastly, this library is instrumented with profiling infrastructure that prints detailed
profiles of execution time; to enables this, compile with `cargo build --features timer`.

## License

`poly-commit` is licensed under either of the following licenses, at your discretion.

 * Apache License Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE) or http://www.apache.org/licenses/LICENSE-2.0)
 * MIT license ([LICENSE-MIT](LICENSE-MIT) or http://opensource.org/licenses/MIT)

Unless you explicitly state otherwise, any contribution submitted for inclusion in `poly-commit` by you shall be dual licensed as above (as defined in the Apache v2 License), without any additional terms or conditions.

[marlin]: https://ia.cr/2019/xxx

## Reference paper

[_<span style="font-variant:small-caps;">Marlin</span>: Preprocessing zkSNARKs with Universal and Updatable SRS_][marlin]
Alessandro Chiesa, Yuncong Hu, Mary Maller, [Pratyush Mishra](https://www.github.com/pratyush), Noah Vesely, [Nicholas Ward](https://www.github.com/npwardberkeley)
*IACR ePrint Report 2019/XXX*

## Acknowledgements

This work was supported by:
an Engineering and Physical Sciences Research Council grant (EP/N028104/1),
a Google Faculty Award,
the RISELab at UC Berkeley,
and
donations from the Ethereum Foundation and the Interchain Foundation.
