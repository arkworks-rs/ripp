<h1 align="center">RIPP (Rust Inner Pairing Products)</h1>

<p align="center">
    <a href="https://travis-ci.org/scipr-lab/ripp"><img src="https://travis-ci.org/scipr-lab/ripp.svg?branch=master"></a>
    <a href="https://github.com/scipr-lab/ripp/blob/master/AUTHORS"><img src="https://img.shields.io/badge/authors-SCIPR%20Lab-orange.svg"></a>
    <a href="https://github.com/scipr-lab/ripp/blob/master/LICENSE-APACHE"><img src="https://img.shields.io/badge/license-APACHE-blue.svg"></a>
    <a href="https://github.com/scipr-lab/ripp/blob/master/LICENSE-MIT"><img src="https://img.shields.io/badge/license-MIT-blue.svg"></a>
    <a href="https://deps.rs/repo/github/scipr-lab/ripp"><img src="https://deps.rs/repo/github/scipr-lab/ripp/status.svg"></a>
</p>


___RIPP___ is a Rust library for proofs about inner pairing products, and applications built atop these. These protocols and applications are described in our paper *"[Proofs for Inner Pairing Products and Applications][ripp]"*

The library currently contains an implementation of our proof system for verifiably outsourcing pairing products. In the future, we intend to implement the other protocols described in our [paper][ripp], along with the polynomial commitment schemes and our protocol for aggregating Groth16 proofs based upon these protocols.

This library is released under the MIT License and the Apache v2 License (see [License](#license)).

**WARNING:** This is an academic proof-of-concept prototype, and in particular has not received careful code review. This implementation is NOT ready for production use.

## Build guide

The library compiles on the `stable` toolchain of the Rust compiler. To install the latest version of Rust, first install `rustup` by following the instructions [here](https://rustup.rs/), or via your platform's package manager. Once `rustup` is installed, install the Rust toolchain by invoking:
```bash
rustup install stable
```

After that, use `cargo`, the standard Rust build tool, to build the library:
```bash
git clone https://github.com/scipr-lab/ripp.git
cd ripp
cargo build --release
```

This library comes with unit tests for each of the provided crates. Run the tests with:
```bash
cargo test
``` 

Lastly, the library comes with benchmarks.
```bash
cargo bench
cargo run --release --example groth16_aggregation
cargo run --release --example scaling-ipp
```

## License

RIPP is licensed under either of the following licenses, at your discretion.

 * Apache License Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE) or http://www.apache.org/licenses/LICENSE-2.0)
 * MIT license ([LICENSE-MIT](LICENSE-MIT) or http://opensource.org/licenses/MIT)

Unless you explicitly state otherwise, any contribution submitted for inclusion in RIPP by you shall be dual licensed as above (as defined in the Apache v2 License), without any additional terms or conditions.

[ripp]: https://eprint.iacr.org/2019/1177

## Reference paper

[_Proofs for Inner Pairing Products and Applications_][ripp]    
[Benedikt Bünz](https://www.github.com/bbuenz), Mary Maller, [Pratyush Mishra](https://www.github.com/pratyush), [Psi Vesely](https://www.github.com/psivesely)    
*IACR ePrint Report 2019/1177*
