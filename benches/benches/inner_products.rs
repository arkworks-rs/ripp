use ark_bls12_381::Bls12_381;
use ark_ec::PairingEngine;
use ark_ff::UniformRand;
use ark_inner_products::{InnerProduct, MultiexponentiationInnerProduct, PairingInnerProduct};

use ark_std::rand::{rngs::StdRng, Rng, SeedableRng};

use std::time::Instant;

fn bench_inner_product<IP: InnerProduct, R: Rng>(rng: &mut R, len: usize)
where
    IP::LeftMessage: UniformRand,
    IP::RightMessage: UniformRand,
{
    let mut l = Vec::new();
    let mut r = Vec::new();
    for _ in 0..len {
        l.push(<IP::LeftMessage>::rand(rng));
        r.push(<IP::RightMessage>::rand(rng));
    }
    let start = Instant::now();
    IP::inner_product(&l, &r).unwrap();
    let bench = start.elapsed().as_millis();
    println!("\t time: {} ms", bench);
}

fn main() {
    const LEN: usize = 16;
    let mut rng = StdRng::seed_from_u64(0u64);
    println!("Benchmarking inner products with vector length: {}", LEN);

    println!("1) Pairing inner product...");
    bench_inner_product::<PairingInnerProduct<Bls12_381>, StdRng>(&mut rng, LEN);

    println!("2) Multiexponentiation G1 inner product...");
    bench_inner_product::<
        MultiexponentiationInnerProduct<<Bls12_381 as PairingEngine>::G1Projective>,
        StdRng,
    >(&mut rng, LEN);

    println!("3) Multiexponentiation G2 inner product...");
    bench_inner_product::<
        MultiexponentiationInnerProduct<<Bls12_381 as PairingEngine>::G2Projective>,
        StdRng,
    >(&mut rng, LEN);
}
