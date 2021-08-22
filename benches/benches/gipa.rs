use ark_bls12_381::Bls12_381;
use ark_dh_commitments::{
    afgho16::{AFGHOCommitmentG1, AFGHOCommitmentG2},
    identity::IdentityCommitment,
    pedersen::PedersenCommitment,
    DoublyHomomorphicCommitment,
};
use ark_ec::PairingEngine;
use ark_ff::UniformRand;
use ark_inner_products::{
    ExtensionFieldElement, InnerProduct, MultiexponentiationInnerProduct, PairingInnerProduct,
};
use ark_ip_proofs::gipa::GIPA;

use ark_std::rand::{rngs::StdRng, Rng, SeedableRng};

use merlin::Transcript;
use std::{ops::MulAssign, time::Instant};

fn bench_gipa<IP, LMC, RMC, IPC, R: Rng>(rng: &mut R, len: usize)
where
    IP: InnerProduct<
        LeftMessage = LMC::Message,
        RightMessage = RMC::Message,
        Output = IPC::Message,
    >,
    LMC: DoublyHomomorphicCommitment,
    RMC: DoublyHomomorphicCommitment<Scalar = LMC::Scalar>,
    IPC: DoublyHomomorphicCommitment<Scalar = LMC::Scalar>,
    RMC::Message: MulAssign<LMC::Scalar>,
    IPC::Message: MulAssign<LMC::Scalar>,
    RMC::Key: MulAssign<LMC::Scalar>,
    IPC::Key: MulAssign<LMC::Scalar>,
    RMC::Output: MulAssign<LMC::Scalar>,
    IPC::Output: MulAssign<LMC::Scalar>,
    IP::LeftMessage: UniformRand,
    IP::RightMessage: UniformRand,
{
    let mut l = Vec::new();
    let mut r = Vec::new();
    for _ in 0..len {
        l.push(<IP::LeftMessage>::rand(rng));
        r.push(<IP::RightMessage>::rand(rng));
    }

    let (ck_l, ck_r, ck_t) = GIPA::<IP, LMC, RMC, IPC>::setup(rng, len).unwrap();
    let com_l = LMC::commit(&ck_l, &l).unwrap();
    let com_r = RMC::commit(&ck_r, &r).unwrap();
    let t = vec![IP::inner_product(&l, &r).unwrap()];
    let com_t = IPC::commit(&vec![ck_t.clone()], &t).unwrap();
    let mut proof_transcript = Transcript::new(b"GIPA-bench");
    let mut start = Instant::now();
    let proof = GIPA::<IP, LMC, RMC, IPC>::prove(
        &mut proof_transcript,
        (&l, &r, &t[0]),
        (&ck_l, &ck_r, &ck_t),
        (&com_l, &com_r, &com_t),
    )
    .unwrap();
    let mut bench = start.elapsed().as_millis();
    println!("\t proving time: {} ms", bench);

    let mut verif_transcript = Transcript::new(b"GIPA-bench");
    start = Instant::now();
    GIPA::<IP, LMC, RMC, IPC>::verify(
        &mut verif_transcript,
        (&ck_l, &ck_r, &ck_t),
        (&com_l, &com_r, &com_t),
        &proof,
    )
    .unwrap();
    bench = start.elapsed().as_millis();
    println!("\t verification time: {} ms", bench);
}

fn main() {
    const LEN: usize = 16;
    type GC1 = AFGHOCommitmentG1<Bls12_381>;
    type GC2 = AFGHOCommitmentG2<Bls12_381>;
    type SC1 = PedersenCommitment<<Bls12_381 as PairingEngine>::G1Projective>;
    let mut rng = StdRng::seed_from_u64(0u64);

    println!("Benchmarking GIPA with vector length: {}", LEN);

    println!("1) Pairing inner product...");
    bench_gipa::<
        PairingInnerProduct<Bls12_381>,
        GC1,
        GC2,
        IdentityCommitment<ExtensionFieldElement<Bls12_381>, <Bls12_381 as PairingEngine>::Fr>,
        StdRng,
    >(&mut rng, LEN);

    println!("2) Multiexponentiation G1 inner product...");
    bench_gipa::<
        MultiexponentiationInnerProduct<<Bls12_381 as PairingEngine>::G1Projective>,
        GC1,
        SC1,
        IdentityCommitment<
            <Bls12_381 as PairingEngine>::G1Projective,
            <Bls12_381 as PairingEngine>::Fr,
        >,
        StdRng,
    >(&mut rng, LEN);
}
