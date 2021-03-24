use ark_bls12_381::Bls12_381;
use ark_dh_commitments::{
    afgho16::{AFGHOCommitmentG1, AFGHOCommitmentG2},
    identity::IdentityCommitment,
    pedersen::PedersenCommitment,
    DoublyHomomorphicCommitment,
};
use ark_ec::{group::Group, PairingEngine};
use ark_ff::{Field, UniformRand};
use ark_inner_products::{
    ExtensionFieldElement, InnerProduct, MultiexponentiationInnerProduct, PairingInnerProduct,
};
use ark_ip_proofs::tipa::{
    structured_scalar_message::{structured_scalar_power, TIPAWithSSM},
    TIPACompatibleSetup, TIPA,
};

use ark_std::rand::{rngs::StdRng, Rng, SeedableRng};
use blake2::Blake2b;
use digest::Digest;

use std::{ops::MulAssign, time::Instant};

fn bench_tipa<IP, LMC, RMC, IPC, P, D, R: Rng>(rng: &mut R, len: usize)
where
    D: Digest,
    P: PairingEngine,
    IP: InnerProduct<
        LeftMessage = LMC::Message,
        RightMessage = RMC::Message,
        Output = IPC::Message,
    >,
    LMC: DoublyHomomorphicCommitment<Scalar = P::Fr, Key = P::G2Projective> + TIPACompatibleSetup,
    RMC: DoublyHomomorphicCommitment<Scalar = LMC::Scalar, Key = P::G1Projective>
        + TIPACompatibleSetup,
    IPC: DoublyHomomorphicCommitment<Scalar = LMC::Scalar>,
    LMC::Message: MulAssign<P::Fr>,
    RMC::Message: MulAssign<P::Fr>,
    IPC::Message: MulAssign<P::Fr>,
    IPC::Key: MulAssign<P::Fr>,
    LMC::Output: MulAssign<P::Fr>,
    RMC::Output: MulAssign<P::Fr>,
    IPC::Output: MulAssign<P::Fr>,
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

    let (srs, ck_t) = TIPA::<IP, LMC, RMC, IPC, P, D>::setup(rng, len).unwrap();
    let (ck_l, ck_r) = srs.get_commitment_keys();
    let v_srs = srs.get_verifier_key();
    let com_l = LMC::commit(&ck_l, &l).unwrap();
    let com_r = RMC::commit(&ck_r, &r).unwrap();
    let t = vec![IP::inner_product(&l, &r).unwrap()];
    let com_t = IPC::commit(&vec![ck_t.clone()], &t).unwrap();
    let mut start = Instant::now();
    let proof =
        TIPA::<IP, LMC, RMC, IPC, P, D>::prove(&srs, (&l, &r), (&ck_l, &ck_r, &ck_t)).unwrap();
    let mut bench = start.elapsed().as_millis();
    println!("\t proving time: {} ms", bench);
    start = Instant::now();
    TIPA::<IP, LMC, RMC, IPC, P, D>::verify(&v_srs, &ck_t, (&com_l, &com_r, &com_t), &proof)
        .unwrap();
    bench = start.elapsed().as_millis();
    println!("\t verification time: {} ms", bench);
}

fn bench_tipa_srs_shift<IP, LMC, RMC, IPC, P, D, R: Rng>(rng: &mut R, len: usize)
where
    D: Digest,
    P: PairingEngine,
    IP: InnerProduct<
        LeftMessage = LMC::Message,
        RightMessage = RMC::Message,
        Output = IPC::Message,
    >,
    LMC: DoublyHomomorphicCommitment<Scalar = P::Fr, Key = P::G2Projective> + TIPACompatibleSetup,
    RMC: DoublyHomomorphicCommitment<Scalar = LMC::Scalar, Key = P::G1Projective>
        + TIPACompatibleSetup,
    IPC: DoublyHomomorphicCommitment<Scalar = LMC::Scalar>,
    LMC::Message: MulAssign<P::Fr> + Group<ScalarField = P::Fr>,
    RMC::Message: MulAssign<P::Fr> + Group,
    IPC::Message: MulAssign<P::Fr>,
    IPC::Key: MulAssign<P::Fr>,
    LMC::Output: MulAssign<P::Fr>,
    RMC::Output: MulAssign<P::Fr>,
    IPC::Output: MulAssign<P::Fr>,
    IPC::Output: MulAssign<LMC::Scalar>,
    IP::LeftMessage: UniformRand + Group,
    IP::RightMessage: UniformRand + Group,
{
    let mut l = Vec::new();
    let mut r = Vec::new();
    for _ in 0..len {
        l.push(<IP::LeftMessage>::rand(rng));
        r.push(<IP::RightMessage>::rand(rng));
    }

    let (srs, ck_t) = TIPA::<IP, LMC, RMC, IPC, P, D>::setup(rng, len).unwrap();
    let (ck_l, ck_r) = srs.get_commitment_keys();
    let v_srs = srs.get_verifier_key();
    let com_l = LMC::commit(&ck_l, &l).unwrap();
    let com_r = RMC::commit(&ck_r, &r).unwrap();
    let a_scalar = <P::Fr>::rand(rng);
    let r_vec = structured_scalar_power(len, &a_scalar);
    let l_a = l
        .iter()
        .zip(&r_vec)
        .map(|(a, r)| a.mul(r))
        .collect::<Vec<LMC::Message>>();
    let ck_l_a = ck_l
        .iter()
        .zip(&r_vec)
        .map(|(ck, r)| ck.mul(&r.inverse().unwrap()))
        .collect::<Vec<LMC::Key>>();

    let t = vec![IP::inner_product(&l_a, &r).unwrap()];
    let com_t = IPC::commit(&vec![ck_t.clone()], &t).unwrap();
    let mut start = Instant::now();
    let proof = TIPA::<IP, LMC, RMC, IPC, P, D>::prove_with_srs_shift(
        &srs,
        (&l_a, &r),
        (&ck_l_a, &ck_r, &ck_t),
        &a_scalar,
    )
    .unwrap();
    let mut bench = start.elapsed().as_millis();
    println!("\t proving time: {} ms", bench);
    start = Instant::now();
    TIPA::<IP, LMC, RMC, IPC, P, D>::verify_with_srs_shift(
        &v_srs,
        &ck_t,
        (&com_l, &com_r, &com_t),
        &proof,
        &a_scalar,
    )
    .unwrap();
    bench = start.elapsed().as_millis();
    println!("\t verification time: {} ms", bench);
}

fn bench_tipa_ssm<IP, LMC, IPC, P, D, R: Rng>(rng: &mut R, len: usize)
where
    D: Digest,
    P: PairingEngine,
    IP: InnerProduct<LeftMessage = LMC::Message, RightMessage = LMC::Scalar, Output = IPC::Message>,
    LMC: DoublyHomomorphicCommitment<Scalar = P::Fr, Key = P::G2Projective> + TIPACompatibleSetup,
    IPC: DoublyHomomorphicCommitment<Scalar = LMC::Scalar>,
    LMC::Message: MulAssign<P::Fr>,
    IPC::Message: MulAssign<P::Fr>,
    IPC::Key: MulAssign<P::Fr>,
    LMC::Output: MulAssign<P::Fr>,
    IPC::Output: MulAssign<P::Fr>,
    IPC::Output: MulAssign<LMC::Scalar>,
    IP::LeftMessage: UniformRand,
    IP::RightMessage: UniformRand,
{
    let mut l = Vec::new();
    for _ in 0..len {
        l.push(<IP::LeftMessage>::rand(rng));
    }
    let scalar = <P::Fr>::rand(rng);
    let r = structured_scalar_power(len, &scalar);

    let (srs, ck_t) = TIPAWithSSM::<IP, LMC, IPC, P, D>::setup(rng, len).unwrap();
    let (ck_l, _) = srs.get_commitment_keys();
    let v_srs = srs.get_verifier_key();
    let com_l = LMC::commit(&ck_l, &l).unwrap();
    let t = vec![IP::inner_product(&l, &r).unwrap()];
    let com_t = IPC::commit(&vec![ck_t.clone()], &t).unwrap();
    let mut start = Instant::now();
    let proof = TIPAWithSSM::<IP, LMC, IPC, P, D>::prove_with_structured_scalar_message(
        &srs,
        (&l, &r),
        (&ck_l, &ck_t),
    )
    .unwrap();
    let mut bench = start.elapsed().as_millis();
    println!("\t proving time: {} ms", bench);
    start = Instant::now();
    TIPAWithSSM::<IP, LMC, IPC, P, D>::verify_with_structured_scalar_message(
        &v_srs,
        &ck_t,
        (&com_l, &com_t),
        &scalar,
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

    println!("Benchmarking TIPA with vector length: {}", LEN);

    println!("1) Pairing inner product...");
    bench_tipa::<
        PairingInnerProduct<Bls12_381>,
        GC1,
        GC2,
        IdentityCommitment<ExtensionFieldElement<Bls12_381>, <Bls12_381 as PairingEngine>::Fr>,
        Bls12_381,
        Blake2b,
        StdRng,
    >(&mut rng, LEN);

    println!("2) Multiexponentiation G1 inner product...");
    bench_tipa::<
        MultiexponentiationInnerProduct<<Bls12_381 as PairingEngine>::G1Projective>,
        GC1,
        SC1,
        IdentityCommitment<
            <Bls12_381 as PairingEngine>::G1Projective,
            <Bls12_381 as PairingEngine>::Fr,
        >,
        Bls12_381,
        Blake2b,
        StdRng,
    >(&mut rng, LEN);

    println!("3) Pairing inner product with SRS shift...");
    bench_tipa_srs_shift::<
        PairingInnerProduct<Bls12_381>,
        GC1,
        GC2,
        IdentityCommitment<ExtensionFieldElement<Bls12_381>, <Bls12_381 as PairingEngine>::Fr>,
        Bls12_381,
        Blake2b,
        StdRng,
    >(&mut rng, LEN);

    println!("4) Multiexponentiation G1 inner product with structured scalar message...");
    bench_tipa_ssm::<
        MultiexponentiationInnerProduct<<Bls12_381 as PairingEngine>::G1Projective>,
        GC1,
        IdentityCommitment<
            <Bls12_381 as PairingEngine>::G1Projective,
            <Bls12_381 as PairingEngine>::Fr,
        >,
        Bls12_381,
        Blake2b,
        StdRng,
    >(&mut rng, LEN);
}
