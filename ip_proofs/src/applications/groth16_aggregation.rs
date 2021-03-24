use ark_ec::{group::Group, AffineCurve, PairingEngine};
use ark_ff::{to_bytes, Field, One};
use ark_groth16::{Proof, VerifyingKey};

use std::ops::AddAssign;

use ark_std::rand::Rng;
use digest::Digest;

use crate::{
    tipa::{
        structured_scalar_message::{structured_scalar_power, TIPAWithSSM, TIPAWithSSMProof},
        TIPAProof, VerifierSRS, SRS, TIPA,
    },
    Error,
};
use ark_dh_commitments::{
    afgho16::{AFGHOCommitmentG1, AFGHOCommitmentG2},
    identity::{HomomorphicPlaceholderValue, IdentityCommitment, IdentityOutput},
};
use ark_inner_products::{
    ExtensionFieldElement, InnerProduct, MultiexponentiationInnerProduct, PairingInnerProduct,
    ScalarInnerProduct,
};

type PairingInnerProductAB<P, D> = TIPA<
    PairingInnerProduct<P>,
    AFGHOCommitmentG1<P>,
    AFGHOCommitmentG2<P>,
    IdentityCommitment<ExtensionFieldElement<P>, <P as PairingEngine>::Fr>,
    P,
    D,
>;

type PairingInnerProductABProof<P, D> = TIPAProof<
    PairingInnerProduct<P>,
    AFGHOCommitmentG1<P>,
    AFGHOCommitmentG2<P>,
    IdentityCommitment<ExtensionFieldElement<P>, <P as PairingEngine>::Fr>,
    P,
    D,
>;

type MultiExpInnerProductC<P, D> = TIPAWithSSM<
    MultiexponentiationInnerProduct<<P as PairingEngine>::G1Projective>,
    AFGHOCommitmentG1<P>,
    IdentityCommitment<<P as PairingEngine>::G1Projective, <P as PairingEngine>::Fr>,
    P,
    D,
>;

type MultiExpInnerProductCProof<P, D> = TIPAWithSSMProof<
    MultiexponentiationInnerProduct<<P as PairingEngine>::G1Projective>,
    AFGHOCommitmentG1<P>,
    IdentityCommitment<<P as PairingEngine>::G1Projective, <P as PairingEngine>::Fr>,
    P,
    D,
>;

pub struct AggregateProof<P: PairingEngine, D: Digest> {
    com_a: ExtensionFieldElement<P>,
    com_b: ExtensionFieldElement<P>,
    com_c: ExtensionFieldElement<P>,
    ip_ab: ExtensionFieldElement<P>,
    agg_c: P::G1Projective,
    tipa_proof_ab: PairingInnerProductABProof<P, D>,
    tipa_proof_c: MultiExpInnerProductCProof<P, D>,
}

pub fn setup_inner_product<P, D, R: Rng>(rng: &mut R, size: usize) -> Result<SRS<P>, Error>
where
    P: PairingEngine,
    D: Digest,
{
    let (srs, _) = PairingInnerProductAB::<P, D>::setup(rng, size)?;
    Ok(srs)
}

pub fn aggregate_proofs<P, D>(
    ip_srs: &SRS<P>,
    proofs: &[Proof<P>],
) -> Result<AggregateProof<P, D>, Error>
where
    P: PairingEngine,
    D: Digest,
{
    let a = proofs
        .iter()
        .map(|proof| proof.a.into_projective())
        .collect::<Vec<P::G1Projective>>();
    let b = proofs
        .iter()
        .map(|proof| proof.b.into_projective())
        .collect::<Vec<P::G2Projective>>();
    let c = proofs
        .iter()
        .map(|proof| proof.c.into_projective())
        .collect::<Vec<P::G1Projective>>();

    let (ck_1, ck_2) = ip_srs.get_commitment_keys();

    let com_a = PairingInnerProduct::<P>::inner_product(&a, &ck_1)?;
    let com_b = PairingInnerProduct::<P>::inner_product(&ck_2, &b)?;
    let com_c = PairingInnerProduct::<P>::inner_product(&c, &ck_1)?;

    // Random linear combination of proofs
    let mut counter_nonce: usize = 0;
    let r = loop {
        let mut hash_input = Vec::new();
        hash_input.extend_from_slice(&counter_nonce.to_be_bytes()[..]);
        //TODO: Should use CanonicalSerialize instead of ToBytes
        hash_input.extend_from_slice(&to_bytes![com_a, com_b, com_c]?);
        if let Some(r) = <P::Fr>::from_random_bytes(&D::digest(&hash_input)) {
            break r;
        };
        counter_nonce += 1;
    };

    let r_vec = structured_scalar_power(proofs.len(), &r);
    let a_r = a
        .iter()
        .zip(&r_vec)
        .map(|(a, r)| a.mul(r))
        .collect::<Vec<P::G1Projective>>();
    let ip_ab = PairingInnerProduct::<P>::inner_product(&a_r, &b)?;
    let agg_c = MultiexponentiationInnerProduct::<P::G1Projective>::inner_product(&c, &r_vec)?;

    let ck_1_r = ck_1
        .iter()
        .zip(&r_vec)
        .map(|(ck, r)| ck.mul(&r.inverse().unwrap()))
        .collect::<Vec<P::G2Projective>>();

    assert_eq!(
        com_a,
        PairingInnerProduct::<P>::inner_product(&a_r, &ck_1_r)?
    );

    let tipa_proof_ab = PairingInnerProductAB::<P, D>::prove_with_srs_shift(
        &ip_srs,
        (&a_r, &b),
        (&ck_1_r, &ck_2, &HomomorphicPlaceholderValue),
        &r,
    )?;

    let tipa_proof_c = MultiExpInnerProductC::<P, D>::prove_with_structured_scalar_message(
        &ip_srs,
        (&c, &r_vec),
        (&ck_1, &HomomorphicPlaceholderValue),
    )?;

    Ok(AggregateProof {
        com_a,
        com_b,
        com_c,
        ip_ab,
        agg_c,
        tipa_proof_ab,
        tipa_proof_c,
    })
}

pub fn verify_aggregate_proof<P, D>(
    ip_verifier_srs: &VerifierSRS<P>,
    vk: &VerifyingKey<P>,
    public_inputs: &Vec<Vec<P::Fr>>, //TODO: Should use ToConstraintField instead
    proof: &AggregateProof<P, D>,
) -> Result<bool, Error>
where
    P: PairingEngine,
    D: Digest,
{
    // Random linear combination of proofs
    let mut counter_nonce: usize = 0;
    let r = loop {
        let mut hash_input = Vec::new();
        hash_input.extend_from_slice(&counter_nonce.to_be_bytes()[..]);
        //TODO: Should use CanonicalSerialize instead of ToBytes
        hash_input.extend_from_slice(&to_bytes![proof.com_a, proof.com_b, proof.com_c]?);
        if let Some(r) = <P::Fr>::from_random_bytes(&D::digest(&hash_input)) {
            break r;
        };
        counter_nonce += 1;
    };

    // Check TIPA proofs
    let tipa_proof_ab_valid = PairingInnerProductAB::<P, D>::verify_with_srs_shift(
        ip_verifier_srs,
        &HomomorphicPlaceholderValue,
        (
            &proof.com_a,
            &proof.com_b,
            &IdentityOutput(vec![proof.ip_ab.clone()]),
        ),
        &proof.tipa_proof_ab,
        &r,
    )?;
    let tipa_proof_c_valid = MultiExpInnerProductC::<P, D>::verify_with_structured_scalar_message(
        ip_verifier_srs,
        &HomomorphicPlaceholderValue,
        (&proof.com_c, &IdentityOutput(vec![proof.agg_c.clone()])),
        &r,
        &proof.tipa_proof_c,
    )?;

    // Check aggregate pairing product equation

    let r_sum =
        (r.pow(&[public_inputs.len() as u64]) - &<P::Fr>::one()) / &(r.clone() - &<P::Fr>::one());
    let p1 = P::pairing(vk.alpha_g1.into_projective().mul(&r_sum), vk.beta_g2);

    assert_eq!(vk.gamma_abc_g1.len(), public_inputs[0].len() + 1);
    let r_vec = structured_scalar_power(public_inputs.len(), &r);
    let mut g_ic = vk.gamma_abc_g1[0].into_projective().mul(&r_sum);
    for (i, b) in vk.gamma_abc_g1.iter().skip(1).enumerate() {
        g_ic.add_assign(
            &b.into_projective().mul(&ScalarInnerProduct::inner_product(
                &public_inputs
                    .iter()
                    .map(|inputs| inputs[i].clone())
                    .collect::<Vec<P::Fr>>(),
                &r_vec,
            )?),
        );
    }
    let p2 = P::pairing(g_ic, vk.gamma_g2);
    let p3 = P::pairing(proof.agg_c, vk.delta_g2);

    let ppe_valid = proof.ip_ab.0 == (p1 * &p2) * &p3;

    Ok(tipa_proof_ab_valid && tipa_proof_c_valid && ppe_valid)
}
