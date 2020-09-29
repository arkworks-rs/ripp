use algebra::{
    curves::{PairingEngine, ProjectiveCurve},
    groups::Group,
    fields::{Field, PrimeField},
    to_bytes, UniformRand,
};
use digest::Digest;
use num_traits::identities::One;
use std::{
    ops::MulAssign,
    marker::PhantomData,
};
use rand::Rng;

use crate::{
    gipa::{GIPA, GIPAProof},
    tipa::{
        TIPACompatibleSetup, VerifierSRS, SRS,
        prove_commitment_key_kzg_opening,
        verify_commitment_key_g2_kzg_opening,
        structured_generators_scalar_power,
    },
    Error,
};
use dh_commitments::{
    DoublyHomomorphicCommitment,
    identity::{HomomorphicPlaceholderValue},
};
use inner_products::InnerProduct;

//TODO: Properly generalize the non-committed message approach of SIPP and MIPP to GIPA
//TODO: Structured message is a special case of the non-committed message and does not rely on TIPA
//TODO: Can support structured group element messages as well as structured scalar messages

// Use placeholder commitment to commit to vector in clear during GIPA execution
#[derive(Clone)]
struct SSMPlaceholderCommitment<F: PrimeField> {
    _field: PhantomData<F>,
}

impl<F: PrimeField> DoublyHomomorphicCommitment for SSMPlaceholderCommitment<F> {
    type Scalar = F;
    type Message = F;
    type Key = HomomorphicPlaceholderValue;
    type Output = F;

    fn setup<R: Rng>(_rng: &mut R, size: usize) -> Result<Vec<Self::Key>, Error> {
        Ok(vec![HomomorphicPlaceholderValue {}; size])
    }

    fn commit(_k: &[Self::Key], m: &[Self::Message]) -> Result<Self::Output, Error> {
        Ok(m[0].clone())
    }
}

pub struct TIPAWithSSM<IP, LMC, IPC, P, D>
{
    _inner_product: PhantomData<IP>,
    _left_commitment: PhantomData<LMC>,
    _inner_product_commitment: PhantomData<IPC>,
    _pair: PhantomData<P>,
    _digest: PhantomData<D>,
}

pub struct TIPAWithSSMProof<IP, LMC, IPC, P, D>
where
    D: Digest,
    P: PairingEngine,
    IP: InnerProduct<
        LeftMessage = LMC::Message,
        RightMessage = LMC::Scalar,
        Output = IPC::Message,
    >,
    LMC: DoublyHomomorphicCommitment<Scalar = P::Fr, Key = P::G2Projective> + TIPACompatibleSetup,
    IPC: DoublyHomomorphicCommitment<Scalar = LMC::Scalar>,
    IPC::Message: MulAssign<LMC::Scalar>,
    IPC::Key: MulAssign<LMC::Scalar>,
    IPC::Output: MulAssign<LMC::Scalar>,
    LMC::Message: MulAssign<P::Fr>,
    LMC::Output: MulAssign<P::Fr>,
{
    gipa_proof: GIPAProof<IP, LMC, SSMPlaceholderCommitment<LMC::Scalar>, IPC, D>,
    final_ck: LMC::Key,
    final_ck_proof: P::G2Projective,
    _pairing: PhantomData<P>,
}

impl<IP, LMC, IPC, P, D> Clone for TIPAWithSSMProof<IP, LMC, IPC, P, D>
    where
        D: Digest,
        P: PairingEngine,
        IP: InnerProduct<
            LeftMessage = LMC::Message,
            RightMessage = LMC::Scalar,
            Output = IPC::Message,
        >,
        LMC: DoublyHomomorphicCommitment<Scalar = P::Fr, Key = P::G2Projective> + TIPACompatibleSetup,
        IPC: DoublyHomomorphicCommitment<Scalar = LMC::Scalar>,
        IPC::Message: MulAssign<LMC::Scalar>,
        IPC::Key: MulAssign<LMC::Scalar>,
        IPC::Output: MulAssign<LMC::Scalar>,
        LMC::Message: MulAssign<P::Fr>,
        LMC::Output: MulAssign<P::Fr>,
{
    fn clone(&self) -> Self {
        Self {
            gipa_proof: self.gipa_proof.clone(),
            final_ck: self.final_ck.clone(),
            final_ck_proof: self.final_ck_proof.clone(),
            _pairing: PhantomData,
        }
    }
}

impl<IP, LMC, IPC, P, D> TIPAWithSSM<IP, LMC, IPC, P, D>
where
    D: Digest,
    P: PairingEngine,
    IP: InnerProduct<
        LeftMessage = LMC::Message,
        RightMessage = LMC::Scalar,
        Output = IPC::Message,
    >,
    LMC: DoublyHomomorphicCommitment<Scalar = P::Fr, Key = P::G2Projective> + TIPACompatibleSetup,
    IPC: DoublyHomomorphicCommitment<Scalar = LMC::Scalar>,
    IPC::Message: MulAssign<P::Fr>,
    IPC::Key: MulAssign<P::Fr>,
    IPC::Output: MulAssign<P::Fr>,
    LMC::Message: MulAssign<P::Fr>,
    LMC::Output: MulAssign<P::Fr>,
{
    //TODO: Don't need full TIPA SRS since only using one set of powers
    pub fn setup<R: Rng>(rng: &mut R, size: usize) -> Result<(SRS<P>, IPC::Key), Error> {
        let alpha = <P::Fr>::rand(rng);
        let beta = <P::Fr>::rand(rng);
        let g = <P::G1Projective>::prime_subgroup_generator();
        let h = <P::G2Projective>::prime_subgroup_generator();
        Ok((
            SRS {
                g_alpha_powers: structured_generators_scalar_power(2 * size - 1, &g, &alpha),
                h_beta_powers: structured_generators_scalar_power(2 * size - 1, &h, &beta),
                g_beta: <P::G1Projective as Group>::mul(&g, &beta),
                h_alpha: <P::G2Projective as Group>::mul(&h, &alpha),
            },
            IPC::setup(rng, 1)?.pop().unwrap(),
        ))
    }

    pub fn prove_with_structured_scalar_message(
        srs: &SRS<P>,
        values: (&[IP::LeftMessage], &[IP::RightMessage]),
        ck: (&[LMC::Key], &IPC::Key),
    ) -> Result<TIPAWithSSMProof<IP, LMC, IPC, P, D>, Error> {
        // Run GIPA
        let (proof, aux) = <GIPA<IP, LMC, SSMPlaceholderCommitment<P::Fr>, IPC, D>>::prove_with_aux(
            values,
            (ck.0, &vec![HomomorphicPlaceholderValue {}; values.1.len()], &vec![ck.1.clone()]),
        )?;

        // Prove final commitment key is wellformed
        let (ck_a_final, _) = aux.ck_base;
        let transcript = aux.r_transcript;
        let transcript_inverse = transcript.iter().map(|x| x.inverse().unwrap()).collect();

        // KZG challenge point
        let mut counter_nonce: usize = 0;
        let c = loop {
            let mut hash_input = Vec::new();
            hash_input.extend_from_slice(&counter_nonce.to_be_bytes()[..]);
            //TODO: Should use CanonicalSerialize instead of ToBytes
            hash_input.extend_from_slice(&to_bytes![
                transcript.first().unwrap(),
                ck_a_final
            ]?);
            if let Some(c) = LMC::Scalar::from_random_bytes(&D::digest(&hash_input)) {
                break c;
            };
            counter_nonce += 1;
        };

        // Complete KZG proof
        let ck_a_kzg_opening = prove_commitment_key_kzg_opening(
            &srs.h_beta_powers,
            &transcript_inverse,
            &<P::Fr>::one(),
            &c,
        )?;

        Ok(TIPAWithSSMProof {
            gipa_proof: proof,
            final_ck: ck_a_final,
            final_ck_proof: ck_a_kzg_opening,
            _pairing: PhantomData,
        })
    }

    pub fn verify_with_structured_scalar_message(
        v_srs: &VerifierSRS<P>,
        ck_t: &IPC::Key,
        com: (&LMC::Output, &IPC::Output),
        scalar_b: &P::Fr,
        proof: &TIPAWithSSMProof<IP, LMC, IPC, P, D>,
    ) -> Result<bool, Error> {
        let (base_com, transcript) =
            GIPA::verify_recursive_challenge_transcript((com.0, scalar_b, com.1), &proof.gipa_proof)?;
        let transcript_inverse = transcript.iter().map(|x| x.inverse().unwrap()).collect();

        let ck_a_final= &proof.final_ck;
        let ck_a_proof= &proof.final_ck_proof;

        // KZG challenge point
        let mut counter_nonce: usize = 0;
        let c = loop {
            let mut hash_input = Vec::new();
            hash_input.extend_from_slice(&counter_nonce.to_be_bytes()[..]);
            //TODO: Should use CanonicalSerialize instead of ToBytes
            hash_input.extend_from_slice(&to_bytes![
                transcript.first().unwrap(),
                ck_a_final
            ]?);
            if let Some(c) = LMC::Scalar::from_random_bytes(&D::digest(&hash_input)) {
                break c;
            };
            counter_nonce += 1;
        };

        // Check commitment key
        let ck_a_valid = verify_commitment_key_g2_kzg_opening(
            v_srs,
            &ck_a_final,
            &ck_a_proof,
            &transcript_inverse,
            &<P::Fr>::one(),
            &c,
        )?;

        // Compute final scalar
        let mut power_2_b = scalar_b.clone();
        let mut product_form = Vec::new();
        for x in transcript.iter() {
            product_form.push(<P::Fr>::one() + &(x.inverse().unwrap() * &power_2_b));
            power_2_b *= &power_2_b.clone();
        }
        let b_base = product_form.iter().product::<P::Fr>();

        // Verify base inner product commitment
        let (com_a, _, com_t) = base_com;
        let a_base = vec![proof.gipa_proof.r_base.0.clone()];
        let t_base = vec![IP::inner_product(&a_base, &vec![b_base])?];
        let base_valid = LMC::verify(&vec![ck_a_final.clone()], &a_base, &com_a)?
            && IPC::verify(&vec![ck_t.clone()], &t_base, &com_t)?;

        Ok(ck_a_valid && base_valid)
    }
}

pub fn structured_scalar_power<F: Field>(num: usize, s: &F) -> Vec<F> {
    let mut powers = vec![F::one()];
    for i in 1..num {
        powers.push(powers[i - 1] * s);
    }
    powers
}

#[cfg(test)]
mod tests {
    use super::*;
    use algebra::{bls12_381::Bls12_381, curves::PairingEngine, UniformRand};
    use blake2::Blake2b;
    use rand::{rngs::StdRng, SeedableRng};

    use dh_commitments::{
        afgho16::AFGHOCommitmentG1, identity::IdentityCommitment, pedersen::PedersenCommitment,
        random_generators,
    };
    use inner_products::{InnerProduct, MultiexponentiationInnerProduct, ScalarInnerProduct};

    type GC1 = AFGHOCommitmentG1<Bls12_381>;
    type SC2 = PedersenCommitment<<Bls12_381 as PairingEngine>::G2Projective>;

    const TEST_SIZE: usize = 8;

    #[test]
    fn tipassm_multiexponentiation_inner_product_test() {
        type IP = MultiexponentiationInnerProduct<<Bls12_381 as PairingEngine>::G1Projective>;
        type IPC = IdentityCommitment<
            <Bls12_381 as PairingEngine>::G1Projective,
            <Bls12_381 as PairingEngine>::Fr,
        >;
        type MultiExpTIPA = TIPAWithSSM<IP, GC1, IPC, Bls12_381, Blake2b>;

        let mut rng = StdRng::seed_from_u64(0u64);
        let (srs, ck_t) = MultiExpTIPA::setup(&mut rng, TEST_SIZE).unwrap();
        let (ck_a, _) = srs.get_commitment_keys();
        let v_srs = srs.get_verifier_key();
        let m_a = random_generators(&mut rng, TEST_SIZE);
        let b = <<Bls12_381 as PairingEngine>::Fr>::rand(&mut rng);
        let m_b = structured_scalar_power(TEST_SIZE, &b);
        let com_a = GC1::commit(&ck_a, &m_a).unwrap();
        let t = vec![IP::inner_product(&m_a, &m_b).unwrap()];
        let com_t = IPC::commit(&vec![ck_t.clone()], &t).unwrap();

        let proof = MultiExpTIPA::prove_with_structured_scalar_message(
            &srs,
            (&m_a, &m_b),
            (&ck_a, &ck_t),
        )
        .unwrap();

        assert!(MultiExpTIPA::verify_with_structured_scalar_message(
            &v_srs,
            &ck_t,
            (&com_a, &com_t),
            &b,
            &proof
        )
        .unwrap());
    }

    #[test]
    fn scalar_inner_product_test() {
        type IP = ScalarInnerProduct<<Bls12_381 as PairingEngine>::Fr>;
        type IPC =
            IdentityCommitment<<Bls12_381 as PairingEngine>::Fr, <Bls12_381 as PairingEngine>::Fr>;
        type ScalarTIPA = TIPAWithSSM<IP, SC2, IPC, Bls12_381, Blake2b>;

        let mut rng = StdRng::seed_from_u64(0u64);
        let (srs, ck_t) = ScalarTIPA::setup(&mut rng, TEST_SIZE).unwrap();
        let (ck_a, _) = srs.get_commitment_keys();
        let v_srs = srs.get_verifier_key();
        let mut m_a = Vec::new();
        for _ in 0..TEST_SIZE {
            m_a.push(<Bls12_381 as PairingEngine>::Fr::rand(&mut rng));
        }
        let b = <<Bls12_381 as PairingEngine>::Fr>::rand(&mut rng);
        let m_b = structured_scalar_power(TEST_SIZE, &b);
        let com_a = SC2::commit(&ck_a, &m_a).unwrap();
        let t = vec![IP::inner_product(&m_a, &m_b).unwrap()];
        let com_t = IPC::commit(&vec![ck_t.clone()], &t).unwrap();

        let proof = ScalarTIPA::prove_with_structured_scalar_message(
            &srs,
            (&m_a, &m_b),
            (&ck_a, &ck_t),
        )
        .unwrap();

        assert!(ScalarTIPA::verify_with_structured_scalar_message(
            &v_srs,
            &ck_t,
            (&com_a, &com_t),
            &b,
            &proof
        )
        .unwrap());
    }
}
