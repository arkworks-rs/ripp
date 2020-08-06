use algebra::{
    curves::{PairingEngine},
    fields::{Field},
};
use digest::Digest;
use num_traits::identities::{One};
use std::{ops::MulAssign};

use crate::{
    gipa::GIPA,
    tipa::{TIPA, TIPACompatibleSetup, TIPAProof, SRS, VerifierSRS},
    Error,
};
use dh_commitments::{
    DoublyHomomorphicCommitment,
};
use inner_products::{InnerProduct};

//TODO: Properly generalize the non-committed message approach of SIPP and MIPP to GIPA
//TODO: Structured message is a special case of the non-committed message and does not rely on TIPA
//TODO: Can support structured group element messages as well as structured scalar messages

// TODO: Currently implemented by completely reusing TIPA, which means wasted prover and verifier effort on B

pub struct TIPAWithSSMProof<IP, LMC, RMC, IPC, P, D>
    where
        D: Digest,
        P: PairingEngine,
        IP: InnerProduct<
            LeftMessage = LMC::Message,
            RightMessage = RMC::Message,
            Output = IPC::Message,
        >,
        LMC: DoublyHomomorphicCommitment + TIPACompatibleSetup,
        RMC: DoublyHomomorphicCommitment<Scalar = LMC::Scalar, Message = P::Fr> + TIPACompatibleSetup,
        IPC: DoublyHomomorphicCommitment<Scalar = LMC::Scalar>,
        RMC::Message: MulAssign<LMC::Scalar>,
        IPC::Message: MulAssign<LMC::Scalar>,
        RMC::Key: MulAssign<LMC::Scalar>,
        IPC::Key: MulAssign<LMC::Scalar>,
        RMC::Output: MulAssign<LMC::Scalar>,
        IPC::Output: MulAssign<LMC::Scalar>,
{
    tipa_proof: TIPAProof<IP, LMC, RMC, IPC, P, D>,
    com_b: RMC::Output, //TODO: Needed because reusing TIPA
}



impl<IP, LMC, RMC, IPC, P, D> TIPA<IP, LMC, RMC, IPC, P, D>
    where
        D: Digest,
        P: PairingEngine,
        IP: InnerProduct<
            LeftMessage = LMC::Message,
            RightMessage = RMC::Message,
            Output = IPC::Message,
        >,
        LMC: DoublyHomomorphicCommitment<Scalar = P::Fr, Key = P::G2Projective>
        + TIPACompatibleSetup,
        RMC: DoublyHomomorphicCommitment<Scalar = LMC::Scalar, Key = P::G1Projective, Message = P::Fr>
        + TIPACompatibleSetup,
        IPC: DoublyHomomorphicCommitment<Scalar = LMC::Scalar>,
        LMC::Message: MulAssign<P::Fr>,
        RMC::Message: MulAssign<P::Fr>,
        IPC::Message: MulAssign<P::Fr>,
        IPC::Key: MulAssign<P::Fr>,
        LMC::Output: MulAssign<P::Fr>,
        RMC::Output: MulAssign<P::Fr>,
        IPC::Output: MulAssign<P::Fr>,
{

    pub fn prove_with_structured_scalar_message(
        srs: &SRS<P>,
        values: (&[IP::LeftMessage], &[IP::RightMessage]),
        ck: (&[LMC::Key], &[RMC::Key], &IPC::Key),
    ) -> Result<TIPAWithSSMProof<IP, LMC, RMC, IPC, P, D>, Error> {
        Ok(TIPAWithSSMProof {
            tipa_proof: TIPA::prove(srs, values, ck)?,
            com_b: RMC::commit(ck.1, values.1)?,
        })
    }

    pub fn verify_with_structured_scalar_message(
        v_srs: &VerifierSRS<P>,
        ck_t: &IPC::Key,
        com: (&LMC::Output, &IPC::Output),
        scalar_b: &P::Fr,
        proof: &TIPAWithSSMProof<IP, LMC, RMC, IPC, P, D>,
    ) -> Result<bool, Error> {
        let tipa_valid = TIPA::verify(v_srs, ck_t, (com.0, &proof.com_b, com.1), &proof.tipa_proof)?;

        // Check final scalar
        //TODO: repeating gathering of transcript from TIPA verify
        let (_, transcript) = GIPA::verify_recursive_challenge_transcript((com.0, &proof.com_b, com.1), &proof.tipa_proof.gipa_proof)?;
        let mut power_2_b = scalar_b.clone();
        let mut product_form = Vec::new();
        for x in transcript.iter() {
            product_form.push(<P::Fr>::one() + &(x.inverse().unwrap() * &power_2_b));
            power_2_b *= power_2_b;
        }
        let final_b = product_form.iter().product::<P::Fr>();
        let final_b_valid = final_b == proof.tipa_proof.gipa_proof.r_base.1;

        Ok(tipa_valid && final_b_valid)
    }

}

pub fn structured_scalar_power<F: Field>(
    num: usize,
    s: &F,
) -> Vec<F> {
    let mut powers= vec![F::one()];
    for i in 1..num {
        powers.push(powers[i-1] * s);
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
        afgho16::{AFGHOCommitmentG1},
        identity::IdentityCommitment,
        pedersen::PedersenCommitment,
        random_generators,
    };
    use inner_products::{
        InnerProduct, MultiexponentiationInnerProduct,
        ScalarInnerProduct,
    };

    type GC1 = AFGHOCommitmentG1<Bls12_381>;
    type SC1 = PedersenCommitment<<Bls12_381 as PairingEngine>::G1Projective>;
    type SC2 = PedersenCommitment<<Bls12_381 as PairingEngine>::G2Projective>;

    const TEST_SIZE: usize = 8;

    #[test]
    fn multiexponentiation_inner_product_test() {
        type IP = MultiexponentiationInnerProduct<<Bls12_381 as PairingEngine>::G1Projective>;
        type IPC = IdentityCommitment<
            <Bls12_381 as PairingEngine>::G1Projective,
            <Bls12_381 as PairingEngine>::Fr,
        >;
        type MultiExpTIPA = TIPA<IP, GC1, SC1, IPC, Bls12_381, Blake2b>;

        let mut rng = StdRng::seed_from_u64(0u64);
        let (srs, ck_t) = MultiExpTIPA::setup(&mut rng, TEST_SIZE).unwrap();
        let (ck_a, ck_b) = srs.get_commitment_keys();
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
            (&ck_a, &ck_b, &ck_t),
        )
            .unwrap();

        assert!(
            MultiExpTIPA::verify_with_structured_scalar_message(&v_srs, &ck_t, (&com_a, &com_t), &b, &proof).unwrap()
        );
    }

    #[test]
    fn scalar_inner_product_test() {
        type IP = ScalarInnerProduct<<Bls12_381 as PairingEngine>::Fr>;
        type IPC =
        IdentityCommitment<<Bls12_381 as PairingEngine>::Fr, <Bls12_381 as PairingEngine>::Fr>;
        type ScalarTIPA = TIPA<IP, SC2, SC1, IPC, Bls12_381, Blake2b>;

        let mut rng = StdRng::seed_from_u64(0u64);
        let (srs, ck_t) = ScalarTIPA::setup(&mut rng, TEST_SIZE).unwrap();
        let (ck_a, ck_b) = srs.get_commitment_keys();
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
            (&ck_a, &ck_b, &ck_t),
        )
            .unwrap();

        assert!(
            ScalarTIPA::verify_with_structured_scalar_message(&v_srs, &ck_t, (&com_a, &com_t), &b, &proof).unwrap()
        );
    }
}
