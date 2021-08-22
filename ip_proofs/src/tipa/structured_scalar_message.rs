use crate::{
    gipa::{GIPAProof, GIPA},
    tipa::{
        prove_commitment_key_kzg_opening, structured_generators_scalar_power,
        verify_commitment_key_g2_kzg_opening, TIPACompatibleSetup, VerifierSRS, SRS,
    },
    util::TranscriptProtocol,
    Error,
};
use ark_dh_commitments::{identity::HomomorphicPlaceholderValue, DoublyHomomorphicCommitment};
use ark_inner_products::InnerProduct;

use std::{marker::PhantomData, ops::MulAssign};

use ark_ec::{group::Group, PairingEngine, ProjectiveCurve};
use ark_ff::{Field, One, PrimeField, UniformRand, Zero};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, Read, SerializationError, Write};
use ark_std::{cfg_iter, rand::Rng};
use ark_std::{end_timer, start_timer};
use merlin::Transcript;

#[cfg(feature = "parallel")]
use rayon::prelude::*;

const GIPA_SSM_DOMAIN_SEP: &[u8] = b"ip_proofs-v0.3-GIPA_SSM";
const TIPA_SSM_DOMAIN_SEP: &[u8] = b"ip_proofs-v0.3-TIPA_SSM";

//TODO: Properly generalize the non-committed message approach of SIPP and MIPP to GIPA
//TODO: Structured message is a special case of the non-committed message and does not rely on TIPA
//TODO: Can support structured group element messages as well as structured scalar messages

// Use placeholder commitment to commit to vector in clear during GIPA execution
#[derive(Clone)]
pub struct SSMPlaceholderCommitment<F: PrimeField> {
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

    //TODO: Doesn't include message which means scalar b not included in generating challenges
    fn commit(_k: &[Self::Key], _m: &[Self::Message]) -> Result<Self::Output, Error> {
        Ok(F::zero())
    }
}

pub struct GIPAWithSSM<IP, LMC, IPC> {
    _inner_product: PhantomData<IP>,
    _left_commitment: PhantomData<LMC>,
    _inner_product_commitment: PhantomData<IPC>,
}

impl<IP, LMC, IPC> GIPAWithSSM<IP, LMC, IPC>
where
    IP: InnerProduct<LeftMessage = LMC::Message, RightMessage = LMC::Scalar, Output = IPC::Message>,
    LMC: DoublyHomomorphicCommitment,
    IPC: DoublyHomomorphicCommitment<Scalar = LMC::Scalar>,
    IPC::Message: MulAssign<LMC::Scalar>,
    IPC::Key: MulAssign<LMC::Scalar>,
    IPC::Output: MulAssign<LMC::Scalar>,
{
    pub fn setup<R: Rng>(rng: &mut R, size: usize) -> Result<(Vec<LMC::Key>, IPC::Key), Error> {
        Ok((LMC::setup(rng, size)?, IPC::setup(rng, 1)?.pop().unwrap()))
    }

    /// Proves an inner product of `m_a` and `m_b`, where `m_b` is of the form
    pub fn prove_with_structured_scalar_message(
        transcript: &mut Transcript,
        values: (&[IP::LeftMessage], &[IP::RightMessage]),
        ck: (&[LMC::Key], &IPC::Key),
    ) -> Result<GIPAProof<IP, LMC, SSMPlaceholderCommitment<LMC::Scalar>, IPC>, Error> {
        // Domain-separate this protocol
        transcript.append_message(b"dom-sep", GIPA_SSM_DOMAIN_SEP);

        let (proof, _) =
            <GIPA<IP, LMC, SSMPlaceholderCommitment<LMC::Scalar>, IPC>>::prove_with_aux(
                transcript,
                values,
                (
                    ck.0,
                    &vec![HomomorphicPlaceholderValue {}; values.1.len()],
                    &vec![ck.1.clone()],
                ),
            )?;
        Ok(proof)
    }

    pub fn verify_with_structured_scalar_message(
        transcript: &mut Transcript,
        ck: (&[LMC::Key], &IPC::Key),
        com: (&LMC::Output, &IPC::Output),
        scalar_b: &LMC::Scalar,
        proof: &GIPAProof<IP, LMC, SSMPlaceholderCommitment<LMC::Scalar>, IPC>,
    ) -> Result<bool, Error> {
        // Domain-separate this protocol
        transcript.append_message(b"dom-sep", GIPA_SSM_DOMAIN_SEP);

        // Calculate base commitments and recursive transcript
        //TODO: Scalar b not included in generating challenges
        let (base_com, transcript) = GIPA::verify_recursive_challenge_transcript(
            transcript,
            (com.0, &LMC::Scalar::zero(), com.1),
            proof,
        )?;
        // Calculate base commitment keys
        let (ck_a_base, ck_b_base) = GIPA::<IP, LMC, SSMPlaceholderCommitment<LMC::Scalar>, IPC>::_compute_final_commitment_keys(
            (&ck.0, &vec![HomomorphicPlaceholderValue {}; ck.0.len()], &ck.1),
            &transcript,
        )?;
        // Verify base commitment
        let gipa_valid = GIPA::_verify_base_commitment(
            (&ck_a_base, &ck_b_base, &vec![ck.1.clone()]),
            base_com.clone(),
            proof,
        )?;

        // Compute final scalar
        let mut power_2_b = scalar_b.clone();
        let mut product_form = Vec::new();
        for x in transcript.iter() {
            product_form.push(<LMC::Scalar>::one() + &(x.inverse().unwrap() * &power_2_b));
            power_2_b *= &power_2_b.clone();
        }
        let b_base = cfg_iter!(product_form).product::<LMC::Scalar>();

        // Verify base inner product commitment
        let (com_a, _, com_t) = base_com;
        let a_base = vec![proof.r_base.0.clone()];
        let t_base = vec![IP::inner_product(&a_base, &vec![b_base])?];
        let base_valid = LMC::verify(&vec![ck_a_base.clone()], &a_base, &com_a)?
            && IPC::verify(&vec![ck.1.clone()], &t_base, &com_t)?;

        Ok(gipa_valid && base_valid)
    }
}

pub struct TIPAWithSSM<IP, LMC, IPC, P> {
    _inner_product: PhantomData<IP>,
    _left_commitment: PhantomData<LMC>,
    _inner_product_commitment: PhantomData<IPC>,
    _pair: PhantomData<P>,
}

#[derive(CanonicalSerialize, CanonicalDeserialize)]
pub struct TIPAWithSSMProof<IP, LMC, IPC, P>
where
    P: PairingEngine,
    IP: InnerProduct<LeftMessage = LMC::Message, RightMessage = LMC::Scalar, Output = IPC::Message>,
    LMC: DoublyHomomorphicCommitment<Scalar = P::Fr, Key = P::G2Projective> + TIPACompatibleSetup,
    IPC: DoublyHomomorphicCommitment<Scalar = LMC::Scalar>,
    IPC::Message: MulAssign<LMC::Scalar>,
    IPC::Key: MulAssign<LMC::Scalar>,
    IPC::Output: MulAssign<LMC::Scalar>,
    LMC::Message: MulAssign<P::Fr>,
    LMC::Output: MulAssign<P::Fr>,
{
    gipa_proof: GIPAProof<IP, LMC, SSMPlaceholderCommitment<LMC::Scalar>, IPC>,
    final_ck: LMC::Key,
    final_ck_proof: P::G2Projective,
    _pairing: PhantomData<P>,
}

impl<IP, LMC, IPC, P> Clone for TIPAWithSSMProof<IP, LMC, IPC, P>
where
    P: PairingEngine,
    IP: InnerProduct<LeftMessage = LMC::Message, RightMessage = LMC::Scalar, Output = IPC::Message>,
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

impl<IP, LMC, IPC, P> TIPAWithSSM<IP, LMC, IPC, P>
where
    P: PairingEngine,
    IP: InnerProduct<LeftMessage = LMC::Message, RightMessage = LMC::Scalar, Output = IPC::Message>,
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
        transcript: &mut Transcript,
        srs: &SRS<P>,
        values: (&[IP::LeftMessage], &[IP::RightMessage]),
        ck: (&[LMC::Key], &IPC::Key),
    ) -> Result<TIPAWithSSMProof<IP, LMC, IPC, P>, Error> {
        // Domain-separate this protocol
        transcript.append_message(b"dom-sep", TIPA_SSM_DOMAIN_SEP);

        // Run GIPA
        let gipa = start_timer!(|| "GIPA");
        let (proof, aux) = <GIPA<IP, LMC, SSMPlaceholderCommitment<P::Fr>, IPC>>::prove_with_aux(
            transcript,
            values,
            (
                ck.0,
                &vec![HomomorphicPlaceholderValue {}; values.1.len()],
                &vec![ck.1.clone()],
            ),
        )?;
        end_timer!(gipa);

        // Prove final commitment key is wellformed
        let ck_kzg = start_timer!(|| "Prove commitment key");
        let (ck_a_final, _) = aux.ck_base;
        let r_challenges = aux.r_challenges;
        let r_chal_inverses = cfg_iter!(r_challenges)
            .map(|x| x.inverse().unwrap())
            .collect();

        // KZG challenge point
        let kzg_chal = transcript.challenge_scalar(b"kzg_chal");

        // Complete KZG proof
        let ck_a_kzg_opening = prove_commitment_key_kzg_opening(
            &srs.h_beta_powers,
            &r_chal_inverses,
            &<P::Fr>::one(),
            &kzg_chal,
        )?;
        end_timer!(ck_kzg);

        Ok(TIPAWithSSMProof {
            gipa_proof: proof,
            final_ck: ck_a_final,
            final_ck_proof: ck_a_kzg_opening,
            _pairing: PhantomData,
        })
    }

    pub fn verify_with_structured_scalar_message(
        transcript: &mut Transcript,
        v_srs: &VerifierSRS<P>,
        ck_t: &IPC::Key,
        com: (&LMC::Output, &IPC::Output),
        scalar_b: &P::Fr,
        proof: &TIPAWithSSMProof<IP, LMC, IPC, P>,
    ) -> Result<bool, Error> {
        // Domain-separate this protocol
        transcript.append_message(b"dom-sep", TIPA_SSM_DOMAIN_SEP);

        let (base_com, r_challenges) = GIPA::verify_recursive_challenge_transcript(
            transcript,
            (com.0, scalar_b, com.1),
            &proof.gipa_proof,
        )?;
        let r_chal_inverses = cfg_iter!(r_challenges)
            .map(|x| x.inverse().unwrap())
            .collect();

        let ck_a_final = &proof.final_ck;
        let ck_a_proof = &proof.final_ck_proof;

        // KZG challenge point
        let kzg_chal = transcript.challenge_scalar(b"kzg_chal");

        // Check commitment key
        let ck_a_valid = verify_commitment_key_g2_kzg_opening(
            v_srs,
            &ck_a_final,
            &ck_a_proof,
            &r_chal_inverses,
            &<P::Fr>::one(),
            &kzg_chal,
        )?;

        // Compute final scalar
        let mut power_2_b = scalar_b.clone();
        let mut product_form = Vec::new();
        for x in r_challenges.iter() {
            product_form.push(<P::Fr>::one() + &(x.inverse().unwrap() * &power_2_b));
            power_2_b *= &power_2_b.clone();
        }
        let b_base = cfg_iter!(product_form).product::<P::Fr>();

        // Verify base inner product commitment
        let (com_a, _, com_t) = base_com;
        let a_base = vec![proof.gipa_proof.r_base.0.clone()];
        let t_base = vec![IP::inner_product(&a_base, &vec![b_base])?];
        let base_valid = LMC::verify(&vec![ck_a_final.clone()], &a_base, &com_a)?
            && IPC::verify(&vec![ck_t.clone()], &t_base, &com_t)?;

        Ok(ck_a_valid && base_valid)
    }
}

/// Returns `(1, s^1, s^2, ..., s^(num-1))`
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
    use ark_bls12_381::Bls12_381;

    use ark_dh_commitments::{
        afgho16::AFGHOCommitmentG1, identity::IdentityCommitment, pedersen::PedersenCommitment,
        random_generators,
    };
    use ark_inner_products::{InnerProduct, MultiexponentiationInnerProduct, ScalarInnerProduct};

    type GC1 = AFGHOCommitmentG1<Bls12_381>;
    type SC1 = PedersenCommitment<<Bls12_381 as PairingEngine>::G1Projective>;

    const TEST_SIZE: usize = 8;

    #[test]
    fn tipa_ssm_multiexponentiation_inner_product_test() {
        type IP = MultiexponentiationInnerProduct<<Bls12_381 as PairingEngine>::G1Projective>;
        type IPC = IdentityCommitment<
            <Bls12_381 as PairingEngine>::G1Projective,
            <Bls12_381 as PairingEngine>::Fr,
        >;
        type MultiExpTIPA = TIPAWithSSM<IP, GC1, IPC, Bls12_381>;

        let mut rng = ark_std::test_rng();
        let (srs, ck_t) = MultiExpTIPA::setup(&mut rng, TEST_SIZE).unwrap();
        let (ck_a, _) = srs.get_commitment_keys();
        let v_srs = srs.get_verifier_key();
        let m_a = random_generators(&mut rng, TEST_SIZE);
        let b = <<Bls12_381 as PairingEngine>::Fr>::rand(&mut rng);
        let m_b = structured_scalar_power(TEST_SIZE, &b);
        let com_a = GC1::commit(&ck_a, &m_a).unwrap();
        let t = vec![IP::inner_product(&m_a, &m_b).unwrap()];
        let com_t = IPC::commit(&vec![ck_t.clone()], &t).unwrap();

        let mut proof_transcript = Transcript::new(b"TIPA_SSM-test");
        let proof = MultiExpTIPA::prove_with_structured_scalar_message(
            &mut proof_transcript,
            &srs,
            (&m_a, &m_b),
            (&ck_a, &ck_t),
        )
        .unwrap();

        let mut verif_transcript = Transcript::new(b"TIPA_SSM-test");
        assert!(MultiExpTIPA::verify_with_structured_scalar_message(
            &mut verif_transcript,
            &v_srs,
            &ck_t,
            (&com_a, &com_t),
            &b,
            &proof
        )
        .unwrap());
    }

    #[test]
    fn gipa_ssm_scalar_inner_product_test() {
        type IP = ScalarInnerProduct<<Bls12_381 as PairingEngine>::Fr>;
        type IPC =
            IdentityCommitment<<Bls12_381 as PairingEngine>::Fr, <Bls12_381 as PairingEngine>::Fr>;
        type ScalarGIPA = GIPAWithSSM<IP, SC1, IPC>;

        let mut rng = ark_std::test_rng();
        let (ck_a, ck_t) = ScalarGIPA::setup(&mut rng, TEST_SIZE).unwrap();
        let mut m_a = Vec::new();
        for _ in 0..TEST_SIZE {
            m_a.push(<Bls12_381 as PairingEngine>::Fr::rand(&mut rng));
        }
        let b = <<Bls12_381 as PairingEngine>::Fr>::rand(&mut rng);
        let m_b = structured_scalar_power(TEST_SIZE, &b);
        let com_a = SC1::commit(&ck_a, &m_a).unwrap();
        let t = vec![IP::inner_product(&m_a, &m_b).unwrap()];
        let com_t = IPC::commit(&vec![ck_t.clone()], &t).unwrap();

        let mut proof_transcript = Transcript::new(b"TIPA_SSM-test");
        let proof = ScalarGIPA::prove_with_structured_scalar_message(
            &mut proof_transcript,
            (&m_a, &m_b),
            (&ck_a, &ck_t),
        )
        .unwrap();

        let mut verif_transcript = Transcript::new(b"TIPA_SSM-test");
        assert!(ScalarGIPA::verify_with_structured_scalar_message(
            &mut verif_transcript,
            (&ck_a, &ck_t),
            (&com_a, &com_t),
            &b,
            &proof
        )
        .unwrap());
    }
}
