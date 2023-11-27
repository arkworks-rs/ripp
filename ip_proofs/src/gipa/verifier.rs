use ark_ff::{Field, One};
use ark_std::borrow::Cow;
use digest::Digest;

use crate::{
    ip_commitment::{IPCommKey, IPCommitment, Scalar},
    Error, InnerProductArgumentError,
};
use ark_inner_products::InnerProduct;

use super::{data_structures::*, GIPA};

impl<IP, IPC, D> GIPA<IP, IPC, D>
where
    D: Digest,
    IP: InnerProduct,
    IPC: IPCommitment<IP = IP>,
{
    pub fn verify<'a>(
        ck: &IPCommKey<'a, IPC>,
        instance: &Instance<IPC>,
        proof: &Proof<IP, IPC, D>,
    ) -> Result<bool, Error> {
        let Instance {
            size,
            output,
            commitment: com,
            random_challenge: c,
        } = instance;
        if !ck.ck_a.len().is_power_of_two() || ck.ck_a.len() != ck.ck_b.len() {
            // Power of 2 length
            return Err(Box::new(InnerProductArgumentError::MessageLengthInvalid(
                ck.ck_a.len(),
                ck.ck_b.len(),
            )));
        }
        // Calculate base commitment and transcript
        let (base_com, transcript) = Self::_compute_recursive_challenges(com, proof)?;
        // Calculate base commitment keys
        let ck_base = Self::_compute_final_commitment_keys(ck, &transcript)?;
        // Verify base commitment
        Self::_verify_base_commitment(&ck_base, &base_com, proof)
    }

    // Helper function used to calculate recursive challenges from proof execution (transcript in reverse)
    pub fn verify_recursive_challenge_transcript(
        com: &IPC::Commitment,
        proof: &Proof<IP, IPC, D>,
    ) -> Result<(IPC::Commitment, Vec<Scalar<IPC>>), Error> {
        Self::_compute_recursive_challenges(com, proof)
    }

    fn _compute_recursive_challenges(
        com: &IPC::Commitment,
        proof: &Proof<IP, IPC, D>,
    ) -> Result<(IPC::Commitment, Vec<Scalar<IPC>>), Error> {
        let mut com = com.clone();
        let mut r_transcript: Vec<Scalar<IPC>> = Vec::new();
        let default_transcript_entry = Default::default();
        for (com_1, com_2) in proof.r_commitment_steps.iter().rev() {
            // Fiat-Shamir challenge
            let transcript = r_transcript.last().unwrap_or(&default_transcript_entry);
            let (c, c_inv) = Self::compute_challenge(transcript, com_1, com_2)?;
            com = com + (com_1.clone() * c + com_2.clone() * c_inv);
            r_transcript.push(c);
        }
        r_transcript.reverse();
        Ok((com, r_transcript))
    }

    pub(crate) fn _compute_final_commitment_keys<'a>(
        ck: &IPCommKey<'a, IPC>,
        transcript: &[Scalar<IPC>],
    ) -> Result<IPCommKey<'a, IPC>, Error> {
        // Calculate base commitment keys
        assert!(ck.ck_a.len().is_power_of_two());

        let mut ck_a_agg_challenge_exponents = vec![Scalar::<IPC>::one()];
        let mut ck_b_agg_challenge_exponents = vec![Scalar::<IPC>::one()];
        for (i, c) in transcript.iter().enumerate() {
            let c_inv = c.inverse().unwrap();
            for j in 0..(2_usize).pow(i as u32) {
                ck_a_agg_challenge_exponents.push(ck_a_agg_challenge_exponents[j] * &c_inv);
                ck_b_agg_challenge_exponents.push(ck_b_agg_challenge_exponents[j] * c);
            }
        }
        assert_eq!(ck_a_agg_challenge_exponents.len(), ck.ck_a.len());
        let ck_a_base = IPC::left_key_msm(&ck.ck_a, &ck_a_agg_challenge_exponents)?;
        let ck_b_base = IPC::right_key_msm(&ck.ck_b, &ck_b_agg_challenge_exponents)?;

        Ok(IPCommKey::new(
            Cow::Owned(vec![ck_a_base]),
            Cow::Owned(vec![ck_b_base]),
            ck.ck_t.clone(),
        ))
    }

    pub(crate) fn _verify_base_commitment<'a>(
        base_ck: &IPCommKey<'a, IPC>,
        base_com: &IPC::Commitment,
        proof: &Proof<IP, IPC, D>,
    ) -> Result<bool, Error> {
        let a_base = [proof.r_base.0.clone()];
        let b_base = [proof.r_base.1.clone()];
        let t_base = IP::inner_product(&a_base, &b_base)?;
        dbg!(&base_com);

        Ok(IPC::verify(&base_ck, &a_base, &b_base, &t_base, &base_com)?)
    }
}
