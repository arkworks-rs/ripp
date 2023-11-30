use ark_ff::{Field, One};
use ark_std::borrow::Cow;
use digest::Digest;

use crate::{
    ip_commitment::{Commitment, IPCommKey, IPCommitment, Scalar},
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
        vk: &VerifierKey<'a, IPC>,
        instance: &Instance<IPC>,
        proof: &Proof<IP, IPC, D>,
    ) -> Result<bool, Error> {
        // TODO: output should appear in checks somewhere...
        let Instance {
            size,
            output,
            commitment: mut com,
            random_challenge,
        } = instance;
        if !vk.ck.ck_a.len().is_power_of_two() || vk.ck.ck_a.len() != vk.ck.ck_b.len() {
            // Power of 2 length
            return Err(Box::new(InnerProductArgumentError::MessageLengthInvalid(
                vk.ck.ck_a.len(),
                vk.ck.ck_b.len(),
            )));
        }
        // We assume that the initial commitment in the instance *does* not commit to the output.
        // That is, the initial commitment is a commitment to `None`.
        com += IPC::commit_only_ip(&vk.ck, *output)?;

        // Calculate base commitment and transcript
        let (base_com, transcript) = Self::_compute_recursive_challenges(&com, proof)?;
        // Calculate base commitment keys
        let ck_base = Self::_compute_final_commitment_keys(&vk.ck, &transcript, random_challenge)?;
        // Verify base commitment
        Self::_verify_base_commitment(&ck_base, &base_com, proof)
    }

    // Helper function used to calculate recursive challenges from proof execution (transcript in reverse)
    pub fn verify_recursive_challenge_transcript(
        com: &Commitment<IPC>,
        proof: &Proof<IP, IPC, D>,
    ) -> Result<(Commitment<IPC>, Vec<Scalar<IPC>>), Error> {
        Self::_compute_recursive_challenges(com, proof)
    }

    fn _compute_recursive_challenges(
        com: &Commitment<IPC>,
        proof: &Proof<IP, IPC, D>,
    ) -> Result<(Commitment<IPC>, Vec<Scalar<IPC>>), Error> {
        let mut com = com.clone();
        let mut r_transcript: Vec<Scalar<IPC>> = Vec::new();
        let default_transcript_entry = Scalar::<IPC>::one();
        for (com_1, com_2) in proof.r_commitment_steps.iter().rev() {
            // Fiat-Shamir challenge
            let transcript = r_transcript.last().unwrap_or(&default_transcript_entry);
            let (c, c_inv) = Self::compute_challenge(transcript, com_1, com_2)?;
            com += com_1.clone() * c + com_2.clone() * c_inv;
            r_transcript.push(c);
        }
        r_transcript.reverse();
        Ok((com, r_transcript))
    }

    pub(crate) fn _compute_final_commitment_keys<'a>(
        ck: &IPCommKey<'a, IPC>,
        transcript: &[Scalar<IPC>],
        &random_challenge: &Scalar<IPC>,
    ) -> Result<IPCommKey<'a, IPC>, Error> {
        // Calculate base commitment keys
        assert!(ck.ck_a.len().is_power_of_two());

        let mut ck_a_poly_coeffs = vec![Scalar::<IPC>::one()];
        let mut ck_b_poly_coeffs = vec![Scalar::<IPC>::one()];
        for (i, c) in transcript.iter().enumerate() {
            let r = random_challenge.pow([(2_u64).pow(i as u32)]);
            let c_inv = (*c * r).inverse().unwrap();
            for j in 0..((2_usize).pow(i as u32)) {
                ck_a_poly_coeffs.push(ck_a_poly_coeffs[j] * &c_inv);
                ck_b_poly_coeffs.push(ck_b_poly_coeffs[j] * *c);
            }
        }
        assert_eq!(ck_a_poly_coeffs.len(), ck.ck_a.len());

        let ck_a_base = IPC::left_key_msm(&ck.ck_a, &ck_a_poly_coeffs)?;
        let ck_b_base = IPC::right_key_msm(&ck.ck_b, &ck_b_poly_coeffs)?;

        Ok(IPCommKey::new(
            Cow::Owned(vec![ck_a_base]),
            Cow::Owned(vec![ck_b_base]),
            ck.ck_t.clone(),
        ))
    }

    pub(crate) fn _verify_base_commitment<'a>(
        base_ck: &IPCommKey<'a, IPC>,
        base_com: &Commitment<IPC>,
        proof: &Proof<IP, IPC, D>,
    ) -> Result<bool, Error> {
        let l = [proof.r_base.0.clone()];
        let r = [proof.r_base.1.clone()];
        Ok(IPC::verify(&base_ck, &l, &r, &base_com)?)
    }
}
