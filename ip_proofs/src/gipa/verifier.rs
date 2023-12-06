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
    pub fn verify(
        vk: &VerifierKey<IPC>,
        instance: &Instance<IPC>,
        proof: &Proof<IP, IPC>,
    ) -> Result<bool, Error> {
        let Instance {
            output,
            commitment: mut com,
            twist,
            ..
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
        let (final_commitment, challenges) =
            Self::verify_recursive_challenge_transcript(&com, proof)?;
        // Calculate base commitment keys
        let final_ck = Self::compute_final_ck(&vk.ck, &challenges, twist)?;
        // Verify base commitment
        Self::verify_final_commitment(&final_ck, &final_commitment, proof)
    }

    // Helper function used to calculate recursive challenges from proof execution (transcript in reverse)
    pub fn verify_recursive_challenge_transcript(
        com: &Commitment<IPC>,
        proof: &Proof<IP, IPC>,
    ) -> Result<(Commitment<IPC>, Vec<Scalar<IPC>>), Error> {
        Self::_compute_recursive_challenges(com, proof)
    }

    fn _compute_recursive_challenges(
        com: &Commitment<IPC>,
        proof: &Proof<IP, IPC>,
    ) -> Result<(Commitment<IPC>, Vec<Scalar<IPC>>), Error> {
        let mut com = com.clone();
        let mut challenges: Vec<Scalar<IPC>> = Vec::new();
        let default_transcript_entry = Scalar::<IPC>::one();
        for (com_1, com_2) in proof.commitments.iter().rev() {
            // Fiat-Shamir challenge
            let transcript = challenges.last().unwrap_or(&default_transcript_entry);
            let (c, c_inv) = Self::compute_challenge(transcript, com_1, com_2)?;
            com += com_1.clone() * c + com_2.clone() * c_inv;
            challenges.push(c);
        }
        challenges.reverse();
        Ok((com, challenges))
    }

    pub(crate) fn compute_final_ck<'a>(
        ck: &IPCommKey<'a, IPC>,
        challenges: &[Scalar<IPC>],
        &twist: &Scalar<IPC>,
    ) -> Result<IPCommKey<'a, IPC>, Error> {
        // Calculate base commitment keys
        assert!(ck.ck_a.len().is_power_of_two());

        let mut ck_a_poly_coeffs = vec![Scalar::<IPC>::one()];
        let mut ck_b_poly_coeffs = vec![Scalar::<IPC>::one()];
        let mut challenges_inv = challenges.to_vec();
        challenges_inv.push(twist);
        ark_ff::batch_inversion(&mut challenges_inv);
        let twist_inv = challenges_inv.pop().unwrap();

        for (i, (&c, c_inv)) in challenges.iter().zip(challenges_inv).enumerate() {
            let c_inv = c_inv * twist_inv.pow([(2_u64).pow(i as u32)]);
            for j in 0..((2_usize).pow(i as u32)) {
                ck_a_poly_coeffs.push(ck_a_poly_coeffs[j] * c_inv);
                ck_b_poly_coeffs.push(ck_b_poly_coeffs[j] * c);
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

    pub(crate) fn verify_final_commitment<'a>(
        final_ck: &IPCommKey<'a, IPC>,
        final_commitment: &Commitment<IPC>,
        proof: &Proof<IP, IPC>,
    ) -> Result<bool, Error> {
        let l = [proof.final_msg.0];
        let r = [proof.final_msg.1];
        Ok(IPC::verify(&final_ck, &l, &r, &final_commitment)?)
    }
}
