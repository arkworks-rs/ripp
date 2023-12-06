use ark_ec::pairing::Pairing;
use digest::Digest;

use crate::{
    gipa::{Instance, GIPA},
    ip_commitment::{IPCommKey, IPCommitment},
    Error,
};

use super::*;

impl<P, D> TIPA<P, D>
where
    D: Digest,
    P: Pairing,
{
    pub fn verify<'a>(
        vk: &VerifierKey<P>,
        instance: &Instance<IPC<P>>,
        proof: &Proof<P>,
    ) -> Result<bool, Error> {
        let Instance {
            output,
            commitment: mut com,
            twist,
            ..
        } = instance;

        // We assume that the initial commitment in the instance *does* not commit to the output.
        // That is, the initial commitment is a commitment to `None`.
        com += IPC::commit_only_ip(&vk.ck_for_ip, *output)?;
        let (final_commitment, challenges) =
            GIPA::<_, _, D>::verify_recursive_challenge_transcript(&com, &proof.gipa_proof)?;

        let kzg_point = Self::compute_kzg_challenge(&proof.final_ck, &challenges)?;

        let final_ck_valid = TIPPCommitment::verify_final_ck(
            vk,
            &proof.final_ck,
            &challenges,
            *twist,
            kzg_point,
            &proof.final_ck_proof,
        );
        let final_ck: IPCommKey<_> = proof.final_ck.into();

        // Verify final inner product commitment
        let final_cm_valid = GIPA::<_, _, D>::verify_final_commitment(
            &final_ck,
            &final_commitment,
            &proof.gipa_proof,
        )?;

        Ok(final_ck_valid & final_cm_valid)
    }
}
