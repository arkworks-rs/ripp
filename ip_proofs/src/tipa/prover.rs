use ark_ec::pairing::Pairing;
use ark_std::{end_timer, start_timer};
use digest::Digest;

use crate::{
    gipa::{Instance, Witness, GIPA},
    Error,
};

use super::*;

impl<P, D> TIPA<P, D>
where
    D: Digest,
    P: Pairing,
{
    pub fn prove<'a>(
        pk: &ProverKey<P>,
        instance: &Instance<TIPPCommitment<P>>,
        witness: &Witness<PairingInnerProduct<P>>,
    ) -> Result<Proof<P>, Error> {
        let prover_time = start_timer!(|| "TIPA::Prove");
        let (gipa_proof, aux) = <GIPA<IP<P>, IPC<P>, D>>::prove_helper(&pk.pk, instance, witness)?;
        let final_ck = aux.final_ck;

        let kzg_point = Self::compute_kzg_challenge(&final_ck, &aux.challenges)?;

        let final_ck_time = start_timer!(|| "FinalCK Correctness Proof");
        let final_ck_proof =
            TIPPCommitment::prove_final_ck(&pk, &aux.challenges, instance.twist, kzg_point);
        end_timer!(final_ck_time);

        let gipa_verify_time = start_timer!(|| "GIPA::Verify");
        debug_assert!(GIPA::<IP<P>, IPC<P>, D>::verify(
            &pk.pk.vk(),
            instance,
            &gipa_proof
        )?);
        end_timer!(gipa_verify_time);
        let final_ck_consistency_time = start_timer!(|| "FinalCK Consistency Check");
        debug_assert_eq!(
            final_ck,
            GIPA::<IP<P>, IPC<P>, D>::compute_final_ck(
                &pk.pk.ck,
                &aux.challenges,
                &instance.twist
            )?
            .try_into()
            .unwrap()
        );
        end_timer!(final_ck_consistency_time);

        let proof = Proof {
            gipa_proof,
            final_ck: final_ck.into(),
            final_ck_proof,
        };
        debug_assert!(Self::verify(&pk.vk(), instance, &proof).unwrap());
        end_timer!(prover_time);

        Ok(proof)
    }
}
