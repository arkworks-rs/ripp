use ark_ec::pairing::Pairing;
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
        let (gipa_proof, aux) = <GIPA<IP<P>, IPC<P>, D>>::prove_helper(&pk.pk, instance, witness)?;
        let final_ck = aux.final_ck;
        // Prove final commitment keys are wellformed
        let kzg_point = Self::compute_kzg_challenge(&final_ck, &aux.challenges)?;

        let final_ck_proof =
            TIPPCommitment::prove_final_ck(&pk, &aux.challenges, instance.twist, kzg_point);

        Ok(Proof {
            gipa_proof,
            final_ck: final_ck.into(),
            final_ck_proof,
        })
    }
}
