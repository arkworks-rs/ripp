use ark_dh_commitments::{
    identity::{IdentityCommitment, IdentityOutput, PlaceholderKey},
    DoublyHomomorphicCommitment, Error,
};
use ark_ec::{
    pairing::{Pairing, PairingOutput},
    CurveGroup, VariableBaseMSM,
};
use ark_inner_products::{multi_pairing, PairingInnerProduct};
use ark_std::{end_timer, marker::PhantomData, rand::Rng, start_timer};

use crate::{
    ip_commitment::{
        Commitment, FinalIPCommKey, IPCommKey, IPCommitment, LeftMessage, OutputMessage,
        RightMessage,
    },
    tipa::{ProverKey, VerifierKey},
};

mod data_structures;
pub use data_structures::*;

use self::kzg::verify::{verify_left_key, verify_right_key};

pub mod kzg;

pub struct TIPPCommitment<E: Pairing>(PhantomData<E>);

impl<E: Pairing> IPCommitment for TIPPCommitment<E> {
    type IP = PairingInnerProduct<E>;
    type LeftKey = LeftKey<E>;
    type RightKey = RightKey<E>;
    type IPKey = PlaceholderKey;
    type LeftRightCommitment = TIPPCommOutput<E>;
    type OutputCommitment = IdentityOutput<PairingOutput<E>>;

    fn setup<'a>(size: usize, rng: impl Rng) -> Result<IPCommKey<'a, Self>, Error> {
        let srs = GenericSRS::sample(size, rng);
        Ok(srs.specialize_to_ck(size))
    }

    fn commit_with_ip<'a>(
        ck: &IPCommKey<'a, Self>,
        l: &[LeftMessage<Self>],
        r: &[RightMessage<Self>],
        ip: impl Into<Option<OutputMessage<Self>>>,
    ) -> Result<Commitment<Self>, Error> {
        let (v1, v2): (Vec<_>, Vec<_>) = ck
            .ck_a
            .iter()
            .map(|LeftKey { v_1, v_2 }| (v_1, v_2))
            .unzip();
        let (w1, w2): (Vec<_>, Vec<_>) = ck
            .ck_b
            .iter()
            .map(|RightKey { w_1, w_2 }| (w_1, w_2))
            .unzip();
        let com_l1 = multi_pairing(l, &v1).ok_or("cfg_multi_pairing failed")?;
        let com_l2 = multi_pairing(l, &v2).ok_or("cfg_multi_pairing failed")?;
        let com_r1 = multi_pairing(&w1, r).ok_or("cfg_multi_pairing failed")?;
        let com_r2 = multi_pairing(&w2, r).ok_or("cfg_multi_pairing failed")?;

        let cm_lr = TIPPCommOutput {
            t: com_l1 + com_r1,
            u: com_l2 + com_r2,
        };
        let cm_ip = ip
            .into()
            .map(|ip| {
                IdentityCommitment::<_, E::ScalarField>::commit(
                    ark_std::slice::from_ref(ck.ck_t.as_ref()),
                    &[ip],
                )
            })
            .transpose()?;

        Ok(Commitment { cm_lr, cm_ip })
    }

    fn left_key_msm<'a>(
        left_keys: &[Self::LeftKey],
        scalars: &[crate::ip_commitment::Scalar<Self>],
    ) -> Result<Self::LeftKey, Error> {
        let (v_1, v_2) = left_keys
            .iter()
            .map(|LeftKey { v_1, v_2 }| (v_1, v_2))
            .unzip::<_, _, Vec<_>, Vec<_>>();

        let v_1 = E::G2::normalize_batch(&v_1);
        let v_2 = E::G2::normalize_batch(&v_2);
        let v_1 = E::G2::msm(&v_1, scalars).unwrap();
        let v_2 = E::G2::msm(&v_2, scalars).unwrap();
        Ok(LeftKey { v_1, v_2 })
    }

    fn right_key_msm<'a>(
        ck: &[Self::RightKey],
        scalars: &[crate::ip_commitment::Scalar<Self>],
    ) -> Result<Self::RightKey, Error> {
        let (w_1, w_2) = ck
            .iter()
            .map(|RightKey { w_1, w_2 }| (w_1, w_2))
            .unzip::<_, _, Vec<_>, Vec<_>>();

        let w_1 = E::G1::normalize_batch(&w_1);
        let w_2 = E::G1::normalize_batch(&w_2);
        let w_1 = E::G1::msm(&w_1, scalars).unwrap();
        let w_2 = E::G1::msm(&w_2, scalars).unwrap();
        Ok(RightKey { w_1, w_2 })
    }
}

impl<E: Pairing> TIPPCommitment<E> {
    pub fn prove_final_ck<'a>(
        pk: &ProverKey<'a, E>,
        challenges: &[E::ScalarField],
        twist: E::ScalarField,
        kzg_point: E::ScalarField,
    ) -> FinalCommKeyProof<E> {
        let mut challenges_inv = challenges.to_vec();
        challenges_inv.push(twist);
        ark_ff::batch_inversion(&mut challenges_inv);
        let twist_inv = challenges_inv.pop().unwrap();

        let left_proof_time = start_timer!(|| "TIPP::LeftProof");
        let left_proof = kzg::prove::prove_left_key(
            &pk.h_alpha_powers,
            &pk.h_beta_powers,
            &challenges_inv,
            twist_inv,
            kzg_point,
        )
        .unwrap();
        end_timer!(left_proof_time);
        let right_proof_time = start_timer!(|| "TIPP::RightProof");
        let right_proof = kzg::prove::prove_right_key(
            &pk.g_alpha_powers,
            &pk.g_beta_powers,
            &challenges,
            kzg_point,
        )
        .unwrap();
        end_timer!(right_proof_time);
        FinalCommKeyProof {
            left_proof,
            right_proof,
        }
    }

    pub fn verify_final_ck(
        vk: &VerifierKey<E>,
        final_ck: &FinalIPCommKey<Self>,
        challenges: &[E::ScalarField],
        twist: E::ScalarField,
        kzg_point: E::ScalarField,
        proof: &FinalCommKeyProof<E>,
    ) -> bool {
        let mut challenges_inv = challenges.to_vec();
        challenges_inv.push(twist);
        ark_ff::batch_inversion(&mut challenges_inv);
        let twist_inv = challenges_inv.pop().unwrap();

        let FinalIPCommKey { ck_a, ck_b, .. } = final_ck;

        let left_is_correct = verify_left_key(
            vk,
            ck_a,
            &proof.left_proof,
            &challenges_inv,
            kzg_point,
            twist_inv,
        );
        let right_is_correct =
            verify_right_key(vk, ck_b, &proof.right_proof, &challenges, kzg_point);
        left_is_correct & right_is_correct
    }
}
