use ark_dh_commitments::{
    identity::{IdentityCommitment, IdentityOutput, PlaceholderKey},
    DoublyHomomorphicCommitment, Error,
};
use ark_ec::{
    pairing::{Pairing, PairingOutput},
    CurveGroup, VariableBaseMSM,
};
use ark_ff::Field;
use ark_inner_products::{multi_pairing, PairingInnerProduct};
use ark_std::{marker::PhantomData, rand::Rng};

use crate::ip_commitment::{
    Commitment, IPCommKey, IPCommitment, LeftMessage, OutputMessage, RightMessage,
};

mod data_structures;
pub use data_structures::*;

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

pub fn structured_scalar_power<F: Field>(num: usize, s: &F) -> Vec<F> {
    let mut powers = vec![F::one()];
    for i in 1..num {
        powers.push(powers[i - 1] * s);
    }
    powers
}
