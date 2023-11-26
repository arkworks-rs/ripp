use ark_dh_commitments::{identity::PlaceholderKey, Error};
use ark_ec::pairing::{Pairing, PairingOutput};
use ark_ff::Field;
use ark_inner_products::{multi_pairing, PairingInnerProduct};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::{
    marker::PhantomData,
    ops::{Add, Mul, MulAssign},
    rand::Rng,
    One,
};
use derivative::Derivative;

use crate::ip_commitment::{IPCommKey, IPCommitment, LeftMessage, OutputMessage, RightMessage};

use super::data_structures::GenericSRS;

pub(crate) struct TIPPCommitment<E: Pairing>(PhantomData<E>);

#[derive(Clone, Debug, Copy, PartialEq, Eq, CanonicalSerialize, CanonicalDeserialize)]
pub(crate) struct TIPPCommOutput<E: Pairing> {
    t: PairingOutput<E>,
    u: PairingOutput<E>,
}

impl<E: Pairing> Add for TIPPCommOutput<E> {
    type Output = Self;
    fn add(mut self, rhs: Self) -> Self::Output {
        self.t += rhs.t;
        self.u += rhs.u;
        self
    }
}

impl<E: Pairing> Mul<E::ScalarField> for TIPPCommOutput<E> {
    type Output = Self;
    fn mul(mut self, rhs: E::ScalarField) -> Self::Output {
        self.t *= rhs;
        self.u *= rhs;
        self
    }
}

impl<E: Pairing> Default for TIPPCommOutput<E> {
    fn default() -> Self {
        let one = E::TargetField::one();
        Self {
            t: PairingOutput(one),
            u: PairingOutput(one),
        }
    }
}

impl<E: Pairing> IPCommitment for TIPPCommitment<E> {
    type IP = PairingInnerProduct<E>;
    type LeftKey = LeftKey<E>;
    type RightKey = RightKey<E>;
    type IPKey = PlaceholderKey;
    type Commitment = TIPPCommOutput<E>;

    fn setup<'a>(size: usize, rng: impl Rng) -> Result<IPCommKey<'a, Self>, Error> {
        let srs = GenericSRS::sample(size, rng);
        let (pk, vk) = srs.specialize(size);
        Ok(pk.ck)
    }

    fn commit<'a>(
        ck: &IPCommKey<'a, Self>,
        l: &[LeftMessage<Self>],
        r: &[RightMessage<Self>],
        _ip: impl Fn() -> OutputMessage<Self>,
    ) -> Result<Self::Commitment, Error> {
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

        Ok(TIPPCommOutput {
            t: com_l1 + com_r1,
            u: com_l2 + com_r2,
        })
    }
}

pub fn structured_scalar_power<F: Field>(num: usize, s: &F) -> Vec<F> {
    let mut powers = vec![F::one()];
    for i in 1..num {
        powers.push(powers[i - 1] * s);
    }
    powers
}

#[derive(Derivative, CanonicalSerialize, CanonicalDeserialize)]
#[derivative(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct LeftKey<E: Pairing> {
    pub v_1: E::G2,
    pub v_2: E::G2,
}

#[derive(Derivative, CanonicalSerialize, CanonicalDeserialize)]
#[derivative(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct RightKey<E: Pairing> {
    pub w_1: E::G1,
    pub w_2: E::G1,
}

#[derive(Derivative, CanonicalSerialize, CanonicalDeserialize)]
#[derivative(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct InnerProductKey<E: Pairing>(PhantomData<E>);

/************************************************/
/************************************************/
/************************************************/

impl<E: Pairing> Add<Self> for LeftKey<E> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self {
            v_1: self.v_1 + rhs.v_1,
            v_2: self.v_2 + rhs.v_2,
        }
    }
}

impl<E: Pairing> MulAssign<E::ScalarField> for LeftKey<E> {
    fn mul_assign(&mut self, rhs: E::ScalarField) {
        self.v_1 *= rhs;
        self.v_2 *= rhs;
    }
}

impl<E: Pairing> Mul<E::ScalarField> for LeftKey<E> {
    type Output = Self;
    fn mul(mut self, rhs: E::ScalarField) -> Self::Output {
        self *= rhs;
        self
    }
}

impl<E: Pairing> Add<Self> for RightKey<E> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self {
            w_1: self.w_1 + rhs.w_1,
            w_2: self.w_2 + rhs.w_2,
        }
    }
}

impl<E: Pairing> MulAssign<E::ScalarField> for RightKey<E> {
    fn mul_assign(&mut self, rhs: E::ScalarField) {
        self.w_1 *= rhs;
        self.w_2 *= rhs;
    }
}

impl<E: Pairing> Mul<E::ScalarField> for RightKey<E> {
    type Output = Self;
    fn mul(mut self, rhs: E::ScalarField) -> Self::Output {
        self *= rhs;
        self
    }
}

impl<E: Pairing> Add<Self> for InnerProductKey<E> {
    type Output = Self;

    fn add(self, _rhs: Self) -> Self::Output {
        self
    }
}

impl<E: Pairing> MulAssign<E::ScalarField> for InnerProductKey<E> {
    fn mul_assign(&mut self, _rhs: E::ScalarField) {}
}

impl<E: Pairing> Mul<E::ScalarField> for InnerProductKey<E> {
    type Output = Self;
    fn mul(mut self, rhs: E::ScalarField) -> Self::Output {
        self *= rhs;
        self
    }
}
