use std::{
    marker::PhantomData,
    ops::{Add, Mul, MulAssign},
};

use ark_dh_commitments::Error;
use ark_ec::pairing::{Pairing, PairingOutput};
use ark_inner_products::PairingInnerProduct;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::rand::Rng;
use derivative::Derivative;

use super::{IPCommKey, IPCommitment, LeftMessage, OutputMessage, RightMessage};

struct TIPPCommitment<E: Pairing>(PhantomData<E>);

impl<E: Pairing> IPCommitment for TIPPCommitment<E> {
    type IP = PairingInnerProduct<E>;
    type LeftKey = LeftKey<E>;
    type RightKey = RightKey<E>;
    type IPKey = InnerProductKey<E>;
    type Commitment = PairingOutput<E>;

    fn setup(size: usize, r: &mut impl Rng) -> Result<IPCommKey<'_, Self>, Error> {
        todo!()
    }

    fn commit<'a>(
        ck: &IPCommKey<'a, Self>,
        l: &[LeftMessage<Self>],
        r: &[RightMessage<Self>],
        ip: &[OutputMessage<Self>],
    ) -> Result<Self::Commitment, Error> {
        todo!()
    }
}

#[derive(Derivative, CanonicalSerialize, CanonicalDeserialize)]
#[derivative(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct LeftKey<E: Pairing> {
    v_1: E::G2,
    v_2: E::G2,
}

#[derive(Derivative, CanonicalSerialize, CanonicalDeserialize)]
#[derivative(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct RightKey<E: Pairing> {
    w_1: E::G1,
    w_2: E::G1,
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

    fn add(self, rhs: Self) -> Self::Output {
        self
    }
}
impl<E: Pairing> MulAssign<E::ScalarField> for InnerProductKey<E> {
    fn mul_assign(&mut self, rhs: E::ScalarField) {}
}

impl<E: Pairing> Mul<E::ScalarField> for InnerProductKey<E> {
    type Output = Self;
    fn mul(mut self, rhs: E::ScalarField) -> Self::Output {
        self *= rhs;
        self
    }
}
