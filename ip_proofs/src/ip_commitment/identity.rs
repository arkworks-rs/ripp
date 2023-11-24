use ark_std::{
    borrow::Cow,
    cfg_into_iter,
    marker::PhantomData,
    ops::{Add, Mul, MulAssign},
    rand::Rng,
};

use ark_inner_products::InnerProduct;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use derivative::Derivative;

#[cfg(feature = "parallel")]
use rayon::prelude::*;

use crate::Error;

use super::{IPCommKey, IPCommitment, LeftMessage, OutputMessage, RightMessage};

#[derive(Copy, Clone)]
pub struct IdentityCommitment<IP: InnerProduct> {
    _ip: PhantomData<IP>,
}

#[derive(CanonicalSerialize, CanonicalDeserialize, Debug, Clone, Copy, Default, Eq, PartialEq)]
pub struct PlaceholderKey;

impl Add for PlaceholderKey {
    type Output = Self;

    fn add(self, _rhs: Self) -> Self {
        self
    }
}

impl<T> Mul<T> for PlaceholderKey {
    type Output = Self;

    fn mul(self, _other: T) -> Self {
        self
    }
}

impl<T> MulAssign<T> for PlaceholderKey {
    fn mul_assign(&mut self, _rhs: T) {}
}

#[derive(Derivative)]
#[derivative(Clone(bound = ""))]
#[derivative(Default(bound = ""))]
#[derivative(PartialEq(bound = ""))]
#[derivative(Eq(bound = ""))]
#[derive(CanonicalSerialize, CanonicalDeserialize)]
pub struct IdentityComm<IP: InnerProduct> {
    left_msg: Vec<IP::LeftMessage>,
    right_msg: Vec<IP::RightMessage>,
    out: IP::Output,
}

impl<IP: InnerProduct> Add for IdentityComm<IP> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        let left_msg = cfg_into_iter!(self.left_msg)
            .zip(rhs.left_msg)
            .map(|(l, r)| l + r)
            .collect();
        let right_msg = cfg_into_iter!(self.right_msg)
            .zip(rhs.right_msg)
            .map(|(l, r)| l + r)
            .collect();
        let out = self.out + rhs.out;

        IdentityComm {
            left_msg,
            right_msg,
            out,
        }
    }
}

impl<IP: InnerProduct> Mul<IP::Scalar> for IdentityComm<IP> {
    type Output = Self;

    fn mul(self, rhs: IP::Scalar) -> Self::Output {
        let left_msg = cfg_into_iter!(self.left_msg).map(|l| l * rhs).collect();
        let right_msg = cfg_into_iter!(self.right_msg).map(|l| l * rhs).collect();
        let out = self.out * rhs;
        IdentityComm {
            left_msg,
            right_msg,
            out,
        }
    }
}

impl<IP: InnerProduct> IPCommitment for IdentityCommitment<IP> {
    type IP = IP;

    type LeftKey = PlaceholderKey;
    type RightKey = PlaceholderKey;
    type IPKey = PlaceholderKey;

    type Commitment = IdentityComm<IP>;

    fn setup<'a>(size: usize, _r: impl Rng) -> Result<IPCommKey<'a, Self>, Error> {
        Ok(IPCommKey {
            ck_a: vec![PlaceholderKey; size].into(),
            ck_b: vec![PlaceholderKey; size].into(),
            ck_t: Cow::Owned(PlaceholderKey),
        })
    }

    fn commit<'a>(
        _ck: &IPCommKey<'a, Self>,
        l: &[LeftMessage<Self>],
        r: &[RightMessage<Self>],
        ip: impl Fn() -> OutputMessage<Self>,
    ) -> Result<Self::Commitment, Error> {
        Ok(IdentityComm {
            left_msg: l.to_vec(),
            right_msg: r.to_vec(),
            out: ip(),
        })
    }
}
