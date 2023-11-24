use std::{
    marker::PhantomData,
    ops::{Add, Mul, MulAssign},
};

use ark_inner_products::InnerProduct;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::rand::Rng;
use derivative::Derivative;

use crate::Error;

use super::{IPCommKey, IPCommitment, LeftMessage, OutputMessage, RightMessage};

#[derive(Copy, Clone)]
pub struct IdentityCommitment<IP: InnerProduct> {
    _ip: PhantomData<IP>,
}

#[derive(CanonicalSerialize, CanonicalDeserialize, Debug, Clone, Default, Eq, PartialEq)]
pub struct PlaceholderKey;

impl Add for PlaceholderKey {
    type Output = Self;

    fn add(self, _rhs: Self) -> Self::Output {
        PlaceholderKey {}
    }
}

impl<T> Mul<T> for PlaceholderKey {
    type Output = PlaceholderKey;

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
pub struct IdentityOutput<IP: InnerProduct>
where
    IP::Output: Default + Eq + Add<Output = IP::Output> + Mul<IP::Scalar, Output = IP::Output>,
{
    left_msg: Vec<IP::LeftMessage>,
    right_msg: Vec<IP::RightMessage>,
    out: Vec<IP::Output>,
}

impl<IP: InnerProduct> Add for IdentityOutput<IP>
where
    IP::Output: Default + Eq + Add<Output = IP::Output> + Mul<IP::Scalar, Output = IP::Output>,
{
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        IdentityOutput {
            left_msg: self
                .left_msg
                .into_iter()
                .zip(rhs.left_msg.into_iter())
                .map(|(l, r)| l + r)
                .collect(),
            right_msg: self
                .right_msg
                .into_iter()
                .zip(rhs.right_msg.into_iter())
                .map(|(l, r)| l + r)
                .collect(),
            out: self
                .out
                .into_iter()
                .zip(rhs.out.into_iter())
                .map(|(l, r)| l + r)
                .collect(),
        }
    }
}

impl<IP: InnerProduct> Mul<IP::Scalar> for IdentityOutput<IP>
where
    IP::Output: Default + Eq + Add<Output = IP::Output> + Mul<IP::Scalar, Output = IP::Output>,
{
    type Output = Self;

    fn mul(self, rhs: IP::Scalar) -> Self::Output {
        IdentityOutput {
            left_msg: self.left_msg.into_iter().map(|l| l * rhs).collect(),
            right_msg: self.right_msg.into_iter().map(|l| l * rhs).collect(),
            out: self.out.into_iter().map(|l| l * rhs).collect(),
        }
    }
}

impl<IP: InnerProduct> IPCommitment for IdentityCommitment<IP>
where
    IP::Output: Default + Eq + Add<Output = IP::Output> + Mul<IP::Scalar, Output = IP::Output>,
{
    type IP = IP;

    type LeftKey = PlaceholderKey;
    type RightKey = PlaceholderKey;
    type IPKey = PlaceholderKey;

    type Commitment = IdentityOutput<IP>;

    fn setup<'a>(size: usize, _r: impl Rng) -> Result<IPCommKey<'a, Self>, Error> {
        Ok(IPCommKey {
            ck_a: vec![PlaceholderKey; size].into(),
            ck_b: vec![PlaceholderKey; size].into(),
            ck_t: vec![PlaceholderKey].into(),
        })
    }

    fn commit<'a>(
        _ck: &IPCommKey<'a, Self>,
        l: &[LeftMessage<Self>],
        r: &[RightMessage<Self>],
        ip: &[OutputMessage<Self>],
    ) -> Result<Self::Commitment, Error> {
        Ok(IdentityOutput {
            left_msg: l.to_vec(),
            right_msg: r.to_vec(),
            out: ip.to_vec(),
        })
    }
}
