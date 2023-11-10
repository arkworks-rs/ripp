use std::{
    convert::TryInto,
    marker::PhantomData,
    ops::{AddAssign, Mul},
};

use ark_dh_commitments::Error;
use ark_std::{
    borrow::Cow,
    cfg_iter, end_timer,
    ops::{Add, MulAssign},
    rand::Rng,
    start_timer,
};

use ark_inner_products::InnerProduct;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};

use derivative::Derivative;
#[cfg(feature = "parallel")]
use rayon::prelude::*;

use crate::mul_helper;

mod mipp;
mod snarkpack;
mod ssm;
mod tipp;

pub type LeftMessage<IPC> = <<IPC as IPCommitment>::IP as InnerProduct>::LeftMessage;
pub type RightMessage<IPC> = <<IPC as IPCommitment>::IP as InnerProduct>::RightMessage;
pub type OutputMessage<IPC> = <<IPC as IPCommitment>::IP as InnerProduct>::Output;
pub type Scalar<IPC> = <<IPC as IPCommitment>::IP as InnerProduct>::Scalar;

pub trait IPCommitment: Sized {
    type IP: InnerProduct;

    type LeftKey: CanonicalSerialize
        + CanonicalDeserialize
        + Clone
        + Default
        + Eq
        + Send
        + Sync
        + Add<Self::LeftKey, Output = Self::LeftKey>
        + Mul<Scalar<Self>, Output = Self::LeftKey>
        + MulAssign<Scalar<Self>>;

    type RightKey: CanonicalSerialize
        + CanonicalDeserialize
        + Clone
        + Default
        + Eq
        + Send
        + Sync
        + Add<Self::RightKey, Output = Self::RightKey>
        + Mul<Scalar<Self>, Output = Self::RightKey>
        + MulAssign<Scalar<Self>>;

    type IPKey: CanonicalSerialize
        + CanonicalDeserialize
        + Clone
        + Default
        + Eq
        + Send
        + Sync
        + Add<Self::IPKey, Output = Self::IPKey>
        + Mul<Scalar<Self>, Output = Self::IPKey>
        + MulAssign<Scalar<Self>>;

    type Commitment: CanonicalSerialize
        + CanonicalDeserialize
        + Default
        + Eq
        + Add<Self::Commitment, Output = Self::Commitment>
        + AddAssign<Self::Commitment>
        + Mul<Scalar<Self>, Output = Self::Commitment>
        + MulAssign<Scalar<Self>>;

    fn setup(size: usize, r: &mut impl Rng) -> Result<IPCommKey<'_, Self>, Error>;

    fn commit<'a>(
        ck: &IPCommKey<'a, Self>,
        l: &[LeftMessage<Self>],
        r: &[RightMessage<Self>],
        ip: &[OutputMessage<Self>],
    ) -> Result<Self::Commitment, Error>;

    fn verify<'a>(
        ck: &IPCommKey<'a, Self>,
        l: &[LeftMessage<Self>],
        r: &[RightMessage<Self>],
        ip: &[OutputMessage<Self>],
        com: &Self::Commitment,
    ) -> Result<bool, Error> {
        Ok(Self::commit(ck, l, r, ip)? == *com)
    }
}

/// A final IP comm key is an IP comm key of length 1
#[derive(Derivative)]
#[derivative(Clone(bound = ""))]
#[derive(CanonicalSerialize, CanonicalDeserialize)]
pub struct FinalIPCommKey<IPC: IPCommitment> {
    pub ck_a: IPC::LeftKey,
    pub ck_b: IPC::RightKey,
    pub ck_t: IPC::IPKey,
}

impl<'a, IPC: IPCommitment> TryInto<FinalIPCommKey<IPC>> for IPCommKey<'a, IPC> {
    type Error = ();

    fn try_into(mut self) -> Result<FinalIPCommKey<IPC>, ()> {
        if self.ck_a.len() != 1 || self.ck_b.len() != 1 || self.ck_t.len() != 1 {
            Err(())
        } else {
            Ok(FinalIPCommKey {
                ck_a: self.ck_a.pop().unwrap(),
                ck_b: self.ck_b.pop().unwrap(),
                ck_t: self.ck_t.pop().unwrap(),
            })
        }
    }
}

impl<'a, 'b, IPC: IPCommitment> Into<IPCommKey<'a, IPC>> for &'b FinalIPCommKey<IPC> {
    fn into(self) -> IPCommKey<'a, IPC> {
        IPCommKey {
            ck_a: vec![self.ck_a].into(),
            ck_b: vec![self.ck_b].into(),
            ck_t: vec![self.ck_t].into(),
        }
    }
}

#[derive(Derivative)]
#[derivative(Clone(bound = ""))]
pub struct IPCommKey<'a, IPC: IPCommitment> {
    pub ck_a: Cow<'a, [IPC::LeftKey]>,
    pub ck_b: Cow<'a, [IPC::RightKey]>,
    pub ck_t: Cow<'a, [IPC::IPKey]>,
}

impl<'a, IPC: IPCommitment> IPCommKey<'a, IPC> {
    pub fn new(
        ck_a: Cow<'a, [IPC::LeftKey]>,
        ck_b: Cow<'a, [IPC::RightKey]>,
        ck_t: Cow<'a, [IPC::IPKey]>,
    ) -> Self {
        Self { ck_a, ck_b, ck_t }
    }

    pub fn split(&self, split: usize) -> (Self, Self) {
        let ck_a_1 = &self.ck_a[..split];
        let ck_a_2 = &self.ck_a[split..];

        let ck_b_1 = &self.ck_b[split..];
        let ck_b_2 = &self.ck_b[..split];
        let ck_1 = IPCommKey::new(
            Cow::Borrowed(ck_a_1),
            Cow::Borrowed(ck_b_1),
            Cow::Borrowed(&self.ck_t[..]),
        );

        let ck_2 = IPCommKey::new(
            Cow::Borrowed(ck_a_2),
            Cow::Borrowed(ck_b_2),
            Cow::Borrowed(&self.ck_t[..]),
        );
        (ck_1, ck_2)
    }

    pub fn fold_into(
        &mut self,
        ck_1: &Self,
        ck_2: &Self,
        c_inv: &Scalar<IPC>,
        c: &Scalar<IPC>,
    ) -> Result<(), Error> {
        let rescale_a = start_timer!(|| "Rescale CK_B");
        let ck_a = cfg_iter!(ck_2.ck_a.as_ref())
            .map(|a| mul_helper(a, &c_inv))
            .zip(ck_1.ck_a.as_ref())
            .map(|(a_1, a_2)| a_1 + a_2.clone())
            .collect::<Vec<IPC::LeftKey>>();
        end_timer!(rescale_a);

        let rescale_b = start_timer!(|| "Rescale CK_A");
        let ck_b = cfg_iter!(ck_1.ck_b.as_ref())
            .map(|b| mul_helper(b, &c))
            .zip(ck_2.ck_b.as_ref())
            .map(|(b_1, b_2)| b_1 + b_2.clone())
            .collect::<Vec<IPC::RightKey>>();
        end_timer!(rescale_b);
        self.ck_a = Cow::Owned(ck_a);
        self.ck_b = Cow::Owned(ck_b);
        Ok(())
    }
}

#[derive(Copy, Clone)]
pub struct IdentityCommitment<IP: InnerProduct> {
    _ip: PhantomData<IP>,
}

#[derive(CanonicalSerialize, CanonicalDeserialize, Clone, Default, Eq, PartialEq)]
pub struct HomomorphicPlaceholderValue;

impl Add for HomomorphicPlaceholderValue {
    type Output = Self;

    fn add(self, _rhs: Self) -> Self::Output {
        HomomorphicPlaceholderValue {}
    }
}

impl<T> Mul<T> for HomomorphicPlaceholderValue {
    type Output = HomomorphicPlaceholderValue;

    fn mul(self, _other: T) -> Self {
        self
    }
}

impl<T> MulAssign<T> for HomomorphicPlaceholderValue {
    fn mul_assign(&mut self, _rhs: T) {}
}

#[derive(Derivative)]
#[derivative(Clone(bound = ""))]
#[derivative(Default(bound = ""))]
#[derivative(PartialEq(bound = ""))]
#[derivative(Eq(bound = ""))]
#[derive(CanonicalSerialize, CanonicalDeserialize)]
pub struct IdentityOutput<IP: InnerProduct>(
    pub Vec<(IP::LeftMessage, IP::RightMessage, IP::Output)>,
)
where
    IP::LeftMessage: Default + Eq,
    IP::RightMessage: Default + Eq,
    IP::Output: Default + Eq;

impl<IP: InnerProduct> Add for IdentityOutput<IP>
where
    IP::LeftMessage: Default + Eq,
    IP::RightMessage: Default + Eq,
    IP::Output: Default + Eq,
{
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        IdentityOutput(
            self.0
                .iter()
                .zip(&rhs.0)
                .map(|((a, b, t), (aa, bb, tt))| (a + aa, b + bb, t + tt))
                .collect::<Vec<_>>(),
        )
    }
}

impl<IP: InnerProduct> AddAssign for IdentityOutput<IP>
where
    IP::LeftMessage: Default + Eq,
    IP::RightMessage: Default + Eq,
    IP::Output: Default + Eq,
{
    fn add_assign(&mut self, other: Self) {
        *self = self + other
    }
}

impl<IP: InnerProduct> MulAssign<IP::Scalar> for IdentityOutput<IP>
where
    IP::LeftMessage: Default + Eq,
    IP::RightMessage: Default + Eq,
    IP::Output: Default + Eq,
{
    fn mul_assign(&mut self, rhs: IP::Scalar) {
        self.0.iter_mut().for_each(|(a, b, t)| {
            a.mul_assign(rhs.clone());
            b.mul_assign(rhs.clone());
            t.mul_assign(rhs.clone());
        })
    }
}

impl<IP: InnerProduct> IPCommitment for IdentityCommitment<IP>
where
    IP::LeftMessage: Default + Eq,
    IP::RightMessage: Default + Eq,
    IP::Output: Default + Eq,
{
    type IP = IP;

    type LeftKey = HomomorphicPlaceholderValue;
    type RightKey = HomomorphicPlaceholderValue;
    type IPKey = HomomorphicPlaceholderValue;

    type Commitment = IdentityOutput<IP>;

    fn setup(size: usize, _r: &mut impl Rng) -> Result<IPCommKey<'_, Self>, Error> {
        Ok(IPCommKey {
            ck_a: vec![HomomorphicPlaceholderValue; size].into(),
            ck_b: vec![HomomorphicPlaceholderValue; size].into(),
            ck_t: vec![HomomorphicPlaceholderValue; size].into(),
        })
    }

    fn commit<'a>(
        _ck: &IPCommKey<'a, Self>,
        l: &[LeftMessage<Self>],
        r: &[RightMessage<Self>],
        ip: &[OutputMessage<Self>],
    ) -> Result<Self::Commitment, Error> {
        Ok(IdentityOutput(
            l.iter().zip(r.iter()).zip(ip.iter()).collect(),
        ))
    }
}
