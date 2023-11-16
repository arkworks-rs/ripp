use std::{
    convert::TryInto,
    marker::PhantomData,
    ops::{Mul},
};

use ark_ec::pairing::{Pairing, PairingOutput};
use ark_dh_commitments::Error;
use ark_std::{
    borrow::Cow,
    cfg_iter, end_timer,
    ops::{Add, MulAssign},
    rand::Rng,
    start_timer,
};
use ark_ff::UniformRand;
use ark_inner_products::{cfg_multi_pairing, InnerProduct, PairingInnerProduct};
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
    + Add<Output=Self::LeftKey>
    + Mul<Scalar<Self>, Output=Self::LeftKey>;

    type RightKey: CanonicalSerialize
    + CanonicalDeserialize
    + Clone
    + Default
    + Eq
    + Send
    + Sync
    + Add<Output=Self::RightKey>
    + Mul<Scalar<Self>, Output=Self::RightKey>;

    type IPKey: CanonicalSerialize
    + CanonicalDeserialize
    + Clone
    + Default
    + Eq
    + Send
    + Sync
    + Add<Output=Self::IPKey>
    + Mul<Scalar<Self>, Output=Self::IPKey>;

    type Commitment: CanonicalSerialize
    + CanonicalDeserialize
    + Clone
    + Default
    + Eq
    + Add<Output=Self::Commitment>
    + Mul<Scalar<Self>, Output=Self::Commitment>;

    fn setup<'a>(size: usize, r: impl Rng) -> Result<IPCommKey<'a, Self>, Error>;

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
                ck_a: self.ck_a.first().unwrap().clone(),
                ck_b: self.ck_b.first().unwrap().clone(),
                ck_t: self.ck_t.first().unwrap().clone(),
            })
        }
    }
}

impl<'a, 'b, IPC: IPCommitment> Into<IPCommKey<'a, IPC>> for &'b FinalIPCommKey<IPC> {
    fn into(self) -> IPCommKey<'a, IPC> {
        IPCommKey {
            ck_a: vec![self.ck_a.clone()].into(),
            ck_b: vec![self.ck_b.clone()].into(),
            ck_t: vec![self.ck_t.clone()].into(),
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

    pub fn to_owned<'b>(&self) -> IPCommKey<'b, IPC> {
        IPCommKey {
            ck_a: Cow::Owned(self.ck_a.to_vec()),
            ck_b: Cow::Owned(self.ck_b.to_vec()),
            ck_t: Cow::Owned(self.ck_t.to_vec()),
        }
    }

    pub fn split(&'a self, split: usize) -> (Self, Self) {
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

    pub fn fold<'b: 'a>(
        ck_1: &Self,
        ck_2: &Self,
        c_inv: &Scalar<IPC>,
        c: &Scalar<IPC>,
    ) -> Result<IPCommKey<'b, IPC>, Error> {
        // We don't do anything to ck_t when folding. These should all have the same ck_t.
        assert!(ck_1.ck_t == ck_2.ck_t);

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

        // TODO: Remove the into_owned here
        Ok(IPCommKey {
            ck_a: Cow::Owned(ck_a),
            ck_b: Cow::Owned(ck_b),
            ck_t: Cow::Owned(ck_1.ck_t.clone().into_owned()),
        })
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
        IP::Output: Default + Eq + Add<Output=IP::Output> + Mul<IP::Scalar, Output=IP::Output>;


impl<IP: InnerProduct> Add for IdentityOutput<IP>
    where
        IP::LeftMessage: Default + Eq,
        IP::RightMessage: Default + Eq,
        IP::Output: Default + Eq + Add<Output=IP::Output> + Mul<IP::Scalar, Output=IP::Output>,
{
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        IdentityOutput(
            self.0
                .iter()
                .zip(&rhs.0)
                .map(|((a, b, t), (aa, bb, tt))| (a.clone() + aa.clone(), b.clone() + bb.clone(), t.clone() + tt.clone()))
                .collect::<Vec<_>>(),
        )
    }
}

impl<IP: InnerProduct> Mul<IP::Scalar> for IdentityOutput<IP>
    where
        IP::LeftMessage: Default + Eq,
        IP::RightMessage: Default + Eq,
        IP::Output: Default + Eq + Add<Output=IP::Output> + Mul<IP::Scalar, Output=IP::Output>,
{
    type Output = Self;

    fn mul(self, rhs: IP::Scalar) -> Self::Output {
        IdentityOutput(
            self.0.into_iter().map(|(a, b, t)| {
                (a * rhs,
                 b * rhs,
                 t * rhs)
            }).collect())
    }
}

impl<IP: InnerProduct> IPCommitment for IdentityCommitment<IP>
    where
        IP::LeftMessage: Default + Eq,
        IP::RightMessage: Default + Eq,
        IP::Output: Default + Eq + Add<Output=IP::Output> + Mul<IP::Scalar, Output=IP::Output>,
{
    type IP = IP;

    type LeftKey = HomomorphicPlaceholderValue;
    type RightKey = HomomorphicPlaceholderValue;
    type IPKey = HomomorphicPlaceholderValue;

    type Commitment = IdentityOutput<IP>;

    fn setup<'a>(size: usize, _r: impl Rng) -> Result<IPCommKey<'a, Self>, Error> {
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
            l.iter().zip(r.iter()).zip(ip.iter()).map(|((a, b), c)| (a.clone(), b.clone(), c.clone())).collect(),
        ))
    }
}

/// Represents a commitment scheme over (G1, G2) which does a pairing commitment on the LHS and a
/// pairing commitment on the RHS, and an identity commitment on the inner product of the two.
pub(crate) struct PairingCommitment<E: Pairing> {
    _marker: PhantomData<E>,
}

#[derive(Clone, CanonicalSerialize, CanonicalDeserialize, PartialEq, Eq)]
pub struct PairingCommOutput<E: Pairing> {
    com_a: PairingOutput<E>,
    com_b: PairingOutput<E>,
    com_t: Vec<PairingOutput<E>>,
}

impl<E: Pairing> Default for PairingCommOutput<E> {
    fn default() -> Self {
        PairingCommOutput {
            com_a: PairingOutput::default(),
            com_b: PairingOutput::default(),
            com_t: vec![PairingOutput::default()],
        }
    }
}

impl<E: Pairing> Add for PairingCommOutput<E> {
    type Output = Self;
    fn add(self, rhs: Self) -> Self::Output {
        PairingCommOutput {
            com_a: self.com_a + rhs.com_a,
            com_b: self.com_b + rhs.com_b,
            com_t: self.com_t.iter().zip(rhs.com_t.iter()).map(|(l, r)| l + r).collect(),
        }
    }
}

impl<E: Pairing> Mul<E::ScalarField> for PairingCommOutput<E> {
    type Output = Self;

    fn mul(self, rhs: E::ScalarField) -> Self::Output {
        PairingCommOutput {
            com_a: self.com_a * rhs,
            com_b: self.com_b * rhs,
            com_t: self.com_t.iter().map(|c| *c * rhs).collect(),
        }
    }
}

impl<E: Pairing> IPCommitment for PairingCommitment<E>
{
    type IP = PairingInnerProduct<E>;

    type LeftKey = E::G2;
    type RightKey = E::G1;
    type IPKey = HomomorphicPlaceholderValue;

    type Commitment = PairingCommOutput<E>;

    fn setup<'a>(size: usize, mut rng: impl Rng) -> Result<IPCommKey<'a, Self>, Error> {
        let random_left_key: Vec<E::G2> = (0..size).map(|_| E::G2::rand(&mut rng)).collect();
        let random_right_key: Vec<E::G1> = (0..size).map(|_| E::G1::rand(&mut rng)).collect();

        Ok(IPCommKey {
            ck_a: random_left_key.into(),
            ck_b: random_right_key.into(),
            ck_t: vec![HomomorphicPlaceholderValue].into(),
        })
    }

    fn commit<'a>(
        ck: &IPCommKey<'a, Self>,
        l: &[LeftMessage<Self>],
        r: &[RightMessage<Self>],
        ip: &[OutputMessage<Self>],
    ) -> Result<Self::Commitment, Error> {
        let com_a = cfg_multi_pairing(l, &ck.ck_a).ok_or("invalid pairing")?;
        let com_b = cfg_multi_pairing(&ck.ck_b, r).ok_or("invalid pairing")?;
        let com_t = ip.to_vec();

        Ok(PairingCommOutput {
            com_a,
            com_b,
            com_t,
        })
    }
}