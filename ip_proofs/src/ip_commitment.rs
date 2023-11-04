use std::borrow::Cow;

use ark_dh_commitments::Error;
use ark_std::{
    cfg_iter, end_timer,
    ops::{Add, MulAssign},
    rand::Rng,
    start_timer,
};

use ark_inner_products::InnerProduct;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};

use crate::mul_helper;

mod mipp;
mod snarkpack;
mod ssm;
mod tipp;

type LeftMessage<IPC> = <<IPC as IPCommitment>::IP as InnerProduct>::LeftMessage;
type RightMessage<IPC> = <<IPC as IPCommitment>::IP as InnerProduct>::RightMessage;
type OutputMessage<IPC> = <<IPC as IPCommitment>::IP as InnerProduct>::Output;

pub trait IPCommitment {
    type IP: InnerProduct;
    type Scalar;
    type LeftKey: CanonicalSerialize
        + CanonicalDeserialize
        + Clone
        + Default
        + Eq
        + Send
        + Sync
        + Add<Self::LeftKey, Output = Self::LeftKey>
        + MulAssign<Self::Scalar>;
    type RightKey: CanonicalSerialize
        + CanonicalDeserialize
        + Clone
        + Default
        + Eq
        + Send
        + Sync
        + Add<Self::RightKey, Output = Self::RightKey>
        + MulAssign<Self::Scalar>;
    type IPKey: CanonicalSerialize
        + CanonicalDeserialize
        + Clone
        + Default
        + Eq
        + Send
        + Sync
        + Add<Self::IPKey, Output = Self::IPKey>
        + MulAssign<Self::Scalar>;

    type Output: CanonicalSerialize
        + CanonicalDeserialize
        + Clone
        + Default
        + Eq
        + Add<Self::Output, Output = Self::Output>
        + MulAssign<Self::Scalar>;

    fn setup(size: usize, r: &mut impl Rng) -> Result<IPCommKey<'_, Self>, Error>;

    fn commit<'a>(
        ck: &IPCommKey<'a, Self>,
        l: &[LeftMessage<Self>],
        r: &[RightMessage<Self>],
        ip: &[OutputMessage<Self>],
    ) -> Result<Self::Output, Error>;

    fn verify<'a>(
        ck: &IPCommKey<'a, Self>,
        l: &[LeftMessage<Self>],
        r: &[RightMessage<Self>],
        ip: &[OutputMessage<Self>],
        com: &Self::Output,
    ) -> Result<bool, Error> {
        Ok(Self::commit(ck, l, r, ip)? == *com)
    }
}

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
        c_inv: &IPC::Scalar,
        c: &IPC::Scalar,
    ) -> Result<(), Error> {
        let rescale_a = start_timer!(|| "Rescale CK_B");
        let ck_a = cfg_iter!(&ck_2.ck_a)
            .map(|a| mul_helper(a, &c_inv))
            .zip(&ck_1.ck_a)
            .map(|(a_1, a_2)| a_1 + a_2.clone())
            .collect::<Vec<IPC::LeftKey>>();
        end_timer!(rescale_a);

        let rescale_b = start_timer!(|| "Rescale CK_A");
        let ck_b = cfg_iter!(&ck_1.ck_b)
            .map(|b| mul_helper(b, &c))
            .zip(ck_2.ck_b)
            .map(|(b_1, b_2)| b_1 + b_2.clone())
            .collect::<Vec<IPC::RightKey>>();
        end_timer!(rescale_b);
        *self.ck_a = Cow::Owned(ck_a);
        *self.ck_b = Cow::Owned(ck_b);
        Ok(())
    }
}
