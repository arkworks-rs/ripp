use ark_dh_commitments::Error;
use ark_inner_products::InnerProduct;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::{
    borrow::Cow,
    cfg_iter,
    convert::TryInto,
    end_timer,
    fmt::Debug,
    ops::{Add, Mul},
    rand::Rng,
    start_timer,
};

use derivative::Derivative;
#[cfg(feature = "parallel")]
use rayon::prelude::*;

pub mod generic;

pub mod mipp;
pub mod pairing;
pub mod scalar;

pub type LeftMessage<IPC> = <<IPC as IPCommitment>::IP as InnerProduct>::LeftMessage;
pub type RightMessage<IPC> = <<IPC as IPCommitment>::IP as InnerProduct>::RightMessage;
pub type OutputMessage<IPC> = <<IPC as IPCommitment>::IP as InnerProduct>::Output;
pub type Scalar<IPC> = <<IPC as IPCommitment>::IP as InnerProduct>::Scalar;

pub trait IPCommitment: Sized {
    type IP: InnerProduct;

    type LeftKey: CanonicalSerialize
        + CanonicalDeserialize
        + Copy
        + Default
        + Debug
        + Eq
        + Send
        + Sync
        + Add<Output = Self::LeftKey>
        + Mul<Scalar<Self>, Output = Self::LeftKey>;

    type RightKey: CanonicalSerialize
        + CanonicalDeserialize
        + Copy
        + Default
        + Debug
        + Eq
        + Send
        + Sync
        + Add<Output = Self::RightKey>
        + Mul<Scalar<Self>, Output = Self::RightKey>;

    type IPKey: CanonicalSerialize
        + CanonicalDeserialize
        + Copy
        + Default
        + Debug
        + Eq
        + Send
        + Sync
        + Add<Output = Self::IPKey>
        + Mul<Scalar<Self>, Output = Self::IPKey>;

    type Commitment: CanonicalSerialize
        + CanonicalDeserialize
        + Clone
        + Debug
        + Default
        + Eq
        + Add<Output = Self::Commitment>
        + Mul<Scalar<Self>, Output = Self::Commitment>;

    fn setup<'a>(size: usize, r: impl Rng) -> Result<IPCommKey<'a, Self>, Error>;

    fn commit<'a>(
        ck: &IPCommKey<'a, Self>,
        l: &[LeftMessage<Self>],
        r: &[RightMessage<Self>],
        ip: impl Fn() -> OutputMessage<Self>,
    ) -> Result<Self::Commitment, Error>;

    fn verify<'a>(
        ck: &IPCommKey<'a, Self>,
        l: &[LeftMessage<Self>],
        r: &[RightMessage<Self>],
        ip: &OutputMessage<Self>,
        com: &Self::Commitment,
    ) -> Result<bool, Error> {
        Ok(dbg!(Self::commit(ck, l, r, || ip.clone())?) == *com)
    }

    fn left_key_msm<'a>(
        ck: &[Self::LeftKey],
        scalars: &[Scalar<Self>],
    ) -> Result<Self::LeftKey, Error> {
        #[cfg(feature = "parallel")]
        let zero = Self::LeftKey::default;

        #[cfg(not(feature = "parallel"))]
        let zero = Self::LeftKey::default();

        Ok(cfg_iter!(ck)
            .zip(scalars)
            .map(|(a, s)| *a * *s)
            .reduce(zero, |a, b| a + b))
    }

    fn right_key_msm<'a>(
        ck: &[Self::RightKey],
        scalars: &[Scalar<Self>],
    ) -> Result<Self::RightKey, Error> {
        #[cfg(feature = "parallel")]
        let zero = Self::RightKey::default;

        #[cfg(not(feature = "parallel"))]
        let zero = Self::RightKey::default();

        Ok(cfg_iter!(ck)
            .zip(scalars)
            .map(|(a, s)| *a * *s)
            .reduce(zero, |a, b| a + b))
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

    fn try_into(self) -> Result<FinalIPCommKey<IPC>, ()> {
        if self.ck_a.len() != 1 || self.ck_b.len() != 1 {
            Err(())
        } else {
            Ok(FinalIPCommKey {
                ck_a: self.ck_a.first().unwrap().clone(),
                ck_b: self.ck_b.first().unwrap().clone(),
                ck_t: self.ck_t.into_owned(),
            })
        }
    }
}

impl<'a, 'b, IPC: IPCommitment> Into<IPCommKey<'a, IPC>> for &'b FinalIPCommKey<IPC> {
    fn into(self) -> IPCommKey<'a, IPC> {
        IPCommKey {
            ck_a: vec![self.ck_a.clone()].into(),
            ck_b: vec![self.ck_b.clone()].into(),
            ck_t: Cow::Owned(self.ck_t.clone()),
        }
    }
}

#[derive(Derivative)]
#[derivative(Clone(bound = ""), Debug(bound = "IPC: IPCommitment"))]
pub struct IPCommKey<'a, IPC: IPCommitment> {
    pub ck_a: Cow<'a, [IPC::LeftKey]>,
    pub ck_b: Cow<'a, [IPC::RightKey]>,
    pub ck_t: Cow<'a, IPC::IPKey>,
}

impl<'a, IPC: IPCommitment> IPCommKey<'a, IPC> {
    pub fn new(
        ck_a: Cow<'a, [IPC::LeftKey]>,
        ck_b: Cow<'a, [IPC::RightKey]>,
        ck_t: Cow<'a, IPC::IPKey>,
    ) -> Self {
        Self { ck_a, ck_b, ck_t }
    }

    pub fn to_owned<'b>(&self) -> IPCommKey<'b, IPC> {
        IPCommKey {
            ck_a: Cow::Owned(self.ck_a.to_vec()),
            ck_b: Cow::Owned(self.ck_b.to_vec()),
            ck_t: Cow::Owned(self.ck_t.clone().into_owned()),
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
            Cow::Borrowed(self.ck_t.as_ref()),
        );

        let ck_2 = IPCommKey::new(
            Cow::Borrowed(ck_a_2),
            Cow::Borrowed(ck_b_2),
            Cow::Borrowed(self.ck_t.as_ref()),
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
        assert_eq!(ck_1.ck_t, ck_2.ck_t);

        let rescale_a = start_timer!(|| "Rescale CK_B");
        let ck_a = cfg_iter!(ck_2.ck_a.as_ref())
            .map(|a| *a * *c_inv)
            .zip(ck_1.ck_a.as_ref())
            .map(|(a_1, a_2)| a_1 + a_2.clone())
            .collect::<Vec<IPC::LeftKey>>();
        end_timer!(rescale_a);

        let rescale_b = start_timer!(|| "Rescale CK_A");
        let ck_b = cfg_iter!(ck_1.ck_b.as_ref())
            .map(|b| *b * *c)
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
