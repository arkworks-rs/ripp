use ark_dh_commitments::Error;
use ark_inner_products::{compute_powers, InnerProduct};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, Valid};
use ark_std::{
    borrow::Cow,
    cfg_iter, cfg_iter_mut, end_timer,
    ops::{Add, AddAssign, Mul},
    start_timer,
};

use derivative::Derivative;
#[cfg(feature = "parallel")]
use rayon::prelude::*;

use super::IPCommitment;

pub type LeftMessage<IPC> = <<IPC as IPCommitment>::IP as InnerProduct>::LeftMessage;
pub type RightMessage<IPC> = <<IPC as IPCommitment>::IP as InnerProduct>::RightMessage;
pub type OutputMessage<IPC> = <<IPC as IPCommitment>::IP as InnerProduct>::Output;
pub type Scalar<IPC> = <<IPC as IPCommitment>::IP as InnerProduct>::Scalar;

#[derive(Derivative)]
#[derivative(Clone(bound = ""), Debug(bound = "IPC: IPCommitment"))]
pub struct IPCommKey<'b, IPC: IPCommitment> {
    pub ck_a: Cow<'b, [IPC::LeftKey]>,
    pub ck_b: Cow<'b, [IPC::RightKey]>,
    pub ck_t: Cow<'b, IPC::IPKey>,
}

impl<'b, IPC: IPCommitment> CanonicalSerialize for IPCommKey<'b, IPC> {
    fn serialize_with_mode<W: std::io::prelude::Write>(
        &self,
        mut writer: W,
        compress: ark_serialize::Compress,
    ) -> Result<(), ark_serialize::SerializationError> {
        self.ck_a.serialize_with_mode(&mut writer, compress)?;
        self.ck_b.serialize_with_mode(&mut writer, compress)?;
        self.ck_t.serialize_with_mode(&mut writer, compress)?;
        Ok(())
    }

    fn serialized_size(&self, compress: ark_serialize::Compress) -> usize {
        self.ck_a.serialized_size(compress)
            + self.ck_b.serialized_size(compress)
            + self.ck_t.serialized_size(compress)
    }
}

impl<'b, IPC: IPCommitment> Valid for IPCommKey<'b, IPC> {
    fn check(&self) -> Result<(), ark_serialize::SerializationError> {
        Valid::batch_check(self.ck_a.iter())?;
        Valid::batch_check(self.ck_b.iter())?;
        self.ck_t.as_ref().check()?;
        Ok(())
    }
}

impl<'b, IPC: IPCommitment> CanonicalDeserialize for IPCommKey<'b, IPC> {
    fn deserialize_with_mode<R: std::io::prelude::Read>(
        mut reader: R,
        compress: ark_serialize::Compress,
        validate: ark_serialize::Validate,
    ) -> Result<Self, ark_serialize::SerializationError> {
        let ck_a = Vec::<IPC::LeftKey>::deserialize_with_mode(&mut reader, compress, validate)?;
        let ck_b = Vec::<IPC::RightKey>::deserialize_with_mode(&mut reader, compress, validate)?;
        let ck_t = IPC::IPKey::deserialize_with_mode(reader, compress, validate)?;
        Ok(Self {
            ck_a: Cow::Owned(ck_a),
            ck_b: Cow::Owned(ck_b),
            ck_t: Cow::Owned(ck_t),
        })
    }
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

    /// Modifies `self` by multiplying `self.ck_a` by powers of `c`.
    pub fn twist_in_place(&mut self, c: Scalar<IPC>) {
        let len = self.ck_b.len();
        cfg_iter_mut!(self.ck_a.to_mut())
            .zip(compute_powers(len, c))
            .for_each(|(a, c)| *a *= c);
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
            .map(|(a_1, a_2)| a_1 + *a_2)
            .collect::<Vec<IPC::LeftKey>>();
        end_timer!(rescale_a);

        let rescale_b = start_timer!(|| "Rescale CK_A");
        let ck_b = cfg_iter!(ck_1.ck_b.as_ref())
            .map(|b| *b * *c)
            .zip(ck_2.ck_b.as_ref())
            .map(|(b_1, b_2)| b_1 + *b_2)
            .collect::<Vec<IPC::RightKey>>();
        end_timer!(rescale_b);

        // TODO: Remove the into_owned here
        Ok(IPCommKey {
            ck_a: Cow::Owned(ck_a),
            ck_b: Cow::Owned(ck_b),
            ck_t: Cow::Owned(ck_1.ck_t.clone().into_owned()),
        })
    }

    pub fn trim_for_only_ip(&self) -> IPCommKey<'static, IPC> {
        IPCommKey {
            ck_a: Cow::Owned(vec![]),
            ck_b: Cow::Owned(vec![]),
            ck_t: Cow::Owned(self.ck_t.clone().into_owned()),
        }
    }
}

/// A final IP comm key is an IP comm key of length 1
#[derive(Derivative)]
#[derivative(
    Clone(bound = ""),
    Copy(bound = ""),
    Debug(bound = "IPC: IPCommitment"),
    PartialEq(bound = "IPC: IPCommitment"),
    Eq(bound = "IPC: IPCommitment")
)]
#[derive(CanonicalSerialize, CanonicalDeserialize)]
pub struct FinalIPCommKey<IPC: IPCommitment> {
    pub ck_a: IPC::LeftKey,
    pub ck_b: IPC::RightKey,
    pub ck_t: IPC::IPKey,
}

impl<'a, IPC: IPCommitment> TryFrom<IPCommKey<'a, IPC>> for FinalIPCommKey<IPC> {
    type Error = ();

    fn try_from(other: IPCommKey<'a, IPC>) -> Result<Self, ()> {
        if other.ck_a.len() != 1 || other.ck_b.len() != 1 {
            Err(())
        } else {
            Ok(Self {
                ck_a: other.ck_a.first().unwrap().clone(),
                ck_b: other.ck_b.first().unwrap().clone(),
                ck_t: other.ck_t.into_owned(),
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

#[derive(Derivative, CanonicalSerialize, CanonicalDeserialize)]
#[derivative(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub struct Commitment<IPC: IPCommitment> {
    pub cm_lr: IPC::LeftRightCommitment,
    pub cm_ip: Option<IPC::OutputCommitment>,
}

impl<IPC: IPCommitment> Add for Commitment<IPC> {
    type Output = Self;
    fn add(mut self, rhs: Self) -> Self::Output {
        self += rhs;
        self
    }
}

impl<IPC: IPCommitment> AddAssign for Commitment<IPC> {
    fn add_assign(&mut self, rhs: Self) {
        self.cm_lr += rhs.cm_lr;
        self.cm_ip = match (self.cm_ip, rhs.cm_ip) {
            (Some(a), Some(b)) => Some(a + b),
            (Some(a), None) | (None, Some(a)) => Some(a),
            (None, None) => None,
        };
    }
}

impl<IPC: IPCommitment> Mul<Scalar<IPC>> for Commitment<IPC> {
    type Output = Self;

    fn mul(mut self, rhs: Scalar<IPC>) -> Self::Output {
        self.cm_lr = self.cm_lr * rhs;
        self.cm_ip.iter_mut().for_each(|cm| *cm = *cm * rhs);
        self
    }
}

impl<'a, IPC: IPCommitment> Into<IPCommKey<'a, IPC>> for FinalIPCommKey<IPC> {
    fn into(self) -> IPCommKey<'a, IPC> {
        IPCommKey {
            ck_a: vec![self.ck_a.clone()].into(),
            ck_b: vec![self.ck_b.clone()].into(),
            ck_t: Cow::Owned(self.ck_t.clone()),
        }
    }
}
