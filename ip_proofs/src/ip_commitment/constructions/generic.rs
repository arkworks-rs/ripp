use ark_dh_commitments::{
    identity::{IdentityCommitment, IdentityOutput, PlaceholderKey},
    DoublyHomomorphicCommitment,
};
/// Generic construction of an inner-product commitment from a doubly-homomorphic commitment.
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::{
    borrow::Cow,
    marker::PhantomData,
    ops::{Add, AddAssign, Mul},
    rand::Rng,
};

use ark_inner_products::InnerProduct;
use derivative::Derivative;

use crate::Error;

use crate::ip_commitment::{Commitment, IPCommKey, IPCommitment, LeftMessage, RightMessage};

/// Represents a commitment scheme over (G1, G2) which does a pairing commitment on the LHS and a
/// pairing commitment on the RHS, and an identity commitment on the inner product of the two.
pub struct GenericCommitment<IP, LC, RC>
where
    IP: InnerProduct,
    LC: DoublyHomomorphicCommitment<Scalar = IP::Scalar, Message = IP::LeftMessage>,
    RC: DoublyHomomorphicCommitment<Scalar = IP::Scalar, Message = IP::RightMessage>,
{
    _marker: PhantomData<(IP, LC, RC)>,
}

#[derive(Derivative, CanonicalSerialize, CanonicalDeserialize)]
#[derivative(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub struct LRComm<IP, LC, RC>
where
    IP: InnerProduct,
    LC: DoublyHomomorphicCommitment<Scalar = IP::Scalar, Message = IP::LeftMessage>,
    RC: DoublyHomomorphicCommitment<Scalar = IP::Scalar, Message = IP::RightMessage>,
{
    #[derivative(Debug(format_with = "ark_std::fmt::Display::fmt"))]
    com_a: LC::Output,
    #[derivative(Debug(format_with = "ark_std::fmt::Display::fmt"))]
    com_b: RC::Output,
    #[derivative(Debug = "ignore")]
    ip: PhantomData<fn(IP) -> IP>,
}

impl<IP, LC, RC> Add for LRComm<IP, LC, RC>
where
    IP: InnerProduct,
    LC: DoublyHomomorphicCommitment<Scalar = IP::Scalar, Message = IP::LeftMessage>,
    RC: DoublyHomomorphicCommitment<Scalar = IP::Scalar, Message = IP::RightMessage>,
{
    type Output = Self;
    fn add(self, rhs: Self) -> Self::Output {
        Self {
            com_a: self.com_a + rhs.com_a,
            com_b: self.com_b + rhs.com_b,
            ip: PhantomData,
        }
    }
}

impl<IP, LC, RC> AddAssign for LRComm<IP, LC, RC>
where
    IP: InnerProduct,
    LC: DoublyHomomorphicCommitment<Scalar = IP::Scalar, Message = IP::LeftMessage>,
    RC: DoublyHomomorphicCommitment<Scalar = IP::Scalar, Message = IP::RightMessage>,
{
    fn add_assign(&mut self, rhs: Self) {
        self.com_a += rhs.com_a;
        self.com_b += rhs.com_b;
    }
}

impl<IP, LC, RC> Mul<IP::Scalar> for LRComm<IP, LC, RC>
where
    IP: InnerProduct,
    LC: DoublyHomomorphicCommitment<Scalar = IP::Scalar, Message = IP::LeftMessage>,
    RC: DoublyHomomorphicCommitment<Scalar = IP::Scalar, Message = IP::RightMessage>,
{
    type Output = Self;

    fn mul(self, rhs: IP::Scalar) -> Self::Output {
        Self {
            com_a: self.com_a * rhs,
            com_b: self.com_b * rhs,
            ip: PhantomData,
        }
    }
}

impl<IP, LC, RC> IPCommitment for GenericCommitment<IP, LC, RC>
where
    IP: InnerProduct,
    LC: DoublyHomomorphicCommitment<Scalar = IP::Scalar, Message = IP::LeftMessage>,
    RC: DoublyHomomorphicCommitment<Scalar = IP::Scalar, Message = IP::RightMessage>,
{
    type IP = IP;

    type LeftKey = LC::Key;
    type RightKey = RC::Key;
    type IPKey = PlaceholderKey;

    type LeftRightCommitment = LRComm<IP, LC, RC>;
    type OutputCommitment = IdentityOutput<IP::Output>;

    fn setup<'a>(size: usize, mut rng: impl Rng) -> Result<IPCommKey<'a, Self>, Error> {
        let left_ck = LC::setup(size, &mut rng)?;
        let right_ck = RC::setup(size, &mut rng)?;
        let ip_ck = PlaceholderKey;

        Ok(IPCommKey {
            ck_a: left_ck.into(),
            ck_b: right_ck.into(),
            ck_t: Cow::Owned(ip_ck),
        })
    }

    fn commit_with_ip<'a>(
        ck: &IPCommKey<'a, Self>,
        l: &[LeftMessage<Self>],
        r: &[RightMessage<Self>],
        ip: impl Into<Option<IP::Output>>,
    ) -> Result<Commitment<Self>, Error> {
        let com_a = LC::commit(&ck.ck_a, l)?;
        let com_b = RC::commit(&ck.ck_b, r)?;
        let com_t = ip
            .into()
            .map(|ip| IdentityCommitment::commit(ark_std::slice::from_ref(ck.ck_t.as_ref()), &[ip]))
            .transpose()?;

        Ok(Commitment {
            cm_lr: LRComm {
                com_a,
                com_b,
                ip: PhantomData,
            },
            cm_ip: com_t,
        })
    }
}
