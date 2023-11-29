use ark_dh_commitments::DoublyHomomorphicCommitment;
/// Generic construction of an inner-product commitment from a doubly-homomorphic commitment.
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::{
    borrow::Cow,
    marker::PhantomData,
    ops::{Add, Mul},
    rand::Rng,
};

use ark_inner_products::InnerProduct;
use derivative::Derivative;

use crate::Error;

use super::{IPCommKey, IPCommitment, LeftMessage, OutputMessage, RightMessage};

/// Represents a commitment scheme over (G1, G2) which does a pairing commitment on the LHS and a
/// pairing commitment on the RHS, and an identity commitment on the inner product of the two.
pub struct GenericCommitment<IP, LC, RC, OC>
where
    IP: InnerProduct,
    LC: DoublyHomomorphicCommitment<Scalar = IP::Scalar, Message = IP::LeftMessage>,
    RC: DoublyHomomorphicCommitment<Scalar = IP::Scalar, Message = IP::RightMessage>,
    OC: DoublyHomomorphicCommitment<Scalar = IP::Scalar, Message = IP::Output>,
{
    _marker: PhantomData<(IP, LC, RC, OC)>,
}

#[derive(Derivative, CanonicalSerialize, CanonicalDeserialize)]
#[derivative(Clone, Debug, PartialEq, Eq, Default)]
pub struct CommOutput<IP, LC, RC, OC>
where
    IP: InnerProduct,
    LC: DoublyHomomorphicCommitment<Scalar = IP::Scalar, Message = IP::LeftMessage>,
    RC: DoublyHomomorphicCommitment<Scalar = IP::Scalar, Message = IP::RightMessage>,
    OC: DoublyHomomorphicCommitment<Scalar = IP::Scalar, Message = IP::Output>,
{
    #[derivative(Debug(format_with = "ark_std::fmt::Display::fmt"))]
    com_a: LC::Output,
    #[derivative(Debug(format_with = "ark_std::fmt::Display::fmt"))]
    com_b: RC::Output,
    #[derivative(Debug(format_with = "ark_std::fmt::Display::fmt"))]
    com_t: OC::Output,
    #[derivative(Debug = "ignore")]
    ip: PhantomData<fn(IP) -> IP>,
}

impl<IP, LC, RC, OC> Add for CommOutput<IP, LC, RC, OC>
where
    IP: InnerProduct,
    LC: DoublyHomomorphicCommitment<Scalar = IP::Scalar, Message = IP::LeftMessage>,
    RC: DoublyHomomorphicCommitment<Scalar = IP::Scalar, Message = IP::RightMessage>,
    OC: DoublyHomomorphicCommitment<Scalar = IP::Scalar, Message = IP::Output>,
{
    type Output = Self;
    fn add(self, rhs: Self) -> Self::Output {
        Self {
            com_a: self.com_a + rhs.com_a,
            com_b: self.com_b + rhs.com_b,
            com_t: self.com_t + rhs.com_t,
            ip: PhantomData,
        }
    }
}
impl<IP, LC, RC, OC> Mul<IP::Scalar> for CommOutput<IP, LC, RC, OC>
where
    IP: InnerProduct,
    LC: DoublyHomomorphicCommitment<Scalar = IP::Scalar, Message = IP::LeftMessage>,
    RC: DoublyHomomorphicCommitment<Scalar = IP::Scalar, Message = IP::RightMessage>,
    OC: DoublyHomomorphicCommitment<Scalar = IP::Scalar, Message = IP::Output>,
{
    type Output = Self;

    fn mul(self, rhs: IP::Scalar) -> Self::Output {
        Self {
            com_a: self.com_a * rhs,
            com_b: self.com_b * rhs,
            com_t: self.com_t * rhs,
            ip: PhantomData,
        }
    }
}

impl<IP, LC, RC, OC> IPCommitment for GenericCommitment<IP, LC, RC, OC>
where
    IP: InnerProduct,
    LC: DoublyHomomorphicCommitment<Scalar = IP::Scalar, Message = IP::LeftMessage>,
    RC: DoublyHomomorphicCommitment<Scalar = IP::Scalar, Message = IP::RightMessage>,
    OC: DoublyHomomorphicCommitment<Scalar = IP::Scalar, Message = IP::Output>,
{
    type IP = IP;

    type LeftKey = LC::Key;
    type RightKey = RC::Key;
    type IPKey = OC::Key;

    type Commitment = CommOutput<IP, LC, RC, OC>;

    fn setup<'a>(size: usize, mut rng: impl Rng) -> Result<IPCommKey<'a, Self>, Error> {
        let left_ck = LC::setup(size, &mut rng)?;
        let right_ck = RC::setup(size, &mut rng)?;
        let ip_ck = OC::setup(1, &mut rng)?[0];

        Ok(IPCommKey {
            ck_a: left_ck.into(),
            ck_b: right_ck.into(),
            ck_t: Cow::Owned(ip_ck),
        })
    }

    fn commit<'a>(
        ck: &IPCommKey<'a, Self>,
        l: &[LeftMessage<Self>],
        r: &[RightMessage<Self>],
        ip: impl Fn() -> OutputMessage<Self>,
    ) -> Result<Self::Commitment, Error> {
        let com_a = LC::commit(&ck.ck_a, l)?;
        let com_b = RC::commit(&ck.ck_b, r)?;
        let com_t = OC::commit(ark_std::slice::from_ref(ck.ck_t.as_ref()), &[ip()])?;

        Ok(CommOutput {
            com_a,
            com_b,
            com_t,
            ip: PhantomData,
        })
    }
}
