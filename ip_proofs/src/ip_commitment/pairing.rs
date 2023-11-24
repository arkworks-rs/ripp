use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::ops::{Add, Mul};
use ark_std::{marker::PhantomData, rand::Rng, UniformRand};

use ark_ec::pairing::{Pairing, PairingOutput};
use ark_inner_products::{multi_pairing, PairingInnerProduct};

use crate::Error;

use super::{
    identity::PlaceholderKey, IPCommKey, IPCommitment, LeftMessage, OutputMessage, RightMessage,
};

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
            com_t: self
                .com_t
                .iter()
                .zip(rhs.com_t.iter())
                .map(|(l, r)| l + r)
                .collect(),
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

impl<E: Pairing> IPCommitment for PairingCommitment<E> {
    type IP = PairingInnerProduct<E>;

    type LeftKey = E::G2;
    type RightKey = E::G1;
    type IPKey = PlaceholderKey;

    type Commitment = PairingCommOutput<E>;

    fn setup<'a>(size: usize, mut rng: impl Rng) -> Result<IPCommKey<'a, Self>, Error> {
        let random_left_key: Vec<E::G2> = (0..size).map(|_| E::G2::rand(&mut rng)).collect();
        let random_right_key: Vec<E::G1> = (0..size).map(|_| E::G1::rand(&mut rng)).collect();

        Ok(IPCommKey {
            ck_a: random_left_key.into(),
            ck_b: random_right_key.into(),
            ck_t: vec![PlaceholderKey].into(),
        })
    }

    fn commit<'a>(
        ck: &IPCommKey<'a, Self>,
        l: &[LeftMessage<Self>],
        r: &[RightMessage<Self>],
        ip: &[OutputMessage<Self>],
    ) -> Result<Self::Commitment, Error> {
        let com_a = multi_pairing(l, &ck.ck_a).ok_or("invalid pairing")?;
        let com_b = multi_pairing(&ck.ck_b, r).ok_or("invalid pairing")?;
        let com_t = ip.to_vec();

        Ok(PairingCommOutput {
            com_a,
            com_b,
            com_t,
        })
    }
}
