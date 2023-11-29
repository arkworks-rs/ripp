use std::ops::Mul;

use ark_ec::Group;
use ark_ff::fields::PrimeField;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::{
    cmp::Eq,
    error::Error as ErrorTrait,
    fmt::{Debug, Display},
    ops::{Add, MulAssign},
    rand::Rng,
};

pub mod afgho16;
pub mod identity;
pub mod pedersen;

pub type Error = Box<dyn ErrorTrait>;

//TODO: support CanonicalSerialize
//TODO: Using MulAssign instead of Mul because Group does not support Mul

pub trait DoublyHomomorphicCommitment: Clone {
    type Scalar: PrimeField;
    type Message: CanonicalSerialize
        + CanonicalDeserialize
        + Copy
        + Debug
        + Display
        + Default
        + Eq
        + Send
        + Sync
        + Add<Self::Message, Output = Self::Message>
        + MulAssign<Self::Scalar>
        + Mul<Self::Scalar, Output = Self::Message>;

    type Key: CanonicalSerialize
        + CanonicalDeserialize
        + Copy
        + Display
        + Debug
        + Default
        + Eq
        + Send
        + Sync
        + Add<Self::Key, Output = Self::Key>
        + MulAssign<Self::Scalar>
        + Mul<Self::Scalar, Output = Self::Key>;
    type Output: CanonicalSerialize
        + CanonicalDeserialize
        + Display
        + Clone
        + Debug
        + Default
        + Eq
        + Add<Self::Output, Output = Self::Output>
        + MulAssign<Self::Scalar>
        + Mul<Self::Scalar, Output = Self::Output>;

    fn setup(size: usize, r: impl Rng) -> Result<Vec<Self::Key>, Error>;

    fn commit(k: &[Self::Key], m: &[Self::Message]) -> Result<Self::Output, Error>;

    fn verify(k: &[Self::Key], m: &[Self::Message], com: &Self::Output) -> Result<bool, Error> {
        Ok(Self::commit(k, m)? == *com)
    }
}

// Helpers for generator commitment keys used by Pedersen and AFGHO16

pub fn random_generators<R: Rng, G: Group>(rng: &mut R, num: usize) -> Vec<G> {
    (0..num).map(|_| G::rand(rng)).collect()
}
