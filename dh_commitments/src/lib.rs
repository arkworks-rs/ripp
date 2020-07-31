use algebra::{
    bytes::ToBytes,
    fields::PrimeField,
    //    serialize::CanonicalSerialize,
    groups::Group,
};
use num_traits::identities::One;
use rand::Rng;
use std::{
    cmp::Eq,
    error::Error as ErrorTrait,
    ops::{Add, MulAssign},
};

pub mod afgho16;
pub mod identity;
pub mod pedersen;

pub type Error = Box<dyn ErrorTrait>;

//TODO: support CanonicalSerialize
//TODO: Using MulAssign instead of Mul because Group does not support Mul

pub trait DoublyHomomorphicCommitment: Clone {
    type Scalar: PrimeField;
    type Message: ToBytes
        + Clone
        + Default
        + Eq
        + Add<Self::Message, Output = Self::Message>
        + MulAssign<Self::Scalar>;
    type Key: ToBytes
        + Clone
        + Default
        + Eq
        + Add<Self::Key, Output = Self::Key>
        + MulAssign<Self::Scalar>;
    type Output: ToBytes
        + Clone
        + Default
        + Eq
        + Add<Self::Output, Output = Self::Output>
        + MulAssign<Self::Scalar>;

    fn setup<R: Rng>(r: &mut R, size: usize) -> Result<Vec<Self::Key>, Error>;

    fn commit(k: &[Self::Key], m: &[Self::Message]) -> Result<Self::Output, Error>;

    fn verify(k: &[Self::Key], m: &[Self::Message], com: &Self::Output) -> Result<bool, Error> {
        Ok(Self::commit(k, m)? == *com)
    }
}

// Helpers for generator commitment keys used by Pedersen and AFGHO16

pub fn random_generators<R: Rng, G: Group>(rng: &mut R, num: usize) -> Vec<G> {
    (0..num).map(|_| G::rand(rng)).collect()
}

pub fn structured_generators_scalar_power<G: Group>(
    num: usize,
    s: &G::ScalarField,
    g: &G,
) -> Vec<G> {
    let mut generators = Vec::new();
    let mut pow_s = G::ScalarField::one();
    for _ in 0..num {
        generators.push(g.mul(&pow_s));
        pow_s *= s;
    }
    generators
}
