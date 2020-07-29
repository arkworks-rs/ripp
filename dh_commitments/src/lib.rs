use algebra::{
    bytes::ToBytes,
    serialize::CanonicalSerialize,
    groups::Group,
    fields::Field,
};
use std::{
    error::Error as ErrorTrait,
    hash::Hash,
    fmt,
    ops::{Add, Deref},
    io::{Result as IoResult, Write},
};
use rand::Rng;
use num_traits::identities::One;

pub mod pedersen;
pub mod afgho16;

pub type Error = Box<dyn ErrorTrait>;


#[derive(Debug)]
pub enum CommitmentError {
    KeyMessageLengthInvalid(usize, usize),
}

impl ErrorTrait for CommitmentError {
    fn source(self: &Self) -> Option<&(dyn ErrorTrait + 'static)> {
        None
    }
}

impl fmt::Display for CommitmentError {
    fn fmt(self: &Self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let msg = match self {
            CommitmentError::KeyMessageLengthInvalid(kl, ml) => format!("key length, message length: {}, {}", kl, ml),
        };
        write!(f, "{}", msg)
    }
}


//TODO: support CanonicalSerialize + Add traits

pub trait DoublyHomomorphicCommitment {
    type Message: ToBytes + Clone + Default + Eq + Add;
    type Key: ToBytes + Clone + Default + Eq + Add;
    type Output: ToBytes + Default + Eq + Add;

    fn setup<R: Rng>(r: &mut R, size: usize) -> Result<Vec<Self::Key>, Error>;

    fn commit(
        k: &[Self::Key],
        m: &[Self::Message],
    ) -> Result<Self::Output, Error>;

    fn verify(
        k: &[Self::Key],
        m: &[Self::Message],
        com: &Self::Output,
    ) -> Result<bool, Error> {
        Ok(Self::commit(k, m)? == *com)
    }
}

// Helpers for generator commitment keys used by Pedersen and AFGHO16

pub fn random_generators<R: Rng, G: Group>(rng: &mut R, num: usize) -> Vec<G> {
    (0..num).map(|_| G::rand(rng)).collect()
}

pub fn structured_generators_scalar_power<G: Group>(num: usize, s: &G::ScalarField, g: &G) -> Vec<G> {
    let mut generators = Vec::new();
    let mut pow_s = G::ScalarField::one();
    for _ in 0..num {
        generators.push(g.mul(&pow_s));
        pow_s *= s;
    }
    generators
}
