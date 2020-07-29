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
    KeyMessageIncompatible(),
}

impl ErrorTrait for CommitmentError {
    fn source(self: &Self) -> Option<&(dyn ErrorTrait + 'static)> {
        None
    }
}

impl fmt::Display for CommitmentError {
    fn fmt(self: &Self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let msg = match self {
            CommitmentError::KeyMessageIncompatible() => "key and message incompatible".to_string(),
        };
        write!(f, "{}", msg)
    }
}


//TODO: support CanonicalSerialize + Add traits

pub trait HomomorphicAdd: Sized {
    fn add(&self, other: &Self) -> Result<Self, Error>;
}

pub trait DoublyHomomorphicCommitment {
    type Message: ToBytes + Clone + Default + Eq + HomomorphicAdd;
    type Key: ToBytes + Clone + Default + Eq + HomomorphicAdd;
    type Output: ToBytes + Default + Eq + HomomorphicAdd;

    fn setup<R: Rng>(r: &mut R) -> Result<Self::Key, Error>;

    fn commit(
        k: &Self::Key,
        m: &Self::Message,
    ) -> Result<Self::Output, Error>;

    fn verify(
        k: &Self::Key,
        m: &Self::Message,
        com: &Self::Output,
    ) -> Result<bool, Error> {
        Ok(Self::commit(k, m)? == *com)
    }
}


// Message spaces for Pedersen and AFGHO16

#[derive(Clone, Eq, PartialEq, Default)]
pub struct GroupElementMessage<G: Group>(Vec<G>);

impl<G: Group> HomomorphicAdd for GroupElementMessage<G> {
    fn add(&self, other: &Self) -> Result<Self, Error> {
        Ok(GroupElementMessage(
            self.0.iter().zip(&other.0)
                .map(|(a, b)| a.add(b))
                .collect()
        ))
    }
}

impl<G: Group> ToBytes for GroupElementMessage<G> {
    fn write<W: Write>(&self, mut writer: W) -> IoResult<()> {
        self.0.write(&mut writer)
    }
}

#[derive(Clone, Eq, PartialEq, Default)]
pub struct ScalarMessage<F: Field>(Vec<F>);

impl<F: Field> HomomorphicAdd for ScalarMessage<F> {
    fn add(&self, other: &Self) -> Result<Self, Error> {
        Ok(ScalarMessage(
            self.0.iter().zip(&other.0)
                .map(|(a, b)| a.add(b))
                .collect()
        ))
    }
}

impl<F: Field> ToBytes for ScalarMessage<F> {
    fn write<W: Write>(&self, mut writer: W) -> IoResult<()> {
        self.0.write(&mut writer)
    }
}

// Output space for Pedersen and AFGHO16

#[derive(Clone, Eq, PartialEq, Default)]
pub struct GroupElementCommitment<G: Group>(G);

impl<G: Group> HomomorphicAdd for GroupElementCommitment<G> {
    fn add(&self, other: &Self) -> Result<Self, Error> {
        Ok(GroupElementCommitment(self.0 + other.0))
    }
}

impl<G: Group> ToBytes for GroupElementCommitment<G> {
    fn write<W: Write>(&self, mut writer: W) -> IoResult<()> {
        self.0.write(&mut writer)
    }
}

//TODO: Extension field can be considered as group?
#[derive(Clone, Eq, PartialEq, Default)]
pub struct ExtensionFieldCommitment<F: Field>(F);

impl<F: Field> HomomorphicAdd for ExtensionFieldCommitment<F> {
    fn add(&self, other: &Self) -> Result<Self, Error> {
        Ok(ExtensionFieldCommitment(self.0 + other.0))
    }
}

impl<F: Field> ToBytes for ExtensionFieldCommitment<F> {
    fn write<W: Write>(&self, mut writer: W) -> IoResult<()> {
        self.0.write(&mut writer)
    }
}


// Helpers for generator commitment keys used by Pedersen and AFGHO16

pub trait GeneratorSetupConfig {
    const SIZE: usize;
}

#[derive(Clone, Eq, PartialEq, Default)]
pub struct GeneratorCommitmentKey<G: Group>(Vec<G>);

impl<G: Group> HomomorphicAdd for GeneratorCommitmentKey<G> {
    fn add(&self, other: &Self) -> Result<Self, Error> {
        Ok(GeneratorCommitmentKey(
            self.0.iter().zip(&other.0)
                .map(|(a, b)| a.add(b))
                .collect()
        ))
    }
}

impl<G: Group> ToBytes for GeneratorCommitmentKey<G> {
    fn write<W: Write>(&self, mut writer: W) -> IoResult<()> {
        self.0.write(&mut writer)
    }
}

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
