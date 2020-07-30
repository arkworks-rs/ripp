use algebra::{
    bytes::ToBytes,
    curves::PairingEngine,
    fields::{Field, PrimeField},
    //    serialize::CanonicalSerialize,
    groups::Group,
};
use num_traits::identities::One;
use rand::Rng;
use std::{
    cmp::{Eq, PartialEq},
    error::Error as ErrorTrait,
    io::{Result as IoResult, Write},
    ops::{Add, MulAssign},
};

pub mod afgho16;
pub mod pedersen;

pub type Error = Box<dyn ErrorTrait>;

//TODO: support CanonicalSerialize

pub trait DoublyHomomorphicCommitment {
    type Scalar: PrimeField;
    type Message: ToBytes + Clone + Default + Eq + Add + MulAssign<Self::Scalar>;
    type Key: ToBytes + Clone + Default + Eq + Add + MulAssign<Self::Scalar>;
    type Output: ToBytes + Default + Eq + Add + MulAssign<Self::Scalar>;

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

// Helper wrapper type around target group commitment output in order to implement MulAssign
//TODO: PairingEngine provides target group GT implementing Group for prime order P::Fr

#[derive(Clone)]
pub struct ExtensionFieldCommitment<P: PairingEngine>(P::Fqk);

impl<P: PairingEngine> Default for ExtensionFieldCommitment<P> {
    fn default() -> Self {
        ExtensionFieldCommitment(<P::Fqk>::default())
    }
}

impl<P: PairingEngine> PartialEq for ExtensionFieldCommitment<P> {
    fn eq(&self, other: &Self) -> bool {
        self.0.eq(&other.0)
    }
}

impl<P: PairingEngine> Eq for ExtensionFieldCommitment<P> {}

impl<P: PairingEngine> MulAssign<P::Fr> for ExtensionFieldCommitment<P> {
    fn mul_assign(&mut self, rhs: P::Fr) {
        *self = ExtensionFieldCommitment(self.0.pow(rhs.into_repr()))
    }
}

impl<P: PairingEngine> Add<Self> for ExtensionFieldCommitment<P> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        ExtensionFieldCommitment(self.0.add(&rhs.0))
    }
}

impl<P: PairingEngine> ToBytes for ExtensionFieldCommitment<P> {
    fn write<W: Write>(&self, mut writer: W) -> IoResult<()> {
        self.0.write(&mut writer)
    }
}
