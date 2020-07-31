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
pub mod identity;
pub mod pedersen;

pub type Error = Box<dyn ErrorTrait>;

//TODO: support CanonicalSerialize
//TODO: Using MulAssign instead of Mul because Group does not support Mul

pub trait DoublyHomomorphicCommitment {
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

// Helper wrapper type around target group commitment output in order to implement MulAssign
//TODO: PairingEngine provides target group GT implementing Group for prime order P::Fr

#[derive(Clone)]
pub struct ExtensionFieldElement<P: PairingEngine>(P::Fqk);

impl<P: PairingEngine> Default for ExtensionFieldElement<P> {
    fn default() -> Self {
        ExtensionFieldElement(<P::Fqk>::default())
    }
}

impl<P: PairingEngine> PartialEq for ExtensionFieldElement<P> {
    fn eq(&self, other: &Self) -> bool {
        self.0.eq(&other.0)
    }
}

impl<P: PairingEngine> Eq for ExtensionFieldElement<P> {}

impl<P: PairingEngine> MulAssign<P::Fr> for ExtensionFieldElement<P> {
    fn mul_assign(&mut self, rhs: P::Fr) {
        *self = ExtensionFieldElement(self.0.pow(rhs.into_repr()))
    }
}

impl<P: PairingEngine> Add<Self> for ExtensionFieldElement<P> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        ExtensionFieldElement(self.0.add(&rhs.0))
    }
}

impl<P: PairingEngine> ToBytes for ExtensionFieldElement<P> {
    fn write<W: Write>(&self, mut writer: W) -> IoResult<()> {
        self.0.write(&mut writer)
    }
}
