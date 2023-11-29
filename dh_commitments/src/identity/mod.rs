use ark_ff::fields::PrimeField;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::{
    fmt::{Debug, Display},
    marker::PhantomData,
    ops::Mul,
    ops::{Add, MulAssign},
    rand::Rng,
};

use crate::{DoublyHomomorphicCommitment, Error};

#[derive(Clone)]
pub struct IdentityCommitment<T, F: PrimeField>(PhantomData<(T, F)>);

#[derive(CanonicalSerialize, CanonicalDeserialize, Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct PlaceholderKey;

impl Display for PlaceholderKey {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "PlaceholderKey")
    }
}

impl Add for PlaceholderKey {
    type Output = Self;

    fn add(self, _rhs: Self) -> Self::Output {
        PlaceholderKey {}
    }
}

impl<T> MulAssign<T> for PlaceholderKey {
    fn mul_assign(&mut self, _rhs: T) {}
}

impl<T> Mul<T> for PlaceholderKey {
    type Output = Self;
    fn mul(self, _rhs: T) -> Self {
        Self
    }
}

#[derive(CanonicalSerialize, CanonicalDeserialize, Clone, Default, Eq, PartialEq, Debug)]
pub struct IdentityOutput<T>(pub Vec<T>)
where
    T: CanonicalSerialize + CanonicalDeserialize + Clone + Default + Eq;

impl<T: Display> Display for IdentityOutput<T>
where
    T: CanonicalSerialize + CanonicalDeserialize + Clone + Default + Eq,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "[")?;
        for (i, t) in self.0.iter().enumerate() {
            write!(f, "{}", t)?;
            if i != self.0.len() - 1 {
                write!(f, ", ")?;
            }
        }
        write!(f, "]")
    }
}

impl<T> Add for IdentityOutput<T>
where
    T: Add<T, Output = T> + CanonicalSerialize + CanonicalDeserialize + Clone + Default + Eq,
{
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        IdentityOutput(
            self.0
                .iter()
                .zip(&rhs.0)
                .map(|(a, b)| a.clone() + b.clone())
                .collect::<Vec<T>>(),
        )
    }
}

impl<T, F> MulAssign<F> for IdentityOutput<T>
where
    T: MulAssign<F> + CanonicalSerialize + CanonicalDeserialize + Clone + Default + Eq,
    F: Clone,
{
    fn mul_assign(&mut self, rhs: F) {
        self.0.iter_mut().for_each(|a| a.mul_assign(rhs.clone()))
    }
}

impl<T, F> Mul<F> for IdentityOutput<T>
where
    T: MulAssign<F> + CanonicalSerialize + CanonicalDeserialize + Clone + Default + Eq,
    F: Clone,
{
    type Output = Self;
    fn mul(mut self, rhs: F) -> Self {
        self *= rhs;
        self
    }
}

impl<T, F> DoublyHomomorphicCommitment for IdentityCommitment<T, F>
where
    T: CanonicalSerialize
        + CanonicalDeserialize
        + Copy
        + Default
        + Eq
        + Add<T, Output = T>
        + Mul<F, Output = T>
        + MulAssign<F>
        + Send
        + Sync
        + Display
        + Debug,
    F: PrimeField,
{
    type Scalar = F;
    type Message = T;
    type Key = PlaceholderKey;
    type Output = IdentityOutput<T>;

    fn setup(size: usize, _: impl Rng) -> Result<Vec<Self::Key>, Error> {
        Ok(vec![PlaceholderKey {}; size])
    }

    fn commit(_k: &[Self::Key], m: &[Self::Message]) -> Result<Self::Output, Error> {
        Ok(IdentityOutput(m.to_vec()))
    }
}
