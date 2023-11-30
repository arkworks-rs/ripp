use ark_ff::fields::PrimeField;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::{
    fmt::{Debug, Display},
    marker::PhantomData,
    ops::{Add, AddAssign, Mul, MulAssign},
    rand::Rng,
};

use crate::{DoublyHomomorphicCommitment, Error};

#[derive(Clone)]
pub struct IdentityCommitment<T, F>(PhantomData<(T, F)>);

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

impl AddAssign for PlaceholderKey {
    fn add_assign(&mut self, _rhs: Self) {}
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

#[derive(CanonicalSerialize, CanonicalDeserialize, Copy, Clone, Default, Eq, PartialEq, Debug)]
pub struct IdentityOutput<T: CanonicalSerialize + CanonicalDeserialize>(pub T);

impl<T: Display> Display for IdentityOutput<T>
where
    T: CanonicalSerialize + CanonicalDeserialize,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl<T> Add for IdentityOutput<T>
where
    T: Add<T, Output = T> + CanonicalSerialize + CanonicalDeserialize,
{
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self(self.0 + rhs.0)
    }
}

impl<T, F> MulAssign<F> for IdentityOutput<T>
where
    T: MulAssign<F> + CanonicalSerialize + CanonicalDeserialize,
{
    fn mul_assign(&mut self, rhs: F) {
        self.0 *= rhs;
    }
}

impl<T> AddAssign<Self> for IdentityOutput<T>
where
    T: AddAssign<T> + CanonicalSerialize + CanonicalDeserialize,
{
    fn add_assign(&mut self, rhs: Self) {
        self.0 += rhs.0;
    }
}

impl<T, F> Mul<F> for IdentityOutput<T>
where
    T: MulAssign<F> + CanonicalSerialize + CanonicalDeserialize,
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
        + AddAssign<T>
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
        Ok(IdentityOutput(m[0]))
    }
}
