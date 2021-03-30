use ark_ff::{bytes::ToBytes, fields::PrimeField};
use ark_std::rand::Rng;
use std::{
    io::{Result as IoResult, Write},
    marker::PhantomData,
    ops::{Add, MulAssign},
};

use crate::{DoublyHomomorphicCommitment, Error};

#[derive(Clone)]
pub struct IdentityCommitment<T, F: PrimeField> {
    _t: PhantomData<T>,
    _field: PhantomData<F>,
}

#[derive(Clone, Default, Eq, PartialEq)]
pub struct HomomorphicPlaceholderValue;

impl ToBytes for HomomorphicPlaceholderValue {
    fn write<W: Write>(&self, _writer: W) -> IoResult<()> {
        Ok(())
    }
}

impl Add for HomomorphicPlaceholderValue {
    type Output = Self;

    fn add(self, _rhs: Self) -> Self::Output {
        HomomorphicPlaceholderValue {}
    }
}

impl<T> MulAssign<T> for HomomorphicPlaceholderValue {
    fn mul_assign(&mut self, _rhs: T) {}
}

#[derive(Clone, Default, Eq, PartialEq)]
pub struct IdentityOutput<T: Clone + Default + Eq>(pub Vec<T>);

impl<T: ToBytes + Clone + Default + Eq> ToBytes for IdentityOutput<T> {
    fn write<W: Write>(&self, mut writer: W) -> IoResult<()> {
        self.0.write(&mut writer)
    }
}

impl<T: Add<T, Output = T> + Clone + Default + Eq> Add for IdentityOutput<T> {
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

impl<T: MulAssign<F> + Clone + Default + Eq, F: Clone> MulAssign<F> for IdentityOutput<T> {
    fn mul_assign(&mut self, rhs: F) {
        self.0.iter_mut().for_each(|a| a.mul_assign(rhs.clone()))
    }
}

impl<T, F> DoublyHomomorphicCommitment for IdentityCommitment<T, F>
where
    T: ToBytes + Clone + Default + Eq + Add<T, Output = T> + MulAssign<F> + Send + Sync,
    F: PrimeField,
{
    type Scalar = F;
    type Message = T;
    type Key = HomomorphicPlaceholderValue;
    type Output = IdentityOutput<T>;

    fn setup<R: Rng>(_rng: &mut R, size: usize) -> Result<Vec<Self::Key>, Error> {
        Ok(vec![HomomorphicPlaceholderValue {}; size])
    }

    fn commit(_k: &[Self::Key], m: &[Self::Message]) -> Result<Self::Output, Error> {
        Ok(IdentityOutput(m.to_vec()))
    }
}
