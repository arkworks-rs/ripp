use ark_ec::{
    pairing::{Pairing, PairingOutput},
    CurveGroup,
};
use ark_ff::Field;
use ark_std::cfg_iter;
use std::{
    error::Error as ErrorTrait,
    fmt::{Display, Formatter, Result as FmtResult},
    marker::PhantomData,
};

#[cfg(feature = "parallel")]
use rayon::prelude::*;

pub type Error = Box<dyn ErrorTrait>;

#[derive(Debug)]
pub enum InnerProductError {
    MessageLengthInvalid(usize, usize),
}

impl ErrorTrait for InnerProductError {
    fn source(self: &Self) -> Option<&(dyn ErrorTrait + 'static)> {
        None
    }
}

impl Display for InnerProductError {
    fn fmt(self: &Self, f: &mut Formatter<'_>) -> FmtResult {
        let msg = match self {
            InnerProductError::MessageLengthInvalid(left, right) => {
                format!("left length, right length: {}, {}", left, right)
            }
        };
        write!(f, "{}", msg)
    }
}

pub trait InnerProduct: Copy {
    type LeftMessage;
    type RightMessage;
    type Output;

    fn inner_product(
        left: &[Self::LeftMessage],
        right: &[Self::RightMessage],
    ) -> Result<Self::Output, Error>;
}

#[derive(Copy, Clone)]
pub struct PairingInnerProduct<P: Pairing> {
    _pair: PhantomData<P>,
}

impl<P: Pairing> InnerProduct for PairingInnerProduct<P> {
    type LeftMessage = P::G1;
    type RightMessage = P::G2;
    type Output = PairingOutput<P>;

    fn inner_product(
        left: &[Self::LeftMessage],
        right: &[Self::RightMessage],
    ) -> Result<Self::Output, Error> {
        if left.len() != right.len() {
            return Err(Box::new(InnerProductError::MessageLengthInvalid(
                left.len(),
                right.len(),
            )));
        };
        let aff_left = P::G1::normalize_batch(left);
        let aff_right = P::G2::normalize_batch(right);
        Ok(P::multi_pairing(&aff_left, &aff_right))
    }
}

#[derive(Copy, Clone)]
pub struct MultiexponentiationInnerProduct<G: CurveGroup> {
    _projective: PhantomData<G>,
}

impl<G: CurveGroup> InnerProduct for MultiexponentiationInnerProduct<G> {
    type LeftMessage = G;
    type RightMessage = G::ScalarField;
    type Output = G;

    fn inner_product(
        left: &[Self::LeftMessage],
        right: &[Self::RightMessage],
    ) -> Result<Self::Output, Error> {
        if left.len() != right.len() {
            return Err(Box::new(InnerProductError::MessageLengthInvalid(
                left.len(),
                right.len(),
            )));
        };

        // Can unwrap because we did the length check above
        Ok(G::msm(&G::normalize_batch(left), &right).unwrap())
    }
}

#[derive(Copy, Clone)]
pub struct ScalarInnerProduct<F: Field> {
    _field: PhantomData<F>,
}

impl<F: Field> InnerProduct for ScalarInnerProduct<F> {
    type LeftMessage = F;
    type RightMessage = F;
    type Output = F;

    fn inner_product(
        left: &[Self::LeftMessage],
        right: &[Self::RightMessage],
    ) -> Result<Self::Output, Error> {
        if left.len() != right.len() {
            return Err(Box::new(InnerProductError::MessageLengthInvalid(
                left.len(),
                right.len(),
            )));
        };
        Ok(cfg_iter!(left).zip(right).map(|(x, y)| *x * y).sum())
    }
}
