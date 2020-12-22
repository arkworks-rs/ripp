use ark_ec::{msm::VariableBaseMSM, PairingEngine, ProjectiveCurve};
use ark_ff::{bytes::ToBytes, Field, PrimeField};
use ark_std::{cfg_into_iter, cfg_iter};
use std::{
    error::Error as ErrorTrait,
    fmt::{Display, Formatter, Result as FmtResult},
    io::{Result as IoResult, Write},
    marker::PhantomData,
    ops::{Add, Mul, MulAssign},
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
pub struct PairingInnerProduct<P: PairingEngine> {
    _pair: PhantomData<P>,
}

impl<P: PairingEngine> InnerProduct for PairingInnerProduct<P> {
    type LeftMessage = P::G1Projective;
    type RightMessage = P::G2Projective;
    type Output = ExtensionFieldElement<P>;

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
        let aff_left = P::G1Projective::batch_normalization_into_affine(left);
        let aff_right = P::G2Projective::batch_normalization_into_affine(right);
        let aff_pairs = cfg_into_iter!(aff_left)
            .zip(aff_right)
            .map(|(a, b)| (P::G1Prepared::from(a), P::G2Prepared::from(b)))
            .collect::<Vec<_>>();
        Ok(ExtensionFieldElement(P::product_of_pairings(&aff_pairs)))
    }
}

#[derive(Copy, Clone)]
pub struct MultiexponentiationInnerProduct<G: ProjectiveCurve> {
    _projective: PhantomData<G>,
}

impl<G: ProjectiveCurve> InnerProduct for MultiexponentiationInnerProduct<G> {
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
        let right_bigints = cfg_iter!(right).map(|b| b.into_repr()).collect::<Vec<_>>();
        Ok(VariableBaseMSM::multi_scalar_mul(
            &G::batch_normalization_into_affine(left),
            &right_bigints,
        ))
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

// Helper wrapper type around target group commitment output in order to implement MulAssign (needed for dh_commitments)
//TODO: PairingEngine provides target group GT implementing Group for prime order P::Fr

#[derive(Clone, Debug)]
pub struct ExtensionFieldElement<P: PairingEngine>(pub P::Fqk);

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
        ExtensionFieldElement(<P::Fqk as Mul>::mul(self.0, rhs.0))
    }
}

impl<P: PairingEngine> ToBytes for ExtensionFieldElement<P> {
    fn write<W: Write>(&self, mut writer: W) -> IoResult<()> {
        self.0.write(&mut writer)
    }
}
