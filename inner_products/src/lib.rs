use ark_ec::{
    pairing::{MillerLoopOutput, Pairing, PairingOutput},
    CurveGroup,
};
use ark_ff::Field;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::{cfg_chunks, cfg_iter};
use std::{
    error::Error as ErrorTrait,
    fmt::{Display, Formatter, Result as FmtResult},
    marker::PhantomData,
    ops::{Add, Mul},
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
    type Scalar: Field;
    type LeftMessage: Clone
        + Send
        + Sync
        + Default
        + Eq
        + CanonicalSerialize
        + CanonicalDeserialize
        + Add<Output = Self::LeftMessage>
        + Mul<Self::Scalar, Output = Self::LeftMessage>;
    type RightMessage: Clone
        + Send
        + Sync
        + Default
        + Eq
        + CanonicalSerialize
        + CanonicalDeserialize
        + Add<Output = Self::RightMessage>
        + Mul<Self::Scalar, Output = Self::RightMessage>;
    type Output: Default + Eq + Clone + Send + Sync + CanonicalSerialize + CanonicalDeserialize;

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
    type Scalar = P::ScalarField;

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

        Ok(multi_pairing(left, right).unwrap())
    }
}

/// Equivalent to `P::multi_pairing`, but with more parallelism (if enabled)
pub fn multi_pairing<P: Pairing>(left: &[P::G1], right: &[P::G2]) -> Option<PairingOutput<P>> {
    // We make the input affine, then convert to prepared. We do this for speed, since the
    // conversion from projective to prepared always goes through affine.
    let aff_left = P::G1::normalize_batch(left);
    let aff_right = P::G2::normalize_batch(right);

    let left = cfg_iter!(aff_left)
        .map(P::G1Prepared::from)
        .collect::<Vec<_>>();
    let right = cfg_iter!(aff_right)
        .map(P::G2Prepared::from)
        .collect::<Vec<_>>();

    // We want to process N chunks in parallel where N is the number of threads available
    #[cfg(feature = "parallel")]
    let num_chunks = rayon::current_num_threads();
    #[cfg(not(feature = "parallel"))]
    let num_chunks = 1;

    let chunk_size = if num_chunks <= left.len() {
        left.len() / num_chunks
    } else {
        // More threads than elements. Just do it all in parallel
        1
    };

    let (left_chunks, right_chunks) = (
        cfg_chunks!(left, chunk_size),
        cfg_chunks!(right, chunk_size),
    );
    // Compute all the (partial) pairings and take the product. We have to take the product over
    // P::TargetField because MillerLoopOutput doesn't impl Product
    let ml_result = left_chunks
        .zip(right_chunks)
        .map(|(aa, bb)| P::multi_miller_loop(aa.iter().cloned(), bb.iter().cloned()).0)
        .product();

    P::final_exponentiation(MillerLoopOutput(ml_result))
}

#[derive(Copy, Clone)]
pub struct MultiexponentiationInnerProduct<G: CurveGroup> {
    _projective: PhantomData<G>,
}

impl<G: CurveGroup> InnerProduct for MultiexponentiationInnerProduct<G> {
    type LeftMessage = G;
    type RightMessage = G::ScalarField;
    type Output = G;
    type Scalar = G::ScalarField;

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
    type Scalar = F;
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
