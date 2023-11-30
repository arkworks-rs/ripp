use ark_dh_commitments::Error;
use ark_inner_products::InnerProduct;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::{
    cfg_iter,
    fmt::Debug,
    ops::{Add, AddAssign, Mul, MulAssign},
    rand::Rng,
};

#[cfg(feature = "parallel")]
use rayon::prelude::*;

pub mod data_structures;
pub use data_structures::*;

pub mod constructions;
pub use constructions::*;

pub trait IPCommitment: Sized {
    /// The inner product this commitment scheme is for.
    type IP: InnerProduct;

    /// The commitment key for the left vector.
    type LeftKey: CanonicalSerialize
        + CanonicalDeserialize
        + Copy
        + Default
        + Debug
        + Eq
        + Send
        + Sync
        + Add<Output = Self::LeftKey>
        + Mul<Scalar<Self>, Output = Self::LeftKey>
        + MulAssign<Scalar<Self>>;

    /// The commitment key for the right vector.
    type RightKey: CanonicalSerialize
        + CanonicalDeserialize
        + Copy
        + Default
        + Debug
        + Eq
        + Send
        + Sync
        + Add<Output = Self::RightKey>
        + Mul<Scalar<Self>, Output = Self::RightKey>
        + MulAssign<Scalar<Self>>;

    /// The commitment key for the inner product `IP::Output`.
    type IPKey: CanonicalSerialize
        + CanonicalDeserialize
        + Copy
        + Default
        + Debug
        + Eq
        + Send
        + Sync
        + Add<Output = Self::IPKey>
        + Mul<Scalar<Self>, Output = Self::IPKey>
        + MulAssign<Scalar<Self>>;

    /// The type of the commitment to the left and right vectors.
    type LeftRightCommitment: CanonicalSerialize
        + CanonicalDeserialize
        + Copy
        + Debug
        + Default
        + Eq
        + Add<Output = Self::LeftRightCommitment>
        + AddAssign<Self::LeftRightCommitment>
        + Mul<Scalar<Self>, Output = Self::LeftRightCommitment>;

    /// The type of the commitment to the inner product (i.e., to `IP::Output`).
    type OutputCommitment: CanonicalSerialize
        + CanonicalDeserialize
        + Copy
        + Debug
        + Default
        + Eq
        + Add<Output = Self::OutputCommitment>
        + AddAssign<Self::OutputCommitment>
        + Mul<Scalar<Self>, Output = Self::OutputCommitment>;

    fn setup<'a>(size: usize, r: impl Rng) -> Result<IPCommKey<'a, Self>, Error>;

    /// Commit to the left and right messages, and the inner product.
    /// `ip` should be `None` if the inner product *should* not be committed to. This
    /// is useful when computing a GIPA for a *twisted* inner product.
    fn commit_with_ip<'a>(
        ck: &IPCommKey<'a, Self>,
        l: &[LeftMessage<Self>],
        r: &[RightMessage<Self>],
        ip: impl Into<Option<OutputMessage<Self>>>,
    ) -> Result<Commitment<Self>, Error>;

    /// Commit to the left and right messages, and the inner product.
    /// The inner product is computed as `IP::inner_product(l, r)`.
    fn commit<'a>(
        ck: &IPCommKey<'a, Self>,
        l: &[LeftMessage<Self>],
        r: &[RightMessage<Self>],
    ) -> Result<Commitment<Self>, Error> {
        let ip = Self::IP::inner_product(l, r)?;
        Self::commit_with_ip(ck, l, r, ip)
    }

    /// Commit *only* to the inner product.
    fn commit_only_ip<'a>(
        ck: &IPCommKey<'a, Self>,
        ip: OutputMessage<Self>,
    ) -> Result<Commitment<Self>, Error> {
        Self::commit_with_ip(ck, &[], &[], ip)
    }

    /// Check that `cm` is indeed a commitment to the left and right messages and their
    /// inner product `IP::inner_product(l, r)`.
    fn verify<'a>(
        ck: &IPCommKey<'a, Self>,
        l: &[LeftMessage<Self>],
        r: &[RightMessage<Self>],
        cm: &Commitment<Self>,
    ) -> Result<bool, Error> {
        Ok(Self::commit(ck, l, r)? == *cm)
    }

    fn left_key_msm<'a>(
        ck: &[Self::LeftKey],
        scalars: &[Scalar<Self>],
    ) -> Result<Self::LeftKey, Error> {
        #[cfg(feature = "parallel")]
        let zero = Self::LeftKey::default;

        #[cfg(not(feature = "parallel"))]
        let zero = Self::LeftKey::default();

        Ok(cfg_iter!(ck)
            .zip(scalars)
            .map(|(a, s)| *a * *s)
            .reduce(zero, |a, b| a + b))
    }

    fn right_key_msm<'a>(
        ck: &[Self::RightKey],
        scalars: &[Scalar<Self>],
    ) -> Result<Self::RightKey, Error> {
        #[cfg(feature = "parallel")]
        let zero = Self::RightKey::default;

        #[cfg(not(feature = "parallel"))]
        let zero = Self::RightKey::default();

        Ok(cfg_iter!(ck)
            .zip(scalars)
            .map(|(a, s)| *a * *s)
            .reduce(zero, |a, b| a + b))
    }
}
