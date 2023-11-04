use std::marker::PhantomData;

use ark_dh_commitments::{identity::IdentityCommitment, pedersen::PedersenCommitment, afgho16::{AFGHOCommitmentG2, AFGHOCommitmentG1}};
use ark_ec::pairing::{Pairing, PairingOutput};
use ark_inner_products::PairingInnerProduct;

struct TIPPCommitment<E: Pairing>(PhantomData<E>);
type GC1<E> = AFGHOCommitmentG1<E>;
type GC2<E> = AFGHOCommitmentG2<E>;
type SC1<E> = PedersenCommitment<<E as Pairing>::G1>;
type SC2<E> = PedersenCommitment<<E as Pairing>::G2>;
type IP<E> = PairingInnerProduct<E>;
type IPC<E> =
    IdentityCommitment<PairingOutput<E>, <E as Pairing>::ScalarField>;



pub trait IPCommitment {
    type IP: InnerProduct;
    type Scalar;
    type LeftKey: CanonicalSerialize
        + CanonicalDeserialize
        + Clone
        + Default
        + Eq
        + Send
        + Sync
        + Add<Self::LeftKey, Output = Self::LeftKey>
        + MulAssign<Self::Scalar>;
    type RightKey: CanonicalSerialize
        + CanonicalDeserialize
        + Clone
        + Default
        + Eq
        + Send
        + Sync
        + Add<Self::RightKey, Output = Self::RightKey>
        + MulAssign<Self::Scalar>;
    type IPKey: CanonicalSerialize
        + CanonicalDeserialize
        + Clone
        + Default
        + Eq
        + Send
        + Sync
        + Add<Self::IPKey, Output = Self::IPKey>
        + MulAssign<Self::Scalar>;

    type Output: CanonicalSerialize
        + CanonicalDeserialize
        + Clone
        + Default
        + Eq
        + Add<Self::Output, Output = Self::Output>
        + MulAssign<Self::Scalar>;
    
    fn setup(size: usize, r: &mut impl Rng) -> Result<IPCommKey<'_, Self>, Error>;

    fn commit<'a>(ck: &IPCommKey<'a, Self>, l: &[LeftMessage<Self>], r: &[RightMessage<Self>], ip: &[OutputMessage<Self>]) -> Result<Self::Output, Error>;

    fn verify<'a>(ck: &IPCommKey<'a, Self>, l: &[LeftMessage<Self>], r: &[RightMessage<Self>], ip: &[OutputMessage<Self>], com: &Self::Output) -> Result<bool, Error> {
        Ok(Self::commit(ck, l, r, ip)? == *com)
    }
}