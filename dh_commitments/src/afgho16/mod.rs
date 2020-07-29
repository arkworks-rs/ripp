use algebra::{
    bytes::ToBytes,
    serialize::CanonicalSerialize,
    groups::Group,
    curves::{PairingEngine, prepare_g1, prepare_g2},
};
use std::marker::PhantomData;
use rand::Rng;

use crate::{
    Error,
    DoublyHomomorphicCommitment,
    random_generators,
};

pub struct AFGHOCommitment<P: PairingEngine> {
    _pair: PhantomData<P>,
}

pub struct AFGHOCommitmentG1<P: PairingEngine>(AFGHOCommitment<P>);
pub struct AFGHOCommitmentG2<P: PairingEngine>(AFGHOCommitment<P>);

impl<P: PairingEngine> DoublyHomomorphicCommitment for AFGHOCommitmentG1<P>
{
    type Message = P::G1Projective;
    type Key = P::G2Projective;
    type Output = P::Fqk;

    fn setup<R: Rng>(rng: &mut R, size: usize) -> Result<Vec<Self::Key>, Error> {
        Ok(random_generators(rng, size))
    }

    fn commit(k: &[Self::Key], m: &[Self::Message]) -> Result<Self::Output, Error> {
        Ok(
                k.iter().zip(m)
                    .map(|(v, a)| P::pairing(a.clone().into(), v.clone().into()))
                    .product()
        )
    }
}

impl<P: PairingEngine> DoublyHomomorphicCommitment for AFGHOCommitmentG2<P>
{
    type Message = P::G2Projective;
    type Key = P::G1Projective;
    type Output = P::Fqk;

    fn setup<R: Rng>(rng: &mut R, size: usize) -> Result<Vec<Self::Key>, Error> {
        Ok(random_generators(rng, size))
    }

    fn commit(k: &[Self::Key], m: &[Self::Message]) -> Result<Self::Output, Error> {
        Ok(
            k.iter().zip(m)
                .map(|(v, a)| P::pairing(v.clone().into(), a.clone().into()))
                .product()
        )
    }
}
