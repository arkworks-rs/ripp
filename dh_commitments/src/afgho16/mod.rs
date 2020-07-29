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
    GeneratorCommitmentKey, GroupElementMessage, ExtensionFieldCommitment,
    random_generators,
};

pub trait AFGHOConfig {
    const SIZE: usize;
}

pub struct AFGHOCommitment<P: PairingEngine, C: AFGHOConfig> {
    _pair: PhantomData<P>,
    _config: PhantomData<C>,
}

pub struct AFGHOCommitmentG1<P: PairingEngine, C: AFGHOConfig>(AFGHOCommitment<P, C>);
pub struct AFGHOCommitmentG2<P: PairingEngine, C: AFGHOConfig>(AFGHOCommitment<P, C>);

impl<P: PairingEngine, C: AFGHOConfig> DoublyHomomorphicCommitment for AFGHOCommitmentG1<P, C>
{
    type Message = GroupElementMessage<P::G1Projective>;
    type Key = GeneratorCommitmentKey<P::G2Projective>;
    type Output = ExtensionFieldCommitment<P::Fqk>;

    fn setup<R: Rng>(rng: &mut R) -> Result<Self::Key, Error> {
        Ok(GeneratorCommitmentKey(random_generators(rng, C::SIZE)))
    }

    fn commit(k: &Self::Key, m: &Self::Message) -> Result<Self::Output, Error> {
        Ok(ExtensionFieldCommitment(
                k.0.iter().zip(&m.0)
                    .map(|(v, a)| P::pairing(a.clone().into(), v.clone().into()))
                    .product()
        ))
    }
}

impl<P: PairingEngine, C: AFGHOConfig> DoublyHomomorphicCommitment for AFGHOCommitmentG2<P, C>
{
    type Message = GroupElementMessage<P::G2Projective>;
    type Key = GeneratorCommitmentKey<P::G1Projective>;
    type Output = ExtensionFieldCommitment<P::Fqk>;

    fn setup<R: Rng>(rng: &mut R) -> Result<Self::Key, Error> {
        Ok(GeneratorCommitmentKey(random_generators(rng, C::SIZE)))
    }

    fn commit(k: &Self::Key, m: &Self::Message) -> Result<Self::Output, Error> {
        Ok(ExtensionFieldCommitment(
            k.0.iter().zip(&m.0)
                .map(|(v, a)| P::pairing(v.clone().into(), a.clone().into()))
                .product()
        ))
    }
}
