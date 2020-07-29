use algebra::{
    bytes::ToBytes,
    serialize::CanonicalSerialize,
    groups::Group,
};
use std::marker::PhantomData;
use rand::Rng;

use crate::{
    Error,
    DoublyHomomorphicCommitment,
    GeneratorCommitmentKey, ScalarMessage, GroupElementCommitment,
    random_generators,
};

pub trait PedersenConfig {
    const SIZE: usize;
}

pub struct PedersenCommitment<G: Group, C: PedersenConfig> {
    _group: PhantomData<G>,
    _config: PhantomData<C>,
}

impl<G: Group, C: PedersenConfig> DoublyHomomorphicCommitment for PedersenCommitment<G, C> {
    type Message = ScalarMessage<G::ScalarField>;
    type Key = GeneratorCommitmentKey<G>;
    type Output = GroupElementCommitment<G>;

    fn setup<R: Rng>(rng: &mut R) -> Result<Self::Key, Error> {
        Ok(GeneratorCommitmentKey(random_generators(rng, C::SIZE)))
    }

    fn commit(k: &Self::Key, m: &Self::Message) -> Result<Self::Output, Error> {
        Ok(GroupElementCommitment(
            k.0.iter().zip(&m.0)
                .map(|(g, x)| g.mul(x))
                .sum()
        ))
    }
}
