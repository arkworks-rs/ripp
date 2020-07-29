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
    random_generators,
};

pub struct PedersenCommitment<G: Group> {
    _group: PhantomData<G>,
}

impl<G: Group> DoublyHomomorphicCommitment for PedersenCommitment<G> {
    type Message = G::ScalarField;
    type Key = G;
    type Output = G;

    fn setup<R: Rng>(rng: &mut R, size: usize) -> Result<Vec<Self::Key>, Error> {
        Ok(random_generators(rng, size))
    }

    fn commit(k: &[Self::Key], m: &[Self::Message]) -> Result<Self::Output, Error> {
        Ok(
            k.iter().zip(m)
                .map(|(g, x)| g.mul(x))
                .sum()
        )
    }
}
