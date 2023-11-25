use ark_ec::CurveGroup;
use ark_std::rand::Rng;
use std::marker::PhantomData;

use crate::{random_generators, DoublyHomomorphicCommitment, Error};

use ark_inner_products::{InnerProduct, MSMInnerProduct};

#[derive(Clone)]
pub struct PedersenCommitment<G: CurveGroup> {
    _group: PhantomData<G>,
}

impl<G: CurveGroup> DoublyHomomorphicCommitment for PedersenCommitment<G> {
    type Scalar = G::ScalarField;
    type Message = G::ScalarField;
    type Key = G;
    type Output = G;

    fn setup(size: usize, mut rng: impl Rng) -> Result<Vec<Self::Key>, Error> {
        Ok(random_generators(&mut rng, size))
    }

    fn commit(k: &[Self::Key], m: &[Self::Message]) -> Result<Self::Output, Error> {
        Ok(MSMInnerProduct::<G>::inner_product(k, m)?)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_ed_on_bls12_381::{EdwardsProjective as JubJub, Fr};
    use ark_ff::UniformRand;

    type C = PedersenCommitment<JubJub>;
    const TEST_SIZE: usize = 8;

    #[test]
    fn pedersen_test() {
        let mut rng = ark_std::test_rng();
        let commit_keys = C::setup(TEST_SIZE, &mut rng).unwrap();
        let mut message = Vec::new();
        let mut wrong_message = Vec::new();
        for _ in 0..TEST_SIZE {
            message.push(Fr::rand(&mut rng));
            wrong_message.push(Fr::rand(&mut rng));
        }
        let com = C::commit(&commit_keys, &message).unwrap();
        assert!(C::verify(&commit_keys, &message, &com).unwrap());
        assert!(!C::verify(&commit_keys, &wrong_message, &com).unwrap());
        message.push(Fr::rand(&mut rng));
        assert!(C::verify(&commit_keys, &message, &com).is_err());
    }
}
