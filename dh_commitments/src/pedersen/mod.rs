use ark_ec::ProjectiveCurve;
use ark_std::rand::Rng;
use std::marker::PhantomData;

use crate::{random_generators, DoublyHomomorphicCommitment, Error};

use ark_inner_products::{InnerProduct, MultiexponentiationInnerProduct};

#[derive(Clone)]
pub struct PedersenCommitment<G: ProjectiveCurve> {
    _group: PhantomData<G>,
}

impl<G: ProjectiveCurve> DoublyHomomorphicCommitment for PedersenCommitment<G> {
    type Scalar = G::ScalarField;
    type Message = G::ScalarField;
    type Key = G;
    type Output = G;

    fn setup<R: Rng>(rng: &mut R, size: usize) -> Result<Vec<Self::Key>, Error> {
        Ok(random_generators(rng, size))
    }

    fn commit(k: &[Self::Key], m: &[Self::Message]) -> Result<Self::Output, Error> {
        Ok(MultiexponentiationInnerProduct::<G>::inner_product(k, m)?)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_ed_on_bls12_381::EdwardsProjective as JubJub;
    use ark_ff::UniformRand;
    use ark_std::rand::{rngs::StdRng, SeedableRng};

    type C = PedersenCommitment<JubJub>;
    const TEST_SIZE: usize = 8;

    #[test]
    fn pedersen_test() {
        let mut rng = StdRng::seed_from_u64(0u64);
        let commit_keys = C::setup(&mut rng, TEST_SIZE).unwrap();
        let mut message = Vec::new();
        let mut wrong_message = Vec::new();
        for _ in 0..TEST_SIZE {
            message.push(<JubJub as ProjectiveCurve>::ScalarField::rand(&mut rng));
            wrong_message.push(<JubJub as ProjectiveCurve>::ScalarField::rand(&mut rng));
        }
        let com = C::commit(&commit_keys, &message).unwrap();
        assert!(C::verify(&commit_keys, &message, &com).unwrap());
        assert!(!C::verify(&commit_keys, &wrong_message, &com).unwrap());
        message.push(<JubJub as ProjectiveCurve>::ScalarField::rand(&mut rng));
        assert!(C::verify(&commit_keys, &message, &com).is_err());
    }
}
