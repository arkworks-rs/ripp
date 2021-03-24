use ark_ec::PairingEngine;
use ark_std::rand::Rng;
use std::marker::PhantomData;

use crate::{random_generators, DoublyHomomorphicCommitment, Error};

use ark_inner_products::{ExtensionFieldElement, InnerProduct, PairingInnerProduct};

#[derive(Clone)]
pub struct AFGHOCommitment<P: PairingEngine> {
    _pair: PhantomData<P>,
}

#[derive(Clone)]
pub struct AFGHOCommitmentG1<P: PairingEngine>(AFGHOCommitment<P>);

#[derive(Clone)]
pub struct AFGHOCommitmentG2<P: PairingEngine>(AFGHOCommitment<P>);

impl<P: PairingEngine> DoublyHomomorphicCommitment for AFGHOCommitmentG1<P> {
    type Scalar = P::Fr;
    type Message = P::G1Projective;
    type Key = P::G2Projective;
    type Output = ExtensionFieldElement<P>;

    fn setup<R: Rng>(rng: &mut R, size: usize) -> Result<Vec<Self::Key>, Error> {
        Ok(random_generators(rng, size))
    }

    fn commit(k: &[Self::Key], m: &[Self::Message]) -> Result<Self::Output, Error> {
        Ok(PairingInnerProduct::<P>::inner_product(m, k)?)
    }
}

impl<P: PairingEngine> DoublyHomomorphicCommitment for AFGHOCommitmentG2<P> {
    type Scalar = P::Fr;
    type Message = P::G2Projective;
    type Key = P::G1Projective;
    type Output = ExtensionFieldElement<P>;

    fn setup<R: Rng>(rng: &mut R, size: usize) -> Result<Vec<Self::Key>, Error> {
        Ok(random_generators(rng, size))
    }

    fn commit(k: &[Self::Key], m: &[Self::Message]) -> Result<Self::Output, Error> {
        Ok(PairingInnerProduct::<P>::inner_product(k, m)?)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_bls12_381::Bls12_381;
    use ark_ff::UniformRand;
    use ark_std::rand::{rngs::StdRng, SeedableRng};

    type C1 = AFGHOCommitmentG1<Bls12_381>;
    type C2 = AFGHOCommitmentG2<Bls12_381>;
    const TEST_SIZE: usize = 8;

    #[test]
    fn afgho_g1_test() {
        let mut rng = StdRng::seed_from_u64(0u64);
        let commit_keys = C1::setup(&mut rng, TEST_SIZE).unwrap();
        let mut message = Vec::new();
        let mut wrong_message = Vec::new();
        for _ in 0..TEST_SIZE {
            message.push(<Bls12_381 as PairingEngine>::G1Projective::rand(&mut rng));
            wrong_message.push(<Bls12_381 as PairingEngine>::G1Projective::rand(&mut rng));
        }
        let com = C1::commit(&commit_keys, &message).unwrap();
        assert!(C1::verify(&commit_keys, &message, &com).unwrap());
        assert!(!C1::verify(&commit_keys, &wrong_message, &com).unwrap());
        message.push(<Bls12_381 as PairingEngine>::G1Projective::rand(&mut rng));
        assert!(C1::verify(&commit_keys, &message, &com).is_err());
    }

    #[test]
    fn afgho_g2_test() {
        let mut rng = StdRng::seed_from_u64(0u64);
        let commit_keys = C2::setup(&mut rng, TEST_SIZE).unwrap();
        let mut message = Vec::new();
        let mut wrong_message = Vec::new();
        for _ in 0..TEST_SIZE {
            message.push(<Bls12_381 as PairingEngine>::G2Projective::rand(&mut rng));
            wrong_message.push(<Bls12_381 as PairingEngine>::G2Projective::rand(&mut rng));
        }
        let com = C2::commit(&commit_keys, &message).unwrap();
        assert!(C2::verify(&commit_keys, &message, &com).unwrap());
        assert!(!C2::verify(&commit_keys, &wrong_message, &com).unwrap());
        message.push(<Bls12_381 as PairingEngine>::G2Projective::rand(&mut rng));
        assert!(C2::verify(&commit_keys, &message, &com).is_err());
    }
}
