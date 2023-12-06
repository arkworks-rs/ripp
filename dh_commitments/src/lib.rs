use ark_ec::Group;
use ark_ff::fields::PrimeField;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::{
    cmp::Eq,
    error::Error as ErrorTrait,
    fmt::{Debug, Display},
    ops::{Add, AddAssign, Mul, MulAssign},
    rand::Rng, cfg_into_iter,
};

pub mod afgho16;
pub mod identity;
pub mod pedersen;

pub type Error = Box<dyn ErrorTrait>;

#[cfg(feature = "parallel")]
use rayon::prelude::*;

//TODO: support CanonicalSerialize
//TODO: Using MulAssign instead of Mul because Group does not support Mul

pub trait DoublyHomomorphicCommitment: Clone {
    type Scalar: PrimeField;
    type Message: CanonicalSerialize
        + CanonicalDeserialize
        + Copy
        + Debug
        + Display
        + Default
        + Eq
        + Send
        + Sync
        + Add<Self::Message, Output = Self::Message>
        + MulAssign<Self::Scalar>
        + Mul<Self::Scalar, Output = Self::Message>;

    type Key: CanonicalSerialize
        + CanonicalDeserialize
        + Copy
        + Display
        + Debug
        + Default
        + Eq
        + Send
        + Sync
        + Add<Self::Key, Output = Self::Key>
        + MulAssign<Self::Scalar>
        + Mul<Self::Scalar, Output = Self::Key>;

    type Output: CanonicalSerialize
        + CanonicalDeserialize
        + Display
        + Copy
        + Debug
        + Default
        + Eq
        + AddAssign<Self::Output>
        + Add<Self::Output, Output = Self::Output>
        + MulAssign<Self::Scalar>
        + Mul<Self::Scalar, Output = Self::Output>;

    fn setup(size: usize, r: impl Rng) -> Result<Vec<Self::Key>, Error>;

    fn commit(k: &[Self::Key], m: &[Self::Message]) -> Result<Self::Output, Error>;

    fn verify(k: &[Self::Key], m: &[Self::Message], com: &Self::Output) -> Result<bool, Error> {
        Ok(Self::commit(k, m)? == *com)
    }
}

// Helpers for generator commitment keys used by Pedersen and AFGHO16

pub fn random_generators<R: Rng, G: Group>(rng: &mut R, num: usize) -> Vec<G> {
    use ark_std::rand::SeedableRng;
    use rand_chacha::ChaChaRng;
    use blake2::Blake2s256;
    use ark_std::UniformRand;
    use blake2::Digest;

    let seed = <[u8; 32]>::rand(rng);

    cfg_into_iter!(0..num).map(|i| {
        let mut hasher = Blake2s256::new();
        hasher.update(seed);
        hasher.update(&(i as u32).to_le_bytes());
        let seed: [u8; 32] = hasher.finalize().as_slice().try_into().unwrap();
        let mut rng = ChaChaRng::from_seed(seed);
        G::rand(&mut rng)
    }).collect()
}
