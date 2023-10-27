use ark_std::rand::{RngCore, SeedableRng};
use digest::{
    generic_array::{typenum::U32, GenericArray},
    Digest,
};
use rand_chacha::ChaChaRng;
use std::marker::PhantomData;

/// A `SeedableRng` that refreshes its seed by hashing together the previous seed
/// and the new seed material.
// TODO: later: re-evaluate decision about ChaChaRng
pub struct FiatShamirRng<D>
where
    D: Digest<OutputSize = U32>,
{
    r: ChaChaRng,
    seed: GenericArray<u8, D::OutputSize>,
    #[doc(hidden)]
    digest: PhantomData<D>,
}

impl<D> RngCore for FiatShamirRng<D>
where
    D: Digest<OutputSize = U32>,
{
    #[inline]
    fn next_u32(&mut self) -> u32 {
        self.r.next_u32()
    }

    #[inline]
    fn next_u64(&mut self) -> u64 {
        self.r.next_u64()
    }

    #[inline]
    fn fill_bytes(&mut self, dest: &mut [u8]) {
        self.r.fill_bytes(dest);
    }

    #[inline]
    fn try_fill_bytes(&mut self, dest: &mut [u8]) -> Result<(), ark_std::rand::Error> {
        Ok(self.r.fill_bytes(dest))
    }
}

impl<D> FiatShamirRng<D>
where
    D: Digest<OutputSize = U32>,
{
    /// Create a new `Self` by initialzing with a fresh seed.
    /// `self.seed = H(self.seed || new_seed)`.
    #[inline]
    pub fn from_seed(seed: &[u8]) -> Self {
        let digested_seed = D::digest(&seed);
        let r = ChaChaRng::from_seed(digested_seed.into());
        Self {
            r,
            seed: digested_seed,
            digest: PhantomData,
        }
    }

    /// Refresh `self.seed` with new material. Achieved by setting
    /// `self.seed = H(self.seed || new_seed)`.
    #[inline]
    pub fn absorb(&mut self, seed: &[u8]) {
        let mut bytes = seed.to_vec();
        bytes.extend_from_slice(&self.seed);
        self.seed = D::digest(&bytes);
        self.r = ChaChaRng::from_seed(self.seed.into());
    }
}
