use ark_ec::pairing::Pairing;
use ark_std::rand::Rng;
use digest::Digest;

use crate::Error;

use super::{
    data_structures::{specialize, GenericSRS},
    *,
};

impl<P, D> TIPA<P, D>
where
    D: Digest,
    P: Pairing,
{
    pub fn setup<'a>(
        size: usize,
        rng: impl Rng,
    ) -> Result<(ProverKey<'a, P>, VerifierKey<P>), Error> {
        let srs = GenericSRS::sample(size, rng);
        Ok(specialize(srs, size))
    }
}
