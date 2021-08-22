use ark_ff::PrimeField;
use ark_serialize::CanonicalSerialize;
use ark_std::rand::{rngs::StdRng, SeedableRng};

// Convenience functions for generateing Fiat-Shamir challenges
pub(crate) trait TranscriptProtocol {
    /// Appends a CanonicalSerialize-able element to the transcript. Panics on serialization error.
    fn append_serializable<S>(&mut self, label: &'static [u8], val: &S)
    where
        S: CanonicalSerialize + ?Sized;

    /// Produces a pseudorandom field element from the current transcript
    fn challenge_scalar<F: PrimeField>(&mut self, label: &'static [u8]) -> F;
}

impl TranscriptProtocol for merlin::Transcript {
    /// Appends a CanonicalSerialize-able element to the transcript. Panics on serialization error.
    fn append_serializable<S>(&mut self, label: &'static [u8], val: &S)
    where
        S: CanonicalSerialize + ?Sized,
    {
        // Serialize the input and give it to the transcript
        let mut buf = Vec::new();
        val.serialize(&mut buf).unwrap();
        self.append_message(label, &buf);
    }

    /// Produces a pseudorandom field element from the current transcript
    fn challenge_scalar<F: PrimeField>(&mut self, label: &'static [u8]) -> F {
        // Fill a buf with random bytes
        let mut buf = <<StdRng as SeedableRng>::Seed as Default>::default();
        self.challenge_bytes(label, &mut buf);

        // Use the buf to make an RNG. Then use that RNG to generate a field element
        let mut rng = StdRng::from_seed(buf);
        let f = F::rand(&mut rng);

        let elems = &[f, f, f, f, f, f, f, f, f, f, f];
        let mut buf = Vec::new();
        elems.serialize(&mut buf).unwrap();

        f
    }
}
