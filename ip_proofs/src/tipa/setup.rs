impl<P, D> TIPA<P, D>
where
    D: Digest,
    P: Pairing,
{
    pub fn setup<'a>(
        size: usize,
        rng: &mut impl Rng,
    ) -> Result<(ProverKey<P>, VerifierKey<P>), Error> {
        let srs = GenericSRS::sample(size, rng);
        Ok(srs.specialize(size))
    }
}