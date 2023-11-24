impl<E: Pairing> IPCommitment for PairingCommitment<E> {
    type IP = PairingInnerProduct<E>;

    type LeftKey = E::G2;
    type RightKey = E::G1;
    type IPKey = HomomorphicPlaceholderValue;

    type Commitment = PairingCommOutput<E>;

    fn setup<'a>(size: usize, mut rng: impl Rng) -> Result<IPCommKey<'a, Self>, Error> {
        let random_left_key: Vec<E::G2> = (0..size).map(|_| E::G2::rand(&mut rng)).collect();
        let random_right_key: Vec<E::G1> = (0..size).map(|_| E::G1::rand(&mut rng)).collect();

        Ok(IPCommKey {
            ck_a: random_left_key.into(),
            ck_b: random_right_key.into(),
            ck_t: vec![HomomorphicPlaceholderValue].into(),
        })
    }

    fn commit<'a>(
        ck: &IPCommKey<'a, Self>,
        l: &[LeftMessage<Self>],
        r: &[RightMessage<Self>],
        ip: &[OutputMessage<Self>],
    ) -> Result<Self::Commitment, Error> {
        let com_a = multi_pairing(l, &ck.ck_a).ok_or("invalid pairing")?;
        let com_b = multi_pairing(&ck.ck_b, r).ok_or("invalid pairing")?;
        let com_t = ip.to_vec();

        Ok(PairingCommOutput {
            com_a,
            com_b,
            com_t,
        })
    }
}
