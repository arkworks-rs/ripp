use ark_ff::{Field, One};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::{
    borrow::Cow, cfg_iter, convert::TryInto, end_timer, marker::PhantomData, rand::Rng, start_timer,
};
use digest::Digest;

use crate::{
    ip_commitment::{FinalIPCommKey, IPCommKey, IPCommitment, Scalar},
    Error, InnerProductArgumentError,
};
use ark_inner_products::InnerProduct;

#[cfg(feature = "parallel")]
use rayon::prelude::*;

pub struct GIPA<IP, IPC, D>
where
    D: Digest,
    IP: InnerProduct,
    IPC: IPCommitment<IP = IP>,
{
    _inner_product: PhantomData<IP>,
    _inner_product_commitment: PhantomData<IPC>,
    _digest: PhantomData<D>,
}

#[derive(CanonicalSerialize, CanonicalDeserialize)]
pub struct GIPAProof<IP, IPC, D>
where
    D: Digest,
    IP: InnerProduct,
    IPC: IPCommitment<IP = IP>,
{
    pub(crate) r_commitment_steps: Vec<(IPC::Commitment, IPC::Commitment)>,
    pub(crate) r_base: (IP::LeftMessage, IP::RightMessage),
    // The fn() is here because PhantomData<T> is Sync iff T is Sync, and these types are not all
    // Sync
    _gipa: PhantomData<fn() -> GIPA<IP, IPC, D>>,
}

#[derive(Clone)]
pub struct GIPAAux<IP, IPC, D>
where
    D: Digest,
    IP: InnerProduct,
    IPC: IPCommitment<IP = IP>,
{
    pub(crate) r_transcript: Vec<Scalar<IPC>>,
    pub(crate) ck_base: FinalIPCommKey<IPC>,
    _gipa: PhantomData<GIPA<IP, IPC, D>>,
}

impl<IP, IPC, D> GIPA<IP, IPC, D>
where
    D: Digest,
    IP: InnerProduct,
    IPC: IPCommitment<IP = IP>,
{
    pub fn setup<'a>(size: usize, rng: impl Rng) -> Result<IPCommKey<'a, IPC>, Error> {
        IPC::setup(size, rng)
    }

    pub fn prove<'a>(
        ck: &IPCommKey<'a, IPC>,
        left: &[IP::LeftMessage],
        right: &[IP::RightMessage],
        output: &IP::Output,
        com: &IPC::Commitment,
    ) -> Result<GIPAProof<IP, IPC, D>, Error> {
        debug_assert_eq!(left.len(), right.len());
        debug_assert_eq!(&IP::inner_product(left, right)?, output);
        if !left.len().is_power_of_two() {
            // Power of 2 length
            return Err(Box::new(InnerProductArgumentError::MessageLengthInvalid(
                left.len(),
                right.len(),
            )));
        }
        if !IPC::verify(ck, left, right, &output, com)? {
            return Err(Box::new(InnerProductArgumentError::InnerProductInvalid));
        }

        let (proof, _) = Self::prove_with_aux(ck, left, right)?;
        Ok(proof)
    }

    pub fn verify<'a>(
        ck: &IPCommKey<'a, IPC>,
        com: &IPC::Commitment,
        proof: &GIPAProof<IP, IPC, D>,
    ) -> Result<bool, Error> {
        if !ck.ck_a.len().is_power_of_two() || ck.ck_a.len() != ck.ck_b.len() {
            // Power of 2 length
            return Err(Box::new(InnerProductArgumentError::MessageLengthInvalid(
                ck.ck_a.len(),
                ck.ck_b.len(),
            )));
        }
        // Calculate base commitment and transcript
        let (base_com, transcript) = Self::_compute_recursive_challenges(com, proof)?;
        // Calculate base commitment keys
        let ck_base = Self::_compute_final_commitment_keys(ck, &transcript)?;
        // Verify base commitment
        Self::_verify_base_commitment(&ck_base, &base_com, proof)
    }

    pub fn prove_with_aux(
        ck: &IPCommKey<IPC>,
        left: &[IP::LeftMessage],
        right: &[IP::RightMessage],
    ) -> Result<(GIPAProof<IP, IPC, D>, GIPAAux<IP, IPC, D>), Error> {
        Self::_prove(ck, left.to_vec(), right.to_vec())
    }

    // Returns vector of recursive commitments and transcripts in reverse order
    fn _prove<'a>(
        ck: &IPCommKey<'a, IPC>,
        mut left: Vec<IP::LeftMessage>,
        mut right: Vec<IP::RightMessage>,
    ) -> Result<(GIPAProof<IP, IPC, D>, GIPAAux<IP, IPC, D>), Error> {
        let mut ck = ck.clone();
        let mut r_commitment_steps = Vec::new();
        let mut r_transcript: Vec<Scalar<IPC>> = Vec::new();
        assert!(left.len().is_power_of_two());
        let (m_base, ck_base) = 'recurse: loop {
            let recurse = start_timer!(|| format!("Recurse round size {}", left.len()));
            if left.len() == 1 {
                // base case
                break 'recurse ((left[0].clone(), right[0].clone()), ck.to_owned());
            } else {
                // recursive step
                // Recurse with problem of half size
                let split = left.len() / 2;

                let left_1 = &left[split..];
                let left_2 = &left[..split];
                let right_1 = &right[..split];
                let right_2 = &right[split..];

                let (ck_l, ck_r) = ck.split(split);

                let cl = start_timer!(|| "Commit L");
                let com_1 = IPC::commit(&ck_l, &left_1, &right_1, || {
                    IP::inner_product(left_1, right_1).unwrap()
                })?;
                end_timer!(cl);
                let cr = start_timer!(|| "Commit R");
                let com_2 = IPC::commit(&ck_r, &left_2, &right_2, || {
                    IP::inner_product(left_2, right_2).unwrap()
                })?;
                end_timer!(cr);

                // Fiat-Shamir challenge
                let default_transcript = Default::default();
                let transcript = r_transcript.last().unwrap_or(&default_transcript);
                let (c, c_inv) = Self::compute_challenge(transcript, &com_1, &com_2)?;

                // Set up values for next step of recursion
                let rescale_m1 = start_timer!(|| "Rescale M1");
                left = cfg_iter!(left_1)
                    .map(|a| *a * c)
                    .zip(left_2)
                    .map(|(a_1, a_2)| a_1 + *a_2)
                    .collect::<Vec<IP::LeftMessage>>();
                end_timer!(rescale_m1);

                let rescale_m2 = start_timer!(|| "Rescale M2");
                right = cfg_iter!(right_2)
                    .map(|b| *b * c_inv)
                    .zip(right_1)
                    .map(|(b_1, b_2)| b_1 + *b_2)
                    .collect::<Vec<IP::RightMessage>>();
                end_timer!(rescale_m2);

                ck = IPCommKey::fold(&ck_l, &ck_r, &c_inv, &c)?;

                r_commitment_steps.push((com_1, com_2));
                r_transcript.push(c);
                end_timer!(recurse);
            }
        };
        r_transcript.reverse();
        r_commitment_steps.reverse();
        Ok((
            GIPAProof {
                r_commitment_steps,
                r_base: m_base,
                _gipa: PhantomData,
            },
            GIPAAux {
                r_transcript,
                ck_base: ck_base.try_into().unwrap(),
                _gipa: PhantomData,
            },
        ))
    }

    // Helper function used to calculate recursive challenges from proof execution (transcript in reverse)
    pub fn verify_recursive_challenge_transcript(
        com: &IPC::Commitment,
        proof: &GIPAProof<IP, IPC, D>,
    ) -> Result<(IPC::Commitment, Vec<Scalar<IPC>>), Error> {
        Self::_compute_recursive_challenges(com, proof)
    }

    fn _compute_recursive_challenges(
        com: &IPC::Commitment,
        proof: &GIPAProof<IP, IPC, D>,
    ) -> Result<(IPC::Commitment, Vec<Scalar<IPC>>), Error> {
        let mut com = com.clone();
        let mut r_transcript: Vec<Scalar<IPC>> = Vec::new();
        let default_transcript_entry = Default::default();
        for (com_1, com_2) in proof.r_commitment_steps.iter().rev() {
            // Fiat-Shamir challenge
            let transcript = r_transcript.last().unwrap_or(&default_transcript_entry);
            let (c, c_inv) = Self::compute_challenge(transcript, com_1, com_2)?;
            com = com + (com_1.clone() * c + com_2.clone() * c_inv);
            r_transcript.push(c);
        }
        r_transcript.reverse();
        Ok((com, r_transcript))
    }

    pub(crate) fn _compute_final_commitment_keys<'a>(
        ck: &IPCommKey<'a, IPC>,
        transcript: &[Scalar<IPC>],
    ) -> Result<IPCommKey<'a, IPC>, Error> {
        // Calculate base commitment keys
        assert!(ck.ck_a.len().is_power_of_two());

        let mut ck_a_agg_challenge_exponents = vec![Scalar::<IPC>::one()];
        let mut ck_b_agg_challenge_exponents = vec![Scalar::<IPC>::one()];
        for (i, c) in transcript.iter().enumerate() {
            let c_inv = c.inverse().unwrap();
            for j in 0..(2_usize).pow(i as u32) {
                ck_a_agg_challenge_exponents.push(ck_a_agg_challenge_exponents[j] * &c_inv);
                ck_b_agg_challenge_exponents.push(ck_b_agg_challenge_exponents[j] * c);
            }
        }
        assert_eq!(ck_a_agg_challenge_exponents.len(), ck.ck_a.len());
        let ck_a_base = IPC::left_key_msm(&ck.ck_a, &ck_a_agg_challenge_exponents)?;
        let ck_b_base = IPC::right_key_msm(&ck.ck_b, &ck_b_agg_challenge_exponents)?;

        Ok(IPCommKey::new(
            Cow::Owned(vec![ck_a_base]),
            Cow::Owned(vec![ck_b_base]),
            ck.ck_t.clone(),
        ))
    }

    pub(crate) fn _verify_base_commitment<'a>(
        base_ck: &IPCommKey<'a, IPC>,
        base_com: &IPC::Commitment,
        proof: &GIPAProof<IP, IPC, D>,
    ) -> Result<bool, Error> {
        let a_base = [proof.r_base.0.clone()];
        let b_base = [proof.r_base.1.clone()];
        let t_base = IP::inner_product(&a_base, &b_base)?;
        dbg!(&base_com);

        Ok(IPC::verify(&base_ck, &a_base, &b_base, &t_base, &base_com)?)
    }

    fn compute_challenge(
        transcript: &Scalar<IPC>,
        com_1: &IPC::Commitment,
        com_2: &IPC::Commitment,
    ) -> Result<(Scalar<IPC>, Scalar<IPC>), Error> {
        let mut counter_nonce = 0u32;
        let mut bytes = Vec::new();
        transcript.serialize_uncompressed(&mut bytes)?;
        com_1.serialize_uncompressed(&mut bytes)?;
        com_2.serialize_uncompressed(&mut bytes)?;
        let (c, c_inv) = 'challenge: loop {
            let mut hash_input = Vec::new();
            hash_input.extend_from_slice(&counter_nonce.to_be_bytes()[..]);
            hash_input.extend_from_slice(&bytes);

            let c = D::digest(&hash_input).as_slice()[..16].try_into().unwrap();
            let c: Scalar<IPC> = u128::from_be_bytes(c).into();
            if let Some(c_inv) = c.inverse() {
                // Optimization for multiexponentiation to rescale G2 elements with 128-bit challenge
                // Swap 'c' and 'c_inv' since can't control bit size of c_inv
                break 'challenge (c_inv, c);
            }
            counter_nonce += 1;
        };

        Ok((c, c_inv))
    }
}

impl<IP, IPC, D> Clone for GIPAProof<IP, IPC, D>
where
    D: Digest,
    IP: InnerProduct,
    IPC: IPCommitment<IP = IP>,
{
    fn clone(&self) -> Self {
        GIPAProof {
            r_commitment_steps: self.r_commitment_steps.clone(),
            r_base: self.r_base.clone(),
            _gipa: PhantomData,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ip_commitment::{
        mipp::MSMCommitment, pairing::PairingCommitment, scalar::ScalarCommitment,
    };

    use ark_bls12_381::Bls12_381;
    use ark_ec::pairing::Pairing;
    use ark_ff::UniformRand;
    use blake2::Blake2b;

    use ark_dh_commitments::random_generators;
    use ark_inner_products::{
        InnerProduct, MSMInnerProduct, PairingInnerProduct, ScalarInnerProduct,
    };

    const TEST_SIZE: usize = 8;

    #[test]
    fn pairing_inner_product_test() {
        type IP = PairingInnerProduct<Bls12_381>;
        type IPC = PairingCommitment<Bls12_381>;
        type PairingGIPA = GIPA<IP, IPC, Blake2b>;

        let mut rng = ark_std::test_rng();
        let ck = PairingGIPA::setup(TEST_SIZE, &mut rng).unwrap();
        let m_a = random_generators(&mut rng, TEST_SIZE);
        let m_b = random_generators(&mut rng, TEST_SIZE);
        let t = IP::inner_product(&m_a, &m_b).unwrap();
        let com = IPC::commit(&ck, &m_a, &m_b, || t).unwrap();

        let proof = PairingGIPA::prove(&ck, &m_a, &m_b, &t, &com).unwrap();

        assert!(PairingGIPA::verify(&ck, &com, &proof).unwrap());
    }

    #[test]
    fn multiexponentiation_inner_product_test() {
        type IP = MSMInnerProduct<<Bls12_381 as Pairing>::G1>;
        type IPC = MSMCommitment<Bls12_381>;
        type MultiExpGIPA = GIPA<IP, IPC, Blake2b>;

        let mut rng = ark_std::test_rng();
        let ck = MultiExpGIPA::setup(TEST_SIZE, &mut rng).unwrap();
        let m_a = random_generators(&mut rng, TEST_SIZE);
        let mut m_b = Vec::new();
        for _ in 0..TEST_SIZE {
            m_b.push(<Bls12_381 as Pairing>::ScalarField::rand(&mut rng));
        }
        let t = IP::inner_product(&m_a, &m_b).unwrap();
        let com = IPC::commit(&ck, &m_a, &m_b, || t).unwrap();

        let proof = MultiExpGIPA::prove(&ck, &m_a, &m_b, &t, &com).unwrap();

        assert!(MultiExpGIPA::verify(&ck, &com, &proof).unwrap());
    }

    #[test]
    fn scalar_inner_product_test() {
        type IP = ScalarInnerProduct<<Bls12_381 as Pairing>::ScalarField>;
        type IPC = ScalarCommitment<Bls12_381>;
        type ScalarGIPA = GIPA<IP, IPC, Blake2b>;

        let mut rng = ark_std::test_rng();
        let ck = ScalarGIPA::setup(TEST_SIZE, &mut rng).unwrap();
        let mut m_a = Vec::new();
        let mut m_b = Vec::new();
        for _ in 0..TEST_SIZE {
            m_a.push(<Bls12_381 as Pairing>::ScalarField::rand(&mut rng));
            m_b.push(<Bls12_381 as Pairing>::ScalarField::rand(&mut rng));
        }

        let t = IP::inner_product(&m_a, &m_b).unwrap();
        let com = IPC::commit(&ck, &m_a, &m_b, || t).unwrap();

        let proof = ScalarGIPA::prove(&ck, &m_a, &m_b, &t, &com).unwrap();

        assert!(ScalarGIPA::verify(&ck, &com, &proof).unwrap());
    }
}
