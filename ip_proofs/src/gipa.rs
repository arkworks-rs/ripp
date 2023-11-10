use ark_ff::{Field, One};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::{
    borrow::Cow, cfg_iter, convert::TryInto, end_timer, marker::PhantomData, rand::Rng, start_timer,
};
use digest::Digest;

use crate::{
    ip_commitment::{FinalIPCommKey, IPCommKey, IPCommitment, Scalar},
    mul_helper, Error, InnerProductArgumentError,
};
use ark_inner_products::InnerProduct;

#[cfg(feature = "parallel")]
use rayon::prelude::*;

pub struct GIPA<IP, IPC, D> {
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

//TODO: Can extend GIPA to support "identity commitments" in addition to "compact commitments", i.e. for SIPP

impl<IP, IPC, D> GIPA<IP, IPC, D>
where
    D: Digest,
    IP: InnerProduct,
    IPC: IPCommitment<IP = IP>,
{
    pub fn setup<'a>(size: usize, rng: &mut impl Rng) -> Result<IPCommKey<'a, IPC>, Error> {
        IPC::setup(size, rng)
    }

    pub fn prove<'a>(
        ck: &IPCommKey<'a, IPC>,
        left: &[IP::LeftMessage],
        right: &[IP::RightMessage],
        output: &IP::Output,
        com: &IPC::Commitment,
    ) -> Result<GIPAProof<IP, IPC, D>, Error> {
        if &IP::inner_product(left, right)? != output {
            return Err(Box::new(InnerProductArgumentError::InnerProductInvalid));
        }
        if left.len().count_ones() != 1 {
            // Power of 2 length
            return Err(Box::new(InnerProductArgumentError::MessageLengthInvalid(
                left.len(),
                right.len(),
            )));
        }
        if !IPC::verify(ck, left, right, core::slice::from_ref(output), com)? {
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
                break 'recurse ((left[0].clone(), right[0].clone()), ck);
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
                let com_1 = IPC::commit(
                    &ck_l,
                    &left_1,
                    &right_1,
                    &[IP::inner_product(left_1, right_1)?],
                )?;
                end_timer!(cl);
                let cr = start_timer!(|| "Commit R");
                let com_2 = IPC::commit(
                    &ck_r,
                    &left_2,
                    &right_2,
                    &[IP::inner_product(left_2, right_2)?],
                )?;
                end_timer!(cr);

                // Fiat-Shamir challenge
                let mut counter_nonce: usize = 0;
                let default_transcript = Default::default();
                let transcript = r_transcript.last().unwrap_or(&default_transcript);
                let (c, c_inv) = 'challenge: loop {
                    let mut hash_input = Vec::new();
                    hash_input.extend_from_slice(&counter_nonce.to_be_bytes()[..]);
                    transcript.serialize_uncompressed(&mut hash_input)?;
                    com_1.serialize_uncompressed(&mut hash_input)?;
                    com_2.serialize_uncompressed(&mut hash_input)?;
                    let c: Scalar<IPC> = u128::from_be_bytes(
                        D::digest(&hash_input).as_slice()[0..16].try_into().unwrap(),
                    )
                    .into();
                    if let Some(c_inv) = c.inverse() {
                        // Optimization for multiexponentiation to rescale G2 elements with 128-bit challenge
                        // Swap 'c' and 'c_inv' since can't control bit size of c_inv
                        break 'challenge (c_inv, c);
                    }
                    counter_nonce += 1;
                };

                // Set up values for next step of recursion
                let rescale_m1 = start_timer!(|| "Rescale M1");
                left = cfg_iter!(left_1)
                    .map(|a| mul_helper(a, &c))
                    .zip(left_2)
                    .map(|(a_1, a_2)| a_1 + a_2.clone())
                    .collect::<Vec<IP::LeftMessage>>();
                end_timer!(rescale_m1);

                let rescale_m2 = start_timer!(|| "Rescale M2");
                right = cfg_iter!(right_2)
                    .map(|b| mul_helper(b, &c_inv))
                    .zip(right_1)
                    .map(|(b_1, b_2)| b_1 + b_2.clone())
                    .collect::<Vec<IP::RightMessage>>();
                end_timer!(rescale_m2);

                ck.fold_into(&ck_l, &ck_r, &c_inv, &c)?;

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
        for (com_1, com_2) in proof.r_commitment_steps.iter().rev() {
            // Fiat-Shamir challenge
            let mut counter_nonce: usize = 0;
            let default_transcript = Default::default();
            let transcript = r_transcript.last().unwrap_or(&default_transcript);
            let (c, c_inv) = 'challenge: loop {
                let mut hash_input = Vec::new();
                hash_input.extend_from_slice(&counter_nonce.to_be_bytes()[..]);
                transcript.serialize_uncompressed(&mut hash_input)?;
                com_1.serialize_uncompressed(&mut hash_input)?;
                com_2.serialize_uncompressed(&mut hash_input)?;
                let c: Scalar<IPC> = u128::from_be_bytes(
                    D::digest(&hash_input).as_slice()[0..16].try_into().unwrap(),
                )
                .into();
                if let Some(c_inv) = c.inverse() {
                    // Optimization for multiexponentiation to rescale G2 elements with 128-bit challenge
                    // Swap 'c' and 'c_inv' since can't control bit size of c_inv
                    break 'challenge (c_inv, c);
                }
                counter_nonce += 1;
            };

            com += com_1 * c + com_2 * c_inv;
            r_transcript.push(c);
        }
        r_transcript.reverse();
        Ok((com, r_transcript))
    }

    pub(crate) fn _compute_final_commitment_keys<'a>(
        ck: &IPCommKey<'a, IPC>,
        transcript: &Vec<Scalar<IPC>>,
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
        //TODO: Optimization: Use VariableMSM multiexponentiation
        let ck_a_base_init = mul_helper(&ck.ck_a[0], &ck_a_agg_challenge_exponents[0]);
        let ck_a_base = ck.ck_a[1..]
            .iter()
            .zip(&ck_a_agg_challenge_exponents[1..])
            .map(|(g, x)| mul_helper(g, &x))
            .fold(ck_a_base_init, |sum, x| sum + x);
        //.reduce(|| ck_a_base_init.clone(), |sum, x| sum + x);
        let ck_b_base_init = mul_helper(&ck.ck_b[0], &ck_b_agg_challenge_exponents[0]);
        let ck_b_base = ck.ck_b[1..]
            .iter()
            .zip(&ck_b_agg_challenge_exponents[1..])
            .map(|(g, x)| mul_helper(g, &x))
            .fold(ck_b_base_init, |sum, x| sum + x);
        //.reduce(|| ck_b_base_init.clone(), |sum, x| sum + x);
        Ok(IPCommKey::new(
            Cow::Owned(vec![ck_a_base]),
            Cow::Owned(vec![ck_b_base]),
            Cow::Owned(ck.ck_t.as_ref().to_vec()),
        ))
    }

    pub(crate) fn _verify_base_commitment<'a>(
        base_ck: &IPCommKey<'a, IPC>,
        base_com: &IPC::Commitment,
        proof: &GIPAProof<IP, IPC, D>,
    ) -> Result<bool, Error> {
        let a_base = [proof.r_base.0.clone()];
        let b_base = [proof.r_base.1.clone()];
        let t_base = [IP::inner_product(&a_base, &b_base)?];

        Ok(IPC::verify(&base_ck, &a_base, &b_base, &t_base, &base_com)?)
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
    use crate::ip_commitment::IdentityCommitment;

    use ark_bls12_381::Bls12_381;
    use ark_ec::pairing::{Pairing, PairingOutput};
    use ark_ff::UniformRand;
    use ark_std::rand::{rngs::StdRng, SeedableRng};
    use blake2::Blake2b;

    use ark_dh_commitments::{
        afgho16::{AFGHOCommitmentG1, AFGHOCommitmentG2},
        pedersen::PedersenCommitment,
        random_generators,
    };
    use ark_inner_products::{
        InnerProduct, MultiexponentiationInnerProduct, PairingInnerProduct, ScalarInnerProduct,
    };

    type GC1 = AFGHOCommitmentG1<Bls12_381>;
    type GC2 = AFGHOCommitmentG2<Bls12_381>;
    type SC1 = PedersenCommitment<<Bls12_381 as Pairing>::G1>;
    type SC2 = PedersenCommitment<<Bls12_381 as Pairing>::G2>;

    const TEST_SIZE: usize = 8;

    #[test]
    fn pairing_inner_product_test() {
        type IP = PairingInnerProduct<Bls12_381>;
        type IPC = IdentityCommitment<IP>;
        type PairingGIPA = GIPA<IP, IPC, Blake2b>;

        let mut rng = StdRng::seed_from_u64(0u64);
        let ck = PairingGIPA::setup(TEST_SIZE, &mut rng).unwrap();
        let m_a = random_generators(&mut rng, TEST_SIZE);
        let m_b = random_generators(&mut rng, TEST_SIZE);
        let com_a = GC1::commit(&ck_a, &m_a).unwrap();
        let com_b = GC2::commit(&ck_b, &m_b).unwrap();
        let t = vec![IP::inner_product(&m_a, &m_b).unwrap()];
        let com_t = IPC::commit(&vec![ck_t.clone()], &t).unwrap();

        let proof = PairingGIPA::prove(
            (&ck_a, &ck_b, &ck_t),
            (&m_a, &m_b, &t[0]),
            (&com_a, &com_b, &com_t),
        )
        .unwrap();

        assert!(
            PairingGIPA::verify((&ck_a, &ck_b, &ck_t), (&com_a, &com_b, &com_t), &proof).unwrap()
        );
    }

    #[test]
    fn multiexponentiation_inner_product_test() {
        type IP = MultiexponentiationInnerProduct<<Bls12_381 as Pairing>::G1>;
        type IPC =
            IdentityCommitment<<Bls12_381 as Pairing>::G1, <Bls12_381 as Pairing>::ScalarField>;
        type MultiExpGIPA = GIPA<IP, IPC, Blake2b>;

        let mut rng = StdRng::seed_from_u64(0u64);
        let (ck_a, ck_b, ck_t) = MultiExpGIPA::setup(&mut rng, TEST_SIZE).unwrap();
        let m_a = random_generators(&mut rng, TEST_SIZE);
        let mut m_b = Vec::new();
        for _ in 0..TEST_SIZE {
            m_b.push(<Bls12_381 as Pairing>::ScalarField::rand(&mut rng));
        }
        let com_a = GC1::commit(&ck_a, &m_a).unwrap();
        let com_b = SC1::commit(&ck_b, &m_b).unwrap();
        let t = vec![IP::inner_product(&m_a, &m_b).unwrap()];
        let com_t = IPC::commit(&vec![ck_t.clone()], &t).unwrap();

        let proof = MultiExpGIPA::prove(
            (&m_a, &m_b, &t[0]),
            (&ck_a, &ck_b, &ck_t),
            (&com_a, &com_b, &com_t),
        )
        .unwrap();

        assert!(
            MultiExpGIPA::verify((&ck_a, &ck_b, &ck_t), (&com_a, &com_b, &com_t), &proof).unwrap()
        );
    }

    #[test]
    fn scalar_inner_product_test() {
        type IP = ScalarInnerProduct<<Bls12_381 as Pairing>::ScalarField>;
        type IPC = IdentityCommitment<
            <Bls12_381 as Pairing>::ScalarField,
            <Bls12_381 as Pairing>::ScalarField,
        >;
        type ScalarGIPA = GIPA<IP, IPC, Blake2b>;

        let mut rng = StdRng::seed_from_u64(0u64);
        let (ck_a, ck_b, ck_t) = ScalarGIPA::setup(&mut rng, TEST_SIZE).unwrap();
        let mut m_a = Vec::new();
        let mut m_b = Vec::new();
        for _ in 0..TEST_SIZE {
            m_a.push(<Bls12_381 as Pairing>::ScalarField::rand(&mut rng));
            m_b.push(<Bls12_381 as Pairing>::ScalarField::rand(&mut rng));
        }
        let com_a = SC2::commit(&ck_a, &m_a).unwrap();
        let com_b = SC2::commit(&ck_b, &m_b).unwrap();
        let t = vec![IP::inner_product(&m_a, &m_b).unwrap()];
        let com_t = IPC::commit(&vec![ck_t.clone()], &t).unwrap();

        let proof = ScalarGIPA::prove(
            (&m_a, &m_b, &t[0]),
            (&ck_a, &ck_b, &ck_t),
            (&com_a, &com_b, &com_t),
        )
        .unwrap();

        assert!(
            ScalarGIPA::verify((&ck_a, &ck_b, &ck_t), (&com_a, &com_b, &com_t), &proof).unwrap()
        );
    }
}
