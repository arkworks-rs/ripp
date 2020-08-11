use algebra::{bytes::ToBytes, fields::Field, to_bytes};
use digest::Digest;
use rand::Rng;
use std::{marker::PhantomData, ops::MulAssign};

use crate::{mul_helper, Error, InnerProductArgumentError, InnerProductCommitmentArgument};
use dh_commitments::DoublyHomomorphicCommitment;
use inner_products::InnerProduct;

pub struct GIPA<IP, LMC, RMC, IPC, D> {
    _inner_product: PhantomData<IP>,
    _left_commitment: PhantomData<LMC>,
    _right_commitment: PhantomData<RMC>,
    _inner_product_commitment: PhantomData<IPC>,
    _digest: PhantomData<D>,
}

pub struct GIPAProof<IP, LMC, RMC, IPC, D>
where
    D: Digest,
    IP: InnerProduct<
        LeftMessage = LMC::Message,
        RightMessage = RMC::Message,
        Output = IPC::Message,
    >,
    LMC: DoublyHomomorphicCommitment,
    RMC: DoublyHomomorphicCommitment<Scalar = LMC::Scalar>,
    IPC: DoublyHomomorphicCommitment<Scalar = LMC::Scalar>,
    RMC::Message: MulAssign<LMC::Scalar>,
    IPC::Message: MulAssign<LMC::Scalar>,
    RMC::Key: MulAssign<LMC::Scalar>,
    IPC::Key: MulAssign<LMC::Scalar>,
    RMC::Output: MulAssign<LMC::Scalar>,
    IPC::Output: MulAssign<LMC::Scalar>,
{
    pub(crate) r_commitment_steps: Vec<(
        (LMC::Output, RMC::Output, IPC::Output),
        (LMC::Output, RMC::Output, IPC::Output),
    )>,
    pub(crate) r_base: (LMC::Message, RMC::Message),
    _gipa: PhantomData<GIPA<IP, LMC, RMC, IPC, D>>,
}

#[derive(Clone)]
pub struct GIPAAux<IP, LMC, RMC, IPC, D>
where
    D: Digest,
    IP: InnerProduct<
        LeftMessage = LMC::Message,
        RightMessage = RMC::Message,
        Output = IPC::Message,
    >,
    LMC: DoublyHomomorphicCommitment,
    RMC: DoublyHomomorphicCommitment<Scalar = LMC::Scalar>,
    IPC: DoublyHomomorphicCommitment<Scalar = LMC::Scalar>,
    RMC::Message: MulAssign<LMC::Scalar>,
    IPC::Message: MulAssign<LMC::Scalar>,
    RMC::Key: MulAssign<LMC::Scalar>,
    IPC::Key: MulAssign<LMC::Scalar>,
    RMC::Output: MulAssign<LMC::Scalar>,
    IPC::Output: MulAssign<LMC::Scalar>,
{
    pub(crate) r_transcript: Vec<LMC::Scalar>,
    pub(crate) ck_base: (LMC::Key, RMC::Key),
    _gipa: PhantomData<GIPA<IP, LMC, RMC, IPC, D>>,
}

//TODO: Can extend GIPA to support "identity commitments" in addition to "compact commitments", i.e. for SIPP

impl<IP, LMC, RMC, IPC, D> InnerProductCommitmentArgument for GIPA<IP, LMC, RMC, IPC, D>
where
    D: Digest,
    IP: InnerProduct<
        LeftMessage = LMC::Message,
        RightMessage = RMC::Message,
        Output = IPC::Message,
    >,
    LMC: DoublyHomomorphicCommitment,
    RMC: DoublyHomomorphicCommitment<Scalar = LMC::Scalar>,
    IPC: DoublyHomomorphicCommitment<Scalar = LMC::Scalar>,
    RMC::Message: MulAssign<LMC::Scalar>,
    IPC::Message: MulAssign<LMC::Scalar>,
    RMC::Key: MulAssign<LMC::Scalar>,
    IPC::Key: MulAssign<LMC::Scalar>,
    RMC::Output: MulAssign<LMC::Scalar>,
    IPC::Output: MulAssign<LMC::Scalar>,
{
    type IP = IP;
    type LMC = LMC;
    type RMC = RMC;
    type IPC = IPC;
    type Proof = GIPAProof<IP, LMC, RMC, IPC, D>;
    type SetupOutput = (Vec<LMC::Key>, Vec<RMC::Key>, IPC::Key);

    fn setup<R: Rng>(
        rng: &mut R,
        size: usize,
    ) -> Result<(Vec<LMC::Key>, Vec<RMC::Key>, IPC::Key), Error> {
        Ok((
            LMC::setup(rng, size)?,
            RMC::setup(rng, size)?,
            IPC::setup(rng, 1)?.pop().unwrap(),
        ))
    }

    fn prove(
        values: (&[IP::LeftMessage], &[IP::RightMessage], &IP::Output),
        ck: (&[LMC::Key], &[RMC::Key], &IPC::Key),
        com: (&LMC::Output, &RMC::Output, &IPC::Output),
    ) -> Result<GIPAProof<IP, LMC, RMC, IPC, D>, Error> {
        if IP::inner_product(values.0, values.1)? != values.2.clone() {
            return Err(Box::new(InnerProductArgumentError::InnerProductInvalid));
        }
        if values.0.len().count_ones() != 1 {
            // Power of 2 length
            return Err(Box::new(InnerProductArgumentError::MessageLengthInvalid(
                values.0.len(),
                values.1.len(),
            )));
        }
        if !(LMC::verify(ck.0, values.0, com.0)?
            && RMC::verify(ck.1, values.1, com.1)?
            && IPC::verify(&vec![ck.2.clone()], &vec![values.2.clone()], com.2)?)
        {
            return Err(Box::new(InnerProductArgumentError::InnerProductInvalid));
        }

        let (proof, _) =
            Self::prove_with_aux((values.0, values.1), (ck.0, ck.1, &vec![ck.2.clone()]))?;
        Ok(proof)
    }

    fn verify(
        ck: (&[LMC::Key], &[RMC::Key], &IPC::Key),
        com: (&LMC::Output, &RMC::Output, &IPC::Output),
        proof: &GIPAProof<IP, LMC, RMC, IPC, D>,
    ) -> Result<bool, Error> {
        if ck.0.len().count_ones() != 1 || ck.0.len() != ck.1.len() {
            // Power of 2 length
            return Err(Box::new(InnerProductArgumentError::MessageLengthInvalid(
                ck.0.len(),
                ck.1.len(),
            )));
        }
        let mut clone = Clone::clone(proof);
        Self::recursive_verify(
            (ck.0, ck.1, &vec![ck.2.clone()]),
            com,
            &mut clone,
            &Default::default(),
        )
    }
}

impl<IP, LMC, RMC, IPC, D> GIPA<IP, LMC, RMC, IPC, D>
where
    D: Digest,
    IP: InnerProduct<
        LeftMessage = LMC::Message,
        RightMessage = RMC::Message,
        Output = IPC::Message,
    >,
    LMC: DoublyHomomorphicCommitment,
    RMC: DoublyHomomorphicCommitment<Scalar = LMC::Scalar>,
    IPC: DoublyHomomorphicCommitment<Scalar = LMC::Scalar>,
    RMC::Message: MulAssign<LMC::Scalar>,
    IPC::Message: MulAssign<LMC::Scalar>,
    RMC::Key: MulAssign<LMC::Scalar>,
    IPC::Key: MulAssign<LMC::Scalar>,
    RMC::Output: MulAssign<LMC::Scalar>,
    IPC::Output: MulAssign<LMC::Scalar>,
{
    pub fn prove_with_aux(
        values: (&[IP::LeftMessage], &[IP::RightMessage]),
        ck: (&[LMC::Key], &[RMC::Key], &[IPC::Key]),
    ) -> Result<
        (
            GIPAProof<IP, LMC, RMC, IPC, D>,
            GIPAAux<IP, LMC, RMC, IPC, D>,
        ),
        Error,
    > {
        let (r_commitment_steps, r_base, r_transcript, ck_base) =
            Self::recursive_prove(values, ck, &Default::default())?;

        Ok((
            GIPAProof {
                r_commitment_steps,
                r_base,
                _gipa: PhantomData,
            },
            GIPAAux {
                r_transcript,
                ck_base,
                _gipa: PhantomData,
            },
        ))
    }

    // Returns vector of recursive commitments and transcripts in reverse order
    fn recursive_prove(
        values: (&[IP::LeftMessage], &[IP::RightMessage]),
        ck: (&[LMC::Key], &[RMC::Key], &[IPC::Key]),
        transcript: &LMC::Scalar,
    ) -> Result<
        (
            Vec<(
                (LMC::Output, RMC::Output, IPC::Output),
                (LMC::Output, RMC::Output, IPC::Output),
            )>,
            (LMC::Message, RMC::Message),
            Vec<LMC::Scalar>,
            (LMC::Key, RMC::Key),
        ),
        Error,
    > {
        let (m_a, m_b) = values;
        let (ck_a, ck_b, ck_t) = ck;
        match m_a.len() {
            1 => Ok((
                Vec::new(),
                (m_a[0].clone(), m_b[0].clone()),
                Vec::new(),
                (ck_a[0].clone(), ck_b[0].clone()),
            )), // base case
            2..=usize::MAX if m_a.len().count_ones() == 1 => {
                // recursive step
                // Recurse with problem of half size
                let split = m_a.len() / 2;

                let m_a_1 = &m_a[split..];
                let m_a_2 = &m_a[..split];
                let ck_a_1 = &ck_a[..split];
                let ck_a_2 = &ck_a[split..];

                let m_b_1 = &m_b[..split];
                let m_b_2 = &m_b[split..];
                let ck_b_1 = &ck_b[split..];
                let ck_b_2 = &ck_b[..split];

                let com_1 = (
                    LMC::commit(ck_a_1, m_a_1)?,
                    RMC::commit(ck_b_1, m_b_1)?,
                    IPC::commit(ck_t, &vec![IP::inner_product(m_a_1, m_b_1)?])?,
                );
                let com_2 = (
                    LMC::commit(ck_a_2, m_a_2)?,
                    RMC::commit(ck_b_2, m_b_2)?,
                    IPC::commit(ck_t, &vec![IP::inner_product(m_a_2, m_b_2)?])?,
                );

                // Fiat-Shamir challenge
                let mut counter_nonce: usize = 0;
                let (c, c_inv) = loop {
                    let mut hash_input = Vec::new();
                    hash_input.extend_from_slice(&counter_nonce.to_be_bytes()[..]);
                    //TODO: Should use CanonicalSerialize instead of ToBytes
                    hash_input.extend_from_slice(&to_bytes![
                        transcript, com_1.0, com_1.1, com_1.2, com_2.0, com_2.1, com_2.2
                    ]?);
                    if let Some(c) = LMC::Scalar::from_random_bytes(&D::digest(&hash_input)) {
                        if let Some(c_inv) = c.inverse() {
                            break (c, c_inv);
                        }
                    };
                    counter_nonce += 1;
                };

                // Set up values for next step of recursion
                let m_a_recurse = m_a_1
                    .iter()
                    .map(|a| mul_helper(a, &c))
                    .zip(m_a_2)
                    .map(|(a_1, a_2)| a_1.clone() + a_2.clone())
                    .collect::<Vec<LMC::Message>>();

                let m_b_recurse = m_b_2
                    .iter()
                    .map(|b| mul_helper(b, &c_inv))
                    .zip(m_b_1)
                    .map(|(b_1, b_2)| b_1.clone() + b_2.clone())
                    .collect::<Vec<RMC::Message>>();

                let ck_a_recurse = ck_a_2
                    .iter()
                    .map(|a| mul_helper(a, &c_inv))
                    .zip(ck_a_1)
                    .map(|(a_1, a_2)| a_1.clone() + a_2.clone())
                    .collect::<Vec<LMC::Key>>();

                let ck_b_recurse = ck_b_1
                    .iter()
                    .map(|b| mul_helper(b, &c))
                    .zip(ck_b_2)
                    .map(|(b_1, b_2)| b_1.clone() + b_2.clone())
                    .collect::<Vec<RMC::Key>>();

                let (mut r_steps, r_base, mut r_transcript, ck_base) = Self::recursive_prove(
                    (&m_a_recurse, &m_b_recurse),
                    (&ck_a_recurse, &ck_b_recurse, ck_t),
                    &c,
                )?;

                r_steps.push((com_1, com_2));
                r_transcript.push(c);
                Ok((r_steps, r_base, r_transcript, ck_base))
            }
            _ => unreachable!(), // If called only on message lengths power of 2
        }
    }

    // Helper function used to calculate recursive challenges from proof execution (transcript in reverse)
    pub fn verify_recursive_challenge_transcript(
        com: (&LMC::Output, &RMC::Output, &IPC::Output),
        proof: &GIPAProof<IP, LMC, RMC, IPC, D>,
    ) -> Result<((LMC::Output, RMC::Output, IPC::Output), Vec<LMC::Scalar>), Error> {
        let mut clone = Clone::clone(proof);
        Self::_verify_recursive_challenges(com, &mut clone, &Default::default())
    }

    fn _verify_recursive_challenges(
        com: (&LMC::Output, &RMC::Output, &IPC::Output),
        proof: &mut GIPAProof<IP, LMC, RMC, IPC, D>,
        transcript: &LMC::Scalar,
    ) -> Result<((LMC::Output, RMC::Output, IPC::Output), Vec<LMC::Scalar>), Error> {
        let (com_a, com_b, com_t) = com;
        match proof.r_commitment_steps.len() {
            0 => Ok(((com_a.clone(), com_b.clone(), com_t.clone()), Vec::new())), // base case
            _ => {
                // recursive step
                // Fiat-Shamir challenge
                let (com_1, com_2) = proof.r_commitment_steps.pop().unwrap();
                let mut counter_nonce: usize = 0;
                let (c, c_inv) = loop {
                    let mut hash_input = Vec::new();
                    hash_input.extend_from_slice(&counter_nonce.to_be_bytes()[..]);
                    hash_input.extend_from_slice(&to_bytes![
                        transcript, com_1.0, com_1.1, com_1.2, com_2.0, com_2.1, com_2.2
                    ]?);
                    if let Some(c) = LMC::Scalar::from_random_bytes(&D::digest(&hash_input)) {
                        if let Some(c_inv) = c.inverse() {
                            break (c, c_inv);
                        }
                    };
                    counter_nonce += 1;
                };

                let com_recurse = (
                    mul_helper(&com_1.0, &c) + com_a.clone() + mul_helper(&com_2.0, &c_inv),
                    mul_helper(&com_1.1, &c) + com_b.clone() + mul_helper(&com_2.1, &c_inv),
                    mul_helper(&com_1.2, &c) + com_t.clone() + mul_helper(&com_2.2, &c_inv),
                );

                let (base_com, mut challenges) = Self::_verify_recursive_challenges(
                    (&com_recurse.0, &com_recurse.1, &com_recurse.2),
                    proof,
                    &c,
                )?;
                challenges.push(c);
                Ok((base_com, challenges))
            }
        }
    }

    fn recursive_verify(
        ck: (&[LMC::Key], &[RMC::Key], &[IPC::Key]),
        com: (&LMC::Output, &RMC::Output, &IPC::Output),
        proof: &mut GIPAProof<IP, LMC, RMC, IPC, D>,
        transcript: &LMC::Scalar,
    ) -> Result<bool, Error> {
        let (ck_a, ck_b, ck_t) = ck;
        let (com_a, com_b, com_t) = com;
        match ck_a.len() {
            1 => {
                let a_base = vec![proof.r_base.0.clone()];
                let b_base = vec![proof.r_base.1.clone()];
                let t_base = vec![IP::inner_product(&a_base, &b_base)?];
                Ok(LMC::verify(ck_a, &a_base, com_a)?
                    && RMC::verify(ck_b, &b_base, com_b)?
                    && IPC::verify(ck_t, &t_base, com_t)?)
            } // base case
            2..=usize::MAX if ck_a.len().count_ones() == 1 => {
                // recursive step
                // Fiat-Shamir challenge
                let (com_1, com_2) = proof.r_commitment_steps.pop().unwrap();
                let mut counter_nonce: usize = 0;
                let (c, c_inv) = loop {
                    let mut hash_input = Vec::new();
                    hash_input.extend_from_slice(&counter_nonce.to_be_bytes()[..]);
                    hash_input.extend_from_slice(&to_bytes![
                        transcript, com_1.0, com_1.1, com_1.2, com_2.0, com_2.1, com_2.2
                    ]?);
                    if let Some(c) = LMC::Scalar::from_random_bytes(&D::digest(&hash_input)) {
                        if let Some(c_inv) = c.inverse() {
                            break (c, c_inv);
                        }
                    };
                    counter_nonce += 1;
                };

                let split = ck_a.len() / 2;
                let ck_a_1 = &ck_a[..split];
                let ck_a_2 = &ck_a[split..];
                let ck_b_1 = &ck_b[split..];
                let ck_b_2 = &ck_b[..split];

                let ck_a_recurse = ck_a_2
                    .iter()
                    .map(|a| mul_helper(a, &c_inv))
                    .zip(ck_a_1)
                    .map(|(a_1, a_2)| a_1.clone() + a_2.clone())
                    .collect::<Vec<LMC::Key>>();

                let ck_b_recurse = ck_b_1
                    .iter()
                    .map(|b| mul_helper(b, &c))
                    .zip(ck_b_2)
                    .map(|(b_1, b_2)| b_1.clone() + b_2.clone())
                    .collect::<Vec<RMC::Key>>();

                let com_recurse = (
                    mul_helper(&com_1.0, &c) + com_a.clone() + mul_helper(&com_2.0, &c_inv),
                    mul_helper(&com_1.1, &c) + com_b.clone() + mul_helper(&com_2.1, &c_inv),
                    mul_helper(&com_1.2, &c) + com_t.clone() + mul_helper(&com_2.2, &c_inv),
                );

                Self::recursive_verify(
                    (&ck_a_recurse, &ck_b_recurse, ck_t),
                    (&com_recurse.0, &com_recurse.1, &com_recurse.2),
                    proof,
                    &c,
                )
            }
            _ => unreachable!(), // If called only on message lengths power of 2
        }
    }
}

impl<IP, LMC, RMC, IPC, D> Clone for GIPAProof<IP, LMC, RMC, IPC, D>
where
    D: Digest,
    IP: InnerProduct<
        LeftMessage = LMC::Message,
        RightMessage = RMC::Message,
        Output = IPC::Message,
    >,
    LMC: DoublyHomomorphicCommitment,
    RMC: DoublyHomomorphicCommitment<Scalar = LMC::Scalar>,
    IPC: DoublyHomomorphicCommitment<Scalar = LMC::Scalar>,
    RMC::Message: MulAssign<LMC::Scalar>,
    IPC::Message: MulAssign<LMC::Scalar>,
    RMC::Key: MulAssign<LMC::Scalar>,
    IPC::Key: MulAssign<LMC::Scalar>,
    RMC::Output: MulAssign<LMC::Scalar>,
    IPC::Output: MulAssign<LMC::Scalar>,
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
    use algebra::{bls12_381::Bls12_381, curves::PairingEngine, UniformRand};
    use blake2::Blake2b;
    use rand::{rngs::StdRng, SeedableRng};

    use dh_commitments::{
        afgho16::{AFGHOCommitmentG1, AFGHOCommitmentG2},
        identity::IdentityCommitment,
        pedersen::PedersenCommitment,
        random_generators,
    };
    use inner_products::{
        ExtensionFieldElement, InnerProduct, MultiexponentiationInnerProduct, PairingInnerProduct,
        ScalarInnerProduct,
    };

    type GC1 = AFGHOCommitmentG1<Bls12_381>;
    type GC2 = AFGHOCommitmentG2<Bls12_381>;
    type SC1 = PedersenCommitment<<Bls12_381 as PairingEngine>::G1Projective>;
    type SC2 = PedersenCommitment<<Bls12_381 as PairingEngine>::G2Projective>;
    const TEST_SIZE: usize = 8;

    #[test]
    fn pairing_inner_product_test() {
        type IP = PairingInnerProduct<Bls12_381>;
        type IPC =
            IdentityCommitment<ExtensionFieldElement<Bls12_381>, <Bls12_381 as PairingEngine>::Fr>;
        type PairingGIPA = GIPA<IP, GC1, GC2, IPC, Blake2b>;

        let mut rng = StdRng::seed_from_u64(0u64);
        let (ck_a, ck_b, ck_t) = PairingGIPA::setup(&mut rng, TEST_SIZE).unwrap();
        let m_a = random_generators(&mut rng, TEST_SIZE);
        let m_b = random_generators(&mut rng, TEST_SIZE);
        let com_a = GC1::commit(&ck_a, &m_a).unwrap();
        let com_b = GC2::commit(&ck_b, &m_b).unwrap();
        let t = vec![IP::inner_product(&m_a, &m_b).unwrap()];
        let com_t = IPC::commit(&vec![ck_t.clone()], &t).unwrap();

        let proof = PairingGIPA::prove(
            (&m_a, &m_b, &t[0]),
            (&ck_a, &ck_b, &ck_t),
            (&com_a, &com_b, &com_t),
        )
        .unwrap();

        assert!(
            PairingGIPA::verify((&ck_a, &ck_b, &ck_t), (&com_a, &com_b, &com_t), &proof,).unwrap()
        );
    }

    #[test]
    fn multiexponentiation_inner_product_test() {
        type IP = MultiexponentiationInnerProduct<<Bls12_381 as PairingEngine>::G1Projective>;
        type IPC = IdentityCommitment<
            <Bls12_381 as PairingEngine>::G1Projective,
            <Bls12_381 as PairingEngine>::Fr,
        >;
        type MultiExpGIPA = GIPA<IP, GC1, SC1, IPC, Blake2b>;

        let mut rng = StdRng::seed_from_u64(0u64);
        let (ck_a, ck_b, ck_t) = MultiExpGIPA::setup(&mut rng, TEST_SIZE).unwrap();
        let m_a = random_generators(&mut rng, TEST_SIZE);
        let mut m_b = Vec::new();
        for _ in 0..TEST_SIZE {
            m_b.push(<Bls12_381 as PairingEngine>::Fr::rand(&mut rng));
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
            MultiExpGIPA::verify((&ck_a, &ck_b, &ck_t), (&com_a, &com_b, &com_t), &proof,).unwrap()
        );
    }

    #[test]
    fn scalar_inner_product_test() {
        type IP = ScalarInnerProduct<<Bls12_381 as PairingEngine>::Fr>;
        type IPC =
            IdentityCommitment<<Bls12_381 as PairingEngine>::Fr, <Bls12_381 as PairingEngine>::Fr>;
        type ScalarGIPA = GIPA<IP, SC2, SC2, IPC, Blake2b>;

        let mut rng = StdRng::seed_from_u64(0u64);
        let (ck_a, ck_b, ck_t) = ScalarGIPA::setup(&mut rng, TEST_SIZE).unwrap();
        let mut m_a = Vec::new();
        let mut m_b = Vec::new();
        for _ in 0..TEST_SIZE {
            m_a.push(<Bls12_381 as PairingEngine>::Fr::rand(&mut rng));
            m_b.push(<Bls12_381 as PairingEngine>::Fr::rand(&mut rng));
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
            ScalarGIPA::verify((&ck_a, &ck_b, &ck_t), (&com_a, &com_b, &com_t), &proof,).unwrap()
        );
    }
}
