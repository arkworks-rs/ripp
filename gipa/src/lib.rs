use algebra::{bytes::ToBytes, fields::Field, to_bytes};
use digest::Digest;
use std::{
    error::Error as ErrorTrait,
    fmt::{Display, Formatter, Result as FmtResult},
    marker::PhantomData,
    ops::MulAssign,
};

use dh_commitments::DoublyHomomorphicCommitment;
use inner_products::InnerProduct;

pub type Error = Box<dyn ErrorTrait>;

pub trait InnerProductCommitmentArgument
where
    <Self::RMC as DoublyHomomorphicCommitment>::Message:
        MulAssign<<Self::LMC as DoublyHomomorphicCommitment>::Scalar>,
    <Self::IPC as DoublyHomomorphicCommitment>::Message:
        MulAssign<<Self::LMC as DoublyHomomorphicCommitment>::Scalar>,
    <Self::RMC as DoublyHomomorphicCommitment>::Key:
        MulAssign<<Self::LMC as DoublyHomomorphicCommitment>::Scalar>,
    <Self::IPC as DoublyHomomorphicCommitment>::Key:
        MulAssign<<Self::LMC as DoublyHomomorphicCommitment>::Scalar>,
    <Self::RMC as DoublyHomomorphicCommitment>::Output:
        MulAssign<<Self::LMC as DoublyHomomorphicCommitment>::Scalar>,
    <Self::IPC as DoublyHomomorphicCommitment>::Output:
        MulAssign<<Self::LMC as DoublyHomomorphicCommitment>::Scalar>,
{
    type IP: InnerProduct<
        LeftMessage = <Self::LMC as DoublyHomomorphicCommitment>::Message,
        RightMessage = <Self::RMC as DoublyHomomorphicCommitment>::Message,
        Output = <Self::IPC as DoublyHomomorphicCommitment>::Message,
    >;
    type LMC: DoublyHomomorphicCommitment;
    type RMC: DoublyHomomorphicCommitment<
        Scalar = <Self::LMC as DoublyHomomorphicCommitment>::Scalar,
    >;
    type IPC: DoublyHomomorphicCommitment<
        Scalar = <Self::LMC as DoublyHomomorphicCommitment>::Scalar,
    >;
    type Proof;

    fn prove(
        values: (
            &[<Self::IP as InnerProduct>::LeftMessage],
            &[<Self::IP as InnerProduct>::RightMessage],
            &<Self::IP as InnerProduct>::Output,
        ),
        ck: (
            &[<Self::LMC as DoublyHomomorphicCommitment>::Key],
            &[<Self::RMC as DoublyHomomorphicCommitment>::Key],
            &<Self::IPC as DoublyHomomorphicCommitment>::Key,
        ),
        com: (
            &<Self::LMC as DoublyHomomorphicCommitment>::Output,
            &<Self::RMC as DoublyHomomorphicCommitment>::Output,
            &<Self::IPC as DoublyHomomorphicCommitment>::Output,
        ),
    ) -> Result<Self::Proof, Error>;

    fn verify(
        ck: (
            &[<Self::LMC as DoublyHomomorphicCommitment>::Key],
            &[<Self::RMC as DoublyHomomorphicCommitment>::Key],
            &<Self::IPC as DoublyHomomorphicCommitment>::Key,
        ),
        com: (
            &<Self::LMC as DoublyHomomorphicCommitment>::Output,
            &<Self::RMC as DoublyHomomorphicCommitment>::Output,
            &<Self::IPC as DoublyHomomorphicCommitment>::Output,
        ),
        proof: &Self::Proof,
    ) -> Result<bool, Error>;
}

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
    r_commitment_steps:
        Vec<(
            (LMC::Output, RMC::Output, IPC::Output),
            (LMC::Output, RMC::Output, IPC::Output),
        )>,
    r_base: (LMC::Message, RMC::Message),
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
    r_transcript: Vec<LMC::Scalar>,
    _gipa: PhantomData<GIPA<IP, LMC, RMC, IPC, D>>,
}


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

    fn prove(
        values: (&[IP::LeftMessage], &[IP::RightMessage], &IP::Output),
        ck: (&[LMC::Key], &[RMC::Key], &IPC::Key),
        com: (&LMC::Output, &RMC::Output, &IPC::Output),
    ) -> Result<GIPAProof<IP, LMC, RMC, IPC, D>, Error> {
        if IP::inner_product(values.0, values.1)? != values.2.clone() {
            return Err(Box::new(InnerProductArgumentError::InnerProductInvalid));
        }
        if values.0.len().count_ones() != 1 {  // Power of 2 length
            return Err(Box::new(InnerProductArgumentError::MessageLengthInvalid(values.0.len(), values.1.len())));
        }
        if !(
            LMC::verify(ck.0, values.0, com.0)?
                && RMC::verify(ck.1, values.1, com.1)?
                && IPC::verify(&vec![ck.2.clone()], &vec![values.2.clone()], com.2)?
        ){
            return Err(Box::new(InnerProductArgumentError::InnerProductInvalid));
        }

        let (proof, _) = Self::prove(
            (values.0, values.1),
            (ck.0, ck.1, &vec![ck.2.clone()]),
        )?;
        Ok(proof)
    }


    fn verify(
        ck: (&[LMC::Key], &[RMC::Key], &IPC::Key),
        com: (&LMC::Output, &RMC::Output, &IPC::Output),
        proof: &GIPAProof<IP, LMC, RMC, IPC, D>,
    ) -> Result<bool, Error> {
        if ck.0.len().count_ones() != 1  || ck.0.len() != ck.1.len() {  // Power of 2 length
            return Err(Box::new(InnerProductArgumentError::MessageLengthInvalid(ck.0.len(), ck.1.len())));
        }

        Self::verify(
            (ck.0, ck.1, &vec![ck.2.clone()]),
            com,
            proof,
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
    pub fn prove(
        values: (&[IP::LeftMessage], &[IP::RightMessage]),
        ck: (&[LMC::Key], &[RMC::Key], &[IPC::Key]),
    ) -> Result<(GIPAProof<IP, LMC, RMC, IPC, D>, GIPAAux<IP, LMC, RMC, IPC, D>), Error> {
        let (mut r_commitment_steps, r_base, mut r_transcript) = Self::recursive_prove(
            values,
            ck,
            &Default::default(),
        )?;
        //r_commitment_steps.reverse();
        //r_transcript.reverse();

        Ok((
            GIPAProof { r_commitment_steps, r_base, _gipa: PhantomData },
            GIPAAux { r_transcript, _gipa: PhantomData },
         ))
    }

    pub fn verify(
        ck: (&[LMC::Key], &[RMC::Key], &[IPC::Key]),
        com: (&LMC::Output, &RMC::Output, &IPC::Output),
        proof: &GIPAProof<IP, LMC, RMC, IPC, D>,
    ) -> Result<bool, Error> {
        let mut clone= Clone::clone(proof);
        Self::recursive_verify(
            ck,
            com,
            &mut clone,
            &Default::default(),
        )
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
        ),
        Error,
    > {
        let (m_a, m_b) = values;
        let (ck_a, ck_b, ck_t) = ck;
        match m_a.len() {
            1 => Ok((Vec::new(), (m_a[0].clone(), m_b[0].clone()), Vec::new())), // base case
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
                let counter_nonce: usize = 0;
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

                let (mut r_steps, r_base, mut r_transcript) = Self::recursive_prove(
                    (&m_a_recurse, &m_b_recurse),
                    (&ck_a_recurse, &ck_b_recurse, ck_t),
                    &c,
                )?;

                r_steps.push((com_1, com_2));
                r_transcript.push(c);
                Ok((r_steps, r_base, r_transcript))
            },
            _ => unreachable!(), // If called only on message lengths power of 2
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
                Ok(
                    LMC::verify(ck_a, &a_base, com_a)?
                        && RMC::verify(ck_b, &b_base, com_b)?
                        && IPC::verify(ck_t, &t_base, com_t)?
                )
            }, // base case
            2..=usize::MAX if ck_a.len().count_ones() == 1 => { // recursive step
                // Fiat-Shamir challenge
                let (com_1, com_2) = proof.r_commitment_steps.pop().unwrap();
                let counter_nonce: usize = 0;
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
            },
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


//TODO: helper function for mul because relying on MulAssign
fn mul_helper<T: MulAssign<F> + Clone, F: Clone>(t: &T, f: &F) -> T {
    let mut clone = t.clone();
    clone.mul_assign(f.clone());
    clone
}

#[derive(Debug)]
pub enum InnerProductArgumentError {
    MessageLengthInvalid(usize, usize),
    InnerProductInvalid,
}

impl ErrorTrait for InnerProductArgumentError {
    fn source(self: &Self) -> Option<&(dyn ErrorTrait + 'static)> {
        None
    }
}

impl Display for InnerProductArgumentError {
    fn fmt(self: &Self, f: &mut Formatter<'_>) -> FmtResult {
        let msg = match self {
            InnerProductArgumentError::MessageLengthInvalid(left, right) => {
                format!("left length, right length: {}, {}", left, right)
            },
            InnerProductArgumentError::InnerProductInvalid => "inner product not sound".to_string(),
        };
        write!(f, "{}", msg)
    }
}
