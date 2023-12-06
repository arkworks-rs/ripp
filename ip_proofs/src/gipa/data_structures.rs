use ark_inner_products::InnerProduct;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, Valid};
use derivative::Derivative;

use crate::ip_commitment::{Commitment, FinalIPCommKey, IPCommKey, IPCommitment, Scalar};

#[derive(Derivative, CanonicalSerialize, CanonicalDeserialize)]
#[derivative(Clone, Debug, PartialEq, Eq)]
/// A "Twisted" GIPA instance
pub struct Instance<IPC: IPCommitment> {
    /// Size of input vectors.
    pub size: usize,
    /// The output of the inner product.
    pub output: <IPC::IP as InnerProduct>::Output,
    /// The commitment to the foregoing.
    pub commitment: Commitment<IPC>,
    /// The challenge used for performing the twist.
    /// Or equivalently, for taking the random linear combination
    pub twist: Scalar<IPC>,
}

#[derive(Derivative, CanonicalSerialize, CanonicalDeserialize)]
#[derivative(Clone, Debug, PartialEq, Eq)]
/// A GIPA witness
pub struct Witness<IP: InnerProduct> {
    /// The left input vector.
    pub left: Vec<IP::LeftMessage>,
    /// The right input vector.
    pub right: Vec<IP::RightMessage>,
}

#[derive(CanonicalSerialize, CanonicalDeserialize, Derivative)]
#[derivative(
    Clone(bound = "IPC: IPCommitment<IP = IP>, IP: InnerProduct"),
    Debug(bound = "IPC: IPCommitment<IP = IP>, IP: InnerProduct"),
    PartialEq(bound = "IPC: IPCommitment<IP = IP>, IP: InnerProduct"),
    Eq(bound = "IPC: IPCommitment<IP = IP>, IP: InnerProduct")
)]
pub struct Proof<IP, IPC>
where
    IP: InnerProduct,
    IPC: IPCommitment<IP = IP>,
{
    pub(crate) commitments: Vec<(Commitment<IPC>, Commitment<IPC>)>,
    pub(crate) final_msg: (IP::LeftMessage, IP::RightMessage),
}

impl<IP, IPC> Proof<IP, IPC>
where
    IP: InnerProduct,
    IPC: IPCommitment<IP = IP>,
{
    pub fn new(
        commitments: Vec<(Commitment<IPC>, Commitment<IPC>)>,
        final_msg: (IP::LeftMessage, IP::RightMessage),
    ) -> Self {
        Self {
            commitments,
            final_msg,
        }
    }
}

#[derive(Clone)]
pub struct Aux<IPC: IPCommitment> {
    pub(crate) challenges: Vec<Scalar<IPC>>,
    pub(crate) final_ck: FinalIPCommKey<IPC>,
}

impl<IPC: IPCommitment> Aux<IPC> {
    pub fn new(challenges: Vec<Scalar<IPC>>, final_ck: FinalIPCommKey<IPC>) -> Self {
        Self {
            challenges,
            final_ck,
        }
    }
}

#[derive(Derivative, CanonicalSerialize)]
#[derivative(Clone(bound = "IPC: IPCommitment"), Debug(bound = "IPC: IPCommitment"))]
pub struct ProverKey<'b, IPC: IPCommitment> {
    pub ck: IPCommKey<'b, IPC>,
}

impl<'b, IPC: IPCommitment> Valid for ProverKey<'b, IPC> {
    fn check(&self) -> Result<(), ark_serialize::SerializationError> {
        self.ck.check()
    }
}

impl<'b, IPC: IPCommitment> CanonicalDeserialize for ProverKey<'b, IPC> {
    fn deserialize_with_mode<R: std::io::prelude::Read>(
        reader: R,
        compress: ark_serialize::Compress,
        validate: ark_serialize::Validate,
    ) -> Result<Self, ark_serialize::SerializationError> {
        let ck = IPCommKey::deserialize_with_mode(reader, compress, validate)?;
        Ok(Self { ck })
    }
}

impl<'a, IPC: IPCommitment> ProverKey<'_, IPC> {
    pub fn vk(&self) -> VerifierKey<'_, IPC> {
        VerifierKey {
            ck: self.ck.clone(),
        }
    }
}

#[derive(Derivative, CanonicalSerialize)]
#[derivative(Clone(bound = "IPC: IPCommitment"), Debug(bound = "IPC: IPCommitment"))]
pub struct VerifierKey<'b, IPC: IPCommitment> {
    pub ck: IPCommKey<'b, IPC>,
}

impl<'b, IPC: IPCommitment> Valid for VerifierKey<'b, IPC> {
    fn check(&self) -> Result<(), ark_serialize::SerializationError> {
        self.ck.check()
    }
}

impl<'b, IPC: IPCommitment> CanonicalDeserialize for VerifierKey<'b, IPC> {
    fn deserialize_with_mode<R: std::io::prelude::Read>(
        reader: R,
        compress: ark_serialize::Compress,
        validate: ark_serialize::Validate,
    ) -> Result<Self, ark_serialize::SerializationError> {
        let ck = IPCommKey::deserialize_with_mode(reader, compress, validate)?;
        Ok(Self { ck })
    }
}
