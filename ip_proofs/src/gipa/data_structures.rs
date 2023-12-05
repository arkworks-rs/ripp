use ark_inner_products::InnerProduct;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
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

#[derive(Derivative)]
#[derivative(Clone(bound = "IPC: IPCommitment"), Debug(bound = "IPC: IPCommitment"))]
pub struct ProverKey<'a, IPC: IPCommitment> {
    pub ck: IPCommKey<'a, IPC>,
}

#[derive(Derivative)]
#[derivative(Clone(bound = "IPC: IPCommitment"), Debug(bound = "IPC: IPCommitment"))]
pub struct VerifierKey<'a, IPC: IPCommitment> {
    pub ck: IPCommKey<'a, IPC>,
}
