use ark_std::marker::PhantomData;

use ark_inner_products::InnerProduct;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use derivative::Derivative;
use digest::Digest;

use super::GIPA;
use crate::ip_commitment::{FinalIPCommKey, IPCommitment, Scalar};

#[derive(Derivative, CanonicalSerialize, CanonicalDeserialize)]
#[derivative(Clone, Debug, PartialEq, Eq)]
/// A "Twisted" GIPA instance
pub struct Instance<IPC: IPCommitment> {
    /// Size of input vectors.
    pub size: usize,
    /// The output of the inner product.
    pub output: <IPC::IP as InnerProduct>::Output,
    /// The commitment to the foregoing.
    pub commitment: IPC::Commitment,
    /// The challenge used for performing the twist.
    pub random_challenge: Scalar<IPC>,
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
#[derivative(Clone, Debug, PartialEq, Eq)]
pub struct Proof<IP, IPC, D>
where
    D: Digest,
    IP: InnerProduct,
    IPC: IPCommitment<IP = IP>,
{
    pub(crate) r_commitment_steps: Vec<(IPC::Commitment, IPC::Commitment)>,
    pub(crate) r_base: (IP::LeftMessage, IP::RightMessage),
    // The fn() is here because PhantomData<T>
    // is Sync iff T is Sync, and these types are not all Sync
    #[derivative(Debug = "ignore")]
    _gipa: PhantomData<fn() -> (IP, IPC, D)>,
}

impl<IP, IPC, D> Proof<IP, IPC, D>
where
    D: Digest,
    IP: InnerProduct,
    IPC: IPCommitment<IP = IP>,
{
    pub fn new(
        r_commitment_steps: Vec<(IPC::Commitment, IPC::Commitment)>,
        r_base: (IP::LeftMessage, IP::RightMessage),
    ) -> Self {
        Self {
            r_commitment_steps,
            r_base,
            _gipa: PhantomData,
        }
    }
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

impl<IP, IPC, D> GIPAAux<IP, IPC, D>
where
    D: Digest,
    IP: InnerProduct,
    IPC: IPCommitment<IP = IP>,
{
    pub fn new(r_transcript: Vec<Scalar<IPC>>, ck_base: FinalIPCommKey<IPC>) -> Self {
        Self {
            r_transcript,
            ck_base,
            _gipa: PhantomData,
        }
    }
}
