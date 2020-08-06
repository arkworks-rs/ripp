use rand::Rng;
use std::{
    error::Error as ErrorTrait,
    fmt::{Display, Formatter, Result as FmtResult},
    ops::MulAssign,
};

use dh_commitments::DoublyHomomorphicCommitment;
use inner_products::InnerProduct;

pub mod gipa;
pub mod tipa;

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
    type SetupOutput;

    fn setup<R: Rng>(rng: &mut R, size: usize) -> Result<Self::SetupOutput, Error>;

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

//TODO: helper function for mul because relying on MulAssign
pub(crate) fn mul_helper<T: MulAssign<F> + Clone, F: Clone>(t: &T, f: &F) -> T {
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
            }
            InnerProductArgumentError::InnerProductInvalid => "inner product not sound".to_string(),
        };
        write!(f, "{}", msg)
    }
}
