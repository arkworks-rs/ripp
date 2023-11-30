use ark_dh_commitments::afgho16::{AFGHOCommitmentG1, AFGHOCommitmentG2};
use ark_inner_products::PairingInnerProduct;

use super::generic::GenericCommitment;

pub type PairingCommitment<E> =
    GenericCommitment<PairingInnerProduct<E>, AFGHOCommitmentG1<E>, AFGHOCommitmentG2<E>>;
