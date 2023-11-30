use ark_dh_commitments::{afgho16::AFGHOCommitmentG1, pedersen::PedersenCommitment};
use ark_ec::pairing::Pairing;
use ark_inner_products::MSMInnerProduct;

use super::generic::GenericCommitment;

pub type MSMCommitment<E> = GenericCommitment<
    MSMInnerProduct<<E as Pairing>::G1>,
    AFGHOCommitmentG1<E>,
    PedersenCommitment<<E as Pairing>::G1>,
>;
