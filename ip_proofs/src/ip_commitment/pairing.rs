use ark_dh_commitments::{
    afgho16::{AFGHOCommitmentG1, AFGHOCommitmentG2},
    identity::IdentityCommitment,
};
use ark_ec::pairing::{Pairing, PairingOutput};
use ark_inner_products::PairingInnerProduct;

use super::generic::GenericCommitment;

pub type PairingCommitment<E> = GenericCommitment<
    PairingInnerProduct<E>,
    AFGHOCommitmentG1<E>,
    AFGHOCommitmentG2<E>,
    IdentityCommitment<PairingOutput<E>, <E as Pairing>::ScalarField>,
>;
