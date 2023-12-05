use ark_dh_commitments::pedersen::PedersenCommitment;
use ark_ec::pairing::Pairing;
use ark_inner_products::ScalarInnerProduct;

use super::generic::GenericCommitment;

pub type ScalarCommitment<E> = GenericCommitment<
    ScalarInnerProduct<<E as Pairing>::ScalarField>,
    PedersenCommitment<<E as Pairing>::G1>,
    PedersenCommitment<<E as Pairing>::G2>,
>;