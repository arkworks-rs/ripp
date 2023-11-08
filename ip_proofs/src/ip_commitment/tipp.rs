use std::marker::PhantomData;

use ark_dh_commitments::{
    afgho16::{AFGHOCommitmentG1, AFGHOCommitmentG2},
    identity::IdentityCommitment,
    pedersen::PedersenCommitment,
};
use ark_ec::pairing::{Pairing, PairingOutput};
use ark_inner_products::PairingInnerProduct;

use super::IPCommitment;

struct TIPPCommitment<E: Pairing>(PhantomData<E>);
type GC1<E> = AFGHOCommitmentG1<E>;
type GC2<E> = AFGHOCommitmentG2<E>;
type SC1<E> = PedersenCommitment<<E as Pairing>::G1>;
type SC2<E> = PedersenCommitment<<E as Pairing>::G2>;
type IP<E> = PairingInnerProduct<E>;
type IPC<E> = IdentityCommitment<PairingOutput<E>, <E as Pairing>::ScalarField>;

impl<E: Pairing> IPCommitment for TIPPCommitment<E> {
    type IP = IP<E>;
    type Scalar = E::ScalarField;
    type LeftKey = GC1<E>;
    type RightKey = GC2<E>;
    type IPKey = IPC<E>;
    type Commitment = PairingOutput<E>;
}
