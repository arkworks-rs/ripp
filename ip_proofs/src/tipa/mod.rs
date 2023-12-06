use ark_ec::pairing::Pairing;
use ark_ff::Field;
use ark_inner_products::PairingInnerProduct;
use ark_serialize::CanonicalSerialize;
use ark_std::marker::PhantomData;
use digest::Digest;

// pub mod structured_scalar_message;
pub mod data_structures;

pub mod prover;
pub mod setup;
pub mod verifier;

use crate::{
    ip_commitment::{snarkpack::TIPPCommitment, FinalIPCommKey},
    Error,
};
pub use data_structures::{Proof, ProverKey, VerifierKey};

type IP<P> = PairingInnerProduct<P>;
type IPC<P> = TIPPCommitment<P>;

//TODO: Could generalize: Don't need TIPA over G1 and G2, would work with G1 and G1 or over different pairing engines
pub trait TIPACompatibleSetup {}

//TODO: May need to add "reverse" MSMInnerProduct to allow for MIP with G2 messages (because TIP hard-coded G1 left and G2 right)
pub struct TIPA<P, D> {
    _pair: PhantomData<P>,
    _digest: PhantomData<D>,
}

impl<P: Pairing, D: Digest> TIPA<P, D> {
    pub fn compute_kzg_challenge<'a>(
        final_ck: &FinalIPCommKey<TIPPCommitment<P>>,
        challenges: &[P::ScalarField],
    ) -> Result<P::ScalarField, Error> {
        // KZG challenge point
        let mut bytes = Vec::new();
        challenges[0].serialize_uncompressed(&mut bytes)?;
        final_ck.ck_a.serialize_uncompressed(&mut bytes)?;
        final_ck.ck_b.serialize_uncompressed(&mut bytes)?;

        let mut counter_nonce: usize = 0;
        let c = loop {
            let mut hash_input = bytes.clone();
            hash_input.extend_from_slice(&counter_nonce.to_be_bytes()[..]);
            if let Some(c) = P::ScalarField::from_random_bytes(&D::digest(&hash_input)) {
                break c;
            };
            counter_nonce += 1;
        };
        Ok(c)
    }
}

#[cfg(test)]
mod tests {
    use crate::{
        gipa::{Instance, Witness},
        ip_commitment::IPCommitment,
    };

    use super::*;
    use ark_bls12_381::{Bls12_381, Fr};
    use ark_std::UniformRand;
    use blake2::Blake2b512;

    use ark_dh_commitments::random_generators;
    use ark_inner_products::{InnerProduct, PairingInnerProduct};

    const TEST_SIZE: usize = 8;

    #[test]
    fn pairing_inner_product_test() {
        type IP = PairingInnerProduct<Bls12_381>;
        type IPC = TIPPCommitment<Bls12_381>;
        type PairingTIPA = TIPA<Bls12_381, Blake2b512>;

        let mut rng = ark_std::test_rng();
        let (pk, vk) = PairingTIPA::setup(TEST_SIZE, &mut rng).unwrap();
        let left = random_generators(&mut rng, TEST_SIZE);
        let right = random_generators(&mut rng, TEST_SIZE);

        let commitment = IPC::commit_with_ip(&pk.pk.ck, &left, &right, None).unwrap();
        let twist = Fr::rand(&mut rng);
        let output = IP::twisted_inner_product(&left, &right, twist).unwrap();

        let instance = Instance {
            size: TEST_SIZE,
            output,
            commitment,
            twist,
        };
        let witness = Witness { left, right };

        let proof = PairingTIPA::prove(&pk, &instance, &witness).unwrap();

        assert!(PairingTIPA::verify(&vk, &instance, &proof).unwrap());
    }
}
