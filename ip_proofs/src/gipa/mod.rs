use ark_ff::Field;
use ark_serialize::CanonicalSerialize;
use ark_std::{convert::TryInto, marker::PhantomData, rand::Rng};
use digest::Digest;

use crate::{
    ip_commitment::{Commitment, IPCommitment, Scalar},
    Error,
};
use ark_inner_products::InnerProduct;

pub mod data_structures;
pub mod prover;
pub mod verifier;

pub use data_structures::*;

pub struct GIPA<IP, IPC, D>
where
    D: Digest,
    IP: InnerProduct,
    IPC: IPCommitment<IP = IP>,
{
    _inner_product: PhantomData<IP>,
    _inner_product_commitment: PhantomData<IPC>,
    _digest: PhantomData<D>,
}

impl<IP, IPC, D> GIPA<IP, IPC, D>
where
    D: Digest,
    IP: InnerProduct,
    IPC: IPCommitment<IP = IP>,
{
    pub fn setup<'a>(
        size: usize,
        rng: impl Rng,
    ) -> Result<(ProverKey<'a, IPC>, VerifierKey<'a, IPC>), Error> {
        let ck = IPC::setup(size, rng)?;
        let pk = ProverKey { ck: ck.clone() };
        let vk = VerifierKey { ck };
        Ok((pk, vk))
    }

    fn compute_challenge(
        transcript: &Scalar<IPC>,
        com_1: &Commitment<IPC>,
        com_2: &Commitment<IPC>,
    ) -> Result<(Scalar<IPC>, Scalar<IPC>), Error> {
        let mut counter_nonce = 0u32;
        let mut bytes = Vec::new();
        transcript.serialize_uncompressed(&mut bytes)?;
        com_1.serialize_uncompressed(&mut bytes)?;
        com_2.serialize_uncompressed(&mut bytes)?;
        let (c, c_inv) = 'challenge: loop {
            let mut hash_input = Vec::new();
            hash_input.extend_from_slice(&counter_nonce.to_be_bytes()[..]);
            hash_input.extend_from_slice(&bytes);

            let c = D::digest(&hash_input).as_slice()[..16].try_into().unwrap();
            let c: Scalar<IPC> = u128::from_be_bytes(c).into();
            if let Some(c_inv) = c.inverse() {
                // Optimization for multiexponentiation to rescale G2 elements with 128-bit challenge
                // Swap 'c' and 'c_inv' since can't control bit size of c_inv
                break 'challenge (c_inv, c);
            }
            counter_nonce += 1;
        };

        // println!("Debugging code; DO NOT USE THIS PRODUCTION!!!");
        // let c = transcript.double();
        // let c_inv = c.inverse().unwrap();
        // println!("(c, c_inv) = ({}, {})", c, c_inv);
        Ok((c, c_inv))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ip_commitment::{
        mipp::MSMCommitment, pairing::PairingCommitment, scalar::ScalarCommitment,
        snarkpack::TIPPCommitment,
    };

    use ark_bls12_381::{Bls12_381, Fr};
    use ark_ec::pairing::Pairing;
    use ark_ff::UniformRand;
    use blake2::Blake2b512;

    use ark_dh_commitments::random_generators;
    use ark_inner_products::{
        InnerProduct, MSMInnerProduct, PairingInnerProduct, ScalarInnerProduct,
    };

    const TEST_SIZE: usize = 4;

    #[test]
    fn pairing_inner_product_test() {
        type IP = PairingInnerProduct<Bls12_381>;
        type IPC = PairingCommitment<Bls12_381>;
        type PairingGIPA = GIPA<IP, IPC, Blake2b512>;

        let mut rng = ark_std::test_rng();
        let (pk, vk) = PairingGIPA::setup(TEST_SIZE, &mut rng).unwrap();
        let left = random_generators(&mut rng, TEST_SIZE);
        let right = random_generators(&mut rng, TEST_SIZE);

        let commitment = IPC::commit_with_ip(&pk.ck, &left, &right, None).unwrap();
        let random_challenge = Fr::rand(&mut rng);
        let output = IP::twisted_inner_product(&left, &right, random_challenge).unwrap();

        let instance = Instance {
            size: TEST_SIZE,
            output,
            commitment,
            random_challenge,
        };
        let witness = Witness { left, right };

        let proof = PairingGIPA::prove(&pk, &instance, &witness).unwrap();

        assert!(PairingGIPA::verify(&vk, &instance, &proof).unwrap());
    }

    #[test]
    fn snarkpack_inner_product_test() {
        type IP = PairingInnerProduct<Bls12_381>;
        type IPC = TIPPCommitment<Bls12_381>;
        type PairingGIPA = GIPA<IP, IPC, Blake2b512>;

        let mut rng = ark_std::test_rng();
        let (pk, vk) = PairingGIPA::setup(TEST_SIZE, &mut rng).unwrap();
        let left = random_generators(&mut rng, TEST_SIZE);
        let right = random_generators(&mut rng, TEST_SIZE);

        let commitment = IPC::commit_with_ip(&pk.ck, &left, &right, None).unwrap();
        let random_challenge = Fr::rand(&mut rng);
        let output = IP::twisted_inner_product(&left, &right, random_challenge).unwrap();

        let instance = Instance {
            size: TEST_SIZE,
            output,
            commitment,
            random_challenge,
        };
        let witness = Witness { left, right };

        let proof = PairingGIPA::prove(&pk, &instance, &witness).unwrap();

        assert!(PairingGIPA::verify(&vk, &instance, &proof).unwrap());
    }

    #[test]
    fn multiexponentiation_inner_product_test() {
        type IP = MSMInnerProduct<<Bls12_381 as Pairing>::G1>;
        type IPC = MSMCommitment<Bls12_381>;
        type MultiExpGIPA = GIPA<IP, IPC, Blake2b512>;

        let mut rng = ark_std::test_rng();
        let (pk, vk) = MultiExpGIPA::setup(TEST_SIZE, &mut rng).unwrap();
        let left = random_generators(&mut rng, TEST_SIZE);
        let mut right = Vec::new();
        for _ in 0..TEST_SIZE {
            right.push(<Bls12_381 as Pairing>::ScalarField::rand(&mut rng));
        }

        let commitment = IPC::commit_with_ip(&pk.ck, &left, &right, None).unwrap();

        let random_challenge = Fr::rand(&mut rng);
        let output = IP::twisted_inner_product(&left, &right, random_challenge).unwrap();

        let instance = Instance {
            size: TEST_SIZE,
            output,
            commitment,
            random_challenge,
        };
        let witness = Witness { left, right };

        let proof = MultiExpGIPA::prove(&pk, &instance, &witness).unwrap();

        assert!(MultiExpGIPA::verify(&vk, &instance, &proof).unwrap());
    }

    #[test]
    fn scalar_inner_product_test() {
        type IP = ScalarInnerProduct<<Bls12_381 as Pairing>::ScalarField>;
        type IPC = ScalarCommitment<Bls12_381>;
        type ScalarGIPA = GIPA<IP, IPC, Blake2b512>;

        let mut rng = ark_std::test_rng();
        let (pk, vk) = ScalarGIPA::setup(TEST_SIZE, &mut rng).unwrap();
        let mut left = Vec::new();
        let mut right = Vec::new();
        for _ in 0..TEST_SIZE {
            left.push(<Bls12_381 as Pairing>::ScalarField::rand(&mut rng));
            right.push(<Bls12_381 as Pairing>::ScalarField::rand(&mut rng));
        }

        // We pass a `None` challenge to `commit_with_ip` to avoid computing a useless
        // inner product
        let commitment = IPC::commit_with_ip(&pk.ck, &left, &right, None).unwrap();
        // Generate random challenge *after* committing to messages
        let random_challenge = Fr::rand(&mut rng);
        // Compute twisted inner product with respect to new challenge.
        let output = IP::twisted_inner_product(&left, &right, random_challenge).unwrap();

        let instance = Instance {
            size: TEST_SIZE,
            output,
            commitment,
            random_challenge,
        };

        let witness = Witness { left, right: right };

        let proof = ScalarGIPA::prove(&pk, &instance, &witness).unwrap();

        assert!(ScalarGIPA::verify(&vk, &instance, &proof).unwrap());
    }
}
