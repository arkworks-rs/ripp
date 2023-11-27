use ark_ff::Field;
use ark_serialize::CanonicalSerialize;
use ark_std::{convert::TryInto, marker::PhantomData, rand::Rng};
use digest::Digest;

use crate::{
    ip_commitment::{IPCommKey, IPCommitment, Scalar},
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
    pub fn setup<'a>(size: usize, rng: impl Rng) -> Result<IPCommKey<'a, IPC>, Error> {
        IPC::setup(size, rng)
    }

    fn compute_challenge(
        transcript: &Scalar<IPC>,
        com_1: &IPC::Commitment,
        com_2: &IPC::Commitment,
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

        Ok((c, c_inv))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ip_commitment::{
        mipp::MSMCommitment, pairing::PairingCommitment, scalar::ScalarCommitment,
    };

    use ark_bls12_381::Bls12_381;
    use ark_ec::pairing::Pairing;
    use ark_ff::UniformRand;
    use blake2::Blake2b;

    use ark_dh_commitments::random_generators;
    use ark_inner_products::{
        InnerProduct, MSMInnerProduct, PairingInnerProduct, ScalarInnerProduct,
    };

    const TEST_SIZE: usize = 8;

    #[test]
    fn pairing_inner_product_test() {
        type IP = PairingInnerProduct<Bls12_381>;
        type IPC = PairingCommitment<Bls12_381>;
        type PairingGIPA = GIPA<IP, IPC, Blake2b>;

        let mut rng = ark_std::test_rng();
        let ck = PairingGIPA::setup(TEST_SIZE, &mut rng).unwrap();
        let m_a = random_generators(&mut rng, TEST_SIZE);
        let m_b = random_generators(&mut rng, TEST_SIZE);
        let t = IP::inner_product(&m_a, &m_b).unwrap();
        let com = IPC::commit(&ck, &m_a, &m_b, || t).unwrap();
        let instance = Instance {
            size: TEST_SIZE,
            output: t,
            commitment: com.clone(),
            random_challenge: <Bls12_381 as Pairing>::ScalarField::rand(&mut rng),
        };
        let witness = Witness {
            left: m_a,
            right: m_b,
        };

        let proof = PairingGIPA::prove(&ck, &instance, &witness).unwrap();

        assert!(PairingGIPA::verify(&ck, &instance, &proof).unwrap());
    }

    #[test]
    fn multiexponentiation_inner_product_test() {
        type IP = MSMInnerProduct<<Bls12_381 as Pairing>::G1>;
        type IPC = MSMCommitment<Bls12_381>;
        type MultiExpGIPA = GIPA<IP, IPC, Blake2b>;

        let mut rng = ark_std::test_rng();
        let ck = MultiExpGIPA::setup(TEST_SIZE, &mut rng).unwrap();
        let m_a = random_generators(&mut rng, TEST_SIZE);
        let mut m_b = Vec::new();
        for _ in 0..TEST_SIZE {
            m_b.push(<Bls12_381 as Pairing>::ScalarField::rand(&mut rng));
        }
        let t = IP::inner_product(&m_a, &m_b).unwrap();
        let com = IPC::commit(&ck, &m_a, &m_b, || t).unwrap();
        let instance = Instance {
            size: TEST_SIZE,
            output: t,
            commitment: com.clone(),
            random_challenge: <Bls12_381 as Pairing>::ScalarField::rand(&mut rng),
        };
        let witness = Witness {
            left: m_a,
            right: m_b,
        };

        let proof = MultiExpGIPA::prove(&ck, &instance, &witness).unwrap();

        assert!(MultiExpGIPA::verify(&ck, &instance, &proof).unwrap());
    }

    #[test]
    fn scalar_inner_product_test() {
        type IP = ScalarInnerProduct<<Bls12_381 as Pairing>::ScalarField>;
        type IPC = ScalarCommitment<Bls12_381>;
        type ScalarGIPA = GIPA<IP, IPC, Blake2b>;

        let mut rng = ark_std::test_rng();
        let ck = ScalarGIPA::setup(TEST_SIZE, &mut rng).unwrap();
        let mut m_a = Vec::new();
        let mut m_b = Vec::new();
        for _ in 0..TEST_SIZE {
            m_a.push(<Bls12_381 as Pairing>::ScalarField::rand(&mut rng));
            m_b.push(<Bls12_381 as Pairing>::ScalarField::rand(&mut rng));
        }

        let t = IP::inner_product(&m_a, &m_b).unwrap();
        let com = IPC::commit(&ck, &m_a, &m_b, || t).unwrap();
        let instance = Instance {
            size: TEST_SIZE,
            output: t,
            commitment: com.clone(),
            random_challenge: <Bls12_381 as Pairing>::ScalarField::rand(&mut rng),
        };
        let witness = Witness {
            left: m_a,
            right: m_b,
        };

        let proof = ScalarGIPA::prove(&ck, &instance, &witness).unwrap();

        assert!(ScalarGIPA::verify(&ck, &instance, &proof).unwrap());
    }
}
