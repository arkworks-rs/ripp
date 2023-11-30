use ark_ec::{pairing::Pairing, scalar_mul::fixed_base::FixedBase, AffineRepr, CurveGroup, Group};
use ark_ff::{Field, One, PrimeField, UniformRand, Zero};
use ark_poly::polynomial::{univariate::DensePolynomial, DenseUVPolynomial};
use ark_serialize::CanonicalSerialize;
use ark_std::rand::Rng;
use ark_std::{end_timer, start_timer};
use digest::Digest;
use itertools::Itertools;
use std::marker::PhantomData;

use crate::{
    gipa::GIPA,
    ip_commitment::{IPCommKey, IPCommitment, Scalar},
    Error,
};
use ark_inner_products::{InnerProduct, PairingInnerProduct};

// pub mod structured_scalar_message;
pub mod data_structures;
pub mod kzg;
pub mod tipp;

pub mod setup;
pub mod prover;
pub mod verifier;

use data_structures::{Proof, ProverKey, VerifierKey};
use tipp::TIPPCommitment;

use self::data_structures::GenericSRS;

type IP<P> = PairingInnerProduct<P>;
type IPC<P> = TIPPCommitment<P>;
type LeftMessage<P> = <IP<P> as InnerProduct>::LeftMessage;
type RightMessage<P> = <IP<P> as InnerProduct>::RightMessage;
type Commitment<P> = <IPC<P> as IPCommitment>::Commitment;

//TODO: Could generalize: Don't need TIPA over G1 and G2, would work with G1 and G1 or over different pairing engines
pub trait TIPACompatibleSetup {}

//TODO: May need to add "reverse" MSMInnerProduct to allow for MIP with G2 messages (because TIP hard-coded G1 left and G2 right)
pub struct TIPA<P, D> {
    _pair: PhantomData<P>,
    _digest: PhantomData<D>,
}

#[cfg(test)]
mod tests {
    use crate::ip_commitment::pairing::PairingCommitment;

    use super::*;
    use ark_bls12_381::Bls12_381;
    use ark_ec::pairing::PairingOutput;
    use ark_std::rand::{rngs::StdRng, SeedableRng};
    use blake2::Blake2b;

    use ark_dh_commitments::{
        afgho16::{AFGHOCommitmentG1, AFGHOCommitmentG2},
        identity::IdentityCommitment,
        pedersen::PedersenCommitment,
        random_generators,
    };
    use ark_inner_products::{
        InnerProduct, MSMInnerProduct, PairingInnerProduct, ScalarInnerProduct,
    };

    pub fn structured_scalar_power<F: Field>(num: usize, s: &F) -> Vec<F> {
        let mut powers = vec![F::one()];
        for i in 1..num {
            powers.push(powers[i - 1] * s);
        }
        powers
    }

    type GC1 = AFGHOCommitmentG1<Bls12_381>;
    type GC2 = AFGHOCommitmentG2<Bls12_381>;
    type SC1 = PedersenCommitment<<Bls12_381 as Pairing>::G1>;
    type SC2 = PedersenCommitment<<Bls12_381 as Pairing>::G2>;

    const TEST_SIZE: usize = 8;

    #[test]
    fn pairing_inner_product_test() {
        type IP = PairingInnerProduct<Bls12_381>;
        type IPC = PairingCommitment<Bls12_381>;
        type PairingTIPA = TIPA<Bls12_381, Blake2b>;

        let mut rng = ark_std::test_rng();
        let (srs, ck_t) = PairingTIPA::setup(TEST_SIZE, &mut rng).unwrap();
        let ck = srs.ck();
        let (ck_a, ck_b) = srs.get_commitment_keys();
        let v_srs = srs.get_verifier_key();
        let m_a = random_generators(&mut rng, TEST_SIZE);
        let m_b = random_generators(&mut rng, TEST_SIZE);
        let com = IPC::commit();
        let t = vec![IP::inner_product(&m_a, &m_b).unwrap()];
        let com_t = IPC::commit(&vec![ck_t.clone()], &t).unwrap();

        let proof = PairingTIPA::prove(&srs, (&m_a, &m_b), (&ck_a, &ck_b, &ck_t)).unwrap();

        assert!(PairingTIPA::verify(&v_srs, &ck_t, (&com_a, &com_b, &com_t), &proof).unwrap());
    }

    #[test]
    fn multiexponentiation_inner_product_test() {
        type IP = MSMInnerProduct<<Bls12_381 as Pairing>::G1>;
        type IPC =
            IdentityCommitment<<Bls12_381 as Pairing>::G1, <Bls12_381 as Pairing>::ScalarField>;
        type MultiExpTIPA = TIPA<Bls12_381, Blake2b>;

        let mut rng = ark_std::test_rng();
        let (srs, ck_t) = MultiExpTIPA::setup(TEST_SIZE, &mut rng).unwrap();
        let (ck_a, ck_b) = srs.get_commitment_keys();
        let v_srs = srs.get_verifier_key();
        let m_a = random_generators(&mut rng, TEST_SIZE);
        let mut m_b = Vec::new();
        for _ in 0..TEST_SIZE {
            m_b.push(<Bls12_381 as Pairing>::ScalarField::rand(&mut rng));
        }
        let com_a = GC1::commit(&ck_a, &m_a).unwrap();
        let com_b = SC1::commit(&ck_b, &m_b).unwrap();
        let t = vec![IP::inner_product(&m_a, &m_b).unwrap()];
        let com_t = IPC::commit(&vec![ck_t.clone()], &t).unwrap();

        let proof = MultiExpTIPA::prove(&srs, (&m_a, &m_b), (&ck_a, &ck_b, &ck_t)).unwrap();

        assert!(MultiExpTIPA::verify(&v_srs, &ck_t, (&com_a, &com_b, &com_t), &proof).unwrap());
    }

    #[test]
    fn scalar_inner_product_test() {
        type IP = ScalarInnerProduct<<Bls12_381 as Pairing>::ScalarField>;
        type IPC = IdentityCommitment<
            <Bls12_381 as Pairing>::ScalarField,
            <Bls12_381 as Pairing>::ScalarField,
        >;
        type ScalarTIPA = TIPA<Bls12_381, Blake2b>;

        let mut rng = StdRng::seed_from_u64(0u64);
        let (srs, ck_t) = ScalarTIPA::setup(TEST_SIZE, &mut rng).unwrap();
        let (ck_a, ck_b) = srs.get_commitment_keys();
        let v_srs = srs.get_verifier_key();
        let mut m_a = Vec::new();
        let mut m_b = Vec::new();
        for _ in 0..TEST_SIZE {
            m_a.push(<Bls12_381 as Pairing>::ScalarField::rand(&mut rng));
            m_b.push(<Bls12_381 as Pairing>::ScalarField::rand(&mut rng));
        }
        let com_a = SC2::commit(&ck_a, &m_a).unwrap();
        let com_b = SC1::commit(&ck_b, &m_b).unwrap();
        let t = vec![IP::inner_product(&m_a, &m_b).unwrap()];
        let com_t = IPC::commit(&vec![ck_t.clone()], &t).unwrap();

        let proof = ScalarTIPA::prove(&srs, (&m_a, &m_b), (&ck_a, &ck_b, &ck_t)).unwrap();

        assert!(ScalarTIPA::verify(&v_srs, &ck_t, (&com_a, &com_b, &com_t), &proof).unwrap());
    }

    #[test]
    fn pairing_inner_product_with_srs_shift_test() {
        type IP = PairingInnerProduct<Bls12_381>;
        type IPC =
            IdentityCommitment<PairingOutput<Bls12_381>, <Bls12_381 as Pairing>::ScalarField>;
        type PairingTIPA = TIPA<Bls12_381, Blake2b>;

        let mut rng = ark_std::test_rng();
        let (srs, ck_t) = PairingTIPA::setup(TEST_SIZE, &mut rng).unwrap();
        let (ck_a, ck_b) = srs.get_commitment_keys();
        let v_srs = srs.get_verifier_key();

        let m_a = random_generators(&mut rng, TEST_SIZE);
        let m_b = random_generators(&mut rng, TEST_SIZE);
        let com_a = GC1::commit(&ck_a, &m_a).unwrap();
        let com_b = GC2::commit(&ck_b, &m_b).unwrap();

        let r_scalar = <<Bls12_381 as Pairing>::ScalarField>::rand(&mut rng);
        let r_vec = structured_scalar_power(TEST_SIZE, &r_scalar);
        let m_a_r = m_a
            .iter()
            .zip(&r_vec)
            .map(|(&a, r)| a * r)
            .collect::<Vec<<Bls12_381 as Pairing>::G1>>();
        let ck_a_r = ck_a
            .iter()
            .zip(&r_vec)
            .map(|(&ck, r)| ck * r.inverse().unwrap())
            .collect::<Vec<<Bls12_381 as Pairing>::G2>>();

        let t = vec![IP::inner_product(&m_a_r, &m_b).unwrap()];
        let com_t = IPC::commit(&vec![ck_t.clone()], &t).unwrap();

        assert_eq!(com_a, IP::inner_product(&m_a_r, &ck_a_r).unwrap());

        let proof = PairingTIPA::prove_with_srs_shift(
            &srs,
            (&m_a_r, &m_b),
            (&ck_a_r, &ck_b, &ck_t),
            &r_scalar,
        )
        .unwrap();

        assert!(PairingTIPA::verify_with_srs_shift(
            &v_srs,
            &ck_t,
            (&com_a, &com_b, &com_t),
            &proof,
            &r_scalar
        )
        .unwrap());
    }
}
