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

impl<P, D> TIPA<P, D>
where
    D: Digest,
    P: Pairing,
{
    pub fn setup<'a>(
        size: usize,
        rng: &mut impl Rng,
    ) -> Result<(ProverKey<P>, VerifierKey<P>), Error> {
        let srs = GenericSRS::sample(size, rng);
        Ok(srs.specialize(size))
    }

    pub fn prove<'a>(
        pk: &ProverKey<P>,
        left: &[LeftMessage<P>],
        right: &[RightMessage<P>],
    ) -> Result<Proof<P, D>, Error> {
        Self::prove_with_srs_shift(pk, left, right, &P::ScalarField::one())
    }

    // Shifts KZG proof for left message by scalar r (used for efficient composition with aggregation protocols)
    // LMC commitment key should already be shifted before being passed as input
    pub fn prove_with_srs_shift<'a>(
        pk: &ProverKey<P>,
        left: &[LeftMessage<P>],
        right: &[RightMessage<P>],
        r_shift: &P::ScalarField,
    ) -> Result<Proof<P, D>, Error> {
        // Run GIPA
        let (proof, aux) = <GIPA<IP<P>, IPC<P>, D>>::prove_with_aux(&pk.ck, left, right)?;

        // Unpack the final commitment key. Make sure it's len 1
        let ck_base = aux.ck_base;

        // Prove final commitment keys are wellformed
        let transcript = aux.r_transcript;
        let transcript_inverse = transcript
            .iter()
            .map(|x| x.inverse().unwrap())
            .collect::<Vec<_>>();
        let r_inverse = r_shift.inverse().unwrap();

        // KZG challenge point
        let mut counter_nonce: usize = 0;
        let c = loop {
            let mut hash_input = Vec::new();
            hash_input.extend_from_slice(&counter_nonce.to_be_bytes()[..]);
            transcript
                .first()
                .unwrap()
                .serialize_uncompressed(&mut hash_input)?;
            ck_base.ck_a.serialize_uncompressed(&mut hash_input)?;
            ck_base.ck_b.serialize_uncompressed(&mut hash_input)?;
            if let Some(c) = Scalar::<IPC<P>>::from_random_bytes(&D::digest(&hash_input)) {
                break c;
            };
            counter_nonce += 1;
        };

        // Complete KZG proofs
        let ck_a_kzg_opening =
            prove_commitment_key_kzg_opening(&pk.h_beta_powers, &transcript_inverse, r_inverse, c)?;
        let ck_b_kzg_opening = prove_commitment_key_kzg_opening(
            &pk.g_alpha_powers,
            &transcript,
            <P::ScalarField>::one(),
            c,
        )?;

        Ok(Proof {
            gipa_proof: proof,
            final_ck: ck_base.into(),
            final_ck_proof: (ck_a_kzg_opening, ck_b_kzg_opening),
            _pair: PhantomData,
        })
    }

    pub fn verify<'a>(
        v_srs: &VerifierKey<P>,
        ck: &IPCommKey<'a, IPC<P>>,
        com: &Commitment<P>,
        proof: &Proof<P, D>,
    ) -> Result<bool, Error> {
        Self::verify_with_srs_shift(v_srs, ck, com, proof, P::ScalarField::one())
    }

    pub fn verify_with_srs_shift<'a>(
        v_srs: &VerifierKey<P>,
        ck: &IPCommKey<'a, IPC<P>>,
        com: &Commitment<P>,
        proof: &Proof<P, D>,
        r_shift: P::ScalarField,
    ) -> Result<bool, Error> {
        let (base_com, transcript) =
            GIPA::verify_recursive_challenge_transcript(com, &proof.gipa_proof)?;
        let transcript_inverse = transcript
            .iter()
            .map(|x| x.inverse().unwrap())
            .collect::<Vec<_>>();

        // Verify commitment keys wellformed
        let ck_final = &proof.final_ck;
        let (ck_a_proof, ck_b_proof) = &proof.final_ck_proof;

        // KZG challenge point
        let mut counter_nonce: usize = 0;
        let c = loop {
            let mut hash_input = Vec::new();
            hash_input.extend_from_slice(&counter_nonce.to_be_bytes()[..]);
            transcript
                .first()
                .unwrap()
                .serialize_uncompressed(&mut hash_input)?;
            ck_final.serialize_uncompressed(&mut hash_input)?;
            if let Some(c) = Scalar::<IPC<P>>::from_random_bytes(&D::digest(&hash_input)) {
                break c;
            };
            counter_nonce += 1;
        };

        let ck_a_valid = verify_commitment_key_g2_kzg_opening(
            v_srs,
            &ck_final.ck_a,
            &ck_a_proof,
            &transcript_inverse,
            r_shift.inverse().unwrap(),
            c,
        )?;
        let ck_b_valid = verify_commitment_key_g1_kzg_opening(
            v_srs,
            &ck_final.ck_b,
            &ck_b_proof,
            &transcript,
            <P::ScalarField>::one(),
            c,
        )?;

        // Verify base inner product commitment
        let (com_a, com_b, com_t) = base_com;
        let a_base = vec![proof.gipa_proof.r_base.0.clone()];
        let b_base = vec![proof.gipa_proof.r_base.1.clone()];
        let t_base = IP::inner_product(&a_base, &b_base)?;
        let base_valid = IPC::verify(&ck_final.into(), &a_base, &b_base, &t_base, com)?;

        Ok(ck_a_valid && ck_b_valid && base_valid)
    }
}

pub fn prove_commitment_key_kzg_opening<G: CurveGroup>(
    srs_powers: &[G::Affine],
    transcript: &[G::ScalarField],
    r_shift: G::ScalarField,
    kzg_challenge: G::ScalarField,
) -> Result<G, Error> {
    let ck_polynomial = DensePolynomial::from_coefficients_slice(
        &polynomial_coefficients_from_transcript(transcript, r_shift),
    );
    assert_eq!(srs_powers.len(), ck_polynomial.coeffs.len());

    let eval = start_timer!(|| "polynomial eval");
    let ck_polynomial_c_eval =
        polynomial_evaluation_product_form_from_transcript(&transcript, kzg_challenge, r_shift);
    end_timer!(eval);

    let quotient = start_timer!(|| "polynomial quotient");
    let quotient_polynomial = &ck_polynomial
        / &(DensePolynomial::from_coefficients_vec(vec![-kzg_challenge, G::ScalarField::one()]));
    end_timer!(quotient);

    let mut quotient_polynomial_coeffs = quotient_polynomial.coeffs;
    quotient_polynomial_coeffs.resize(srs_powers.len(), <G::ScalarField>::zero());

    let multiexp = start_timer!(|| "opening multiexp");
    let opening = G::msm(srs_powers, &quotient_polynomial_coeffs).unwrap();
    end_timer!(multiexp);
    Ok(opening)
}

//TODO: Figure out how to avoid needing two separate methods for verification of opposite groups
pub fn verify_commitment_key_g2_kzg_opening<P: Pairing>(
    v_srs: &VerifierKey<P>,
    ck_final: &P::G2,
    ck_opening: &P::G2,
    transcript: &[P::ScalarField],
    r_shift: P::ScalarField,
    kzg_challenge: P::ScalarField,
) -> Result<bool, Error> {
    let ck_polynomial_c_eval =
        polynomial_evaluation_product_form_from_transcript(transcript, kzg_challenge, r_shift);
    Ok(
        P::pairing(v_srs.g, *ck_final - v_srs.h * ck_polynomial_c_eval)
            == P::pairing(
                v_srs.g_beta.into_group() - v_srs.g * kzg_challenge,
                *ck_opening,
            ),
    )
}

pub fn verify_commitment_key_g1_kzg_opening<P: Pairing>(
    v_srs: &VerifierKey<P>,
    ck_final: &P::G1,
    ck_opening: &P::G1,
    transcript: &[P::ScalarField],
    r_shift: P::ScalarField,
    kzg_challenge: P::ScalarField,
) -> Result<bool, Error> {
    let ck_polynomial_c_eval =
        polynomial_evaluation_product_form_from_transcript(transcript, kzg_challenge, r_shift);
    Ok(
        P::pairing(*ck_final - v_srs.g * ck_polynomial_c_eval, v_srs.h)
            == P::pairing(
                *ck_opening,
                v_srs.h_alpha.into_group() - v_srs.h * kzg_challenge,
            ),
    )
}

pub fn structured_generators_scalar_power<G: CurveGroup>(
    num: usize,
    g: &G,
    s: &G::ScalarField,
) -> Vec<G> {
    assert!(num > 0);
    let mut powers_of_scalar = vec![];
    let mut pow_s = G::ScalarField::one();
    for _ in 0..num {
        powers_of_scalar.push(pow_s);
        pow_s *= s;
    }

    let window_size = FixedBase::get_mul_window_size(num);

    let scalar_bits = G::ScalarField::MODULUS_BIT_SIZE as usize;
    let g_table = FixedBase::get_window_table(scalar_bits, window_size, g.clone());
    let powers_of_g = FixedBase::msm::<G>(scalar_bits, window_size, &g_table, &powers_of_scalar);
    powers_of_g
}

fn polynomial_evaluation_product_form_from_transcript<F: Field>(
    transcript: &[F],
    z: F,
    r_shift: F,
) -> F {
    let mut power_2_zr = (z.clone() * z) * r_shift;
    let mut product_form = Vec::new();
    for x in transcript.iter() {
        product_form.push(F::one() + (x.clone() * &power_2_zr));
        power_2_zr *= power_2_zr;
    }
    product_form.iter().product()
}

fn polynomial_coefficients_from_transcript<F: Field>(transcript: &[F], r_shift: F) -> Vec<F> {
    let mut coefficients = Vec::with_capacity(2usize.pow((transcript.len() + 1) as u32));
    coefficients.push(F::one());
    let mut coefficients = vec![F::one()];
    let mut power_2_r = r_shift.clone();
    for (i, x) in transcript.iter().enumerate() {
        for j in 0..(2_usize).pow(i as u32) {
            coefficients.push(coefficients[j] * &(x.clone() * &power_2_r));
        }
        power_2_r *= power_2_r;
    }
    // Interleave with 0 coefficients
    coefficients
        .into_iter()
        .interleave([F::zero()].into_iter().cycle().take(coefficients.len() - 1))
        .collect()
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
