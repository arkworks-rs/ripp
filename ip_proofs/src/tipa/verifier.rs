use ark_ec::{pairing::Pairing, scalar_mul::fixed_base::FixedBase, AffineRepr, CurveGroup};
use ark_ff::{Field, One, Zero};
use ark_serialize::CanonicalSerialize;
use digest::Digest;

use crate::{
    gipa::GIPA,
    ip_commitment::{IPCommKey, Scalar},
    Error,
};

use super::*;

impl<P, D> TIPA<P, D>
where
    D: Digest,
    P: Pairing,
{
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
