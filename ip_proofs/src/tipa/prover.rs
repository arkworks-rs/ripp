use ark_ec::{pairing::Pairing, CurveGroup};
use ark_ff::{Field, One, Zero};
use ark_poly::polynomial::{univariate::DensePolynomial, DenseUVPolynomial};
use ark_serialize::CanonicalSerialize;
use ark_std::{end_timer, start_timer};
use digest::Digest;

use crate::{gipa::GIPA, ip_commitment::Scalar, Error};

use super::*;

impl<P, D> TIPA<P, D>
where
    D: Digest,
    P: Pairing,
{
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
        })
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
