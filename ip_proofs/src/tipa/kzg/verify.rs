use ark_ec::{
    pairing::{Pairing, PairingOutput},
    AffineRepr, CurveGroup,
};
use ark_ff::Field;
use crossbeam_channel::Sender;

use crate::{pairing_check::PairingCheck, srs::VerifierKey};

use super::{evaluate_ipa_polynomial, EvaluationProof};

/// verify_kzg_opening_g2 takes a KZG opening, the final commitment key, SRS and
/// any shift (in TIPP we shift the v commitment by r^-1) and returns a pairing
/// tuple to check if the opening is correct or not.
pub fn verify_kzg_v<E: Pairing>(
    v_srs: &VerifierKey<E>,
    final_vkey: &(E::G2Affine, E::G2Affine),
    proof: &EvaluationProof<E::G2Affine>,
    challenges: &[E::ScalarField],
    point: E::ScalarField,
    checks: Sender<Option<PairingCheck<E>>>,
) {
    // f_v(z)
    let v_poly_at_z = evaluate_ipa_polynomial(challenges, point, E::ScalarField::ONE);
    // -g such that when we test a pairing equation we only need to check if
    // it's equal 1 at the end:
    // e(a,b) = e(c,d) <=> e(a,b)e(-c,d) = 1
    // e(A,B) = e(C,D) <=> e(A,B)e(-C,D) == 1 <=> e(A,B)e(C,D)^-1 == 1

    let v1 = checks.clone();
    let v2 = checks.clone();

    par! {
        // e(g, C_f * h^{-y}) == e(v1 * g^{-x}, \pi) = 1
        let _check1 = kzg_check_v::<E>(
            v_srs,
            point,
            v_poly_at_z,
            final_vkey.0,
            v_srs.g_alpha,
            proof.0,
            v1,
        ),

        // e(g, C_f * h^{-y}) == e(v2 * g^{-x}, \pi) = 1
        let _check2 = kzg_check_v::<E>(
            v_srs,
            point,
            v_poly_at_z,
            final_vkey.1,
            v_srs.g_beta,
            proof.1,
            v2,
        )
    };
}

fn kzg_check_v<E: Pairing>(
    v_srs: &VerifierKey<E>,
    x: E::ScalarField,
    y: E::ScalarField,
    cf: E::G2Affine,
    vk: E::G1Affine,
    pi: E::G2Affine,
    checks: Sender<Option<PairingCheck<E>>>,
) {
    // KZG Check: e(g, C_f * h^{-y}) = e(vk * g^{-x}, \pi)
    // Transformed, such that
    // e(-g, C_f * h^{-y}) * e(vk * g^{-x}, \pi) = 1

    // C_f - (y * h)
    let b = (cf.into_group() - v_srs.h * y).into();

    // vk - (g * x)
    let c = (vk.into_group() - v_srs.g * x).into_affine();
    let p = PairingCheck::rand(
        &[(v_srs.neg_g, b), (c, pi)],
        PairingOutput(E::TargetField::ONE),
    );
    checks.send(Some(p)).unwrap();
}

/// Similar to verify_kzg_opening_g2 but for g1.
pub fn verify_kzg_w<E: Pairing>(
    v_srs: &VerifierKey<E>,
    final_wkey: &(E::G1Affine, E::G1Affine),
    wkey_opening: &EvaluationProof<E::G1Affine>,
    challenges: &[E::ScalarField],
    r_shift: E::ScalarField,
    kzg_challenge: E::ScalarField,
    checks: Sender<Option<PairingCheck<E>>>,
) {
    // compute in parallel f(z) and z^n and then combine into f_w(z) = z^n * f(z)
    par! {
        let fz = evaluate_ipa_polynomial(challenges, kzg_challenge, r_shift),
        let zn = kzg_challenge.pow(&[v_srs.n as u64])
    };

    let w_poly_at_z = fz * zn;

    let w1 = checks.clone();
    let w2 = checks.clone();
    par! {
        // e(C_f * g^{-y}, h) = e(\pi, w1 * h^{-x})
        let _check1 = kzg_check_w::<E>(
            v_srs,
            kzg_challenge,
            w_poly_at_z,
            final_wkey.0,
            v_srs.h_alpha,
            wkey_opening.0,
            w1,
        ),

        // e(C_f * g^{-y}, h) = e(\pi, w2 * h^{-x})
        let _check2 = kzg_check_w::<E>(
            v_srs,
            kzg_challenge,
            w_poly_at_z,
            final_wkey.1,
            v_srs.h_beta,
            wkey_opening.1,
            w2,
        )
    };
}

fn kzg_check_w<E: Pairing>(
    v_srs: &VerifierKey<E>,
    x: E::ScalarField,
    y: E::ScalarField,
    cf: E::G1Affine,
    wk: E::G2Affine,
    pi: E::G1Affine,
    checks: Sender<Option<PairingCheck<E>>>,
) {
    // KZG Check: e(C_f * g^{-y}, h) = e(\pi, wk * h^{-x})
    // Transformed, such that
    // e(C_f * g^{-y}, -h) * e(\pi, wk * h^{-x}) = 1

    // C_f - (y * g)
    let a = (cf.into_group() - v_srs.g * y).into();

    // wk - (x * h)
    let d = (wk.into_group() - v_srs.h * x).into();
    let p = PairingCheck::rand(
        &[(a, v_srs.neg_h), (pi, d)],
        PairingOutput(E::TargetField::ONE),
    );
    checks.send(Some(p)).unwrap();
}
