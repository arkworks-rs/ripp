use ark_ec::{pairing::Pairing, AffineRepr, CurveGroup};
use ark_ff::{Field, Zero};

use crate::{
    ip_commitment::snarkpack::{kzg::evaluate_ipa_polynomial_shifted, LeftKey, RightKey},
    tipa::VerifierKey,
};

use super::{evaluate_ipa_polynomial, EvaluationProof};

/// verify_kzg_opening_g2 takes a KZG opening, the final commitment key, SRS and
/// any shift (in TIPP we shift the v commitment by r^-1) and returns a pairing
/// tuple to check if the opening is correct or not.
#[must_use]
pub fn verify_left_key<E: Pairing>(
    vk: &VerifierKey<E>,
    // These are the KZG commitments that arise from the final commitment key
    &LeftKey { v_1, v_2 }: &LeftKey<E>,
    &EvaluationProof(proof_1, proof_2): &EvaluationProof<E::G2Affine>,
    challenges_inv: &[E::ScalarField],
    point: E::ScalarField,
) -> bool {
    // f_v(z)
    let v_poly_at_point = evaluate_ipa_polynomial(challenges_inv, point, E::ScalarField::ONE);
    let [v_1, v_2]: [E::G2Affine; 2] = E::G2::normalize_batch(&[v_1, v_2]).try_into().unwrap();

    let check1 = check_left::<E>(vk, vk.g_alpha, point, v_poly_at_point, v_1, proof_1);

    let check2 = check_left::<E>(vk, vk.g_beta, point, v_poly_at_point, v_2, proof_2);
    check1 & check2
}

/// Here `tau` is the SRS secret.
#[must_use]
pub(super) fn check_left<E: Pairing>(
    vk: &VerifierKey<E>,
    tau_g: E::G1Affine,
    point: E::ScalarField,
    eval: E::ScalarField,
    commitment: E::G2Affine,
    evaluation_proof: E::G2Affine,
) -> bool {
    // Let tau be the SRS secret. Then the KZG Verifier check in G2 looks like:
    // e(G, C - eval * H) == e(tau * G - point * G, W)
    // We rewrite this as a multi-pairing:
    // e(-G, C - eval * H) * e(tau * G - point * G, W) == 1
    let lhs_g2 = (commitment.into_group() - vk.h * eval).into_affine();
    let rhs_g1 = (tau_g.into_group() - vk.g * point).into_affine();
    E::multi_pairing([vk.neg_g, rhs_g1], [lhs_g2, evaluation_proof]).is_zero()
}

/// Similar to verify_kzg_opening_g2 but for g1.
#[must_use]
pub fn verify_right_key<E: Pairing>(
    vk: &VerifierKey<E>,
    // These are the KZG commitments that arise from the final commitment key
    &RightKey { w_1, w_2 }: &RightKey<E>,
    &EvaluationProof(proof_1, proof_2): &EvaluationProof<E::G1Affine>,
    challenges: &[E::ScalarField],
    twist_inv: E::ScalarField,
    point: E::ScalarField,
) -> bool {
    // compute f(z) and z^n and then combine into f_w(z) = z^n * f(z)
    let w_poly_at_point = evaluate_ipa_polynomial_shifted(challenges, point, twist_inv);
    let [w_1, w_2]: [E::G1Affine; 2] = E::G1::normalize_batch(&[w_1, w_2]).try_into().unwrap();

    // e(C_f * g^{-y}, h) = e(\pi, w1 * h^{-x})
    let check1 = check_right::<E>(vk, vk.h_alpha, point, w_poly_at_point, w_1, proof_1);

    // e(C_f * g^{-y}, h) = e(\pi, w2 * h^{-x})
    let check2 = check_right::<E>(vk, vk.h_beta, point, w_poly_at_point, w_2, proof_2);
    check1 & check2
}

/// Here `tau` is the SRS secret.
#[must_use]
pub(super) fn check_right<E: Pairing>(
    vk: &VerifierKey<E>,
    tau_h: E::G2Affine,
    point: E::ScalarField,
    eval: E::ScalarField,
    commitment: E::G1Affine,
    evaluation_proof: E::G1Affine,
) -> bool {
    // Let tau be the SRS secret. Then the KZG Verifier check in G1 looks like:
    // e(C - eval * G, H) == e(W, tau * H - point * H)
    // We rewrite this as a multi-pairing:
    // e(-G, C - eval * H) == e(tau * G - point * G, W)
    let lhs_g1 = (commitment.into_group() - vk.g * eval).into_affine();
    let rhs_g2 = (tau_h.into_group() - vk.h * point).into_affine();
    E::multi_pairing([lhs_g1, evaluation_proof], [vk.neg_h, rhs_g2]).is_zero()
}
