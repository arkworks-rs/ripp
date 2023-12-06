use super::{data_structures::EvaluationProof, ipa_polynomial, ipa_polynomial_shifted};
use crate::Error;

use ark_ec::{AffineRepr, VariableBaseMSM};
use ark_ff::Field;
use ark_poly::{univariate::DensePolynomial, DenseUVPolynomial};

pub(crate) fn prove_left_key<G: AffineRepr>(
    h_alpha_powers: &[G],
    h_beta_powers: &[G],
    challenges: &[G::ScalarField],
    point: G::ScalarField,
) -> Result<EvaluationProof<G>, Error> {
    // f_v
    let left_poly = ipa_polynomial(challenges, G::ScalarField::ONE);

    prove_evaluation(h_alpha_powers, h_beta_powers, left_poly, point)
}

pub(crate) fn prove_right_key<G: AffineRepr>(
    g_alpha_powers: &[G],
    g_beta_powers: &[G],
    challenges: &[G::ScalarField],
    twist_inv: G::ScalarField,
    point: G::ScalarField,
) -> Result<EvaluationProof<G>, Error> {
    // this computes f(X) = \prod (1 + x (rX)^{2^j})
    let right_poly = ipa_polynomial_shifted(challenges, twist_inv);
    prove_evaluation(g_alpha_powers, g_beta_powers, right_poly, point)
}

/// Returns the KZG opening proof for the given commitment key. Specifically, it
/// returns $g^{f(alpha) - f(z) / (alpha - z)}$ for $a$ and $b$.
fn prove_evaluation<G: AffineRepr>(
    alpha_powers: &[G], // h^alpha^i
    beta_powers: &[G],  // h^beta^i
    poly: DensePolynomial<G::ScalarField>,
    point: G::ScalarField,
) -> Result<EvaluationProof<G>, Error> {
    assert_eq!(alpha_powers.len(), beta_powers.len());
    assert_eq!(alpha_powers.len(), poly.coeffs().len());

    // f_v(X) - f_v(z) / (X - z)
    let witness_poly =
        &poly / &DensePolynomial::from_coefficients_vec(vec![-point, G::ScalarField::ONE]);

    assert!(witness_poly.coeffs.len() <= alpha_powers.len());
    assert!(witness_poly.coeffs.len() <= beta_powers.len());

    // we do one proof over h^a and one proof over h^b (or g^a and g^b depending
    // on the curve we are on). that's the extra cost of the commitment scheme
    // used which is compatible with Groth16 CRS instead of the scheme
    // defined in [BMMTV19]
    let (proof_1, proof_2) = rayon::join(
        || G::Group::msm_unchecked(&alpha_powers, &witness_poly.coeffs),
        || G::Group::msm_unchecked(&beta_powers, &witness_poly.coeffs),
    );
    #[cfg(debug_assertions)]
    {
        use ark_poly::Polynomial;
        use ark_std::UniformRand;
        let rand = G::ScalarField::rand(&mut ark_std::test_rng());
        debug_assert_eq!(
            witness_poly.evaluate(&rand) * (rand - point),
            poly.evaluate(&rand) - poly.evaluate(&point)
        );
    }

    Ok(EvaluationProof::new(proof_1, proof_2))
}
