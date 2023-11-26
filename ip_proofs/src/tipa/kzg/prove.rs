use super::{data_structures::EvaluationProof, ipa_polynomial};
use crate::Error;

use ark_ec::{AffineRepr, VariableBaseMSM};
use ark_ff::Field;
use ark_poly::{univariate::DensePolynomial, DenseUVPolynomial};

pub(crate) fn prove_commitment_v<G: AffineRepr>(
    h_alpha_powers: &[G],
    h_beta_powers: &[G],
    challenges: &[G::ScalarField],
    point: G::ScalarField,
) -> Result<EvaluationProof<G>, Error> {
    // f_v
    let vkey_poly = ipa_polynomial(challenges, G::ScalarField::ONE);

    prove_evaluation(h_alpha_powers, h_beta_powers, vkey_poly, point)
}

pub(crate) fn prove_commitment_w<G: AffineRepr>(
    g_alpha_powers: &[G],
    g_beta_powers: &[G],
    challenges: &[G::ScalarField],
    r: G::ScalarField,
    point: G::ScalarField,
) -> Result<EvaluationProof<G>, Error> {
    // this computes f(X) = \prod (1 + x (rX)^{2^j})
    let f = ipa_polynomial(challenges, r);
    // this computes f_w(X) = X^n * f(X) - it simply shifts all coefficients to by n
    let fw_coeffs = [vec![G::ScalarField::ZERO; f.len()], f.coeffs].concat();
    let fw = DensePolynomial::from_coefficients_vec(fw_coeffs);

    prove_evaluation(g_alpha_powers, g_beta_powers, fw, point)
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
    let mut witness_poly =
        &poly / &DensePolynomial::from_coefficients_vec(vec![-point, G::ScalarField::ONE]);

    witness_poly
        .coeffs
        .resize(alpha_powers.len(), <G::ScalarField>::ZERO);

    assert_eq!(witness_poly.coeffs.len(), alpha_powers.len());
    assert_eq!(witness_poly.coeffs.len(), beta_powers.len());

    // we do one proof over h^a and one proof over h^b (or g^a and g^b depending
    // on the curve we are on). that's the extra cost of the commitment scheme
    // used which is compatible with Groth16 CRS instead of the scheme
    // defined in [BMMTV19]
    let (a, b) = rayon::join(
        || G::Group::msm(&alpha_powers, &witness_poly.coeffs).expect("msm for a failed!"),
        || G::Group::msm(&beta_powers, &witness_poly.coeffs).expect("msm for b failed!"),
    );
    Ok(EvaluationProof::new(a, b))
}
