// Adapted from https://github.com/nikkolasg/snarkpack
use ark_ff::Field;

mod data_structures;
mod prove;
mod verify;

use ark_poly::{univariate::DensePolynomial, DenseUVPolynomial};
pub use data_structures::*;
pub(crate) use prove::*;
pub(crate) use verify::*;

/// It returns the evaluation of the polynomial $\prod (1 + x_{l-j}(rX)^{2j}$ at
/// the point z, where transcript contains the reversed order of all challenges (the x).
/// THe challenges must be in reversed order for the correct evaluation of the
/// polynomial in O(logn)
pub(super) fn evaluate_ipa_polynomial<F: Field>(challenges: &[F], z: F, r: F) -> F {
    // this is the term (rz) that will get squared at each step to produce the
    // $(rz)^{2j}$ of the formula
    let mut evaluation_point = z * r;

    let one = F::one();

    let mut res = one;
    for x in challenges {
        res *= one + *x * &evaluation_point;
        evaluation_point.square_in_place();
    }

    res
}

// Compute the coefficients of the polynomial $\prod_{j=0}^{l-1} (1 + x_{l-j}(rX)^{2j})$
// It does this in logarithmic time directly; here is an example with 2
// challenges:
//
//     We wish to compute $(1+x_1ra)(1+x_0(ra)^2) = 1 +  x_1ra + x_0(ra)^2 + x_0x_1(ra)^3$
//     Algorithm: $c_{-1} = [1]$; $c_j = c_{i-1} \| (x_{l-j} * c_{i-1})$; $r = r*r$
//     $c_0 = c_{-1} \| (x_1 * r * c_{-1}) = [1] \| [rx_1] = [1, rx_1]$, $r = r^2$
//     $c_1 = c_0 \| (x_0 * r^2c_0) = [1, rx_1] \| [x_0r^2, x_0x_1r^3] = [1, x_1r, x_0r^2, x_0x_1r^3]$
//     which is equivalent to $f(a) = 1 + x_1ra + x_0(ra)^2 + x_0x_1r^2a^3$
//
// This method expects the coefficients in reverse order so transcript[i] =
// x_{l-j}.
// f(Y) = Y^n * \prod (1 + x_{l-j-1} (r_shiftY^{2^j}))
fn ipa_polynomial<F: Field>(challenges: &[F], r: F) -> DensePolynomial<F> {
    let mut coefficients = vec![F::one()];
    let mut power_2_r = r;

    for (i, x) in challenges.iter().enumerate() {
        let n = coefficients.len();
        if i > 0 {
            power_2_r = power_2_r.square();
        }
        for j in 0..n {
            let coeff = coefficients[j] * &(*x * &power_2_r);
            coefficients.push(coeff);
        }
    }

    DensePolynomial::from_coefficients_vec(coefficients)
}

#[test]
fn ipa_polynomial_consistency() {
    use ark_bls12_381::Fr;
    use ark_ff::UniformRand;
    use ark_poly::Polynomial;
    use ark_std::test_rng;

    let mut rng = test_rng();
    let r = Fr::rand(&mut rng);
    let evaluation_point = Fr::rand(&mut rng);

    let challenges = vec![Fr::rand(&mut rng), Fr::rand(&mut rng)];
    let poly = ipa_polynomial(&challenges, r);
    let direct_eval = poly.evaluate(&evaluation_point);
    let fast_eval = evaluate_ipa_polynomial(&challenges, evaluation_point, r);
    assert_eq!(fast_eval, direct_eval);
}
