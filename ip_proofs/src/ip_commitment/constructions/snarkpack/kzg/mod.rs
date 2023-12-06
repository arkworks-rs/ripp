// Adapted from https://github.com/nikkolasg/snarkpack
use ark_ff::Field;
use ark_poly::{univariate::DensePolynomial, DenseUVPolynomial};

mod data_structures;
pub(crate) mod prove;
pub(crate) mod verify;

pub use data_structures::*;

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

pub(super) fn evaluate_ipa_polynomial_shifted<F: Field>(challenges: &[F], z: F, r: F) -> F {
    let e = evaluate_ipa_polynomial(challenges, z, r);
    e * z.pow([2u64.pow(challenges.len() as u32)])
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

fn ipa_polynomial_shifted<F: Field>(challenges: &[F], r: F) -> DensePolynomial<F> {
    let f = ipa_polynomial(challenges, r);
    let shifted_coeffs = [vec![F::ZERO; f.len()], f.coeffs].concat();
    DensePolynomial::from_coefficients_vec(shifted_coeffs)
}

#[cfg(test)]
mod tests {
    use super::evaluate_ipa_polynomial;
    use super::ipa_polynomial;
    use crate::{
        ip_commitment::snarkpack::{
            kzg::{evaluate_ipa_polynomial_shifted, ipa_polynomial_shifted, verify::check_left},
            LeftKey, RightKey,
        },
        tipa::data_structures::{specialize, GenericSRS},
    };
    use ark_bls12_381::{Bls12_381, Fr, G1Projective as G1, G2Projective as G2};
    use ark_ec::{pairing::Pairing, CurveGroup, VariableBaseMSM};
    use ark_ff::Field;
    use ark_poly::{univariate::DensePolynomial, DenseUVPolynomial, Polynomial};
    use ark_std::{test_rng, UniformRand};

    #[test]
    fn ipa_polynomial_consistency() {
        let mut rng = test_rng();
        let r = Fr::rand(&mut rng);
        let evaluation_point = Fr::rand(&mut rng);

        let challenges = vec![Fr::rand(&mut rng), Fr::rand(&mut rng)];
        let poly = ipa_polynomial(&challenges, r);
        let direct_eval = poly.evaluate(&evaluation_point);
        let fast_eval = evaluate_ipa_polynomial(&challenges, evaluation_point, r);
        assert_eq!(fast_eval, direct_eval);
    }

    #[test]
    fn prove_verify_left() {
        use crate::ip_commitment::snarkpack::kzg::prove::prove_left_key;
        use crate::ip_commitment::snarkpack::kzg::verify::verify_left_key;

        let mut rng = test_rng();
        let point = Fr::rand(&mut rng);

        let challenges = vec![Fr::ONE, Fr::rand(&mut rng), Fr::rand(&mut rng)];
        let poly = ipa_polynomial(&challenges, Fr::ONE);
        let eval = poly.evaluate(&point);
        let fast_eval = evaluate_ipa_polynomial(&challenges, point, Fr::ONE);
        assert_eq!(fast_eval, eval);

        let srs = GenericSRS::<Bls12_381>::sample(16, &mut rng);
        let (pk, vk) = specialize(srs, 8);
        let (v_1, v_2) = {
            (
                G2::msm_unchecked(&pk.h_alpha_powers, &poly.coeffs),
                G2::msm_unchecked(&pk.h_beta_powers, &poly.coeffs),
            )
        };

        let proof =
            prove_left_key(&pk.h_alpha_powers, &pk.h_beta_powers, &challenges, point).unwrap();
        assert!(verify_left_key(
            &vk,
            &LeftKey { v_1, v_2 },
            &proof,
            &challenges,
            point
        ));
    }

    #[test]
    fn prove_verify_right() {
        use crate::ip_commitment::snarkpack::kzg::prove::prove_right_key;
        use crate::ip_commitment::snarkpack::kzg::verify::verify_right_key;

        let mut rng = test_rng();
        let twist = Fr::rand(&mut rng);
        let twist_inv = twist.inverse().unwrap();
        let point = Fr::rand(&mut rng);

        let challenges = vec![Fr::ONE, Fr::rand(&mut rng), Fr::rand(&mut rng)];
        let poly = ipa_polynomial_shifted(&challenges, twist_inv);
        let eval = poly.evaluate(&point);
        let fast_eval = evaluate_ipa_polynomial_shifted(&challenges, point, twist_inv);
        assert_eq!(fast_eval, eval);

        let srs = GenericSRS::<Bls12_381>::sample(16, &mut rng);
        let (pk, vk) = specialize(srs, 8);
        let (w_1, w_2) = {
            (
                G1::msm_unchecked(&pk.g_alpha_powers, &poly.coeffs),
                G1::msm_unchecked(&pk.g_beta_powers, &poly.coeffs),
            )
        };

        let proof = prove_right_key(
            &pk.g_alpha_powers,
            &pk.g_beta_powers,
            &challenges,
            twist_inv,
            point,
        )
        .unwrap();
        assert!(verify_right_key(
            &vk,
            &RightKey { w_1, w_2 },
            &proof,
            &challenges,
            twist_inv,
            point
        ));
    }

    #[test]
    fn prove_verify_normal() {
        let mut rng = test_rng();
        let point = Fr::rand(&mut rng);

        let poly = DensePolynomial::rand(5, &mut rng);
        let eval = poly.evaluate(&point);

        let srs = GenericSRS::<Bls12_381>::sample(16, &mut rng);
        let (pk, vk) = specialize(srs, 8);
        let (v_1, v_2) = {
            (
                G2::msm_unchecked(&pk.h_alpha_powers, &poly.coeffs).into_affine(),
                G2::msm_unchecked(&pk.h_beta_powers, &poly.coeffs).into_affine(),
            )
        };
        let witness_poly = &poly / &DensePolynomial::from_coefficients_vec(vec![-point, Fr::ONE]);

        let test_point = Fr::rand(&mut rng);
        // Check that the witness polynomial is correct
        assert_eq!(
            poly.evaluate(&test_point) - eval,
            witness_poly.evaluate(&test_point) * (test_point - point)
        );

        let (proof_1, proof_2) = rayon::join(
            || G2::msm_unchecked(&pk.h_alpha_powers, &witness_poly.coeffs).into_affine(),
            || G2::msm_unchecked(&pk.h_beta_powers, &witness_poly.coeffs).into_affine(),
        );
        let check_1 = check_left(&vk, vk.g_alpha, point, eval, v_1, proof_1);
        let check_2 = check_left(&vk, vk.g_beta, point, eval, v_2, proof_2);
        assert!(check_1);
        assert!(check_2);
    }

    #[test]
    fn test_srs() {
        let mut rng = test_rng();
        let srs = GenericSRS::<Bls12_381>::sample(16, &mut rng);
        let (pk, vk) = specialize(srs, 16);
        let g = pk.g_alpha_powers[0];
        let h = pk.h_alpha_powers[0];
        let alpha_g = pk.g_alpha_powers[1];
        let alpha_h = pk.h_alpha_powers[1];
        let beta_g = pk.g_beta_powers[1];
        let beta_h = pk.h_beta_powers[1];
        for (g, alpha_g) in pk.g_alpha_powers[..=15].iter().zip(&pk.g_alpha_powers[1..]) {
            assert_eq!(
                Bls12_381::pairing(g, alpha_h),
                Bls12_381::pairing(alpha_g, h)
            );
        }
        for (g, beta_g) in pk.g_beta_powers[..=15].iter().zip(&pk.g_beta_powers[1..]) {
            assert_eq!(Bls12_381::pairing(g, beta_h), Bls12_381::pairing(beta_g, h));
        }
        for (h, alpha_h) in pk.h_alpha_powers[..=15].iter().zip(&pk.h_alpha_powers[1..]) {
            assert_eq!(
                Bls12_381::pairing(g, alpha_h),
                Bls12_381::pairing(alpha_g, h)
            );
        }
        for (h, beta_h) in pk.h_beta_powers[..=15].iter().zip(&pk.h_beta_powers[1..]) {
            assert_eq!(Bls12_381::pairing(g, beta_h), Bls12_381::pairing(beta_g, h));
        }
        assert_eq!(
            Bls12_381::pairing(g, vk.h_alpha),
            Bls12_381::pairing(alpha_g, vk.h)
        );
    }
}
