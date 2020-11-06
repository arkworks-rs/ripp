//! A crate for inner pairing product arguments/proofs.
#![deny(warnings, unused, missing_docs)]
#![forbid(unsafe_code)]

use ark_ec::{msm::VariableBaseMSM, AffineCurve, PairingEngine, ProjectiveCurve};
use ark_ff::{to_bytes, Field, One, PrimeField, UniformRand};
use digest::Digest;
use rayon::prelude::*;
use std::marker::PhantomData;

/// Fiat-Shamir Rng
pub mod rng;

use rng::FiatShamirRng;

/// SIPP is a inner-pairing product proof that allows a verifier to check an
/// inner-pairing product over `n` elements with only a single pairing.
pub struct SIPP<E: PairingEngine, D: Digest> {
    _engine: PhantomData<E>,
    _digest: PhantomData<D>,
}

/// `Proof` contains the GT elements produced by the prover.
// TODO(psi): why not just make Proof an alias since there's only one field?
pub struct Proof<E: PairingEngine> {
    gt_elems: Vec<(E::Fqk, E::Fqk)>,
}

impl<E: PairingEngine, D: Digest> SIPP<E, D> {
    /// Produce a proof of the inner pairing product.
    pub fn prove(
        a: &[E::G1Affine],
        b: &[E::G2Affine],
        r: &[E::Fr],
        value: E::Fqk,
    ) -> Result<Proof<E>, ()> {
        assert_eq!(a.len(), b.len());
        // Ensure the order of the input vectors is a power of 2
        assert_eq!(a.len().count_ones(), 1);
        let mut length = a.len();
        assert_eq!(length, b.len());
        assert_eq!(length.count_ones(), 1);
        let mut proof_vec = Vec::new();
        // TODO(psi): should we also input a succinct bilinear group description to the rng?
        let mut rng = FiatShamirRng::<D>::from_seed(&to_bytes![a, b, r, value].unwrap());
        let a = a
            .into_par_iter()
            .zip(r)
            .map(|(a, r)| a.mul(*r))
            .collect::<Vec<_>>();
        let mut a = E::G1Projective::batch_normalization_into_affine(&a);
        let mut b = b.to_vec();

        while length != 1 {
            length /= 2;
            let a_l = &a[..length];
            let a_r = &a[length..];

            let b_l = &b[..length];
            let b_r = &b[length..];

            let z_l = product_of_pairings::<E>(a_r, b_l);
            let z_r = product_of_pairings::<E>(a_l, b_r);
            proof_vec.push((z_l, z_r));
            rng.absorb(&to_bytes![z_l, z_r].unwrap());
            let x: E::Fr = u128::rand(&mut rng).into();

            let a_proj = a_l
                .par_iter()
                .zip(a_r)
                .map(|(a_l, a_r)| {
                    let mut temp = a_r.mul(x);
                    temp.add_assign_mixed(a_l);
                    temp
                })
                .collect::<Vec<_>>();
            a = E::G1Projective::batch_normalization_into_affine(&a_proj);

            let x_inv = x.inverse().unwrap();
            let b_proj = b_l
                .par_iter()
                .zip(b_r)
                .map(|(b_l, b_r)| {
                    let mut temp = b_r.mul(x_inv);
                    temp.add_assign_mixed(b_l);
                    temp
                })
                .collect::<Vec<_>>();
            b = E::G2Projective::batch_normalization_into_affine(&b_proj);
        }

        Ok(Proof {
            gt_elems: proof_vec,
        })
    }

    /// Verify an inner-pairing-product proof.
    pub fn verify(
        a: &[E::G1Affine],
        b: &[E::G2Affine],
        r: &[E::Fr],
        claimed_value: E::Fqk,
        proof: &Proof<E>,
    ) -> Result<bool, ()> {
        // Ensure the order of the input vectors is a power of 2
        let length = a.len();
        assert_eq!(length.count_ones(), 1);
        assert!(length >= 2);
        assert_eq!(length, b.len());
        // Ensure there are the correct number of proof elements
        let proof_len = proof.gt_elems.len();
        assert_eq!(proof_len as f32, f32::log2(length as f32));

        // TODO(psi): should we also input a succinct bilinear group description to the rng?
        let mut rng = FiatShamirRng::<D>::from_seed(&to_bytes![a, b, r, claimed_value].unwrap());

        let x_s = proof
            .gt_elems
            .iter()
            .map(|(z_l, z_r)| {
                rng.absorb(&to_bytes![z_l, z_r].unwrap());
                let x: E::Fr = u128::rand(&mut rng).into();
                x
            })
            .collect::<Vec<_>>();

        let mut x_invs = x_s.clone();
        ark_ff::batch_inversion(&mut x_invs);

        let z_prime = claimed_value
            * &proof
                .gt_elems
                .par_iter()
                .zip(&x_s)
                .zip(&x_invs)
                .map(|(((z_l, z_r), x), x_inv)| {
                    z_l.pow(x.into_repr()) * &z_r.pow(x_inv.into_repr())
                })
                .reduce(|| E::Fqk::one(), |a, b| a * &b);

        let mut s: Vec<E::Fr> = vec![E::Fr::one(); length];
        let mut s_invs: Vec<E::Fr> = vec![E::Fr::one(); length];
        // TODO(psi): batch verify
        for (j, (x, x_inv)) in x_s.into_iter().zip(x_invs).enumerate() {
            for i in 0..length {
                if i & (1 << (proof_len - j - 1)) != 0 {
                    s[i] *= &x;
                    s_invs[i] *= &x_inv;
                }
            }
        }

        let s = s
            .into_iter()
            .zip(r)
            .map(|(x, r)| (x * r).into_repr())
            .collect::<Vec<_>>();
        let s_invs = s_invs
            .iter()
            .map(|x_inv| x_inv.into_repr())
            .collect::<Vec<_>>();

        let a_prime = VariableBaseMSM::multi_scalar_mul(&a, &s);
        let b_prime = VariableBaseMSM::multi_scalar_mul(&b, &s_invs);

        let accept = E::pairing(a_prime, b_prime) == z_prime;

        Ok(accept)
    }
}

/// Compute the product of pairings of `r_i * a_i` and `b_i`.
pub fn product_of_pairings_with_coeffs<E: PairingEngine>(
    a: &[E::G1Affine],
    b: &[E::G2Affine],
    r: &[E::Fr],
) -> E::Fqk {
    let a = a
        .into_par_iter()
        .zip(r)
        .map(|(a, r)| a.mul(*r))
        .collect::<Vec<_>>();
    let a = E::G1Projective::batch_normalization_into_affine(&a);
    let elements = a
        .par_iter()
        .zip(b)
        .map(|(a, b)| (E::G1Prepared::from(*a), E::G2Prepared::from(*b)))
        .collect::<Vec<_>>();
    let num_chunks = elements.len() / rayon::current_num_threads();
    let num_chunks = if num_chunks == 0 {
        elements.len()
    } else {
        num_chunks
    };
    let ml_result = elements
        .par_chunks(num_chunks)
        .map(E::miller_loop)
        .product();
    E::final_exponentiation(&ml_result).unwrap()
}

/// Compute the product of pairings of `a` and `b`.
#[must_use]
pub fn product_of_pairings<E: PairingEngine>(a: &[E::G1Affine], b: &[E::G2Affine]) -> E::Fqk {
    let r = vec![E::Fr::one(); a.len()];
    product_of_pairings_with_coeffs::<E>(a, b, &r)
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_bls12_377::{Bls12_377, Fr, G1Projective, G2Projective};
    use blake2::Blake2s;

    #[test]
    fn prove_and_verify_base_case() {
        let mut rng = FiatShamirRng::<Blake2s>::from_seed(&to_bytes![b"falafel"].unwrap());
        let mut a = Vec::with_capacity(32);
        let mut b = Vec::with_capacity(32);
        let mut r = Vec::with_capacity(32);
        for _ in 0..32 {
            a.push(G1Projective::rand(&mut rng).into_affine());
            b.push(G2Projective::rand(&mut rng).into_affine());
            r.push(Fr::rand(&mut rng));
        }

        let z = product_of_pairings_with_coeffs::<Bls12_377>(&a, &b, &r);

        let proof = SIPP::<Bls12_377, Blake2s>::prove(&a, &b, &r, z);
        assert!(proof.is_ok());
        let proof = proof.unwrap();

        let accept = SIPP::<Bls12_377, Blake2s>::verify(&a, &b, &r, z, &proof);
        assert!(accept.is_ok());
        assert!(accept.unwrap());
    }
}
