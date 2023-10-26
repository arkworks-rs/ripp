//! A crate for inner pairing product arguments/proofs.
#![deny(warnings, unused, missing_docs)]
#![forbid(unsafe_code)]

use std::marker::PhantomData;

use ark_ec::{
    pairing::{MillerLoopOutput, Pairing, PairingOutput},
    scalar_mul::variable_base::VariableBaseMSM,
    CurveGroup,
};
use ark_ff::{Field, One, UniformRand};
use ark_serialize::CanonicalSerialize;
use ark_std::Zero;
use digest::{generic_array::typenum::U32, Digest};
use rayon::prelude::*;

/// Fiat-Shamir Rng
pub mod rng;

use rng::FiatShamirRng;

/// SIPP is a inner-pairing product proof that allows a verifier to check an
/// inner-pairing product over `n` elements with only a single pairing.
pub struct SIPP<E: Pairing, D: Digest> {
    _engine: PhantomData<E>,
    _digest: PhantomData<D>,
}

/// `Proof` contains the GT elements produced by the prover.
// TODO(psi): why not just make Proof an alias since there's only one field?
pub struct Proof<E: Pairing> {
    gt_elems: Vec<(PairingOutput<E>, PairingOutput<E>)>,
}

impl<E, D> SIPP<E, D>
where
    E: Pairing,
    D: Digest<OutputSize = U32>,
{
    /// Produce a proof of the inner pairing product.
    pub fn prove(
        a: &[E::G1Affine],
        b: &[E::G2Affine],
        r: &[E::ScalarField],
        value: PairingOutput<E>,
    ) -> Result<Proof<E>, ()> {
        assert_eq!(a.len(), b.len());
        // Ensure the order of the input vectors is a power of 2
        assert_eq!(a.len().count_ones(), 1);
        let mut length = a.len();
        assert_eq!(length, b.len());
        assert_eq!(length.count_ones(), 1);
        let mut proof_vec = Vec::new();
        // TODO(psi): should we also input a succinct bilinear group description to the rng?
        let mut rng = {
            let mut seed = Vec::new();
            (a, b, r, value).serialize_uncompressed(&mut seed).unwrap();
            FiatShamirRng::<D>::from_seed(&seed)
        };
        let a = a
            .into_par_iter()
            .zip(r)
            .map(|(&a, r)| a * r)
            .collect::<Vec<_>>();
        let mut a = E::G1::normalize_batch(&a);
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
            {
                let mut buf = Vec::new();
                (z_l, z_r).serialize_uncompressed(&mut buf).unwrap();
                rng.absorb(&buf);
            }
            let x: E::ScalarField = u128::rand(&mut rng).into();

            let a_proj = a_l
                .par_iter()
                .zip(a_r)
                .map(|(a_l, &a_r)| a_r * x + a_l)
                .collect::<Vec<_>>();
            a = E::G1::normalize_batch(&a_proj);

            let x_inv = x.inverse().unwrap();
            let b_proj = b_l
                .par_iter()
                .zip(b_r)
                .map(|(b_l, &b_r)| b_r * x_inv + b_l)
                .collect::<Vec<_>>();
            b = E::G2::normalize_batch(&b_proj);
        }

        Ok(Proof {
            gt_elems: proof_vec,
        })
    }

    /// Verify an inner-pairing-product proof.
    pub fn verify(
        a: &[E::G1Affine],
        b: &[E::G2Affine],
        r: &[E::ScalarField],
        claimed_value: PairingOutput<E>,
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
        let mut rng = {
            let mut seed = Vec::new();
            (a, b, r, claimed_value)
                .serialize_uncompressed(&mut seed)
                .unwrap();
            FiatShamirRng::<D>::from_seed(&seed)
        };

        let x_s = proof
            .gt_elems
            .iter()
            .map(|(z_l, z_r)| {
                {
                    let mut buf = Vec::new();
                    (*z_l, *z_r).serialize_uncompressed(&mut buf).unwrap();
                    rng.absorb(&buf);
                }
                let x: E::ScalarField = u128::rand(&mut rng).into();
                x
            })
            .collect::<Vec<_>>();

        let mut x_invs = x_s.clone();
        ark_ff::batch_inversion(&mut x_invs);

        let z_prime = claimed_value
            + proof
                .gt_elems
                .par_iter()
                .zip(&x_s)
                .zip(&x_invs)
                .map(|(((z_l, z_r), x), x_inv)| (*z_l * x) + (*z_r * x_inv))
                .reduce(|| PairingOutput::<E>::zero(), |a, b| a + b);

        let mut s: Vec<E::ScalarField> = vec![E::ScalarField::one(); length];
        let mut s_invs: Vec<E::ScalarField> = vec![E::ScalarField::one(); length];
        // TODO(psi): batch verify
        for (j, (x, x_inv)) in x_s.into_iter().zip(x_invs).enumerate() {
            for i in 0..length {
                if i & (1 << (proof_len - j - 1)) != 0 {
                    s[i] *= &x;
                    s_invs[i] *= &x_inv;
                }
            }
        }

        let s = s.into_iter().zip(r).map(|(x, r)| x * r).collect::<Vec<_>>();

        let a_prime = E::G1::msm(&a, &s).unwrap();
        let b_prime = E::G2::msm(&b, &s_invs).unwrap();

        let accept = E::pairing(a_prime, b_prime) == z_prime;

        Ok(accept)
    }
}

/// Compute the product of pairings of `r_i * a_i` and `b_i`.
pub fn product_of_pairings_with_coeffs<E: Pairing>(
    a: &[E::G1Affine],
    b: &[E::G2Affine],
    r: &[E::ScalarField],
) -> PairingOutput<E> {
    let a = a
        .into_par_iter()
        .zip(r)
        .map(|(&a, r)| a * r)
        .collect::<Vec<_>>();
    let a = E::G1::normalize_batch(&a);

    let a = a.par_iter().map(E::G1Prepared::from).collect::<Vec<_>>();
    let b = b.par_iter().map(E::G2Prepared::from).collect::<Vec<_>>();

    // We want to process N chunks in parallel where N is the number of threads available
    let num_chunks = rayon::current_num_threads();
    let chunk_size = if num_chunks <= a.len() {
        a.len() / num_chunks
    } else {
        // More threads than elements. Just do it all in parallel
        1
    };

    // Compute all the (partial) pairings and take the product. We have to take the product over
    // P::TargetField because MillerLoopOutput doesn't impl Product
    let ml_result = a
        .par_chunks(chunk_size)
        .zip(b.par_chunks(chunk_size))
        .map(|(aa, bb)| E::multi_miller_loop(aa.iter().cloned(), bb.iter().cloned()).0)
        .product();

    E::final_exponentiation(MillerLoopOutput(ml_result)).unwrap()
}

/// Compute the product of pairings of `a` and `b`.
#[must_use]
pub fn product_of_pairings<E: Pairing>(a: &[E::G1Affine], b: &[E::G2Affine]) -> PairingOutput<E> {
    let r = vec![E::ScalarField::one(); a.len()];
    product_of_pairings_with_coeffs::<E>(a, b, &r)
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_bls12_377::{Bls12_377, Fr, G1Projective, G2Projective};
    use blake2::Blake2s;

    #[test]
    fn prove_and_verify_base_case() {
        let mut rng = FiatShamirRng::<Blake2s>::from_seed(b"falafel");
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
