use ark_dh_commitments::identity::PlaceholderKey;
use ark_ec::{pairing::Pairing, scalar_mul::fixed_base::FixedBase, AffineRepr, CurveGroup, Group};
use ark_ff::PrimeField;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::{borrow::Cow, marker::PhantomData, rand::Rng, One, UniformRand};
use derivative::Derivative;
use digest::Digest;

use crate::{
    gipa::Proof as GIPAProof,
    ip_commitment::{FinalIPCommKey, IPCommKey},
    tipa::tipp::{LeftKey, RightKey},
};

use super::{IP, IPC};

/// A generic and reusable powers-of-tau SRS for TIPP
#[derive(Clone, Debug)]
pub struct GenericSRS<E: Pairing> {
    /// $\{a^i G\}_{i=0}^{N}$
    pub g_alpha_powers: Vec<E::G1Affine>,
    /// $\{a^i H\}_{i=0}^{N}$
    pub h_alpha_powers: Vec<E::G2Affine>,
    /// $\{b^i G\}_{i=n}^{N}$
    pub g_beta_powers: Vec<E::G1Affine>,
    /// $\{b^i H\}_{i=0}^{N}$
    pub h_beta_powers: Vec<E::G2Affine>,
}

impl<E: Pairing> GenericSRS<E> {
    pub fn sample(size: usize, mut rng: impl Rng) -> Self {
        let alpha = E::ScalarField::rand(&mut rng);
        let beta = E::ScalarField::rand(&mut rng);
        let g = E::G1::generator();
        let h = E::G2::generator();

        let mut g_alpha_powers = Vec::new();
        let mut g_beta_powers = Vec::new();
        let mut h_alpha_powers = Vec::new();
        let mut h_beta_powers = Vec::new();
        rayon::scope(|s| {
            let alpha = &alpha;
            let h = &h;
            let g = &g;
            let beta = &beta;
            let g_alpha_powers = &mut g_alpha_powers;
            s.spawn(move |_| {
                *g_alpha_powers = structured_generators_scalar_power(2 * size, g, alpha);
            });
            let g_beta_powers = &mut g_beta_powers;
            s.spawn(move |_| {
                *g_beta_powers = structured_generators_scalar_power(2 * size, g, beta);
            });

            let h_alpha_powers = &mut h_alpha_powers;
            s.spawn(move |_| {
                *h_alpha_powers = structured_generators_scalar_power(2 * size, h, alpha);
            });

            let h_beta_powers = &mut h_beta_powers;
            s.spawn(move |_| {
                *h_beta_powers = structured_generators_scalar_power(2 * size, h, beta);
            });
        });

        debug_assert_eq!(h_alpha_powers[0], E::G2Affine::generator());
        debug_assert_eq!(h_beta_powers[0], E::G2Affine::generator());
        debug_assert_eq!(g_alpha_powers[0], E::G1Affine::generator());
        debug_assert_eq!(g_beta_powers[0], E::G1Affine::generator());
        GenericSRS {
            g_alpha_powers,
            g_beta_powers,
            h_alpha_powers,
            h_beta_powers,
        }
    }

    /// Returns [`ProverKey`] and [`VerifierKey`] specialized to a specific number of
    /// proofs.
    ///
    /// # Panics
    /// * If the number of proofs is not a power of two, it
    /// * If the number of proofs is more than half the SRS size.
    pub(crate) fn specialize<'b>(&self, num_proofs: usize) -> (ProverKey<'b, E>, VerifierKey<E>) {
        assert!(num_proofs.is_power_of_two());
        let supported_size = num_proofs;
        let tn = 2 * num_proofs; // size of the CRS we need
        assert!(self.g_alpha_powers.len() >= tn);
        assert!(self.h_alpha_powers.len() >= tn);
        assert!(self.g_beta_powers.len() >= tn);
        assert!(self.h_beta_powers.len() >= tn);
        let n = num_proofs;
        // when doing the KZG opening we need _all_ coefficients from 0
        // to 2n-1 because the polynomial is of degree 2n-1.
        let g_low = 0;
        let g_up = tn;
        let h_low = 0;
        let h_up = h_low + n;
        // TODO  precompute window
        let g_alpha_powers = self.g_alpha_powers[g_low..g_up].to_vec();
        let g_beta_powers = self.g_beta_powers[g_low..g_up].to_vec();
        let h_alpha_powers = self.h_alpha_powers[h_low..h_up].to_vec();
        let h_beta_powers = self.h_beta_powers[h_low..h_up].to_vec();

        println!(
            "\nPROVER SRS -- nun_proofs {}, tn {}, alpha_power_table {}\n",
            num_proofs,
            tn,
            g_alpha_powers.len()
        );

        let v1 = self.h_alpha_powers[h_low..h_up].to_vec();
        let v2 = self.h_beta_powers[h_low..h_up].to_vec();
        let ck_a: Vec<_> = v1
            .into_iter()
            .zip(v2.into_iter())
            .map(|(v_1, v_2)| LeftKey {
                v_1: v_1.into(),
                v_2: v_2.into(),
            })
            .collect();

        // however, here we only need the "right" shifted bases for the
        // commitment scheme.
        let w1 = self.g_alpha_powers[n..g_up].to_vec();
        let w2 = self.g_beta_powers[n..g_up].to_vec();
        let ck_b: Vec<_> = w1
            .into_iter()
            .zip(w2.into_iter())
            .map(|(w_1, w_2)| RightKey {
                w_1: w_1.into(),
                w_2: w_2.into(),
            })
            .collect();

        let pk = ProverKey {
            supported_size,
            ck: IPCommKey::new(
                Cow::Owned(ck_a),
                Cow::Owned(ck_b),
                Cow::Owned(PlaceholderKey),
            ),
            g_alpha_powers,
            g_beta_powers,
            h_alpha_powers,
            h_beta_powers,
            g_alpha: self.g_alpha_powers[0].clone(),
            g_beta: self.g_beta_powers[0].clone(),
            h_alpha: self.h_alpha_powers[0].clone(),
            h_beta: self.h_beta_powers[0].clone(),
        };
        let vk = pk.vk();
        (pk, vk)
    }
}

pub(crate) fn structured_generators_scalar_power<G: CurveGroup>(
    num: usize,
    g: &G,
    s: &G::ScalarField,
) -> Vec<G::Affine> {
    assert!(num > 0);
    let mut powers_of_scalar = Vec::with_capacity(num);
    let mut pow_s = G::ScalarField::one();
    for _ in 0..num {
        powers_of_scalar.push(pow_s);
        pow_s *= s;
    }
    let scalar_bits = G::ScalarField::MODULUS_BIT_SIZE as usize;
    let window_size = FixedBase::get_mul_window_size(num);
    let g_table = FixedBase::get_window_table::<G>(scalar_bits, window_size, g.clone());
    let powers_of_g =
        FixedBase::msm::<G>(scalar_bits, window_size, &g_table, &powers_of_scalar[..]);
    G::normalize_batch(&powers_of_g)
}

#[derive(Clone)]
pub struct ProverKey<'a, P: Pairing> {
    pub supported_size: usize,
    pub ck: IPCommKey<'a, IPC<P>>,
    pub g_alpha_powers: Vec<P::G1Affine>,
    pub g_beta_powers: Vec<P::G1Affine>,
    pub h_alpha_powers: Vec<P::G2Affine>,
    pub h_beta_powers: Vec<P::G2Affine>,
    pub g_alpha: P::G1Affine,
    pub g_beta: P::G1Affine,
    pub h_alpha: P::G2Affine,
    pub h_beta: P::G2Affine,
}

#[derive(Clone)]
pub struct VerifierKey<P: Pairing> {
    pub supported_size: usize,
    pub g: P::G1Affine,
    pub h: P::G2Affine,
    pub neg_g: P::G1Affine,
    pub neg_h: P::G2Affine,
    pub g_alpha: P::G1Affine,
    pub g_beta: P::G1Affine,
    pub h_alpha: P::G2Affine,
    pub h_beta: P::G2Affine,
}

//TODO: Change SRS to return reference iterator - requires changes to TIPA and GIPA signatures
impl<P: Pairing> ProverKey<'_, P> {
    pub fn vk(&self) -> VerifierKey<P> {
        use core::ops::Neg;
        let g = self.g_alpha_powers[0].clone();
        let h = self.h_alpha_powers[0].clone();
        let neg_g = g.into_group().neg().into_affine();
        let neg_h = h.into_group().neg().into_affine();
        VerifierKey {
            supported_size: self.supported_size,
            g,
            h,
            neg_g,
            neg_h,
            g_alpha: self.g_alpha,
            g_beta: self.g_beta,
            h_alpha: self.h_alpha,
            h_beta: self.h_beta,
        }
    }
}

#[derive(CanonicalSerialize, CanonicalDeserialize, Derivative)]
#[derivative(Clone)]
pub struct Proof<P, D>
where
    D: Digest,
    P: Pairing,
{
    pub gipa_proof: GIPAProof<IP<P>, IPC<P>, D>,
    pub final_ck: FinalIPCommKey<IPC<P>>,
    pub final_ck_proof: (P::G2, P::G1),
    pub _pair: PhantomData<P>,
}
