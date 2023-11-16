use std::{
    marker::PhantomData,
    ops::{Add, Mul, MulAssign},
};

use ark_dh_commitments::Error;
use ark_ec::{AffineRepr, CurveGroup, Group, pairing::{Pairing, PairingOutput}};
use ark_ec::scalar_mul::fixed_base::FixedBase;
use ark_ff::{Field, PrimeField, UniformRand};
use ark_inner_products::{PairingInnerProduct, cfg_multi_pairing};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::{One, rand::Rng};
use derivative::Derivative;

use super::{IPCommKey, IPCommitment, LeftMessage, OutputMessage, RightMessage, HomomorphicPlaceholderValue};

/// A generic and reusable powers-of-tau SRS for TIPP
#[derive(Clone, Debug)]
pub struct GenericSRS<E: Pairing> {
    /// $\{g^a^i\}_{i=0}^{N}$
    pub g_alpha_powers: Vec<E::G1Affine>,
    /// $\{h^a^i\}_{i=0}^{N}$
    pub h_alpha_powers: Vec<E::G2Affine>,
    /// $\{g^b^i\}_{i=n}^{N}$
    pub g_beta_powers: Vec<E::G1Affine>,
    /// $\{h^b^i\}_{i=0}^{N}$
    pub h_beta_powers: Vec<E::G2Affine>,
}

/// Generates a SRS of the given size. It must NOT be used in production, only
/// in testing, as this is insecure given we know the secret exponent of the SRS.
pub fn setup_fake_srs<E: Pairing>(mut rng: impl Rng, size: usize) -> GenericSRS<E> {
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

    debug_assert!(h_alpha_powers[0] == E::G2Affine::generator());
    debug_assert!(h_beta_powers[0] == E::G2Affine::generator());
    debug_assert!(g_alpha_powers[0] == E::G1Affine::generator());
    debug_assert!(g_beta_powers[0] == E::G1Affine::generator());
    GenericSRS {
        g_alpha_powers,
        g_beta_powers,
        h_alpha_powers,
        h_beta_powers,
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
        pow_s.mul_assign(s);
    }
    let scalar_bits = G::ScalarField::MODULUS_BIT_SIZE as usize;
    let window_size = FixedBase::get_mul_window_size(num);
    let g_table = FixedBase::get_window_table::<G>(scalar_bits, window_size, g.clone());
    let powers_of_g = FixedBase::msm::<G>(
        //let powers_of_g = msm::fixed_base::multi_scalar_mul::<G>(
        scalar_bits,
        window_size,
        &g_table,
        &powers_of_scalar[..],
    );
    powers_of_g.into_iter().map(|v| v.into_affine()).collect()
}

impl<E: Pairing> GenericSRS<E> {
    /// specializes returns the prover and verifier SRS for a specific number of
    /// proofs to aggregate. The number of proofs MUST BE a power of two, it
    /// panics otherwise. The number of proofs must be inferior to half of the
    /// size of the generic srs otherwise it panics.
    pub(crate) fn specialize<'b>(&self, num_proofs: usize) -> IPCommKey<'b, TIPPCommitment<E>> {
        assert!(num_proofs.is_power_of_two());
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
        let ck_a: Vec<_> = v1.into_iter().zip(v2.into_iter()).map(|(v_1, v_2)| LeftKey { v_1: v_1.into(), v_2: v_2.into() }).collect();

        // however, here we only need the "right" shifted bases for the
        // commitment scheme.
        let w1 = self.g_alpha_powers[n..g_up].to_vec();
        let w2 = self.g_beta_powers[n..g_up].to_vec();
        let ck_b: Vec<_> = w1.into_iter().zip(w2.into_iter()).map(|(w_1, w_2)| RightKey { w_1: w_1.into(), w_2: w_2.into() }).collect();

        IPCommKey {
            ck_a: ck_a.into(),
            ck_b: ck_b.into(),
            ck_t: vec![HomomorphicPlaceholderValue].into(),
        }
    }
}

pub(crate) struct TIPPCommitment<E: Pairing>(PhantomData<E>);

#[derive(Clone, PartialEq, Eq, CanonicalSerialize, CanonicalDeserialize)]
pub(crate) struct TIPPCommOutput<E: Pairing>(PairingOutput<E>, PairingOutput<E>);

impl<E: Pairing> Add for TIPPCommOutput<E> {
    type Output = Self;
    fn add(self, rhs: Self) -> Self::Output {
        TIPPCommOutput(self.0 + rhs.0, self.1 + rhs.1)
    }
}

impl<E: Pairing> Mul<E::ScalarField> for TIPPCommOutput<E> {
    type Output = Self;
    fn mul(self, rhs: E::ScalarField) -> Self::Output {
        TIPPCommOutput(self.0 * rhs, self.1 * rhs)
    }
}

impl<E: Pairing> Default for TIPPCommOutput<E> {
    fn default() -> Self {
        TIPPCommOutput(PairingOutput::default(), PairingOutput::default())
    }
}

impl<E: Pairing> IPCommitment for TIPPCommitment<E> {
    type IP = PairingInnerProduct<E>;
    type LeftKey = LeftKey<E>;
    type RightKey = RightKey<E>;
    type IPKey = HomomorphicPlaceholderValue;
    type Commitment = TIPPCommOutput<E>;

    fn setup<'a>(size: usize, mut rng: impl Rng) -> Result<IPCommKey<'a, Self>, Error> {
        let srs = setup_fake_srs(rng, size);
        Ok(srs.specialize(size))
    }

    fn commit<'a>(
        ck: &IPCommKey<'a, Self>,
        l: &[LeftMessage<Self>],
        r: &[RightMessage<Self>],
        ip: &[OutputMessage<Self>],
    ) -> Result<Self::Commitment, Error> {
        let (v1, v2): (Vec<_>, Vec<_>) = ck.ck_a.iter().map(|LeftKey { v_1, v_2 }| (v_1, v_2)).unzip();
        let (w1, w2): (Vec<_>, Vec<_>) = ck.ck_b.iter().map(|RightKey { w_1, w_2 }| (w_1, w_2)).unzip();
        let com_l1 = cfg_multi_pairing(l, &v1).ok_or("cfg_multi_pairing failed")?;
        let com_l2 = cfg_multi_pairing(l, &v2).ok_or("cfg_multi_pairing failed")?;
        let com_r1 = cfg_multi_pairing(&w1, r).ok_or("cfg_multi_pairing failed")?;
        let com_r2 = cfg_multi_pairing(&w2, r).ok_or("cfg_multi_pairing failed")?;

        Ok(TIPPCommOutput(com_l1 + com_r1, com_l2 + com_r2))
    }
}


pub fn structured_scalar_power<F: Field>(num: usize, s: &F) -> Vec<F> {
    let mut powers = vec![F::one()];
    for i in 1..num {
        powers.push(powers[i - 1] * s);
    }
    powers
}

#[derive(Derivative, CanonicalSerialize, CanonicalDeserialize)]
#[derivative(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct LeftKey<E: Pairing> {
    v_1: E::G2,
    v_2: E::G2,
}

#[derive(Derivative, CanonicalSerialize, CanonicalDeserialize)]
#[derivative(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct RightKey<E: Pairing> {
    w_1: E::G1,
    w_2: E::G1,
}

#[derive(Derivative, CanonicalSerialize, CanonicalDeserialize)]
#[derivative(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct InnerProductKey<E: Pairing>(PhantomData<E>);

/************************************************/
/************************************************/
/************************************************/

impl<E: Pairing> Add<Self> for LeftKey<E> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self {
            v_1: self.v_1 + rhs.v_1,
            v_2: self.v_2 + rhs.v_2,
        }
    }
}

impl<E: Pairing> MulAssign<E::ScalarField> for LeftKey<E> {
    fn mul_assign(&mut self, rhs: E::ScalarField) {
        self.v_1 *= rhs;
        self.v_2 *= rhs;
    }
}

impl<E: Pairing> Mul<E::ScalarField> for LeftKey<E> {
    type Output = Self;
    fn mul(mut self, rhs: E::ScalarField) -> Self::Output {
        self *= rhs;
        self
    }
}

impl<E: Pairing> Add<Self> for RightKey<E> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self {
            w_1: self.w_1 + rhs.w_1,
            w_2: self.w_2 + rhs.w_2,
        }
    }
}

impl<E: Pairing> MulAssign<E::ScalarField> for RightKey<E> {
    fn mul_assign(&mut self, rhs: E::ScalarField) {
        self.w_1 *= rhs;
        self.w_2 *= rhs;
    }
}

impl<E: Pairing> Mul<E::ScalarField> for RightKey<E> {
    type Output = Self;
    fn mul(mut self, rhs: E::ScalarField) -> Self::Output {
        self *= rhs;
        self
    }
}

impl<E: Pairing> Add<Self> for InnerProductKey<E> {
    type Output = Self;

    fn add(self, _rhs: Self) -> Self::Output {
        self
    }
}

impl<E: Pairing> MulAssign<E::ScalarField> for InnerProductKey<E> {
    fn mul_assign(&mut self, _rhs: E::ScalarField) {}
}

impl<E: Pairing> Mul<E::ScalarField> for InnerProductKey<E> {
    type Output = Self;
    fn mul(mut self, rhs: E::ScalarField) -> Self::Output {
        self *= rhs;
        self
    }
}
