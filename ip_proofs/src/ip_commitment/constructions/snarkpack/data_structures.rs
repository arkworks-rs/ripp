use ark_dh_commitments::identity::PlaceholderKey;
use ark_ec::{
    pairing::{Pairing, PairingOutput},
    scalar_mul::fixed_base::FixedBase,
    AffineRepr, CurveGroup, Group,
};
use ark_ff::PrimeField;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::{
    borrow::Cow,
    marker::PhantomData,
    ops::{Add, AddAssign, Mul, MulAssign},
    rand::Rng,
    One, UniformRand,
};
use derivative::Derivative;

use super::kzg::EvaluationProof;

use crate::ip_commitment::IPCommKey;

#[derive(Clone, Debug, Copy, PartialEq, Eq, CanonicalSerialize, CanonicalDeserialize)]
pub struct TIPPCommOutput<E: Pairing> {
    pub t: PairingOutput<E>,
    pub u: PairingOutput<E>,
}

impl<E: Pairing> Add for TIPPCommOutput<E> {
    type Output = Self;
    fn add(mut self, rhs: Self) -> Self::Output {
        self.t += rhs.t;
        self.u += rhs.u;
        self
    }
}
impl<E: Pairing> AddAssign for TIPPCommOutput<E> {
    fn add_assign(&mut self, rhs: Self) {
        self.t += rhs.t;
        self.u += rhs.u;
    }
}

impl<E: Pairing> Mul<E::ScalarField> for TIPPCommOutput<E> {
    type Output = Self;
    fn mul(mut self, rhs: E::ScalarField) -> Self::Output {
        self.t *= rhs;
        self.u *= rhs;
        self
    }
}

impl<E: Pairing> Default for TIPPCommOutput<E> {
    fn default() -> Self {
        let one = E::TargetField::one();
        Self {
            t: PairingOutput(one),
            u: PairingOutput(one),
        }
    }
}

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

    pub fn specialize_to_ck<'b>(
        &self,
        num_proofs: usize,
    ) -> IPCommKey<'b, super::TIPPCommitment<E>> {
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

        IPCommKey::new(
            Cow::Owned(ck_a),
            Cow::Owned(ck_b),
            Cow::Owned(PlaceholderKey),
        )
    }
}

pub(crate) fn structured_generators_scalar_power<G: CurveGroup>(
    num: usize,
    g: G,
    s: G::ScalarField,
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
    let g_table = FixedBase::get_window_table::<G>(scalar_bits, window_size, g);
    let powers_of_g =
        FixedBase::msm::<G>(scalar_bits, window_size, &g_table, &powers_of_scalar[..]);
    G::normalize_batch(&powers_of_g)
}

/************************************************/
/************************************************/
/************************************************/

#[derive(Derivative, CanonicalSerialize, CanonicalDeserialize)]
#[derivative(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct LeftKey<E: Pairing> {
    pub v_1: E::G2,
    pub v_2: E::G2,
}

#[derive(Derivative, CanonicalSerialize, CanonicalDeserialize)]
#[derivative(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct RightKey<E: Pairing> {
    pub w_1: E::G1,
    pub w_2: E::G1,
}

#[derive(Derivative, CanonicalSerialize, CanonicalDeserialize)]
#[derivative(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct InnerProductKey<E: Pairing>(PhantomData<E>);

/// Proof of correctness of the final commitment key.
#[derive(Derivative, CanonicalSerialize, CanonicalDeserialize)]
#[derivative(
    Clone(bound = "E: Pairing"),
    Copy(bound = "E: Pairing"),
    PartialEq(bound = "E: Pairing"),
    Eq(bound = "E: Pairing"),
    Debug(bound = "E: Pairing")
)]
pub struct FinalCommKeyProof<E: Pairing> {
    pub left_proof: EvaluationProof<E::G2Affine>,
    pub right_proof: EvaluationProof<E::G1Affine>,
}

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
