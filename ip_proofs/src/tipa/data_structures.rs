use ark_dh_commitments::identity::PlaceholderKey;
use ark_ec::{pairing::Pairing, AffineRepr, CurveGroup};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::borrow::Cow;
use derivative::Derivative;

use crate::{
    gipa::Proof as GIPAProof,
    ip_commitment::snarkpack::{FinalCommKeyProof, LeftKey, RightKey, TIPPCommitment},
    ip_commitment::{FinalIPCommKey, IPCommKey},
};

pub use crate::ip_commitment::snarkpack::GenericSRS;

use super::{IP, IPC};

#[derive(Clone, CanonicalDeserialize, CanonicalSerialize)]
pub struct ProverKey<'b, P: Pairing> {
    pub supported_size: usize,
    pub pk: crate::gipa::ProverKey<'b, TIPPCommitment<P>>,
    pub g_alpha_powers: Vec<P::G1Affine>,
    pub g_beta_powers: Vec<P::G1Affine>,
    pub h_alpha_powers: Vec<P::G2Affine>,
    pub h_beta_powers: Vec<P::G2Affine>,
    pub g_alpha: P::G1Affine,
    pub g_beta: P::G1Affine,
    pub h_alpha: P::G2Affine,
    pub h_beta: P::G2Affine,
}

//TODO: Change SRS to return reference iterator - requires changes to TIPA and GIPA signatures
impl<P: Pairing> ProverKey<'_, P> {
    pub fn vk(&self) -> VerifierKey<P> {
        use core::ops::Neg;
        let g = P::G1Affine::generator();
        let h = P::G2Affine::generator();
        let neg_g = g.into_group().neg().into_affine();
        let neg_h = h.into_group().neg().into_affine();
        let ck_for_ip = self.pk.ck.trim_for_only_ip();
        VerifierKey {
            supported_size: self.supported_size,
            ck_for_ip,
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

#[derive(Clone, CanonicalDeserialize, CanonicalSerialize)]
pub struct VerifierKey<P: Pairing> {
    pub supported_size: usize,
    pub ck_for_ip: IPCommKey<'static, TIPPCommitment<P>>,
    pub g: P::G1Affine,
    pub h: P::G2Affine,
    pub neg_g: P::G1Affine,
    pub neg_h: P::G2Affine,
    pub g_alpha: P::G1Affine,
    pub g_beta: P::G1Affine,
    pub h_alpha: P::G2Affine,
    pub h_beta: P::G2Affine,
}

/// Returns [`ProverKey`] and [`VerifierKey`] specialized to a specific number of
/// proofs.
///
/// # Panics
/// * If the number of proofs is not a power of two, it
/// * If the number of proofs is more than half the SRS size.
pub fn specialize<'a, E: Pairing>(
    srs: GenericSRS<E>,
    num_proofs: usize,
) -> (ProverKey<'a, E>, VerifierKey<E>) {
    assert!(num_proofs.is_power_of_two());
    let supported_size = num_proofs;
    let tn = 2 * num_proofs; // size of the CRS we need
    assert!(srs.g_alpha_powers.len() >= tn);
    assert!(srs.h_alpha_powers.len() >= tn);
    assert!(srs.g_beta_powers.len() >= tn);
    assert!(srs.h_beta_powers.len() >= tn);
    let n = num_proofs;
    // when doing the KZG opening we need _all_ coefficients from 0
    // to 2n-1 because the polynomial is of degree 2n-1.
    let g_low = 0;
    let g_up = tn;
    let h_low = 0;
    let h_up = h_low + n;
    // TODO  precompute window
    let g_alpha_powers = srs.g_alpha_powers[g_low..g_up].to_vec();
    let g_beta_powers = srs.g_beta_powers[g_low..g_up].to_vec();
    let h_alpha_powers = srs.h_alpha_powers[h_low..h_up].to_vec();
    let h_beta_powers = srs.h_beta_powers[h_low..h_up].to_vec();

    let ck_a: Vec<_> = h_alpha_powers
        .iter()
        .zip(&h_beta_powers)
        .map(|(&v_1, &v_2)| LeftKey {
            v_1: v_1.into(),
            v_2: v_2.into(),
        })
        .collect();
    assert_eq!(ck_a.len(), supported_size);

    // however, here we only need the "right" shifted bases for the
    // commitment scheme.
    let ck_b: Vec<_> = g_alpha_powers[n..g_up]
        .iter()
        .zip(&g_beta_powers[n..g_up])
        .map(|(&w_1, &w_2)| RightKey {
            w_1: w_1.into(),
            w_2: w_2.into(),
        })
        .collect();
    assert_eq!(ck_b.len(), supported_size);

    let pk = ProverKey {
        supported_size,
        pk: crate::gipa::ProverKey {
            ck: IPCommKey::new(
                Cow::Owned(ck_a),
                Cow::Owned(ck_b),
                Cow::Owned(PlaceholderKey),
            ),
        },
        g_alpha_powers,
        g_beta_powers,
        h_alpha_powers,
        h_beta_powers,
        g_alpha: srs.g_alpha_powers[1].clone(),
        g_beta: srs.g_beta_powers[1].clone(),
        h_alpha: srs.h_alpha_powers[1].clone(),
        h_beta: srs.h_beta_powers[1].clone(),
    };
    let vk = pk.vk();
    (pk, vk)
}

#[derive(CanonicalSerialize, CanonicalDeserialize, Derivative)]
#[derivative(Clone(bound = "P: Pairing"))]
pub struct Proof<P: Pairing> {
    pub gipa_proof: GIPAProof<IP<P>, IPC<P>>,
    pub final_ck: FinalIPCommKey<IPC<P>>,
    pub final_ck_proof: FinalCommKeyProof<P>,
}
