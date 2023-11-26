use ark_ec::{AffineRepr, CurveGroup};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};

/// OpeningProof represents the KZG evaluation proof for the SRS used in our scheme.
#[derive(Clone, Debug, PartialEq, CanonicalSerialize, CanonicalDeserialize)]
pub struct EvaluationProof<G: AffineRepr>(pub G, pub G);

impl<G: AffineRepr> EvaluationProof<G> {
    pub fn new(a: G::Group, b: G::Group) -> Self {
        let s = [a, b];
        let s = G::Group::normalize_batch(&s);

        EvaluationProof(s[0], s[1])
    }
}
