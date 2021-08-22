use crate::{
    gipa::GIPAProof,
    tipa::structured_scalar_message::{
        structured_scalar_power, GIPAWithSSM, SSMPlaceholderCommitment,
    },
    Error,
};

use std::marker::PhantomData;

use ark_dh_commitments::{
    afgho16::AFGHOCommitmentG1,
    identity::{HomomorphicPlaceholderValue, IdentityCommitment, IdentityOutput},
    pedersen::PedersenCommitment,
    DoublyHomomorphicCommitment,
};
use ark_ec::PairingEngine;
use ark_ff::{Field, Zero};
use ark_inner_products::{
    ExtensionFieldElement, MultiexponentiationInnerProduct, ScalarInnerProduct,
};
use ark_poly::polynomial::{
    univariate::DensePolynomial as UnivariatePolynomial, Polynomial, UVPolynomial,
};
use ark_std::rand::Rng;
use ark_std::{end_timer, start_timer};
use merlin::Transcript;

const UNIVARIATE_DOMAIN_SEP: &[u8] = b"ip_proofs-Univariate_KZG_transparent";
const BIVARIATE_DOMAIN_SEP: &[u8] = b"ip_proofs-Bivariate_KZG_transparent";

type PolynomialEvaluationSecondTierIPA<P> = GIPAWithSSM<
    MultiexponentiationInnerProduct<<P as PairingEngine>::G1Projective>,
    AFGHOCommitmentG1<P>,
    IdentityCommitment<<P as PairingEngine>::G1Projective, <P as PairingEngine>::Fr>,
>;

type PolynomialEvaluationSecondTierIPAProof<P> = GIPAProof<
    MultiexponentiationInnerProduct<<P as PairingEngine>::G1Projective>,
    AFGHOCommitmentG1<P>,
    SSMPlaceholderCommitment<<P as PairingEngine>::Fr>,
    IdentityCommitment<<P as PairingEngine>::G1Projective, <P as PairingEngine>::Fr>,
>;

type PolynomialEvaluationFirstTierIPA<P> = GIPAWithSSM<
    ScalarInnerProduct<<P as PairingEngine>::Fr>,
    PedersenCommitment<<P as PairingEngine>::G1Projective>,
    IdentityCommitment<<P as PairingEngine>::Fr, <P as PairingEngine>::Fr>,
>;

type PolynomialEvaluationFirstTierIPAProof<P> = GIPAProof<
    ScalarInnerProduct<<P as PairingEngine>::Fr>,
    PedersenCommitment<<P as PairingEngine>::G1Projective>,
    SSMPlaceholderCommitment<<P as PairingEngine>::Fr>,
    IdentityCommitment<<P as PairingEngine>::Fr, <P as PairingEngine>::Fr>,
>;

pub struct BivariatePolynomial<F: Field> {
    y_polynomials: Vec<UnivariatePolynomial<F>>,
}

impl<F: Field> BivariatePolynomial<F> {
    pub fn evaluate(&self, point: &(F, F)) -> F {
        let (x, y) = point;
        let mut point_x_powers = vec![];
        let mut cur = F::one();
        for _ in 0..(self.y_polynomials.len()) {
            point_x_powers.push(cur);
            cur *= x;
        }
        point_x_powers
            .iter()
            .zip(&self.y_polynomials)
            .map(|(x_power, y_polynomial)| x_power.clone() * y_polynomial.evaluate(&y))
            .sum()
    }
}

pub struct OpeningProof<P: PairingEngine> {
    second_tier_ip_proof: PolynomialEvaluationSecondTierIPAProof<P>,
    y_eval_comm: P::G1Projective,
    first_tier_ip_proof: PolynomialEvaluationFirstTierIPAProof<P>,
}

pub struct BivariatePolynomialCommitment<P: PairingEngine> {
    _pairing: PhantomData<P>,
}

impl<P: PairingEngine> BivariatePolynomialCommitment<P> {
    pub fn setup<R: Rng>(
        rng: &mut R,
        x_degree: usize,
        y_degree: usize,
    ) -> Result<(Vec<P::G1Projective>, Vec<P::G2Projective>), Error> {
        let first_tier_ck = PolynomialEvaluationFirstTierIPA::<P>::setup(rng, y_degree + 1)?.0;
        let second_tier_ck = PolynomialEvaluationSecondTierIPA::<P>::setup(rng, x_degree + 1)?.0;
        Ok((first_tier_ck, second_tier_ck))
    }

    pub fn commit(
        ck: &(Vec<P::G1Projective>, Vec<P::G2Projective>),
        bivariate_polynomial: &BivariatePolynomial<P::Fr>,
    ) -> Result<(ExtensionFieldElement<P>, Vec<P::G1Projective>), Error> {
        let (first_tier_ck, second_tier_ck) = ck;
        assert!(second_tier_ck.len() >= bivariate_polynomial.y_polynomials.len());

        // Create first-tier commitments to Y polynomials
        let y_polynomial_coms = bivariate_polynomial
            .y_polynomials
            .iter()
            .chain(vec![UnivariatePolynomial::zero()].iter().cycle())
            .take(second_tier_ck.len())
            .map(|y_polynomial| {
                let mut coeffs = y_polynomial.coeffs.to_vec();
                assert!(first_tier_ck.len() >= coeffs.len());
                coeffs.resize(first_tier_ck.len(), <P::Fr>::zero());
                PedersenCommitment::<<P as PairingEngine>::G1Projective>::commit(
                    first_tier_ck,
                    &coeffs,
                )
            })
            .collect::<Result<Vec<P::G1Projective>, Error>>()?;

        // Create AFGHO commitment to Y polynomial commitments
        Ok((
            AFGHOCommitmentG1::<P>::commit(&second_tier_ck, &y_polynomial_coms)?,
            y_polynomial_coms,
        ))
    }

    pub fn open(
        transcript: &mut Transcript,
        ck: &(Vec<P::G1Projective>, Vec<P::G2Projective>),
        bivariate_polynomial: &BivariatePolynomial<P::Fr>,
        y_polynomial_comms: &Vec<P::G1Projective>,
        point: &(P::Fr, P::Fr),
    ) -> Result<OpeningProof<P>, Error> {
        // Domain-separate this protocol
        transcript.append_message(b"dom-sep", BIVARIATE_DOMAIN_SEP);

        let (x, y) = point;
        let (first_tier_ck, second_tier_ck) = ck;
        assert!(second_tier_ck.len() >= bivariate_polynomial.y_polynomials.len());

        let precomp_time = start_timer!(|| "Computing coefficients and Pedersen commitment");
        let powers_of_x = structured_scalar_power(second_tier_ck.len(), x);

        let coeffs = bivariate_polynomial
            .y_polynomials
            .iter()
            .chain(vec![UnivariatePolynomial::zero()].iter().cycle())
            .take(second_tier_ck.len())
            .map(|y_polynomial| {
                let mut c = y_polynomial.coeffs.to_vec();
                c.resize(first_tier_ck.len(), <P::Fr>::zero());
                c
            })
            .collect::<Vec<Vec<P::Fr>>>();
        let y_eval_coeffs = (0..first_tier_ck.len())
            .map(|j| {
                (0..second_tier_ck.len())
                    .map(|i| powers_of_x[i].clone() * &coeffs[i][j])
                    .sum()
            })
            .collect::<Vec<P::Fr>>();
        let y_eval_comm = PedersenCommitment::<<P as PairingEngine>::G1Projective>::commit(
            first_tier_ck,
            &y_eval_coeffs,
        )?;
        end_timer!(precomp_time);

        let ipa_time = start_timer!(|| "Computing second tier IPA opening proof");
        let second_tier_ip_proof =
            PolynomialEvaluationSecondTierIPA::<P>::prove_with_structured_scalar_message(
                transcript,
                (y_polynomial_comms, &powers_of_x),
                (second_tier_ck, &HomomorphicPlaceholderValue),
            )?;
        end_timer!(ipa_time);
        let first_tier_ipa_time = start_timer!(|| "Computing first tier IPA opening proof");

        let powers_of_y = structured_scalar_power(first_tier_ck.len(), y);
        let first_tier_ip_proof =
            PolynomialEvaluationFirstTierIPA::<P>::prove_with_structured_scalar_message(
                transcript,
                (&y_eval_coeffs, &powers_of_y),
                (first_tier_ck, &HomomorphicPlaceholderValue),
            )?;
        end_timer!(first_tier_ipa_time);

        Ok(OpeningProof {
            second_tier_ip_proof,
            y_eval_comm,
            first_tier_ip_proof,
        })
    }

    pub fn verify(
        transcript: &mut Transcript,
        ck: &(Vec<P::G1Projective>, Vec<P::G2Projective>),
        com: &ExtensionFieldElement<P>,
        point: &(P::Fr, P::Fr),
        eval: &P::Fr,
        proof: &OpeningProof<P>,
    ) -> Result<bool, Error> {
        // Domain-separate this protocol
        transcript.append_message(b"dom-sep", BIVARIATE_DOMAIN_SEP);

        let (first_tier_ck, second_tier_ck) = ck;
        let (x, y) = point;
        let second_tier_ip_proof_valid =
            PolynomialEvaluationSecondTierIPA::<P>::verify_with_structured_scalar_message(
                transcript,
                (second_tier_ck, &HomomorphicPlaceholderValue),
                (com, &IdentityOutput(vec![proof.y_eval_comm.clone()])),
                x,
                &proof.second_tier_ip_proof,
            )?;
        let first_tier_ip_proof_valid =
            PolynomialEvaluationFirstTierIPA::<P>::verify_with_structured_scalar_message(
                transcript,
                (first_tier_ck, &HomomorphicPlaceholderValue),
                (&proof.y_eval_comm, &IdentityOutput(vec![eval.clone()])),
                y,
                &proof.first_tier_ip_proof,
            )?;
        Ok(second_tier_ip_proof_valid && first_tier_ip_proof_valid)
    }
}

pub struct UnivariatePolynomialCommitment<P: PairingEngine> {
    _pairing: PhantomData<P>,
}

impl<P: PairingEngine> UnivariatePolynomialCommitment<P> {
    fn bivariate_degrees(univariate_degree: usize) -> (usize, usize) {
        //(((univariate_degree + 1) as f64).sqrt().ceil() as usize).next_power_of_two() - 1;
        let sqrt = (((univariate_degree + 1) as f64).sqrt().ceil() as usize).next_power_of_two();
        // Skew split between bivariate degrees to account for scalar IPA being less expensive than MIPP
        let skew_factor = if sqrt >= 8 { 4_usize } else { sqrt / 2 };
        (sqrt / skew_factor - 1, sqrt * skew_factor - 1)
    }

    fn parse_bivariate_degrees_from_ck(
        ck: &(Vec<P::G1Projective>, Vec<P::G2Projective>),
    ) -> (usize, usize) {
        let x_degree = ck.1.len() - 1;
        let y_degree = ck.0.len() - 1;
        (x_degree, y_degree)
    }

    fn bivariate_form(
        bivariate_degrees: (usize, usize),
        polynomial: &UnivariatePolynomial<P::Fr>,
    ) -> BivariatePolynomial<P::Fr> {
        let (x_degree, y_degree) = bivariate_degrees;
        let default_zero = vec![P::Fr::zero()];
        let mut coeff_iter = polynomial
            .coeffs
            .iter()
            .chain(default_zero.iter().cycle())
            .take((x_degree + 1) * (y_degree + 1));

        let mut y_polynomials = Vec::new();
        for _ in 0..x_degree + 1 {
            let mut y_polynomial_coeffs = vec![];
            for _ in 0..y_degree + 1 {
                y_polynomial_coeffs.push(Clone::clone(coeff_iter.next().unwrap()))
            }
            y_polynomials.push(UnivariatePolynomial::from_coefficients_slice(
                &y_polynomial_coeffs,
            ));
        }
        BivariatePolynomial { y_polynomials }
    }

    pub fn setup<R: Rng>(
        rng: &mut R,
        degree: usize,
    ) -> Result<(Vec<P::G1Projective>, Vec<P::G2Projective>), Error> {
        let (x_degree, y_degree) = Self::bivariate_degrees(degree);
        BivariatePolynomialCommitment::<P>::setup(rng, x_degree, y_degree)
    }

    pub fn commit(
        ck: &(Vec<P::G1Projective>, Vec<P::G2Projective>),
        polynomial: &UnivariatePolynomial<P::Fr>,
    ) -> Result<(ExtensionFieldElement<P>, Vec<P::G1Projective>), Error> {
        let bivariate_degrees = Self::parse_bivariate_degrees_from_ck(ck);
        BivariatePolynomialCommitment::<P>::commit(
            ck,
            &Self::bivariate_form(bivariate_degrees, polynomial),
        )
    }

    pub fn open(
        transcript: &mut Transcript,
        ck: &(Vec<P::G1Projective>, Vec<P::G2Projective>),
        polynomial: &UnivariatePolynomial<P::Fr>,
        y_polynomial_comms: &Vec<P::G1Projective>,
        point: &P::Fr,
    ) -> Result<OpeningProof<P>, Error> {
        // Domain-separate this protocol
        transcript.append_message(b"dom-sep", UNIVARIATE_DOMAIN_SEP);

        let (x_degree, y_degree) = Self::parse_bivariate_degrees_from_ck(ck);
        let y = point.clone();
        let x = point.pow(&vec![(y_degree + 1) as u64]);
        BivariatePolynomialCommitment::open(
            transcript,
            ck,
            &Self::bivariate_form((x_degree, y_degree), polynomial),
            y_polynomial_comms,
            &(x, y),
        )
    }

    pub fn verify(
        transcript: &mut Transcript,
        ck: &(Vec<P::G1Projective>, Vec<P::G2Projective>),
        com: &ExtensionFieldElement<P>,
        point: &P::Fr,
        eval: &P::Fr,
        proof: &OpeningProof<P>,
    ) -> Result<bool, Error> {
        // Domain-separate this protocol
        transcript.append_message(b"dom-sep", UNIVARIATE_DOMAIN_SEP);

        let (_, y_degree) = Self::parse_bivariate_degrees_from_ck(ck);
        let y = point.clone();
        let x = y.pow(&vec![(y_degree + 1) as u64]);
        BivariatePolynomialCommitment::verify(transcript, ck, com, &(x, y), eval, proof)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_bls12_381::Bls12_381;
    use ark_ec::PairingEngine;
    use ark_ff::UniformRand;

    const BIVARIATE_X_DEGREE: usize = 7;
    const BIVARIATE_Y_DEGREE: usize = 7;
    //const UNIVARIATE_DEGREE: usize = 56;
    const UNIVARIATE_DEGREE: usize = 65535;
    //const UNIVARIATE_DEGREE: usize = 1048575;

    type TestBivariatePolyCommitment = BivariatePolynomialCommitment<Bls12_381>;
    type TestUnivariatePolyCommitment = UnivariatePolynomialCommitment<Bls12_381>;

    #[test]
    fn transparent_bivariate_poly_commit_test() {
        let mut rng = ark_std::test_rng();

        let ck =
            TestBivariatePolyCommitment::setup(&mut rng, BIVARIATE_X_DEGREE, BIVARIATE_Y_DEGREE)
                .unwrap();

        let mut y_polynomials = Vec::new();
        for _ in 0..BIVARIATE_X_DEGREE + 1 {
            let mut y_polynomial_coeffs = vec![];
            for _ in 0..BIVARIATE_Y_DEGREE + 1 {
                y_polynomial_coeffs.push(<Bls12_381 as PairingEngine>::Fr::rand(&mut rng));
            }
            y_polynomials.push(UnivariatePolynomial::from_coefficients_slice(
                &y_polynomial_coeffs,
            ));
        }
        let bivariate_polynomial = BivariatePolynomial { y_polynomials };

        // Commit to polynomial
        let (com, y_polynomial_comms) =
            TestBivariatePolyCommitment::commit(&ck, &bivariate_polynomial).unwrap();

        // Evaluate at challenge point
        let point = (UniformRand::rand(&mut rng), UniformRand::rand(&mut rng));
        let mut open_transcript = Transcript::new(b"Transparent-bivariate-test");
        let eval_proof = TestBivariatePolyCommitment::open(
            &mut open_transcript,
            &ck,
            &bivariate_polynomial,
            &y_polynomial_comms,
            &point,
        )
        .unwrap();
        let eval = bivariate_polynomial.evaluate(&point);

        // Verify proof
        let mut verif_transcript = Transcript::new(b"Transparent-bivariate-test");
        assert!(TestBivariatePolyCommitment::verify(
            &mut verif_transcript,
            &ck,
            &com,
            &point,
            &eval,
            &eval_proof
        )
        .unwrap());
    }

    // `cargo test transparent_univariate_poly_commit_test --release --features print-trace -- --ignored --nocapture`
    #[ignore]
    #[test]
    fn transparent_univariate_poly_commit_test() {
        let mut rng = ark_std::test_rng();

        let ck = TestUnivariatePolyCommitment::setup(&mut rng, UNIVARIATE_DEGREE).unwrap();

        let mut polynomial_coeffs = vec![];
        for _ in 0..UNIVARIATE_DEGREE + 1 {
            polynomial_coeffs.push(<Bls12_381 as PairingEngine>::Fr::rand(&mut rng));
        }
        let polynomial = UnivariatePolynomial::from_coefficients_slice(&polynomial_coeffs);

        // Commit to polynomial
        let (com, y_polynomial_comms) =
            TestUnivariatePolyCommitment::commit(&ck, &polynomial).unwrap();

        // Evaluate at challenge point
        let point = UniformRand::rand(&mut rng);
        let mut open_transcript = Transcript::new(b"Transparent-univariate-test");
        let eval_proof = TestUnivariatePolyCommitment::open(
            &mut open_transcript,
            &ck,
            &polynomial,
            &y_polynomial_comms,
            &point,
        )
        .unwrap();
        let eval = polynomial.evaluate(&point);

        // Verify proof
        let mut verif_transcript = Transcript::new(b"Transparent-univariate-test");
        assert!(TestUnivariatePolyCommitment::verify(
            &mut verif_transcript,
            &ck,
            &com,
            &point,
            &eval,
            &eval_proof
        )
        .unwrap());
    }
}
