use ark_bls12_377::{
    constraints::PairingVar as BLS12PairingVar, Bls12_377, Fr as BLS12Fr,
    FrParameters as BLS12FrParameters,
};
use ark_bw6_761::{Fr as BW6Fr, BW6_761};
use ark_crypto_primitives::{
    prf::{
        blake2s::{constraints::Blake2sGadget, Blake2s},
        constraints::PRFGadget,
        PRF,
    },
    snark::*,
};
use ark_ec::{AffineCurve, PairingEngine, ProjectiveCurve};
use ark_ff::{
    biginteger::BigInteger, FftParameters, One, PrimeField, ToConstraintField, UniformRand,
};
use ark_groth16::{constraints::*, Groth16, PreparedVerifyingKey, Proof, VerifyingKey};
use ark_r1cs_std::prelude::*;
use ark_relations::r1cs::{ConstraintSynthesizer, ConstraintSystemRef, SynthesisError};

use ark_ip_proofs::applications::groth16_aggregation::{
    aggregate_proofs, setup_inner_product, verify_aggregate_proof,
};

use ark_std::rand::{rngs::StdRng, SeedableRng};
use blake2::Blake2b;
use csv::Writer;

use std::{io::stdout, time::Instant};

#[derive(Default)]
struct SingleBlake2SCircuit {
    input: [u8; 32],
    output: [u8; 32],
}

impl<F: PrimeField> ConstraintSynthesizer<F> for SingleBlake2SCircuit {
    fn generate_constraints(self, cs: ConstraintSystemRef<F>) -> Result<(), SynthesisError> {
        let seed = UInt8::constant_vec(&[0; 32]);
        let input = UInt8::new_witness_vec(cs.clone(), &self.input)?;
        let hash = <Blake2sGadget as PRFGadget<_, F>>::OutputVar::new_variable(
            cs.clone(),
            || Ok(self.output.clone()),
            AllocationMode::Input,
        )?;
        hash.enforce_equal(&<Blake2sGadget as PRFGadget<_, F>>::evaluate(
            &seed, &input,
        )?)?;
        Ok(())
    }
}

#[derive(Clone)]
struct ManyBlake2SCircuit {
    input: Vec<[u8; 32]>,
    output: Vec<[u8; 32]>,
}

impl<F: PrimeField> ConstraintSynthesizer<F> for ManyBlake2SCircuit {
    fn generate_constraints(self, cs: ConstraintSystemRef<F>) -> Result<(), SynthesisError> {
        let seed = UInt8::constant_vec(&[0; 32]);

        for (hash_input, hash_output) in self.input.iter().zip(&self.output) {
            let input = UInt8::new_witness_vec(cs.clone(), hash_input)?;
            let hash = <Blake2sGadget as PRFGadget<_, F>>::OutputVar::new_variable(
                cs.clone(),
                || Ok(hash_output.clone()),
                AllocationMode::Input,
            )?;
            hash.enforce_equal(&<Blake2sGadget as PRFGadget<_, F>>::evaluate(
                &seed, &input,
            )?)?;
        }
        Ok(())
    }
}

#[derive(Clone)]
struct AggregateBlake2SCircuitVerificationCircuit {
    hash_outputs: Vec<[u8; 32]>,
    proofs: Vec<Proof<Bls12_377>>,
    vk: VerifyingKey<Bls12_377>,
}

impl ConstraintSynthesizer<BW6Fr> for AggregateBlake2SCircuitVerificationCircuit {
    fn generate_constraints(self, cs: ConstraintSystemRef<BW6Fr>) -> Result<(), SynthesisError> {
        let input_gadgets = self
            .hash_outputs
            .iter()
            .map(|h| h.to_field_elements())
            .collect::<Option<Vec<Vec<BLS12Fr>>>>()
            .unwrap()
            .iter()
            .map(|h_as_bls_fr| {
                let h_as_bls_fr_bytes = h_as_bls_fr
                    .iter()
                    .map(|bls_fr| {
                        bls_fr
                            .into_repr()
                            .as_ref()
                            .iter()
                            .map(|bls_fr_int| bls_fr_int.to_le_bytes().to_vec())
                            .collect::<Vec<Vec<u8>>>()
                    })
                    .collect::<Vec<Vec<Vec<u8>>>>()
                    .iter()
                    .flatten()
                    .flatten()
                    .cloned()
                    .collect::<Vec<u8>>();
                UInt8::new_input_vec(cs.clone(), &h_as_bls_fr_bytes)
            })
            .collect::<Result<Vec<Vec<UInt8<BW6Fr>>>, SynthesisError>>()?
            // Allocated as BLS12-377 field elements byte representation packed together to BW6-761 field elements
            // Now split BW6-761 byte representation back to iterator over BLS12-377 field element byte representations
            .iter()
            .map(|h_as_bls_fr_bytes| {
                let bls_field_element_size_in_bytes =
                    <BLS12FrParameters as FftParameters>::BigInt::NUM_LIMBS * 8;
                h_as_bls_fr_bytes
                    .chunks(bls_field_element_size_in_bytes)
                    .map(|bls_field_element_chunk| bls_field_element_chunk.to_vec())
                    .collect::<Vec<Vec<UInt8<BW6Fr>>>>()
            })
            .collect::<Vec<Vec<Vec<UInt8<BW6Fr>>>>>();

        let vk_gadget =
            VerifyingKeyVar::<Bls12_377, BLS12PairingVar>::new_constant(cs.clone(), &self.vk)?;

        let proof_gadgets =
            self.proofs
                .iter()
                .map(|proof| {
                    ProofVar::<Bls12_377, BLS12PairingVar>::new_witness(cs.clone(), || {
                        Ok(proof.clone())
                    })
                })
                .collect::<Result<Vec<ProofVar<Bls12_377, BLS12PairingVar>>, SynthesisError>>()?;

        assert_eq!(input_gadgets.len(), proof_gadgets.len());

        for (input_gadget, proof_gadget) in input_gadgets.iter().zip(&proof_gadgets) {
            let input = input_gadget
                .iter()
                .map(|bytes| {
                    bytes
                        .iter()
                        .flat_map(|b| b.to_bits_le().unwrap())
                        .collect::<Vec<_>>()
                })
                .collect::<Vec<_>>();
            let input = BooleanInputVar::new(input);
            <Groth16VerifierGadget<Bls12_377, BLS12PairingVar>>::verify(
                &vk_gadget,
                &input,
                proof_gadget,
            )?
            .enforce_equal(&Boolean::TRUE)?;
        }

        Ok(())
    }
}

struct AggregateBlake2SCircuitVerificationCircuitInput {
    hash_outputs: Vec<[u8; 32]>,
}

impl ToConstraintField<BW6Fr> for AggregateBlake2SCircuitVerificationCircuitInput {
    fn to_field_elements(&self) -> Option<Vec<BW6Fr>> {
        let mut fr_elements: Vec<BW6Fr> = vec![];

        for h_as_bls_fr_bytes in self
            .hash_outputs
            .iter()
            .map(|h| h.to_field_elements())
            .collect::<Option<Vec<Vec<BLS12Fr>>>>()?
            .iter()
            .map(|h_as_bls_fr| {
                h_as_bls_fr
                    .iter()
                    .map(|bls_fr| {
                        bls_fr
                            .into_repr()
                            .as_ref()
                            .iter()
                            .map(|bls_fr_int| bls_fr_int.to_le_bytes().to_vec())
                            .collect::<Vec<Vec<u8>>>()
                    })
                    .collect::<Vec<Vec<Vec<u8>>>>()
                    .iter()
                    .flatten()
                    .flatten()
                    .cloned()
                    .collect::<Vec<u8>>()
            })
        {
            fr_elements.extend_from_slice(&h_as_bls_fr_bytes.to_field_elements()?);
        }

        Some(fr_elements)
    }
}

fn main() {
    let mut args: Vec<String> = std::env::args().collect();
    if args.last().unwrap() == "--bench" {
        args.pop();
    }
    let temp = if args.len() < 2 || args[1] == "-h" || args[1] == "--help" {
        println!("Usage: ``cargo bench --bench groth16_aggregation -- <num_trials> <num_proofs> <bench_rec=(true/false)> <generate_all_proofs=(true/false)> <monolithic_proof=(true/false)>``");
        return;
    } else {
        let bench_recursion = match args.get(3).unwrap_or(&"true".to_string()).as_ref() {
            "true" => true,
            "false" => false,
            _ => panic!("<bench_rec> should be true/false"),
        };
        let generate_all_proofs = match args.get(4).unwrap_or(&"true".to_string()).as_ref() {
            "true" => true,
            "false" => false,
            _ => panic!("<generate_all_proofs> should be true/false"),
        };
        let monolithic_proof = match args.get(4).unwrap_or(&"true".to_string()).as_ref() {
            "true" => true,
            "false" => false,
            _ => panic!("<monolithic_proof> should be true/false"),
        };
        (
            String::from(args[1].clone())
                .parse()
                .expect("<num_trials> should be integer"),
            String::from(args[2].clone())
                .parse()
                .expect("<num_proofs> should be integer"),
            bench_recursion,
            generate_all_proofs,
            monolithic_proof,
        )
    };
    let (num_trials, num_proofs, bench_recursion, generate_all_proofs, generate_monolithic_proof): (usize, usize, _, _, _)  = temp;

    let mut csv_writer = Writer::from_writer(stdout());
    csv_writer
        .write_record(&["trial", "num_proofs", "scheme", "function", "time"])
        .unwrap();
    csv_writer.flush().unwrap();
    let mut start;
    let mut time;
    let mut rng = StdRng::seed_from_u64(0u64);

    // Compute hashes
    let mut hash_inputs = vec![];
    let mut hash_outputs = vec![];
    for _ in 0..num_proofs {
        let hash_input = UniformRand::rand(&mut rng);
        hash_outputs.push(Blake2s::evaluate(&[0; 32], &hash_input).unwrap());
        hash_inputs.push(hash_input);
    }

    {
        // Prove individual proofs
        start = Instant::now();
        let hash_circuit_parameters =
            Groth16::<Bls12_377>::setup(SingleBlake2SCircuit::default(), &mut rng).unwrap();
        time = start.elapsed().as_millis();
        if generate_all_proofs {
            csv_writer
                .write_record(&[
                    1.to_string(),
                    num_proofs.to_string(),
                    "single_circuit".to_string(),
                    "setup".to_string(),
                    time.to_string(),
                ])
                .unwrap();
            csv_writer.flush().unwrap();
        }

        start = Instant::now();
        let mut proofs = vec![];
        for i in 0..num_proofs {
            proofs.push(
                Groth16::<Bls12_377>::prove(
                    &hash_circuit_parameters.0,
                    SingleBlake2SCircuit {
                        input: hash_inputs[i].clone(),
                        output: hash_outputs[i].clone(),
                    },
                    &mut rng,
                )
                .unwrap(),
            );
            if !generate_all_proofs {
                break;
            }
            //assert!(Groth16::<Bls12_377, SingleBlake2SCircuit, [u8; 32]>::verify(&hash_circuit_parameters.1, &hash_outputs[i], &proofs[i]).unwrap());
        }
        time = start.elapsed().as_millis();

        if !generate_all_proofs {
            proofs = vec![proofs[0].clone(); num_proofs];
            hash_inputs = vec![hash_inputs[0]; num_proofs];
            hash_outputs = vec![hash_outputs[0]; num_proofs];
        } else {
            csv_writer
                .write_record(&[
                    1.to_string(),
                    num_proofs.to_string(),
                    "single_circuit".to_string(),
                    "prove".to_string(),
                    time.to_string(),
                ])
                .unwrap();
            csv_writer.flush().unwrap();
        }
        {
            start = Instant::now();
            let result = batch_verify_proof(
                &ark_groth16::prepare_verifying_key(&hash_circuit_parameters.0.vk),
                &hash_outputs
                    .iter()
                    .map(|h| h.to_field_elements())
                    .collect::<Option<Vec<Vec<<Bls12_377 as PairingEngine>::Fr>>>>()
                    .unwrap(),
                &proofs,
            )
            .unwrap();
            time = start.elapsed().as_millis();
            csv_writer
                .write_record(&[
                    1.to_string(),
                    num_proofs.to_string(),
                    "single_circuit".to_string(),
                    "verify".to_string(),
                    time.to_string(),
                ])
                .unwrap();
            csv_writer.flush().unwrap();
            assert!(result);
        }

        // Benchmark aggregation via IPA
        {
            start = Instant::now();
            let srs = setup_inner_product::<Bls12_377, Blake2b, _>(&mut rng, num_proofs).unwrap();
            time = start.elapsed().as_millis();
            csv_writer
                .write_record(&[
                    1.to_string(),
                    num_proofs.to_string(),
                    "ipa".to_string(),
                    "setup".to_string(),
                    time.to_string(),
                ])
                .unwrap();
            csv_writer.flush().unwrap();
            let v_srs = srs.get_verifier_key();

            for i in 1..=num_trials {
                start = Instant::now();
                let aggregate_proof =
                    aggregate_proofs::<Bls12_377, Blake2b>(&srs, &proofs).unwrap();
                time = start.elapsed().as_millis();
                csv_writer
                    .write_record(&[
                        i.to_string(),
                        num_proofs.to_string(),
                        "ipa".to_string(),
                        "aggregate".to_string(),
                        time.to_string(),
                    ])
                    .unwrap();
                csv_writer.flush().unwrap();

                start = Instant::now();
                let result = verify_aggregate_proof(
                    &v_srs,
                    &hash_circuit_parameters.0.vk,
                    &hash_outputs
                        .iter()
                        .map(|h| h.to_field_elements())
                        .collect::<Option<Vec<Vec<<Bls12_377 as PairingEngine>::Fr>>>>()
                        .unwrap(),
                    &aggregate_proof,
                )
                .unwrap();
                time = start.elapsed().as_millis();
                csv_writer
                    .write_record(&[
                        i.to_string(),
                        num_proofs.to_string(),
                        "ipa".to_string(),
                        "verify".to_string(),
                        time.to_string(),
                    ])
                    .unwrap();
                csv_writer.flush().unwrap();
                assert!(result);
            }
        }

        // Benchmark aggregation via one-layer recursion
        if bench_recursion {
            let agg_verification_circuit = AggregateBlake2SCircuitVerificationCircuit {
                hash_outputs: hash_outputs.clone(),
                proofs: proofs.clone(),
                vk: hash_circuit_parameters.0.vk.clone(),
            };
            let agg_verifier_input = AggregateBlake2SCircuitVerificationCircuitInput {
                hash_outputs: hash_outputs.clone(),
            };
            start = Instant::now();
            let agg_verification_circuit_parameters =
                Groth16::<BW6_761>::setup(agg_verification_circuit.clone(), &mut rng).unwrap();
            time = start.elapsed().as_millis();
            csv_writer
                .write_record(&[
                    1.to_string(),
                    num_proofs.to_string(),
                    "olr".to_string(),
                    "setup".to_string(),
                    time.to_string(),
                ])
                .unwrap();
            csv_writer.flush().unwrap();

            for i in 1..=num_trials {
                start = Instant::now();
                let aggregate_proof = Groth16::<BW6_761>::prove(
                    &agg_verification_circuit_parameters.0,
                    agg_verification_circuit.clone(),
                    &mut rng,
                )
                .unwrap();
                time = start.elapsed().as_millis();
                csv_writer
                    .write_record(&[
                        i.to_string(),
                        num_proofs.to_string(),
                        "olr".to_string(),
                        "aggregate".to_string(),
                        time.to_string(),
                    ])
                    .unwrap();
                csv_writer.flush().unwrap();

                start = Instant::now();
                let result = Groth16::<BW6_761>::verify(
                    &agg_verification_circuit_parameters.1,
                    &agg_verifier_input.to_field_elements().unwrap(),
                    &aggregate_proof,
                )
                .unwrap();
                time = start.elapsed().as_millis();
                csv_writer
                    .write_record(&[
                        i.to_string(),
                        num_proofs.to_string(),
                        "olr".to_string(),
                        "verify".to_string(),
                        time.to_string(),
                    ])
                    .unwrap();
                csv_writer.flush().unwrap();
                assert!(result);
            }
        }
    }

    // Benchmark complete circuit
    if generate_monolithic_proof {
        let circuit = ManyBlake2SCircuit {
            input: hash_inputs.clone(),
            output: hash_outputs.clone(),
        };
        start = Instant::now();
        let circuit_parameters = Groth16::<Bls12_377>::setup(circuit.clone(), &mut rng).unwrap();
        time = start.elapsed().as_millis();
        csv_writer
            .write_record(&[
                1.to_string(),
                num_proofs.to_string(),
                "complete_circuit".to_string(),
                "setup".to_string(),
                time.to_string(),
            ])
            .unwrap();
        csv_writer.flush().unwrap();

        start = Instant::now();
        let proof = Groth16::<Bls12_377>::prove(&circuit_parameters.0, circuit, &mut rng).unwrap();
        time = start.elapsed().as_millis();
        csv_writer
            .write_record(&[
                1.to_string(),
                num_proofs.to_string(),
                "complete_circuit".to_string(),
                "prove".to_string(),
                time.to_string(),
            ])
            .unwrap();
        csv_writer.flush().unwrap();

        start = Instant::now();
        let result = Groth16::<Bls12_377>::verify(
            &circuit_parameters.1,
            &hash_outputs
                .iter()
                .flat_map(|h| h.to_field_elements().unwrap())
                .collect::<Vec<BLS12Fr>>(),
            &proof,
        )
        .unwrap();
        time = start.elapsed().as_millis();
        csv_writer
            .write_record(&[
                1.to_string(),
                num_proofs.to_string(),
                "complete_circuit".to_string(),
                "verify".to_string(),
                time.to_string(),
            ])
            .unwrap();
        csv_writer.flush().unwrap();
        assert!(result);
    }
}

pub fn batch_verify_proof<E: PairingEngine>(
    pvk: &PreparedVerifyingKey<E>,
    public_inputs: &[Vec<E::Fr>],
    proofs: &[Proof<E>],
) -> Result<bool, SynthesisError> {
    let mut rng = StdRng::seed_from_u64(0u64);
    let mut r_powers = Vec::with_capacity(proofs.len());
    for _ in 0..proofs.len() {
        let challenge: E::Fr = u128::rand(&mut rng).into();
        r_powers.push(challenge);
    }

    let combined_inputs = public_inputs
        .iter()
        .zip(&r_powers)
        .map(|(input, r)| {
            let mut g_ic = pvk.vk.gamma_abc_g1[0].into_projective();
            for (i, b) in input.iter().zip(pvk.vk.gamma_abc_g1.iter().skip(1)) {
                g_ic += &b.mul(i.into_repr());
            }
            g_ic.mul(r.into_repr())
        })
        .sum::<E::G1Projective>()
        .into_affine();

    let combined_proof_a_s = proofs
        .iter()
        .zip(&r_powers)
        .map(|(proof, r)| proof.a.mul(*r))
        .collect::<Vec<_>>();
    let combined_proof_a_s = E::G1Projective::batch_normalization_into_affine(&combined_proof_a_s);
    let ml_inputs = proofs
        .iter()
        .zip(&combined_proof_a_s)
        .map(|(proof, a)| ((*a).into(), proof.b.into()))
        .collect::<Vec<_>>();
    let a_r_times_b = E::miller_loop(ml_inputs.iter());

    let combined_c_s = proofs
        .iter()
        .zip(&r_powers)
        .map(|(proof, r)| proof.c.mul(*r))
        .sum::<E::G1Projective>()
        .into_affine();

    let sum_of_rs = (&r_powers).iter().copied().sum::<E::Fr>();
    let combined_alpha = (-pvk.vk.alpha_g1.mul(sum_of_rs)).into_affine();
    let qap = E::miller_loop(
        [
            (combined_alpha.into(), pvk.vk.beta_g2.into()),
            (combined_inputs.into(), pvk.gamma_g2_neg_pc.clone()),
            (combined_c_s.into(), pvk.delta_g2_neg_pc.clone()),
        ]
        .iter(),
    );

    let test =
        E::final_exponentiation(&(qap * &a_r_times_b)).ok_or(SynthesisError::UnexpectedIdentity)?;

    Ok(test == E::Fqk::one())
}
