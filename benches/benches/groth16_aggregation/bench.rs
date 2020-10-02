use algebra::{
    fields::PrimeField,
    curves::PairingEngine,
    bls12_377::{Bls12_377, Fr as BLS12Fr},
    bw6_761::{BW6_761, Fr as BW6Fr, Fq as BW6Fq},
    UniformRand, ToConstraintField,
};
use zexe_cp::{
    prf::{PRF, constraints::PRFGadget, blake2s::{Blake2s, constraints::Blake2sGadget}},
    nizk::{groth16::{Groth16, constraints::Groth16VerifierGadget}, NIZK, constraints::NIZKVerifierGadget},
};
use r1cs_core::{ConstraintSynthesizer, ConstraintSystemRef, SynthesisError};
use r1cs_std::{
    prelude::*,
    bls12_377::PairingVar as BLS12PairingVar,
};
use groth16::{Proof, VerifyingKey};

use ip_proofs::applications::groth16_aggregation::{
    aggregate_proofs, verify_aggregate_proof, setup_inner_product,
};
use ip_proofs::Error;

use rand::{rngs::StdRng, SeedableRng};
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
        hash.enforce_equal(&<Blake2sGadget as PRFGadget<_, F>>::evaluate(&seed, &input)?)?;
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
    fn generate_constraints(
        self,
        cs: ConstraintSystemRef<BW6Fr>,
    ) -> Result<(), SynthesisError> {
        let input_gadgets = self.hash_outputs.iter()
            .map(|h| UInt8::new_input_vec(cs.clone(), h))
            .collect::<Result<Vec<Vec<UInt8<BW6Fr>>>, SynthesisError>>()?;
        let vk_gadget = <Groth16VerifierGadget<Bls12_377, BLS12PairingVar> as NIZKVerifierGadget<Groth16<Bls12_377, SingleBlake2SCircuit, [u8; 32]>, BW6Fr>>::VerificationKeyVar::new_variable(
            cs.clone(),
            || Ok(self.vk.clone()),
            AllocationMode::Input,
        )?;
        let proof_gadgets = self.proofs.iter()
            .map(|proof| {
                <Groth16VerifierGadget<Bls12_377, BLS12PairingVar> as NIZKVerifierGadget<Groth16<Bls12_377, SingleBlake2SCircuit, [u8; 32]>, BW6Fr>>::ProofVar::new_variable(
                    cs.clone(),
                    || Ok(proof.clone()),
                    AllocationMode::Witness,
                )
            }).collect::<Result<Vec<
            <Groth16VerifierGadget<Bls12_377, BLS12PairingVar> as NIZKVerifierGadget<Groth16<Bls12_377, SingleBlake2SCircuit, [u8; 32]>, BW6Fr>>::ProofVar>, SynthesisError>>()?;

        assert_eq!(input_gadgets.len(), proof_gadgets.len());

        for (input_gadget, proof_gadget) in input_gadgets.iter().zip(&proof_gadgets) {
            <Groth16VerifierGadget<Bls12_377, BLS12PairingVar> as NIZKVerifierGadget<Groth16<Bls12_377, SingleBlake2SCircuit, [u8; 32]>, BW6Fr>>::verify(
                &vk_gadget,
                input_gadget,
                proof_gadget,
            )?.enforce_equal(&Boolean::TRUE)?;
        }

        Ok(())
    }
}

struct AggregateBlake2SCircuitVerificationCircuitInput {
    hash_outputs: Vec<[u8; 32]>,
    vk: VerifyingKey<Bls12_377>,
}

impl ToConstraintField<BW6Fr> for AggregateBlake2SCircuitVerificationCircuitInput {
    fn to_field_elements(&self) -> Result<Vec<BW6Fr>, Error> {
        let mut fr_elements = vec![];
        for h in self.hash_outputs.iter() {
            fr_elements.extend_from_slice(&h.to_field_elements()?);
        }
        let VerifyingKey {
            alpha_g1,
            beta_g2,
            gamma_g2,
            delta_g2,
            gamma_abc_g1,
        } = &self.vk;
        fr_elements.extend_from_slice(&alpha_g1.to_field_elements()?);
        fr_elements.extend_from_slice(&beta_g2.to_field_elements()?);
        fr_elements.extend_from_slice(&gamma_g2.to_field_elements()?);
        fr_elements.extend_from_slice(&delta_g2.to_field_elements()?);
        for g1 in gamma_abc_g1.iter() {
            fr_elements.extend_from_slice(&g1.to_field_elements()?);
        }
        Ok(fr_elements)
    }
}

fn main() {
    let mut args: Vec<String> = std::env::args().collect();
    if args.last().unwrap() == "--bench" {
        args.pop();
    }
    let (num_trials, num_proofs): (usize, usize) = if args.len() < 2 || args[1] == "-h" || args[1] == "--help" {
        println!("Usage: ``cargo bench --bench groth16_aggregation -- <num_trials> <num_proofs>``");
        return
    } else {
        (
            String::from(args[1].clone()).parse().expect("<num_trials> should be integer"),
            String::from(args[2].clone()).parse().expect("<num_proofs> should be integer"),
        )
    };

    let mut csv_writer = Writer::from_writer(stdout());
    csv_writer.write_record(&["trial", "num_proofs", "scheme", "function", "time"]).unwrap();
    csv_writer.flush().unwrap();
    let mut start;
    let mut time;
    let mut rng = StdRng::seed_from_u64(0u64);

    // Prove individual proofs
    start = Instant::now();
    let hash_circuit_parameters = Groth16::<Bls12_377, SingleBlake2SCircuit, [u8; 32]>::setup(
        SingleBlake2SCircuit::default(), &mut rng,
    ).unwrap();
    time = start.elapsed().as_millis();
    csv_writer.write_record(&[1.to_string(), num_proofs.to_string(), "ipa_and_olr".to_string(), "individual_proof_setup".to_string(), time.to_string()]).unwrap();
    csv_writer.flush().unwrap();

    let mut hash_inputs = vec![];
    let mut hash_outputs = vec![];
    for _ in 0..num_proofs {
        let hash_input = UniformRand::rand(&mut rng);
        hash_outputs.push(Blake2s::evaluate(&[0; 32], &hash_input).unwrap());
        hash_inputs.push(hash_input);
    }
    start = Instant::now();
    let mut proofs = vec![];
    for i in 0..num_proofs {
        proofs.push(Groth16::<Bls12_377, SingleBlake2SCircuit, [u8; 32]>::prove(
            &hash_circuit_parameters.0,
            SingleBlake2SCircuit {
                input: hash_inputs[i].clone(),
                output: hash_outputs[i].clone(),
            },
            &mut rng,
        ).unwrap());
        //assert!(Groth16::<Bls12_377, SingleBlake2SCircuit, [u8; 32]>::verify(&hash_circuit_parameters.1, &hash_outputs[i], &proofs[i]).unwrap());
    }
    time = start.elapsed().as_millis();
    csv_writer.write_record(&[1.to_string(), num_proofs.to_string(), "ipa_and_olr".to_string(), "individual_proof_prove".to_string(), time.to_string()]).unwrap();
    csv_writer.flush().unwrap();

    // Benchmark aggregation via IPA
    {
        start = Instant::now();
        let srs = setup_inner_product::<Bls12_377, Blake2b, _>(&mut rng, num_proofs).unwrap();
        time = start.elapsed().as_millis();
        csv_writer.write_record(&[1.to_string(), num_proofs.to_string(), "ipa".to_string(), "setup".to_string(), time.to_string()]).unwrap();
        csv_writer.flush().unwrap();
        let v_srs = srs.get_verifier_key();

        for i in 1..=num_trials {
            start = Instant::now();
            let aggregate_proof = aggregate_proofs::<Bls12_377, Blake2b>(&srs, &proofs).unwrap();
            time = start.elapsed().as_millis();
            csv_writer.write_record(&[i.to_string(), num_proofs.to_string(), "ipa".to_string(), "aggregate".to_string(), time.to_string()]).unwrap();
            csv_writer.flush().unwrap();

            start = Instant::now();
            let result = verify_aggregate_proof(
                &v_srs,
                &hash_circuit_parameters.0.vk,
                &hash_outputs.iter().map(|h| h.to_field_elements()).collect::<Result<Vec<Vec<<Bls12_377 as PairingEngine>::Fr>>, Error>>().unwrap(),
                &aggregate_proof,
            ).unwrap();
            time = start.elapsed().as_millis();
            csv_writer.write_record(&[i.to_string(), num_proofs.to_string(), "ipa".to_string(), "verify".to_string(), time.to_string()]).unwrap();
            csv_writer.flush().unwrap();
            assert!(result);
        }
    }
    // Benchmark aggregation via one-layer recursion
    {
        let agg_verification_circuit = AggregateBlake2SCircuitVerificationCircuit {
            hash_outputs: hash_outputs.clone(),
            proofs: proofs.clone(),
            vk: hash_circuit_parameters.0.vk.clone(),
        };
        let agg_verififier_input = AggregateBlake2SCircuitVerificationCircuitInput {
            hash_outputs: hash_outputs.clone(),
            vk: hash_circuit_parameters.0.vk.clone(),
        };
        start = Instant::now();
        let agg_verification_circuit_parameters = Groth16::<BW6_761, AggregateBlake2SCircuitVerificationCircuit, AggregateBlake2SCircuitVerificationCircuitInput>::setup(
            agg_verification_circuit.clone(), &mut rng,
        ).unwrap();
        time = start.elapsed().as_millis();
        csv_writer.write_record(&[1.to_string(), num_proofs.to_string(), "olr".to_string(), "setup".to_string(), time.to_string()]).unwrap();
        csv_writer.flush().unwrap();

        for i in 1..=num_trials {
            start = Instant::now();
            let aggregate_proof = Groth16::<BW6_761, AggregateBlake2SCircuitVerificationCircuit, AggregateBlake2SCircuitVerificationCircuitInput>::prove(
                &agg_verification_circuit_parameters.0,
                agg_verification_circuit.clone(),
                &mut rng,
            ).unwrap();
            time = start.elapsed().as_millis();
            csv_writer.write_record(&[i.to_string(), num_proofs.to_string(), "olr".to_string(), "aggregate".to_string(), time.to_string()]).unwrap();
            csv_writer.flush().unwrap();

            start = Instant::now();
            let result = Groth16::<BW6_761, AggregateBlake2SCircuitVerificationCircuit, AggregateBlake2SCircuitVerificationCircuitInput>::verify(
                &agg_verification_circuit_parameters.1,
                &agg_verififier_input,
                &aggregate_proof,
            ).unwrap();
            time = start.elapsed().as_millis();
            csv_writer.write_record(&[i.to_string(), num_proofs.to_string(), "olr".to_string(), "verify".to_string(), time.to_string()]).unwrap();
            csv_writer.flush().unwrap();
            assert!(result);
        }
    }
}
