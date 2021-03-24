use ark_ip_proofs::applications::groth16_aggregation::{
    aggregate_proofs, setup_inner_product, verify_aggregate_proof,
};

use std::time::Instant;

use ark_bls12_381::{Bls12_381, Fr};
use ark_crypto_primitives::snark::*;
use ark_ff::UniformRand;
use ark_groth16::Groth16;
use ark_r1cs_std::{fields::fp::FpVar, prelude::*};
use ark_relations::r1cs::{ConstraintSynthesizer, ConstraintSystemRef, SynthesisError};

use ark_std::rand::{rngs::StdRng, SeedableRng};
use blake2::Blake2b;

#[derive(Clone)]
struct TestCircuit {
    public_inputs: Vec<Fr>,
    witness_input: Fr,
    public_sum: Fr,
}

impl ConstraintSynthesizer<Fr> for TestCircuit {
    fn generate_constraints(self, cs: ConstraintSystemRef<Fr>) -> Result<(), SynthesisError> {
        let input_variables =
            Vec::<FpVar<Fr>>::new_input(cs.clone(), || Ok(self.public_inputs.clone()))?;
        let sum = FpVar::new_input(cs.clone(), || Ok(&self.public_sum))?;
        let witness = FpVar::new_witness(cs.clone(), || Ok(&self.witness_input))?;

        let mut computed_sum = witness;
        for x in &input_variables {
            computed_sum += x;
        }

        sum.enforce_equal(&computed_sum)?;

        Ok(())
    }
}

fn main() {
    const NUM_PUBLIC_INPUTS: usize = 4;
    const NUM_PROOFS_TO_AGGREGATE: usize = 1024;
    let mut rng = StdRng::seed_from_u64(0u64);

    // Generate parameters for Groth16
    let test_circuit = TestCircuit {
        public_inputs: vec![Default::default(); NUM_PUBLIC_INPUTS],
        public_sum: Default::default(),
        witness_input: Default::default(),
    };
    let parameters = Groth16::<Bls12_381>::setup(test_circuit, &mut rng).unwrap();

    // Generate parameters for inner product aggregation
    let srs = setup_inner_product::<_, Blake2b, _>(&mut rng, NUM_PROOFS_TO_AGGREGATE).unwrap();

    // Generate proofs
    println!("Generating {} Groth16 proofs...", NUM_PROOFS_TO_AGGREGATE);
    let mut start = Instant::now();
    let mut proofs = Vec::new();
    let mut statements = Vec::new();
    for _ in 0..NUM_PROOFS_TO_AGGREGATE {
        // Generate random inputs to sum together
        let mut public_inputs = Vec::new();
        let mut statement = Vec::new();
        for i in 0..NUM_PUBLIC_INPUTS {
            public_inputs.push(Fr::rand(&mut rng));
            statement.push(public_inputs[i].clone());
        }
        let w = Fr::rand(&mut rng);
        let sum: Fr = w.clone() + &public_inputs.iter().sum();
        statement.push(sum.clone());
        let circuit = TestCircuit {
            public_inputs: public_inputs,
            public_sum: sum,
            witness_input: w,
        };

        let proof = Groth16::<Bls12_381>::prove(&parameters.0, circuit.clone(), &mut rng).unwrap();
        proofs.push(proof);
        statements.push(statement);

        //let result = Groth16::<Bls12_381, TestCircuit, [Fr]>::verify(&parameters.1, &statement, &proof).unwrap();
        //assert!(result);
    }
    let generation_time = start.elapsed().as_millis();

    // Aggregate proofs using inner product proofs
    start = Instant::now();
    println!("Aggregating {} Groth16 proofs...", NUM_PROOFS_TO_AGGREGATE);
    let aggregate_proof = aggregate_proofs::<Bls12_381, Blake2b>(&srs, &proofs).unwrap();
    let prover_time = start.elapsed().as_millis();

    println!("Verifying aggregated proof...");
    start = Instant::now();
    let result = verify_aggregate_proof(
        &srs.get_verifier_key(),
        &parameters.0.vk,
        &statements,
        &aggregate_proof,
    )
    .unwrap();
    let verifier_time = start.elapsed().as_millis();
    assert!(result);

    println!("Proof generation time: {} ms", generation_time);
    println!("Proof aggregation time: {} ms", prover_time);
    println!("Proof verification time: {} ms", verifier_time);
}
