use dh_commitments::{
    afgho16::{AFGHOCommitmentG1, AFGHOCommitmentG2},
    identity::IdentityCommitment,
};
use groth16_aggregation::{aggregate_proofs, verify_aggregate_proof};
use inner_products::{ExtensionFieldElement, PairingInnerProduct};
use ip_proofs::tipa::TIPA;

use std::time::Instant;

use algebra::{
    bls12_381::{Bls12_381, Fr},
    UniformRand,
};
use r1cs_core::{ConstraintSynthesizer, ConstraintSystem, SynthesisError};
use r1cs_std::{
    alloc::AllocGadget,
    eq::EqGadget,
    fields::{fp::FpGadget, FieldGadget},
};
use zexe_cp::nizk::{groth16::Groth16, NIZK};

use blake2::Blake2b;
use rand::{rngs::StdRng, SeedableRng};

mod groth16_aggregation;

#[derive(Clone)]
struct TestCircuit {
    public_inputs: Vec<Fr>,
    witness_input: Fr,
    public_sum: Fr,
}

impl ConstraintSynthesizer<Fr> for TestCircuit {
    fn generate_constraints<CS: ConstraintSystem<Fr>>(
        self,
        cs: &mut CS,
    ) -> Result<(), SynthesisError> {
        let input_variables =
            Vec::<FpGadget<Fr>>::alloc_input(&mut cs.ns(|| "public_inputs"), || {
                Ok(self.public_inputs.clone())
            })?;
        let sum = <FpGadget<Fr>>::alloc_input(&mut cs.ns(|| "sum_input"), || Ok(&self.public_sum))?;
        let witness = <FpGadget<Fr>>::alloc(&mut cs.ns(|| "witness"), || Ok(&self.witness_input))?;

        let mut computed_sum = witness;
        for (i, x) in input_variables.iter().enumerate() {
            computed_sum = computed_sum.add(&mut cs.ns(|| format!("comp_sum_{}", i)), x)?;
        }

        sum.enforce_equal(&mut cs.ns(|| "check_sum"), &computed_sum)?;

        Ok(())
    }
}

fn main() {
    const NUM_PUBLIC_INPUTS: usize = 4;
    const NUM_PROOFS_TO_AGGREGATE: usize = 16;
    let mut rng = StdRng::seed_from_u64(0u64);
    let mut generation_time = 0.0;
    let mut prover_time = 0.0;
    let mut verifier_time = 0.0;

    // Generate parameters for Groth16
    let test_circuit = TestCircuit {
        public_inputs: vec![Default::default(); NUM_PUBLIC_INPUTS],
        public_sum: Default::default(),
        witness_input: Default::default(),
    };
    let parameters =
        Groth16::<Bls12_381, TestCircuit, [Fr]>::setup(test_circuit, &mut rng).unwrap();

    // Generate parameters for inner product aggregation
    type IP = PairingInnerProduct<Bls12_381>;
    type GC1 = AFGHOCommitmentG1<Bls12_381>;
    type GC2 = AFGHOCommitmentG2<Bls12_381>;
    type IPC = IdentityCommitment<ExtensionFieldElement<Bls12_381>, Fr>;
    type PairingTIPA = TIPA<IP, GC1, GC2, IPC, Bls12_381, Blake2b>;
    let (srs, _) = PairingTIPA::setup(&mut rng, NUM_PROOFS_TO_AGGREGATE).unwrap();

    // Generate proofs
    println!("Generating {} Groth16 proofs...", NUM_PROOFS_TO_AGGREGATE);
    let start = Instant::now();
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

        let proof = Groth16::<Bls12_381, TestCircuit, [Fr]>::prove(
            &parameters.0,
            circuit.clone(),
            &mut rng,
        )
        .unwrap();
        proofs.push(proof);
        statements.push(statement);

        //let result = Groth16::<Bls12_381, TestCircuit, [Fr]>::verify(&parameters.1, &statement, &proof).unwrap();
        //assert!(result);
    }
    generation_time += (start.elapsed().as_millis() as f64) / 1_000.0;

    // Aggregate proofs using inner product proofs
    let start = Instant::now();
    println!("Aggregating {} Groth16 proofs...", NUM_PROOFS_TO_AGGREGATE);
    let aggregate_proof = aggregate_proofs::<Bls12_381, Blake2b>(&srs, &proofs).unwrap();
    prover_time += (start.elapsed().as_millis() as f64) / 1_000.0;

    println!("Verifying aggregated proof...");
    let start = Instant::now();
    let result = verify_aggregate_proof(
        &srs.get_verifier_key(),
        &parameters.0.vk,
        &statements,
        &aggregate_proof,
    )
    .unwrap();
    verifier_time += (start.elapsed().as_millis() as f64) / 1_000.0;
    assert!(result);

    println!("Proof generation time: {:?} seconds", generation_time);
    println!("Proof aggregation time: {:?} seconds", prover_time);
    println!("Proof verification time: {:?} seconds", verifier_time);
}
