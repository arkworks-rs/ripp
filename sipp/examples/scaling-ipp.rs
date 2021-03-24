// For benchmarking
use ark_bls12_377::*;
use ark_ec::ProjectiveCurve;
use ark_ff::UniformRand;
use ark_sipp::{rng::FiatShamirRng, SIPP};
use ark_std::rand::seq::SliceRandom;
use blake2::Blake2s;
use std::time::Instant;

type ExampleSIPP = SIPP<Bls12_377, Blake2s>;
use serde::Serialize;

#[derive(Debug, Serialize)]
struct ProfileData {
    size: usize,
    direct: f64,
    prover: f64,
    verifier: f64,
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 2 || args[1] == "-h" || args[1] == "--help" {
        println!("\nHelp: Invoke this as <program> <log_min_num_inputs> <log_max_num_inputs> <path_to_output_dir>\n");
    }
    let min_num_inputs: usize = String::from(args[1].clone())
        .parse()
        .expect("<log_min_num_constraints> should be integer");
    let max_num_inputs: usize = String::from(args[2].clone())
        .parse()
        .expect("<log_max_num_constraints> should be integer");
    let output_directory = String::from(args[3].clone());

    let num_threads: usize =
        std::env::var("RAYON_NUM_THREADS").map_or(rayon::current_num_threads(), |n| {
            n.parse()
                .expect("The environment variable `RAYON_NUM_THREADS` must be an integer")
        });

    let mut rng = FiatShamirRng::<Blake2s>::from_seed(b"falafel");
    let g = G1Projective::rand(&mut rng);
    let h = G2Projective::rand(&mut rng);
    let mut a_s = Vec::new();
    let mut b_s = Vec::new();
    for _ in 0..(1 << max_num_inputs) {
        a_s.push(g.double());
        b_s.push(h.double());
    }
    let mut a_s = ProjectiveCurve::batch_normalization_into_affine(&a_s);
    let mut b_s = ProjectiveCurve::batch_normalization_into_affine(&b_s);
    let r_s = vec![Fr::rand(&mut rng); a_s.len()];

    let output_file = output_directory + &format!("/ipp-{}-threads.csv", num_threads);
    let mut wtr = csv::Writer::from_writer(std::fs::File::create(output_file).unwrap());
    for i in min_num_inputs..=max_num_inputs {
        let m = 1 << i;
        let num_iters = match i {
            1..=5 => 20,
            7..=10 => 10,
            11..=13 => 5,
            _ => 1,
        };
        let a_s = &mut a_s[..m];
        let b_s = &mut b_s[..m];

        let mut direct_time = 0.0;
        let mut prover_time = 0.0;
        let mut verifier_time = 0.0;
        for _ in 0..num_iters {
            a_s.shuffle(&mut rng);
            b_s.shuffle(&mut rng);
            let start = Instant::now();
            let z = ark_sipp::product_of_pairings_with_coeffs::<Bls12_377>(a_s, b_s, &r_s);
            direct_time += (start.elapsed().as_millis() as f64) / 1_000.0;

            let start = Instant::now();
            let proof = ExampleSIPP::prove(a_s, b_s, &r_s, z.clone()).unwrap();
            prover_time += (start.elapsed().as_millis() as f64) / 1_000.0;

            let start = Instant::now();
            assert!(ExampleSIPP::verify(a_s, b_s, &r_s, z.clone(), &proof).unwrap());
            verifier_time += (start.elapsed().as_millis() as f64) / 1_000.0;
        }
        let num_iters = num_iters as f64;
        println!(
            "=== Benchmarking SIPP over Bls12-377 with {} input(s) and {} thread(s) ====",
            m, num_threads,
        );
        println!("Direct time: {:?} seconds", direct_time / num_iters);
        println!("Prover time: {:?} seconds", prover_time / num_iters);
        println!("Verifier time: {:?} seconds", verifier_time / num_iters);
        println!();
        let d = ProfileData {
            size: m,
            direct: direct_time / num_iters,
            prover: prover_time / num_iters,
            verifier: verifier_time / num_iters,
        };
        wtr.serialize(d).unwrap();
    }
    wtr.flush().unwrap();
}
