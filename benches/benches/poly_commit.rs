use ark_bls12_381::Bls12_381;
use ark_ec::PairingEngine;
use ark_ff::UniformRand;
use ark_ip_proofs::applications::poly_commit::{
    transparent::UnivariatePolynomialCommitment as TransparentIPA,
    UnivariatePolynomialCommitment as IPA, KZG,
};
use ark_poly::polynomial::{
    univariate::DensePolynomial as UnivariatePolynomial, Polynomial, UVPolynomial,
};

use ark_std::rand::{rngs::StdRng, SeedableRng};
use csv::Writer;

use blake2::Blake2b;
use std::{
    io::stdout,
    time::{Duration, Instant},
};

fn main() {
    let mut args: Vec<String> = std::env::args().collect();
    if args.last().unwrap() == "--bench" {
        args.pop();
    }
    let (num_trials, num_data_points): (usize, usize) = if args.len() < 2
        || args[1] == "-h"
        || args[1] == "--help"
    {
        println!("Usage: ``cargo bench --bench poly_commit -- <num_trials> <num_data_points>``");
        return;
    } else {
        (
            String::from(args[1].clone())
                .parse()
                .expect("<num_trials> should be integer"),
            String::from(args[2].clone())
                .parse()
                .expect("<num_data_points> should be integer"),
        )
    };

    let mut csv_writer = Writer::from_writer(stdout());
    csv_writer
        .write_record(&["trial", "scheme", "function", "degree", "time"])
        .unwrap();
    csv_writer.flush().unwrap();
    let mut start;
    let mut time;

    for degree in (0..num_data_points).map(|i| 4_usize.pow((i + 1) as u32) - 1) {
        // Benchmark KZG
        {
            let mut rng = StdRng::seed_from_u64(0u64);
            start = Instant::now();
            let (g_alpha_powers, v_srs) = KZG::<Bls12_381>::setup(&mut rng, degree).unwrap();
            time = start.elapsed().as_millis();
            csv_writer
                .write_record(&[
                    1.to_string(),
                    "kzg".to_string(),
                    "setup".to_string(),
                    degree.to_string(),
                    time.to_string(),
                ])
                .unwrap();
            csv_writer.flush().unwrap();
            for i in 1..num_trials + 1 {
                let polynomial = UnivariatePolynomial::rand(degree, &mut rng);
                let point = <Bls12_381 as PairingEngine>::Fr::rand(&mut rng);
                let eval = polynomial.evaluate(&point);

                // Commit
                start = Instant::now();
                let com = KZG::<Bls12_381>::commit(&g_alpha_powers, &polynomial).unwrap();
                time = start.elapsed().as_millis();
                csv_writer
                    .write_record(&[
                        i.to_string(),
                        "kzg".to_string(),
                        "commit".to_string(),
                        degree.to_string(),
                        time.to_string(),
                    ])
                    .unwrap();

                // Open
                start = Instant::now();
                let proof = KZG::<Bls12_381>::open(&g_alpha_powers, &polynomial, &point).unwrap();
                time = start.elapsed().as_millis();
                csv_writer
                    .write_record(&[
                        i.to_string(),
                        "kzg".to_string(),
                        "open".to_string(),
                        degree.to_string(),
                        time.to_string(),
                    ])
                    .unwrap();

                // Verify
                std::thread::sleep(Duration::from_millis(5000));
                start = Instant::now();
                for _ in 0..50 {
                    let is_valid =
                        KZG::<Bls12_381>::verify(&v_srs, &com, &point, &eval, &proof).unwrap();
                    assert!(is_valid);
                }
                time = start.elapsed().as_millis() / 50;
                csv_writer
                    .write_record(&[
                        i.to_string(),
                        "kzg".to_string(),
                        "verify".to_string(),
                        degree.to_string(),
                        time.to_string(),
                    ])
                    .unwrap();
                csv_writer.flush().unwrap();
            }
        }

        // Benchmark IPA
        {
            let mut rng = StdRng::seed_from_u64(0u64);
            start = Instant::now();
            let srs = IPA::<Bls12_381, Blake2b>::setup(&mut rng, degree).unwrap();
            let v_srs = srs.0.get_verifier_key();
            time = start.elapsed().as_millis();
            csv_writer
                .write_record(&[
                    1.to_string(),
                    "ipa".to_string(),
                    "setup".to_string(),
                    degree.to_string(),
                    time.to_string(),
                ])
                .unwrap();
            csv_writer.flush().unwrap();
            for i in 1..num_trials + 1 {
                let polynomial = UnivariatePolynomial::rand(degree, &mut rng);
                let point = <Bls12_381 as PairingEngine>::Fr::rand(&mut rng);
                let eval = polynomial.evaluate(&point);

                // Commit
                start = Instant::now();
                let (com, prover_aux) =
                    IPA::<Bls12_381, Blake2b>::commit(&srs, &polynomial).unwrap();
                time = start.elapsed().as_millis();
                csv_writer
                    .write_record(&[
                        i.to_string(),
                        "ipa".to_string(),
                        "commit".to_string(),
                        degree.to_string(),
                        time.to_string(),
                    ])
                    .unwrap();

                // Open
                start = Instant::now();
                let proof = IPA::<Bls12_381, Blake2b>::open(&srs, &polynomial, &prover_aux, &point)
                    .unwrap();
                time = start.elapsed().as_millis();
                csv_writer
                    .write_record(&[
                        i.to_string(),
                        "ipa".to_string(),
                        "open".to_string(),
                        degree.to_string(),
                        time.to_string(),
                    ])
                    .unwrap();

                // Verify
                std::thread::sleep(Duration::from_millis(5000));
                start = Instant::now();
                for _ in 0..50 {
                    let is_valid = IPA::<Bls12_381, Blake2b>::verify(
                        &v_srs, degree, &com, &point, &eval, &proof,
                    )
                    .unwrap();
                    assert!(is_valid);
                }
                time = start.elapsed().as_millis() / 50;
                csv_writer
                    .write_record(&[
                        i.to_string(),
                        "ipa".to_string(),
                        "verify".to_string(),
                        degree.to_string(),
                        time.to_string(),
                    ])
                    .unwrap();
                csv_writer.flush().unwrap();
            }

            // Benchmark transparent IPA
            {
                let mut rng = StdRng::seed_from_u64(0u64);
                start = Instant::now();
                let ck = TransparentIPA::<Bls12_381, Blake2b>::setup(&mut rng, degree).unwrap();
                time = start.elapsed().as_millis();
                csv_writer
                    .write_record(&[
                        1.to_string(),
                        "transparent_ipa".to_string(),
                        "setup".to_string(),
                        degree.to_string(),
                        time.to_string(),
                    ])
                    .unwrap();
                csv_writer.flush().unwrap();
                for i in 1..num_trials + 1 {
                    let polynomial = UnivariatePolynomial::rand(degree, &mut rng);
                    let point = <Bls12_381 as PairingEngine>::Fr::rand(&mut rng);
                    let eval = polynomial.evaluate(&point);

                    // Commit
                    start = Instant::now();
                    let (com, prover_aux) =
                        TransparentIPA::<Bls12_381, Blake2b>::commit(&ck, &polynomial).unwrap();
                    time = start.elapsed().as_millis();
                    csv_writer
                        .write_record(&[
                            i.to_string(),
                            "transparent_ipa".to_string(),
                            "commit".to_string(),
                            degree.to_string(),
                            time.to_string(),
                        ])
                        .unwrap();

                    // Open
                    start = Instant::now();
                    let proof = TransparentIPA::<Bls12_381, Blake2b>::open(
                        &ck,
                        &polynomial,
                        &prover_aux,
                        &point,
                    )
                    .unwrap();
                    time = start.elapsed().as_millis();
                    csv_writer
                        .write_record(&[
                            i.to_string(),
                            "transparent_ipa".to_string(),
                            "open".to_string(),
                            degree.to_string(),
                            time.to_string(),
                        ])
                        .unwrap();

                    // Verify
                    std::thread::sleep(Duration::from_millis(5000));
                    start = Instant::now();
                    for _ in 0..50 {
                        let is_valid = TransparentIPA::<Bls12_381, Blake2b>::verify(
                            &ck, &com, &point, &eval, &proof,
                        )
                        .unwrap();
                        assert!(is_valid);
                    }
                    time = start.elapsed().as_millis() / 50;
                    csv_writer
                        .write_record(&[
                            i.to_string(),
                            "transparent_ipa".to_string(),
                            "verify".to_string(),
                            degree.to_string(),
                            time.to_string(),
                        ])
                        .unwrap();
                    csv_writer.flush().unwrap();
                }
            }
        }
    }
}
