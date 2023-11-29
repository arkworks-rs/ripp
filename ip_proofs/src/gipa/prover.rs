use ark_ff::Field;
use ark_std::{cfg_iter, cfg_iter_mut, convert::TryInto, end_timer, start_timer};
use digest::Digest;
use num_traits::One;

use crate::{
    ip_commitment::{IPCommKey, IPCommitment, Scalar},
    Error, InnerProductArgumentError,
};
use ark_inner_products::{compute_powers, InnerProduct};

#[cfg(feature = "parallel")]
use rayon::prelude::*;

use super::{data_structures::*, GIPA};

impl<IP, IPC, D> GIPA<IP, IPC, D>
where
    D: Digest,
    IP: InnerProduct,
    IPC: IPCommitment<IP = IP>,
{
    pub fn prove<'a>(
        ck: &IPCommKey<'a, IPC>,
        instance: &Instance<IPC>,
        witness: &Witness<IP>,
    ) -> Result<Proof<IP, IPC, D>, Error> {
        let Witness { left, right } = witness;
        let Instance {
            size,
            output,
            commitment: com,
            random_challenge: c,
        } = instance;

        if !left.len().is_power_of_two() {
            // Power of 2 length
            return Err(Box::new(InnerProductArgumentError::MessageLengthInvalid(
                left.len(),
                right.len(),
            )));
        }

        debug_assert_eq!(
            left.len(),
            right.len(),
            "left and right vectors are of unequal length"
        );
        debug_assert_eq!(
            &IP::twisted_inner_product(left, right, *c)?,
            output,
            "invalid witness"
        );

        let mut left = left.to_vec();
        let mut ck = ck.clone();
        if !c.is_one() {
            let powers_of_c = compute_powers(*size, *c);
            cfg_iter_mut!(left)
                .zip(&powers_of_c)
                .for_each(|(l, c)| *l *= *c);
            ck.twist_in_place(c.inverse().unwrap());
        }

        // if !IPC::verify(&ck, &left, right, &output, com)? {
        //     return Err(Box::new(InnerProductArgumentError::InnerProductInvalid));
        // }

        let (proof, _) = Self::prove_with_aux(&ck, &left, right)?;
        Ok(proof)
    }

    pub fn prove_with_aux(
        ck: &IPCommKey<IPC>,
        left: &[IP::LeftMessage],
        right: &[IP::RightMessage],
    ) -> Result<(Proof<IP, IPC, D>, GIPAAux<IP, IPC, D>), Error> {
        Self::_prove(ck, left.to_vec(), right.to_vec())
    }

    // Returns vector of recursive commitments and transcripts in reverse order
    fn _prove<'a>(
        ck: &IPCommKey<'a, IPC>,
        mut left: Vec<IP::LeftMessage>,
        mut right: Vec<IP::RightMessage>,
    ) -> Result<(Proof<IP, IPC, D>, GIPAAux<IP, IPC, D>), Error> {
        let mut ck = ck.clone();
        let mut r_commitment_steps = Vec::new();
        let mut r_transcript: Vec<Scalar<IPC>> = Vec::new();
        assert!(left.len().is_power_of_two());
        let (m_base, ck_base) = 'recurse: loop {
            let recurse = start_timer!(|| format!("Recurse round size {}", left.len()));
            if left.len() == 1 {
                // base case
                break 'recurse ((left[0].clone(), right[0].clone()), ck.to_owned());
            } else {
                // recursive step
                // Recurse with problem of half size
                let split = left.len() / 2;

                let left_1 = &left[split..];
                let left_2 = &left[..split];
                let right_1 = &right[..split];
                let right_2 = &right[split..];

                let (ck_l, ck_r) = ck.split(split);

                let cl = start_timer!(|| "Commit L");
                let com_1 = IPC::commit(&ck_l, &left_1, &right_1, || {
                    IP::inner_product(left_1, right_1).unwrap()
                })?;
                end_timer!(cl);
                let cr = start_timer!(|| "Commit R");
                let com_2 = IPC::commit(&ck_r, &left_2, &right_2, || {
                    IP::inner_product(left_2, right_2).unwrap()
                })?;
                end_timer!(cr);

                // Fiat-Shamir challenge
                let default_transcript = Default::default();
                let transcript = r_transcript.last().unwrap_or(&default_transcript);
                let (c, c_inv) = Self::compute_challenge(transcript, &com_1, &com_2)?;

                // Set up values for next step of recursion
                let rescale_m1 = start_timer!(|| "Rescale M1");
                left = cfg_iter!(left_1)
                    .map(|a| *a * c)
                    .zip(left_2)
                    .map(|(a_1, a_2)| a_1 + *a_2)
                    .collect::<Vec<IP::LeftMessage>>();
                end_timer!(rescale_m1);

                let rescale_m2 = start_timer!(|| "Rescale M2");
                right = cfg_iter!(right_2)
                    .map(|b| *b * c_inv)
                    .zip(right_1)
                    .map(|(b_1, b_2)| b_1 + *b_2)
                    .collect::<Vec<IP::RightMessage>>();
                end_timer!(rescale_m2);

                ck = IPCommKey::fold(&ck_l, &ck_r, &c_inv, &c)?;

                r_commitment_steps.push((com_1, com_2));
                r_transcript.push(c);
                end_timer!(recurse);
            }
        };
        r_transcript.reverse();
        r_commitment_steps.reverse();
        Ok((
            Proof::new(r_commitment_steps, m_base),
            GIPAAux::new(r_transcript, ck_base.try_into().unwrap()),
        ))
    }
}
