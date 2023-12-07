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
        pk: &ProverKey<'a, IPC>,
        instance: &Instance<IPC>,
        witness: &Witness<IP>,
    ) -> Result<Proof<IP, IPC>, Error> {
        let (proof, _) = Self::prove_helper(pk, instance, witness)?;
        Ok(proof)
    }

    /// Used as a subroutine in more complex GIPA protocols that
    /// support succinct verification of the final commitment key.
    pub fn prove_helper<'a>(
        pk: &ProverKey<'a, IPC>,
        instance: &Instance<IPC>,
        witness: &Witness<IP>,
    ) -> Result<(Proof<IP, IPC>, Aux<IPC>), Error> {
        let prover_time = start_timer!(|| "GIPA::Prove");
        let Witness { left, right } = witness;
        let Instance {
            size,
            output,
            commitment: mut com,
            twist,
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
            &IP::twisted_inner_product(left, right, *twist)?,
            output,
            "invalid witness"
        );
        debug_assert!(com.cm_ip.is_none());

        let mut left = left.to_vec();
        let mut right = right.to_vec();
        let mut ck = pk.ck.clone();
        let twist_time = start_timer!(|| "Twist");
        if !twist.is_one() {
            cfg_iter_mut!(left)
                .zip(compute_powers(*size, *twist))
                .for_each(|(l, c)| *l *= c);
            ck.twist_in_place(twist.inverse().unwrap());
        }
        end_timer!(twist_time);

        com += IPC::commit_only_ip(&ck, *output)?;

        #[cfg(debug_assertions)]
        if !IPC::verify(&ck, &left, &right, &com)? {
            return Err(Box::new(InnerProductArgumentError::InnerProductInvalid));
        }
        let mut ck = ck.clone();
        let mut commitments = Vec::new();
        let mut challenges: Vec<Scalar<IPC>> = Vec::new();
        assert!(left.len().is_power_of_two());

        let loop_time = start_timer!(|| "Loop");
        let (final_msg, final_ck) = 'recurse: loop {
            let recurse = start_timer!(|| format!("Recurse round size {}", left.len()));
            if left.len() == 1 {
                // base case
                end_timer!(recurse);
                break 'recurse ((left[0], right[0]), ck.to_owned());
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
                let com_1 = IPC::commit(&ck_l, &left_1, &right_1)?;
                end_timer!(cl);
                let cr = start_timer!(|| "Commit R");
                let com_2 = IPC::commit(&ck_r, &left_2, &right_2)?;
                end_timer!(cr);

                // Fiat-Shamir challenge
                let default_transcript = Scalar::<IPC>::one();
                let prev_challenge = challenges.last().unwrap_or(&default_transcript);
                let (c, c_inv) = Self::compute_challenge(prev_challenge, &com_1, &com_2)?;

                // Set up values for next step of recursion
                let rescale_m1 = start_timer!(|| "Rescale M1");
                left = cfg_iter!(left_1)
                    .zip(left_2)
                    .map(|(a_1, a_2)| *a_1 * c + *a_2)
                    .collect::<Vec<IP::LeftMessage>>();
                end_timer!(rescale_m1);

                let rescale_m2 = start_timer!(|| "Rescale M2");
                right = cfg_iter!(right_2)
                    .zip(right_1)
                    .map(|(b_1, b_2)| *b_1 * c_inv + *b_2)
                    .collect::<Vec<IP::RightMessage>>();
                end_timer!(rescale_m2);

                let fold_time = start_timer!(|| "Fold");
                ck = IPCommKey::fold(&ck_l, &ck_r, &c_inv, &c)?;
                end_timer!(fold_time);

                commitments.push((com_1, com_2));
                challenges.push(c);
                end_timer!(recurse);
            }
        };
        end_timer!(loop_time);
        challenges.reverse();
        commitments.reverse();
        end_timer!(prover_time);
        Ok((
            Proof::new(commitments, final_msg),
            Aux::new(challenges, final_ck.try_into().unwrap()),
        ))
    }
}
