use std::{collections::BinaryHeap, path::Path};

use bio::{Pwm, Sequence};
use itertools::Itertools;

mod bio;
mod fasta_reader;

struct RankedPwm(Pwm, f64, Vec<usize>);

impl PartialEq for RankedPwm {
    fn eq(&self, other: &Self) -> bool {
        self.1 == other.1
    }
}

impl Eq for RankedPwm {}

impl PartialOrd for RankedPwm {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for RankedPwm {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.1.partial_cmp(&other.1).unwrap()
    }
}

fn main() {
    let path = Path::new("assets/test_input.fa");
    let seqs = fasta_reader::read_fasta(path);
    // .into_iter()
    // .take(500)
    // .collect::<Vec<_>>();

    println!("Input file: {:?} read", path);
    println!("Sequences: {}", seqs.len());

    let mut priority_queue = BinaryHeap::new();

    let only_take_top_score = 1000;
    let longer_than = 6;

    get_all_shift_pwms(&seqs[0], &seqs[1], longer_than)
        .into_iter()
        .for_each(|pwm| {
            let score = pwm.get_custom_score();
            priority_queue.push(RankedPwm(pwm, score, vec![0, 1]));
        });

    // for RankedPwm(pwm, score, indicies) in priority_queue.into_iter().take(3) {
    //     println!("Score: {}", score);
    //     println!("{:?}", pwm.history);
    // }
    //

    let mut top_results = Vec::new();

    while let Some(RankedPwm(pwm, score, indicies)) = priority_queue.pop() {
        println!("Score: {:.2}, Depth: {}", score, indicies.len());

        if indicies.len() >= seqs.len() {
            top_results.push((pwm, score));

            if top_results.len() >= only_take_top_score {
                break;
            }

            continue;
        }

        let next_seq = &seqs[indicies.len()];

        get_all_shift_pwms_with_pwm(&pwm, next_seq, longer_than)
            .into_iter()
            .for_each(|new_pwm| {
                let new_score = new_pwm.get_custom_score();
                let mut new_indicies = indicies.clone();
                new_indicies.push(indicies.len());
                priority_queue.push(RankedPwm(new_pwm, new_score, new_indicies));
            });
    }

    for (pwm, score) in top_results
        .into_iter()
        .sorted_by(|a, b| a.1.partial_cmp(&b.1).unwrap())
        .take(5)
    {
        println!("Score: {:.2}", score);
        println!("{:?}", pwm.history);
    }
}

fn get_all_shift_pwms_with_pwm<'a>(pwm: &'a Pwm, seq: &'a Sequence, min_len: usize) -> Vec<Pwm> {
    let mut pwms = Vec::new();
    let min_len = min_len as i32;

    // if seq_1.len() > seq_2.len() {
    //     panic!("Hmmm...")
    // }

    let mut pwm_shift = -(pwm.len() as i32) + min_len;

    // assert!(seq_1.len() <= seq_2.len());

    loop {
        if pwm_shift > seq.len() as i32 - min_len {
            break;
        }

        let pwm_from = (-pwm_shift).max(0) as usize; // Cut off the start if seq_1_shift is negative
                                                     // seq_1_to only needs to be smaller if it is at the end of the str2
        let last_point_of_pwm = pwm_shift + pwm.len() as i32;
        let overflow = (seq.len() as i32 - last_point_of_pwm).min(0);
        let pwm_to = (pwm.len() as i32 + overflow) as usize;

        let seq_from = (pwm_shift).max(0) as usize;
        let seq_to = (pwm_shift + pwm_to as i32) as usize;

        let mut pwm_clone = pwm.slice(pwm_from..pwm_to);

        pwm_clone.add_sequence_to_pwm(&seq.slice(seq_from..seq_to));

        pwms.push(pwm_clone);

        pwm_shift += 1;
    }

    pwms
}
fn get_all_shift_pwms<'a>(
    mut seq_1: &'a Sequence,
    mut seq_2: &'a Sequence,
    min_len: usize,
) -> Vec<Pwm> {
    let mut pwms = Vec::new();
    let min_len = min_len as i32;

    if seq_1.len() > seq_2.len() {
        std::mem::swap(&mut seq_1, &mut seq_2);
    }

    let mut seq_1_shift = -(seq_1.len() as i32) + min_len;

    assert!(seq_1.len() <= seq_2.len());

    loop {
        if seq_1_shift > seq_2.len() as i32 - min_len {
            break;
        }

        let seq_1_from = (-seq_1_shift).max(0) as usize; // Cut off the start if seq_1_shift is negative
                                                         // seq_1_to only needs to be smaller if it is at the end of the str2
        let last_point_of_seq_1 = seq_1_shift + seq_1.len() as i32;
        let overflow = (seq_2.len() as i32 - last_point_of_seq_1).min(0);
        let seq_1_to = (seq_1.len() as i32 + overflow) as usize;

        let seq_2_from = (seq_1_shift).max(0) as usize;
        let seq_2_to = (seq_1_shift + seq_1_to as i32) as usize;

        assert!(
            seq_1_from < seq_1_to,
            "seq_1_from: {}, seq_1_to: {}",
            seq_1_from,
            seq_1_to
        );
        assert!(
            seq_2_from < seq_2_to,
            "seq_2_from: {}, seq_2_to: {}",
            seq_2_from,
            seq_2_to
        );

        pwms.push(Pwm::new_pfm(&[
            seq_1.slice(seq_1_from..seq_1_to),
            seq_2.slice(seq_2_from..seq_2_to),
        ]));

        seq_1_shift += 1;
    }

    pwms
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_all_shift_pwms() {
        let seq_1 = Sequence::from("ACGT");
        let seq_2 = Sequence::from("ACGTA");
        let min_len = 2;

        let res = get_all_shift_pwms(&seq_1, &seq_2, min_len);

        assert_eq!(
            res[0],
            Pwm::new_pfm(&[Sequence::from("GT"), Sequence::from("AC")])
        );
        assert_eq!(
            res[1],
            Pwm::new_pfm(&[Sequence::from("CGT"), Sequence::from("ACG")])
        );
        assert_eq!(
            res[2],
            Pwm::new_pfm(&[Sequence::from("ACGT"), Sequence::from("ACGT")])
        );
        assert_eq!(
            res[3],
            Pwm::new_pfm(&[Sequence::from("ACGT"), Sequence::from("CGTA")])
        );
        assert_eq!(
            res[4],
            Pwm::new_pfm(&[Sequence::from("ACG"), Sequence::from("GTA")])
        );
        assert_eq!(
            res[5],
            Pwm::new_pfm(&[Sequence::from("AC"), Sequence::from("TA")])
        );
    }

    #[test]
    fn test_get_all_shift_pwms_with_pwm_longer_seq() {
        let pwm = Pwm::new_pfm(&[Sequence::from("ACGT"), Sequence::from("ACGT")]);
        let seq = Sequence::from("ACGTA");

        let res = get_all_shift_pwms_with_pwm(&pwm, &seq, 2);

        assert_eq!(res.len(), 6);
    }

    #[test]
    fn test_get_all_shift_pwms_with_pwm_shorter_seq() {
        let pwm = Pwm::new_pfm(&[Sequence::from("ACGT"), Sequence::from("ACGT")]);
        let seq = Sequence::from("ACG");

        let res = get_all_shift_pwms_with_pwm(&pwm, &seq, 2);

        assert_eq!(res.len(), 4);
    }
}
