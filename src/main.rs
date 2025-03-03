use std::{collections::BinaryHeap, path::Path};

use bio::{Pwm, Sequence};

mod bio;
mod fasta_reader;

struct RankedPwm(Pwm, f64);

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
    let path = Path::new("assets/output.fa");
    let seqs = fasta_reader::read_fasta(path);

    let mut priority_queue = BinaryHeap::new();

    get_all_shift_pwms(&seqs[0], &seqs[1], 5)
        .into_iter()
        .for_each(|pwm| {
            let score = pwm.get_custom_score();
            priority_queue.push(RankedPwm(pwm, score));
        });

    for RankedPwm(pwm, score) in priority_queue.into_iter().take(3) {
        println!("Score: {}", score);
        println!("{:?}", pwm);
    }
}

fn get_all_shift_pwms<'a>(
    mut seq_1: &'a Sequence,
    mut seq_2: &'a Sequence,
    min_len: usize,
) -> Vec<Pwm> {
    let mut pwms = Vec::new();
    let min_len = min_len as i32;

    let mut seq_1_shift = -(seq_1.len() as i32) + min_len;

    if seq_1.len() > seq_2.len() {
        std::mem::swap(&mut seq_1, &mut seq_2);
    }

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
}
