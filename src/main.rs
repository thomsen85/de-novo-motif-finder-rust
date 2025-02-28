use std::path::Path;

use fasta_reader::Base;
use itertools::Itertools;

mod fasta_reader;

fn main() {
    let path = Path::new("assets/output.fa");
    let sequences = fasta_reader::read_fasta(path);
    dbg!(sequences.len());

    // Just limiting for testing
    let sequences = sequences.iter().take(4).collect::<Vec<_>>();

    let pwms = get_all_slided_pwm(sequences[0], sequences[1], 5);

    dbg!(pwms.len());
    let a = pwms
        .iter()
        .map(|pwm| (pwm, get_custom_score(pwm)))
        .sorted_by(|a, b| a.1.partial_cmp(&b.1).unwrap())
        .collect::<Vec<_>>();

    dbg!(a.first().unwrap());
    dbg!(a.last().unwrap());
}

fn get_all_slided_pwm(seq1: &[Base], seq2: &[Base], min_len: usize) -> Vec<Pwm> {
    let seq1_len = seq1.len();
    let seq2_len = seq2.len();

    let mut pwms = Vec::new();

    for i in min_len..seq1_len.min(seq2_len) {
        pwms.push(make_pwm(&seq1[seq1_len - i..], &seq2[..i]));
    }

    pwms
}

type Pwm = Vec<Vec<i32>>;

fn make_pwm(seq1: &[Base], seq2: &[Base]) -> Pwm {
    assert_eq!(seq1.len(), seq2.len());

    let mut pwm = vec![vec![0; 4]; seq1.len()];
    for (i, (base1, base2)) in seq1.iter().zip(seq2.iter()).enumerate() {
        pwm[i][base1.to_index()] += 1;
        pwm[i][base2.to_index()] += 1;
    }
    pwm
}

fn get_consensus_score(pwm: &Pwm) -> i32 {
    pwm.iter().map(|row| row.iter().max().unwrap()).sum()
}

fn get_custom_score(pwm: &Pwm) -> f64 {
    // Award big similaries and penalize big differences
    const POWF: f64 = 1.5;
    let pwm_len = pwm.len() as f64;

    pwm.iter()
        .map(|row| {
            row.iter()
                .map(|&count| (count as f64).powf(POWF))
                .sum::<f64>()
        })
        .sum::<f64>()
        / pwm_len
}

fn transpose2<T>(v: Vec<Vec<T>>) -> Vec<Vec<T>> {
    assert!(!v.is_empty());
    let len = v[0].len();
    let mut iters: Vec<_> = v.into_iter().map(|n| n.into_iter()).collect();
    (0..len)
        .map(|_| {
            iters
                .iter_mut()
                .map(|n| n.next().unwrap())
                .collect::<Vec<T>>()
        })
        .collect()
}
