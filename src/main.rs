use std::path::Path;

use bio::{Pwm, Sequence};

mod bio;
mod fasta_reader;

fn main() {
    let path = Path::new("assets/output.fa");
    let sequences = fasta_reader::read_fasta(path);
}

fn get_pwms(sequences: &[Sequence], min_len: usize) -> Vec<Pwm> {
    let mut pwms = Vec::new();

    todo!();

    pwms
}

fn get_all_shift_pwms(seq1: &Sequence, seq2: &Sequence, min_len: usize) -> Vec<Pwm> {
    let mut pwms = Vec::new();

    
    let seq_1_shift = seq1.len() - min_len;
    for

    pwms
}
