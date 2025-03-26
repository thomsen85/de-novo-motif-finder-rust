use std::{collections::BinaryHeap, io::Write, path::Path};

use bio::{Pwm, Sequence};
use clap::Parser;
use itertools::Itertools;
use plot::clear_cache;

mod bio;
mod fasta_reader;
mod plot;

struct RankedPfm(Pwm, f64, usize);

impl PartialEq for RankedPfm {
    fn eq(&self, other: &Self) -> bool {
        self.1 == other.1
    }
}

impl Eq for RankedPfm {}

impl PartialOrd for RankedPfm {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for RankedPfm {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.1.partial_cmp(&other.1).unwrap()
    }
}

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    /// Name of the person to greet
    #[arg(short, long)]
    input_file: String,
}

fn main() {
    let args = Args::parse();

    let path = Path::new(&args.input_file);

    if !path.exists() {
        panic!("File does not exist: {:?}", path);
    }

    let seqs = fasta_reader::read_fasta(path);

    println!("Input file: {:?}", path);
    println!("Sequences: {}", seqs.len());
    println!(
        "Average sequence length: {:.2}",
        seqs.iter().map(|x| x.len()).sum::<usize>() as f64 / seqs.len() as f64
    );
    print!("====================");

    let mut priority_queue = BinaryHeap::new();
    let max_priority_queue_size = 1000;
    let shrinked_priority_queue_size = 50;

    let only_take_top_score = 3;
    let pfm_min_length = 6;

    get_all_shift_pfms(&seqs[0], &seqs[1], pfm_min_length)
        .into_iter()
        .for_each(|pfm| {
            let score = pfm.clone().get_custom_score();
            priority_queue.push(RankedPfm(pfm, score, 2));
        });

    println!("Starting search...");
    let mut top_results = Vec::new();
    // let mut seen = HashSet::new();

    while let Some(RankedPfm(pfm, score, indicies)) = priority_queue.pop() {
        print!(
            "\rScore: {:.2}, Depth: {}, Size: {}",
            score,
            indicies,
            priority_queue.len()
        );

        std::io::stdout().flush().unwrap();

        // let consensus = pfm.get_consensus_string();
        // if seen.contains(&(consensus.clone(), indicies.clone())) {
        //     continue;
        // }
        // seen.insert((consensus, indicies.clone()));

        if indicies >= seqs.len() {
            top_results.push((pfm, score));

            if top_results.len() >= only_take_top_score {
                break;
            }
            println!("Score: {:.2}", score);

            continue;
        }

        let next_seq = &seqs[indicies];

        get_all_shift_pfms_with_pfm(&pfm, next_seq, pfm_min_length)
            .into_iter()
            .map(|pfm| {
                let new_score = pfm.get_custom_score();
                RankedPfm(pfm, new_score, indicies + 1)
            })
            .sorted_by(|a, b| b.1.partial_cmp(&a.1).unwrap())
            .take(3)
            .for_each(|x| {
                priority_queue.push(x);
            });

        if priority_queue.len() > max_priority_queue_size {
            let mut new_priority_queue = BinaryHeap::new();
            for _ in 0..shrinked_priority_queue_size {
                new_priority_queue.push(priority_queue.pop().unwrap());
            }
            priority_queue = new_priority_queue;
        }
    }

    for (pwm, score) in top_results
        .into_iter()
        .sorted_by(|a, b| b.1.partial_cmp(&a.1).unwrap())
        .take(3)
    {
        println!("Score: {:.2}", score);
        println!("{:?}", pwm.get_consensus_string());
        println!("{:?}", pwm.matrix);
        let pwm = pwm.pfm_to_pwm();
        plot::plot_pwm(
            &format!("{}.png", pwm.clone().get_consensus_string()),
            &pwm,
            score,
        )
        .unwrap();
    }

    clear_cache();
}

fn get_all_shift_pfms_with_pfm<'a>(pfm: &'a Pwm, seq: &'a Sequence, min_len: usize) -> Vec<Pwm> {
    assert!(!pfm.is_pwm);

    let mut pwms = Vec::new();
    let min_len = min_len as i32;

    let mut pwm_shift = -(pfm.len() as i32) + min_len;

    loop {
        if pwm_shift > seq.len() as i32 - min_len {
            break;
        }

        let pwm_from = (-pwm_shift).max(0) as usize; // Cut off the start if seq_1_shift is negative
                                                     // seq_1_to only needs to be smaller if it is at the end of the str2
        let last_point_of_pwm = pwm_shift + pfm.len() as i32;
        let overflow = (seq.len() as i32 - last_point_of_pwm).min(0);
        let pwm_to = (pfm.len() as i32 + overflow) as usize;

        let seq_from = (pwm_shift).max(0) as usize;
        let seq_to = (pwm_shift + pwm_to as i32) as usize;

        let mut pwm_clone = pfm.slice(pwm_from..pwm_to);

        pwm_clone.add_sequence_to_pwm(&seq.slice(seq_from..seq_to));

        pwms.push(pwm_clone);

        pwm_shift += 1;
    }

    pwms
}

fn get_all_shift_pfms<'a>(
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

        let res = get_all_shift_pfms(&seq_1, &seq_2, min_len);

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

        let res = get_all_shift_pfms_with_pfm(&pwm, &seq, 2);

        assert_eq!(res.len(), 6);
    }

    #[test]
    fn test_get_all_shift_pwms_with_pwm_shorter_seq() {
        let pwm = Pwm::new_pfm(&[Sequence::from("ACGT"), Sequence::from("ACGT")]);
        let seq = Sequence::from("ACG");

        let res = get_all_shift_pfms_with_pfm(&pwm, &seq, 2);

        assert_eq!(res.len(), 4);
    }
}
