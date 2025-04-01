use std::{
    collections::{BinaryHeap, HashMap},
    io::Write,
};

use itertools::Itertools;

use crate::{
    datastructures::{pfm::Pfm, pwm::Pwm, sequence::Sequence},
    plot,
};

struct RankedPfm(Pfm, f64, usize);

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

pub fn motif_finder(seqs: Vec<Sequence>, plot_logos: bool, hits: usize) {
    let mut priority_queue = BinaryHeap::new();
    let max_priority_queue_size = 1_000_000;
    let shrinked_priority_queue_size = 50;

    let only_take_top_score = 3;
    let pfm_min_length = 6;

    let threshold = 0.5;
    let max_gap = Some(2);

    get_all_shift_pfms(&seqs[0], &seqs[1], pfm_min_length)
        .into_iter()
        .flat_map(|x| extraxt_high_interest_area(x, threshold, pfm_min_length, max_gap))
        .for_each(|pfm| {
            let score = pfm.clone().get_custom_score();
            priority_queue.push(RankedPfm(pfm, score, 2));
        });

    println!("Starting search...");
    let mut top_results = HashMap::new();

    while let Some(RankedPfm(pfm, score, indicies)) = priority_queue.pop() {
        if top_results.contains_key(&pfm.get_consensus_string()) {
            continue;
        }
        println!(
            "Score: {:.2}, Depth: {}, Size: {}, String: {:?}",
            score,
            indicies,
            priority_queue.len(),
            pfm.get_consensus_string()
        );

        std::io::stdout().flush().unwrap();

        if indicies >= seqs.len() {
            top_results
                .entry(pfm.get_consensus_string())
                .or_insert((pfm, score));

            if top_results.len() >= hits {
                println!("Found all hits!");
                break;
            }

            continue;
        }

        let next_seq = &seqs[indicies];

        get_all_shift_pfms_with_pfm(&pfm, next_seq, pfm_min_length)
            .into_iter()
            .flat_map(|x| extraxt_high_interest_area(x, threshold, pfm_min_length, max_gap))
            .map(|pfm| {
                let new_score = pfm.get_custom_score();
                RankedPfm(pfm, new_score, indicies + 1)
            })
            .sorted_by(|a, b| b.1.partial_cmp(&a.1).unwrap())
            .take(only_take_top_score)
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

    for (consensus_string, (pwm, score)) in top_results
        .into_iter()
        .sorted_by(|a, b| b.1 .1.partial_cmp(&a.1 .1).unwrap())
        .take(3)
    {
        println!("Score: {:.2}", score);
        println!("{:?}", consensus_string);
        println!("{:?}", pwm.matrix);

        if plot_logos {
            let pwm = Pwm::pfm_into_ppm(pwm);
            plot::plot_pwm(&format!("{}.png", consensus_string), &pwm, score).unwrap();
        }
    }

    plot::clear_cache();
}

fn get_all_shift_pfms_with_pfm<'a>(pfm: &'a Pfm, seq: &'a Sequence, min_len: usize) -> Vec<Pfm> {
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

        pwm_clone.add_sequence(&seq.slice(seq_from..seq_to));

        pwms.push(pwm_clone);

        pwm_shift += 1;
    }

    pwms
}

fn get_all_shift_pfms<'a>(
    mut seq_1: &'a Sequence,
    mut seq_2: &'a Sequence,
    min_len: usize,
) -> Vec<Pfm> {
    let mut pfms = Vec::new();
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

        pfms.push(Pfm::from_sequences(&[
            seq_1.slice(seq_1_from..seq_1_to),
            seq_2.slice(seq_2_from..seq_2_to),
        ]));

        seq_1_shift += 1;
    }

    pfms
}

fn extraxt_high_interest_area(
    pfm: Pfm,
    threshold: f64,
    min_len: usize,
    max_gap: Option<usize>,
) -> Vec<Pfm> {
    let mut pfm_copy = pfm.clone();
    pfm_copy.additive_smoothing(1);

    // dbg!(&pfm_copy);
    let pwm = Pwm::from(pfm_copy);

    let mut high_interest_areas = Vec::new();
    // Find the high interest areas
    // - Starts giving the first high interest area when a positive score is given, backtracks with
    // max gap to find the start of the high interest area, but also adds the first step, including
    // the chance if ealy pfms has faults

    let mut p = 0;
    let mut start_points = Vec::from([0]);
    let mut gap_streak = 0;
    // dbg!(&pwm);
    loop {
        if p >= pwm.matrix.len() {
            break;
        }

        let max_col_value = *pwm.matrix[p]
            .iter()
            .max_by(|&&a, &b| a.partial_cmp(b).unwrap())
            .unwrap();

        if max_col_value > threshold {
            gap_streak = 0;
        } else {
            gap_streak += 1;
        }

        if gap_streak > max_gap.unwrap_or(0) {
            start_points = Vec::new();
        }

        if max_col_value > threshold && start_points.is_empty() {
            start_points.push(p);

            // Add padding to the start incase there is an error in front
            let saturated_padded_start_point = p.saturating_sub(max_gap.unwrap_or(0));
            if saturated_padded_start_point != p {
                start_points.push(saturated_padded_start_point);
            }
        }

        if max_col_value < threshold && !start_points.is_empty() {
            for &start_point in start_points.iter() {
                // Dont add to small sequences
                if p - start_point < min_len {
                    continue;
                }

                high_interest_areas.push(pfm.slice(start_point..p));

                if let Some(max_gap) = max_gap {
                    let saturated_padded_end_point = (p + max_gap).min(pwm.matrix.len() - 1);

                    if saturated_padded_end_point != p {
                        high_interest_areas
                            .push(pfm.slice(start_point..saturated_padded_end_point));
                    }
                }
            }
        }

        p += 1;
    }

    if !start_points.is_empty() {
        for &start_point in start_points.iter() {
            // Dont add to small sequences
            if p - start_point < min_len {
                continue;
            }

            high_interest_areas.push(pfm.slice(start_point..p));

            if let Some(max_gap) = max_gap {
                let saturated_padded_end_point = (p + max_gap).min(pwm.matrix.len() - 1);

                if saturated_padded_end_point != p {
                    high_interest_areas.push(pfm.slice(start_point..saturated_padded_end_point));
                }
            }
        }
    }

    high_interest_areas
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_all_shift_pfms() {
        let seq_1 = Sequence::from("ACGT");
        let seq_2 = Sequence::from("ACGTA");
        let min_len = 2;

        let res = get_all_shift_pfms(&seq_1, &seq_2, min_len);

        assert_eq!(
            res[0],
            Pfm::from_sequences(&[Sequence::from("GT"), Sequence::from("AC")])
        );
        assert_eq!(
            res[1],
            Pfm::from_sequences(&[Sequence::from("CGT"), Sequence::from("ACG")])
        );
        assert_eq!(
            res[2],
            Pfm::from_sequences(&[Sequence::from("ACGT"), Sequence::from("ACGT")])
        );
        assert_eq!(
            res[3],
            Pfm::from_sequences(&[Sequence::from("ACGT"), Sequence::from("CGTA")])
        );
        assert_eq!(
            res[4],
            Pfm::from_sequences(&[Sequence::from("ACG"), Sequence::from("GTA")])
        );
        assert_eq!(
            res[5],
            Pfm::from_sequences(&[Sequence::from("AC"), Sequence::from("TA")])
        );
    }

    #[test]
    fn test_get_all_shift_pfms_with_pfm_longer_seq() {
        let pfm = Pfm::from_sequences(&[Sequence::from("ACGT"), Sequence::from("ACGT")]);
        let seq = Sequence::from("ACGTA");

        let res = get_all_shift_pfms_with_pfm(&pfm, &seq, 2);

        assert_eq!(res.len(), 6);
    }

    #[test]
    fn test_get_all_shift_pfms_with_pfm_shorter_seq() {
        let pfm = Pfm::from_sequences(&[Sequence::from("ACGT"), Sequence::from("ACGT")]);
        let seq = Sequence::from("ACG");

        let res = get_all_shift_pfms_with_pfm(&pfm, &seq, 2);

        assert_eq!(res.len(), 4);
    }

    #[test]
    fn test_get_all_shift_extract_high_interest_area_and_sort_by_custom_score() {
        let pfms = get_all_shift_pfms(
            &Sequence::from("ACGTACGTAACCGGTTACGTACGT"),
            &Sequence::from("TGCATGCAAACCGGTTTGCATGCA"),
            4,
        );

        let clipped_pfms = pfms
            .into_iter()
            .flat_map(|x| extraxt_high_interest_area(x, 0.5, 4, Some(1)))
            .sorted_by(|a, b| {
                b.get_custom_score()
                    .partial_cmp(&a.get_custom_score())
                    .unwrap_or_else(|| {
                        panic!("{:?} {:?}", a.get_custom_score(), b.get_custom_score())
                    })
            })
            .collect::<Vec<_>>();

        assert_eq!(
            clipped_pfms[0].get_consensus_string(),
            String::from("AACCGGTT")
        );
    }

    #[test]
    fn test_extract_high_interest_area_on_already_high_interest_area() {
        let pfm = Pfm::from_sequences(&[Sequence::from("ACGT"), Sequence::from("ACGT")]);

        let res = extraxt_high_interest_area(pfm, 0.5, 2, Some(1));

        let best = res
            .into_iter()
            .sorted_by(|a, b| {
                b.get_custom_score()
                    .partial_cmp(&a.get_custom_score())
                    .unwrap_or_else(|| {
                        panic!("{:?} {:?}", a.get_custom_score(), b.get_custom_score())
                    })
            })
            .next()
            .unwrap()
            .get_consensus_string();

        assert_eq!(best, "ACGT");
    }
}
