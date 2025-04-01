// use std::{
//     hash::Hash,
//     ops::{Index, RangeBounds},
//     slice::SliceIndex,
// };
//
// use itertools::Itertools;
//
// #[derive(Debug, Clone, PartialEq)]
// pub struct Pwm {
//     pub matrix: Vec<[f64; 4]>,
//     pub total_sequences: usize,
//     pub is_pwm: bool,
// }
//
// impl Pwm {
//     pub const POWF: f64 = 0.8;
//
//     pub fn get_custom_score(&self) -> f64 {
//         // let amount_of_sequences = self.history.len() as f64;
//         // let pfm_len = self.len() as f64;
//         //
//         // self.matrix
//         //     .iter()
//         //     .map(|row| row.iter().map(|&count| count.powf(Self::POWF)).sum::<f64>())
//         //     .sum::<f64>()
//         //     / (self.len() as f64 * 0.7)
//         //     / (self.history.len() as f64 * 0.7)
//         //
//         // (self.clone().pfm_to_pwm().get_consensus_score() / (self.len() as f64 * 0.90))
//         //     + (self.history.len() / 1000) as f64
//
//         // let amount_of_significant_bases = self
//         //     .matrix
//         //     .iter()
//         //     .filter(|row| row.iter().any(|&p| p > 0.8))
//         //     .count() as f64;
//         //
//         // amount_of_significant_bases + self.history.len() as f64
//         //
//         // assert!(!self.is_pwm);
//         // self.clone().get_consensus_score() / (self.len() as f64)
//
//         // let mut score = 0.0;
//         // for row in self.clone().pfm_to_pwm().matrix.iter() {
//         //     let max = row
//         //         .iter()
//         //         .max_by(|a, b| a.partial_cmp(b).unwrap())
//         //         .unwrap()
//         //         .powf(Self::POWF);
//         //     score += max;
//         // }
//         // score
//
//         assert!(!self.is_pwm);
//         // Dirichlet Prior or (Additive smoothing is the same when it is uniform distribution)
//         let additive_strength = 2.;
//         let threshold = 0.7;
//
//         self.matrix
//             .iter()
//             .map(|row| {
//                 let max = row
//                     .iter()
//                     .map(|f| {
//                         assert!(*f == 0. || *f >= 1.0, "F was: {}", f);
//                         (f + additive_strength) / (self.total_sequences as f64 + additive_strength)
//                     })
//                     .max_by(|a, b| a.partial_cmp(b).unwrap())
//                     .unwrap();
//                 if max > threshold {
//                     max * (self.total_sequences as f64 * 0.02)
//                 } else {
//                     0.0
//                 }
//             })
//             .sum::<f64>()
//     }
//
//     pub fn get_consensus_string(&self) -> String {
//         self.matrix
//             .iter()
//             .map(|row| {
//                 let max = row
//                     .iter()
//                     .position_max_by(|a, b| a.partial_cmp(b).unwrap())
//                     .unwrap();
//
//                 match max {
//                     0 => 'A',
//                     1 => 'C',
//                     2 => 'G',
//                     3 => 'T',
//                     _ => unreachable!(),
//                 }
//             })
//             .collect()
//     }
//
//     pub fn get_consensus_score(&self) -> f64 {
//         assert!(!self.is_pwm);
//
//         let mut score = 0.0;
//         for row in self.matrix.iter() {
//             let max = row.iter().max_by(|a, b| a.partial_cmp(b).unwrap()).unwrap();
//             score += max;
//         }
//         score
//     }
// }
//
// // #[cfg(test)]
// // mod tests {
// //     use super::*;
// //
// //     #[test]
// //     fn test_pfm_to_pwm() {
// //         let seqs = vec![
// //             Sequence::from("ACGT"),
// //             Sequence::from("ACGT"),
// //             Sequence::from("ACGT"),
// //             Sequence::from("AAAA"),
// //         ];
// //
// //         let pfm = Pwm::new_pfm(&seqs);
// //         let pwm = pfm.pfm_to_pwm();
// //
// //         assert_eq!(pwm.matrix[0][0], 1.);
// //         assert_eq!(pwm.matrix[1][0], 0.25);
// //         assert_eq!(pwm.matrix[2][0], 0.25);
// //         assert_eq!(pwm.matrix[3][0], 0.25);
// //         assert_eq!(pwm.matrix[1][1], 0.75);
// //     }
// //
// //     #[test]
// //     fn test_consensus_score() {
// //         let seqs = vec![
// //             Sequence::from("ACGT"),
// //             Sequence::from("ACGT"),
// //             Sequence::from("ACGT"),
// //             Sequence::from("AAAA"),
// //         ];
// //
// //         let pfm = Pwm::new_pfm(&seqs);
// //         let score = pfm.get_consensus_score();
// //
// //         assert_eq!(score, 4. + 3. + 3. + 3.);
// //     }
// // }
