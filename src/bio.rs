use std::{
    ops::{Index, RangeBounds},
    slice::SliceIndex,
};

use itertools::Itertools;

#[derive(Debug, Clone, PartialEq)]
pub struct Pwm {
    pub matrix: Vec<[f64; 4]>,
    pub total_sequences: usize,
    pub is_pwm: bool,
}

impl Pwm {
    pub const POWF: f64 = 0.8;

    pub fn new_pfm(seqs: &[Sequence]) -> Self {
        assert!(seqs.iter().map(|a| a.len()).all_equal());

        let mut matrix = vec![[0.0; 4]; seqs[0].len()];

        for seq in seqs {
            for (i, base) in seq.bases.iter().enumerate() {
                matrix[i][base.to_index()] += 1.0;
            }
        }

        Self {
            matrix,
            total_sequences: seqs.len(),
            is_pwm: false,
        }
    }

    pub fn get_custom_score(&self) -> f64 {
        // let amount_of_sequences = self.history.len() as f64;
        // let pfm_len = self.len() as f64;
        //
        // self.matrix
        //     .iter()
        //     .map(|row| row.iter().map(|&count| count.powf(Self::POWF)).sum::<f64>())
        //     .sum::<f64>()
        //     / (self.len() as f64 * 0.7)
        //     / (self.history.len() as f64 * 0.7)
        //
        // (self.clone().pfm_to_pwm().get_consensus_score() / (self.len() as f64 * 0.90))
        //     + (self.history.len() / 1000) as f64

        // let amount_of_significant_bases = self
        //     .matrix
        //     .iter()
        //     .filter(|row| row.iter().any(|&p| p > 0.8))
        //     .count() as f64;
        //
        // amount_of_significant_bases + self.history.len() as f64
        //
        // assert!(!self.is_pwm);
        // self.clone().get_consensus_score() / (self.len() as f64)

        // let mut score = 0.0;
        // for row in self.clone().pfm_to_pwm().matrix.iter() {
        //     let max = row
        //         .iter()
        //         .max_by(|a, b| a.partial_cmp(b).unwrap())
        //         .unwrap()
        //         .powf(Self::POWF);
        //     score += max;
        // }
        // score

        assert!(!self.is_pwm);
        // Dirichlet Prior or (Additive smoothing is the same when it is uniform distribution)
        let additive_strength = 2.;
        let threshold = 0.7;

        self.matrix
            .iter()
            .map(|row| {
                let max = row
                    .iter()
                    .map(|f| {
                        assert!(*f == 0. || *f >= 1.0, "F was: {}", f);
                        (f + additive_strength) / (self.total_sequences as f64 + additive_strength)
                    })
                    .max_by(|a, b| a.partial_cmp(b).unwrap())
                    .unwrap();
                if max > threshold {
                    max * (self.total_sequences as f64 * 0.02)
                } else {
                    0.0
                }
            })
            .sum::<f64>()
    }

    pub fn pfm_to_pwm(self) -> Self {
        let mut matrix = self.matrix;

        for row in matrix.iter_mut() {
            let sum = row.iter().sum::<f64>();
            for p in row.iter_mut() {
                *p /= sum;
            }
        }

        Self {
            matrix,
            total_sequences: self.total_sequences,
            is_pwm: true,
        }
    }

    /// Assums ACTG is uniform. Must give pfm
    pub fn kullback_leibler_divergence(&self) -> f64 {
        assert!(!self.is_pwm);
        let uniform = 1.0 / 4.0; // Assuming a uniform distribution of ACTG
        let mut kl_divergence = 0.0;
        let pseudo_count = 0.0001; // To avoid log(0)

        for row in self.matrix.iter() {
            for &p in row.iter() {
                let p = p + pseudo_count;

                kl_divergence += p * (p / uniform).ln();
            }
        }

        kl_divergence
    }

    pub fn add_sequence_to_pwm(&mut self, seq: &Sequence) {
        assert_eq!(self.len(), seq.len());
        self.total_sequences += 1;

        for (i, base) in seq.bases.iter().enumerate() {
            self.matrix[i][base.to_index()] += 1.0;
        }
    }

    pub fn len(&self) -> usize {
        self.matrix.len()
    }

    pub fn slice(
        &self,
        range: impl RangeBounds<usize>
            + SliceIndex<[[f64; 4]], Output = [[f64; 4]]>
            + SliceIndex<[Base], Output = [Base]>
            + Clone,
    ) -> Self {
        let matrix = self.matrix[range.clone()].to_vec();

        Self {
            matrix,
            is_pwm: self.is_pwm,
            total_sequences: self.total_sequences,
        }
    }

    pub fn get_consensus_string(&self) -> String {
        self.matrix
            .iter()
            .map(|row| {
                let max = row
                    .iter()
                    .position_max_by(|a, b| a.partial_cmp(b).unwrap())
                    .unwrap();

                match max {
                    0 => 'A',
                    1 => 'C',
                    2 => 'G',
                    3 => 'T',
                    _ => unreachable!(),
                }
            })
            .collect()
    }

    pub fn get_consensus_score(&self) -> f64 {
        assert!(!self.is_pwm);

        let mut score = 0.0;
        for row in self.matrix.iter() {
            let max = row.iter().max_by(|a, b| a.partial_cmp(b).unwrap()).unwrap();
            score += max;
        }
        score
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Base {
    A,
    C,
    G,
    T,
}

impl TryFrom<char> for Base {
    type Error = String;
    fn try_from(c: char) -> Result<Self, Self::Error> {
        match c {
            'A' => Ok(Base::A),
            'C' => Ok(Base::C),
            'G' => Ok(Base::G),
            'T' => Ok(Base::T),
            _ => Err(format!("Invalid base: {}", c)),
        }
    }
}

impl Base {
    pub fn to_index(self) -> usize {
        match self {
            Base::A => 0,
            Base::C => 1,
            Base::G => 2,
            Base::T => 3,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Sequence {
    pub bases: Vec<Base>,
}

impl Sequence {
    pub fn len(&self) -> usize {
        self.bases.len()
    }

    pub fn slice(
        &self,
        range: impl RangeBounds<usize> + SliceIndex<[Base], Output = [Base]>,
    ) -> Self {
        Self {
            bases: self.bases.as_slice()[range].to_vec(),
        }
    }
}

impl FromIterator<Base> for Sequence {
    fn from_iter<T: IntoIterator<Item = Base>>(iter: T) -> Self {
        Sequence {
            bases: iter.into_iter().collect(),
        }
    }
}

impl From<&str> for Sequence {
    fn from(s: &str) -> Self {
        s.chars().filter_map(|c| Base::try_from(c).ok()).collect()
    }
}

impl Index<usize> for Sequence {
    type Output = Base;
    fn index(&self, index: usize) -> &Self::Output {
        &self.bases[index]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pfm_to_pwm() {
        let seqs = vec![
            Sequence::from("ACGT"),
            Sequence::from("ACGT"),
            Sequence::from("ACGT"),
            Sequence::from("AAAA"),
        ];

        let pfm = Pwm::new_pfm(&seqs);
        let pwm = pfm.pfm_to_pwm();

        assert_eq!(pwm.matrix[0][0], 1.);
        assert_eq!(pwm.matrix[1][0], 0.25);
        assert_eq!(pwm.matrix[2][0], 0.25);
        assert_eq!(pwm.matrix[3][0], 0.25);
        assert_eq!(pwm.matrix[1][1], 0.75);
    }

    #[test]
    fn test_consensus_score() {
        let seqs = vec![
            Sequence::from("ACGT"),
            Sequence::from("ACGT"),
            Sequence::from("ACGT"),
            Sequence::from("AAAA"),
        ];

        let pfm = Pwm::new_pfm(&seqs);
        let score = pfm.get_consensus_score();

        assert_eq!(score, 4. + 3. + 3. + 3.);
    }
}
