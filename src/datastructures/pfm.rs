#![allow(dead_code)]

use itertools::Itertools;

use super::{base_matrix::BaseMatrix, sequence::Sequence};

pub type Pfm = BaseMatrix<usize>;

impl Pfm {
    pub fn from_sequences(seqs: &[Sequence]) -> Self {
        assert!(seqs.iter().map(|a| a.len()).all_equal());

        let mut matrix = vec![[0; 4]; seqs[0].len()];

        for seq in seqs {
            for (i, base) in seq.bases.iter().enumerate() {
                matrix[i][base.to_index()] += 1;
            }
        }

        Self {
            matrix,
            sample_size: seqs.len(),
        }
    }

    pub fn add_sequence(&mut self, seq: &Sequence) {
        assert_eq!(self.len(), seq.len());

        self.sample_size += 1;

        for (i, base) in seq.bases.iter().enumerate() {
            self.matrix[i][base.to_index()] += 1;
        }
    }

    pub fn additive_smoothing(&mut self, pseudocount: usize) {
        for row in self.matrix.iter_mut() {
            for count in row.iter_mut() {
                *count += pseudocount;
            }
        }

        self.sample_size += pseudocount * 4;
    }

    /// Assums ACTG is uniform. Must give pfm
    pub fn kullback_leibler_divergence(&self) -> f64 {
        let uniform = 1.0 / 4.0; // Assuming a uniform distribution of ACTG
        let mut kl_divergence = 0.0;
        let pseudo_count = 0.0001; // To avoid log(0)

        for row in self.matrix.iter() {
            for &p in row.iter() {
                let p = p as f64 + pseudo_count;

                kl_divergence += p * (p / uniform).ln();
            }
        }

        kl_divergence
    }

    //https://www.maths.usyd.edu.au/u/uri/my_papers/2006_Evalue_finders_RecombRG_draft.pdf
    pub fn get_custom_score(&self) -> f64 {
        let smoothing = 0.001;
        self.matrix
            .iter()
            .map(|row| {
                row.iter()
                    .map(|&val| {
                        (val as f64 + smoothing)
                            * (((val as f64 + smoothing) / self.sample_size as f64) / 0.25).log2()
                    })
                    .sum::<f64>()
            })
            .sum()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_custom_score_1() {
        let pfm_best = Pfm::from_sequences(&["ACGT".into(), "ACGT".into()]);
        let pfm_not_best = Pfm::from_sequences(&["AACGT".into(), "CACGT".into()]);

        assert!(pfm_best.get_custom_score() > pfm_not_best.get_custom_score());
    }

    #[test]
    fn test_get_custom_score_2() {
        let pfm_best = Pfm::from_sequences(&["CCGTTT".into(), "ACGTTT".into()]);
        let pfm_not_best = Pfm::from_sequences(&["CCGT".into(), "ACGT".into()]);

        assert!(
            pfm_best.get_custom_score() > pfm_not_best.get_custom_score(),
            "{:?} > {:?}",
            pfm_best.get_custom_score(),
            pfm_not_best.get_custom_score()
        );
    }
}
