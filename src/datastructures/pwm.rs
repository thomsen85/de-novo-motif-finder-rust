use super::{base_matrix::BaseMatrix, pfm::Pfm};

pub type Pwm = BaseMatrix<f64>;

impl Pwm {
    /// Tecnically, this is a Ppm, but it's easier to just call it a Pwm
    pub fn pfm_into_ppm(pfm: Pfm) -> Pwm {
        let matrix = pfm
            .matrix
            .iter()
            .map(|row| row.map(|count| count as f64 / pfm.sample_size as f64))
            .collect::<Vec<[f64; 4]>>();

        Self {
            matrix,
            sample_size: pfm.sample_size,
        }
    }
}

impl From<Pfm> for Pwm {
    /// Converts a Pfm to a Pwm using log-odds scoring
    fn from(pfm: Pfm) -> Self {
        let mut matrix = vec![[0.; 4]; pfm.matrix.len()];
        let background = 0.25;

        for (i, row) in pfm.matrix.iter().enumerate() {
            for (j, count) in row.iter().enumerate() {
                let p = *count as f64 / pfm.sample_size as f64;
                let adjusted_p = if p == 0.0 { 1e-6 } else { p };
                matrix[i][j] = (adjusted_p / background).log2();
            }
        }

        Self {
            matrix,
            sample_size: pfm.sample_size,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
}
