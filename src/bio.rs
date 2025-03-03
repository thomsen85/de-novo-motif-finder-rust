use std::{
    ops::{Index, RangeBounds},
    slice::SliceIndex,
};

use itertools::Itertools;

#[derive(Debug, PartialEq)]
pub struct Pwm {
    matrix: Vec<[f64; 4]>,
}

impl Pwm {
    pub const POWF: f64 = 1.5;

    pub fn new_pfm(seqs: &[Sequence]) -> Self {
        assert!(seqs.iter().map(|a| a.len()).all_equal());

        let mut matrix = vec![[0.0; 4]; seqs[0].len()];

        for seq in seqs {
            for (i, base) in seq.bases.iter().enumerate() {
                matrix[i][base.to_index()] += 1.0;
            }
        }

        Self { matrix }
    }

    pub fn get_custom_score(&self) -> f64 {
        self.matrix
            .iter()
            .map(|row| row.iter().map(|&count| count.powf(Self::POWF)).sum::<f64>())
            .sum::<f64>()
            / self.len() as f64
    }

    pub fn add_sequence_to_pwm(&mut self, seq: &Sequence) {
        assert_eq!(self.len(), seq.len());
        todo!()
    }

    pub fn len(&self) -> usize {
        self.matrix.len()
    }
}

#[derive(Debug, Clone, Copy)]
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

#[derive(Debug, Clone)]
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

// impl<T: RangeBounds<usize>> Index<T> for Sequence {
//     type Output = Sequence;
//     fn index(&self, index: T) -> &Self::Output {
//         &Self {
//             bases: self.bases[index],
//         }
//     }
// }
