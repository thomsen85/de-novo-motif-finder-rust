#![allow(dead_code)]

use std::{ops::RangeBounds, slice::SliceIndex};

use itertools::Itertools;

use super::{base::Base, sequence::Sequence};

#[derive(Debug, Clone, PartialEq)]
pub struct BaseMatrix<T> {
    pub matrix: Vec<[T; 4]>,

    /// The number of sequences in used to create the matrix.
    pub sample_size: usize,
}

impl<T> BaseMatrix<T>
where
    T: Clone + PartialOrd,
{
    pub fn slice(
        &self,
        range: impl RangeBounds<usize> + SliceIndex<[[T; 4]], Output = [[T; 4]]> + Clone,
    ) -> Self {
        let matrix = self.matrix[range].to_vec();

        Self {
            matrix,
            sample_size: self.sample_size,
        }
    }

    pub fn len(&self) -> usize {
        self.matrix.len()
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

    pub fn get_consensus_sequence(&self) -> Sequence {
        self.matrix
            .iter()
            .map(|row| {
                let max = row
                    .iter()
                    .position_max_by(|a, b| a.partial_cmp(b).unwrap())
                    .unwrap();

                Base::try_from(max).unwrap()
            })
            .collect()
    }
}
