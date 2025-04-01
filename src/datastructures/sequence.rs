use std::{
    fmt::Debug,
    ops::{Index, RangeBounds},
    slice::SliceIndex,
};

use super::base::Base;

#[derive(Clone, PartialEq)]
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

impl Debug for Sequence {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(
            f,
            "{}",
            self.bases
                .iter()
                .map(|&b| char::from(b))
                .collect::<String>()
        )
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

impl From<Vec<Base>> for Sequence {
    fn from(bases: Vec<Base>) -> Self {
        Sequence { bases }
    }
}

impl Index<usize> for Sequence {
    type Output = Base;
    fn index(&self, index: usize) -> &Self::Output {
        &self.bases[index]
    }
}
