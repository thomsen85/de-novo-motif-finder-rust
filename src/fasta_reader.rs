use std::{
    fs::File,
    io::{BufReader, Read},
    path::Path,
};

use itertools::Itertools;

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
    pub fn to_index(&self) -> usize {
        match self {
            Base::A => 0,
            Base::C => 1,
            Base::G => 2,
            Base::T => 3,
        }
    }
}

pub fn read_fasta(path: &Path) -> Vec<Vec<Base>> {
    let file = File::open(path)
        .unwrap_or_else(|_| panic!("Could not opne file at path: {}", path.to_str().unwrap()));

    let mut reader = BufReader::new(file);
    let mut file_contens = String::new();
    reader.read_to_string(&mut file_contens).unwrap();

    file_contens
        .split("\n")
        .chunks(2)
        .into_iter()
        .map(|chunk| {
            chunk
                .last()
                .unwrap()
                .chars()
                .filter_map(|c| Base::try_from(c.to_ascii_uppercase()).ok())
                .collect()
        })
        .collect()
}
