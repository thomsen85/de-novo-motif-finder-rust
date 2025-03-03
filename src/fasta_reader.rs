use std::{
    fs::File,
    io::{BufReader, Read},
    path::Path,
};

use itertools::Itertools;

use crate::bio::{Base, Sequence};

pub fn read_fasta(path: &Path) -> Vec<Sequence> {
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
