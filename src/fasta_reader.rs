use std::path::Path;

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

pub fn read_fasta(path: &Path) -> Vec<Vec<Base>> {
    todo!();
}
