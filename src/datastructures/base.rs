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

impl TryFrom<usize> for Base {
    type Error = String;

    fn try_from(i: usize) -> Result<Self, Self::Error> {
        match i {
            0 => Ok(Base::A),
            1 => Ok(Base::C),
            2 => Ok(Base::G),
            3 => Ok(Base::T),
            _ => Err(format!("Invalid base index: {}", i)),
        }
    }
}

impl From<Base> for char {
    fn from(base: Base) -> char {
        match base {
            Base::A => 'A',
            Base::C => 'C',
            Base::G => 'G',
            Base::T => 'T',
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
