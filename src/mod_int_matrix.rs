use crate::mod_int::ModInt;
use std::fmt;
use std::fmt::Formatter;

pub struct ModIntMatrix {
    entries: Vec<Vec<ModInt>>,
    pub dim: usize,
    pub modulus: u64,
}

impl ModIntMatrix {
    pub fn zero(dim: usize, modulus: u64) -> Self {
        ModIntMatrix {
            entries: vec![vec![ModInt::zero(modulus); dim]; dim],
            dim,
            modulus,
        }
    }

    pub fn new(entries: Vec<Vec<ModInt>>, dim: usize, modulus: u64) -> Self {
        ModIntMatrix {
            entries,
            dim,
            modulus,
        }
    }

    pub fn add(&self, other: &Self) -> Self {
        assert_eq!(self.dim, other.dim);
        assert_eq!(self.modulus, other.modulus);
        let dim = self.dim;
        let modulus = self.modulus;

        let mut entries = vec![vec![ModInt::zero(modulus); dim]; dim];
        for i in 0..dim {
            for j in 0..dim {
                entries[i][j] = self.entries[i][j] + other.entries[i][j];
            }
        }

        ModIntMatrix {
            entries,
            dim: self.dim,
            modulus: self.modulus,
        }
    }

    pub fn mul(&self, other: &Self) -> Self {
        assert_eq!(self.dim, other.dim);
        assert_eq!(self.modulus, other.modulus);
        let dim = self.dim;
        let modulus = self.modulus;
        let mut entries = vec![vec![ModInt::zero(modulus); dim]; dim];
        for i in 0..dim {
            for j in 0..dim {
                for k in 0..dim {
                    entries[i][j] = entries[i][j] + (self.entries[i][k] * other.entries[k][j]);
                }
            }
        }

        ModIntMatrix {
            entries,
            dim: self.dim,
            modulus: self.modulus,
        }
    }
}

impl fmt::Display for ModIntMatrix {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        let mut str = String::new();
        for i in 0..self.dim {
            for j in 0..self.dim {
                str.push_str(&format!("{} ", self.entries[i][j]));
            }
            str.push('\n');
        }
        
        write!(f, "{}", str)
    }
}
