use crate::mod_int::ModInt;
use crate::mod_int_vector::ModIntVector;
use std::fmt;
use std::fmt::Formatter;

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
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

    pub fn id(dim: usize, modulus: u64) -> Self {
        let mut entries = vec![vec![ModInt::zero(modulus); dim]; dim];
        for i in 0..dim {
            entries[i][i] = ModInt::new(1, modulus);
        }

        ModIntMatrix::new(entries, dim, modulus)
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

    pub fn left_mul(&self, col_vec: &ModIntVector) -> ModIntVector {
        assert_eq!(self.dim, col_vec.dim);
        assert_eq!(self.modulus, col_vec.modulus);
        let mut entries = vec![ModInt::zero(self.modulus); self.dim];
        for i in 0..self.dim {
            for j in 0..self.dim {
                entries[i] = entries[i] + (self.entries[i][j] * col_vec.entries[j]);
            }
        }

        ModIntVector::new(entries)
    }

    pub fn right_mul(&self, row_vec: &ModIntVector) -> ModIntVector {
        assert_eq!(self.dim, row_vec.dim);
        assert_eq!(self.modulus, row_vec.modulus);
        let mut entries = vec![ModInt::zero(self.modulus); self.dim];
        for j in 0..self.dim {
            for i in 0..self.dim {
                entries[j] = entries[j] + (row_vec.entries[i] * self.entries[i][j]);
            }
        }

        ModIntVector::new(entries)
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
