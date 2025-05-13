use crate::dfao::ConstantTerm;
use crate::laurent_poly::LaurentPoly;
use crate::mod_int::ModInt;

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct ModIntVector {
    pub entries: Vec<ModInt>,
    pub dim: usize,
    pub modulus: u64,
    pub is_row: bool,
}

impl ModIntVector {
    pub fn new_row(entries: Vec<ModInt>) -> Self {
        assert_eq!(entries.len() % 2, 1);
        let dim = entries.len();
        let modulus = entries.first().unwrap().modulus;
        entries
            .iter()
            .for_each(|entry| assert_eq!(entry.modulus, modulus));
        ModIntVector {
            entries,
            dim,
            modulus,
            is_row: true,
        }
    }

    pub fn new_col(entries: Vec<ModInt>) -> Self {
        assert_eq!(entries.len() % 2, 1);
        let dim = entries.len();
        let modulus = entries.first().unwrap().modulus;
        entries
            .iter()
            .for_each(|entry| assert_eq!(entry.modulus, modulus));
        ModIntVector {
            entries,
            dim,
            modulus,
            is_row: false,
        }
    }

    pub fn from_poly(poly: &LaurentPoly, dim: usize) -> Self {
        assert_eq!(dim % 2, 1);
        let mut row_vec = vec![ModInt::zero(poly.modulus); dim];
        let max_deg = (dim - 1) / 2;
        for k in 0..dim {
            row_vec[k] = poly.get_coefficient(&(k as i64 - max_deg as i64));
        }

        Self::new_row(row_vec)
    }

    pub fn max_deg(&self) -> usize {
        (self.dim - 1) / 2
    }

    pub fn dot(&self, other: &Self) -> ModInt {
        assert_eq!(self.dim, other.dim);
        assert_eq!(self.modulus, other.modulus);
        let mut sum = ModInt::zero(self.modulus);
        for i in 0..self.dim {
            sum = sum + (self.entries[i] * other.entries[i]);
        }
        sum
    }
}

impl ConstantTerm for ModIntVector {
    fn constant_term(&self) -> ModInt {
        self.entries[self.max_deg()]
    }

    fn modulus(&self) -> u64 {
        self.modulus
    }
}
