use crate::dfao::ConstantTerm;
use crate::mod_int::ModInt;

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct ModIntVector {
    pub entries: Vec<ModInt>,
    pub dim: usize,
    pub modulus: u64,
}

impl ModIntVector {
    pub fn new(entries: Vec<ModInt>) -> Self {
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
        }
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
