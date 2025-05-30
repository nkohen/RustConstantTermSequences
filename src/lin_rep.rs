use crate::laurent_poly::LaurentPoly;
use crate::mod_int::ModInt;
use crate::mod_int_matrix::ModIntMatrix;
use crate::mod_int_vector::ModIntVector;

pub struct LinRep {
    pub row_vec: ModIntVector,
    pub mat_func: Vec<ModIntMatrix>,
    pub col_vec: ModIntVector,
    pub rank: usize,
    pub modulus: u64,
}

impl LinRep {
    fn compute_mat_for_poly(poly: &LaurentPoly, max_deg: usize) -> ModIntMatrix {
        let dim = 2 * max_deg + 1;
        let modulus = poly.modulus;
        let mut entries = vec![vec![ModInt::zero(modulus); dim]; dim];
        let deg = max_deg as i64;

        for i in 0..dim {
            for j in 0..dim {
                let index = (deg - i as i64) - ((deg - j as i64) * modulus as i64);
                entries[i][j] = poly.get_coefficient(&index);
            }
        }

        ModIntMatrix::new(entries, 2 * max_deg + 1, modulus)
    }

    pub fn for_ct_sequence(P: &LaurentPoly, Q: &LaurentPoly) -> Self {
        assert_eq!(P.modulus, Q.modulus);
        let modulus = P.modulus;
        let max_deg = std::cmp::max(P.degree() - 1, Q.degree()) as usize;
        let mut poly = LaurentPoly::one(modulus);
        let mut mats: Vec<ModIntMatrix> = vec![];
        for _ in 0..modulus {
            mats.push(Self::compute_mat_for_poly(&poly, max_deg));
            poly = poly.mul(P);
        }

        let mut col_vec = vec![ModInt::zero(modulus); 2 * max_deg + 1];
        col_vec[max_deg] = ModInt::new(1, modulus);

        LinRep {
            row_vec: ModIntVector::from_poly(Q, 2 * max_deg + 1),
            mat_func: mats,
            col_vec: ModIntVector::new_col(col_vec),
            rank: 2 * max_deg + 1,
            modulus,
        }
    }

    pub fn compute_functional(&self, n: u64) -> ModIntVector {
        let p = self.modulus;
        let mut n = n;
        let mut digits: Vec<ModInt> = Vec::new();
        while n > 0 {
            let r = n % p;
            digits.push(ModInt::new(r, p));
            n = (n - r) / p;
        }
        digits.reverse();

        let mut functional: ModIntVector = self.col_vec.clone();
        for digit in digits {
            functional = self.mat_func[digit.value as usize].left_mul(&functional);
        }

        functional
    }

    pub fn compute(&self, n: u64) -> ModInt {
        let functional = self.compute_functional(n);
        self.row_vec.dot(&functional)
    }
}
