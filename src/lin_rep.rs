use crate::laurent_poly::LaurentPoly;
use crate::mod_int::ModInt;
use crate::mod_int_matrix::ModIntMatrix;
use crate::mod_int_vector::ModIntVector;
use std::collections::HashMap;

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

    // BFS the reachable states of the machine whose start state is `start` and whose
    // transition on digit d is step(state, &mat_func[d]). Returns the list of distinct
    // reachable states (deduped by their entry values), or None if the count exceeds bound.
    fn reachable_states<F: Fn(&ModIntVector, &ModIntMatrix) -> ModIntVector>(
        &self,
        start: ModIntVector,
        bound: usize,
        step: F,
    ) -> Option<Vec<ModIntVector>> {
        let p = self.modulus;
        let key = |v: &ModIntVector| -> Vec<u64> { v.entries.iter().map(|e| e.value).collect() };
        let mut states: Vec<ModIntVector> = vec![start.clone()];
        let mut seen: HashMap<Vec<u64>, ()> = HashMap::new();
        seen.insert(key(&start), ());
        let mut k = 0usize;
        while k < states.len() {
            for d in 0..p as usize {
                let ns = step(&states[k], &self.mat_func[d]);
                let kk = key(&ns);
                if !seen.contains_key(&kk) {
                    seen.insert(kk, ());
                    states.push(ns);
                    if states.len() > bound {
                        return None;
                    }
                }
            }
            k += 1;
        }
        Some(states)
    }

    // Basis of the span of `rows` over GF(p) by Gaussian elimination; its length is the rank.
    fn row_basis(rows: &[Vec<u64>], dim: usize, p: u64) -> Vec<Vec<u64>> {
        let mut basis: Vec<Vec<u64>> = Vec::new();
        let mut pivot_col: Vec<usize> = Vec::new();
        for r in rows {
            let mut v = r.clone();
            for (bi, b) in basis.iter().enumerate() {
                let c = pivot_col[bi];
                if v[c] != 0 {
                    let f = v[c]; // b[c] == 1 (normalized)
                    for k in 0..dim {
                        v[k] = (v[k] + (p - f * b[k] % p) % p) % p;
                    }
                }
            }
            if let Some(c) = (0..dim).find(|&k| v[k] != 0) {
                let inv = ModInt::new(v[c], p).inv().value;
                for k in 0..dim {
                    v[k] = v[k] * inv % p;
                }
                basis.push(v);
                pivot_col.push(c);
            }
            if basis.len() == dim {
                break; // full rank; no point scanning further
            }
        }
        basis
    }

    // Rank of the matrix `m` over GF(p) by Gaussian elimination.
    fn mat_rank(mut m: Vec<Vec<u64>>, p: u64) -> usize {
        let rows = m.len();
        if rows == 0 {
            return 0;
        }
        let cols = m[0].len();
        if cols == 0 {
            return 0;
        }
        let mut rank = 0usize;
        let mut row = 0usize;
        for col in 0..cols {
            let sel = (row..rows).find(|&r| m[r][col] != 0);
            if let Some(s) = sel {
                m.swap(row, s);
                let inv = ModInt::new(m[row][col], p).inv().value;
                for k in 0..cols {
                    m[row][k] = m[row][k] * inv % p;
                }
                for r in 0..rows {
                    if r != row && m[r][col] != 0 {
                        let f = m[r][col];
                        for k in 0..cols {
                            m[r][k] = (m[r][k] + (p - f * m[row][k] % p) % p) % p;
                        }
                    }
                }
                row += 1;
                rank += 1;
                if row == rows {
                    break;
                }
            }
        }
        rank
    }

    // The Hankel rank = dimension of a MINIMAL linear representation of the sequence.
    //
    // The `rank` field of this struct is the AMBIENT dimension (2*max_deg + 1); the Hankel
    // rank is generally smaller. It equals rank(RowMat * ColMat), where RowMat's rows are
    // the reachable forward (lsd) row states (start = row_vec, step s -> s * mat) and
    // ColMat's columns are the reachable reverse (msd) column states (start = col_vec,
    // step s -> mat * s). It is <= min(N_lsd, N_msd), the two minimal-automaton state counts.
    //
    // The product is computed via row/column-space bases so the cost is
    // O((R + C)*dim^2 + dim^3) rather than materializing the full R x C product. `bound`
    // caps the reachable-state BFS in each direction; None is returned if it is exceeded.
    pub fn hankel_rank(&self, bound: usize) -> Option<usize> {
        let p = self.modulus;
        let dim = self.row_vec.dim;
        let lsd = self.reachable_states(self.row_vec.clone(), bound, |s, m| m.right_mul(s))?;
        let msd = self.reachable_states(self.col_vec.clone(), bound, |s, m| m.left_mul(s))?;
        let to_u64 = |s: &ModIntVector| -> Vec<u64> { s.entries.iter().map(|e| e.value).collect() };
        let lsd_rows: Vec<Vec<u64>> = lsd.iter().map(to_u64).collect();
        let msd_rows: Vec<Vec<u64>> = msd.iter().map(to_u64).collect();
        let ba = Self::row_basis(&lsd_rows, dim, p); // basis of lsd-state span (rows of RowMat)
        let bb = Self::row_basis(&msd_rows, dim, p); // basis of msd-state span (cols of ColMat)
        let g: Vec<Vec<u64>> = ba
            .iter()
            .map(|a| {
                bb.iter()
                    .map(|b| {
                        let mut acc = 0u64;
                        for k in 0..dim {
                            acc += a[k] * b[k];
                        }
                        acc % p
                    })
                    .collect()
            })
            .collect();
        Some(Self::mat_rank(g, p))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // P = t^-1 + 1 + t, Q = 1 - t^2: the Motzkin functional ct[P^n (1-t^2)].
    fn motzkin(p: u64) -> (LaurentPoly, LaurentPoly) {
        let big = LaurentPoly::from_vec(vec![(-1, 1), (0, 1), (1, 1)], p);
        let q = LaurentPoly::from_vec(vec![(0, 1), (2, p - 1)], p);
        (big, q)
    }

    fn poly(terms: &[(i64, i64)], p: u64) -> LaurentPoly {
        let v: Vec<(i64, u64)> = terms
            .iter()
            .map(|&(e, c)| (e, ((c % p as i64 + p as i64) % p as i64) as u64))
            .collect();
        LaurentPoly::from_vec(v, p)
    }

    #[test]
    fn motzkin_hankel_rank_is_three() {
        // The census records rank=3 for Motzkin at p=3,5,7 (and the ambient is 5, so this
        // is strictly smaller than 2*max_deg+1 -- it is not re-reporting the ambient dim).
        for p in [3u64, 5, 7] {
            let (big, q) = motzkin(p);
            let rep = LinRep::for_ct_sequence(&big, &q);
            assert_eq!(rep.rank, 5, "ambient dim should be 5 at p={p}");
            assert_eq!(
                rep.hankel_rank(100000),
                Some(3),
                "Motzkin Hankel rank should be 3 at p={p}"
            );
            assert!(
                rep.hankel_rank(100000).unwrap() < rep.rank,
                "Hankel rank must be strictly below ambient at p={p}"
            );
        }
    }

    #[test]
    fn asymmetric_anchors_match_census() {
        // Anchors lifted from experiments/lsd-msd-state-census smoke (main.rs).
        // (P-terms, Q-terms, p, expected Hankel rank)
        let cases: &[(&[(i64, i64)], &[(i64, i64)], u64, usize)] = &[
            (&[(-1, 1), (2, 1)], &[(0, 1), (2, -1)], 11, 3),
            (&[(-1, 1), (2, 1)], &[(0, 1), (2, -1)], 13, 3),
            (&[(-2, 1), (0, 1), (1, 1)], &[(0, 1), (2, -1)], 13, 4),
            (&[(-1, 1), (0, 1), (1, 1), (2, 1)], &[(1, 1)], 3, 2),
        ];
        for (pt, qt, p, expected) in cases {
            let big = poly(pt, *p);
            let q = poly(qt, *p);
            let rep = LinRep::for_ct_sequence(&big, &q);
            let r = rep.hankel_rank(1_000_000).unwrap();
            assert_eq!(r, *expected, "Hankel rank mismatch at p={p}");
            assert!(r <= rep.rank, "Hankel rank exceeds ambient at p={p}");
        }
    }
}
