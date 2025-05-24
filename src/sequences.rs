use crate::laurent_poly::LaurentPoly;
use crate::mod_int::ModInt;

pub fn constant_term_slow(P: &LaurentPoly, Q: &LaurentPoly, n: &u64) -> ModInt {
    P.pow(n).mul(Q).constant_term()
}

pub fn constant_term_reduce(P: &LaurentPoly, Q: &LaurentPoly, n: &u64) -> LaurentPoly {
    let mut Q = Q.clone();
    let mut n = n.clone();
    let p = P.modulus;

    while n > 0 {
        let r = n % p;
        n = (n - r) / p;
        Q = P.pow(&r).mul(&Q).lambda_reduce()
    }

    Q
}

pub fn constant_term(P: &LaurentPoly, Q: &LaurentPoly, n: &u64) -> ModInt {
    constant_term_reduce(P, Q, n).constant_term()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_constant_term() {
        let primes = vec![2, 3, 5, 7, 11, 13];
        let mots: Vec<u64> = vec![
            1, 1, 2, 4, 9, 21, 51, 127, 323, 835, 2188, 5798, 15511, 41835, 113634, 310572, 853467,
            2356779, 6536382, 18199284, 50852019, 142547559, 400763223, 1129760415, 3192727797,
        ];

        for p in primes {
            let P = LaurentPoly::from_string("x + 1 + x^-1", p);
            let Q = LaurentPoly::from_string("1 - x^2", p);
            for n in 0..mots.len() as u64 {
                let mot_n = constant_term(&P, &Q, &n);
                assert_eq!(mot_n, constant_term_slow(&P, &Q, &n));
                assert_eq!(mot_n, ModInt::new(mots[n as usize], p));
            }
        }
    }
}
