use crate::laurent_poly::LaurentPoly;
use crate::mod_int::ModInt;

pub fn constant_term_slow(P: &LaurentPoly, Q: &LaurentPoly, n: &u64) -> ModInt {
    P.pow(n).mul(Q).constant_term()
}

pub fn constant_term(P: &LaurentPoly, Q: &LaurentPoly, n: &u64) -> ModInt {
    let mut Q = Q.clone();
    let mut n = n.clone();
    let p = P.modulus;

    while n > 0 {
        let r = n % p;
        n = (n - r) / p;
        Q = P.pow(&r).mul(&Q).lambda_reduce()
    }

    Q.constant_term()
}
