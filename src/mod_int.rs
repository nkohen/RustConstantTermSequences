use std::fmt;

#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
pub struct ModInt {
    pub value: u64,
    pub modulus: u64,
}

impl ModInt {
    pub fn new(value: u64, modulus: u64) -> ModInt {
        ModInt {
            value: value % modulus,
            modulus,
        }
    }

    pub fn from_i64(value: i64, modulus: u64) -> Self {
        let mut val = value % (modulus as i64);
        if val < 0 {
            val += modulus as i64;
        }
        ModInt {
            value: val as u64,
            modulus,
        }
    }

    pub fn zero(modulus: u64) -> Self {
        ModInt::new(0, modulus)
    }

    pub fn pow(self, exp: u64) -> Self {
        let mut base = self.value;
        let mut exp = exp;
        let mut result = 1;
        let m = self.modulus;
        while exp > 0 {
            if exp % 2 == 1 {
                result = result * base % m;
            }
            base = base * base % m;
            exp /= 2;
        }
        ModInt::new(result, self.modulus)
    }

    pub fn inv(self) -> Self {
        // Fermat's little theorem (modulus must be prime)
        self.pow(self.modulus - 2)
    }

    pub fn get_digits(n: u64, base: u64) -> Vec<ModInt> {
        let mut n = n;
        let mut digits: Vec<ModInt> = Vec::new();
        while n > 0 {
            let r = n % base;
            digits.push(ModInt::new(r, base));
            n = (n - r) / base;
        }
        digits
    }
}

use std::ops::{Add, Div, Mul, Sub};

impl Add for ModInt {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        assert_eq!(self.modulus, rhs.modulus);
        ModInt::new((self.value + rhs.value) % self.modulus, self.modulus)
    }
}

impl Sub for ModInt {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        assert_eq!(self.modulus, rhs.modulus);
        ModInt::from_i64(self.value as i64 - rhs.value as i64, self.modulus)
    }
}

impl Mul for ModInt {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        assert_eq!(self.modulus, rhs.modulus);
        ModInt::new((self.value * rhs.value) % self.modulus, self.modulus)
    }
}

impl Div for ModInt {
    type Output = Self;
    fn div(self, rhs: Self) -> Self {
        self * rhs.inv()
    }
}

impl fmt::Display for ModInt {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.value)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new() {
        assert_eq!(ModInt::new(42, 10).value, 2);
        assert_eq!(ModInt::new(11, 5), ModInt::new(1, 5));
        assert_ne!(ModInt::new(0, 2), ModInt::new(0, 3));
    }

    #[test]
    fn test_from_i64() {
        assert_eq!(ModInt::from_i64(-23, 10), ModInt::from_i64(7, 10));
        assert_eq!(ModInt::from_i64(11, 5), ModInt::from_i64(1, 5));
        assert_ne!(ModInt::from_i64(0, 2), ModInt::from_i64(0, 3));
    }

    #[test]
    fn test_zero() {
        assert_eq!(ModInt::zero(42).value, 0);
        assert_eq!(ModInt::zero(11), ModInt::zero(11));
        assert_ne!(ModInt::zero(11), ModInt::zero(12));
    }

    #[test]
    fn test_pow() {
        for n in 1..10 {
            for k in 2..5 {
                for i in 0..10 {
                    assert_eq!(ModInt::new(k, n).pow(i as u64), ModInt::new(k.pow(i), n))
                }
            }
        }

        assert_eq!(ModInt::new(23, 42).pow(2300), ModInt::new(25, 42));
    }

    #[test]
    fn test_inv() {
        let primes = vec![2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
        for p in primes {
            for n in 1..p {
                let num = ModInt::new(n, p);
                assert_eq!(num * num.inv(), ModInt::new(1, p));
            }
        }

        assert_eq!(ModInt::new(87, 101).inv(), ModInt::new(36, 101));
    }

    #[test]
    fn test_get_digits() {
        let digits: Vec<ModInt> = vec![4, 7, 1, 1, 2, 2]
            .iter()
            .map(|n: &u64| ModInt::new(n.clone(), 8))
            .collect();
        assert_eq!(ModInt::get_digits(74364, 8), digits);

        let digits: Vec<ModInt> = vec![2, 2, 4, 0, 4, 0, 4, 0, 2]
            .iter()
            .map(|n: &u64| ModInt::new(n.clone(), 5))
            .collect();
        assert_eq!(ModInt::get_digits(846362, 5), digits);
    }

    #[test]
    fn test_add() {
        for i in 0..100 {
            for j in 0..100 {
                for m in 1..100 {
                    assert_eq!(ModInt::new(i, m) + ModInt::new(j, m), ModInt::new(i + j, m));
                }
            }
        }
    }

    #[test]
    #[should_panic]
    fn test_add_panic() {
        let _ = ModInt::new(1, 5) + ModInt::new(2, 7);
    }

    #[test]
    fn test_mul() {
        for i in 0..100 {
            for j in 0..100 {
                for m in 1..100 {
                    assert_eq!(ModInt::new(i, m) * ModInt::new(j, m), ModInt::new(i * j, m));
                }
            }
        }
    }

    #[test]
    #[should_panic]
    fn test_mul_panic() {
        let _ = ModInt::new(1, 5) * ModInt::new(2, 7);
    }

    #[test]
    fn test_div() {
        let primes = vec![2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];

        for p in primes {
            for i in 0..100 {
                for j in 1..p {
                    assert_eq!(
                        (ModInt::new(i, p) / ModInt::new(j, p)) * ModInt::new(j, p),
                        ModInt::new(i, p)
                    );
                }
            }
        }
    }

    #[test]
    #[should_panic]
    fn test_div_panic() {
        let _ = ModInt::new(1, 5) / ModInt::new(2, 7);
    }
}
