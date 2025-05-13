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
