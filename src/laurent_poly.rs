use crate::mod_int::ModInt;
use std::collections::BTreeMap;
use std::fmt;

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct LaurentPoly {
    terms: BTreeMap<i64, ModInt>, // exponent -> coefficient
    pub modulus: u64,
}

impl LaurentPoly {
    pub fn new(terms: BTreeMap<i64, ModInt>, modulus: u64) -> LaurentPoly {
        LaurentPoly { terms, modulus }
    }

    pub fn zero(modulus: u64) -> Self {
        LaurentPoly::new(BTreeMap::new(), modulus)
    }

    pub fn one(modulus: u64) -> Self {
        let mut poly = LaurentPoly::zero(modulus);

        poly.terms.insert(0, ModInt::new(1, modulus));

        poly
    }

    pub fn from_vec(vec: Vec<(i64, u64)>, modulus: u64) -> Self {
        let mut poly = LaurentPoly::zero(modulus);
        for (exp, coeff) in vec {
            let m = ModInt::new(coeff, modulus);
            if m.value != 0 {
                poly.terms.insert(exp, m);
            }
        }
        poly
    }

    pub fn from_string(s: &str, modulus: u64) -> Self {
        let mut terms = BTreeMap::new();
        let mut tokens = Vec::new();
        let mut i = 0;
        let chars: Vec<char> = s.chars().collect();
        let mut current = String::new();

        while i < chars.len() {
            if (chars[i] == '+' || chars[i] == '-') && i > 0 && chars[i - 1] != '^' {
                if !current.trim().is_empty() {
                    tokens.push(current.trim().to_string());
                }
                current = chars[i].to_string();
            } else {
                current.push(chars[i]);
            }
            i += 1;
        }
        if !current.trim().is_empty() {
            tokens.push(current.trim().to_string());
        }

        for token in tokens {
            let (coeff, exp) = if let Some(idx) = token.find('x') {
                let coeff_str = token[..idx].trim();
                let coeff = match coeff_str {
                    "" | "+" => 1,
                    "-" => -1,
                    _ => coeff_str
                        .trim_start_matches('+')
                        .chars()
                        .filter(|c| !c.is_whitespace())
                        .collect::<String>()
                        .parse::<i64>()
                        .unwrap(),
                };

                let exp = if let Some(exp_str) = token.get((idx + 1)..) {
                    if exp_str.starts_with('^') {
                        exp_str[1..].parse::<i64>().unwrap()
                    } else {
                        1
                    }
                } else {
                    1
                };

                (coeff, exp)
            } else {
                // Constant term — must strip leading '+' if present
                let coeff_str = token.trim_start_matches('+').trim();
                let coeff = coeff_str
                    .chars()
                    .filter(|c| !c.is_whitespace())
                    .collect::<String>()
                    .parse::<i64>()
                    .unwrap();
                (coeff, 0)
            };

            let modint = ModInt::from_i64(coeff, modulus);
            if modint.value != 0 {
                terms
                    .entry(exp)
                    .and_modify(|c| *c = *c + modint)
                    .or_insert(modint);
            }
        }

        LaurentPoly::new(terms, modulus)
    }

    pub fn add(&self, other: &Self) -> Self {
        assert_eq!(self.modulus, other.modulus);
        let mut result = self.terms.clone();
        for (&exp, &coeff) in &other.terms {
            result
                .entry(exp)
                .and_modify(|c| *c = *c + coeff)
                .or_insert(coeff);
        }
        result.retain(|_, &mut c| c.value != 0);
        LaurentPoly::new(result, self.modulus)
    }

    pub fn mul(&self, other: &Self) -> Self {
        assert_eq!(self.modulus, other.modulus);
        let mut result = BTreeMap::new();
        for (&e1, &c1) in &self.terms {
            for (&e2, &c2) in &other.terms {
                let exp = e1 + e2;
                let coeff = c1 * c2;
                result
                    .entry(exp)
                    .and_modify(|c| *c = *c + coeff)
                    .or_insert(coeff);
            }
        }
        result.retain(|_, &mut c| c.value != 0);
        LaurentPoly::new(result, self.modulus)
    }

    pub fn pow(&self, exponent: &u64) -> Self {
        let mut i = exponent.clone();
        let mut poly = LaurentPoly::one(self.modulus);

        while i > 0 {
            i -= 1;
            poly = poly.mul(&self);
        }

        poly
    }

    pub fn get_coefficient(&self, index: &i64) -> ModInt {
        self.terms
            .get(index)
            .unwrap_or(&ModInt::zero(self.modulus))
            .clone()
    }

    pub fn constant_term(&self) -> ModInt {
        self.get_coefficient(&0)
    }

    pub fn lambda_reduce(&self) -> Self {
        let new_terms: BTreeMap<i64, ModInt> = self
            .terms
            .iter()
            .filter(|&(exp, _)| exp % (self.modulus as i64) == 0)
            .map(|(exp, coef)| (exp.clone() / (self.modulus as i64), coef.clone()))
            .collect();
        LaurentPoly::new(new_terms, self.modulus)
    }

    pub fn degree(&self) -> u64 {
        self.terms
            .keys()
            .max_by_key(|e| e.abs())
            .unwrap()
            .unsigned_abs()
    }
}

impl fmt::Display for LaurentPoly {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.terms.is_empty() {
            return write!(f, "0");
        }

        let mut parts = Vec::new();

        for (&exp, coeff) in self.terms.iter().rev() {
            if coeff.value == 0 {
                continue;
            }

            let coeff_str = if coeff.value == 1 && exp != 0 {
                String::new()
            } else if coeff.value == self.modulus - 1 && exp != 0 {
                "-".to_string()
            } else {
                format!("{}", coeff)
            };

            let term = match exp {
                0 => format!("{}", coeff),
                1 => format!("{}x", coeff_str),
                _ => format!("{}x^{}", coeff_str, exp),
            };

            parts.push(term);
        }

        write!(f, "{}", parts.join(" + ").replace("+ -", "- "))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_vec_and_string() {
        assert_eq!(LaurentPoly::from_vec(vec![(0, 1)], 5), LaurentPoly::one(5));
        assert_eq!(LaurentPoly::from_string("1", 5), LaurentPoly::one(5));
        assert_ne!(LaurentPoly::from_vec(vec![(0, 1)], 3), LaurentPoly::one(5));

        assert_eq!(
            LaurentPoly::from_vec(vec![(-1, 1), (0, 1), (1, 1)], 10),
            LaurentPoly::from_string("x^-1 + 1 + x", 10)
        );
        assert_eq!(
            LaurentPoly::from_vec(vec![(-1, 1), (0, 2), (2, 3)], 10),
            LaurentPoly::from_string("x^-1 + 2 + 3x^2", 10)
        );
        assert_eq!(
            LaurentPoly::from_vec(vec![(-1, 9), (0, 2), (2, 3)], 10),
            LaurentPoly::from_string("-x^-1 + 2 + 3x^2", 10)
        );
        debug_assert_eq!(
            LaurentPoly::from_vec(vec![(-1, 9), (0, 2), (2, 7)], 10),
            LaurentPoly::from_string("-x^-1 + 2 - 3x^2", 10)
        );
        assert_eq!(
            LaurentPoly::from_vec(vec![(-1, 9), (0, 8), (2, 7)], 10),
            LaurentPoly::from_string("-x^-1 - 2 - 3x^2", 10)
        );

        assert_eq!(
            LaurentPoly::from_string("-x^-1 + 2 - 3x^2", 10),
            LaurentPoly::from_string("2 - x^-1 - 3x^2", 10)
        );
        assert_eq!(
            LaurentPoly::from_string("2 - x^-1 - 3x^2", 10),
            LaurentPoly::from_string("2 - 3x^2 - x^-1", 10)
        );
        assert_eq!(
            LaurentPoly::from_string("2 - 3x^2 - x^-1", 10),
            LaurentPoly::from_string("-x^-1 - 3x^2 + 2", 10)
        );
        assert_eq!(
            LaurentPoly::from_string("-x^-1 - 3x^2 + 2", 10),
            LaurentPoly::from_string("-3x^2 - x^-1 + 2", 10)
        );

        assert_eq!(
            LaurentPoly::from_string("-3x^2 - x^-1 + 2", 10),
            LaurentPoly::from_string("-3x^2 + 2 - x^-1", 10)
        );
    }

    #[test]
    fn test_add() {
        for i in 0..10 {
            for j in 0..10 {
                for k in -10..=10 {
                    for m in 2..6 {
                        assert_eq!(
                            LaurentPoly::from_vec(vec![(k, i)], m)
                                .add(&LaurentPoly::from_vec(vec![(k, j)], m)),
                            LaurentPoly::from_vec(vec![(k, i + j)], m)
                        );
                    }
                }
            }
        }

        assert_eq!(
            LaurentPoly::from_string("2x", 10).add(&LaurentPoly::from_string("-2x^-1", 10)),
            LaurentPoly::from_string("-2x^-1 + 2x", 10)
        );
        assert_eq!(
            LaurentPoly::from_string("2x", 10).add(&LaurentPoly::from_string("3x+1-2x^-1", 10)),
            LaurentPoly::from_string("-2x^-1 + 1 + 5x", 10)
        );
    }

    #[test]
    #[should_panic]
    fn test_add_panic() {
        LaurentPoly::from_string("2x", 10).add(&LaurentPoly::from_string("2x", 9));
    }

    #[test]
    fn test_mul() {
        for i in 0..10 {
            for j in 0..10 {
                for k in -10..=10 {
                    for m in 2..6 {
                        assert_eq!(
                            LaurentPoly::from_vec(vec![(k, i)], m)
                                .mul(&LaurentPoly::from_vec(vec![(k, j)], m)),
                            LaurentPoly::from_vec(vec![(2 * k, i * j)], m)
                        );
                    }
                }
            }
        }

        assert_eq!(
            LaurentPoly::from_string("5x^3 + 3x^2 + 7 + x^-1 + 10x^-100", 11)
                .mul(&LaurentPoly::from_string("4x^100 + x^2 + 3 + 5x^-3", 11)),
            LaurentPoly::from_string(
                "9x^103 + x^102 + 6x^100 + 4x^99 + 5x^5 + 3x^4 + 4x^3 + 5x^2 + x + 9 + 7x^-1 + 2x^-3 + 5x^-4 - x^-98 + 8x^-100 + 6x^-103",
                11
            )
        );
    }

    #[test]
    #[should_panic]
    fn test_mul_panic() {
        LaurentPoly::from_string("2x", 10).mul(&LaurentPoly::from_string("2x", 9));
    }

    #[test]
    fn test_pow() {
        assert_eq!(
            LaurentPoly::from_string("2x + 1", 10).pow(&3),
            LaurentPoly::from_string("8x^3 + 2x^2 + 6x + 1", 10)
        );
        assert_eq!(
            LaurentPoly::from_string("2x + 1", 10).pow(&0),
            LaurentPoly::from_string("1", 10)
        );
        assert_eq!(
            LaurentPoly::from_string("5x^3 + 3x^2 + 7 + x^-1 + 10x^-100", 11).pow(&3),
            LaurentPoly::from_string(
                "4x^9 + 5x^8 + 3x^7 + 2x^6 + x^5 + 4x^4 + 3x^3 + 2x^2 + 9x + 4x^-1 - x^-2 + x^-3 + 2x^-94 + 9x^-95 + 6x^-96 - x^-97 + 9x^-98 + 4x^-99 + 7x^-100 + 2x^-101 + 8x^-102 + 4x^-197 + 9x^-198 - x^-200 + 3x^-201 - x^-300",
                11
            )
        );
    }

    #[test]
    fn test_get_coefficient() {
        let poly = LaurentPoly::from_string("5x^3 + 3x^2 + 7 + x^-1 + 10x^-5", 11);
        let coeffs: Vec<ModInt> = vec![5, 3, 0, 7, 1, 0, 0, 0, 10]
            .iter()
            .map(|x: &u64| ModInt::new(x.clone(), 11))
            .collect();
        for i in 0..8 {
            assert_eq!(poly.get_coefficient(&(i + 4)), ModInt::zero(11));
            assert_eq!(poly.get_coefficient(&(3 - i)), coeffs[i as usize]);
        }
    }

    #[test]
    fn test_constant_term() {
        assert_eq!(
            LaurentPoly::from_string("5x^3 + 3x^2 + 7 + x^-1 + 10x^-5", 11).constant_term(),
            ModInt::new(7, 11)
        );
        assert_eq!(
            LaurentPoly::from_string("5x^3 + 3x^2 + 7 + x^-1 + 10x^-100", 11)
                .pow(&3)
                .constant_term(),
            ModInt::zero(11)
        );
        assert_eq!(LaurentPoly::zero(11).constant_term(), ModInt::zero(11));
        assert_eq!(LaurentPoly::one(11).constant_term(), ModInt::new(1, 11));
    }

    #[test]
    fn test_lambda_reduce() {
        assert_eq!(
            LaurentPoly::from_string("5x^3 + 3x^2 + 7 + x^-1 + 10x^-110", 11).lambda_reduce(),
            LaurentPoly::from_string("7 - x^-10", 11)
        );
        assert_eq!(
            LaurentPoly::from_string("4x^22 + x^2 + 3 + 4x^-5 + 5x^-11", 11).lambda_reduce(),
            LaurentPoly::from_string("4x^2 + 3 + 5x^-1", 11)
        );
    }

    #[test]
    fn test_degree() {
        assert_eq!(
            LaurentPoly::from_string("5x^3 + 3x^2 + 7 + x^-1 + 10x^-110", 11).degree(),
            110
        );
        assert_eq!(
            LaurentPoly::from_string("5x^3 + 3x^2 + 7 + x^-1", 11).degree(),
            3
        );
        assert_eq!(LaurentPoly::from_string("x + 1 + x^-1", 11).degree(), 1);
        assert_eq!(LaurentPoly::from_string("7", 11).degree(), 0);
    }

    #[test]
    fn test_to_from_string() {
        let poly = LaurentPoly::from_string("5x^3 + 3x^2 + 7 + x^-1 + 10x^-110", 11);
        assert_eq!(
            LaurentPoly::from_string(format!("{}", poly).as_str(), 11),
            poly
        );
    }
}
