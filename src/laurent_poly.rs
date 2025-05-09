use crate::mod_int::ModInt;
use std::collections::BTreeMap;
use std::fmt;

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct LaurentPoly {
    terms: BTreeMap<i64, ModInt>, // exponent -> coefficient
    pub modulus: u64
}

impl LaurentPoly {
    
    pub fn new(terms: BTreeMap<i64, ModInt>, modulus: u64) -> LaurentPoly {
        LaurentPoly {
            terms,
            modulus
        }
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
                        .trim()
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
                let coeff = coeff_str.trim().parse::<i64>().unwrap();
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
