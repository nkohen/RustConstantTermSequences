use crate::laurent_poly::LaurentPoly;
use crate::mod_int::ModInt;
use std::collections::HashMap;
use std::hash::Hash;

pub struct DFAO<A: Eq + Hash, S: Clone + Eq + Hash> {
    pub states: Vec<S>, // states[0] is the initial state
    pub transitions: HashMap<(S, A), S>,
}

impl<A: Eq + Hash, S: Clone + Eq + Hash> DFAO<A, S> {
    pub fn compute<O>(&self, input: Vec<A>, output_func: &dyn Fn(S) -> O) -> O {
        let mut state = self.states.get(0).unwrap().clone();
        for character in input {
            state = self
                .transitions
                .get(&(state.clone(), character))
                .unwrap()
                .clone();
        }
        output_func(state)
    }
}

impl DFAO<ModInt, LaurentPoly> {
    pub fn poly_auto(P: &LaurentPoly, Q: &LaurentPoly) -> Self {
        assert_eq!(P.modulus, Q.modulus);
        let p = P.modulus;
        let mut states: Vec<LaurentPoly> = Vec::new();
        states.push(Q.clone());
        let mut k = 0;
        let mut transitions = HashMap::new();

        while k < states.len() {
            let current_state = states.get(k).unwrap().clone();
            for i in 0..p {
                let new_state = P.pow(&i).mul(&current_state).lambda_reduce();
                let mut new_state_index = states.len();
                for j in 0..states.len() {
                    if states.get(j).unwrap() == &new_state {
                        new_state_index = j;
                        break;
                    }
                }

                if new_state_index == states.len() {
                    states.push(new_state);
                }

                transitions.insert(
                    (current_state.clone(), ModInt::new(i, p)),
                    states.get(new_state_index).unwrap().clone(),
                );
            }
            k += 1;
        }

        DFAO {
            states,
            transitions,
        }
    }

    pub fn compute_ct(&self, n: u64) -> ModInt {
        let p = self.states.get(0).unwrap().modulus;
        let mut n = n;
        let mut digits: Vec<ModInt> = Vec::new();
        while n > 0 {
            let r = n % p;
            digits.push(ModInt::new(r, p));
            n = (n - r) / p;
        }

        fn output_func(poly: LaurentPoly) -> ModInt {
            poly.constant_term()
        }

        self.compute(digits, &output_func)
    }

    pub fn serialize(&self) -> String {
        let p = self.states.get(0).unwrap().modulus;
        let mut str = String::from(format!("lsd_{p}\n\n"));
        for k in 0..self.states.len() {
            let current_state = self.states.get(k).unwrap();
            str.push_str(&format!("{k} {}\n", current_state.constant_term()));
            for i in 0..p {
                let index = ModInt::new(i, p);
                let to_state = self
                    .transitions
                    .get(&(current_state.clone(), index))
                    .unwrap();
                let to_index = self
                    .states
                    .iter()
                    .position(|state| state == to_state)
                    .unwrap();
                str.push_str(&format!("{i} -> {to_index}\n"));
            }
            str.push('\n');
        }
        str
    }
}
