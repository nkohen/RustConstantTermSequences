use crate::laurent_poly::LaurentPoly;
use crate::lin_rep::LinRep;
use crate::mod_int::ModInt;
use crate::mod_int_vector::ModIntVector;
use graphviz_rust::cmd::{CommandArg, Format};
use graphviz_rust::dot_structures::Graph;
use graphviz_rust::printer::PrinterContext;
use graphviz_rust::{exec, parse};
use std::collections::HashMap;
use std::fmt::Display;
use std::hash::Hash;

pub struct DFAO<A: Eq + Hash, S: Clone + Eq + Hash> {
    pub states: Vec<S>, // states[0] is the initial state
    pub transitions: HashMap<(S, A), S>,
}

pub trait ConstantTerm {
    fn constant_term(&self) -> ModInt;
    fn modulus(&self) -> u64;
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
}

impl DFAO<ModInt, ModIntVector> {
    // TODO: Fix this function, then make its reverse version, then make constructions which terminate on property
    pub fn lin_rep_to_machine(P: &LaurentPoly, Q: &LaurentPoly) -> Self {
        assert_eq!(P.modulus, Q.modulus);
        let p = P.modulus;
        let lin_rep = LinRep::for_ct_sequence(P, Q);
        let mut states: Vec<ModIntVector> = Vec::new();
        states.push(lin_rep.row_vec);
        let mut k = 0;
        let mut transitions = HashMap::new();

        while k < states.len() {
            let current_state = states.get(k).unwrap().clone();
            for i in 0..p {
                let new_state = lin_rep.mat_func[i as usize].right_mul(&current_state);
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
}

impl<S: ConstantTerm + Clone + Eq + Hash> DFAO<ModInt, S> {
    pub fn compute_ct(&self, n: u64) -> ModInt {
        let p = self.states.get(0).unwrap().modulus();
        let mut n = n;
        let mut digits: Vec<ModInt> = Vec::new();
        while n > 0 {
            let r = n % p;
            digits.push(ModInt::new(r, p));
            n = (n - r) / p;
        }

        fn output_func<S: ConstantTerm>(poly: S) -> ModInt {
            poly.constant_term()
        }

        self.compute(digits, &output_func)
    }

    pub fn serialize(&self) -> String {
        let p = self.states.get(0).unwrap().modulus();
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

impl<S: ConstantTerm + Clone + Eq + Hash + Display> DFAO<ModInt, S> {
    pub fn to_graphviz(&self) -> String {
        let p = self.states.get(0).unwrap().modulus();
        let mut str = String::from("digraph G {\nrankdir = LR;\nnode [shape = point ]; qi\n");
        let mut index_map: HashMap<S, usize> = HashMap::new();

        for i in 0..self.states.len() {
            let state = self.states.get(i).unwrap();
            index_map.insert(state.clone(), i);
            let shape = if self.states[i].constant_term() == ModInt::zero(p) {
                "circle"
            } else {
                "doublecircle"
            };
            str.push_str(&format!(
                "node [shape = {shape}, label=\"{state}\", fontsize=12]{i};\n"
            ));
        }

        str.push_str("qi -> 0;\n");

        for i in 0..self.states.len() {
            let state = self.states.get(i).unwrap();
            let mut transitions_map: HashMap<usize, Vec<String>> = HashMap::new();
            for j in 0..p {
                let to_state = self
                    .transitions
                    .get(&(state.clone(), ModInt::new(j, p)))
                    .unwrap();
                let to_index = index_map.get(to_state).unwrap();
                if transitions_map.contains_key(to_index) {
                    transitions_map
                        .get_mut(to_index)
                        .unwrap()
                        .push(j.to_string());
                } else {
                    transitions_map.insert(to_index.clone(), vec![j.to_string()]);
                }
            }
            for to_index in transitions_map.keys() {
                str.push_str(&format!(
                    "{i} -> {to_index} [label = \"{}\"];\n",
                    transitions_map.get(to_index).unwrap().join(", ")
                ));
            }
        }

        str.push('}');
        str
    }

    pub fn save_png(&self, filename: &str) -> () {
        let g: Graph = parse(&self.to_graphviz()).unwrap();
        let _ = exec(
            g,
            &mut PrinterContext::default(),
            vec![Format::Png.into(), CommandArg::Output(filename.to_string())],
        );
    }
}
