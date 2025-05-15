use crate::laurent_poly::LaurentPoly;
use crate::lin_rep::LinRep;
use crate::mod_int::ModInt;
use crate::mod_int_vector::ModIntVector;
use either::{Either, Left, Right};
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

impl<A: Eq + Hash, S: Clone + Eq + Hash> DFAO<A, S> {
    pub fn evaluate<F, O>(&self, input: Vec<A>, output_func: F) -> O
    where
        F: FnOnce(&S) -> O,
    {
        let mut state = self.states.get(0).unwrap();
        for character in input {
            state = self.transitions.get(&(state.clone(), character)).unwrap();
        }
        output_func(state)
    }
}

impl<S: Clone + Eq + Hash> DFAO<ModInt, S> {
    pub fn from_reduction_rules_until_prop<F1, F2>(
        initial_state: &S,
        modulus: u64,
        reduction_rule: F1,
        stop_prop: F2,
    ) -> Either<Self, Vec<u64>>
    where
        F1: Fn(&S, ModInt) -> S,
        F2: Fn(&S) -> bool,
    {
        if stop_prop(initial_state) {
            return Right(vec![]);
        }

        let mut states_with_paths: Vec<(S, Vec<u64>)> = Vec::new();
        states_with_paths.push((initial_state.clone(), vec![]));
        let mut k = 0;
        let mut transitions = HashMap::new();

        while k < states_with_paths.len() {
            let (current_state, path) = states_with_paths.get(k).unwrap().clone();
            for i in 0..modulus {
                let new_state = reduction_rule(&current_state, ModInt::new(i, modulus));
                let mut new_state_index = states_with_paths.len();
                for j in 0..states_with_paths.len() {
                    let (state_j, _) = states_with_paths.get(j).unwrap();
                    if state_j == &new_state {
                        new_state_index = j;
                        break;
                    }
                }

                if new_state_index == states_with_paths.len() {
                    let mut new_path = path.clone();
                    new_path.push(i);

                    if stop_prop(&new_state) {
                        return Right(new_path);
                    }
                    states_with_paths.push((new_state, new_path));
                }

                transitions.insert(
                    (current_state.clone(), ModInt::new(i, modulus)),
                    match states_with_paths.get(new_state_index).unwrap() {
                        (state, _) => state.clone(),
                    },
                );
            }
            k += 1;
        }

        let states = states_with_paths
            .iter()
            .map(|(state, _)| state.clone())
            .collect();

        Left(DFAO {
            states,
            transitions,
        })
    }

    pub fn from_reduction_rules<F>(initial_state: &S, modulus: u64, reduction_rule: F) -> Self
    where
        F: Fn(&S, ModInt) -> S,
    {
        Self::from_reduction_rules_until_prop(initial_state, modulus, reduction_rule, |_| false)
            .unwrap_left()
    }

    pub fn compute_lsd<F, O>(&self, n: u64, modulus: u64, output_func: F) -> O
    where
        F: Fn(&S) -> O,
    {
        let digits = ModInt::get_digits(n, modulus);
        self.evaluate(digits, output_func)
    }

    pub fn compute_msd<F, O>(&self, n: u64, modulus: u64, output_func: F) -> O
    where
        F: Fn(&S) -> O,
    {
        let mut digits = ModInt::get_digits(n, modulus);
        digits.reverse();
        self.evaluate(digits, output_func)
    }
}

impl DFAO<ModInt, LaurentPoly> {
    pub fn poly_auto(P: &LaurentPoly, Q: &LaurentPoly) -> Self {
        assert_eq!(P.modulus, Q.modulus);
        let p = P.modulus;
        DFAO::from_reduction_rules(&Q, p, |state, i| P.pow(&i.value).mul(state).lambda_reduce())
    }

    pub fn poly_auto_fail_on_zero(P: &LaurentPoly, Q: &LaurentPoly) -> Option<Self> {
        assert_eq!(P.modulus, Q.modulus);
        let p = P.modulus;
        match Self::from_reduction_rules_until_prop(
            &Q,
            p,
            |state, i| P.pow(&i.value).mul(state).lambda_reduce(),
            |state| state.constant_term() == ModInt::zero(p),
        ) {
            Left(machine) => Some(machine),
            Right(_) => None,
        }
    }

    pub fn compute_ct(&self, n: u64) -> ModInt {
        let p = self.states.first().unwrap().modulus;
        DFAO::compute_lsd(&self, n, p, |poly| poly.constant_term())
    }
}

impl DFAO<ModInt, ModIntVector> {
    pub fn lin_rep_machine(P: &LaurentPoly, Q: &LaurentPoly) -> Self {
        let lin_rep = LinRep::for_ct_sequence(P, Q);
        assert_eq!(P.modulus, Q.modulus);
        let p = P.modulus;

        DFAO::from_reduction_rules(&lin_rep.row_vec, p, |state, i| {
            lin_rep.mat_func[i.value as usize].right_mul(&state)
        })
    }

    pub fn lin_rep_reverse_machine(P: &LaurentPoly, Q: &LaurentPoly) -> Self {
        assert_eq!(P.modulus, Q.modulus);
        let p = P.modulus;
        let lin_rep = LinRep::for_ct_sequence(P, Q);

        DFAO::from_reduction_rules(&lin_rep.col_vec, p, |state, i| {
            lin_rep.mat_func[i.value as usize].left_mul(&state)
        })
    }

    pub fn compute_ct(&self, n: u64) -> ModInt {
        let p = self.states.first().unwrap().modulus;
        DFAO::compute_lsd(&self, n, p, |poly| poly.constant_term())
    }

    pub fn compute_ct_reverse(&self, n: u64, Q: &LaurentPoly) -> ModInt {
        let p = self.states.first().unwrap().modulus;

        self.compute_msd(n, p, |state| {
            state.dot(&ModIntVector::from_poly(Q, state.dim))
        })
    }

    pub fn compute_shortest_prop<F>(P: &LaurentPoly, Q: &LaurentPoly, prop: F) -> Option<u64>
    where
        F: Fn(&ModIntVector) -> bool,
    {
        assert_eq!(P.modulus, Q.modulus);
        let p = P.modulus;
        let lin_rep = LinRep::for_ct_sequence(P, Q);
        match DFAO::from_reduction_rules_until_prop(
            &lin_rep.col_vec,
            p,
            |state, i| lin_rep.mat_func[i.value as usize].left_mul(&state),
            prop,
        ) {
            Left(_) => None,
            Right(digits) => {
                let mut result = 0;
                for digit in digits {
                    result *= p;
                    result += digit;
                }
                Some(result)
            }
        }
    }

    pub fn compute_shortest_ct_prop<F>(P: &LaurentPoly, Q: &LaurentPoly, prop: F) -> Option<u64>
    where
        F: Fn(&ModInt) -> bool,
    {
        Self::compute_shortest_prop(P, Q, |vec| {
            prop(&vec.dot(&ModIntVector::from_poly(Q, vec.dim)))
        })
    }

    // TODO: Test this
    pub fn compute_shortest_zero(P: &LaurentPoly, Q: &LaurentPoly) -> Option<u64> {
        Self::compute_shortest_ct_prop(P, Q, |value| value == &ModInt::zero(P.modulus))
    }
}

impl<S: Clone + Eq + Hash> DFAO<ModInt, S> {
    pub fn serialize<F>(&self, p: u64, output_func: F) -> String
    where
        F: Fn(&S) -> ModInt,
    {
        let mut str = String::from(format!("lsd_{p}\n\n"));
        for k in 0..self.states.len() {
            let current_state = self.states.get(k).unwrap();
            str.push_str(&format!("{k} {}\n", output_func(&current_state)));
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

impl<S: Clone + Eq + Hash + Display> DFAO<ModInt, S> {
    pub fn to_graphviz(&self, p: u64) -> String {
        let mut str = String::from("digraph G {\nrankdir = LR;\nnode [shape = point ]; qi\n");
        let mut index_map: HashMap<S, usize> = HashMap::new();

        for i in 0..self.states.len() {
            let state = self.states.get(i).unwrap();
            index_map.insert(state.clone(), i);
            str.push_str(&format!(
                "node [shape = circle, label=\"{state}\", fontsize=12]{i};\n"
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

    pub fn save_png(&self, p: u64, filename: &str) -> () {
        let g: Graph = parse(&self.to_graphviz(p)).unwrap();
        let _ = exec(
            g,
            &mut PrinterContext::default(),
            vec![Format::Png.into(), CommandArg::Output(filename.to_string())],
        );
    }
}
