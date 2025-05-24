use crate::laurent_poly::LaurentPoly;
use crate::lin_rep::LinRep;
use crate::mod_int::ModInt;
use crate::mod_int_vector::ModIntVector;
use crate::sequences::constant_term_reduce;
use either::{Either, Left, Right};
use graphviz_rust::cmd::{CommandArg, Format};
use graphviz_rust::dot_structures::Graph;
use graphviz_rust::printer::PrinterContext;
use graphviz_rust::{exec, parse};
use std::collections::HashMap;
use std::fmt::Display;
use std::hash::Hash;
use std::sync::atomic::AtomicBool;
use std::sync::{Arc, mpsc};
use std::thread;
use std::time::Duration;

#[derive(Debug)]
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
        state_bound: usize,
        cancel_flag_opt: Option<Arc<AtomicBool>>,
    ) -> Result<Either<Self, Vec<u64>>, String>
    where
        F1: Fn(&S, ModInt) -> S,
        F2: Fn(&S) -> bool,
    {
        if stop_prop(initial_state) {
            return Ok(Right(vec![]));
        }

        let mut states_with_paths: Vec<(S, Vec<u64>)> = Vec::new();
        states_with_paths.push((initial_state.clone(), vec![]));
        let mut k = 0;
        let mut transitions = HashMap::new();

        while k < states_with_paths.len() {
            if let Some(cancel_flag) = &cancel_flag_opt {
                if cancel_flag.load(std::sync::atomic::Ordering::Relaxed) {
                    return Err("Process was cancelled".to_string());
                }
            }

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
                        return Ok(Right(new_path));
                    }
                    states_with_paths.push((new_state, new_path));
                    if states_with_paths.len() > state_bound {
                        return Err(format!("Number of states exceeded {}.", state_bound));
                    }
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

        Ok(Left(DFAO {
            states,
            transitions,
        }))
    }

    pub fn from_reduction_rules<F>(
        initial_state: &S,
        modulus: u64,
        reduction_rule: F,
        state_bound: usize,
        cancel_flag_opt: Option<Arc<AtomicBool>>,
    ) -> Result<Self, String>
    where
        F: Fn(&S, ModInt) -> S,
    {
        Self::from_reduction_rules_until_prop(
            initial_state,
            modulus,
            reduction_rule,
            |_| false,
            state_bound,
            cancel_flag_opt,
        )
        .map(|machine| machine.unwrap_left())
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
    pub fn poly_auto(P: &LaurentPoly, Q: &LaurentPoly, state_bound: usize) -> Result<Self, String> {
        assert_eq!(P.modulus, Q.modulus);
        let p = P.modulus;
        DFAO::from_reduction_rules(
            &Q,
            p,
            |state, i| P.pow(&i.value).mul(state).lambda_reduce(),
            state_bound,
            None,
        )
    }

    pub fn poly_auto_fail_on_prop<F>(
        P: &LaurentPoly,
        Q: &LaurentPoly,
        prop: F,
        state_bound: usize,
        cancel_flag_opt: Option<Arc<AtomicBool>>,
    ) -> Result<Option<Self>, String>
    where
        F: Fn(&LaurentPoly) -> bool,
    {
        assert_eq!(P.modulus, Q.modulus);
        Self::from_reduction_rules_until_prop(
            &Q,
            P.modulus,
            |state, i| P.pow(&i.value).mul(state).lambda_reduce(),
            prop,
            state_bound,
            cancel_flag_opt,
        )
        .map(|machine_or_zero| match machine_or_zero {
            Left(machine) => Some(machine),
            Right(_) => None,
        })
    }

    pub fn poly_auto_fail_on_zero(
        P: &LaurentPoly,
        Q: &LaurentPoly,
        state_bound: usize,
    ) -> Result<Option<Self>, String> {
        Self::poly_auto_fail_on_prop(
            P,
            Q,
            |state| state.constant_term() == ModInt::zero(P.modulus),
            state_bound,
            None,
        )
    }

    pub fn compute_ct(&self, n: u64) -> ModInt {
        let p = self.states.first().unwrap().modulus;
        self.compute_lsd(n, p, |poly| poly.constant_term())
    }
}

impl DFAO<ModInt, ModIntVector> {
    pub fn lin_rep_machine(
        P: &LaurentPoly,
        Q: &LaurentPoly,
        state_bound: usize,
    ) -> Result<Self, String> {
        let lin_rep = LinRep::for_ct_sequence(P, Q);
        assert_eq!(P.modulus, Q.modulus);
        let p = P.modulus;

        DFAO::from_reduction_rules(
            &lin_rep.row_vec,
            p,
            |state, i| lin_rep.mat_func[i.value as usize].right_mul(&state),
            state_bound,
            None,
        )
    }

    pub fn lin_rep_reverse_machine(
        P: &LaurentPoly,
        Q: &LaurentPoly,
        state_bound: usize,
    ) -> Result<Self, String> {
        assert_eq!(P.modulus, Q.modulus);
        let p = P.modulus;
        let lin_rep = LinRep::for_ct_sequence(P, Q);

        DFAO::from_reduction_rules(
            &lin_rep.col_vec,
            p,
            |state, i| lin_rep.mat_func[i.value as usize].left_mul(&state),
            state_bound,
            None,
        )
    }

    pub fn compute_ct(&self, n: u64) -> ModInt {
        let p = self.states.first().unwrap().modulus;
        self.compute_lsd(n, p, |poly| poly.constant_term())
    }

    pub fn compute_ct_reverse(&self, n: u64, Q: &LaurentPoly) -> ModInt {
        let p = self.states.first().unwrap().modulus;

        self.compute_msd(n, p, |state| {
            state.dot(&ModIntVector::from_poly(Q, state.dim))
        })
    }

    // This function will not terminate without using cancel_flag if prop is never satisfied
    pub fn compute_shortest_poly_prop_directly<F>(
        P: &LaurentPoly,
        Q: &LaurentPoly,
        prop: F,
        cancel_flag: Arc<AtomicBool>,
    ) -> Option<Option<u64>>
    where
        F: Fn(&LaurentPoly) -> bool,
    {
        assert_eq!(P.modulus, Q.modulus);
        let mut n = 0;
        loop {
            if cancel_flag.load(std::sync::atomic::Ordering::Relaxed) {
                return None;
            }

            if prop(&constant_term_reduce(&P, &Q, &n)) {
                return Some(Some(n));
            }

            n += 1;
        }
    }

    pub fn compute_shortest_poly_prop<F>(
        P: &LaurentPoly,
        Q: &LaurentPoly,
        prop: F,
        state_bound: usize,
        cancel_flag_opt: Option<Arc<AtomicBool>>,
    ) -> Result<Option<u64>, String>
    where
        F: Fn(&LaurentPoly) -> bool + Send + Sync + 'static,
    {
        let cancel_flag = cancel_flag_opt.unwrap_or(Arc::new(AtomicBool::new(false)));
        let (tx, rx) = mpsc::channel();
        let prop = Arc::new(prop);
        let P = Arc::new(P.clone());
        let Q = Arc::new(Q.clone());

        // This thread directly compute the sequence of polynomials
        // This thread returns the value if there is a value for which prop is satisfied
        {
            let tx = tx.clone();
            let P = Arc::clone(&P);
            let Q = Arc::clone(&Q);
            let prop = Arc::clone(&prop);
            let flag = Arc::clone(&cancel_flag);

            thread::spawn(move || {
                if let Some(result) =
                    Self::compute_shortest_poly_prop_directly(&P, &Q, |v| prop(v), flag)
                {
                    let _ = tx.send(Ok(result));
                }
            });
        }

        // This thread computes the lsd-DFAO, but only returns if it is completed and no value satisfies prop, or state_bound was exceeded
        // Having a shortest-path to a lsd-DFAO state does not guarantee that this is the smallest value
        {
            let tx = tx.clone();
            let P = Arc::clone(&P);
            let Q = Arc::clone(&Q);
            let prop = Arc::clone(&prop);
            let flag = Arc::clone(&cancel_flag);

            thread::spawn(move || {
                match DFAO::poly_auto_fail_on_prop(
                    &P,
                    &Q,
                    |poly| prop(poly),
                    state_bound,
                    Some(flag),
                ) {
                    Ok(Some(_)) => {
                        let _ = tx.send(Ok(None));
                    }
                    Ok(None) => {} // The other thread handles this case
                    Err(err_msg) => {
                        if err_msg != "Process was cancelled" {
                            let _ = tx.send(Err(err_msg));
                        };
                    }
                }
            });
        }

        // This thread makes this method return if cancel_flag is set to true
        {
            let tx = tx.clone();
            let flag = Arc::clone(&cancel_flag);
            thread::spawn(move || {
                loop {
                    thread::sleep(Duration::from_secs(2));
                    if flag.load(std::sync::atomic::Ordering::Relaxed) {
                        let _ = tx.send(Err("Process was cancelled".to_string()));
                        break;
                    }
                }
            });
        }

        // Receive the first result
        let result = rx.recv().unwrap();

        // Signal cancellation to the other threads
        cancel_flag.store(true, std::sync::atomic::Ordering::Relaxed);

        result
    }

    // This function will not terminate without using cancel_flag if prop is never satisfied
    pub fn compute_shortest_functional_prop_directly<F>(
        lin_rep: &LinRep,
        prop: F,
        cancel_flag: Arc<AtomicBool>,
    ) -> Option<Option<u64>>
    where
        F: Fn(&ModIntVector) -> bool,
    {
        let mut n = 0;
        loop {
            if cancel_flag.load(std::sync::atomic::Ordering::Relaxed) {
                return None;
            }

            if prop(&lin_rep.compute_functional(n)) {
                return Some(Some(n));
            }

            n += 1;
        }
    }

    pub fn compute_shortest_prop_using_msd_dfao<F>(
        lin_rep: &LinRep,
        prop: F,
        state_bound: usize,
        cancel_flag_opt: Option<Arc<AtomicBool>>,
    ) -> Result<Option<u64>, String>
    where
        F: Fn(&ModIntVector) -> bool,
    {
        let p = lin_rep.modulus;
        DFAO::from_reduction_rules_until_prop(
            &lin_rep.col_vec,
            p,
            |state, i| lin_rep.mat_func[i.value as usize].left_mul(&state),
            prop,
            state_bound,
            cancel_flag_opt,
        )
        .map(|machine_or_first| match machine_or_first {
            Left(_) => None,
            Right(first_digits) => {
                let mut result = 0;
                for digit in first_digits {
                    result *= p;
                    result += digit;
                }
                Some(result)
            }
        })
    }

    pub fn compute_shortest_functional_prop<F>(
        P: &LaurentPoly,
        Q: &LaurentPoly,
        prop: F,
        state_bound: usize,
        cancel_flag_opt: Option<Arc<AtomicBool>>,
    ) -> Result<Option<u64>, String>
    where
        F: Fn(&ModIntVector) -> bool + Send + Sync + 'static,
    {
        let lin_rep = Arc::new(LinRep::for_ct_sequence(P, Q));
        let cancel_flag = cancel_flag_opt.unwrap_or(Arc::new(AtomicBool::new(false)));
        let (tx, rx) = mpsc::channel();
        let prop = Arc::new(prop);

        // This thread computes values of the sequence of functionals directly and checks prop
        // This thread usually completes first if there is a value for which prop is satisfied
        {
            let tx = tx.clone();
            let lin_rep = Arc::clone(&lin_rep);
            let prop = Arc::clone(&prop);
            let flag = Arc::clone(&cancel_flag);

            thread::spawn(move || {
                if let Some(result) =
                    Self::compute_shortest_functional_prop_directly(&lin_rep, |v| prop(v), flag)
                {
                    let _ = tx.send(Ok(result));
                }
            });
        }

        // This thread computes the msd-DFAO and returns the value of the shortest path to a prop-satisfying state
        {
            let tx = tx.clone();
            let lin_rep = Arc::clone(&lin_rep);
            let prop = Arc::clone(&prop);
            let flag = Arc::clone(&cancel_flag);

            thread::spawn(move || {
                let maybe_result = Self::compute_shortest_prop_using_msd_dfao(
                    &lin_rep,
                    |v| prop(v),
                    state_bound,
                    Some(flag),
                );
                if let Ok(result) = maybe_result {
                    let _ = tx.send(Ok(result));
                } else if let Err(err_msg) = maybe_result {
                    if err_msg != "Process was cancelled" {
                        let _ = tx.send(Err(err_msg));
                    }
                }
            });
        }

        // This thread makes this method return if cancel_flag is set to true
        {
            let tx = tx.clone();
            let flag = Arc::clone(&cancel_flag);
            thread::spawn(move || {
               loop {
                   thread::sleep(Duration::from_secs(2));
                   if flag.load(std::sync::atomic::Ordering::Relaxed) {
                       let _ = tx.send(Err("Process was cancelled".to_string()));
                       break;
                   } 
               }
            });
        }

        // Receive the first result
        let result = rx.recv().unwrap();

        // Signal cancellation to the other threads
        cancel_flag.store(true, std::sync::atomic::Ordering::Relaxed);

        result
    }

    pub fn compute_shortest_ct_prop<F>(
        P: &LaurentPoly,
        Q: &LaurentPoly,
        prop: F,
        state_bound: usize,
    ) -> Result<Option<u64>, String>
    where
        F: Fn(&ModInt) -> bool + Send + Sync + 'static,
    {
        let cancel_flag = Arc::new(AtomicBool::new(false));
        let (tx, rx) = mpsc::channel();
        let prop = Arc::new(prop);
        let P = Arc::new(P.clone());
        let Q = Arc::new(Q.clone());

        // This thread finds the first value satisfying prop using msd-first
        {
            let tx = tx.clone();
            let P = Arc::clone(&P);
            let Q = Arc::clone(&Q);
            let prop = Arc::clone(&prop);
            let flag = Arc::clone(&cancel_flag);

            thread::spawn(move || {
                let result = Self::compute_shortest_functional_prop(
                    &P,
                    &Q,
                    move |v| prop(&v.constant_term()),
                    state_bound,
                    Some(flag),
                );
                let _ = tx.send(result);
            });
        }

        // This thread finds the first value satisfying prop using lsd-first
        {
            let tx = tx.clone();
            let P = Arc::clone(&P);
            let Q = Arc::clone(&Q);
            let prop = Arc::clone(&prop);
            let flag = Arc::clone(&cancel_flag);

            thread::spawn(move || {
                let result = Self::compute_shortest_poly_prop(
                    &P,
                    &Q,
                    move |poly| prop(&poly.constant_term()),
                    state_bound,
                    Some(flag),
                );
                let _ = tx.send(result);
            });
        }

        // Receive the first result
        let result = rx.recv().unwrap();

        // Signal cancellation to the other threads
        cancel_flag.store(true, std::sync::atomic::Ordering::Relaxed);

        result
    }

    pub fn compute_shortest_zero(
        P: &LaurentPoly,
        Q: &LaurentPoly,
        state_bound: usize,
    ) -> Result<Option<u64>, String> {
        let P = P.clone();
        Self::compute_shortest_ct_prop(
            &P,
            Q,
            move |value| value == &ModInt::zero(P.modulus),
            state_bound,
        )
    }

    pub fn compute_shortest_non_zero(
        P: &LaurentPoly,
        Q: &LaurentPoly,
        state_bound: usize,
    ) -> Result<Option<u64>, String> {
        let P = P.clone();
        Self::compute_shortest_ct_prop(
            &P,
            Q,
            move |value| value != &ModInt::zero(P.modulus),
            state_bound,
        )
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_poly_auto() {
        let primes = vec![2, 3, 5, 7, 11, 13];
        let mots: Vec<u64> = vec![
            1, 1, 2, 4, 9, 21, 51, 127, 323, 835, 2188, 5798, 15511, 41835, 113634, 310572, 853467,
            2356779, 6536382, 18199284, 50852019, 142547559, 400763223, 1129760415, 3192727797,
        ];

        for p in primes {
            let P = LaurentPoly::from_string("x + 1 + x^-1", p);
            let Q = LaurentPoly::from_string("1 - x^2", p);
            let dfao = DFAO::poly_auto(&P, &Q, 10000).unwrap();
            for n in 0..mots.len() as u64 {
                assert_eq!(dfao.compute_ct(n), ModInt::new(mots[n as usize], p));
            }
        }

        let P = LaurentPoly::from_string("x + 1 + x^-1", 11);
        let Q = LaurentPoly::from_string("1 - x^2", 11);
        let dfao = DFAO::poly_auto(&P, &Q, 15);
        assert_eq!(dfao.unwrap_err(), "Number of states exceeded 15.");
    }

    #[test]
    fn test_poly_auto_fail_on_zero() {
        let primes = vec![2, 3, 5, 7, 11, 13, 17, 19, 23, 29];
        let no_zero_primes = vec![2, 5, 11, 13, 23, 29];

        for p in primes {
            let P = LaurentPoly::from_string("x + 1 + x^-1", p);
            let Q = LaurentPoly::one(p);
            let dfao_opt = DFAO::poly_auto_fail_on_zero(&P, &Q, 10000).unwrap();

            if no_zero_primes.contains(&p) {
                let dfao = dfao_opt.unwrap();
                for n in 0..25 {
                    assert_eq!(dfao.compute_ct(n), P.pow(&n).constant_term());
                }
            } else {
                assert!(dfao_opt.is_none());
            }
        }
    }

    #[test]
    fn test_lin_rep_machine() {
        let primes = vec![2, 3, 5, 7, 11, 13];
        let mots: Vec<u64> = vec![
            1, 1, 2, 4, 9, 21, 51, 127, 323, 835, 2188, 5798, 15511, 41835, 113634, 310572, 853467,
            2356779, 6536382, 18199284, 50852019, 142547559, 400763223, 1129760415, 3192727797,
        ];

        for p in primes {
            let P = LaurentPoly::from_string("x + 1 + x^-1", p);
            let Q = LaurentPoly::from_string("1 - x^2", p);
            let dfao = DFAO::lin_rep_machine(&P, &Q, 10000).unwrap();
            for n in 0..mots.len() as u64 {
                assert_eq!(dfao.compute_ct(n), ModInt::new(mots[n as usize], p));
            }

            let dfao_poly_auto = DFAO::poly_auto(&P, &Q, 10000).unwrap();
            assert_eq!(
                dfao.serialize(p, |state| state.constant_term()),
                dfao_poly_auto.serialize(p, |state| state.constant_term())
            );
        }

        let P = LaurentPoly::from_string("x + 1 + x^-1", 11);
        let Q = LaurentPoly::from_string("1 - x^2", 11);
        let dfao = DFAO::lin_rep_machine(&P, &Q, 15);
        assert_eq!(dfao.unwrap_err(), "Number of states exceeded 15.");
    }

    #[test]
    fn test_lin_rep_reverse_machine() {
        let primes = vec![2, 3, 5, 7, 11, 13];
        let mots: Vec<u64> = vec![
            1, 1, 2, 4, 9, 21, 51, 127, 323, 835, 2188, 5798, 15511, 41835, 113634, 310572, 853467,
            2356779, 6536382, 18199284, 50852019, 142547559, 400763223, 1129760415, 3192727797,
        ];

        for p in primes {
            let P = LaurentPoly::from_string("x + 1 + x^-1", p);
            let Q = LaurentPoly::from_string("1 - x^2", p);
            let dfao = DFAO::lin_rep_reverse_machine(&P, &Q, 10000).unwrap();
            for n in 0..mots.len() as u64 {
                assert_eq!(
                    dfao.compute_ct_reverse(n, &Q),
                    ModInt::new(mots[n as usize], p)
                );
            }
        }

        let P = LaurentPoly::from_string("x + 1 + x^-1", 11);
        let Q = LaurentPoly::from_string("1 - x^2", 11);
        let dfao = DFAO::lin_rep_reverse_machine(&P, &Q, 279);
        assert_eq!(dfao.unwrap_err(), "Number of states exceeded 279.");
    }

    #[test]
    fn test_compute_shortest_poly_prop_directly() {
        for (poly_str, p, first_zero) in [
            ("x + 1 + x^-1", 3, 2),
            ("x + 1 + x^-2", 5, 39),
            ("x + 1 + x^-7", 5, 14),
            ("x + 1 + x^-48", 5, 49),
        ]
        .iter()
        {
            assert_eq!(
                DFAO::compute_shortest_poly_prop_directly(
                    &LaurentPoly::from_string(poly_str, *p),
                    &LaurentPoly::one(*p),
                    |v| v.constant_term() == ModInt::zero(*p),
                    Arc::new(AtomicBool::new(false))
                )
                .unwrap()
                .unwrap(),
                *first_zero
            );
        }
    }

    #[test]
    fn test_compute_shortest_prop_using_msd_dfao() {
        for (poly_str, p, first_zero) in [
            ("x + 1 + x^-1", 3, 2),
            ("x + 1 + x^-2", 5, 39),
            ("x + 1 + x^-7", 5, 14),
            ("x + 1 + x^-48", 5, 49),
        ]
            .iter()
        {
            let lin_rep = LinRep::for_ct_sequence(
                &LaurentPoly::from_string(poly_str, *p),
                &LaurentPoly::one(*p),
            );
            assert_eq!(
                DFAO::compute_shortest_prop_using_msd_dfao(
                    &lin_rep,
                    |v| v.constant_term() == ModInt::zero(*p),
                    10000,
                    None
                )
                    .unwrap()
                    .unwrap(),
                *first_zero
            );
        }
        
        let lin_rep = LinRep::for_ct_sequence(
            &LaurentPoly::from_string("x + 1 + x^-11", 2),
            &LaurentPoly::one(2),
        );
        assert_eq!(
            DFAO::compute_shortest_prop_using_msd_dfao(
                &lin_rep,
                |v| v.constant_term() == ModInt::zero(2),
                10000,
                None
            )
            .unwrap(),
            None
        );
    }

    #[test]
    fn test_compute_shortest_zero() {
        for (poly_str, p, first_zero) in [
            ("x + 1 + x^-1", 3, 2),
            ("x + 1 + x^-2", 5, 39),
            ("x + 1 + x^-7", 5, 14),
            ("x + 1 + x^-48", 5, 49),
        ]
            .iter()
        {
            assert_eq!(
                DFAO::compute_shortest_zero(
                    &LaurentPoly::from_string(poly_str, *p),
                    &LaurentPoly::one(*p),
                    10000
                )
                    .unwrap()
                    .unwrap(),
                *first_zero
            );
        }
        
        assert_eq!(
            DFAO::compute_shortest_zero(
                &LaurentPoly::from_string("x + 1 + x^-1001", 2),
                &LaurentPoly::one(2),
                10000
            )
            .unwrap(),
            None
        );
    }
}
