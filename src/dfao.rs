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
        // Membership index (state -> its position in `states_with_paths`) so the dedup is an
        // O(1) hash lookup instead of an O(N) linear scan; mirrors `build_kernel`'s `index`
        // and `forward_reachable`'s `seen`. The per-state `path` stays in `states_with_paths`.
        let mut state_index: HashMap<S, usize> = HashMap::new();
        state_index.insert(initial_state.clone(), 0);
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
                let new_state_index = match state_index.get(&new_state) {
                    Some(&j) => j,
                    None => {
                        let mut new_path = path.clone();
                        new_path.push(i);

                        if stop_prop(&new_state) {
                            return Ok(Right(new_path));
                        }
                        let j = states_with_paths.len();
                        state_index.insert(new_state.clone(), j);
                        states_with_paths.push((new_state, new_path));
                        if states_with_paths.len() > state_bound {
                            return Err(format!("Number of states exceeded {}.", state_bound));
                        }
                        j
                    }
                };

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
        // P^0..P^{p-1} do not depend on the state; compute them once and index by digit,
        // as `LinRep::for_ct_sequence` does with its transition matrices.
        let powers = Self::poly_powers(P, p);
        DFAO::from_reduction_rules(
            &Q,
            p,
            |state, i| powers[i.value as usize].mul(state).lambda_reduce(),
            state_bound,
            None,
        )
    }

    /// The **raw Rowland–Zeilberger automaton** of `ct[P^n Q] mod p`, with **no minimization
    /// applied**.
    ///
    /// "No minimization" means specifically:
    /// - the Moore partition-refinement minimizer (`DFAO::minimize`) is **not** run, and
    /// - no kernel / minimal-dual reduction (`minimization::minimal_dual` etc.) is applied.
    ///
    /// The states are the distinct reduced Laurent polynomials reachable from `Q` under the RZ
    /// transition `u |-> Lambda_p(P^d * u)` (digit `d`), and its output on a state is that
    /// polynomial's constant term. The BFS construction does identify two reachable states when
    /// they are the *same reduced polynomial* — this construction-time merge is intrinsic to the
    /// RZ machine itself, **not** a minimization step; the returned machine can still have more
    /// states than the Myhill–Nerode minimal one.
    ///
    /// This is the explicit, supported entry point for callers that need the un-minimized RZ
    /// machine (e.g. to compare it against its own minimization). It is a thin alias of
    /// `poly_auto`, which already builds exactly this machine. `state_bound` caps the BFS;
    /// `Err` is returned if the raw state count exceeds it.
    ///
    /// To obtain the minimized machine, call `.minimize(p, |s| s.constant_term())` on the result
    /// (that pair — raw then minimize — is equivalent to the standard minimized path).
    pub fn rz_machine(
        P: &LaurentPoly,
        Q: &LaurentPoly,
        state_bound: usize,
    ) -> Result<Self, String> {
        Self::poly_auto(P, Q, state_bound)
    }

    /// `[P^0, P^1, ..., P^{p-1}]` built incrementally (one multiply per power).
    fn poly_powers(P: &LaurentPoly, p: u64) -> Vec<LaurentPoly> {
        let mut powers = Vec::with_capacity(p as usize);
        let mut acc = LaurentPoly::one(P.modulus);
        for _ in 0..p {
            powers.push(acc.clone());
            acc = acc.mul(P);
        }
        powers
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
        // Hoist the per-digit powers P^0..P^{p-1} out of the per-state reduction rule.
        let powers = Self::poly_powers(P, P.modulus);
        Self::from_reduction_rules_until_prop(
            &Q,
            P.modulus,
            |state, i| powers[i.value as usize].mul(state).lambda_reduce(),
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
                    thread::sleep(Duration::from_millis(100));
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
                   thread::sleep(Duration::from_millis(100));
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
    /// Returns the minimal DFAO equivalent to `self` via Moore partition refinement
    /// (the Myhill-Nerode quotient), using `output_func` to label states. The machine is
    /// assumed to contain only reachable states; `from_reduction_rules` builds exactly that
    /// (BFS over the reachable state space), so the constructors above satisfy the contract.
    ///
    /// The result keeps the same `DFAO<ModInt, S>` shape (one representative state per block,
    /// `states[0]` still the initial state) so it can be serialized/composed like any other
    /// machine. Two states are merged iff they have equal output and, recursively, equal
    /// behavior on every digit (`O: Eq + Hash` is the output type that drives the partition).
    pub fn minimize<F, O>(&self, modulus: u64, output_func: F) -> DFAO<ModInt, S>
    where
        F: Fn(&S) -> O,
        O: Eq + Hash,
    {
        let n = self.states.len();
        let p = modulus as usize;

        // Flatten to next[s][digit] = target state index.
        let index_of: HashMap<&S, usize> =
            self.states.iter().enumerate().map(|(i, s)| (s, i)).collect();
        let mut next = vec![vec![0usize; p]; n];
        for (s, state) in self.states.iter().enumerate() {
            for d in 0..p {
                let to = self
                    .transitions
                    .get(&(state.clone(), ModInt::new(d as u64, modulus)))
                    .expect("transition missing");
                next[s][d] = *index_of.get(to).expect("transition target not in states");
            }
        }

        // Initial partition: group states by output value.
        let mut out_label: HashMap<O, u32> = HashMap::new();
        let mut block = vec![0u32; n];
        for s in 0..n {
            let label = out_label.len() as u32;
            block[s] = *out_label.entry(output_func(&self.states[s])).or_insert(label);
        }

        // Refine on (own block, blocks of successors) until the block count stabilizes.
        loop {
            let mut sig: HashMap<Vec<u32>, u32> = HashMap::new();
            let mut new_block = vec![0u32; n];
            for s in 0..n {
                let mut key = Vec::with_capacity(1 + p);
                key.push(block[s]);
                for d in 0..p {
                    key.push(block[next[s][d]]);
                }
                let label = sig.len() as u32;
                new_block[s] = *sig.entry(key).or_insert(label);
            }
            let old_count = block.iter().collect::<std::collections::HashSet<_>>().len();
            if sig.len() == old_count {
                break;
            }
            block = new_block;
        }

        // Build the quotient: the lowest-indexed original state of each block is its
        // representative, blocks are re-indexed in order of first appearance (so the block
        // containing the original initial state stays at index 0).
        let mut rep: Vec<usize> = Vec::new(); // new index -> representative original state
        let mut new_index: HashMap<u32, usize> = HashMap::new();
        for s in 0..n {
            new_index.entry(block[s]).or_insert_with(|| {
                rep.push(s);
                rep.len() - 1
            });
        }

        let new_states: Vec<S> = rep.iter().map(|&s| self.states[s].clone()).collect();
        let mut new_transitions: HashMap<(S, ModInt), S> = HashMap::new();
        for (new_i, &s) in rep.iter().enumerate() {
            for d in 0..p {
                let target_new = new_index[&block[next[s][d]]];
                new_transitions.insert(
                    (new_states[new_i].clone(), ModInt::new(d as u64, modulus)),
                    new_states[target_new].clone(),
                );
            }
        }

        DFAO {
            states: new_states,
            transitions: new_transitions,
        }
    }

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
    fn test_rz_machine_raw_reproduces_sequence() {
        // The raw (un-minimized) RZ machine must compute the exact ct[P^n Q] values, and must be
        // identical to poly_auto (it is a documented alias for the raw RZ construction).
        let primes = vec![2, 3, 5, 7, 11, 13];
        let mots: Vec<u64> = vec![
            1, 1, 2, 4, 9, 21, 51, 127, 323, 835, 2188, 5798, 15511, 41835, 113634, 310572, 853467,
            2356779, 6536382, 18199284, 50852019, 142547559, 400763223, 1129760415, 3192727797,
        ];

        for p in primes {
            let P = LaurentPoly::from_string("x + 1 + x^-1", p);
            let Q = LaurentPoly::from_string("1 - x^2", p);
            let raw = DFAO::rz_machine(&P, &Q, 10000).unwrap();
            for n in 0..mots.len() as u64 {
                assert_eq!(raw.compute_ct(n), ModInt::new(mots[n as usize], p));
            }
            // rz_machine is the raw RZ construction == poly_auto (same states & transitions).
            let pa = DFAO::poly_auto(&P, &Q, 10000).unwrap();
            assert_eq!(
                raw.serialize(p, |s| s.constant_term()),
                pa.serialize(p, |s| s.constant_term()),
                "rz_machine must equal poly_auto at p={p}"
            );
        }
    }

    #[test]
    fn test_rz_machine_minimize_matches_standard() {
        // raw + minimize == standard: minimizing the raw RZ machine must yield the same state
        // count as the standard minimized lsd path (lin_rep_machine.minimize, which the census
        // anchors at p=5 -> 10). lin_rep_machine and poly_auto are the same construction
        // (serialization-identical, per test_lin_rep_machine), so the minimized counts agree.
        let expected: [(u64, usize); 3] = [(3, 6), (5, 10), (7, 10)];
        for (p, exp_min) in expected {
            let P = LaurentPoly::from_string("x + 1 + x^-1", p);
            let Q = LaurentPoly::from_string("1 - x^2", p);

            let raw = DFAO::rz_machine(&P, &Q, 100000).unwrap();
            let raw_min = raw.minimize(p, |s: &LaurentPoly| s.constant_term());
            assert_eq!(
                raw_min.states.len(),
                exp_min,
                "raw+minimize count at p={p}"
            );

            // The standard minimized path (lin_rep forward machine, then minimize) must agree.
            let std = DFAO::lin_rep_machine(&P, &Q, 100000).unwrap();
            let std_min = std.minimize(p, |s: &ModIntVector| s.constant_term());
            assert_eq!(
                raw_min.states.len(),
                std_min.states.len(),
                "raw+minimize != standard minimized count at p={p}"
            );

            // And the minimized raw machine still computes the sequence.
            for n in 0..100u64 {
                assert_eq!(raw_min.compute_ct(n), raw.compute_ct(n), "raw_min value at n={n}, p={p}");
            }
        }
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
    fn test_minimize_parity_and_equivalence() {
        // Motzkin ct[P^n (1-x^2)] mod p: minimized lsd / msd state counts must match the
        // census known-answer anchor (p=5 -> 10 lsd / 40 msd). The lin_rep machines are the
        // same construction the experiment driver minimizes, so this is a direct parity check.
        let expected: [(u64, usize, usize); 3] = [(3, 6, 6), (5, 10, 40), (7, 10, 97)];
        for (p, exp_lsd, exp_msd) in expected {
            let pp = LaurentPoly::from_string("x + 1 + x^-1", p);
            let qq = LaurentPoly::from_string("1 - x^2", p);

            let fwd = DFAO::lin_rep_machine(&pp, &qq, 100000).unwrap();
            let fmin = fwd.minimize(p, |s: &ModIntVector| s.constant_term());
            assert_eq!(fmin.states.len(), exp_lsd, "lsd min count at p={p}");

            let rev = DFAO::lin_rep_reverse_machine(&pp, &qq, 100000).unwrap();
            let qv = ModIntVector::from_poly(&qq, rev.states[0].dim);
            let rmin = rev.minimize(p, |s: &ModIntVector| s.dot(&qv));
            assert_eq!(rmin.states.len(), exp_msd, "msd min count at p={p}");

            // Equivalence: minimized lsd machine computes the same sequence as the original.
            for n in 0..200u64 {
                assert_eq!(
                    fmin.compute_ct(n),
                    fwd.compute_ct(n),
                    "minimized lsd value mismatch at n={n}, p={p}"
                );
            }

            // Idempotence: re-minimizing a minimal machine is a no-op (same count).
            let fmin2 = fmin.minimize(p, |s: &ModIntVector| s.constant_term());
            assert_eq!(fmin2.states.len(), fmin.states.len());
        }
    }

    #[test]
    fn test_minimize_collapses_duplicates() {
        // Build a DFAO over ModInt states (output = the state value mod 2) with deliberately
        // Nerode-equivalent states, and check the minimizer merges exactly them.
        let p = 2u64;
        // states 0 and 1 are equivalent: same output (both 0), and each sends every digit into
        // {0,1}; state 2 has distinct output 1 and self-loops. Expect 3 -> 2. (State labels live
        // under modulus 3 so the value-2 label stays distinct from value 0; the alphabet is p=2.)
        let s = |v: u64| ModInt::new(v, 3);
        let states = vec![s(0), s(1), s(2)];
        let mut transitions: HashMap<(ModInt, ModInt), ModInt> = HashMap::new();
        for d in 0..p {
            transitions.insert((s(0), ModInt::new(d, p)), s(1));
            transitions.insert((s(1), ModInt::new(d, p)), s(0));
            transitions.insert((s(2), ModInt::new(d, p)), s(2));
        }
        let machine = DFAO { states, transitions };
        // output: state 2 -> 1, states 0/1 -> 0
        let mm = machine.minimize(p, |st: &ModInt| if st.value == 2 { 1u64 } else { 0u64 });
        assert_eq!(mm.states.len(), 2, "equivalent states should collapse to 2");

        // A machine whose states differ by output must NOT collapse.
        let states2 = vec![s(0), s(1)];
        let mut transitions2: HashMap<(ModInt, ModInt), ModInt> = HashMap::new();
        for d in 0..p {
            transitions2.insert((s(0), ModInt::new(d, p)), s(0));
            transitions2.insert((s(1), ModInt::new(d, p)), s(1));
        }
        let machine2 = DFAO { states: states2, transitions: transitions2 };
        let mm2 = machine2.minimize(p, |st: &ModInt| st.value);
        assert_eq!(mm2.states.len(), 2, "distinct-output states must not merge");
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
