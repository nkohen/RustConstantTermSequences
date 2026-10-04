//! mod-p^k constant-term automaton via the **exact value-signature p-kernel** build.
//!
//! This is a deliberately *separate* code path from the GF(p) builders in `dfao.rs`
//! (`poly_auto`, `lin_rep_machine`): those rely on the Frobenius identity
//! `P^p ≡ P(t^p)` baked into `LaurentPoly::lambda_reduce`, which is valid over the
//! field GF(p) but **false over Z/p^k Z** (`p^k` is not prime for `k ≥ 2`, so the
//! ring is not a field and the Frobenius endomorphism collapses). There is therefore
//! no reduction-rule on polynomial/vector states; instead the automaton is built from
//! the **exact integer values** `a(n) = ct[P^n Q] mod p^k` and its states are keyed by
//! a *value signature* sampled from those values.
//!
//! Pipeline (correct-by-construction; no Frobenius assumption):
//!   1. `mod_pk_values` — exact `a(n)` for `n in 0..=cap` by `w <- P·w`, `a(n)=[t^0]w`.
//!   2. `build_kernel`  — BFS the p-kernel `{ n ↦ a(p^j·n + i) }`, dedupe by value
//!      signature; output of a state is the exact `a(i)`. The automaton reads `n`
//!      **LSD-first** (δ((j,i),d) = (j+1, i + d·p^j)).
//!   3. `DFAO::minimize` — the *shared* Moore minimizer from `dfao.rs` (brief 04); this
//!      module does **not** carry its own minimizer (see reconciliation note below).
//!
//! ## Gotchas (do not ignore)
//! - **Validity window.** The automaton is valid only up to index `cap`. A child state
//!   whose offset exceeds `cap`, or whose signature cannot be sampled to `MIN_SAMPLES`,
//!   sets `incomplete = true`. Size `cap` so that `cap >> p·win` or the signature is
//!   under-sampled and the build is unsound. A caller MUST check `incomplete` (and is
//!   strongly advised to `validate` against the exact values) before trusting the result.
//! - **Outputs are exact values mod p^k (`u64`), not GF(p) residues.** Do not coerce
//!   them to a GF(p) residue unless you are deliberately extracting a residue class.
//! - **One minimizer.** The originating experiment carried a bespoke Moore pass that
//!   returned `(count, block_map)` so it could rebuild the quotient by hand. That is
//!   redundant here: `DFAO::minimize` (brief 04) already performs the partition
//!   refinement *and* returns the rebuilt minimal `DFAO`, so the block map never needs
//!   to escape. We reuse it unchanged with `output_func = |s| s.out`.

use crate::laurent_poly::LaurentPoly;
use crate::mod_int::ModInt;
use crate::dfao::DFAO;
use std::collections::HashMap;

/// Minimum number of value samples a kernel-state signature must contain before two
/// states are considered distinguishable. A signature with fewer samples is treated as
/// under-sampled and flags the build `incomplete`.
pub const MIN_SAMPLES: usize = 24;

/// A state of the kernel automaton. Identity is the value signature (carried implicitly
/// by the BFS dedup); the fields below are what we keep per surviving state. `out` is the
/// exact value `a(off) mod p^k` emitted at this state. `stride`/`off` are retained for
/// debugging / provenance (the kernel rep is `n ↦ a(stride·n + off)`).
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct KernelState {
    pub stride: u64, // p^j
    pub off: u64,    // i
    pub out: u64,    // exact a(off) mod p^k
}

/// The result of a mod-p^k kernel build, before minimization.
pub struct KernelAuto {
    /// Number of (reachable) kernel states.
    pub n: usize,
    /// `trans[s][d]` = next state index on digit `d ∈ [0, p)` (LSD-first).
    pub trans: Vec<Vec<u32>>,
    /// `out[s]` = exact `a(i)` for the kernel rep of state `s` (NOT reduced to GF(p)).
    pub out: Vec<u64>,
    /// Prime base `p`.
    pub p: usize,
    /// `true` if the value `cap` was exceeded / a signature was under-sampled while
    /// building — in which case the automaton is only a candidate and must not be
    /// trusted without enlarging `cap`.
    pub incomplete: bool,
}

/// Exact value sweep: `a(n) = ct[P^n Q] mod m`, for `n = 0..=cap`, where `m = P.modulus`.
///
/// Incremental: `w <- P·w` with exact integer multiplication reduced mod `m` each step
/// (this is `LaurentPoly::mul`, which carries the modulus). **No Frobenius / lambda
/// reduction** is used — that is the whole point of the mod-p^k path.
pub fn mod_pk_values(p_poly: &LaurentPoly, q_poly: &LaurentPoly, cap: usize) -> Vec<u64> {
    assert_eq!(p_poly.modulus, q_poly.modulus, "P and Q moduli must match");
    let mut a = Vec::with_capacity(cap + 1);
    let mut w = q_poly.clone();
    for _ in 0..=cap {
        a.push(w.constant_term().value);
        w = w.mul(p_poly);
    }
    a
}

/// Build the p-kernel automaton from an exact value table `a` (as produced by
/// [`mod_pk_values`]). `p` is the prime base; `win` is the signature window length.
///
/// A kernel rep `(stride = p^j, off = i)` is identified by the signature
/// `[ a(stride·n + off) : n = 0..win ]` (truncated where indices exceed `cap`). Two reps
/// with equal signature are the same state. The transition `δ((j,i), d) = (j+1, i+d·p^j)`
/// reads `n` LSD-first. If a child offset exceeds `cap`, or a signature cannot reach
/// `MIN_SAMPLES`, `incomplete` is set.
pub fn build_kernel(a: &[u64], p: usize, win: usize) -> KernelAuto {
    let cap = a.len() - 1;
    // signature of kernel rep (stride, off): sample a(stride·n + off) for as many n as fit.
    let sig = |stride: u64, off: u64| -> Option<Vec<u64>> {
        let mut v = Vec::with_capacity(win);
        for nn in 0..win as u64 {
            let idx = stride.checked_mul(nn).and_then(|x| x.checked_add(off));
            match idx {
                Some(ix) if (ix as usize) <= cap => v.push(a[ix as usize]),
                _ => break,
            }
        }
        if v.len() >= MIN_SAMPLES.min(win) {
            Some(v)
        } else {
            None
        }
    };

    let mut index: HashMap<Vec<u64>, u32> = HashMap::new();
    let mut states: Vec<(u64, u64)> = Vec::new(); // (stride = p^j, off = i)
    let mut trans: Vec<Vec<u32>> = Vec::new();
    let mut out: Vec<u64> = Vec::new();
    let mut incomplete = false;

    let s0 = sig(1, 0).expect("cap too small for the start signature");
    index.insert(s0, 0);
    states.push((1, 0));
    out.push(a[0]);

    let mut k = 0usize;
    while k < states.len() {
        let (stride, off) = states[k];
        let mut row = vec![0u32; p];
        for d in 0..p as u64 {
            let cstride = stride.saturating_mul(p as u64);
            let coff = off + d * stride;
            let id = match sig(cstride, coff) {
                Some(s) => {
                    if let Some(&i) = index.get(&s) {
                        i
                    } else {
                        let i = states.len() as u32;
                        index.insert(s, i);
                        states.push((cstride, coff));
                        out.push(if (coff as usize) <= cap {
                            a[coff as usize]
                        } else {
                            incomplete = true;
                            0
                        });
                        i
                    }
                }
                None => {
                    // Under-sampled child: flag and self-loop-ish fallback to keep the
                    // automaton total. The `incomplete` flag is the caller's signal that
                    // this state is not trustworthy.
                    incomplete = true;
                    k as u32
                }
            };
            row[d as usize] = id;
        }
        trans.push(row);
        k += 1;
    }

    KernelAuto {
        n: states.len(),
        trans,
        out,
        p,
        incomplete,
    }
}

/// Convert a freshly-built [`KernelAuto`] into the generic [`DFAO`] over `ModInt` digits
/// and [`KernelState`] states, so it composes with the rest of the library (serialize,
/// minimize, …). State outputs are the exact values via [`KernelState::out`].
fn kernel_to_dfao(auto: &KernelAuto) -> DFAO<ModInt, KernelState> {
    let p = auto.p as u64;
    let states: Vec<KernelState> = (0..auto.n)
        .map(|s| KernelState {
            // stride/off provenance is not retained past the build; only `out` drives
            // behavior, so we use the index as a stable, distinct placeholder for the
            // first two fields (they do not participate in minimization — `out` does).
            stride: s as u64,
            off: s as u64,
            out: auto.out[s],
        })
        .collect();

    let mut transitions: HashMap<(KernelState, ModInt), KernelState> = HashMap::new();
    for s in 0..auto.n {
        for d in 0..auto.p {
            let to = auto.trans[s][d] as usize;
            transitions.insert(
                (states[s].clone(), ModInt::new(d as u64, p)),
                states[to].clone(),
            );
        }
    }

    DFAO {
        states,
        transitions,
    }
}

/// The full mod-p^k pipeline: exact values → p-kernel automaton → shared Moore minimizer.
///
/// Returns the **minimized** [`DFAO`] (states keyed by exact value, LSD-first), together
/// with the `incomplete` flag. Build params:
/// - `p_poly`, `q_poly` — `P`, `Q` with modulus `p^k` (use [`LaurentPoly::from_vec`] /
///   [`LaurentPoly::from_string`] over modulus `p^k`).
/// - `p` — the prime base (the alphabet size; `p_poly.modulus` must equal `p^k`).
/// - `k` — the exponent; checked against `p_poly.modulus == p^k`.
/// - `cap` — exact-value window. **Must be `>> p·win`** (see module gotchas); the build
///   sets `incomplete` if it under-samples.
/// - `win` — signature window length for kernel-state identity (e.g. 96).
///
/// `Err` on a modulus mismatch. A non-`incomplete` result still warrants a [`validate`]
/// pass before being trusted at scale.
pub fn mod_pk_kernel(
    p_poly: &LaurentPoly,
    q_poly: &LaurentPoly,
    p: u64,
    k: u32,
    cap: usize,
    win: usize,
) -> Result<(DFAO<ModInt, KernelState>, bool), String> {
    let m = p.checked_pow(k).ok_or("p^k overflows u64")?;
    if p_poly.modulus != m || q_poly.modulus != m {
        return Err(format!(
            "modulus mismatch: P/Q must be built over p^k = {m}, got P={} Q={}",
            p_poly.modulus, q_poly.modulus
        ));
    }
    let a = mod_pk_values(p_poly, q_poly, cap);
    let auto = build_kernel(&a, p as usize, win);
    let incomplete = auto.incomplete;
    let dfao = kernel_to_dfao(&auto);
    let minimized = dfao.minimize(p, |s: &KernelState| s.out);
    Ok((minimized, incomplete))
}

/// Evaluate a (possibly minimized) kernel DFAO on `n`, LSD-first, returning the exact
/// output value. Mirrors the kernel transition order `δ(s,d)` reading the least
/// significant digit first.
pub fn dfao_value(auto: &DFAO<ModInt, KernelState>, n: u64, p: u64) -> u64 {
    auto.compute_lsd(n, p, |s: &KernelState| s.out)
}

/// Replay validation: confirm the DFAO reproduces the exact `a(n)` on `0..=cap2`.
/// Returns `Err((n, got, want))` at the first disagreement, `Ok(())` if all agree.
pub fn validate(
    auto: &DFAO<ModInt, KernelState>,
    a: &[u64],
    p: u64,
    cap2: usize,
) -> Result<(), (usize, u64, u64)> {
    let lim = cap2.min(a.len() - 1);
    for n in 0..=lim {
        let got = dfao_value(auto, n as u64, p);
        if got != a[n] {
            return Err((n, got, a[n]));
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn poly(terms: &[(i64, i64)], m: u64) -> LaurentPoly {
        let v: Vec<(i64, u64)> = terms
            .iter()
            .map(|&(e, c)| (e, ((c % m as i64 + m as i64) % m as i64) as u64))
            .collect();
        LaurentPoly::from_vec(v, m)
    }

    // Build (P over modulus p^k, Q over modulus p^k) from term lists at a chosen (p,k).
    fn build(
        p_terms: &[(i64, i64)],
        q_terms: &[(i64, i64)],
        p: u64,
        k: u32,
        cap: usize,
        win: usize,
    ) -> (usize, bool, bool) {
        let m = p.pow(k);
        let pp = poly(p_terms, m);
        let qq = poly(q_terms, m);
        let a = mod_pk_values(&pp, &qq, cap);
        let (dfao, incomplete) = mod_pk_kernel(&pp, &qq, p, k, cap, win).unwrap();
        let valid = validate(&dfao, &a, p, cap).is_ok();
        (dfao.states.len(), incomplete, valid)
    }

    #[test]
    fn modulus_mismatch_is_an_error() {
        // P built over the prime, not p^k = 4 -> must be rejected.
        let pp = poly(&[(-1, 1), (1, 1)], 2);
        let qq = poly(&[(0, 1)], 2);
        let err = mod_pk_kernel(&pp, &qq, 2, 2, 1024, 96).unwrap_err();
        assert!(err.contains("modulus mismatch"), "got: {err}");
    }

    #[test]
    fn anchor_powers_of_two_mod_4() {
        // P = x + x^-1, Q = 1, mod 4. Proven: a(2m) = C(2m,m); C(2m,m) ≡ 2 (mod 4)
        // iff m is a power of two, so S_2 = {2^{j+1}} (powers of two). The minimal
        // mod-4 automaton has 4 states (recorded experiment value).
        let cap = 1 << 14;
        let win = 96;
        let (nmin, incomplete, valid) =
            build(&[(-1, 1), (1, 1)], &[(0, 1)], 2, 2, cap, win);
        assert_eq!(nmin, 4, "anchor mod 4 minimized states");
        assert!(!incomplete, "anchor kernel must close within cap");
        assert!(valid, "anchor DFAO must reproduce exact a(n)");

        // S_2 == powers of two, read directly off the exact values.
        let m = 4u64;
        let a = mod_pk_values(&poly(&[(-1, 1), (1, 1)], m), &poly(&[(0, 1)], m), cap);
        let s2: Vec<usize> = (0..=cap).filter(|&n| a[n] == 2).collect();
        let pw: Vec<usize> = (1..).map(|j| 1usize << j).take_while(|&x| x <= cap).collect();
        assert_eq!(s2, pw, "S_2 must be exactly the powers of two");
    }

    #[test]
    fn reproduces_experiment_state_counts() {
        // Recorded minimized state counts from experiments/mod-pk-sparseness
        // (output.txt / output-sweep.txt). The library build must match exactly on the
        // closed (incomplete == false, valid == true) cases.
        // (p, k, P-terms, Q-terms, expected MIN, cap). cap matches the experiment's
        // window (2^14) where the kernel needs it to close; smaller caps suffice for the
        // small-p cases and keep the test cheap. The per-case `!incomplete && valid` guard
        // catches any under-sampling. Each cap is well above p·win.
        let win = 96;
        type Case = (u64, u32, Vec<(i64, i64)>, Vec<(i64, i64)>, usize, usize);
        let cases: Vec<Case> = vec![
            // mod 4 (p=2, k=2)
            (2, 2, vec![(-1, 1), (1, 1)], vec![(0, 1)], 4, 1 << 12),
            (2, 2, vec![(-1, 1), (0, 1), (1, 1)], vec![(0, 1)], 3, 1 << 12),
            (2, 2, vec![(-2, 1), (0, 1), (2, 1)], vec![(0, 1)], 3, 1 << 12),
            (2, 2, vec![(-2, 1), (-1, 1), (0, 1), (1, 1), (2, 1)], vec![(0, 1)], 4, 1 << 12),
            // mod 8 (p=2, k=3)
            (2, 3, vec![(-1, 1), (1, 1)], vec![(0, 1)], 7, 1 << 12),
            (2, 3, vec![(-1, 1), (0, 1), (1, 1)], vec![(0, 1)], 9, 1 << 12),
            (2, 3, vec![(-2, 1), (0, 1), (2, 1)], vec![(0, 1)], 9, 1 << 12),
            // p=3, k=1
            (3, 1, vec![(-1, 1), (1, 1)], vec![(0, 1)], 3, 1 << 12),
            (3, 1, vec![(-1, 1), (1, 1), (0, 1)], vec![(0, 1)], 2, 1 << 12),
            (3, 1, vec![(-2, 1), (0, 1), (2, 1)], vec![(0, 1)], 2, 1 << 12),
            // p=5, k=1 (the closed/valid cases from the sweep; need the full 2^14 window).
            (5, 1, vec![(-1, 1), (0, 1), (1, 1)], vec![(0, 1)], 4, 1 << 14),
            (5, 1, vec![(-1, 1), (0, 4), (1, 1)], vec![(0, 1)], 4, 1 << 14),
            (5, 1, vec![(-2, 1), (0, 1), (2, 1)], vec![(0, 1)], 4, 1 << 14),
        ];
        for (p, k, pt, qt, exp, cap) in cases {
            let (nmin, incomplete, valid) = build(&pt, &qt, p, k, cap, win);
            assert!(
                !incomplete && valid,
                "case p={p} k={k} P={pt:?} expected to close & validate (inc={incomplete} valid={valid})"
            );
            assert_eq!(
                nmin, exp,
                "minimized state count mismatch for p={p} k={k} P={pt:?}"
            );
        }
    }

    #[test]
    fn k1_agrees_with_gf_p_poly_auto() {
        // At k=1 the mod-p^k build is over the field GF(p), so it must agree with the
        // GF(p) poly_auto minimal automaton for the same sequence: same minimized state
        // count AND same value sequence.
        let cap = 1 << 13;
        let win = 96;
        let mut closed = 0usize;
        for p in [3u64, 5, 7] {
            // Use a recurrent (non-degenerate) P so the GF(p) automaton is interesting:
            // P = x + 1 + x^-1, Q = 1.
            let pt = vec![(-1i64, 1i64), (0, 1), (1, 1)];
            let qt = vec![(0i64, 1i64)];

            // mod-p^k path at k=1.
            let pp = poly(&pt, p);
            let qq = poly(&qt, p);
            let (mk, incomplete) = mod_pk_kernel(&pp, &qq, p, 1, cap, win).unwrap();
            // Closure is a `cap`-vs-`p·win` artifact (larger p needs a larger cap to
            // sample the signature); only assert agreement on closed kernels, and require
            // a coverage floor so the cross-check is not vacuous.
            if incomplete {
                continue;
            }
            closed += 1;

            // GF(p) path: poly_auto, minimized with the shared minimizer.
            let gf = DFAO::poly_auto(&pp, &qq, 100000).unwrap();
            let gf_min = gf.minimize(p, |s: &LaurentPoly| s.constant_term());

            assert_eq!(
                mk.states.len(),
                gf_min.states.len(),
                "k=1 mod-p^k vs GF(p) minimal state count must agree for p={p}"
            );

            // Same sequence of values on a prefix.
            for n in 0..300u64 {
                let v_mk = dfao_value(&mk, n, p);
                let v_gf = gf_min.compute_lsd(n, p, |s: &LaurentPoly| s.constant_term().value);
                assert_eq!(v_mk, v_gf, "k=1 value mismatch at n={n}, p={p}");
            }
        }
        assert!(
            closed >= 2,
            "k=1 cross-check vacuous: only {closed} kernels closed at cap={cap}"
        );
    }

    #[test]
    fn incomplete_window_is_surfaced() {
        // An aggressively small cap (relative to p·win) under-samples the signature, so
        // the build must report `incomplete = true` rather than return a silently-wrong
        // automaton. p=3, k=2 at the recorded cap closes; at a tiny cap it must not.
        let win = 96;
        let p = 3u64;
        let k = 2u32;
        let m = p.pow(k);
        let pp = poly(&[(-1, 1), (1, 1)], m);
        let qq = poly(&[(0, 1)], m);
        let (_dfao, incomplete) = mod_pk_kernel(&pp, &qq, p, k, 200, win).unwrap();
        assert!(incomplete, "tiny cap must surface incomplete=true");
    }

    #[test]
    fn mod_pk_values_matches_central_binomial() {
        // Independent closed-form anchor for the exact value sweep (no automaton involved):
        // for P = x + x^-1 and Q = 1, ct[P^n] = C(n, n/2), i.e. a(2m) = C(2m, m) and a(odd) = 0.
        // C(2m,m) is accumulated in u128 so it never overflows for the n-range used; this checks
        // mod_pk_values across several (p, k), including the non-field rings Z/p^k Z (k >= 2)
        // whose whole point is that no Frobenius reduction is applied.
        fn central_binomial(m: u64, modulus: u64) -> u64 {
            // C(2m, m) via the exact integer recurrence C(2m,i+1) = C(2m,i)*(2m-i)/(i+1);
            // each partial product is itself a binomial coefficient, so the division is exact.
            let mut c: u128 = 1;
            for i in 0..m {
                c = c * ((2 * m - i) as u128) / ((i + 1) as u128);
            }
            (c % modulus as u128) as u64
        }
        let cap = 80usize; // 2m <= 80 keeps C(2m,m) well within u128
        let cases: &[(u64, u32)] = &[(2, 2), (2, 3), (3, 2), (5, 1)];
        for &(p, k) in cases {
            let m = p.pow(k);
            let a = mod_pk_values(&poly(&[(-1, 1), (1, 1)], m), &poly(&[(0, 1)], m), cap);
            for n in 0..=cap {
                let want = if n % 2 == 0 {
                    central_binomial((n / 2) as u64, m)
                } else {
                    0
                };
                assert_eq!(a[n], want, "mod_pk_values mismatch at n={n}, p={p}, k={k}");
            }
            // Discriminating: the sweep is neither all-zero nor constant (a vacuous sweep of all
            // 1s or all 0s would pass a weaker check but fails this one).
            assert!(
                a.iter().any(|&v| v != 0 && v != a[0]),
                "sweep is suspiciously trivial at p={p}, k={k}"
            );
        }
    }
}
