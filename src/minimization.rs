//! Minimization of the Rowland-Zeilberger linear representation of ct[P^n Q] mod p.
//!
//! Port of the Sage `Minimization` module (branch `upstream/rz-minimization` of the
//! ConstantTermSequences library). Algorithm 4.2 of the RZ-minimization note: build the
//! minimal forward (lsd) linear representation by a single linear quotient by the
//! forward-invisible kernel
//!
//! ```text
//! K = { u : ct[P^n u] = 0  for all n }   (= W^perp, Thm 3.2),
//! ```
//!
//! WITHOUT a reverse / observability pass. The forward machine has row-vector states `u`
//! over the ambient space of dimension `2*max_deg + 1` (entry order = `ModIntVector::from_poly`:
//! coefficient of `t^(-max_deg)` first), with transition `u |-> u * M_d = Lambda_p(P^d * u)`
//! (`mat_func[d].right_mul(u)` from `lin_rep.rs`) and output `u |-> ct[u]`.
//!
//! `K` is gamma-invariant (`u * M_k` stays in `K` for `k in [0, p)`) and a gamma-invariant slab
//! `K_0 <= K` is available in closed form with NO reverse pass from the §7 structural sources
//! (`kernel_seeds`): exactness/Euler, symmetry (palindromic `P` only), and support/Newton-polytope
//! gaps. Seeding `K` with these and gamma-closing recovers (a subspace of) `K`; when `K_0 = K`
//! the seeded quotient is already minimal (Cor 6.1). With INCOMPLETE seeds (`K_0 ⊊ K`) the result
//! is a correct but possibly NON-minimal automaton -- this port preserves that semantics and does
//! not claim minimality it cannot guarantee.
//!
//! The §7.4 syzygy/Jacobian source is not implemented (redundant for univariate trinomials), as
//! in the Sage reference.
//!
//! All GF(p) rank / span / dimension work routes through the single Gaussian-elimination routine
//! already in `lin_rep.rs` (`LinRep::row_basis` / `LinRep::mat_rank`, made crate-visible here);
//! the only extra linear-algebra primitive is `solve_coords` (express a vector in a given basis),
//! which those two do not provide.

use crate::laurent_poly::LaurentPoly;
use crate::lin_rep::LinRep;
use crate::mod_int::ModInt;
use crate::mod_int_matrix::ModIntMatrix;
use std::collections::HashMap;

/// A minimized forward linear representation `(v', {M_d'}, w')` of dimension `dim`
/// (the Schuetzenberger quotient of `lin_rep(P, Q, p)` by the kernel `K`).
///
/// `dim` is the minimal forward state-space dimension (= the Hankel rank) precisely when
/// `K = W^perp`, i.e. when the seeding is complete (Cor 6.1); with incomplete seeds it is a
/// correct but possibly non-minimal automaton.
pub struct MinimizedRep {
    pub row_vec: Vec<u64>,           // v', length `dim`
    pub mat_func: Vec<Vec<Vec<u64>>>, // {M_d'}, p matrices, each `dim` x `dim`
    pub col_vec: Vec<u64>,           // w', length `dim`
    pub dim: usize,
    pub modulus: u64,
}

// ---- §7 closed-form kernel seeds (brief 02) ----------------------------------------------

/// True iff `P(t) = P(1/t)` as Laurent polynomials over GF(p) (the monomial symmetry
/// sigma: t -> 1/t). Mirrors `is_palindromic` in the Sage reference.
pub fn is_palindromic(poly: &LaurentPoly, p: u64) -> bool {
    let reflected = reflect(poly, p);
    poly == &reflected
}

/// P(1/t): reflect every exponent.
fn reflect(poly: &LaurentPoly, p: u64) -> LaurentPoly {
    let mut terms = Vec::new();
    for e in lowest_exp(poly)..=highest_exp(poly) {
        let c = poly.get_coefficient(&e);
        if c.value != 0 {
            terms.push((-e, c.value));
        }
    }
    LaurentPoly::from_vec(terms, p)
}

// Lowest / highest exponent present (treating the zero poly as the single exponent 0).
fn highest_exp(poly: &LaurentPoly) -> i64 {
    let mut hi: Option<i64> = None;
    for e in -(scan_bound(poly))..=scan_bound(poly) {
        if poly.get_coefficient(&e).value != 0 {
            hi = Some(e);
        }
    }
    hi.unwrap_or(0)
}
fn lowest_exp(poly: &LaurentPoly) -> i64 {
    for e in -(scan_bound(poly))..=scan_bound(poly) {
        if poly.get_coefficient(&e).value != 0 {
            return e;
        }
    }
    0
}
// A safe scanning window: the magnitude of the largest exponent (LaurentPoly::degree is max |e|),
// or 0 for the zero polynomial (degree() would panic on an empty poly).
fn scan_bound(poly: &LaurentPoly) -> i64 {
    if poly == &LaurentPoly::zero(poly.modulus) {
        0
    } else {
        poly.degree() as i64
    }
}

/// `Lambda_p(poly)`: keep terms whose exponent is divisible by `p`, divide those exponents by `p`.
/// `LaurentPoly::lambda_reduce` already implements this.
fn lambda(poly: &LaurentPoly) -> LaurentPoly {
    poly.lambda_reduce()
}

/// The Euler operator `D P = t * dP/dt` (kills constant terms, so `ct[P^n D P] = 0`).
fn euler_derivative(poly: &LaurentPoly, p: u64) -> LaurentPoly {
    let mut terms = Vec::new();
    for e in lowest_exp(poly)..=highest_exp(poly) {
        let c = poly.get_coefficient(&e);
        if c.value != 0 {
            // t * d/dt (c t^e) = e c t^e
            let scaled = ModInt::from_i64(e, p) * c;
            if scaled.value != 0 {
                terms.push((e, scaled.value));
            }
        }
    }
    LaurentPoly::from_vec(terms, p)
}

/// Length-`(2*max_deg+1)` row of a Laurent polynomial (mass outside `[-max_deg, max_deg]` dropped),
/// in the `ModIntVector::from_poly` entry order. Returns `None` if the row is all-zero.
fn poly_to_row(poly: &LaurentPoly, max_deg: usize, p: u64) -> Option<Vec<u64>> {
    let dim = 2 * max_deg + 1;
    let mut row = vec![0u64; dim];
    let md = max_deg as i64;
    let mut nonzero = false;
    for k in 0..dim {
        let e = k as i64 - md;
        let c = poly.get_coefficient(&e).value % p;
        row[k] = c;
        if c != 0 {
            nonzero = true;
        }
    }
    if nonzero {
        Some(row)
    } else {
        None
    }
}

/// §7.1 exactness/Euler seeds for the FULL in-ambient module `{Lambda_p(P^s * DP)}`, `DP = t dP/dt`.
/// Seeds every `Lambda_p(P^s * DP)`, `s >= 0`, that lands inside `[-max_deg, max_deg]`; gamma-closure
/// recomposes these into the in-ambient shadow of the whole module. (A single `Lambda_p(DP)` seed is
/// strictly weaker -- `Lambda_p` does not commute with multiplication by `P^s`.)
pub fn seeds_exactness(poly: &LaurentPoly, p: u64, max_deg: usize) -> Vec<Vec<u64>> {
    let dp = euler_derivative(poly, p);
    let md = max_deg as i64;
    let s_max = 4 * max_deg + 4 * (p as usize) + 8; // window: enough to expose every in-ambient member
    let mut seeds = Vec::new();
    let mut ps = LaurentPoly::one(p);
    for _ in 0..=s_max {
        let g = lambda(&ps.mul(&dp));
        // include only if g already lives inside [-max_deg, max_deg]
        let mut inside = true;
        if g != LaurentPoly::zero(p) {
            for e in lowest_exp(&g)..=highest_exp(&g) {
                if g.get_coefficient(&e).value != 0 && (e < -md || e > md) {
                    inside = false;
                    break;
                }
            }
        }
        if inside {
            if let Some(row) = poly_to_row(&g, max_deg, p) {
                seeds.push(row);
            }
        }
        ps = ps.mul(poly);
    }
    seeds
}

/// §7.2 symmetry seeds `t^j - t^(-j)`, `j in [1, max_deg]` (the whole antisymmetric subspace lies
/// in `K` for the reflection sigma: t -> 1/t). Returns an empty vector unless `P` is palindromic.
pub fn seeds_symmetry(poly: &LaurentPoly, p: u64, max_deg: usize) -> Vec<Vec<u64>> {
    if !is_palindromic(poly, p) {
        return Vec::new();
    }
    let mut seeds = Vec::new();
    for j in 1..=(max_deg as i64) {
        let g = LaurentPoly::from_vec(vec![(j, 1), (-j, p - 1)], p);
        if let Some(row) = poly_to_row(&g, max_deg, p) {
            seeds.push(row);
        }
    }
    seeds
}

/// §7.3 support/Newton-polytope seeds: `t^j` with `-j` NOT in the numerical semigroup generated by
/// `supp(P)` (then `ct[P^n t^j] = 0` for all n). The semigroup is BFS-closed within the window
/// `[-2*max_deg, 2*max_deg]` so it cannot run away.
pub fn seeds_support(poly: &LaurentPoly, p: u64, max_deg: usize) -> Vec<Vec<u64>> {
    let md = max_deg as i64;
    // support exponents
    let mut supp: Vec<i64> = Vec::new();
    for e in lowest_exp(poly)..=highest_exp(poly) {
        if poly.get_coefficient(&e).value != 0 {
            supp.push(e);
        }
    }
    // BFS closure of supp under addition, restricted to [-2*max_deg, 2*max_deg].
    let mut semigroup: std::collections::HashSet<i64> = std::collections::HashSet::new();
    semigroup.insert(0);
    let mut frontier: Vec<i64> = vec![0];
    while !frontier.is_empty() {
        let mut next = Vec::new();
        for &s in &frontier {
            for &g in &supp {
                let u = s + g;
                if u >= -2 * md && u <= 2 * md && !semigroup.contains(&u) {
                    semigroup.insert(u);
                    next.push(u);
                }
            }
        }
        frontier = next;
    }
    let mut seeds = Vec::new();
    for j in -md..=md {
        if !semigroup.contains(&(-j)) {
            let g = LaurentPoly::from_vec(vec![(j, 1)], p);
            if let Some(row) = poly_to_row(&g, max_deg, p) {
                seeds.push(row);
            }
        }
    }
    seeds
}

/// Closed-form gamma-closure seeds for `K` (brief 02): the union of the applicable §7 sources.
/// Symmetry is auto-selected ONLY when `P` is palindromic; exactness and support are always
/// included. `max_deg = None` defaults to the natural ambient `LinRep::for_ct_sequence(P, 1)` half-width.
pub fn kernel_seeds(poly: &LaurentPoly, p: u64, max_deg: Option<usize>) -> Vec<Vec<u64>> {
    let max_deg = max_deg.unwrap_or_else(|| {
        let one = LaurentPoly::one(p);
        let rep = LinRep::for_ct_sequence(poly, &one);
        (rep.rank - 1) / 2
    });
    let mut seeds = seeds_exactness(poly, p, max_deg);
    seeds.extend(seeds_symmetry(poly, p, max_deg));
    seeds.extend(seeds_support(poly, p, max_deg));
    seeds
}

// ---- gamma-closure + Schuetzenberger quotient (brief 01) ---------------------------------

/// `u * M`: the row action (`M.right_mul(u)`), expressed on plain `Vec<u64>` rows. Uses the raw
/// `right_mul_u64` path so it skips the `ModInt`/`ModIntVector` round-trip.
fn row_times_mat(row: &[u64], mat: &ModIntMatrix, _p: u64) -> Vec<u64> {
    mat.right_mul_u64(row)
}

/// An incrementally-maintained reduced row basis over GF(p): each stored vector is normalized
/// (leading pivot = 1) and reduced against the earlier ones, so a new candidate's independence
/// is tested by one reduction pass. This lets `gamma_closure` reduce only the *new* rows each
/// round instead of re-basing the whole generator set.
struct RunningBasis {
    rows: Vec<Vec<u64>>,
    pivot_col: Vec<usize>,
    dim: usize,
    p: u64,
    inv_tab: Vec<u64>,
}

impl RunningBasis {
    fn new(dim: usize, p: u64) -> Self {
        RunningBasis {
            rows: Vec::new(),
            pivot_col: Vec::new(),
            dim,
            p,
            inv_tab: LinRep::inverse_table(p),
        }
    }

    /// Reduce `v` against the current basis and, if independent, normalize and append it.
    /// Returns true iff `v` extended the span (i.e. the basis grew).
    fn insert(&mut self, mut v: Vec<u64>) -> bool {
        let p = self.p;
        for (bi, b) in self.rows.iter().enumerate() {
            let c = self.pivot_col[bi];
            if v[c] != 0 {
                let f = v[c]; // b[c] == 1 (normalized)
                for k in 0..self.dim {
                    v[k] = (v[k] + (p - f * b[k] % p) % p) % p;
                }
            }
        }
        if let Some(c) = (0..self.dim).find(|&k| v[k] != 0) {
            let inv = self.inv_tab[v[c] as usize];
            for k in 0..self.dim {
                v[k] = v[k] * inv % p;
            }
            self.rows.push(v);
            self.pivot_col.push(c);
            true
        } else {
            false
        }
    }

    fn len(&self) -> usize {
        self.rows.len()
    }
}

/// Smallest gamma-invariant subspace of GF(p)^dim containing `seeds`: close under
/// right-multiplication by each `mats[k]` (the action `u |-> Lambda_p(P^k u)`) until the dimension
/// stabilizes. Returns a basis of `K`. `iter_cap` guards against non-termination (returns `Err`
/// if the dimension has not stabilized within the cap).
///
/// Incremental: a running reduced basis is kept and only the rows added in the previous round are
/// multiplied by each `mats[k]` and reduced against it, instead of re-basing the whole generator
/// set every iteration. The resulting span is identical (closure is monotone), so the kernel and
/// every downstream count are unchanged.
pub fn gamma_closure(
    seeds: &[Vec<u64>],
    mats: &[ModIntMatrix],
    p: u64,
    dim: usize,
    iter_cap: usize,
) -> Result<Vec<Vec<u64>>, String> {
    let mut basis = RunningBasis::new(dim, p);
    for s in seeds {
        basis.insert(s.clone());
    }
    // Frontier = the rows added in the previous round (initially the whole seed basis); only these
    // need their gamma-images examined this round.
    let mut frontier: Vec<Vec<u64>> = basis.rows.clone();
    for _ in 0..iter_cap {
        let mut next_frontier: Vec<Vec<u64>> = Vec::new();
        for b in &frontier {
            for mat in mats.iter() {
                let img = row_times_mat(b, mat, p);
                if basis.insert(img.clone()) {
                    next_frontier.push(img);
                }
            }
        }
        if next_frontier.is_empty() {
            return Ok(basis.rows);
        }
        frontier = next_frontier;
    }
    Err(format!("gamma-closure did not stabilize within iter_cap={iter_cap}"))
}

/// Build the transition matrices `{M_d = Lambda_p(P^d * .)}` at an arbitrary ambient half-width
/// `max_deg`, mirroring `LinRep::compute_mat_for_poly` (allows the padded-ambient C004 setup).
fn build_mats(poly: &LaurentPoly, p: u64, max_deg: usize) -> Vec<ModIntMatrix> {
    let dim = 2 * max_deg + 1;
    let deg = max_deg as i64;
    let mut mats = Vec::new();
    let mut power = LaurentPoly::one(p);
    // Window of exponents the (i, j) loop indexes, identical to LinRep::compute_mat_for_poly.
    let lo = -deg - deg * p as i64;
    let hi = deg + deg * p as i64;
    for _ in 0..p {
        // Dense coefficient array once per power, then O(1) cell reads (vs O(log t) BTreeMap).
        let mut dense = vec![ModInt::zero(p); (hi - lo + 1) as usize];
        for e in lo..=hi {
            dense[(e - lo) as usize] = power.get_coefficient(&e);
        }
        let mut entries = vec![vec![ModInt::zero(p); dim]; dim];
        for i in 0..dim {
            for j in 0..dim {
                let index = (deg - i as i64) - ((deg - j as i64) * p as i64);
                entries[i][j] = dense[(index - lo) as usize];
            }
        }
        mats.push(ModIntMatrix::new(entries, dim, p));
        power = power.mul(poly);
    }
    mats
}

/// The forward-invisible kernel `K` seeded by `seeds` (default `kernel_seeds`) and gamma-closed.
/// Returns a basis of `K` in the ambient `2*max_deg + 1`.
pub fn forward_kernel(
    poly: &LaurentPoly,
    q: &LaurentPoly,
    p: u64,
    seeds: Option<Vec<Vec<u64>>>,
    max_deg: Option<usize>,
    iter_cap: usize,
) -> Result<Vec<Vec<u64>>, String> {
    let max_deg = max_deg.unwrap_or_else(|| {
        let rep = LinRep::for_ct_sequence(poly, q);
        (rep.rank - 1) / 2
    });
    let dim = 2 * max_deg + 1;
    let mats = build_mats(poly, p, max_deg);
    let seeds = seeds.unwrap_or_else(|| kernel_seeds(poly, p, Some(max_deg)));
    gamma_closure(&seeds, &mats, p, dim, iter_cap)
}

/// Forward-reachable row space `U` = span of all `v * M(word)` (the row span of the cheap forward
/// machine). Returns a basis of `U`. `lin_rep` already produces only reachable states.
fn forward_reachable(start: &[u64], mats: &[ModIntMatrix], p: u64, dim: usize) -> Vec<Vec<u64>> {
    let key = |r: &Vec<u64>| -> Vec<u64> { r.clone() };
    let mut states: Vec<Vec<u64>> = vec![start.to_vec()];
    let mut seen: HashMap<Vec<u64>, ()> = HashMap::new();
    seen.insert(key(&states[0]), ());
    let mut k = 0;
    while k < states.len() {
        for mat in mats.iter() {
            let ns = row_times_mat(&states[k], mat, p);
            if !seen.contains_key(&ns) {
                seen.insert(ns.clone(), ());
                states.push(ns);
            }
        }
        k += 1;
    }
    let _ = p;
    LinRep::row_basis(&states, dim, p)
}

/// Express `x` (assumed in span of the ordered `basis`) as coordinates in that basis: returns
/// `coords` with `x = sum coords[i] * basis[i]` over GF(p), or `None` if `x` is not in the span.
/// This is the one linear-algebra primitive `row_basis`/`mat_rank` do not provide (a solve, not a
/// rank). Augmented Gaussian elimination on `[basis | I]`.
fn solve_coords(x: &[u64], basis: &[Vec<u64>], p: u64) -> Option<Vec<u64>> {
    let n = basis.len();
    let dim = x.len();
    if n == 0 {
        return if x.iter().all(|&c| c == 0) {
            Some(vec![])
        } else {
            None
        };
    }
    // Build matrix M whose rows are the basis vectors, augmented with the identity so the row
    // operations that reduce M also record the basis-combination used. We instead transpose the
    // problem: solve basis^T * coords = x. Do plain Gaussian elimination on the augmented system.
    // Columns: dim equations, n unknowns. Build A (dim x n) and rhs x (dim).
    let mut a = vec![vec![0u64; n]; dim];
    for j in 0..n {
        for i in 0..dim {
            a[i][j] = basis[j][i] % p;
        }
    }
    let mut b: Vec<u64> = x.iter().map(|&c| c % p).collect();
    let inv_tab = LinRep::inverse_table(p);
    let mut where_pivot = vec![usize::MAX; n]; // which row pivots each unknown column
    let mut row = 0usize;
    for col in 0..n {
        // find pivot at or below `row`
        let mut sel = None;
        for r in row..dim {
            if a[r][col] != 0 {
                sel = Some(r);
                break;
            }
        }
        let s = match sel {
            Some(s) => s,
            None => continue,
        };
        a.swap(s, row);
        b.swap(s, row);
        let inv = inv_tab[a[row][col] as usize];
        for c in 0..n {
            a[row][c] = a[row][c] * inv % p;
        }
        b[row] = b[row] * inv % p;
        for r in 0..dim {
            if r != row && a[r][col] != 0 {
                let f = a[r][col];
                for c in 0..n {
                    a[r][c] = (a[r][c] + (p - f * a[row][c] % p) % p) % p;
                }
                b[r] = (b[r] + (p - f * b[row] % p) % p) % p;
            }
        }
        where_pivot[col] = row;
        row += 1;
    }
    // consistency: any all-zero row of A must have zero rhs
    for r in 0..dim {
        if a[r].iter().all(|&c| c == 0) && b[r] != 0 {
            return None;
        }
    }
    let mut coords = vec![0u64; n];
    for col in 0..n {
        if where_pivot[col] != usize::MAX {
            coords[col] = b[where_pivot[col]];
        }
    }
    Some(coords)
}

/// Schuetzenberger quotient: the minimal forward linear representation on the reachable space `U`
/// modulo `U cap K`. Returns `(v', {M_d'}, w')` of dimension `dim(U/(U cap K))` = the minimal forward
/// state-space dimension when `K = W^perp`. `kernel` must be `M_d`-invariant (guaranteed by
/// `gamma_closure`).
pub fn quotient_rep(
    v: &[u64],
    mats: &[ModIntMatrix],
    w: &[u64],
    kernel: &[Vec<u64>],
    p: u64,
    dim: usize,
) -> MinimizedRep {
    let u_basis = forward_reachable(v, mats, p, dim);
    // U cap K via U + K dimension count, then build an explicit basis of U cap K.
    // C = basis of (U cap K): vectors of U that lie in span(K). We find them by reducing each
    // U-basis vector against K and collecting a maximal independent set whose K-combination exists.
    // Simpler & robust: pick a basis `c_basis` of U∩K, then complete to a basis of U; the added
    // vectors are the coset representatives.
    let inter = intersect(&u_basis, kernel, p);
    // complete `inter` to a basis of U: greedily add u-basis vectors independent of current span.
    // A single running reduced basis (seeded with `inter`) tests each candidate's independence in
    // one reduction pass, instead of re-basing the whole growing `spanning` set twice per vector.
    let mut spanning = RunningBasis::new(dim, p);
    for iv in &inter {
        spanning.insert(iv.clone());
    }
    let mut reps: Vec<Vec<u64>> = Vec::new();
    for ub in &u_basis {
        if spanning.insert(ub.clone()) {
            reps.push(ub.clone());
        }
    }
    let r = reps.len();

    // quotient coordinates of a vector `x in U`: reduce x against [inter ; reps], read off the
    // coefficients on `reps` (coefficients on `inter` are killed in the quotient).
    let combined: Vec<Vec<u64>> = inter.iter().chain(reps.iter()).cloned().collect();
    let n_inter = inter.len();
    let coords = |x: &[u64]| -> Vec<u64> {
        match solve_coords(x, &combined, p) {
            Some(c) => c[n_inter..].to_vec(),
            None => vec![0u64; r], // x outside U (should not happen for reachable images)
        }
    };

    let row_vec = coords(v);
    let col_vec: Vec<u64> = reps
        .iter()
        .map(|rep| {
            // ct[rep] = entry at exponent 0 = index max_deg = (dim-1)/2; equivalently dot with w.
            let mut acc = 0u64;
            for k in 0..dim {
                acc += rep[k] * w[k] % p;
            }
            acc % p
        })
        .collect();
    let mut mat_func = Vec::new();
    for mat in mats.iter() {
        let mut m = vec![vec![0u64; r]; r];
        for (i, rep) in reps.iter().enumerate() {
            let img = row_times_mat(rep, mat, p);
            let c = coords(&img);
            for j in 0..r {
                m[i][j] = c[j];
            }
        }
        mat_func.push(m);
    }

    MinimizedRep {
        row_vec,
        mat_func,
        col_vec,
        dim: r,
        modulus: p,
    }
}

/// Basis of `span(a) cap span(b)` over GF(p). Uses the Zassenhaus-style identity
/// `dim(A cap B) = dim(A) + dim(B) - dim(A + B)` to find the dimension, and recovers an explicit
/// basis by reducing each `a`-vector against `b` to keep those expressible in `B`.
fn intersect(a: &[Vec<u64>], b: &[Vec<u64>], p: u64) -> Vec<Vec<u64>> {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let dim = a[0].len();
    let ba = LinRep::row_basis(a, dim, p);
    let bb = LinRep::row_basis(b, dim, p);
    let mut union = ba.clone();
    union.extend(bb.clone());
    let dim_union = LinRep::row_basis(&union, dim, p).len();
    let target = ba.len() + bb.len() - dim_union; // dim(A cap B)
    if target == 0 {
        return Vec::new();
    }
    // Recover the intersection: a vector in A is in A cap B iff it lies in span(B). Search the span
    // of A for `target` independent such vectors. Iterate over A's basis and their pairwise sums is
    // not guaranteed complete, so instead solve directly: x in A cap B <=> x = A^T s = B^T u. Build
    // the kernel of [A^T | -B^T] and read off the A-part. We do this via a small null-space solve.
    intersection_via_nullspace(&ba, &bb, dim, p, target)
}

/// Explicit basis of `span(A) cap span(B)`: null-space of `[A^T | B^T]` gives `(s, u)` with
/// `A^T s = B^T u` ... we want `A^T s = B^T u`, i.e. coefficients with `sum s_i a_i = sum u_j b_j`.
/// Build the matrix `M = [a_1 ... a_k | b_1 ... b_l]` (columns are the basis vectors) and find its
/// null space; for each null vector `(s, u)`, `sum s_i a_i` is in the intersection.
fn intersection_via_nullspace(
    ba: &[Vec<u64>],
    bb: &[Vec<u64>],
    dim: usize,
    p: u64,
    target: usize,
) -> Vec<Vec<u64>> {
    let k = ba.len();
    let l = bb.len();
    let ncols = k + l;
    // M is dim x ncols: column c<k is ba[c], column c>=k is -bb[c-k] (so M*(s,u)^T = 0 means
    // sum s_i a_i = sum u_j b_j).
    let mut m = vec![vec![0u64; ncols]; dim];
    for i in 0..dim {
        for c in 0..k {
            m[i][c] = ba[c][i] % p;
        }
        for c in 0..l {
            m[i][k + c] = (p - bb[c][i] % p) % p;
        }
    }
    let null = null_space(&m, dim, ncols, p);
    // For each null vector (s,u), the intersection element is sum_{i<k} s_i * ba[i].
    let mut inter_vecs: Vec<Vec<u64>> = Vec::new();
    for nv in &null {
        let mut x = vec![0u64; dim];
        for i in 0..k {
            let s = nv[i];
            if s != 0 {
                for d in 0..dim {
                    x[d] = (x[d] + s * ba[i][d]) % p;
                }
            }
        }
        if x.iter().any(|&c| c != 0) {
            inter_vecs.push(x);
        }
    }
    let basis = LinRep::row_basis(&inter_vecs, dim, p);
    debug_assert_eq!(basis.len(), target, "intersection dimension mismatch");
    let _ = target;
    basis
}

/// Null space of an `rows x cols` matrix over GF(p): a basis of `{ y : M y = 0 }`.
fn null_space(matrix: &[Vec<u64>], rows: usize, cols: usize, p: u64) -> Vec<Vec<u64>> {
    let mut m: Vec<Vec<u64>> = matrix.iter().map(|r| r.iter().map(|&c| c % p).collect()).collect();
    let inv_tab = LinRep::inverse_table(p);
    let mut pivot_col_of_row: Vec<Option<usize>> = vec![None; rows];
    let mut is_pivot_col = vec![false; cols];
    let mut row = 0usize;
    for col in 0..cols {
        if row >= rows {
            break;
        }
        let mut sel = None;
        for r in row..rows {
            if m[r][col] != 0 {
                sel = Some(r);
                break;
            }
        }
        let s = match sel {
            Some(s) => s,
            None => continue,
        };
        m.swap(s, row);
        let inv = inv_tab[m[row][col] as usize];
        for c in 0..cols {
            m[row][c] = m[row][c] * inv % p;
        }
        for r in 0..rows {
            if r != row && m[r][col] != 0 {
                let f = m[r][col];
                for c in 0..cols {
                    m[r][c] = (m[r][c] + (p - f * m[row][c] % p) % p) % p;
                }
            }
        }
        pivot_col_of_row[row] = Some(col);
        is_pivot_col[col] = true;
        row += 1;
    }
    // free columns parametrize the null space
    let mut basis = Vec::new();
    for free in 0..cols {
        if is_pivot_col[free] {
            continue;
        }
        let mut y = vec![0u64; cols];
        y[free] = 1;
        for r in 0..rows {
            if let Some(pc) = pivot_col_of_row[r] {
                // pivot row: y[pc] = -m[r][free]
                y[pc] = (p - m[r][free] % p) % p;
            }
        }
        basis.push(y);
    }
    basis
}

/// Brief 01 / Algorithm 4.2: build the forward-invisible kernel `K` seeded by `seeds`
/// (default `kernel_seeds(P, p)`), gamma-close it, and return `(K_basis, minimized_rep)` where
/// `minimized_rep` is the Schuetzenberger quotient of `lin_rep(P, Q, p)` by `K` on the reachable
/// forward space. Its dimension equals the minimal forward state-space dimension (= the Hankel rank)
/// precisely when `K = W^perp`, i.e. when the seeding is complete (Cor 6.1); with incomplete seeds it
/// is a CORRECT but possibly NON-minimal automaton. No reverse / observability pass is run.
///
/// `max_deg = None` uses the natural ambient of `lin_rep(P, Q, p)`; pass `Some(m)` to seed into a
/// padded ambient half-width `m` (the C004 setup).
pub fn minimal_dual(
    poly: &LaurentPoly,
    q: &LaurentPoly,
    p: u64,
    seeds: Option<Vec<Vec<u64>>>,
    max_deg: Option<usize>,
    iter_cap: usize,
) -> Result<(Vec<Vec<u64>>, MinimizedRep), String> {
    let max_deg = max_deg.unwrap_or_else(|| {
        let rep = LinRep::for_ct_sequence(poly, q);
        (rep.rank - 1) / 2
    });
    let dim = 2 * max_deg + 1;
    let mats = build_mats(poly, p, max_deg);
    // v = Q as a row, w = unit at exponent 0 (ct functional), in this ambient.
    let v: Vec<u64> = poly_to_row(q, max_deg, p).unwrap_or_else(|| vec![0u64; dim]);
    let mut w = vec![0u64; dim];
    w[max_deg] = 1 % p;
    let seeds = seeds.unwrap_or_else(|| kernel_seeds(poly, p, Some(max_deg)));
    let kernel = gamma_closure(&seeds, &mats, p, dim, iter_cap)?;
    let minimized = quotient_rep(&v, &mats, &w, &kernel, p, dim);
    Ok((kernel, minimized))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn poly(terms: &[(i64, i64)], p: u64) -> LaurentPoly {
        let v: Vec<(i64, u64)> = terms
            .iter()
            .map(|&(e, c)| (e, ((c % p as i64 + p as i64) % p as i64) as u64))
            .collect();
        LaurentPoly::from_vec(v, p)
    }

    // ---- symmetric C004 family: P = t^-1 + 1 + t, seeded into PADDED ambient half-width m. ----
    // Reproduces the Sage reference (p, m, dimK) table exactly.
    #[test]
    fn c004_symmetric_dimk_table() {
        let cases: &[(u64, usize, usize)] = &[
            (2, 2, 2),
            (3, 3, 3),
            (5, 3, 3),
            (7, 4, 4),
            (11, 3, 3),
            (13, 2, 2),
        ];
        for &(p, m, expected_dimk) in cases {
            let p_poly = poly(&[(-1, 1), (0, 1), (1, 1)], p);
            assert!(is_palindromic(&p_poly, p), "P must be palindromic at p={p}");
            let kernel = forward_kernel(&p_poly, &LaurentPoly::one(p), p, None, Some(m), 1000)
                .expect("gamma-closure should stabilize");
            assert_eq!(
                kernel.len(),
                expected_dimk,
                "dimK mismatch at p={p}, m={m}"
            );
        }
    }

    // ---- non-symmetric basket: minimal_dual minrank == hankel_rank, symmetry seed omitted. ----
    #[test]
    fn nonsymmetric_minrank_matches_hankel() {
        // (P-terms, p, expected hankel rank). Sage-confirmed 12/12.
        let cases: &[(&[(i64, i64)], u64, usize)] = &[
            (&[(-1, 1), (0, 1), (2, 1)], 5, 2),
            (&[(-1, 1), (0, 1), (2, 1)], 7, 2),
            (&[(-1, 1), (0, 1), (2, 1)], 11, 2),
            (&[(-1, 1), (0, 1), (1, 2)], 5, 1),
            (&[(-1, 1), (0, 1), (1, 2)], 7, 1),
            (&[(-1, 1), (0, 1), (1, 2)], 11, 1),
            (&[(-1, 1), (2, 1)], 5, 2),
            (&[(-1, 1), (2, 1)], 7, 2),
            (&[(-1, 1), (2, 1)], 13, 2),
            (&[(-2, 1), (0, 1), (1, 1)], 5, 2),
            (&[(-2, 1), (0, 1), (1, 1)], 7, 2),
            (&[(-1, 1), (0, 2), (2, 3)], 7, 2),
        ];
        let mut agree = 0;
        for &(pt, p, expected) in cases {
            let p_poly = poly(pt, p);
            let q = LaurentPoly::one(p);
            // symmetry seed must be omitted for these non-palindromic P
            assert!(
                !is_palindromic(&p_poly, p),
                "test basket P should be non-symmetric (p={p})"
            );
            assert!(
                seeds_symmetry(&p_poly, p, 4).is_empty(),
                "symmetry seeds must be omitted for non-palindromic P (p={p})"
            );

            let rep = LinRep::for_ct_sequence(&p_poly, &q);
            let hankel = rep.hankel_rank(1_000_000).expect("hankel rank");
            assert_eq!(hankel, expected, "hankel-rank anchor mismatch at p={p}");

            let (_k, minimized) =
                minimal_dual(&p_poly, &q, p, None, None, 1000).expect("minimal_dual");
            assert_eq!(
                minimized.dim, hankel,
                "minimal_dual minrank != hankel_rank at p={p}, P={pt:?}"
            );
            agree += 1;
        }
        assert_eq!(agree, 12, "expected 12/12 agreement");
    }

    // ---- the three brief anchors, spelled out explicitly. ----
    #[test]
    fn brief_anchors() {
        // P = t^-1 + 1 + t^2, p=7 -> 2
        let (_k, m1) = minimal_dual(
            &poly(&[(-1, 1), (0, 1), (2, 1)], 7),
            &LaurentPoly::one(7),
            7,
            None,
            None,
            1000,
        )
        .unwrap();
        assert_eq!(m1.dim, 2);
        // P = t^-1 + 1 + 2t, p=5 -> 1
        let (_k, m2) = minimal_dual(
            &poly(&[(-1, 1), (0, 1), (1, 2)], 5),
            &LaurentPoly::one(5),
            5,
            None,
            None,
            1000,
        )
        .unwrap();
        assert_eq!(m2.dim, 1);
        // P = t^-1 + t^2, p=5 -> 2
        let (_k, m3) = minimal_dual(
            &poly(&[(-1, 1), (2, 1)], 5),
            &LaurentPoly::one(5),
            5,
            None,
            None,
            1000,
        )
        .unwrap();
        assert_eq!(m3.dim, 2);
    }

    // ---- palindromicity detection ----
    #[test]
    fn palindromic_detection() {
        assert!(is_palindromic(&poly(&[(-1, 1), (0, 1), (1, 1)], 5), 5));
        assert!(is_palindromic(&poly(&[(-1, 3), (0, 2), (1, 3)], 7), 7));
        assert!(!is_palindromic(&poly(&[(-1, 1), (0, 1), (2, 1)], 5), 5));
        assert!(!is_palindromic(&poly(&[(-1, 1), (0, 1), (1, 2)], 5), 5));
    }

    // ---- gamma-closure termination + iter_cap guard ----
    #[test]
    fn gamma_closure_terminates_and_cap_fires() {
        let p = 5u64;
        let p_poly = poly(&[(-1, 1), (0, 1), (1, 1)], p);
        let m = 3usize;
        let dim = 2 * m + 1;
        let mats = build_mats(&p_poly, p, m);
        let seeds = kernel_seeds(&p_poly, p, Some(m));
        // converges within the cap
        let k = gamma_closure(&seeds, &mats, p, dim, 1000).expect("should stabilize");
        assert!(k.len() > 0 && k.len() <= dim);
        // iter_cap = 0 must fail when seeds are not already closed (one expansion step is needed).
        // Use bare exactness seeds (not yet gamma-closed) so a step is genuinely required.
        let raw = seeds_exactness(&p_poly, p, m);
        let pre = LinRep::row_basis(&raw, dim, p).len();
        // confirm a closure step actually grows the space here, so cap=1-with-no-progress is a real guard
        let closed = gamma_closure(&raw, &mats, p, dim, 1000).unwrap().len();
        if closed > pre {
            assert!(
                gamma_closure(&raw, &mats, p, dim, 0).is_err(),
                "iter_cap=0 must fail before any expansion"
            );
        }
    }

    // ---- per-source seed sanity: exactness seeds are nonempty & live in-ambient ----
    #[test]
    fn seed_sources_basic() {
        let p = 7u64;
        let p_poly = poly(&[(-1, 1), (0, 1), (2, 1)], p);
        let ex = seeds_exactness(&p_poly, p, 3);
        assert!(!ex.is_empty(), "exactness seeds should be nonempty");
        for s in &ex {
            assert_eq!(s.len(), 2 * 3 + 1);
        }
        // support seeds: for a trinomial with full support there may be none; just ensure it runs.
        let _ = seeds_support(&p_poly, p, 3);
    }

    // ---- the minimized rep must COMPUTE the sequence, not merely have the right dimension. ----
    #[test]
    fn minimal_dual_rep_reproduces_sequence() {
        // Evaluate (v', {M_d'}, w') as a forward lsd row machine -- state u starts at v', reads n
        // least-significant digit first with u |-> u * M_d' (row times matrix), output u . w' --
        // and compare to the ambient LinRep on a range of n. A wrong quotient (e.g. collapsing
        // non-equivalent cosets) would still have a plausible `dim` but break the values here.
        // Covers non-symmetric P (Q=1) and symmetric Motzkin P with Q = 1 - x^2.
        fn eval(rep: &MinimizedRep, n: u64) -> u64 {
            let p = rep.modulus;
            let mut state = rep.row_vec.clone();
            let mut nn = n;
            let mut digits = Vec::new();
            while nn > 0 {
                digits.push((nn % p) as usize);
                nn /= p;
            }
            for d in digits {
                // lsd first
                let mat = &rep.mat_func[d];
                let mut ns = vec![0u64; rep.dim];
                for j in 0..rep.dim {
                    let mut acc = 0u64;
                    for i in 0..rep.dim {
                        acc += state[i] * mat[i][j] % p;
                    }
                    ns[j] = acc % p;
                }
                state = ns;
            }
            let mut acc = 0u64;
            for i in 0..rep.dim {
                acc += state[i] * rep.col_vec[i] % p;
            }
            acc % p
        }

        // (P-terms, Q-terms, p)
        let cases: &[(&[(i64, i64)], &[(i64, i64)], u64)] = &[
            (&[(-1, 1), (0, 1), (2, 1)], &[(0, 1)], 7), // non-symmetric, Q=1
            (&[(-1, 1), (2, 1)], &[(0, 1)], 5),         // non-symmetric, Q=1
            (&[(-1, 1), (0, 1), (1, 1)], &[(0, 1), (2, -1)], 5), // Motzkin P (palindromic), Q=1-x^2
        ];
        for &(pt, qt, p) in cases {
            let pp = poly(pt, p);
            let qq = poly(qt, p);
            let (_k, rep) = minimal_dual(&pp, &qq, p, None, None, 1000).unwrap();
            let ambient = LinRep::for_ct_sequence(&pp, &qq);
            for n in 0..120u64 {
                assert_eq!(
                    eval(&rep, n),
                    ambient.compute(n).value,
                    "minimal_dual rep value mismatch at n={n}, p={p}, P={pt:?}"
                );
            }
        }
    }
}
