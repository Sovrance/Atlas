# Candidate-A Convergence and the Scalar E1 Canary (ATLAS-RH-ENG-004)

**Scope:** the P0→P1 bridge. One canonical pole primitive, one promotion
predicate, and exactly one certificate regenerated and released — the scalar
cell entry on `[log 3, log 4]`.

**No RH proof claim is made anywhere in this work order.** The claim scope is
`finite_dimensional_weil_compression`.

---

## 0. What changed

| | before | after |
|---|---|---|
| pole formula | two copies, one of them the rejected block | `src/pole.py`, Candidate A only |
| rejected block | live in `finite_weil` and `weil_fourier_jets` | archival in `src/rejected_pole.py`, CI-guarded |
| promotion rules | inline in `pir_bridge` | `src/promotion.py`, one predicate |
| normalization id | recomputed from code | read from the adjudication artifact |
| SymPy | optional; degraded the artifact when missing | required; fails the job |
| cross-check name | "four-way" | `three_way_internal_crosscheck` |
| scalar E1 | samples + `PENDING_FULL_CELL_COVER` | uniform rigorous bound, `PROMOTED` |

Two independent copies of the rejected Candidate B existed, not one. The pole
assembly in `finite_weil.g0_even_block` was known; the second was
`weil_fourier_jets._dL_g0_even`, which differentiated the rejected block. The
import scan in `tests/test_production_imports.py` found it.

---

## 1. The canonical pole primitive (§1)

`src/pole.py` implements Candidate A and nothing else:

```
G0_ij = E_i^+ E_j^- + E_i^- E_j^+,     E_i^± = ∫_0^L h_i(x) e^{±x/2} dx
```

API: `laplace_plus`, `laplace_minus`, `pole_gram_entry`, `pole_gram_matrix`,
`pole_gram_entry_dL`. Carrier-generic (float / mpmath / sympy / Arb); on an
interval carrier the result is an outward enclosure, and no midpoint is ever
taken to narrow a value.

`normalization.pole_entry` and `finite_weil.g0_even_block` delegate here. The
rejected block lives in `src/rejected_pole.py`, which production may not import —
CI fails on an import *or* on the scale appearing in executable code (string
literals are excluded, because the quarantine reason legitimately quotes it).

---

## 2. Exact invariants (§2)

Checked at the adjudication points `log 3`, `1.1059498113`, `1.20`, `log 4` **and**
at the out-of-cell point `L = 3.5` — a calibration fitted at one `L` is precisely
what WO-RH-17 rejected, so no invariant is allowed to hold only on the research
cell.

* midpoint parity `E^- = ± e^{-L/2} E^+`;
* even/odd pole cross terms vanish (the pole matrix is parity block diagonal);
* odd pivot `G0[q1,q1] = -8 A(L)^2` with `A = L cosh(L/4) - 4 sinh(L/4)`;
* agreement with the real-space route `∫_0^L K_ij(a) 2cosh(a/2) da`, using the
  same `K_ij` as the prime block;
* Candidate B's ratio is `(√3/2)cosh(L/2)`, equal to 1 only at `L = log 3`.

---

## 3. The scalar canary (§5)

### 3.1 The assembled entry

For `h = 1` on `[0, L]`, `G00(L) = G0(L) - Gp(L) + Ginf(L)` with

```
G0(L)   = 16(cosh(L/2) - 1)                                   (Candidate A, exact)
Gp(L)   = Σ_{q=p^k, log q < L} (2 log p / √q)(L - log q)        (exact)
Ginf(L) = (1/π) ∫_0^∞ h_+(t) (2 - 2cos(Lt))/t² dt
```

### 3.2 Why the cell is convex — and why that is Candidate A's doing

This is the load-bearing result. Take the three components in turn.

`Gp` is piecewise linear in `L` with breakpoints exactly at `log q`. The prime
powers below `e^L` on this cell are `{2, 3}`, and `log 3`, `log 4` are the cell's
own endpoints, so **no breakpoint is interior** and `Gp'' = 0` there.

The pole term is a closed form, so `G0'' = 4 cosh(L/2)` exactly.

For the archimedean term, start from the digamma series

```
Re ψ(1/4 + it/2) = -γ + Σ_{n≥0} [ 1/(n+1) - a_n/(a_n² + t²/4) ],    a_n = n + 1/4,
```

and the standard transform `∫_0^∞ cos(Lt) · a/(a² + t²/4) dt = π e^{-2aL}`. The
constant terms contribute only at `L = 0`, so for `L > 0`

```
Ginf''(L) = (2/π) ∫_0^∞ h_+(t) cos(Lt) dt
          = -2 Σ_{n≥0} e^{-2(n+1/4)L}
          = -2 e^{-L/2}/(1 - e^{-2L})
          = -e^{L/2}/sinh(L).
```

Summing, with `r = e^L`:

```
G00''(L) = 4cosh(L/2) - e^{L/2}/sinh(L)
         = 2(r+1)/√r - 2r^{3/2}/(r²-1)
         = 2(r³ - r - 1)/(√r (r² - 1)).
```

That right-hand side is **exactly** the repository's E0-algebraic curvature
`W00''` (`core.scalar_curvature`), already proved positive on the cell by
`scalar.w00_second_positive_on_r_interval`: on `r ∈ [3,4]`, `r³ - r - 1 ≥ 23 > 0`
and `√r(r²-1) > 0`. Each step is verified numerically in
`tests/test_promotion_and_canary.py::ScalarCanaryMath`.

**This is independent evidence for WO-RH-17.** The `4cosh(L/2)` term is produced
by Candidate A's pole and by nothing else. Candidate B multiplies the pole by
`(√3/2)cosh(L/2)`, so its second derivative differs and it cannot reproduce the
certified E0 curvature. The adjudication was decided on the pole's *scale*; the
scalar cell's *geometry* independently agrees.

So `G00` is convex on the cell, and a uniform bound needs only point evaluations —
no interval-`L` quadrature of an oscillatory integrand, which is what makes this
tractable at all.

### 3.3 The quadrature

`Ginf` is evaluated on `[0, T]`, `T = 2·10⁵`, by Arb's rigorous adaptive
integrator `acb.integral`, after continuing the integrand analytically:

* `h_+(z) = (ψ(1/4+iz/2) + ψ(1/4-iz/2))/2 - log π` — equal to `Re ψ(...) - log π`
  on the real axis but, unlike `Re`, analytic, which `acb.integral` requires;
* `(2 - 2cos u)/u²` via an alternating series with an explicit remainder ball near
  `u = 0`, where the quotient form is a `0/0` ball.

Integrating `[0, T]` in one call exhausts the evaluation budget and returns a
non-finite ball. Decade panels converge every time; per-panel radii (~10⁻³⁰) are
recorded in the certificate as the subdivision statistics.

The composite-trapezoid path already in the tree is not adequate here: its
rigorous remainder uses a single global `M2` over `[0, T]` and returns a radius of
~2·10⁴ — useless. That is why this path exists rather than reusing it.

### 3.4 The discarded tail (Lemma T)

```
0 ≤ R_T(L) ≤ (4/π)(h_+(T) + κ(T))/T,     κ(T) = 1 + 1.3/T
```

*Sign.* Every summand of the series above increases with `t`, so `h_+` is
increasing; `h_+(T) > 0` is checked with Arb. With `2 - 2cos ≥ 0` the tail is
therefore **non-negative**, and dropping it can only lower the certified bound —
the direction the claim needs.

*Size (Lemma A′).* `h_+'(t) = (t/2) Σ_n a_n/(a_n²+c²)²` with `c = t/2`, all terms
positive. The summand `f(a) = a/(a²+c²)²` is unimodal with peak
`f(c/√3) = 9/(16√3 c³)`, so sampling at spacing 1 gives
`Σ f(a_n) ≤ ∫_0^∞ f + max f = 1/(2c²) + 9/(16√3 c³)`, hence

```
t h_+'(t) ≤ 1 + 1.29904/t.
```

Note this exceeds 1: numerically `t h_+'(t) = 1.0601` at `t = 2`, approaching 1
from below only for large `t`. The simpler claim `t h_+' ≤ 1` is **false** near
`t = 2`, which is why the constant is carried. Then
`h_+(t) ≤ h_+(T) + κ log(t/T)` and `∫_T^∞ log(t/T)/t² dt = 1/T` give Lemma T.

At `T = 2·10⁵` the bound is `7.24·10⁻⁵`.

### 3.5 From points to a uniform bound

At any interior grid point `p`, convexity puts the tangent below `G00` everywhere,
and brackets `G00'(p)` between the neighbouring secants — both computable from the
point enclosures alone:

```
s⁻ = (lo(p) - hi(p_prev))/(p - p_prev) ≤ G00'(p) ≤ (hi(p_next) - lo(p))/(p_next - p) = s⁺
G00(L) ≥ lo(p) + min over s ∈ {s⁻, s⁺}, L ∈ {a, b} of s·(L - p)
```

Every grid point yields an independently valid global bound; the certificate takes
their maximum. Correctness does not depend on the grid being well chosen — only
the tightness does.

### 3.6 Result

```
G00(L) ≥ 0.0696…  for every L ∈ [log 3, log 4]
```

anchored at `L ≈ 1.2728`, with `T = 2·10⁵`, 200-bit precision, an 18-point grid
refined near the minimizer, and 6 quadrature panels per point.

---

## 4. Regression review (§5)

The historical notebook figure is `0.0753795566117244`. The recovered enclosure at
`L ≈ 1.2784` is `[0.075349745, 0.075422117]` — the notebook value falls **inside**
it.

Under the rejected Candidate B this entry never dropped below ≈0.1276 on the cell
(it ran 0.1276 → 0.4171, monotone), so the notebook minimum was unreachable and
the repository had recorded the regenerated scale as "~0.12–0.13, differs from the
notebook". Candidate A recovers it. The interior minimum that the `L`-dependent
rescaling had flattened away is back.

This is **regression evidence only, never an acceptance constant** — the gate is
`certified_lower_bound > 0`. It is recorded because a normalization that recovers
a previously unreachable historical value is worth noting, not because any number
was fitted to it.

---

## 5. Promotion (§3, §9)

`src/promotion.py` is the single predicate. It refuses a certificate that is
quarantined; that is rigorous but carries no `normalization_certificate_id`; whose
id does not match the active one; whose recorded source hashes no longer match
disk; or that declares E1 without `hard_constraints_certified`.

The active id is read from `normalization_adjudication.json`, never inferred from
a filename, and cross-checked against the content id computed from the frozen
definition. A disagreement is a §14 stop condition, not something to paper over.

The dependency-hash rule earned its keep during this work order: editing
`src/normalization.py` after certifying immediately invalidated the certificate
and the runner re-quarantined it, which is exactly the intended behaviour.

PIR export: 7 promoted (the six prior facts plus the scalar), 4 refused — the
other disputed E1 artifacts.

---

## 6. What stays quarantined

`e1_degree1_log3_log4`, `e1_degree2_compact_log3_log4`, `e1_fourier_T84_points`,
`e1_fourier_T84_uniform_degree2` remain `QUARANTINED_NORMALIZATION_ADJUDICATION`
with their `prior_state` intact. ENG-005 recovers them.

The scalar's release is recorded in `normalization.RELEASED_CERTIFICATES`. It
stays in `QUARANTINED_CERTIFICATES` on purpose: the legacy `certify_scalar_e1.py`
still exists, so an *unauthorised* write must still fail closed. Only a body that
carries `quarantine_released` **and** passes the promotion predicate is left alone.

---

## 7. Limitations (stated, not hidden)

* The bound is uniform over the closed cell but the certificate covers **one
  cell**. Nothing here says anything about other cells or about RH.
* `T = 2·10⁵` with a non-negative discarded tail. The bound is therefore a genuine
  lower bound on the untruncated `G00`, but the *upper* half of each enclosure
  carries the `7.24·10⁻⁵` tail slack.
* Convexity rests on the `Ginf''` identity in §3.2. It is derived, and each step
  is verified numerically, but it is not machine-checked.
* The certified bound `0.0696` is conservative relative to the observed grid
  minimum `0.07535`: the gap is the convexity slope bracket, which shrinks with a
  finer grid near the minimizer and a larger `T`. Tightening it was not needed for
  the gate.
* The cross-check remains three-way internal. `ConnesCvSProjectedProvider`
  quantifies no projection or truncation error and reports `NOT_COMPARABLE`.
