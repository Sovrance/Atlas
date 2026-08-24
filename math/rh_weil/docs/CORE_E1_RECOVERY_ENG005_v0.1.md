# Core E1 Recovery and Candidate-A T=84 Reconstruction (ATLAS-RH-ENG-005)

**Baseline:** ENG-004, merged. **Scope:** recover degree-1 and compact degree-2,
then rebuild the direct-Fourier T=84 topology and E1 from scratch under
Candidate A.

**No RH proof claim.** Claim scope is `finite_dimensional_weil_compression`.

---

## 1. The blocker, and what removed it

ENG-004 certified the scalar cell by convexity, because the archimedean term

```
Ginf_ij(L) = (1/π) ∫₀^∞ h₊(t) Re(conj(H_i) H_j) dt
```

is an oscillatory half-line integral. At a *point* `L` that is fine. On an
`L`-**interval** it is hopeless: `cos(Lt)` has argument width `rad(L)·t`, so a box
wide enough to be worth taking destroys the enclosure long before `t` reaches the
cutoff. ENG-004 escaped via a closed form for `G00''`. That escape does not
generalise — `E2 = G00·Gbb − G0b²` is a product difference with no such form.

The way out is to stop integrating over frequency. Using
`Re(conj(H_i)H_j) = ∫₀^L K_ij(u)cos(tu)du` — the **same** `K_ij` as the prime
block — together with the transform proved in §2 and the digamma series:

```
Ginf_ij(L) = (K(0)/2)·h₊(0) + ∫₀^L [K(0) − K(u)]·w(u) du + K(0)·S(L)

w(u) = e^{−u/2}/(1 − e^{−2u}),   S(L) = Σ_{n≥0} e^{−(2n+1/2)L}/(2n+1/2)
```

Every piece is benign on an interval: `K` is an exact polynomial, the integral is
over the *compact* `[0, L]` with no oscillation, and `S` converges geometrically
(ratio `e^{−2L} < 0.09` on this cell).

**Why the naive interchange fails.** Swapping sum and integral directly gives
`−∫₀^L K(u)w(u)du`, which **diverges**: `w(u) ~ 1/(2u)` while `K(0) ≠ 0`. The
constant part of `h₊` contributes a delta at `u = 0` that the naive swap silently
drops. Keeping the `F(0)` terms produces the `K(0) − K(u)` numerator, which
vanishes linearly at the origin and cancels the `1/u` exactly.

**Two independent routes.** The frequency-space route survives as a cross-check.
They agree to 2e-16 on the fast-decaying entries `(one,b)` and `(b,b)`, and on the
slow-decaying `(one,one)` and `(q1,q1)` they differ by exactly the expected
`T`-truncation tail, with the right sign. The real-space route is also ~1000×
faster (0.0s against 36–64s per entry).

---

## 2. Reproducible `Ginf''` (§1)

`src/curvature_derivation.py` re-derives ENG-004's load-bearing identity with
SymPy on every run. Five steps are exact symbolic identities:

1. `Re[1/(n+1/4+it/2)] = a_n/(a_n²+t²/4)`, `a_n = n+1/4`
2. `∫₀^∞ cos(Lt)·a/(a²+t²/4) dt = π e^{−2aL}`
3. the geometric sum, as an exact **finite** partial-sum identity plus the
   elementary `N → ∞` limit
4. `−2e^{−L/2}/(1−e^{−2L}) = −e^{L/2}/sinh(L)`
5. `4cosh(L/2) − e^{L/2}/sinh(L) = 2(r³−r−1)/(√r(r²−1))`, `r = e^L`

Step 3 is split deliberately. SymPy will not evaluate the infinite sum (it does
not assume `|e^{−2L}| < 1`), and asserting the closed form without the ratio
condition would be exactly the unearned step §1 warns against. Three steps needed
`rewrite(exp)` to compare, since SymPy returns the transform in sinh/cosh form.

**What is not proved.** The termwise interchange is *not* machine-verified and is
not claimed to be. The integral is only conditionally convergent (`h₊` grows like
`log t`), so the exchange needs an Abel/Cesàro argument. It is recorded under
`analytic_hypotheses`, the report is classified **E0 for the algebra only**, and
an independent convergent-series regression — sharing no code with the SymPy path
— checks the closed form numerically. §1 forbids overclassifying this as stronger
than the actual proof technology, and it is not.

---

## 3. Canonical panel integration (§2)

`src/rigorous_integration.py` is now the one rigorous quadrature. `acb.integral`
over a whole range in a single call exhausts its evaluation budget on this
integrand and returns a **non-finite ball** — sound but useless. Panels converge
every time, so panels are the canonical path, not a fallback.

The schedule is a pure function of `T` (T=84 gets the fixed dyadic split the spec
names), is **validated to cover `[0, T]` exactly**, and is recorded in every
certificate. A non-finite panel raises rather than returning an infinite
enclosure, since a bound derived from one is vacuously true.

That coverage check immediately caught a live bug in ENG-004's private schedule:
its decade edges were zipped pairwise and filtered by `hi > lo`, so for any `T`
below its last hard-coded edge it integrated *past* `T` — `T = 20000` was
integrated over `[0, 100000]`. The shipped canary used `T = 200000`, where the
schedule happens to be exact, so the merged ENG-004 certificate is unaffected;
nothing had checked that.

The old global composite-trapezoid path may not emit E1 for this program. Its
rigorous remainder uses a single global `M2` and returns a radius ~2e4 on this
integrand — six orders of magnitude larger than the quantity being bounded.

---

## 4. Tail lemma (§3)

`κ(T) = 1 + 1.3/T`, from a unimodal sum-vs-integral bound giving
`t·h₊'(t) ≤ 1 + 1.29904/t`.

The tempting `t·h₊'(t) ≤ 1` is **false** near the low end. A rigorous enclosure at
`t = 2` gives `[1.06010377, 1.06010377]`, so the bound provably fails there. It is
now regression-guarded (`invalid_assumption_is_rejected`) so it cannot be quietly
re-adopted, and `κ` is carried as a function of the tail domain rather than
replaced by 1.

---

## 5. Degree-1 and compact degree-2 (§4/§5)

Both cutoff-free, both uniform on the closed cell, by adaptive interval cover:

* `O1(L) = G[q1,q1] ≥ 1.4978e-2`, basis `q1 = x − L/2`,
  `K_q1q1(a;L) = (L−a)(L²−2La−2a²)/6`
* `E2(L) = G00·Gbb − G0b² ≥ 2.0130e-6`, basis `{1, b}`, `b = x(L−x)`

Neither approaches zero, so §15's stop condition does not fire. `E2`'s true
minimum is ~4.6e-6 at the left endpoint, so the cover is asked to clear 2e-6
rather than merely 0 — a bound thirty times below the minimum would be a weak
certificate for a quantity this small.

`D2 = E2 + L²·G00·O1` and `det(G_deg≤2) = O1·E2` verified; both follow from the
pole and prime blocks being parity block diagonal.

**Two implementation traps, both caught before use.** Integrating to `sup(L)` on
an `L`-ball is wrong — the limit itself depends on `L`, so for `L' < sup(L)` the
true integral stops earlier and the enclosure includes mass the real integral
never sees; fixed by substituting `u = L·s` onto the fixed `[0,1]`. And centring
only the archimedean term is not enough: `G0`, `−Gp` and `Ginf` are individually
of order 0.1 and cancel heavily (`Gbb ~ 3.6e-5`), so evaluating them separately on
a box makes their widths add while the true variation stays tiny
(`d/dL Gbb ~ 1e-3`, four orders below the per-block slopes). The mean-value form
is applied to the whole assembled entry.

---

## 6. T=84: a different object (§6/§7/§9/§10)

The T=84 matrix is **not** an approximation of the cutoff-free entries. Its
archimedean term stops at `T = 84` by definition, so the frequency route *is* the
definition there — and it is affordable, because 8 panels suffice and `t ≤ 84`
keeps interval-`L` usable.

### Exact support-length jets (§9)

No finite differences in any E1 path:

```
d_L^n H₀  = (it)^{n−1} e^{itL}                        (n ≥ 1)
d_L   Hb  = ∫₀^L x e^{itx} dx
d_L² Hb  = L e^{itL}
d_L^n Hb  = e^{itL}[L(it)^{n−2} + (n−2)(it)^{n−3}]    (n ≥ 3)
d_L   Hq1 = (L/2)e^{itL} − H₀/2
d_L² Hq1 = (L/2)(it)e^{itL}
```

with binomial convolution for the Gram entries. The pole jets come from
`pole.pole_gram_entry_dL` / `_d2L` (closed forms; every basis element is linear in
`L`, so `d²_L E^± ` reduces to boundary terms) and the prime jets from exact
kernel coefficient expansions. Candidate-B derivative code is unreachable.

Finite differences appear **only as a check** on those jets — and earned their
keep: they caught `d²/dL²(L³/6)` written as `L/2` instead of `L`, which threw
`d²O1` off by exactly `Σ_q w_q(L − log q) = 0.7028` at `L = 1.25`. Nothing else in
the pipeline would have noticed.

### Topology, chosen not assumed (§8)

The previous Candidate-B monotonicity topology is **not** reused; the superseded
scan is preserved under `certificates/history/` as rejected-normalization
provenance. The fresh scan reports stationary points and curvature changes of
`E2` as *apparent* grid features — E3 evidence, never a warrant.

Two independent warrants are then produced, and the certificates carry both.

**Warrant 1 — exhaustive cover.** A branch-and-bound over a box partition of the
closed cell that assumes **no** topology at all: every box has to clear the
target on its own. This is the fallback and the tie-breaker. If it ever
disagreed with the scan, it wins.

**Warrant 2 — interior minimum.** The exhaustive cover proves a bound but says
nothing about *where* the minimum is, and it pays for that: it has to resolve a
very flat minimum with boxes narrow enough to separate `E2 ~ 3.4e-6` from zero
everywhere at once. The interior-minimum argument locates the minimiser instead:

1. **Isolate.** Certified bisection on `sign(E2')` to a bracket `[a, b]` with
   `E2'(a) < 0 < E2'(b)`, both signs from rigorous Arb enclosures. Point
   enclosures of `E2'` come back with width ~3e-13, so the bracket closes to
   ~5e-10 and `L*` is pinned to nine digits.
2. **Uniqueness and minimality.** `E2'' > 0` on a window `W` around `[a, b]`, so
   `E2'` is strictly increasing there: `L*` is the *only* critical point in `W`
   and it is a strict local minimum.
3. **Basin bound.** `E2 >= m` on `W` by interval cover. `W` is short, so this is
   nearly a point evaluation and `m` is within a hair of the true minimum.
4. **No lower values elsewhere.** `E2' < 0` on `[log 3, W_lo]` and `E2' > 0` on
   `[W_hi, W_hi + 0.03]`, so `E2` falls strictly onto `W` from the left endpoint
   and rises strictly off it, and no critical point exists in between.

The headline bound is `min(m, cover bound outside the governed interval)`,
floored at the plain whole-cell cover bound so §8's second warrant can only
sharpen the number, never weaken it.

**Why the derivative argument stops at 0.03 rather than reaching `log 4`.** The
`E2''` enclosure on a box of radius `r` carries a dependency blow-up of ~`300 r`
on top of the true `|E2''| ~ 1.4e-3` — and unlike `E2'`, it cannot be centred,
because centring needs an exact third jet and `pole.pole_gram_entry_d2L` is the
last closed form available (a finite-difference third jet is exactly what §9
forbids in an E1 path). So the centred `E2'` form separates from zero only while
`300 r^2 < |E2'|`; out where `|E2'| ~ 4e-5` that forces `r < 1.2e-4`, about 1200
boxes to reach `log 4` — more than covering `E2` directly costs, and each box
dearer. The honest split is: derivative signs near the minimiser, where they are
cheap and where the direct cover is weakest, and the direct cover further out,
where `E2` has climbed well clear of `m`. The certificate names the interval each
warrant governs, and they abut exactly.

**Centring the derivative was not optional.** Evaluated raw on an `L`-box of
radius `r`, `E2'` comes back with halfwidth ~`25 r` — at `r = 1e-4` that is
`1.1e-3`, two orders above the `|E2'| <= 1e-5` the sign conditions have to
resolve. No sign could be certified anywhere near `L*`. The same mean-value form
the entries already use — exact point value at the midpoint plus `r` times an
enclosure of the exact second jet — replaces `25 r` with `r(|E2''| + O(r))`, or
`3e-6` instead of `1.1e-3` at that radius.

---

## 7. Limitations (stated, not hidden)

* Everything covers **one cell**, `[log 3, log 4]`. Nothing here says anything
  about other cells or about RH.
* The real-space archimedean form is E0 algebra **conditional on** the interchange
  hypothesis of §2. The frequency route is an independent numerical check, not a
  proof of the interchange.
* `E2` is genuinely small (~1e-6 to 1e-5) because the basis `{1, b}` is nearly
  degenerate on this cell. The bound is positive and uniform, but it is not a
  well-conditioned quantity, and the cover has to work to separate it.
* The T=84 certificates depend on the truncation being part of the *definition*.
  They say nothing about the untruncated matrix.
* The scan is a grid. Features between grid points are not excluded by it — which
  is precisely why the uniform certificate does not rest on it.
* The interior-minimum argument governs `[log 3, W_hi + 0.03]`, not the whole
  cell. Beyond that band the only warrant is the exhaustive cover. The two
  intervals abut exactly and are both recorded, but the sharper `E2 >= m`
  statement is a statement about the band, and the headline bound is the weaker
  of the two where they meet.
* `E2''` is certified raw, not centred, because centring it needs an exact third
  jet the pole primitive does not provide. That is why the curvature window is
  1e-4 wide rather than something more comfortable: it is the widest interval on
  which the raw `E2''` enclosure still separates from zero at reasonable cost.
* Uniqueness of the critical point is proved on the curvature window and, off it,
  only out to the ends of the sign-certified regions. Nothing here excludes a
  critical point beyond the band — only values below the bound, which the
  exhaustive cover excludes directly.
* Integrator settings (140-bit precision, `rel_tol` 1e-25) were chosen for speed;
  they leave the enclosures ~1e-26, far tighter than the ~1e-10 the claims need,
  but they are a deliberate accuracy/time trade and are recorded per certificate.
