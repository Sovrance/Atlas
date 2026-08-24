# RH / Weil Positivity Program (Research Notebook V2)

This directory integrates the RH Research Notebook V2 into Atlas as the reproducible source behind the Constant Atlas positivity-verifier example (`G00 > 0`, Schur pivot positivity).

## Scope and claim boundary

This program studies finite-dimensional polynomial compressions of the localized Weil quadratic form. **It does not prove the Riemann Hypothesis.** A finite block certificate is evidence only for that stated block, interval, normalization, and cutoff/representation.

Accepted assembly convention:

`G = G^0 - G^p + G^infinity`.

Current certified/research cell: `L in [log(3), log(4)]`.

## Evidence policy

- E0 / SOUND: exact algebraic identities independently re-derived/tested.
- E1 / SOUND: interval-certified numerical statements with reproducible interval engine and certificate.
- E3 / HEURISTIC: floating scans and localization of candidate minima.
- Spectral/quantum analogies: diagnostics only; never proof warrants.

Imported notebook claims are intentionally marked `IMPORTED_PENDING_REGENERATION` until an Atlas agent reproduces them from the code in this repository. Do not promote them merely because they appeared in the prior notebook transcript.

## Normalization is frozen (WO-RH-17, P0) and now lives in one module (ENG-004)

The pole block is **adjudicated**. The adopted form is derived from the explicit
formula,

`G0_ij = Fhat_ij(i/2) + Fhat_ij(-i/2) = E_i^+ E_j^- + E_i^- E_j^+`,  `E_i^± = ∫_0^L h_i(x) e^{±x/2} dx`,

and the legacy even block `(√3/2)(v₊v₊ᵀ+v₋v₋ᵀ)` is **rejected**: it equals the
adopted pole times `(√3/2)cosh(L/2)`, a factor that is 1 only at `L = log 3` — a
calibration fitted at one test point (+8.25 % at `L = log 4`). The odd pivot
`−8A²` was already correct and is unchanged.

`src/pole.py` is the **single** implementation (`laplace_plus`, `laplace_minus`,
`pole_gram_entry`, `pole_gram_matrix`, `pole_gram_entry_dL`). Every production
path routes through it. The rejected candidate has been moved out of production
into `src/rejected_pole.py`, which is archival: `tests/test_production_imports.py`
fails CI if anything under `src/` imports it or spells its scale in executable
code. Two copies of the rejected block existed before this change — the pole
assembly and, separately, its `L`-derivative in the jet module.

See [`docs/NORMALIZATION_ADJUDICATION_v0.1.md`](docs/NORMALIZATION_ADJUDICATION_v0.1.md),
`certificates/normalization_adjudication.json` and
`certificates/normalization_crosscheck.json`. That comparison is a
**three-way internal cross-check** (explicit formula, compact real space, direct
Fourier); the external Connes/CvS provider quantifies no projection or truncation
error, so it reports `NOT_COMPARABLE` and never certifies.

Certificates that depended on the rejected block are
`QUARANTINED_NORMALIZATION_ADJUDICATION` — preserved, not deleted, not relabelled —
and the promotion predicate refuses them. The quarantine is enforced at the point
of write, so re-running a legacy certifier cannot erase it. They must be
**regenerated** under the frozen normalization, never reinterpreted.

## Scalar E1 canary (ATLAS-RH-ENG-004)

The scalar cell entry is the first artifact regenerated under Candidate A and the
only one ENG-004 releases from quarantine. `certificates/e1_scalar_log3_log4.json`
carries a **uniform** rigorous lower bound over `[log 3, log 4]`, computed with
python-flint/Arb — there is no mpmath path that may emit E1.

The bound rests on convexity, which is algebraic rather than numeric. `Gp` is
piecewise linear with breakpoints exactly at the cell endpoints, the pole
contributes `G0'' = 4cosh(L/2)`, and the archimedean term's cosine transform gives
`Ginf'' = −e^{L/2}/sinh L`, so

`G00''(L) = 4cosh(L/2) − e^{L/2}/sinh(L) = 2(r³−r−1)/(√r(r²−1))`,  `r = eᴸ`,

which is exactly the repository's E0 curvature `W00''`, already proved positive on
the cell. That identity is itself independent evidence for the adjudication: the
`4cosh(L/2)` term is produced by Candidate A's pole and by nothing else, so the
rejected block cannot reproduce the certified curvature.

Regression note: the historical notebook minimum `0.0753795566…` — which the
rejected calibration had put out of reach, holding this entry above ≈0.1276 on the
cell — falls **inside** the recovered enclosure. It is recorded as regression
evidence only; the acceptance gate is the positive bound, never a fitted constant.

Degree-1, compact degree-2 and both T=84 E1 artifacts remain quarantined. ENG-005
recovers them.

```bash
python3 scripts/derive_normalization.py           # derivation + adjudication certificate (requires SymPy)
python3 scripts/run_normalization_crosscheck.py   # three-way internal cross-check (--no-arch to skip the slow route)
python3 scripts/quarantine_normalization.py       # idempotent quarantine
python3 scripts/run_rigorous_scalar.py            # the rigorous scalar path, in order
```

## Core E1 recovery (ATLAS-RH-ENG-005)

ENG-004 recovered the scalar cell. ENG-005 recovers the rest of the Candidate-A
chain and rebuilds the T=84 topology from scratch.

**Degree-1 and compact degree-2** are cutoff-free — no frequency truncation to
bound away — because the archimedean term now has an exact real-space form:

`Ginf_ij(L) = (K(0)/2)h₊(0) + ∫₀^L [K(0)−K(u)]·w(u)du + K(0)·S(L)`

with `w(u) = e^{−u/2}/(1−e^{−2u})` and `S(L) = Σₙ e^{−(2n+1/2)L}/(2n+1/2)`, using
the same `K_ij` as the prime block. This matters because the frequency-space
definition is an oscillatory half-line integral: fine at a point, useless on an
`L`-interval, which is why ENG-004 had to reach for convexity. The real-space form
has no oscillation, so ordinary interval subdivision works — and it is ~1000×
faster. The two routes agree exactly on the fast-decaying entries and differ on
the slow-decaying ones by precisely the expected truncation tail, so each
cross-checks the other.

The naive interchange gives `−∫₀^L K(u)w(u)du`, which diverges; the constant part
of `h₊` contributes a delta at the origin that the naive swap drops. Keeping it
produces the `K(0)−K(u)` numerator, which vanishes linearly and cancels the `1/u`.

**T=84** is a different object, not an approximation of these: its archimedean
term stops at `T = 84` by definition, so the frequency route *is* the definition
there. Its topology was rescanned fresh under Candidate A — the rejected
Candidate-B monotonicity topology is not reused, and the superseded scan is kept
under `certificates/history/` as provenance.

All `L`-derivatives are **exact support-length jets**, never finite differences:
`d_L^n H₀ = (it)^{n−1}e^{itL}`, `d_L Hb = ∫₀^L x e^{itx}dx`, `d_L² Hb = L e^{itL}`,
with binomial convolution for the Gram entries. Finite differences appear only as
a *check* on those jets — and earned their place: they caught a `d²/dL²(L³/6)`
coefficient written as `L/2` instead of `L`, which threw `d²O1` off by 0.70.

Uniform bounds come from one shared adaptive interval cover
(`src/interval_cover.py`), which assumes no topology at all — §8 forbids
precommitting to convexity or monotonicity, so each certificate records the
topology it actually established.

At `T = 84` a second, independent warrant sits alongside that cover: instead of
merely bounding `E2`, it **locates the minimiser**. Certified bisection on
`sign(E2')` pins `L*` to a nine-digit interval; `E2'' > 0` on a window around it
makes `L*` the unique critical point there and a strict minimum; an interval
cover of that short window bounds `E2` almost as tightly as a point evaluation;
and certified derivative signs (`E2' < 0` to the left, `E2' > 0` through a band
to the right) show nothing outside the window can be lower. The headline bound is
the sharper of the two warrants, floored at what the plain cover alone proves, so
the second warrant can only improve the number.

The derivative argument governs a band, not the whole cell: the `E2''` enclosure
on a box of radius `r` carries a ~`300 r` dependency blow-up that cannot be
centred away (that would need an exact third jet, and a finite-difference one is
forbidden in an E1 path), so past the band it is cheaper and no less rigorous to
bound `E2` directly. Each certificate names the interval its warrant governs, and
the intervals abut exactly.

```bash
python3 scripts/run_rigorous_chain.py      # the whole chain, in canonical order
```

### Rigorous dependencies

SymPy and python-flint are **required** for the rigorous/research path, not
optional extras. A missing one now fails the job rather than degrading an artifact
to weaker evidence and dirtying the tree:

```bash
pip install -r requirements-rigorous.txt
```

`scripts/run_rh_weil_suite.py` is the fast path. Passing it does **not** mean the
rigorous certificates are current — it does not re-derive the scalar canary. Read
`scripts/run_rigorous_scalar.py`'s exit code before believing an E1 claim is fresh.

## Inertia, rank-trace and spectral moments (ATLAS-RH-ENG-006)

Through ENG-005 every block was asked one question — is it positive? That is a
single bit, and when the answer is no the run ends with nothing. ENG-006 adds
three channels so an indefinite block still yields a rigorous result:

* **inertia** (`inertia/`) — the full signature `(n₊, n₋, n₀)` by interval
  Hermitian LDL congruence, checked against an independent oracle that reads the
  inertia off the characteristic polynomial by Descartes' rule of signs (exact,
  because a symmetric matrix is real-rooted);
* **rank–trace** (`ranktrace/`) — `rank(P) ≥ 2 tr(P) + 4 tr(Q) − 4b − ‖P+Q‖²_HS`
  with its four hypotheses enforced, not assumed;
* **spectral moments** (`moments/`) — `m₁..m₄` as traces of matrix powers, fed to
  the existing Atlas B1 truncated-moment solver rather than a second one.

Three things these refuse to do. An interval result never claims an exact zero,
so singular inputs come back `INCONCLUSIVE` rather than `n₀ = 1`. A rank–trace
call with any unverified hypothesis produces no number at all. And "the moments
force PSD" comes back `INSUFFICIENT_INFORMATION` even for a positive definite
matrix, because PSD-ness of a truncated localizing matrix is necessary and not
sufficient — only the refuting direction is available from four moments.

**The odd degree-3 block** `[[G_q1q1, G_q1b3], [G_q1b3, G_b3b3]]` is the first
live workload, and it came out **positive definite**: inertia `(2,0,0)` uniformly
on `[log 3, log 4]`, one stratum and no transition regions, with
`O1 ≥ 1.5331e-02` and `det ≥ 1.0731e-06`. The work order did not require this —
an inertia stratification would have counted as success — and the same machinery
would have produced one had any part of the cell been indefinite.

Its two kernels are re-derived from the basis with SymPy on every run:

`K_q1b3 = (L−a)²(L³ + 2L²a − 12La² − 6a³)/60` and
`K_b3b3 = (L−a)³(L⁴ + 3L³a − 15L²a² − 18La³ − 6a⁴)/420`.

Every active prime-shift block on the cell is indefinite — determinant negative,
inertia `(1,1,0)` — which is kept as a regression test. It is why the assembled
entry has to be bounded as a whole and why termwise PSD domination is unavailable.

Positivity and inertia are distinct content kinds. An inertia certificate never
satisfies a consumer requiring PSD, even when its signature is `(2,0,0)`; the
degree-3 artifact answers such a consumer as a *positivity* certificate carrying
certified bounds, with the inertia object nested inside it still refusing.

```bash
python3 scripts/certify_degree3.py            # scan, E1 result, moments
python3 scripts/ci_inertia.py --gate fast     # exact gates, no python-flint needed
python3 scripts/ci_inertia.py --gate rigorous # interval gates, python-flint required
python3 scripts/report_information_comparison.py
```

## Current priority

ENG-005 recovered the E1 chain and ENG-006 delivered the inertia/rank-trace/moment
channels and the degree-3 pilot. Next is ENG-007: formalize the stabilized theorem
boundary in Lean — Sylvester inertia under congruence, the 2x2/3x3 criteria, the
rank-trace theorem, certificate semantics, and selected degree-3 exact identities.
Only after that should the program widen to additional prime-power cells or higher
polynomial degree at scale.

**Executed in-repo so far (see `certificates/work_order_status.json`):**
- WO-RH-01/03/04 — exact identities, f1 audit, bubble block (E0)
- WO-RH-02 — scalar cell `[log 3, log 4]` algebraic `W00''>0` (E0)
- WO-RH-05 — stable `H0`/`Hb` + L-jets; E3 energy-probe scan only (interval E1 **open**)
- WO-RH-06/07 — regenerated E0 certificates + dedicated runner (does not expand root CI)

Run:

```bash
python math/rh_weil/scripts/run_rh_weil_suite.py
```

See `AGENT_INSTRUCTIONS.md` and `notebook/RH_RESEARCH_NOTEBOOK_V2_INTEGRATION.md`.

Optional external oracle (not required at runtime): `external/` wraps
`connes-cvs` for independent cross-checks of shared Weil ingredients. See
`external/CONNES_CVS_MAPPING.md` and `external/PROVENANCE.md`. Do not compare
Galerkin matrix entries to Atlas polynomial Gram blocks by index.

Regenerate the external cross-validation certificate (research env):

```bash
pip install 'connes-cvs==0.3.1' python-flint mpmath
python math/rh_weil/scripts/run_connes_cvs_crosschecks.py
```

