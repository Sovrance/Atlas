# The First Higher-Dimensional Certified Weil Block (ATLAS-RH-ENG-008)

**Baseline:** ENG-007, merged (`015519e`). **Scope:** convert the ENG-007 3×3
even pilot into a rigorously certified spectral result on `L ∈ [log 3, log 4]`.

**No RH proof claim.** Claim scope is `finite_dimensional_weil_compression`.
Everything below concerns one finite block, on one interval, under one
normalization.

---

## 1. The result

The 3×3 even Weil block over the frozen basis `{1, b, b²}`,

`e0 = 1`,  `e1 = b = x(L−x)`,  `e2 = b² = x²(L−x)²`,

is **positive definite for every `L` in `[log 3, log 4]`**: inertia `(3, 0, 0)`,
one stratum, no transition regions.

Two independent routes reach that conclusion and both are kept.

| route | what it does | outcome |
|---|---|---|
| interval LDL* congruence | eliminates the preconditioned block, stratifies the cell by signature | `(3,0,0)`, 2920 boxes, max depth 3 |
| Sylvester's criterion | three separate adaptive interval covers of the leading principal minors | `Δ1, Δ2, Δ3` all uniformly positive |

They share the assembly and nothing after it. Agreement is therefore evidence;
disagreement would have been a stop condition, and the rigorous chain fails if
the two artifacts ever report different signatures.

## 2. Why this block and not another

Writing `u = x − L/2`, the even sector about the midpoint is `span{1, u², u⁴}`
and the three basis elements supply it in order — `b = L²/4 − u²` contributes
`u²`, and `b²` contributes `u⁴`.

The odd sector cannot be extended at degree 3 at all. Through degree 3 it is
exactly `span{u, u³}`, and `q1³ = (L²/4)q1 − b3` already lies in it. That is why
ENG-006 stopped at a 2×2 odd block, and why ENG-007's pilot chose the even
sector.

The choice also matters for continuity with what was already certified: this
block's second leading minor **is** the `E2 = G00·Gbb − G0b²` that ENG-004 and
ENG-005 certified. At `L = log 3` the 3×3 assembly reproduces it to every digit,
`4.5617133218400482562565030758111264343806077379e-6`, which is the
reconciliation §WO-RH-51 requires.

### A disagreement that was not one

The ENG-007 pilot preview reported different numbers for the same entries, and
that looked like the stop condition "independent assemblies disagree". It is
not. The pilot computed the **T = 84 truncated** object; this work computes the
**cutoff-free** one. ENG-005 recorded that those are different objects rather
than approximations of each other, and the two sets of numbers differ exactly as
that implies.

## 3. The preconditioner

The raw block spans ten orders of magnitude — entries `O(1e-1)` down to a third
leading minor around `1.4e-11` — and the elimination separates badly on it.
§WO-RH-47 permits a diagonal preconditioner provided the congruence is certified
and the mathematical claim is unchanged.

`D = diag(2^{-e₀}, 2^{-e₁}, 2^{-e₂})` with `eₖ = round(log₂ √(G_kk))`, **frozen
for the whole cell** at `(-2, -6, -9)`.

Three properties, each load-bearing:

* **Exactly invertible.** Diagonal, nonzero dyadic entries,
  `det D = 2^{-Σeₖ}`. No enclosure is involved, so invertibility is not
  something the certificate asks a reader to accept on numerical grounds.
* **Exactly applied.** Scaling an Arb ball by a power of two shifts its exponent
  and changes nothing else — not the midpoint's significand, not the radius'. So
  `DᵀGD` adds **no width at all**. A general Jacobi scaling by `1/√(G_kk)` would
  have to round, and would inflate every entry it touched.
* **Frozen, not per box.** Choosing the exponents from each box's own diagonal
  would give slightly better scaling and would mean the certified numbers
  described a different matrix on every box. The sign conclusion would survive —
  every admissible `D` is a positive diagonal scaling — but a uniform numerical
  bound would not be a statement about any single object.

It buys ten orders of magnitude: the rescaled minors sit at `O(1)` where the raw
third minor is `~1e-11`. That this costs nothing mathematically is a **theorem**,
not an assumption: `AtlasRH.posIndexAtLeast_congruence_iff` and
`AtlasRH.rank_congruence`, proved in ENG-007 against a pinned Mathlib, say
congruence by an invertible matrix preserves the positive index, the negative
index and the rank.

## 4. Derived algebra, not hand tables

Two pieces of machinery were generalized before the block could be assembled at
all, and both replaced hand-written per-element tables.

**The overlap kernels** (`src/basis_algebra.py`). `K_ij(a; L)` is now computed as
an exact bivariate polynomial in `(a, L)` from the one genuinely primitive thing
— each basis element's monomial coefficients as exact polynomials in `L`. The
three runtime forms (closed form, coefficients in `a`, and `d/dL` of those) are
read off that single table. Previously they were three hand-written tables in
three modules, one entry per pair, each raising `KeyError` for an unknown pair.

**The second `L`-derivative** (`src/pole.py`). `_laplace_d2L` was a hand-written
closed form per element. It now evaluates the general expression

`F'' = H'(L)e^{sL/2} + (s/2)H(L)e^{sL/2} + h_L(L;L)e^{sL/2} + ∫₀^L h_LL(x;L)e^{sx/2}dx`

from the coefficient tables.

<!-- docs-check: superseded-quote start -->
The ENG-007 record said `pole.py` "drops the second integral because
`d²_L h = 0`". That was too strong: `b3` is quadratic in `L` and its integral
term was already carried. The actual defect was narrower and worse — the
provider was hand-specialized per element and simply raised `KeyError` for
anything not in its table, and `b2` was not in the table.
<!-- docs-check: superseded-quote end -->

Both generalizations were checked two ways. They reproduce every retired table
entry **exactly**, in exact rational arithmetic at several `(a, L)`; and each is
independently verified against direct symbolic integration or differentiation of
its own defining integral. A generalization that faithfully reproduced a wrong
table would pass the first check and fail the second.

### What the refactor cost, and what it did not

The refactor made every E1 certificate stale by construction, so the whole chain
was regenerated. The scalar, T=84 uniform and degree-3 bounds came back
**bit-identical**. The degree-1 and degree-2 bounds agree with their predecessors
to ten and nine significant figures respectively, and the reason is understood
rather than assumed:

* the exact rational kernel coefficients are **provably identical** — the tests
  check them in exact arithmetic at several `(a, L)`;
* the covers are identical — same box counts (24 and 178) and same maximum
  depths (0 and 3);
* the pole side is **bit-identical**, verified by running the old and new
  `pole_gram_entry`, `_dL` and `_d2L` side by side;
* what differs is the association order when those exact coefficients are
  evaluated on an Arb ball. Midpoints agree to 47 digits; several of the new
  radii are *tighter* than the old ones.

That is a rounding difference with a cause, not a change in the mathematics.

### Two regressions the refactor did introduce, and how they were caught

Both were the same kind of error and neither was a mathematical one. In each
case the generalization computed the *right number* and a *wider interval*, and
the width is what a certificate is made of.

#### The endpoint, and cancellation on the carrier

The first version of the generalization computed `h(L; L)` by substituting
`x = L` on the carrier. `b`, `b3` and `b2` all vanish identically there, but the
cancellation then happens *in interval arithmetic*: on a box of radius `1e-2`,
`b(L; L) = L·L − L·L` came back as a ball of radius `2.2e-2` instead of exact
zero. That width propagated through every derivative bound built on it, and it
moved the degree-1 bound down and the degree-2 bound up by 6%.

It was caught because the certified numbers are checked against the README by
`scripts/check_docs.py`, so a moved bound fails a gate rather than passing
quietly.

The fix is to do the cancellation in exact rational arithmetic and let the
carrier see only the simplified polynomial — an identically zero quantity is an
empty polynomial and evaluates to exact zero. `basis_algebra.endpoint_poly` and
its two derivatives now hold those, and `pole.py` evaluates them.

#### The kernel, and a factor that had been there all along

Every retired closed form displayed a factor of `(L − a)^m` — `m = 1` for
`K_one,one` and `K_q1,q1`, rising to `5` for `K_b2,b2`. That is not decoration:
each kernel integrates over `[0, L − a]`, so the factor is structural. The
derived engine produced the same polynomial but in *expanded* form, and
evaluated it by Horner in `a`.

On an exact carrier those are the same number. On a ball they are not. Horner
treats each occurrence of `a` as an independent quantity, so the correlation
that makes `(L − a)^5` small when `a` is near `L` is thrown away term by term.
The measured cost, on the box the certificates actually bind at:

| quantity | cost of evaluating expanded instead of factored |
|---|---|
| `K_q1,q1` enclosure radius | 3× wider |
| `K_b3,b3` enclosure radius | 12× wider |
| prime block, `q1q1` | 3× wider |
| prime block, `b3b3` | 48× wider |
| degree-3 determinant bound | 26% lower |
| 3×3 third minor `Δ3` | 73% lower |

The fix is to recover the factorization rather than to special-case it: repeated
exact synthetic division of the bivariate kernel by `(a − L)` in `a`, at the root
`a = L`, until the remainder no longer vanishes. `basis_algebra.kernel_factored`
returns the multiplicity and the exact quotient coefficients, and `kernel_value`
evaluates the quotient by Horner and multiplies by `(L − a)^m` once. That
restores every radius bit-for-bit, and the degree-3 determinant enclosure back to
`[1.077648488215e-06, 1.804795874074e-06]` — identical to the pre-refactor run.

The multiplicities are pinned as literals in `tests/test_kernel_algebra.py`, and
one test measures the width ratio on a ball directly, so the factorization
silently ceasing to be *used* fails a gate rather than quietly costing 73% of a
bound.

The same tests record why the multiplicity claim is restricted to same-parity
pairs: a cross-parity kernel is the **zero polynomial** in `(a, L)`, exactly —
the parity block structure appearing at the kernel level, and the reason there
is no `(L − a)` factor to find in it.

The general lesson is worth keeping: *a generalization that is mathematically
correct can still be numerically much worse, and interval arithmetic is where
that shows up.* Reproducing a table's values is not the same as reproducing its
tightness.

## 5. The independent cross-check

A cross-check that shares code is a re-run. `src/independent_even3.py` imports
none of the modules it checks — a test asserts that by parsing its imports — and
reaches the same formula by a different route: SymPy integrates the overlap
kernels straight from their definition where the rigorous path uses exact
bivariate integer arithmetic, and the pole Laplace transforms are evaluated by
quadrature where the rigorous path uses a closed form.

Worst relative difference across five points and all six entries: **1.56e-13**.

It is E3 and promotes nothing. mpmath never certifies in this program.

## 6. What the channels actually said

This is the first block where the four channels *could* disagree, because it is
the first where the determinant does not fix the spectrum. A positive
determinant on a 3×3 is consistent with `(3,0,0)` and with `(1,2,0)` — two
negative eigenvalues multiply to a positive contribution — so det-only reporting
would not have distinguished a positive definite block from one with a
two-dimensional negative subspace. That ambiguity does not exist at 2×2.

| channel | outcome |
|---|---|
| positivity | certified, and no longer a complete characterization |
| inertia | `(3,0,0)`; would have survived a negative outcome, which positivity would not |
| moments `m₁..m₄` | **INSUFFICIENT_INFORMATION** — they constrain the inertia and do not force it |
| rank–trace | `rank ≥ 1` against a true rank of 3 |

The moment result is the substantive change from ENG-006. At `n = 2` the map
from a spectrum to `(m₁, m₂)` is injective, so the moments recovered the inertia
exactly and the channel looked stronger than it is. At `n = 3` that map is not
injective. The ENG-006 finding was an artefact of the dimension; this one is the
general behaviour.

Rank–trace got *weaker* with dimension — `1` out of `2` at degree 3, `1` out of
`3` here. The inequality is tight at projections and this block's eigenvalues are
nowhere near 1 (its trace is around `0.1` across the whole cell). Recorded as a
weak result rather than tuned until it looked better.

## 7. Formal layer

Five theorems added, no `sorry`, axioms exactly `{propext, Classical.choice,
Quot.sound}`; the manifest now carries fifteen.

* `pd_three_by_three_certificate` — positive certified lower bounds on the three
  leading minors imply positive definiteness;
* `diagonal_congruence_preserves_pd` — a nonzero diagonal preconditioner
  preserves definiteness;
* `preconditioned_certificate3` — the two composed, which is the implication the
  runtime actually performs;
* `diagonal_congruence_preserves_index` and `_rank` — the signature, not just the
  yes/no, because the runtime reads a signature and would have needed these had
  the block been indefinite.

The 3×3 certificate states its payload as bounds on the minors rather than as
entry enclosures. At 2×2 the worst corner of a determinant is findable in closed
form, so that certificate carries enclosures; at 3×3 it is not, and encoding
interval arithmetic in Lean would be the wrong division of labour. §WO-RH-53 is
explicit that Lean verifies the implication, not the Arb arithmetic.

The numeric warrant stays E1 and the implication warrant is FORMAL, in separate
fields, as ENG-007 established.

## 8. A baseline defect found and fixed

Running every gate on the baseline before touching anything — step 1 of the
runbook — found that `ci_inertia.py --gate fast` had been **red since PR #11
merged**. It asserted the PIR content-kind set against a frozen literal of the
four ENG-006 kinds, and ENG-007 added two more without updating it. Separately,
`WEIL_DEGREE3_POSITIVITY_CERTIFICATE` was being emitted and promoted while
appearing in no declared list at all.

The fix is a single registry (`src/content_kinds.py`) with an explicit
`psd_licensable` answer per kind, and a gate that checks a *property* rather than
an equality: every kind PIR publishes must be declared, and the predicate that
actually decides must return the declared answer.

That property check immediately found a second, latent hole.
`satisfies_psd_requirement` refused inertia kinds and trusted every other kind
that said `psd_claim`, so a kind invented later was licensed by omission — a
`FORMAL_THEOREM_CERTIFICATE` claiming PSD would have passed, and it proves an
implication and asserts no number. The predicate is now default-deny.

## 9. Certified numbers

Read from the certificates by `scripts/check_docs.py`, not transcribed.

| quantity | bound (raw) | bound (preconditioned) |
|---|---|---|
| `Δ1` | `≥ 0.07537591825740127` | `≥ 1.2060146921184203` |
| `Δ2` | `≥ 3.335516179674528e-06` | `≥ 0.21859638835114986` |
| `Δ3` | `≥ 2.4352136989119354e-14` | `≥ 0.00041836652782391813` |

The `Δ3` bound is conservative: pointwise determinants across the cell run from
`1.4e-11` to `6.5e-11`, so the certified lower bound sits well below the true
minimum. That gap is the interval widening at the binding box, not the block's
own scale, and it is the clearest measure of what conditioning still costs at
this dimension.

## 10. For ENG-009

`certificates/eng009_structural_diagnostics.json` compares all five certified
blocks and records three candidate invariants, each with the falsifier that would
kill it. **No infinite-dimensional theorem is inferred, suggested or implied.**
Every row is a finite block on one cell under one normalization, and the patterns
across them are observations about five certificates.
