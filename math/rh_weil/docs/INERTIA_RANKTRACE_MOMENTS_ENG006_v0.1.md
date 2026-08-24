# Weil Inertia, Rank–Trace, Spectral Moments, and the Degree-3 Pilot (ATLAS-RH-ENG-006)

**Baseline:** ENG-005, merged (`f9e52e0`). **Scope:** build certified inertia,
rank–trace and spectral-moment channels, then run the odd degree-3 Weil block
through them as the first live workload.

**No RH proof claim.** Claim scope is `finite_dimensional_weil_compression`.
Everything below concerns one cell, `L ∈ [log 3, log 4]`.

---

## 1. Why stop treating positivity as the only observable

Through ENG-005 the program asked one question of every block: is it positive?
That is a single bit, and when the answer is no the run ends with nothing. The
three channels here are all attempts to have something left in that case:

* **inertia** — the full signature `(n₊, n₋, n₀)`, so an indefinite block still
  yields a rigorous result;
* **rank–trace** — a finite-dimensional lower bound on rank from traces and a
  Hilbert–Schmidt norm;
* **spectral moments** — `m₁..m₄` and what they do and do not pin down.

The degree-3 block was expected to be the case that needed them. It turned out
not to be — §5 below — but the channels are what established that, and they are
what makes the negative case survivable next time.

---

## 2. The inertia engine (§3/§4)

`inertia/` certifies `Inertia(A) = (n₊, n₋, n₀)` by congruence reduction.
Sylvester's law is what makes it rigorous: symmetric elimination and symmetric
permutation are both congruences, so the signature of the pivots *is* the
signature of `A`.

**Soundness on intervals.** Run the elimination on a matrix of balls and record
the pivot order. For any point matrix inside the box, the same pivot order gives
point pivots inside the corresponding pivot balls, because interval arithmetic
encloses the point computation step by step. So if every pivot ball has a
determined sign, every point matrix in the box has that signature.

**Two pivot kinds, not one.** Diagonal pivots alone cannot handle `[[0,1],[1,0]]`
— inertia `(1,1,0)`, no usable diagonal entry anywhere. When no diagonal has a
determined sign the engine looks for a symmetric 2×2 sub-block with *definitely
negative* determinant, which contributes exactly `(1,1)` whatever its diagonal
does. Without it, ordinary indefinite matrices would report `INCONCLUSIVE` for no
better reason than pivot order.

**Exact zero is not numerical zero.** The interval path never reports a nonzero
`n₀`. A ball containing zero is not a proof that the entry is zero, and a
zero-radius ball produced by arithmetic is not one either. Singular and
near-singular inputs come back `INCONCLUSIVE`; exact zero multiplicity is
reported only on the exact rational path, where a vanishing pivot with a
vanishing row is structural.

**An independent oracle.** The engine is checked against inertia read off the
characteristic polynomial by Descartes' rule of signs — exact, because a
symmetric matrix is real-rooted, which is precisely the condition that turns
Descartes from a bound into an equality. It shares no code with the elimination.
They agree on every case and across 50 exact rational congruences per signature.

**Stratification.** `certify_inertia_family` subdivides a parameter cell until
each piece carries a determined signature, merges neighbours that agree, and
reports the leftovers as `INCONCLUSIVE_TRANSITION_REGION`. That leftover is the
honest part: a signature can only change where the matrix is singular, and an
interval enclosure can never certify singularity, so the cells bracketing a
crossing will never resolve however far they are split. Their width is the
resolution to which the transition has been located. Strata and transitions tile
the cell exactly.

---

## 3. Rank–trace (§5)

```
rank(P) ≥ 2 tr(P) + 4 tr(Q) − 4b − ‖P+Q‖²_HS
```

under `P ⪰ 0`, `Q` Hermitian with at most `b` positive directions, all terms in
the theorem's normalization.

The hypotheses are objects, not comments. An inequality is only as good as the
conditions attached to it, and those are exactly what gets lost when a bound is
carried between settings, so the engine refuses to produce a number unless every
hypothesis arrives with a verification status. There is no partial credit: one
unverified hypothesis makes the certificate `INCONCLUSIVE` with no number to
quote. The theorem id and normalization tag are compared, not assumed.

**The normalization is not cosmetic.** The inequality is not scale-free — rank is
scale-invariant, `tr(P)` is degree 1 and the HS term degree 2 — so it holds only
under a normalization. Equality at a projection (`tr P = ‖P‖²_HS = r`, giving
`rhs = r = rank P`) identifies that normalization as *spectrum in `[0,1]`*, which
is checkable, so it is checked rather than declared.

A right-hand side that comes out non-positive is reported as `trivial`: true, but
saying nothing, since rank is a non-negative integer anyway. A null result must
not be able to masquerade as a finding.

---

## 4. Spectral moments and the B1 adapter (§6)

`m_k = tr(G^k)` are enclosed as traces of matrix powers — never via an eigenvalue
solver, which §14.4 forbids from supporting an E1 claim.

The eigenvalue counting measure `μ = Σ δ_{λᵢ}` has `m₀ = n`, so asking what
spectra are consistent with given moments *is* a truncated moment problem, and
Atlas already has an exact engine for it in `b1_moment_solver`. The exact path
calls B1's Hankel/PSD/rank/flatness routines directly. B1 is exact-rational and
the moments of a Weil block are enclosures, so interval inputs go to the ENG-006
inertia engine for the same definiteness questions. Neither is a second solver.

**The asymmetry that matters.** A truncated localizing matrix failing to be PSD
is conclusive — no representing measure on `[0,∞)`, so a negative eigenvalue
exists. A truncated localizing matrix *being* PSD proves nothing on its own:
sufficiency needs a flat extension (Curto–Fialkow) that four moments do not
supply. So "the moments force PSD" comes back `INSUFFICIENT_INFORMATION` even for
a matrix that is in fact positive definite. `INSUFFICIENT_INFORMATION` is a
certified outcome, not a failure to compute.

Two-sided eigenvalue bounds come from Wolkowicz–Styan on `m₁, m₂`. At `n = 2`
these collapse to equalities, so the first two moments determine the spectrum
outright — and the odd degree-3 block is 2×2.

---

## 5. The degree-3 block (§7/§8/§9/§10)

Midpoint parity splits the basis: `1` and `b = x(L−x)` are even about `x = L/2`,
`q1 = x − L/2` and `b3 = x(L−x)(x−L/2)` are odd. The Gram is block diagonal and
the odd block is the 2×2 `[[G_q1q1, G_q1b3], [G_q1b3, G_b3b3]]`.

Both kernels are re-derived from the basis with SymPy on every run:

```
K_q1b3 = (L−a)² (L³ + 2L²a − 12La² − 6a³) / 60
K_b3b3 = (L−a)³ (L⁴ + 3L³a − 15L²a² − 18La³ − 6a⁴) / 420
```

A transcription error in either would have been invisible: the assembled block
would still be a smooth positive-looking 2×2 and every certificate built on it
would be internally consistent and wrong.

**The convention had to be pinned.** The repository's `K_ij` is the *symmetrized*
overlap, covering both the `+log q` and `−log q` shifts the explicit formula sums
over. The one-sided integral gives exactly half — a silent factor-2 error through
the whole prime block. The test re-derives ENG-005's `K_q1q1` as well, whose value
is already fixed by merged work, so the convention is checked against something
that cannot move.

**Prime-shift indefiniteness, preserved.** Every active shift block on the cell
has negative determinant — inertia `(1,1,0)`. §7 forbids termwise PSD domination
and this is why: the assembled entry is positive while every individual term of
`Gp` is indefinite. It is kept as a regression test rather than worked around.

**Cross-validation.** Real-space against the frequency route at `L = 1.2`: the
fast-decaying `(q1,b3)` and `(b3,b3)` agree to ~11 digits with the truncation tail
shrinking ~45000× from `T = 84` to `T = 4000`; slow-decaying `(q1,q1)` differs by
the expected tail. Same pattern ENG-005 documented, so §16's disagreement stop
condition does not fire.

### The result (§9, Outcome A)

The block is **positive definite on the whole closed cell**:

```
Inertia(G_odd,3)(L) = (2, 0, 0)   for every L in [log 3, log 4]
```

by interval LDL congruence with adaptive subdivision — one stratum, no
transition regions, 11024 boxes at maximum subdivision depth 1. §9's two named
quantities, certified by exhaustive cover:

| quantity | certified uniform lower bound | cover |
|---|---|---|
| `O1 = G[q1,q1]` | `1.5331267702267608e-02` | 8192 boxes, depth 0 |
| `det(G_odd,3)` | `1.073120529992708e-06` | 32768 boxes, depth 0 |

The determinant bound is deliberately sized. A cover accepts a box the moment
its lower end clears the target, so a cover sized only to decide positivity
reports whatever the worst box happened to give — at 8192 boxes that was
`9.6e-11`, four orders below the grid minimum of `1.407e-6`. True, but a weak
certificate for the quantity it bounds, which is the point ENG-005 made about
`E2`. Halving the box radius halves the enclosure width, so the reported bound
uses a finer cover and lands within a factor of 1.3 of the truth.

**This is not the outcome the work order was built for.** §0 says the goal is
"not prove degree 3 positive" and §17 counts an inertia stratification as
success. The block turned out definite anyway. The stratification machinery is
still what proved it, and would have delivered Outcome B unchanged had any part
of the cell been indefinite.

### Positivity is a different content kind from inertia (§11)

§11 requires that an inertia certificate never satisfy a consumer demanding PSD.
That rule binds the *inertia* artifact, and it does: the nested stratification
here has signature `(2,0,0)` and still refuses, because "I know the signature"
must not be silently read as "it is positive".

The outer artifact is a different claim — the block is positive definite, with
certified positive lower bounds on `O1` and `det` — so it carries a positivity
content kind and answers a PSD consumer honestly. Both are checked: the same
signature is accepted on the positivity certificate and refused on the inertia
object nested inside it.

---

## 6. What each channel was worth (§10)

Five questions, answered from the emitted certificates rather than from prose,
so the report cannot drift from the artifacts it describes.

**Did positivity alone characterize the block?** Yes — but only because it is
2×2, where `(trace > 0, det > 0)` *is* the signature. There is no gap between
"positive definite" and "inertia known" at this size. That is a fact about the
dimension, not about the method.

**What did inertia add?** Beyond positivity, nothing here. What it delivered was
the route: it proved the signature constant over the closed cell rather than at
sample points, and it is what would have produced a usable result had the block
been indefinite. A channel that adds nothing when the easy answer holds is
behaving correctly; it earns its place on the blocks where the easy answer does
not.

**Did the moments recover the inertia?** Yes, exactly, and again because `n = 2`.
The Wolkowicz–Styan inequalities become equalities there, so `m₁` and `m₂` invert
to the spectrum outright and the signature is a consequence of two moments rather
than a separate observation. `m₃` and `m₄` were not needed. This is the sharpest
the moment channel can ever be and it does not survive to larger blocks.

**Were the moments insufficient?** Also yes, by the general truncated-moment
route — and the two answers are not in conflict, because they are different
questions. Asking whether the moments *force* a non-negative spectrum via the
localizing matrix returns `INSUFFICIENT_INFORMATION`, correctly: only the
refuting direction is available from four moments. The `n = 2` inversion answers
an easier question and is what settled the inertia.

**Did rank–trace give a nontrivial bound?** Yes but weak: `rank ≥ 1` against a
true rank of 2. The right-hand side is positive, so it is not the vacuous
`rank ≥ 0`, but it is off by a factor of two — and structurally, not through
tuning. The inequality is tight at projections, where every eigenvalue is 0 or 1;
this block's eigenvalues run from `3.3e-05` to `4.2e-02`, four orders below 1, so
each contributes about `2λ` to the right-hand side instead of the 1 it would
contribute at a projection. A Weil Gram block is not the kind of operator this
inequality is sharp for. Recorded as a null-ish result per §10 rather than tuned
until it looked better.

---

## 7. Limitations (stated, not hidden)

* One cell, `[log 3, log 4]`, and one 2×2 block. Nothing here says anything about
  other cells, higher degree, or RH.
* Being 2×2 does a lot of work. The inertia follows from trace and determinant
  with no elimination, the eigenvalues have a closed form, and `m₁, m₂` invert to
  the spectrum. None of that survives to larger blocks, where the moment channel
  degrades to bounds and the inertia engine has to earn its pivots.
* The rank–trace bound is weak here, and weak for a structural reason: the
  inequality is tight at projections and this block's eigenvalues are four orders
  of magnitude below 1.
* `INSUFFICIENT_INFORMATION` from the localizing-matrix route is a statement
  about what four moments can decide, not about the block.
* The inertia engine's interval path can never certify a zero eigenvalue. A
  genuinely singular family would be reported as a transition region, not as
  `n₀ = 1`.
