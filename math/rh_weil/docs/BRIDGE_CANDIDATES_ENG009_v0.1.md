# Finite-to-Infinite Bridge Candidates (ATLAS-RH-ENG-009 §WO-RH-64)

**Status: conjectural. Evidence class E3 throughout. Nothing in this note is
proved, promoted, or a warrant for anything.** Its purpose is the opposite of
a claim: to state exactly which theorem is *missing* for each route from the
finite certified results to an infinite statement, so that no future document
can blur the line by accident. No RH proof claim is made here or anywhere in
this program; `rh_proof_claim: false`.

The finite side that actually exists: for each certified block, uniform
positive definiteness on `[log 3, log 4]`, inertia `(n, 0, 0)`, and a
certified generalized gap `λmin(G, M) ≥ λ_n > 0` against the exact `L²` Gram
metric. Five blocks, dimensions 1–3, one cell.

For each candidate below: the finite quantity Atlas certifies, the infinite
object it would need to become, the exact missing theorem, what Atlas
establishes today, what it does not, and the falsifier that would kill the
route.

---

## B1. Uniform generalized gap

* **Finite quantity.** `λ_n = inf_L λmin(G_n(L), M_n(L))` for the certified
  ladder of nested even (and odd) blocks.
* **Conjectured infinite object.** `inf_n λ_n > 0` — a uniform coercivity
  constant for the Weil form restricted to the full polynomial test family,
  in the `L²` metric.
* **Missing theorem.** That `λ_n` is bounded below away from zero as
  `n → ∞`, plus the identification of the limit form as the Weil distribution
  applied to the closure of the test family.
* **Atlas establishes.** `λ_1..λ_3` (even), `λ_1..λ_2` (odd), certified;
  a decay of ~3 orders over the ladder — i.e. the data currently point *away*
  from an easy uniform bound.
* **Atlas does not.** Any statement about `n ≥ 4`; any limit object.
* **Falsifier.** Certified `λ_n → 0` at any fitted decaying rate (the E3
  models' own falsifiers apply); or a single certified negative eigenvalue
  anywhere in the family, which would end the positivity program itself.

## B2. Interlacing across nested subspaces

* **Finite quantity.** The generalized spectra of nested blocks
  `G_n ⊂ G_{n+1}` (same cell, same metric family).
* **Conjectured infinite object.** Cauchy-interlacing-type monotone control:
  `λmin` of the pencil is non-increasing in `n` and its decrement is governed
  by the new element's component orthogonal to the old span.
* **Missing theorem.** Generalized-eigenvalue interlacing for the pencil
  sequence with *varying* metric (`M_n` is the leading principal block of
  `M_{n+1}` — true here by construction, which makes standard interlacing
  applicable; the missing part is the quantitative decrement bound).
* **Atlas establishes.** The premise `M_n = (M_{n+1})_{[1..n]}` exactly, by
  the monomial form; monotone decrease of the certified even-family gaps,
  consistent with interlacing.
* **Atlas does not.** Any proved decrement bound; anything at `n = 4` until
  ENG-010 certifies it.
* **Falsifier.** A certified `λmin` *increase* under adding an element would
  contradict interlacing's premise being applied correctly (it cannot
  increase; observing an increase means an assembly error, and is a stop
  condition, not a discovery).

## B3. Coercivity of the limiting quadratic form

* **Finite quantity.** `vᵀG_n(L)v ≥ λ_n·vᵀM_n(L)v` for all `v`, uniformly in
  `L` on the cell.
* **Conjectured infinite object.** The Weil quadratic form is coercive
  (`≥ λ‖·‖²_{L²}`) on a dense subspace of the relevant test-function space.
* **Missing theorem.** Density/completeness of the polynomial family (see
  B6), plus a limit interchange: coercivity constants surviving the passage
  from the finite ladder to its closure — precisely the kind of interchange
  ENG-005 records as an analytic hypothesis, never machine-verified.
* **Atlas establishes.** The finite inequality with certified constants.
* **Atlas does not.** Any function-space statement; the cell `[log 3, log 4]`
  is one bounded window, not a support family exhausting anything.
* **Falsifier.** `λ_n → 0` (as B1); or the density premise failing (B6).

## B4. Schur-complement recurrences

* **Finite quantity.** The certified leading minors and the LDL* pivots the
  inertia engine produces per block.
* **Conjectured infinite object.** A recurrence for the pivot sequence as
  `n` grows (a Jacobi/orthogonal-polynomial-type three-term structure) whose
  asymptotics decide the gap's limit.
* **Missing theorem.** That the Weil Gram family, in a suitable orthogonal
  basis for `M`, has bounded-coefficient recurrences — none of which is
  currently even conjectured precisely.
* **Atlas establishes.** The pivots at n ≤ 3, certified.
* **Atlas does not.** Any recurrence; the basis change that would exhibit one.
* **Falsifier.** Computed pivots at n = 4, 5 failing every bounded-coefficient
  ansatz.

## B5. Operator convergence

* **Finite quantity.** Compressions `P_n W P_n` of the Weil form to the
  n-dimensional test spaces.
* **Conjectured infinite object.** Strong (or norm-resolvent) convergence of
  the compressions to a self-adjoint operator whose positivity is RH-relevant
  (the Weil-positivity ⟺ RH equivalence supplies the last step *only* for
  the full test family, not a polynomial slice on one cell).
* **Missing theorem.** Convergence itself, and the identification of the
  limit; both are far outside anything Atlas computes.
* **Atlas establishes.** Nothing beyond the finite compressions.
* **Atlas does not.** Everything else; this is the widest gap in the list.
* **Falsifier.** Not falsifiable inside the program at all today — which is
  the honest reason it stays a note, not a work order.

## B6. Density / completeness of the polynomial test family

* **Finite quantity.** The basis ladder `{1, b, b², b³, …}` (even) and
  `{q1, b3, b²q1, …}` (odd) on `[0, L]`.
* **Conjectured infinite object.** Completeness of the family, as `L` ranges
  over supports and the degree grows, in the test space in which Weil
  positivity is equivalent to RH.
* **Missing theorem.** A density statement for compactly supported polynomial
  slices in that space, with quantitative approximation rates (the rates are
  what B3's interchange would consume).
* **Atlas establishes.** Exact linear independence (the metric's positive
  minors are exactly that) and the parity decomposition.
* **Atlas does not.** Density in any function space; any approximation rate.
* **Falsifier.** A test function provably outside every closure of the family
  whose Weil functional is negative while all polynomial slices stay positive
  — this would show the family's positivity is strictly weaker than RH.

---

## What this note is not

It is not a roadmap claiming RH is `k` theorems away. Several of the missing
theorems above are of unknown difficulty, B5 in particular. The single
concrete next step the finite program can take — measuring which decay law
the gap actually follows — is ENG-010's, and its outcome feeds B1/B2 whether
it confirms or falsifies the current fits.
