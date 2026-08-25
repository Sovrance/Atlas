# Spectral Scaling Laws and Margin Diagnostics (ATLAS-RH-ENG-009)

**Scope.** Finite-dimensional Weil compression only. No RH proof claim is made
anywhere in this document or by any artifact it describes; `rh_proof_claim:
false`.
Baseline: Atlas main after the ENG-008 merge (`d7c47b7`), plus the regenerated
certificate set.

ENG-008 ended with a certified positive definite 3×3 block whose raw third
minor is bounded below by `6.451586222238981e-15` — next to entries of order
`1e-1`. The question this work order answers is whether that number measures
anything: *is the shrinking margin an intrinsic spectral gap closing, or an
artifact of the coordinates the block happens to be written in?*

## 1. The answer

**Mostly coordinates — but not only.** Three numbers, all certified, carry the
whole verdict:

| quantity, scalar block → 3×3 even block | orders of magnitude lost |
|---|---|
| raw determinant lower bound | ~13.0 |
| determinant of the exact reference metric `M` over the same span | ~6.9 |
| congruence-invariant generalized gap `λmin(G, M)` | ~3.2 |

The raw determinant's collapse is spectacular, but seven of its thirteen orders
are already present in the determinant of the *perfectly healthy* reference
metric — a Hankel moment matrix, positive definite for every `L > 0`, whose
minors on the even basis are exactly `1`, `1/180`, `1/7938000`. Determinants of
Gram matrices of increasingly similar functions collapse; that is what
determinants do, and it has nothing to do with the spectrum approaching
indefiniteness.

What remains after the coordinates are quotiented out is the generalized gap,
and it tells a different story: a decay of ~3 orders across the ladder, slow
and so far unclassified. Both halves are the result. The determinant's
collapse is not distance-to-failure; the gap's decay is real and is exactly
what ENG-010 is set up to measure next.

The dataset records this as
`MOSTLY_COORDINATE_DRIVEN_BUT_THE_GAP_ALSO_DECAYS`, and §Anti-overclaim is
enforced in both directions: no claim that the spectrum is collapsing, and no
claim that it is stable.

## 2. The reference metric (§WO-RH-58)

The pencil needs a positive reference form on the same basis, and the choice
is the ordinary `L²` Gram matrix on the support interval:

    M_ij(L) = ∫₀^L h_i(x; L) h_j(x; L) dx.

Three facts, all derived exactly in `src/reference_metric.py` rather than
asserted:

1. **Exact.** Every entry is an exact polynomial in `L` with rational
   coefficients, derived from the same primitive table
   (`basis_algebra.BASIS_L_POLY`) the kernels come from.
2. **A monomial.** Every basis element is homogeneous in `(x, L)` — degrees
   0, 1, 2, 3, 4, (6 for the prepared `bcube`) — so
   `M_ij(L) = m_ij·L^(dᵢ+dⱼ+1)` with `m_ij` a single rational. Equivalently
   `M(L) = D(L)ᵀ M(1) D(L)` with `D(L) = diag(L^(dᵢ+1/2))`, invertible for
   every `L > 0`.
3. **E0-positive.** By that congruence, `M(L)` is PD for all `L > 0` iff the
   constant rational matrix `M(1)` is — and that is decided by exact rational
   Sylvester minors. No interval arithmetic, no floating point, no covers.

`e0_eng009_reference_metric.json` is, as a result, the cheapest certificate in
the program: its warrant is exact arithmetic end to end.

## 3. The generalized gap (§WO-RH-58)

Everything rigorous reduces to *shifted positivity*: if `G − λM ⪰ 0` then
`vᵀGv ≥ λ·vᵀMv` for every `v`, i.e. every generalized eigenvalue of the pencil
`(G, M)` is at least `λ` — stated and certified without any eigensolver. The
scout that proposes `λ` is float Sylvester bisection (E3, recorded as such);
the warrant is adaptive interval covers of the leading minors of the exactly
preconditioned shifted block, exactly as ENG-008 certified positivity itself.
Upper bounds are certified Rayleigh quotients of rational witness vectors at
the scouted bottleneck.

Certified enclosures for `inf_L λmin(G, M)(L)` on `[log 3, log 4]`:

| block | dim | raw det ≥ | `λmin(G, M)` enclosure |
|---|---|---|---|
| scalar `{1}` | 1 | `6.96e-2` | `[0.052666425704956055, 0.058518285644676]` |
| odd `{q1}` | 1 | `1.50e-2` | `[0.12558043003082275, 0.1395338743236466]` |
| even `{1, b}` | 2 | `2.07e-6` | `[0.0005867481231689453, 0.0006520075621795992]` |
| odd `{q1, b3}` | 2 | `1.07e-6` | `[0.004970192909240723, 0.005522444641828629]` |
| even `{1, b, b²}` | 3 | `6.45e-15` | `[3.606081008911133e-05, 4.007732435284577e-05]` |

Two structural observations the raw column could never have shown:

* the odd family is *healthier* than the even family at equal dimension — its
  n = 2 gap is an order of magnitude larger, while the raw determinants of the
  two n = 2 blocks are indistinguishable;
* the even family's gap decays by factors of ~90 then ~16 — decelerating,
  where the raw determinant accelerates (×3·10⁴ then ×3·10⁸).

The enclosure widths are dominated by the deliberate 10% slack between the
scouted crossing and the certified shift, not by interval arithmetic.

## 4. Invariance, proved not narrated (§WO-RH-57/63)

The claim that licenses cross-basis comparison is
`det(SᵀGS − λ·SᵀMS) = det(S)²·det(G − λM)`: the pencil's roots do not move
when both forms are transformed together. This is checked as an exact
polynomial identity in `λ` (all coefficients) in `tests/test_generalized_gap.py`,
and proved in Lean:

* `AtlasRH.rayleigh_lower_of_shifted_psd` — shifted PSD ⟹ the Rayleigh bound;
* `AtlasRH.congruence_sub_smul` and `shifted_posDef_congruence` /
  `shifted_posDef_of_congruence` — the certified gap transports across any
  invertible simultaneous congruence, in both directions, with the same `λ`;
* `AtlasRH.gap_of_preconditioned_certificate3` — the composition the runtime
  actually performs: preconditioned 3×3 minor bounds ⟹ the generalized
  Rayleigh bound for the original pencil.

All three (with their statement-layer twins) enter the comparator and the
manifest: **18 theorems**, no `sorry`, the three standard axioms.

By contrast, *raw* eigenvalues are not invariant — the same test that proves
the pencil identity exhibits the trace moving under the same congruence. The
program's rule going forward: **raw determinant magnitude is never quoted as a
distance-to-failure observable**; a margin is quoted against a named metric or
not at all.

## 5. Scaling models and their falsifiers (§WO-RH-60)

`e3_eng009_scaling_models.json` fits exponential and power-law decay per
parity family — E3, `EXPLORATORY_NEVER_PROMOTED`, over the certified
enclosures' midpoints. With three even points the two models genuinely
diverge at n = 4:

| model (even family) | predicted `λmin` at n = 4 |
|---|---|
| exponential `C·ρⁿ` | `~7.8e-7` |
| power law `C·n^(−p)` | `~6.1e-6` |

The float preview of the prepared 4×4 block scouts its midpoint gap at
`~4.5e-5` — *above both predictions*, hinting the gap may be decelerating
faster than either model. That is E3 twice over and decides nothing; it makes
the ENG-010 measurement genuinely discriminating, which is the property
§WO-RH-62 selects for. Every model record carries the explicit falsifier
interval its rejection requires.

## 6. The next block (§WO-RH-62)

**Selected: `even4 = {1, b, b², b³}`.** The even family is where the collapse
phenomenon lives, where the models disagree, and where a fourth point falsifies
at least one of them. The odd-3 candidate is recorded with its scores and the
reason it lost; the mixed 5×5 is rejected as a repackaging (parity makes it
block-diagonal); T=84 is rejected because finite-T and cutoff-free families
are not one scaling sequence.

Preparation done here (E0/E3 only, per the order):

* `bcube = b³ = x³(L−x)³` added to the primitive table — kernels, endpoint
  polynomials, derivative machinery and reference metric all *derive*; the
  generic exact tests (including direct symbolic integration) sweep it, and
  its kernel multiplicities `(L−a)^m` with `m = 4, 5, 6, 7` are pinned;
* `e3_eng010_even4_preview.json` — float conditioning (raw minors down to
  `6.4e-19`), proposed frozen preconditioner exponents `[2, 6, 10, 11]`, and
  the scouted pencil gap.

No E1 certificate for the 4×4 exists, deliberately.

## 7. The channel comparison (§WO-RH-61)

The dataset scores seven channels on the seven criteria the order names. The
selected primary diagnostic vector:

    (inertia, generalized_λmin_lower vs L²-Gram, trace/n, m2/n,
     conditioning, moments m1..m4, rank_trace_bound)

The generalized gap earns its place as the only channel that is simultaneously
basis-invariant, graded rather than boolean, certifiable without an
eigensolver, meaningful after positivity fails, comparable across dimension,
and already formalized. Inertia stays first because it is the claim the others
qualify.

## 8. Anti-overclaim, restated

* Raw determinant shrinkage is not intrinsic spectral collapse — and raw
  determinant magnitude is not a basis-invariant distance-to-RH signal.
* Five finite blocks do not establish asymptotics; the E3 fits are falsifiable
  plans, not facts.
* Every generalized gap names its reference metric
  (`l2_gram_on_support`, `e0_eng009_reference_metric.json`).
* Finite-T and cutoff-free families are never pooled.
* Any RH/off-line-zero conclusion still requires a proved finite-to-infinite
  bridge — see `docs/BRIDGE_CANDIDATES_ENG009_v0.1.md`, which is conjectural
  and says so.
* `rh_proof_claim: false`, everywhere.
