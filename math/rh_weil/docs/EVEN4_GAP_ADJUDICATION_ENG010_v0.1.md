# The 4×4 Even Block and the Model Adjudication (ATLAS-RH-ENG-010)

**Scope.** Finite-dimensional Weil compression only. No RH proof claim is made
anywhere in this document or by any artifact it describes; `rh_proof_claim:
false`.

ENG-010 is the first work order in this program designed as a **prediction
test**. ENG-009 preregistered two exploratory decay models for the even-family
generalized gap, disagreeing by ~8× about dimension 4; this order certified
the block and adjudicated them — against the preregistered artifact, before
any refit, with the artifact's content hash pinned as a tripwire.

## 1. The certified results

> The 4×4 even Weil block `G[{1, b, b², b³}]` is **positive definite** for
> every `L ∈ [log 3, log 4]`, inertia `(4, 0, 0)`, and its generalized gap
> against the exact L² reference metric satisfies
> **`λmin(G, M) ∈ [1.90734863e-06, 2.415977410246834e-06]`**,
> with the bottleneck at the right cell edge `L → log 4`.

Definiteness is certified twice, by routes sharing the assembly and nothing
after it:

| route | result | cost |
|---|---|---|
| interval LDL* congruence, stratified | `(4, 0, 0)`, one stratum, no transition regions | 26,854 boxes, depth 6 |
| Sylvester, four adaptive covers | `Δ1 ≥ 0.0753`, `Δ2 ≥ 1.4363e-07`, `Δ3 ≥ 6.3384e-15`, `Δ4 ≥ 5.2788569975917046e-24` (raw) | up to 82,380 boxes, depth 5 |

The gap's lower bound is a shifted-positivity certificate: all four leading
minors of the frozen dyadic congruence `D(G − λM)D` positive over the cell at
the exact dyadic `λ = 1.90734863e-06` (the fourth needed 92,658 boxes, depth
6). The upper bound is a certified Rayleigh quotient of a rational witness at
the bottleneck. The implication chain — minor bounds ⟹ shifted PD ⟹ the
generalized Rayleigh bound for the *original* pencil — is proved in Lean
(`preconditioned_gap_certificate4`), as is the 4×4 Sylvester criterion itself
via a division-free Jacobi completed-square identity.

## 2. The scout was wrong, and the order predicted it

ENG-009's E3 preview scouted the 4×4 gap at the cell **midpoint**: `4.5e-5`.
The certified bottleneck sits at the **right cell edge**, twenty times lower.
`λ*(L)` is far from flat — it rises to `9.2e-5` at `L ≈ 1.17` and collapses to
`2.4e-6` at `log 4`. §12's instruction — *do not present the E3 scout as the
certified n=4 result* — was exactly the right discipline: an E3 number quoted
as the answer would have been off by a factor of twenty and on the wrong side
of both model predictions.

Why the edge? That is now a named ENG-011 question. What ENG-010 records is
the fact: the pencil weakens toward `log 4`, where a new prime power (4 = 2²)
enters the prime block at the cell boundary.

## 3. The adjudication (§WO-RH-71)

Preregistered predictions (ENG-009, fitted on n = 1..3, tolerance ×5):

| model | predicted `λmin(4)` | falsifier window | certified enclosure | verdict |
|---|---|---|---|---|
| A: exponential `C·ρⁿ` | `7.494e-07` | `[1.499e-07, 3.747e-06]` | `[1.907e-06, 2.416e-06]` | **NOT FALSIFIED** |
| B: power law `C·n^(−p)` | `5.912e-06` | `[1.182e-06, 2.956e-05]` | same | **NOT FALSIFIED** |

**Verdict: `NEITHER_FALSIFIED`** — expected outcome E of the work order. The
certified enclosure is *tight* (width ~25%); the models survived because
their ×5 windows overlap precisely where the result fell, between the two
point predictions (geometric mean of `7.5e-7` and `5.9e-6` is `2.1e-6`; the
truth is `~2.4e-6`). The experiment was decisive about the *number* and
indecisive between these two *models* — a statement about the preregistered
tolerance, not about the block, and recorded as exactly that.

The verdict was recorded before any refit; the adjudication artifact carries
the preregistered file's content hash and its verification.

## 4. What n = 4 adds to the ENG-009 verdict

The determinant story continues on schedule: the raw `Δ4` uniform bound is
`5.28e-24` — five more orders of collapse beyond the 3×3 — while the
invariant gap fell by ~15× (`[3.61e-5, 4.01e-5] → [1.91e-6, 2.42e-6]`). The
reference metric's own determinant loses another seven orders at n = 4
(`det M(1) = 1/88104560544000`). Raw determinant magnitude remains a
statement about coordinates; the ENG-009 rule — a margin is quoted against a
named metric or not at all — stands.

The even-family gap ladder, all certified:

| n | basis | `λmin(G, M)` enclosure |
|---|---|---|
| 1 | `{1}` | `[5.267e-02, 5.852e-02]` |
| 2 | `{1, b}` | `[5.867e-04, 6.520e-04]` |
| 3 | `{1, b, b²}` | `[3.606e-05, 4.008e-05]` |
| 4 | `{1, b, b², b³}` | `[1.907e-06, 2.416e-06]` |

Successive decay factors ~90, ~16, ~16: no longer the deceleration the first
three points suggested, and not the sharp exponential either.

## 5. After the adjudication: refits and ENG-011 (§WO-RH-71/75)

Refitted on all four points (E3, `EXPLORATORY_NEVER_PROMOTED`), the two
models diverge by **11.9×** at n = 5: exponential predicts `5.66e-08`, power
law `6.71e-07`. The `even5` block `{1, b, b², b³, b⁴}` is therefore selected
for ENG-011 — a genuinely discriminating fourth test, with the additional
question of whether the bottleneck stays pinned at the right edge.
`b⁴` derives from the same primitive table when needed; **no E1 work on
n = 5 is started here**, per §WO-RH-75.

## 6. The information channels at n = 4 (§WO-RH-72)

The comparison report records: the moments still do not force the inertia
(already false at n = 3, worse here); the moments say nothing useful about
the pencil; rank–trace weakens again (`rank ≥ 1` against a true rank of 4);
and the generalized gap is more clearly the primary margin than at n = 3 —
it measured the real decay and located a bottleneck no raw quantity showed.

## 7. Anti-overclaim, restated

- n = 4 does not establish asymptotics; four points do not select a law.
- Surviving a ×5 falsifier window is not confirmation.
- A positive 4×4 block is finite-dimensional evidence only.
- The gap is always quoted against the exact L² reference metric
  (`e0_eng010_even4_reference_metric.json`).
- Finite-T remains a separate family.
- Any RH/off-line-zero conclusion still requires a proved finite-to-infinite
  bridge (`docs/BRIDGE_CANDIDATES_ENG009_v0.1.md`, still conjectural).
- `rh_proof_claim: false`.
