# The 5×5 Even Block and the Migrating Bottleneck (ATLAS-RH-ENG-011)

**Scope.** Finite-dimensional Weil compression only. No RH proof claim is made
anywhere in this document or by any artifact it describes; `rh_proof_claim:
false`.

ENG-011 certified the 5×5 even block, explained the right-edge bottleneck
ENG-010 found, decomposed the new `b⁴` direction's Schur channel, and ran the
second preregistered adjudication — which returned the one verdict nobody
preregisters hoping for, and the honest one.

## 1. The certified results

> The 5×5 even Weil block `G[{1, b, b², b³, b⁴}]` is **positive definite** for
> every `L ∈ [log 3, log 4]`, inertia `(5, 0, 0)`, and its generalized gap
> against the exact L² reference metric satisfies
> **`λmin(G, M) ∈ [1.9033905118703842e-07, 2.5221165652583946e-07]`** —
> with the bottleneck **returned to the interior**, certified inside
> `[1.11, 1.30]` (the E3 scout locates it near `L ≈ 1.173`).

Both heavy warrants are parallel interval LDL* sweeps of the frozen-dyadic
preconditioned matrices (exponents `(−2, −6, −10, −10, −13)`, spanning ~3×10⁶
of raw diagonal ratio exactly):

* the **unshifted sweep** stratifies the cell — one stratum, signature
  `(5, 0, 0)`, no transition regions (the inertia certificate);
* the **shifted sweep** proves `G − λM ≻ 0` at the exact dyadic
  `λ = 1.9033905118703842e-07` on every box — which is simultaneously the
  uniform gap lower bound (`generalized_rayleigh`), the second positivity
  route (`shifted_positivity_transfer`: `G = (G − λM) + λM` with `M ≻ 0` by
  the E0 metric certificate), and, restricted by `nested_gap_regression`,
  a re-proof of every smaller block's bound at this `λ`.

The bottleneck classification is rigorous, not scouted: side sweeps prove
`λ* ≥ 4.0e-07` on `[log 3, 1.11]` and `[1.30, log 4]`, endpoint certificates
prove `λ*(log 3) ≥ 6.4e-07` and `λ*(log 4) ≥ 1.28e-06`, and the certified
Rayleigh witness at `L = 1.173` sits below all three — so the infimum over the
compact cell is attained strictly inside `[1.11, 1.30]`.

## 2. The bottleneck migrated — and the n=4 edge was real (§WO-RH-78)

ENG-010 found the n = 4 bottleneck at the right cell edge. Verdict now:
**`INTERIOR_MINIMUM_RETURNS`**. Two certified facts, not one grid picture:

* for the frozen n = 4 bottleneck witness, `d/dL [vᵀ(G₄ − λ₄M₄)v]` is
  **certified negative** at every sampled point of `[1.30, log 4]` — the
  n = 4 quotient genuinely falls into the edge, with the pole term driving
  the decrease and the prime term partially cancelling it;
* at n = 5 the infimum provably moved inside `[1.11, 1.30]` (§1).

So the right edge was not structural: it was a property of the n = 4 pencil's
weakest direction. The observation that `q = 4` joins the prime set exactly at
`L = log 4` is recorded as an observation — the neighboring cell stays out of
scope (§Anti-overclaim).

## 3. The Schur channel of `b⁴` (§WO-RH-79)

Writing the pencil in nested form `[[A₄(λ), c], [cᵀ, d]]` and bounding
`S₅ = d − cᵀA₄⁻¹c` by certified interval solves (Arb's verified LU — never a
floating inversion): the `b⁴` column's coupling to `span{1, b, b², b³}`
consumes **98.4–99.8%** of its preconditioned diagonal across the cell, with
`S₅` smallest (`≈ 9.4e-04` preconditioned) at `L ≈ 1.173` — exactly where the
pencil bottleneck sits. The n = 5 gap loss is therefore **introduced by the
new direction**, not inherited: the leading 4×4 pencil keeps its certified
ENG-010 gap (`shifted_shift_monotone` licenses it at the smaller `λ₅` with no
new cover). The semantic payload is formal: `schur_witness_block` proves the
bordered form positive from a witness solve and a positive residual.

## 4. The adjudication: `TOLERANCE_TOO_WIDE` (§WO-RH-83)

The frozen post-ENG-010 refits (artifact content hash pinned before any n = 5
E1 number existed) predicted:

| model | predicted `λmin(5)` | ×5 window | enclosure vs window |
|---|---|---|---|
| A: exponential refit | `5.657e-08` | `[1.131e-08, 2.828e-07]` | inside |
| B: power-law refit | `6.713e-07` | `[1.343e-07, 3.356e-06]` | inside |

**Verdict: `TOLERANCE_TOO_WIDE`** — for the *second consecutive dimension*, a
tight certified enclosure (width ~25%) lands in the overlap of both windows,
between the two point predictions. Two decisive numbers in a row have now
failed to separate the models only because ×5 windows overlap exactly where
the results keep falling. The verdict says that, rather than pretending the
models were tested harder than they were. Recorded before any refit, per the
order.

The certified even-family gap ladder now:

| n | `λmin(G, M)` enclosure | bottleneck |
|---|---|---|
| 1 | `[5.267e-02, 5.852e-02]` | — |
| 2 | `[5.867e-04, 6.520e-04]` | — |
| 3 | `[3.606e-05, 4.008e-05]` | — |
| 4 | `[1.907e-06, 2.416e-06]` | right edge `log 4` |
| 5 | `[1.903e-07, 2.522e-07]` | interior `≈ 1.173` |

Successive decay ratios ~89, 16.3, 17.6, 9.4 — still no one-parameter law.

## 5. Nesting is now a theorem, and it held (§WO-RH-84)

`M₄` and `G₄` are exactly the leading blocks of `M₅` and `G₅`, so generalized
Cauchy interlacing demands `λmin(5) ≤ λmin(4)`. The certified intervals
confirm it strictly (`2.522e-07 < 1.907e-06`), and the direction used is now
formal: `nested_gap_regression` proves that a shifted-PSD certificate at n = 5
restricts along any index embedding — certified gap lower bounds regress
upward through nesting, the formal regression bound §WO-RH-84 asked for.

The rest of the information comparison repeats the pattern more strongly:
the raw determinant fell to `~1e-29` pointwise (ten more orders of coordinate
collapse; the invariant gap moved by one), rank–trace weakens again, moments
constrain neither the inertia nor the pencil, and the generalized gap alone
expressed the bottleneck's migration.

## 6. Lean (§WO-RH-85)

Four new comparator-audited theorems — `shifted_positivity_transfer`,
`shifted_shift_monotone`, `nested_gap_regression`, `schur_witness_block` —
bring the manifest to **25 theorems**, no `sorry`, three standard axioms. The
5×5 PD implication itself rides the generic transcript/congruence theorems
(the sanctioned "generic theorem" route): nothing about dimension 5 needed a
new criterion.

## 7. ENG-012 (§WO-RH-87)

Selection: **adjudication reform + bottleneck dynamics**. Two consecutive
`TOLERANCE_TOO_WIDE`-class outcomes mean the n = 6 test must be preregistered
at a defensibly tighter tolerance (chosen by a stated power analysis, so the
two refits' diverging predictions cannot both contain one tight enclosure) —
and the bottleneck location `L*(n)` is promoted to a first-class certified
observable with its own model, since its migration (edge at n = 4, interior at
n = 5) is now the most structured unexplained phenomenon in the family. The
n = 6 element `b⁵` derives from the primitive table on demand; **no n = 6 E1
work is started here**.

## 8. Anti-overclaim, restated

Five positive finite blocks do not imply RH. A decreasing gap does not prove
convergence to zero. The `log 4` observation is not globally special until
neighboring cells are studied. Finite nested-subspace interlacing is not the
finite-to-infinite bridge. `rh_proof_claim: false`.
