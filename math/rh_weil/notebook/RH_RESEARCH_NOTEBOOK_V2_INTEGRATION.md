---
status: CURRENT
work_order: ATLAS-RH-ENG-007
---

# RH Research Notebook V2 → Atlas: integration history and current state

**Claim boundary: finite-dimensional Weil-form research. No RH proof claim.**

This document has two halves, and the split is the point. The first records what the
imported notebook claimed, as history. The second records what Atlas has actually certified
in-repo. Conflating them is how a research input gets mistaken for a warrant.

---

## Part 1 — Historical notebook checkpoint (evidence, not instructions)

> **HISTORICAL.** Everything in Part 1 describes the state at import time. None of it was
> promoted by the integration commit, and several of its "next targets" are now closed.
> It is retained because WO-RH-17 forbids deleting contrary evidence.

### Accepted research trajectory (as imported)

1. Scalar f0 verifier: prime-power cell splitting, strict convexity above log(plastic
   constant), downward derivative jumps at breakpoints.
2. f1 audit: midpoint-odd basis `q1 = x − L/2`; exact overlap formula; active prime shifts
   strengthen the odd channel under `G = G0 − Gp + Ginf`.
3. Even degree-2 block: bubble basis `b = x(L−x)`; exact prime kernels; midpoint parity
   factorization.
4. Structural bounds: endpoint-jet filtration gives omitted-tail scaling
   `O(log(T)/T^(2r+1))` when `r` endpoint jets vanish.
5. Fourier cross-check: `T = 84` selected; stable entire low-frequency transforms and
   support-length Taylor jets developed.

### What the notebook reported

A hybrid uniform `T=84` degree-2 certificate, and later independent direct-Fourier point
certificates at `log 3`, an interior bottleneck near `L ≈ 1.10595`, and `log 4`. Its
remaining gap was interval sign coverage of `E2,84''` and `E2,84'`.

**Atlas policy applied at import:** none of those numerical claims was promoted. They were
research inputs pending in-repo regeneration. `certificates/imported_notebook_state.json`
still carries `IMPORTED_PENDING_REGENERATION` and
`hard_constraints_certified: false`, and the suite asserts that on every run.

### One notebook reading that did not survive

The notebook's "monotone `E2' > 0`" reading was a **Candidate-B artifact** and is not
reused. ENG-005 rescanned the `T=84` topology fresh under Candidate A: `E2'` changes sign
once, at `L* ≈ 1.10595`, and the minimum is interior. This is a worked example of why
imported numbers stay pending until regenerated.

---

## Part 2 — Current Atlas state (ENG-006 merged)

Authoritative sources: [`../README.md`](../README.md) and the generated
[`../certificates/work_order_status.json`](../certificates/work_order_status.json).

| Object | Domain | Result | Warrant |
|---|---|---|---|
| scalar | `[log 3, log 4]` | uniform positive lower bound | E1 |
| degree-1 odd | cell | positive | E1 |
| degree-2 even (compact) | cell | positive | E1 |
| `T=84` degree-2 | cell | positive, point + uniform | E1 |
| degree-3 odd | cell | inertia `(2,0,0)`, positive | E1 |
| rank–trace | degree-3 sample points | rank lower bound, nontrivial | E1 |
| spectral moments | degree-3 points | `m1..m4`, mixed conclusive/insufficient | E1 |

The direct-Fourier uniform check Part 1 lists as missing was **closed by ENG-005**, by two
independent warrants: an exhaustive interval cover assuming no topology, and an
interior-minimum argument locating `L*`. The headline bound is floored at what the cover
alone proves.

Degree 3 is **implemented and certified**, not future work. The odd block `{q1, b3}` with
`b3 = x(L−x)(x − L/2)` carries a rigorous inertia result of `(2,0,0)` on the current cell.

## Part 3 — Current frontier

ENG-007: the formal theorem boundary in Lean — see [`../formal/README.md`](../formal/README.md).
The next mathematical target is ENG-008, the first genuinely `>2`-dimensional parity block,
prepared under WO-RH-46. The degree-3 odd block is still `2×2`, so determinant and trace
already encode most of its spectral behaviour; inertia and moments only start earning their
keep at `3×3`.
