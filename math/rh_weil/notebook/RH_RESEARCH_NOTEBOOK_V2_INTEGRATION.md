# RH Research Notebook V2 — integration checkpoint, and what became of it

Two things live in this file, kept apart on purpose. The first half is the
**historical** notebook checkpoint as it was imported: the research trajectory it
established, and the numerical claims it reported. The second half is the
**current** repository state, which supersedes those claims entirely — not
because they were wrong in spirit, but because none of them was regenerated
under the normalization that is now in force until an Atlas agent did it here.

**No RH proof claim is made** by the notebook, by this repository, or by
anything below. Everything concerns finite-dimensional polynomial compressions of
the localized Weil quadratic form on one cell.

---

## Part 1 — Historical: the imported notebook checkpoint

This is the state the notebook reached, recorded verbatim in substance. It is
provenance, not instruction.

### Research trajectory (accepted, and still the shape of the program)

1. **Scalar `f0` verifier** — prime-power cell splitting, strict convexity above
   `log` of the plastic constant, downward derivative jumps at breakpoints.
2. **`f1` audit** — midpoint-odd basis `q1 = x − L/2`; exact overlap formula;
   active prime shifts strengthen the odd channel under `G = G0 − Gp + Ginf`.
3. **Even degree-2 block** — bubble basis `b = x(L − x)`; exact prime kernels;
   midpoint parity factorization.
4. **Structural bounds** — endpoint-jet filtration gives omitted-tail scaling
   `O(log T / T^{2r+1})` when `r` endpoint jets vanish.
5. **Fourier cross-check** — `T = 84` selected; stable entire low-frequency
   transforms and support-length Taylor jets developed.

That trajectory survived contact with the repository. Every item above is still
the structure of the current implementation.

### Numerical claims the notebook reported

<!-- docs-check: superseded-quote start -->
The notebook reported a hybrid uniform `T = 84` degree-2 certificate, and later
independently reproduced direct-Fourier point certificates at `log 3`, an
interior bottleneck near `L = 1.10595`, and `log 4`. Its final direct-Fourier
work described the remaining gap as interval sign coverage of `E2,84''` and
`E2,84'`, and named as the next target: certify `E2,84'' > 0` on
`[log 3, 1.20]`, `E2,84' > 0` on `[1.20, log 4]`, and one interval point ball
near `L = 1.1059498113`, then infer uniform `E2,84(L) > 0`.

It also proposed the degree-3 basis `q1 = x − L/2`, `b3 = x(L−x)(x−L/2)`, and
instructed that degree 3 not be started until the independent degree-2 Fourier
certificate was reproducibly closed.
<!-- docs-check: superseded-quote end -->

**None of those numbers was promoted by the integration commit**, and none is
promoted now. They were research inputs, marked `IMPORTED_PENDING_REGENERATION`
until regenerated in-repo. That policy is unchanged.

---

## Part 2 — Current: what the repository actually holds

### The normalization the notebook used was later rejected

WO-RH-17 adjudicated the pole block. The even outer-product form
`(√3/2)(v₊v₊ᵀ + v₋v₋ᵀ)` was **rejected**: it equals the explicit-formula pole
times `(√3/2)cosh(L/2)`, a factor equal to 1 only at `L = log 3` — a calibration
fitted at a single test point, and +8.25 % off at `L = log 4`. Certificates that
depended on it were quarantined, preserved rather than deleted, and had to be
**regenerated** rather than reinterpreted.

This is why the notebook's numbers could not simply be adopted. A number computed
under a rejected normalization is not a number that needs re-checking; it is a
number about a different object.

### Atlas-native recovery, ENG-004 through ENG-006

| Notebook item | Repository state now | Warrant |
|---|---|---|
| scalar `f0` positivity | uniform rigorous lower bound on `[log 3, log 4]`, Arb | E1, promoted |
| degree-1 odd | cutoff-free uniform bound | E1, promoted |
| compact degree-2 even | cutoff-free uniform bound | E1, promoted |
| `T = 84` points | regenerated under Candidate A at three points | E1, promoted |
| `T = 84` topology | **rescanned from scratch.** The notebook's monotonicity reading was a Candidate-B artifact and is not reused: `E2'` changes sign exactly once, and the minimum is interior | E3 scan → E1 result |
| `T = 84` uniform | **closed.** Two independent warrants: an interval cover assuming no topology, and an interior-minimum argument that isolates `L*` to a nine-digit interval and certifies the derivative signs on either side | E1, promoted |
| degree-3 basis `{q1, b3}` | implemented, kernels E0-verified, block certified **positive definite**, inertia `(2,0,0)` uniformly on the cell | E1, promoted |

The exact numbers are in
[`../README.md`](../README.md) §3, read from the certificates by
`scripts/check_docs.py` rather than transcribed.

### Historical notebook evidence versus current repository certificates

The distinction is enforced, not merely stated:

* an imported claim carries `IMPORTED_PENDING_REGENERATION` and the promotion
  predicate refuses it;
* a regenerated certificate binds the active normalization id and the source
  hashes of everything it depends on, and is refused the moment either drifts;
* the notebook's historical minimum for the scalar entry, `0.0753795566…`, is
  recorded as *regression evidence only* — it falls inside the recovered
  enclosure, which is reassuring and is not the acceptance gate. The gate is the
  positive bound, never a fitted constant.

### Current research frontier

ENG-007 (current) formalizes the finite theorem boundary in Lean 4 / Mathlib and
makes the documentation mechanically unable to fall behind the implementation.
ENG-008 takes the 3×3 even parity block `{1, b, b²}` prepared by WO-RH-46, where
the determinant is no longer the whole story and inertia and moments have room to
add information a 2×2 does not have.

See [`../AGENT_INSTRUCTIONS.md`](../AGENT_INSTRUCTIONS.md) for live instructions
and [`../certificates/work_order_status.json`](../certificates/work_order_status.json)
for the authoritative per-order state.
