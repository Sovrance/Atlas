---
status: CURRENT
work_order: ATLAS-RH-ENG-010
supersedes: docs/history/agent-instructions-initial-integration.md
---

# RH / Weil — live agent instructions

**Read this before touching `math/rh_weil/`.** It describes the repository as it
is now, not as it was at integration. The original integration work order is
preserved at
[`docs/history/agent-instructions-initial-integration.md`](docs/history/agent-instructions-initial-integration.md)
and is **historical**:

<!-- docs-check: superseded-quote start -->
it says WO-RH-05 is open and degree 3 must not start, and both of those are
false today.
<!-- docs-check: superseded-quote end -->

**No RH proof claim is made** by this program. Every artifact is scoped to
`finite_dimensional_weil_compression`; a finite block certificate is evidence for
that block, interval, normalization and cutoff, and for nothing else.

Machine-readable status is
[`certificates/work_order_status.json`](certificates/work_order_status.json).
When this file and that file disagree, that file wins and this one is a bug —
`scripts/check_docs.py` fails CI on exactly that disagreement.

## Non-negotiable rules

These have not changed and are not up for renegotiation by a future work order
that finds them inconvenient.

1. Never emit or imply `RH PROVED` from finite polynomial blocks. Every artifact
   carries `rh_proof_claim: false` and `claim_scope:
   finite_dimensional_weil_compression`.
2. The sign convention is `G = G0 - Gp + Ginf` under the frozen Candidate-A
   normalization (WO-RH-17). The rejected `(sqrt(3)/2)` calibration is archival
   in `src/rejected_pole.py` and production may not import it.
3. `uncertified != pass`. Imported or historical numbers stay pending until
   regenerated locally.
4. Exact identities (E0), interval certificates (E1) and floating scans (E3) are
   separate evidence classes. A scan never promotes an E1 claim.
5. Every E1 claim records interval, precision, algorithm, source hashes and
   outward-rounded bounds. A certificate whose source hashes are stale is
   refused by the promotion predicate, not merely flagged.
6. Spectral/quantum material is diagnostic only.
7. Do not modify unrelated Atlas certificates or benchmark outputs.

## What is already done

| Work order | State |
|---|---|
| WO-RH-01…04 | exact identities, f1 audit, even block — **E0** |
| WO-RH-05, 10–15 | **recovered by ENG-005** (were quarantined by WO-RH-17) |
| WO-RH-08 | **done by ENG-006** — odd degree-3 implemented and certified |
| WO-RH-17/18 | normalization adjudicated, Candidate A adopted |
| WO-RH-28…36 | inertia, rank–trace, moments, degree-3 pilot — ENG-006 |
| WO-RH-37…46 | **done by ENG-007** — formal boundary and documentation truth pass |
| WO-RH-47…55 | **current work** — the 3×3 even block `{1, b, b²}`, certified positive definite |

See the README's *Current certified results* table for the numbers, and
`work_order_status.json` for the authoritative per-order state.

## How to run things

```bash
python3 scripts/run_rh_weil_suite.py             # fast path — does NOT re-derive E1
python3 scripts/run_rigorous_chain.py --release  # the real chain, in canonical order
python3 scripts/ci_inertia.py --gate fast        # exact gates, no python-flint needed
python3 scripts/ci_inertia.py --gate rigorous    # interval gates, python-flint required
python3 scripts/check_docs.py                    # documentation truth gate
```

Passing the fast suite does **not** mean the E1 certificates are current. Read
the rigorous chain's exit code before believing any E1 claim is fresh.

## Working rules for an agent

* **Regenerate, don't hand-edit.** Certificates are build outputs. If you change
  a file listed in a certifier's `DEPENDENCIES`, that certificate is stale by
  construction and the chain will refuse it until you re-run the certifier.
* **A stop condition is a result.** If an interval cover cannot separate, if two
  independent implementations disagree, or if a hypothesis cannot be discharged,
  report it. Do not widen a tolerance until the red goes away.
* **INCONCLUSIVE is a valid output.** So is a weak bound. Both are preferable to
  a number whose derivation you cannot defend.
* **Keep the evidence classes apart.** The most common way this program could go
  wrong is a floating diagnostic quietly becoming a warrant.

## Current frontier

ENG-008: the 3×3 even Weil block `G[{1, b, b²}]` is certified **positive
definite** on `[log 3, log 4]`, inertia `(3, 0, 0)`, by two independent routes.
This is the first block in the program where the determinant does not fix the
spectrum, and the first where the four channels give genuinely different answers:
the moments no longer force the inertia, rank–trace got weaker rather than
stronger, and conditioning became necessary. See
[`docs/HIGHER_DIMENSIONAL_BLOCK_ENG008_v0.1.md`](docs/HIGHER_DIMENSIONAL_BLOCK_ENG008_v0.1.md).

Two working rules this block added, both worth keeping:

* **Derive, do not tabulate.** The overlap kernels and the `L`-derivative
  machinery are now computed from the basis coefficients. Both replaced
  hand-written per-element tables that raised `KeyError` for anything new, and
  both were verified to reproduce every retired entry exactly *and* to agree with
  independent symbolic integration. Reproducing a table is not the same as being
  right; check both.
* **A preconditioner is a claim about a congruence.** If you rescale, the
  rescaling has to be exactly invertible, exactly applied, and covered by a
  theorem that says the inertia is unchanged. Powers of two give all three for
  free; a general Jacobi scaling gives none of them.

ENG-009 answered what those margins measure: the determinant collapse is
mostly the coordinate system (the healthy exact reference metric loses seven
of the thirteen orders by itself), and the congruence-invariant generalized
gap `λmin(G, M)` — certified for all five blocks against the exact `L²` Gram
metric — is the primary margin observable from now on. Raw determinant
magnitude is never quoted as distance-to-failure. See
[`docs/GENERALIZED_GAP_ENG009_v0.1.md`](docs/GENERALIZED_GAP_ENG009_v0.1.md).

One working rule ENG-009 added:

* **A margin is quoted against a named metric or not at all.** Anything that
  moves by `det(S)²` under a change of basis is a statement about the
  coordinates. The pencil `(G, M)` with the exact `L²` reference metric is the
  invariant object, its shifted-positivity certificates are eigensolver-free,
  and the implication is proved in Lean (`generalized_rayleigh`,
  `generalized_pencil_congruence`).

ENG-010 certified that block — positive definite, inertia `(4, 0, 0)`, two
routes — and enclosed its gap: `[1.907e-06, 2.416e-06]`, bottleneck at the
right cell edge, twenty times below the E3 midpoint scout. The preregistered
models were adjudicated `NEITHER_FALSIFIED` against the pinned artifact,
before any refit. Two working rules this order added:

* **A scout is a place to look, never a number to quote.** The E3 preview's
  midpoint gap was off by 20× because the bottleneck sits at the cell edge.
  Certified covers are the only source of quoted numbers, and the docs gate
  holds every quoted number to a certificate path.
* **Preregister, pin, adjudicate, then refit.** The adjudication script
  verifies the ENG-009 artifact's content hash and refuses to run against
  anything else; the verdict lands in its own certificate before any model
  is refitted. The tests pin the hash too.

Next, ENG-011 targets the even 5×5 `{1, b, b², b³, b⁴}` (refits diverge
11.9× at n = 5, `certificates/eng011_target_selection.json`), plus the
structural question of the edge bottleneck. No E1 work on n = 5 is launched
before ENG-010 is interpreted; no infinite-dimensional limit is inferred from
any of it.
