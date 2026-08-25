---
status: CURRENT
work_order: ATLAS-RH-ENG-007
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
| WO-RH-37…46 | **current work** — formal boundary and documentation truth pass |

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

ENG-007: make the finite theorem boundary formally replayable in Lean, and make
the documentation mechanically unable to fall behind the implementation. See
[`docs/FORMAL_BOUNDARY_ENG007_v0.1.md`](docs/FORMAL_BOUNDARY_ENG007_v0.1.md).

After that, ENG-008 attacks the first genuinely >2-dimensional parity block,
where inertia and moments add information that a 2x2 determinant does not.
