---
status: CURRENT
work_order: ATLAS-RH-ENG-007
latest_completed: ATLAS-RH-ENG-006
---

# RH / Weil — live agent instructions

**Read this before touching anything under `math/rh_weil/`.** This file is the live
entrypoint. The original integration work order is preserved, and superseded, at
[`docs/history/agent-instructions-initial-integration.md`](docs/history/agent-instructions-initial-integration.md);
do not take instructions from it.

## Claim boundary

This program studies **finite-dimensional** Weil-form compressions. It does **not** prove the
Riemann Hypothesis and no artifact here may imply that it does. `rh_proof_claim` is `false`
in every certificate, and it stays false. Finite positivity on a cell is a statement about
that cell.

## Non-negotiable rules

1. Never emit or imply `RH PROVED` from finite polynomial blocks.
2. Keep the sign convention `G = G0 - Gp + Ginf` under an explicit normalization audit.
   The adopted pole is Candidate A; the `(sqrt(3)/2)` outer product was **rejected** by
   WO-RH-17 and is archival only in `src/rejected_pole.py`. Production must not import it.
3. `uncertified != pass`. Imported or historical numbers stay pending until regenerated
   in-repo under the active normalization.
4. Exact identities (E0), rigorous interval certificates (E1), and floating scans (E2/E3)
   are different evidence classes. A missing rigorous dependency must **fail** a job, never
   silently degrade the artifact to a weaker class.
5. Every E1 claim records interval, precision, algorithm/version, source hashes, and
   outward-rounded bounds.
6. Spectral/quantum material is diagnostic only.
7. An **inertia** certificate is not a **positivity** certificate. A consumer that requires
   PSD is not satisfied by a signature, and ENG-006 keeps these as distinct content kinds
   deliberately.
8. Do not delete contrary evidence. Supersede it, label it, and link it.

## Current state (ENG-006 merged)

Everything below is certified in-repo; see the table in [`README.md`](README.md) and the
machine-readable [`certificates/work_order_status.json`](certificates/work_order_status.json).

- Normalization adjudicated and frozen (WO-RH-17); Candidate A active.
- Scalar cell `[log 3, log 4]`: uniform rigorous lower bound, E1, PROMOTED.
- Degree-1 odd and compact degree-2 even blocks: E1, PROMOTED.
- T=84 degree-2, point and uniform: E1, PROMOTED.
- Inertia, rank–trace and spectral-moment engines: operational.
- Odd degree-3 block: E1, positive, inertia `(2,0,0)`.

**WO-RH-05 is closed.** **Degree 3 is done.** If you find a document telling you otherwise,
it is stale — fix it, and check why `scripts/check_docs.py` did not catch it.

## Current work order: ENG-007

Formalize the stable finite theorem boundary in Lean and keep documentation mechanically
truthful. See [`formal/README.md`](formal/README.md).

The formal layer proves **implications**, not numerics. Arb produces enclosures; Lean proves
those enclosures imply the advertised finite conclusion. A formal theorem never converts an
E1 interval result into a FORMAL one.

## Canonical commands

```bash
python3 math/rh_weil/scripts/run_rh_weil_suite.py          # fast suite + status regen
python3 math/rh_weil/scripts/run_rigorous_chain.py --release   # full E1 chain
python3 math/rh_weil/scripts/ci_inertia.py --gate fast     # no python-flint needed
python3 math/rh_weil/scripts/ci_inertia.py --gate rigorous # python-flint required
python3 scripts/check_docs.py                              # documentation truth gate
python3 scripts/check_formal_manifest.py                   # formal manifest + axiom audit
cd math/rh_weil/formal && lake build                       # AtlasRH + comparator
```

## Before you commit

- `scripts/check_docs.py` passes. Stale live documentation is a correctness defect: the next
  agent will act on it.
- `scripts/check_formal_manifest.py` passes — no `sorry`, no project axioms, no statement
  drift.
- Existing rigorous certificates still validate; a changed certified number needs a stated
  dependency change, not a shrug.
