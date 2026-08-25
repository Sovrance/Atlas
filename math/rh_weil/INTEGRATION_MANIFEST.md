# Integration Manifest

What the original integration commit added to Atlas, and what has been added
since. **No RH proof claim** is made by any of it; the claim scope throughout is
`finite_dimensional_weil_compression`.

For live status read [`AGENT_INSTRUCTIONS.md`](AGENT_INSTRUCTIONS.md) and
[`certificates/work_order_status.json`](certificates/work_order_status.json).
This file records *what exists*, not *what is true*.

## Original integration (WO-RH-01…07)

- `math/rh_weil/README.md`, `AGENT_INSTRUCTIONS.md`, `INTEGRATION_MANIFEST.md`, `SHA256SUMS.txt`
- `math/rh_weil/theory/FORMULAS.md`
- `math/rh_weil/src/` — `core.py`, `cells.py`, `scalar.py`, `fourier.py`, `certificate_io.py`, `mpmath_core.py`
- `math/rh_weil/tests/` — exact, scalar, Fourier, Connes contract/XC
- `math/rh_weil/certificates/` — E0 regenerated, E3 Fourier scan, imported claims left pending
- `math/rh_weil/external/` — the optional Connes–CvS oracle
- `math/rh_weil/scripts/run_rh_weil_suite.py`, `run_connes_cvs_crosschecks.py`
- `docs/rh-weil-integration-v0.1.md`

## Added since

**ENG-003/004 — normalization and the scalar canary.** `src/pole.py` as the
single Candidate-A pole primitive; `src/rejected_pole.py` as the archival home of
the rejected calibration; `src/promotion.py` and `src/normalization.py` for
binding and quarantine; `scripts/derive_normalization.py`,
`run_normalization_crosscheck.py`, `quarantine_normalization.py`,
`certify_scalar_canary.py`, `run_rigorous_scalar.py`.

**ENG-005 — core E1 recovery.** `src/archimedean_realspace.py` (the exact
real-space archimedean form), `src/interval_cover.py`, `src/interval_backend.py`,
`src/e1_cutoff_free.py`, `src/e1_t84.py`, `src/fourier_jets.py`,
`src/weil_fourier_jets.py`, `src/curvature_derivation.py`,
`src/rigorous_integration.py`; `scripts/certify_cutoff_free_e1.py`,
`certify_t84_e1.py`, `certify_degree1_e1.py`, `certify_degree2_compact_e1.py`,
`certify_fourier_T84_points.py`, `certify_fourier_T84_uniform.py`,
`run_rigorous_chain.py`.

**ENG-006 — inertia, rank–trace, moments.** The `inertia/`, `ranktrace/` and
`moments/` packages beside `src/`; `src/degree3.py`;
`scripts/certify_degree3.py`, `ci_inertia.py`, `report_information_comparison.py`.

**ENG-007 — the formal boundary and the documentation gate.** `formal/` (the
Lean 4 project: `AtlasRH/`, `comparator/`, `manifests/`); `src/formal_evidence.py`;
`src/pilot3.py`; `scripts/check_formal_manifest.py`, `certify_formal_boundary.py`,
`ci_formal.py`, `check_docs.py`, `preview_pilot3.py`; `docs_status.json`;
`external/zeta23/`.

## Existing files intentionally not changed

No Atlas root certificates, PIR schemas or `ci/run_all_certified.py` were
expanded. RH has dedicated runners only (WO-RH-07), now four gates:

```bash
python3 math/rh_weil/scripts/run_rh_weil_suite.py        # fast
python3 math/rh_weil/scripts/run_rigorous_chain.py       # rigorous
python3 math/rh_weil/scripts/ci_inertia.py --gate fast   # inertia
python3 math/rh_weil/scripts/ci_formal.py                # rh-formal
python3 math/rh_weil/scripts/check_docs.py               # rh-docs
```

Exact/algebraic tests are stdlib. Interval work needs `python-flint`; symbolic
derivations need `sympy`; both are required, not optional (see
`requirements-rigorous.txt`). The Lean layer needs a toolchain on `PATH` or
`ATLAS_LEAN_BIN` set, and degrades to its offline checks without one.
