# Integration Manifest

## Added paths

- `math/rh_weil/README.md`
- `math/rh_weil/AGENT_INSTRUCTIONS.md`
- `math/rh_weil/INTEGRATION_MANIFEST.md`
- `math/rh_weil/SHA256SUMS.txt`
- `math/rh_weil/theory/FORMULAS.md`
- `math/rh_weil/src/core.py`
- `math/rh_weil/src/mpmath_core.py`
- `math/rh_weil/src/__init__.py`
- `math/rh_weil/tests/test_exact_identities.py`
- `math/rh_weil/tests/test_connes_cvs_adapter_contract.py`
- `math/rh_weil/tests/test_connes_cvs_crosschecks.py`
- `math/rh_weil/certificates/imported_notebook_state.json`
- `math/rh_weil/certificates/external/connes_cvs_crossvalidation_v0.1.json`
- `math/rh_weil/notebook/RH_RESEARCH_NOTEBOOK_V2_INTEGRATION.md`
- `math/rh_weil/external/` (Connes–CvS adapter, crosschecks XC-01..06, provenance)
- `math/rh_weil/scripts/run_connes_cvs_crosschecks.py`
- `docs/rh-weil-integration-v0.1.md`

## Existing files intentionally not changed

No existing Atlas certificate, PIR schema, canonical document, CI runner, or benchmark implementation was overwritten. Integration is additive so a remote coding agent can review and promote the module in controlled stages.

## Local validation performed during packaging

`python math/rh_weil/tests/test_exact_identities.py` -> 5 tests PASS.
`python math/rh_weil/tests/test_connes_cvs_adapter_contract.py` -> 2 tests PASS (stdlib; optional `connes-cvs` not required).
`python math/rh_weil/tests/test_connes_cvs_crosschecks.py` -> PASS (skips XC suite if oracle absent; with oracle: XC-01..03 at dps 40 and 80).
`python math/rh_weil/scripts/run_connes_cvs_crosschecks.py` -> ACCEPTANCE GATE PASSED; external certificate written with `rh_proof_claim: false`.

Upstream fast tests (source tip): `pytest tests/test_validation.py tests/test_operator_hardening.py` -> 98 passed.

Exact/algebraic tests are not interval certificates. External XC results are E3 diagnostics (XC-06 external finite-matrix E1 only).
