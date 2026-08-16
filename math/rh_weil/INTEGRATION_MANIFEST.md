# Integration Manifest

## Added paths

- `math/rh_weil/README.md`
- `math/rh_weil/AGENT_INSTRUCTIONS.md`
- `math/rh_weil/INTEGRATION_MANIFEST.md`
- `math/rh_weil/SHA256SUMS.txt`
- `math/rh_weil/theory/FORMULAS.md`
- `math/rh_weil/src/core.py`
- `math/rh_weil/src/__init__.py`
- `math/rh_weil/tests/test_exact_identities.py`
- `math/rh_weil/certificates/imported_notebook_state.json`
- `math/rh_weil/notebook/RH_RESEARCH_NOTEBOOK_V2_INTEGRATION.md`
- `docs/rh-weil-integration-v0.1.md`

## Existing files intentionally not changed

No existing Atlas certificate, PIR schema, canonical document, CI runner, or benchmark implementation was overwritten. Integration is additive so a remote coding agent can review and promote the module in controlled stages.

## Local validation performed during packaging

`python math/rh_weil/tests/test_exact_identities.py` -> 5 tests PASS.

These tests validate algebraic implementation consistency only. They are not interval certificates.
