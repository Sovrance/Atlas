---
status: HISTORICAL
---

> **HISTORICAL.** Record of the original notebook->Atlas import. Paths and contents have moved since;
> the live inventory is the tree itself and `certificates/work_order_status.json`.

# Integration Manifest

## Added paths

- `math/rh_weil/README.md`
- `math/rh_weil/AGENT_INSTRUCTIONS.md`
- `math/rh_weil/INTEGRATION_MANIFEST.md`
- `math/rh_weil/SHA256SUMS.txt`
- `math/rh_weil/theory/FORMULAS.md`
- `math/rh_weil/src/core.py`
- `math/rh_weil/src/cells.py`
- `math/rh_weil/src/scalar.py`
- `math/rh_weil/src/fourier.py`
- `math/rh_weil/src/certificate_io.py`
- `math/rh_weil/src/mpmath_core.py`
- `math/rh_weil/src/__init__.py`
- `math/rh_weil/tests/` (exact, scalar, Fourier, Connes contract/XC)
- `math/rh_weil/certificates/` (E0 regenerated + E3 Fourier scan; imported pending unchanged)
- `math/rh_weil/external/` (optional Connes–CvS oracle)
- `math/rh_weil/scripts/run_rh_weil_suite.py`
- `math/rh_weil/scripts/run_connes_cvs_crosschecks.py`
- `docs/rh-weil-integration-v0.1.md`

## Existing files intentionally not changed

No existing Atlas root certificates, PIR schemas, or `ci/run_all_certified.py` were
expanded. RH has a dedicated runner only (WO-RH-07).

## Work-order status

| Order | Status |
|---|---|
| WO-RH-01 | done (exact identities) |
| WO-RH-02 | done (E0 scalar cell [log3,log4]) |
| WO-RH-03 | done (f1 / midpoint-odd) |
| WO-RH-04 | done (bubble algebraic block) |
| WO-RH-05 | partial — stable H0/Hb + L-jets + E3 probe scan; **interval E1 coverage still open** |
| WO-RH-06 | partial — E0 certs regenerated; imported notebook not promoted |
| WO-RH-07 | done — `scripts/run_rh_weil_suite.py` |
| WO-RH-08 | blocked until WO-RH-05 E1 closes |

## Local validation

```bash
python math/rh_weil/scripts/run_rh_weil_suite.py
```

Exact/algebraic tests are stdlib. Fourier forms and Connes XC need `mpmath`
(and optionally `connes-cvs`). No RH proof claim is emitted.
