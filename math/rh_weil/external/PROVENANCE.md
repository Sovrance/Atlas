# Connes-CvS provenance (WO-CVS-01)


**No RH proof claim is made** by this layer or by anything it feeds. The
external oracle is a cross-check on shared ingredients; it certifies nothing.

## Reviewed package

| Field | Value |
|---|---|
| Package | `connes-cvs==0.3.1` |
| License | MIT |
| Repository | https://github.com/akivag613/connes-cvs- |
| Tip SHA at integration | `8ce0fc791ed9c9ca6f4ba512322720b4be80421b` |
| Role | Optional external cross-validation oracle (not canonical Atlas proof engine) |

## Local research environment (integration host)

Recorded while executing WO-CVS-01 on the Atlas agent host. Re-run
`python math/rh_weil/scripts/run_connes_cvs_crosschecks.py` to refresh the
certificate's live dependency block.

| Dependency | Version / note |
|---|---|
| Python | 3.14.4 |
| mpmath | 1.4.1 |
| python-flint | 0.9.0 |
| connes-cvs | 0.3.1 |

## Upstream fast tests

Wheel installs do not ship tests. From upstream source tip matching the
integration SHA window:

```text
pytest tests/test_validation.py tests/test_operator_hardening.py
→ 98 passed
```

## Install (research / cross-check env only)

```bash
pip install 'connes-cvs==0.3.1' python-flint mpmath
```

Do not add these to Atlas's required runtime dependencies. Core exact-identity
tests remain stdlib-only.
