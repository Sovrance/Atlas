# External cross-validation certificates

Artifacts here are **not** Atlas E0/E1 RH certificates.

- Upstream floats → E3 diagnostic
- Upstream Arb residual on a supplied finite matrix → external E1, finite-matrix scope only
- `rh_proof_claim` must remain `false`

Regenerate:

```bash
python math/rh_weil/scripts/run_connes_cvs_crosschecks.py
```
