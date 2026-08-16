# RH / Weil Positivity Program (Research Notebook V2)

This directory integrates the RH Research Notebook V2 into Atlas as the reproducible source behind the Constant Atlas positivity-verifier example (`G00 > 0`, Schur pivot positivity).

## Scope and claim boundary

This program studies finite-dimensional polynomial compressions of the localized Weil quadratic form. **It does not prove the Riemann Hypothesis.** A finite block certificate is evidence only for that stated block, interval, normalization, and cutoff/representation.

Accepted assembly convention:

`G = G^0 - G^p + G^infinity`.

Current certified/research cell: `L in [log(3), log(4)]`.

## Evidence policy

- E0 / SOUND: exact algebraic identities independently re-derived/tested.
- E1 / SOUND: interval-certified numerical statements with reproducible interval engine and certificate.
- E3 / HEURISTIC: floating scans and localization of candidate minima.
- Spectral/quantum analogies: diagnostics only; never proof warrants.

Imported notebook claims are intentionally marked `IMPORTED_PENDING_REGENERATION` until an Atlas agent reproduces them from the code in this repository. Do not promote them merely because they appeared in the prior notebook transcript.

## Current priority

Reproduce the fully direct Fourier-side uniform degree-2 certificate at `T=84` by interval coverage of the centered-even determinant, then proceed to the midpoint-odd degree-3 block.

**Executed in-repo so far (see `certificates/work_order_status.json`):**
- WO-RH-01/03/04 — exact identities, f1 audit, bubble block (E0)
- WO-RH-02 — scalar cell `[log 3, log 4]` algebraic `W00''>0` (E0)
- WO-RH-05 — stable `H0`/`Hb` + L-jets; E3 energy-probe scan only (interval E1 **open**)
- WO-RH-06/07 — regenerated E0 certificates + dedicated runner (does not expand root CI)

Run:

```bash
python math/rh_weil/scripts/run_rh_weil_suite.py
```

See `AGENT_INSTRUCTIONS.md` and `notebook/RH_RESEARCH_NOTEBOOK_V2_INTEGRATION.md`.

Optional external oracle (not required at runtime): `external/` wraps
`connes-cvs` for independent cross-checks of shared Weil ingredients. See
`external/CONNES_CVS_MAPPING.md` and `external/PROVENANCE.md`. Do not compare
Galerkin matrix entries to Atlas polynomial Gram blocks by index.

Regenerate the external cross-validation certificate (research env):

```bash
pip install 'connes-cvs==0.3.1' python-flint mpmath
python math/rh_weil/scripts/run_connes_cvs_crosschecks.py
```

