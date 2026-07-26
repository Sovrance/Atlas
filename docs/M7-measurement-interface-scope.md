# M7 — Measurement-Interface Layer: scope

**Status:** DONE (5/5). Implemented as `measurement_interface/` +
`tests/test_m7.py`, committed `certificates/m7_certificate.json`. The last
unbuilt *supporting machine* of roadmap-v2 (§0, §3 P3). M7 makes SPEC §2
**layer 3** — the measurement interface 𝖬 — executable, and is the gate the
roadmap requires *before any new cosmology / biology / BEC claim*.

## Why M7, and what it generalizes

The decoding model is `y = 𝖬[𝒰(x)] + ε`: measured data `y` is the unknown rule
system `𝒰` seen through an apparatus/preprocessing map `𝖬` plus noise `ε`
(Maudlin R14, `docs/notes/measurement-interface-maudlin.md`). Layers 1–2 (M1
canonicalization, M2 generation, the B-suite verifiers) all live *above* 𝖬 and
silently assume it. M7 is the machine that models 𝖬 itself and decides when a
structural conclusion is actually a statement about the apparatus.

It generalizes two existing, one-off results into a reusable engine:

- **B9-T3** — the "independent calibration route" requirement: a recovered
  quantity is only trustworthy if a *second, independent* calibration route
  recovers it within tolerance. M7 turns this into a general
  **injection–recovery + dual-route** primitive.
- **The S4 `APPARATUS_LIMITED` lesson** — a rank/threshold verdict that changes
  within the detector's declared resolution/bandwidth (rank plateaus at cutoffs)
  is not a structural fact. M7 turns this into a general **apparatus-limited
  gate** that sweeps the declared resolution window and refuses a structural
  verdict that is not stable across it.

## Contract (SPEC-locked)

- **Operationalizes SPEC §2 layer 3.** Every M7 result declares the apparatus
  transfer function, calibration route, noise family, and selection/censoring it
  assumes, and lists them as `m_layer_stipulations` (mandatory ≥ v0.3). M7
  certificates additionally carry a first-class **`calibration_route`** field.
- **Emits SPEC §4 verdicts with a stated cause.** Primary: `APPARATUS_LIMITED`
  (the verdict would change within the apparatus's declared
  resolution/bandwidth → no structural conclusion permitted). Also
  `OBSERVATIONALLY_EQUIVALENT` under a declared 𝖬 (§5) where two latents are
  indistinguishable after the interface, and a plain pass-through when the
  structural verdict is stable across the whole declared window.
- **Honesty (the M1-T5 discipline).** `APPARATUS_LIMITED` is *not* a synonym for
  "noisy": an ample apparatus must pass the structural verdict through
  unchanged. The refusal must fire on a genuine cutoff-dependence and *only*
  then.
- **No layer-1 promotion from M7 alone.** An interface model is layer-3 work; it
  can *demote* (to `APPARATUS_LIMITED`) or *license* a downstream claim, never by
  itself supply layer-1 evidence.

## Engine surface (`measurement_interface/engine.py`)

Stdlib + numpy, matching the M1/M2 machine convention (package = `__init__.py`
re-export + `engine.py`; results returned as plain dicts/tuples; no certificate
logic in the engine). Planned public API:

- `apply_transfer(signal, transfer)` — apply a declared apparatus transfer
  function (gain/bandwidth/affine) to a latent signal.
- `add_noise(y, noise, rng)` / noise-family helpers — the declared `ε`.
- `injection_recovery(latent, transfer, noise, estimator, ...)` → recovery dict
  with a certified interval (generalized B9-T3, single route).
- `dual_route_agreement(latent, route_a, route_b, tol)` → whether two independent
  calibration routes agree within tolerance (B9-T3 core), with the disagreement
  reported honestly when they don't.
- `censoring_correct(sample, selection, model)` — a selection/censoring model and
  the bias it removes (naive vs interface-aware estimator).
- `apparatus_limited_gate(verdict_fn, resolution_grid, declared_window)` → the
  central gate: sweep the structural verdict across the apparatus resolution
  window; return `PASS(verdict)` if stable, `APPARATUS_LIMITED(cause)` if it
  flips inside the declared window.

## Acceptance tests (`tests/test_m7.py`, T1–T5)

Framework-free script harness (matches `tests/test_m1.py`), writes
`certificates/m7_certificate.json` via `b1_moment_solver.certificate.save_certificate`.

- **T1 — injection–recovery (B9-T3 generalized).** Inject a known latent through a
  declared transfer + noise; recover it inside a certified interval. Positive.
- **T2 — independent calibration routes.** Two independent routes recover the same
  latent within tolerance → agree; a miscalibrated route disagrees and the
  disagreement is surfaced (calibration = an M-layer stipulation, made explicit).
- **T3 — `APPARATUS_LIMITED` trigger (S4 generalized).** A rank/threshold verdict
  stable above a detector cutoff but flipping at/below it → the gate returns
  `APPARATUS_LIMITED`, refusing the structural conclusion. The refusal fires.
- **T4 — censoring/selection.** A threshold-censored estimator is biased; the
  declared censoring model removes the bias (interface-aware ≫ naive), or, if
  uncorrectable at the declared resolution, the gate returns `APPARATUS_LIMITED`.
- **T5 — negative / no-collapse.** With an ample apparatus (verdict stable across
  the whole declared window) the gate does **not** return `APPARATUS_LIMITED`; it
  passes the structural verdict through. Guards the honesty invariant.

## Definition of done

- `measurement_interface/` machine + `tests/test_m7.py` green (T1–T5) +
  committed `certificates/m7_certificate.json` (class NUMERICAL-DISCOVERY, with
  `calibration_route` + `m_layer_stipulations`).
- README M-machines section + roadmap-v2 P3 updated (M7 done); SPEC §2/§4
  cross-reference to the executable layer-3 machine.
- Deferred to **M7-b** (recorded, not silently dropped): continuous/functional
  transfer families (deconvolution with regularization), full
  injection-recovery over real B4/B10 posteriors, and coupling M7 into the PIR
  `measurement_interface` fact field as an analyzer pass.
