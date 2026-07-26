"""M7 -- the measurement-interface layer ('the declared 𝖬' as executable code).

SPEC §2 layer 3 made runnable. The decoding model is ``y = 𝖬[𝒰(x)] + ε``:
measured data ``y`` is the unknown rule system ``𝒰`` seen through an apparatus /
preprocessing map ``𝖬`` plus noise ``ε`` (Maudlin R14,
``docs/notes/measurement-interface-maudlin.md``). M1/M2 and the B-suite verifiers
all live *above* 𝖬 and silently assume it; M7 models 𝖬 itself and decides when a
"structural" conclusion is really a statement about the apparatus.

It generalizes two one-off results into reusable primitives:

* **B9-T3** (``b9_circuit/estimate.py``, ``pir/domains/circuit_semantics.py``):
  a recovered quantity is only promotable if an *independent* calibration route
  confirms it. A self-consistent route (scale inferred from the same data) is
  fooled by a hidden bath and is non-promotable. Here: :func:`dual_route_recovery`
  with the same FDT-residual audit (threshold 0.15) and the spectroscopy/Gibbs
  distinction.
* **S4** (``s4_ds/kernel.py``): a rank verdict defined by an ε-cutoff plateaus at
  the detector resolution. Here: :func:`effective_rank` (mirrors the S4 kernel)
  and :func:`apparatus_limited_gate`, which sweeps the declared resolution window
  and returns SPEC §4 ``APPARATUS_LIMITED`` when the structural verdict is not
  stable across it.

stdlib + numpy only, matching the M1/M2 machine convention (engine returns plain
dicts; no certificate logic here).
"""

from __future__ import annotations

import math
from typing import Callable, Dict, List, Sequence, Tuple

import numpy as np

# The FDT-audit rejection threshold, carried over verbatim from B9
# (``b9_circuit/estimate.py``: fdt_resid < 0.15 == EQUILIBRIUM_CONSISTENT).
FDT_THRESHOLD = 0.15

# Calibration-route promotability, aligned with
# ``pir/domains/circuit_semantics.py``: an independent channel is promotable; a
# self-consistent (data-inferred) route is not.
PROMOTABLE_ROUTES: Dict[str, bool] = {
    "cal:spectroscopy_route": True,   # independent calibration channel
    "cal:gibbs_route": False,         # self-consistent; NOT promotable evidence
}


def is_promotable_route(route: str) -> bool:
    """Whether a named calibration route may supply promotable evidence."""
    if route not in PROMOTABLE_ROUTES:
        raise ValueError(f"unknown calibration route {route!r}")
    return PROMOTABLE_ROUTES[route]


# --------------------------------------------------------------------------- #
# 1. Apparatus transfer function + noise (the '𝖬[·] + ε' primitives)
# --------------------------------------------------------------------------- #

def apply_transfer(signal, transfer: Dict) -> np.ndarray:
    """Apply a declared affine apparatus transfer ``y = gain * x + offset``.

    ``transfer`` is the declared 𝖬: ``{"gain": g, "offset": b}`` (both floats or
    broadcastable arrays). Deliberately invertible so a *declared* interface can
    be undone in recovery; a non-invertible interface (gain 0) is an apparatus
    that erases the latent and is rejected here rather than silently divided by.
    """
    g = np.asarray(transfer.get("gain", 1.0), float)
    if np.any(g == 0.0):
        raise ValueError("degenerate transfer (gain 0): apparatus erases the latent")
    b = np.asarray(transfer.get("offset", 0.0), float)
    return g * np.asarray(signal, float) + b


def invert_transfer(measurement, transfer: Dict) -> np.ndarray:
    """Undo a declared affine transfer: ``x_hat = (y - offset) / gain``."""
    g = np.asarray(transfer.get("gain", 1.0), float)
    if np.any(g == 0.0):
        raise ValueError("degenerate transfer (gain 0): not invertible")
    b = np.asarray(transfer.get("offset", 0.0), float)
    return (np.asarray(measurement, float) - b) / g


# --------------------------------------------------------------------------- #
# 2. Injection–recovery with a certified interval (B9-T3, single route)
# --------------------------------------------------------------------------- #

def injection_recovery(latent: float, transfer: Dict, sigma: float,
                       n_samples: int = 400, n_trials: int = 400,
                       z: float = 1.96, seed: int = 0) -> Dict:
    """Inject a known latent through a declared 𝖬 + Gaussian ε, recover it.

    Each trial draws ``n_samples`` noisy measurements ``y = 𝖬[latent] + ε``,
    inverts the *declared* transfer, and forms a z-interval estimate of the
    latent. Returns the empirical coverage of the certified interval (should
    track the nominal ``2Φ(z)−1``) and the bias — the audit that a declared
    interface plus a stated noise family recovers ground truth without bias.
    """
    rng = np.random.default_rng(seed)
    nominal = 2.0 * 0.5 * (1.0 + math.erf(z / math.sqrt(2.0))) - 1.0
    covered = 0
    ests: List[float] = []
    widths: List[float] = []
    for _ in range(n_trials):
        clean = apply_transfer(np.full(n_samples, float(latent)), transfer)
        y = clean + rng.normal(0.0, sigma, n_samples)
        xhat = invert_transfer(y, transfer)
        est = float(xhat.mean())
        se = float(xhat.std(ddof=1) / math.sqrt(n_samples))
        lo, hi = est - z * se, est + z * se
        covered += int(lo <= latent <= hi)
        ests.append(est)
        widths.append(hi - lo)
    return {
        "latent": float(latent),
        "mean_estimate": float(np.mean(ests)),
        "bias": float(np.mean(ests) - latent),
        "coverage": covered / n_trials,
        "nominal_coverage": float(nominal),
        "mean_interval_width": float(np.mean(widths)),
        "n_trials": n_trials,
        "n_samples": n_samples,
    }


# --------------------------------------------------------------------------- #
# 3. Dual calibration routes: independent vs self-consistent (B9-T3 core)
# --------------------------------------------------------------------------- #

def _fdt_residual(cov_hat: np.ndarray, response_target: np.ndarray) -> float:
    """Relative FDT residual ``||Σ̂ − target|| / ||target||`` (B9 convention)."""
    return float(np.linalg.norm(cov_hat - response_target)
                 / np.linalg.norm(response_target))


def dual_route_recovery(structure: np.ndarray, scale_true: float,
                        hidden_bath: np.ndarray, n: int = 20000,
                        seed: int = 1) -> Dict:
    """Audit two calibration routes on data carrying a hidden bath.

    Fluctuation–dissipation ties the stationary covariance ``Σ = scale·K`` to a
    separately-declared response. Data is generated with a *hidden* extra bath
    (``Σ_true = scale·K + hidden_bath``) — an undeclared apparatus contribution.

    * **Spectroscopy route** (independent channel): the scale is supplied
      externally as ``scale_true``; the target response is ``scale_true·K``. The
      hidden bath shows up as a large FDT residual → the fraud is DETECTED.
    * **Gibbs route** (self-consistent): the response is *read off the same
      family of data* (a second self-consistent covariance estimate), so FDT is
      satisfied by construction up to sampling noise and the residual stays small
      → the route is FOOLED (the Voss–Webb ambiguity), and is non-promotable.

    Mirrors ``b9_circuit/estimate.py`` and its T3. Returns both residuals, both
    verdicts, and the promotable route.
    """
    K = np.asarray(structure, float)
    d = K.shape[0]
    bath = np.asarray(hidden_bath, float)
    sigma_true = scale_true * K + bath
    # PSD Cholesky factor to synthesize correlated samples.
    L = np.linalg.cholesky(sigma_true)
    rng = np.random.default_rng(seed)
    Z1 = (L @ rng.standard_normal((d, n))).T
    Z2 = (L @ rng.standard_normal((d, n))).T
    cov_hat = np.cov(Z1.T)

    # Independent (spectroscopy) route: response supplied by an external channel
    # as the *declared* equilibrium expectation scale_true·K. The hidden bath
    # makes the observed covariance disagree -> the fraud is detected.
    spec_target = scale_true * K
    spec_resid = _fdt_residual(cov_hat, spec_target)

    # Self-consistent (Gibbs) route: the response is a second covariance estimate
    # of the *same* data family, so it tracks cov_hat by construction (only
    # sampling noise remains) and FDT appears satisfied -> the route is fooled.
    gibbs_target = np.cov(Z2.T)
    scale_gibbs = float(np.trace(cov_hat @ np.linalg.inv(K)) / d)
    gibbs_resid = _fdt_residual(cov_hat, gibbs_target)

    def verdict(r: float) -> str:
        return ("EQUILIBRIUM_CONSISTENT" if r < FDT_THRESHOLD
                else "EQUILIBRIUM_REJECTED (hidden bath / undeclared apparatus)")

    return {
        "spectroscopy_route": "cal:spectroscopy_route",
        "gibbs_route": "cal:gibbs_route",
        "spec_residual": spec_resid,
        "gibbs_residual": gibbs_resid,
        "spec_verdict": verdict(spec_resid),
        "gibbs_verdict": verdict(gibbs_resid),
        "scale_true": float(scale_true),
        "scale_gibbs": scale_gibbs,
        "gibbs_more_forgiving": gibbs_resid < spec_resid,
        "promotable_route": "cal:spectroscopy_route",
        "gibbs_promotable": is_promotable_route("cal:gibbs_route"),
    }


# --------------------------------------------------------------------------- #
# 4. Effective rank at a resolution cutoff + the APPARATUS_LIMITED gate (S4)
# --------------------------------------------------------------------------- #

def effective_rank(matrix: np.ndarray, eps: float) -> int:
    """# eigenvalues above ``eps · λ_max`` — the ε-rank of ``s4_ds/kernel.py``.

    Mirrors S4 exactly: the positive spectrum is kept *before* thresholding, so a
    noisy/indefinite Gram cannot inflate the rank by counting large-magnitude
    negative eigenvalues as resolved modes (they are apparatus/indefiniteness
    artifacts, not S4 ε-rank modes).
    """
    lam = np.linalg.eigvalsh(np.asarray(matrix, float))
    lam = lam[lam > 0.0]
    if lam.size == 0:
        return 0
    return int(np.sum(lam > eps * lam.max()))


def apparatus_limited_gate(verdict_fn: Callable[[float], object],
                           resolutions: Sequence[float],
                           declared_window: Tuple[float, float]) -> Dict:
    """Sweep a structural verdict across the declared resolution/bandwidth.

    ``verdict_fn(r)`` returns the structural verdict at apparatus resolution
    ``r`` (e.g. an ε-cutoff or a bandwidth). ``declared_window`` is the
    ``(lo, hi)`` range the apparatus actually reaches. If the verdict is constant
    across every resolution in that window, the structural conclusion is
    permitted and passed through; if it changes *inside* the window, no
    structural conclusion is permitted and the gate returns SPEC §4
    ``APPARATUS_LIMITED`` naming the bounding threshold (the resolution at which
    the verdict first flips). Honest by construction: an ample apparatus (stable
    verdict) never returns APPARATUS_LIMITED.
    """
    lo, hi = declared_window
    grid = sorted(float(r) for r in resolutions)
    in_window = [r for r in grid if lo <= r <= hi]
    if not in_window:
        raise ValueError("declared_window contains no sampled resolution")
    vals = [(r, verdict_fn(r)) for r in in_window]
    distinct = {v for _, v in vals}
    if len(distinct) == 1:
        return {
            "verdict": "PASS",
            "apparatus_limited": False,
            "structural_verdict": next(iter(distinct)),
            "declared_window": [lo, hi],
            "resolutions_tested": len(vals),
            "cause": None,
        }
    # find the first resolution at which the verdict departs from the finest.
    boundary = None
    for (r_prev, v_prev), (r, v) in zip(vals, vals[1:]):
        if v != v_prev:
            boundary = r
            break
    return {
        "verdict": "APPARATUS_LIMITED",
        "apparatus_limited": True,
        "cause": "structural verdict changes within the declared "
                 "resolution/bandwidth",
        "bounding_threshold": boundary,
        "declared_window": [lo, hi],
        "resolutions_tested": len(vals),
        "verdicts_in_window": {f"{r:.3g}": (bool(v) if isinstance(v, (bool, np.bool_))
                                            else v) for r, v in vals},
    }


# --------------------------------------------------------------------------- #
# 5. Selection / censoring correction (declared truncation model)
# --------------------------------------------------------------------------- #

def _std_normal_pdf(x: float) -> float:
    return math.exp(-0.5 * x * x) / math.sqrt(2.0 * math.pi)


def _std_normal_cdf(x: float) -> float:
    return 0.5 * (1.0 + math.erf(x / math.sqrt(2.0)))


def _std_normal_sf(x: float) -> float:
    """Upper-tail survival ``1 − Φ(x)`` via ``erfc`` — stable deep in the tail
    (``erfc`` stays accurate to ~1e-300, where ``1 − Φ`` from ``erf`` underflows
    to exactly 0)."""
    return 0.5 * math.erfc(x / math.sqrt(2.0))


def _inverse_mills_upper(alpha: float) -> float:
    """Upper-tail inverse Mills ratio ``φ(α)/(1−Φ(α))``.

    For large ``α`` both ``φ`` and the survival function underflow, so the ratio
    is taken from its asymptotic series ``α + 1/α − 2/α³`` (→ ``α`` as ``α→∞``),
    which keeps :func:`truncated_normal_mean_mle` monotone and identifiable far
    into the heavily-censored regime instead of flattening to the truncation
    floor.
    """
    if alpha > 25.0:
        return alpha + 1.0 / alpha - 2.0 / (alpha ** 3)
    q = _std_normal_sf(alpha)
    if q < 1e-300:
        return alpha + 1.0 / alpha - 2.0 / (alpha ** 3)
    return _std_normal_pdf(alpha) / q


def truncated_normal_mean_mle(observed_mean: float, threshold: float,
                              sigma: float) -> float:
    """Recover μ of a normal from a *left-truncated* sample mean (x > threshold).

    For ``x ~ N(μ, σ²)`` observed only when ``x > t``,
    ``E[x | x > t] = μ + σ·φ(α)/(1−Φ(α))`` with ``α = (t−μ)/σ`` (inverse Mills
    ratio). This is monotone increasing in μ, so a bisection on μ inverts the
    declared censoring model. The inverse Mills ratio uses a tail-stable survival
    function (:func:`_inverse_mills_upper`) so heavily-censored samples stay
    identifiable instead of collapsing to the truncation floor. The naive sample
    mean (ignoring the interface) is biased high; this removes the bias.
    """
    def trunc_mean(mu: float) -> float:
        alpha = (threshold - mu) / sigma
        return mu + sigma * _inverse_mills_upper(alpha)

    # Wide bracket: trunc_mean → threshold as μ → −∞ and exceeds any observed
    # mean (> threshold) at μ = observed_mean, so the root is bracketed.
    lo, hi = threshold - 60.0 * sigma, observed_mean
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if trunc_mean(mid) < observed_mean:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def censoring_audit(mu_true: float, sigma: float, threshold: float,
                    n: int = 40000, seed: int = 2) -> Dict:
    """Naive vs interface-aware estimate of μ under threshold censoring.

    Draw ``x ~ N(μ, σ²)`` but observe only ``x > threshold`` (a detector
    selection). The naive mean of the survivors is biased; the declared
    truncation model recovers μ.
    """
    rng = np.random.default_rng(seed)
    x = rng.normal(mu_true, sigma, n)
    obs = x[x > threshold]
    naive = float(obs.mean())
    corrected = float(truncated_normal_mean_mle(naive, threshold, sigma))
    return {
        "mu_true": float(mu_true),
        "threshold": float(threshold),
        "n_observed": int(obs.size),
        "selection_fraction": float(obs.size) / n,
        "naive_estimate": naive,
        "naive_error": abs(naive - mu_true),
        "corrected_estimate": corrected,
        "corrected_error": abs(corrected - mu_true),
    }
