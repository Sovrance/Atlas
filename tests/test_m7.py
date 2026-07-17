"""M7 measurement-interface tests -- 'the declared 𝖬' as executable code.

T1  INJECTION-RECOVERY (B9-T3 generalized, single route): a known latent pushed
    through a declared affine transfer + Gaussian noise is recovered without
    bias, and the certified z-interval attains its nominal coverage. A declared
    interface + stated noise family recovers ground truth.
T2  INDEPENDENT vs SELF-CONSISTENT CALIBRATION (B9-T3 core): on data carrying a
    hidden bath, the independent spectroscopy route REJECTS via the FDT audit
    (residual > 0.15) while the self-consistent Gibbs route is FOOLED
    (residual < 0.15) and is flagged non-promotable -- the Voss-Webb ambiguity.
    'Self-consistent' is never promotable evidence.
T3  APPARATUS_LIMITED TRIGGER (S4 generalized): a 'rank >= 3' structural verdict
    on a Gram with a near-floor eigenvalue flips inside the declared resolution
    window; the gate returns SPEC-§4 APPARATUS_LIMITED naming the bounding
    threshold -- no structural conclusion permitted at the detector cutoff.
T4  SELECTION / CENSORING: a threshold-censored sample mean is biased; the
    declared truncation model removes the bias (interface-aware error << naive
    error). Ignoring the interface is the error, not the noise.
T5  NO-COLLAPSE (honesty): with an ample apparatus (verdict stable across the
    whole declared window) the gate does NOT cry APPARATUS_LIMITED -- it passes
    the structural verdict (rank >= 3 = True) through. APPARATUS_LIMITED is not a
    synonym for 'noisy'.
"""

import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np

from measurement_interface.engine import (
    FDT_THRESHOLD,
    apparatus_limited_gate,
    censoring_audit,
    dual_route_recovery,
    effective_rank,
    injection_recovery,
    is_promotable_route,
)
from b1_moment_solver.certificate import save_certificate


def _gram(eigs, seed):
    """A symmetric PSD matrix with the given spectrum, in a random basis."""
    Q, _ = np.linalg.qr(np.random.default_rng(seed).standard_normal((len(eigs), len(eigs))))
    return Q @ np.diag(np.asarray(eigs, float)) @ Q.T


def t1_injection_recovery():
    r = injection_recovery(3.5, {"gain": 2.0, "offset": -1.0}, sigma=0.8,
                           n_samples=400, n_trials=400, seed=0)
    assert abs(r["bias"]) < 0.02, ("recovery biased", r["bias"])
    # empirical coverage within Monte-Carlo tolerance of nominal 95%
    assert abs(r["coverage"] - r["nominal_coverage"]) < 0.03, ("coverage off", r)
    print(f"T1 injection-recovery: PASS  bias={r['bias']:.4f} "
          f"coverage={r['coverage']:.3f} (nominal {r['nominal_coverage']:.3f})")
    return r


def t2_dual_calibration_route():
    K = np.array([[1.0, 0.3], [0.3, 1.0]])
    bath = np.diag([0.9, 0.0])  # undeclared apparatus contribution on channel 0
    r = dual_route_recovery(K, scale_true=1.0, hidden_bath=bath, n=20000, seed=1)
    # independent route detects the hidden bath; self-consistent route is fooled
    assert r["spec_residual"] > FDT_THRESHOLD, ("spectroscopy failed to detect", r)
    assert r["gibbs_residual"] < FDT_THRESHOLD, ("gibbs should be fooled", r)
    # the exact B9-T3 shape: the Gibbs route is strictly more forgiving
    assert r["gibbs_residual"] < 0.6 * r["spec_residual"], ("not more forgiving", r)
    assert not is_promotable_route("cal:gibbs_route")
    assert is_promotable_route("cal:spectroscopy_route")
    print(f"T2 dual calibration route: PASS  spec={r['spec_residual']:.3f} REJECTED, "
          f"gibbs={r['gibbs_residual']:.3f} fooled (non-promotable)")
    return r


def t3_apparatus_limited_trigger():
    # a Gram whose 3rd eigenvalue sits inside the declared resolution window
    G = _gram([1.0, 0.5, 2e-3, 1e-4], seed=5)
    resolutions = [1e-2, 5e-3, 1e-3, 5e-4, 1e-4]
    verdict_fn = lambda eps: bool(effective_rank(G, eps) >= 3)
    gate = apparatus_limited_gate(verdict_fn, resolutions, declared_window=(5e-4, 1e-2))
    assert gate["apparatus_limited"] is True, ("should be apparatus-limited", gate)
    assert gate["verdict"] == "APPARATUS_LIMITED"
    assert gate["bounding_threshold"] is not None, ("no threshold named", gate)
    print(f"T3 APPARATUS_LIMITED trigger: PASS  verdict flips at eps="
          f"{gate['bounding_threshold']:.1e} within the declared window")
    return gate


def t4_censoring_correction():
    r = censoring_audit(mu_true=0.0, sigma=1.0, threshold=0.5, n=40000, seed=2)
    assert r["naive_error"] > 0.5, ("naive should be badly biased", r)
    assert r["corrected_error"] < 0.05, ("declared model should correct it", r)
    assert r["corrected_error"] < 0.1 * r["naive_error"], ("insufficient correction", r)
    print(f"T4 censoring correction: PASS  naive_err={r['naive_error']:.3f} -> "
          f"corrected_err={r['corrected_error']:.4f} (sel {r['selection_fraction']:.2f})")
    return r


def t5_no_collapse():
    # same Gram, but an AMPLE apparatus whose window sits entirely in the
    # fine-resolution regime where rank>=3 holds stably.
    G = _gram([1.0, 0.5, 2e-3, 1e-4], seed=5)
    resolutions = [1e-2, 5e-3, 1e-3, 5e-4, 1e-4]
    verdict_fn = lambda eps: bool(effective_rank(G, eps) >= 3)
    gate = apparatus_limited_gate(verdict_fn, resolutions, declared_window=(1e-4, 1e-3))
    assert gate["apparatus_limited"] is False, ("must not be apparatus-limited", gate)
    assert gate["verdict"] == "PASS"
    assert gate["structural_verdict"] is True, ("ample apparatus licenses rank>=3", gate)
    print("T5 no-collapse (honesty): PASS  stable verdict passes through "
          "(rank>=3 = True), no false APPARATUS_LIMITED")
    return gate


if __name__ == "__main__":
    r1 = t1_injection_recovery()
    r2 = t2_dual_calibration_route()
    r3 = t3_apparatus_limited_trigger()
    r4 = t4_censoring_correction()
    r5 = t5_no_collapse()
    cert = {
        "certificate_version": "0.3",
        "certificate_class": "NUMERICAL-DISCOVERY (measurement-interface layer; "
                             "apparatus/calibration/censoring audits, E3)",
        "problem": "M7: measurement-interface layer (SPEC §2 layer 3 as "
                   "executable code)",
        "timestamp_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "calibration_route": "cal:spectroscopy_route (independent channel; the "
                             "self-consistent cal:gibbs_route is audited as "
                             "non-promotable in T2)",
        "m_layer_stipulations": [
            "transfer function is the declared affine apparatus map y = gain*x + "
            "offset; a degenerate (gain 0) apparatus is rejected, not inverted",
            "noise family is additive i.i.d. Gaussian with a stated sigma; the "
            "certified interval is a normal z-interval (nominal 95%)",
            "FDT rejection threshold 0.15 carried over verbatim from B9; a "
            "self-consistent (data-inferred) calibration route is non-promotable "
            "by construction (Voss-Webb ambiguity), aligned with "
            "pir/domains/circuit_semantics.CALIBRATION_ROUTES",
            "APPARATUS_LIMITED is defined relative to the DECLARED "
            "resolution/bandwidth window; eps-rank uses eps * lambda_max exactly "
            "as s4_ds/kernel.py",
            "censoring is a hard left-truncation x > threshold with a KNOWN sigma; "
            "the correction is the truncated-normal (inverse Mills ratio) MLE. "
            "continuous/functional transfers, calibration priors, and recovery "
            "over real B4/B10 posteriors are deferred to M7-b",
        ],
        "headline": "the measurement interface as code: a declared apparatus "
                    "recovers a latent without bias (T1), an independent "
                    "calibration route detects a hidden bath that the "
                    "self-consistent route is fooled by and cannot promote (T2), "
                    "a rank verdict that flips inside the detector's resolution "
                    "window is refused as APPARATUS_LIMITED (T3), declared "
                    "censoring is corrected (T4), and an ample apparatus is NOT "
                    "cried apparatus-limited (T5).",
        "results": {"T1": r1, "T2": r2, "T3": r3, "T4": r4, "T5": r5},
    }
    save_certificate(cert, os.path.join(os.path.dirname(__file__), "..",
                                        "certificates", "m7_certificate.json"))
    print("\nCertificate written: certificates/m7_certificate.json")
    print("ALL M7 TESTS PASS")
