"""Scalar cutoff-free verifier path (WO-RH-09).

Uses algebraic W00'' (E0) plus Arb sampling. Absolute E1 lower bound requires
a calibrated continuous antiderivative / full explicit-formula assembly —
status remains pending until that calibration is first-principles closed.
"""
from __future__ import annotations

from typing import Any, Dict

import core
import scalar
from interval_backend import (
    backend_info,
    is_definitely_positive,
    require_flint,
    set_precision_bits,
)


def w00_second_arb(L, arb):
    r = L.exp()
    return 2 * (r**3 - r - 1) / (r.sqrt() * (r**2 - 1))


def certify_scalar_e1(precision_bits: int = 256) -> Dict[str, Any]:
    """Structural + sample curvature positivity; not a full G00 E1 yet."""
    _, arb, _, _ = require_flint()
    set_precision_bits(precision_bits)
    report = scalar.verify_scalar_cell().to_dict()
    import math

    samples = [math.log(3) + 1e-9, 1.1059498113, 1.20, math.log(4) - 1e-9]
    rows = []
    all_pos = True
    min_lower = None
    for L in samples:
        w = w00_second_arb(arb(L), arb)
        low = float(w.lower())
        rows.append({"L": L, "W00_second_lower": low, "W00_second_upper": float(w.upper())})
        all_pos = all_pos and is_definitely_positive(w)
        min_lower = low if min_lower is None else min(min_lower, low)
    return {
        "certificate_version": "0.2",
        "program": "RH/Weil scalar cell numerical path",
        "work_order": "WO-RH-09",
        "evidence_class": "E0_PLUS_SAMPLES",
        "status": "STRUCTURAL_E0_RETAINED_G00_E1_PENDING_CALIBRATION",
        "hard_constraints_certified": bool(report.get("w00_second_positive")) and all_pos,
        "rh_proof_claim": False,
        "normalization": core.NORMALIZATION,
        "domain": {"L_interval": ["log(3)", "log(4)"]},
        "structural_report": report,
        "curvature_samples": rows,
        "min_sample_W00_second_lower": min_lower,
        "backend": backend_info(precision_bits).to_dict(),
        "notebook_regression_target_G00": "0.0753795566117244",
        "note": (
            "W00'' positivity remains E0 algebraic. Full G00(L) E1 lower bound "
            "awaits first-principles absolute calibration of the continuous "
            "antiderivative (notebook value is regression-only)."
        ),
    }
