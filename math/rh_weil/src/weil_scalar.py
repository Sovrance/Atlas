"""Scalar cutoff-free verifier path (WO-RH-09).

Absolute G00 uses the even (1,1) entry of the true Weil assembly
G = G0 − Gp + Ginf_∞ with even pole outer-product and an archimedean tail
majorant. Notebook 0.075… is regression-only (regenerated scale differs).
"""
from __future__ import annotations

from typing import Any, Dict

import core
import scalar
from finite_weil import (
    g0_even_block,
    gp_even_block,
    ginf_even_block_quad,
    archimedean_tail_even,
)
from interval_backend import (
    backend_info,
    is_definitely_positive,
    require_flint,
    set_precision_bits,
)


def w00_second_arb(L, arb):
    r = L.exp()
    return 2 * (r**3 - r - 1) / (r.sqrt() * (r**2 - 1))


def assemble_G00_cutoff_free(
    L,
    precision_bits: int = 256,
    T_cap: int = 2000,
    n: int = 8000,
) -> Dict[str, Any]:
    """Cutoff-free scalar entry = even G00 with tail-enclosed Ginf."""
    _, arb, _, _ = require_flint()
    set_precision_bits(precision_bits)
    L_a = arb(L)
    g000, _, _ = g0_even_block(L_a, arb)
    gp00, _, _ = gp_even_block(L_a, arb)
    gi00, _, _ = ginf_even_block_quad(L_a, T_cap, arb, n=n)
    t00, _, _ = archimedean_tail_even(L_a, float(T_cap), arb)
    gi = gi00 + arb(0, float(t00.upper()))
    G00 = g000 - gp00 + gi
    return {
        "G00": G00,
        "G0": g000,
        "Gp": gp00,
        "Ginf_trunc": gi00,
        "tail_majorant": t00,
        "T_cap": T_cap,
        "normalization": core.NORMALIZATION,
        "rh_proof_claim": False,
    }


def certify_scalar_e1(precision_bits: int = 256) -> Dict[str, Any]:
    """Structural W00'' E0 + regenerated absolute G00 samples with tail bound."""
    _, arb, _, _ = require_flint()
    set_precision_bits(precision_bits)
    report = scalar.verify_scalar_cell().to_dict()
    import math

    samples = [math.log(3) + 1e-9, 1.1059498113, 1.20, math.log(4) - 1e-9]
    rows = []
    all_curv_pos = True
    all_g00_pos = True
    min_g00 = None
    min_curv = None
    for L in samples:
        w = w00_second_arb(arb(L), arb)
        a = assemble_G00_cutoff_free(L, precision_bits=precision_bits)
        g00 = a["G00"]
        rows.append(
            {
                "L": L,
                "W00_second_lower": float(w.lower()),
                "W00_second_upper": float(w.upper()),
                "G00_lower": float(g00.lower()),
                "G00_upper": float(g00.upper()),
                "G00_definitely_positive": is_definitely_positive(g00),
            }
        )
        all_curv_pos = all_curv_pos and is_definitely_positive(w)
        all_g00_pos = all_g00_pos and is_definitely_positive(g00)
        min_curv = (
            float(w.lower())
            if min_curv is None
            else min(min_curv, float(w.lower()))
        )
        min_g00 = (
            float(g00.lower())
            if min_g00 is None
            else min(min_g00, float(g00.lower()))
        )

    # Full cell E1 needs interval cover of the minimizer; samples + convexity
    # give a provisional lower bound once G00>0 on a fine net and W00''>0.
    notebook_tgt = 0.0753795566117244
    notebook_in_min = min_g00 is not None and min_g00 <= notebook_tgt <= (
        rows[1]["G00_upper"] if rows else notebook_tgt
    )

    return {
        "certificate_version": "0.3",
        "program": "RH/Weil scalar cell numerical path",
        "work_order": "WO-RH-09",
        "evidence_class": "E1_SAMPLES_PLUS_E0_CURVATURE" if all_g00_pos else "E3",
        "status": (
            "ABSOLUTE_G00_REGENERATED_PENDING_FULL_CELL_COVER"
            if all_g00_pos
            else "G00_ASSEMBLY_FAILED"
        ),
        "hard_constraints_certified": bool(report.get("w00_second_positive"))
        and all_curv_pos
        and all_g00_pos,
        "rh_proof_claim": False,
        "normalization": core.NORMALIZATION,
        "domain": {"L_interval": ["log(3)", "log(4)"]},
        "structural_report": report,
        "samples": rows,
        "min_sample_W00_second_lower": min_curv,
        "min_sample_G00_lower": min_g00,
        "all_sample_G00_definitely_positive": all_g00_pos,
        "backend": backend_info(precision_bits).to_dict(),
        "notebook_regression_target_G00": str(notebook_tgt),
        "notebook_target_in_regenerated_range": notebook_in_min,
        "note": (
            "W00''>0 remains E0 algebraic. Absolute G00 uses even pole "
            "outer-product + truncated Fourier + archimedean tail majorant. "
            "Regenerated scale (~0.12–0.13 near the cell minimizer) differs from "
            "the notebook 0.075… figure; notebook value is not a warrant. "
            "Promote to full-cell E1 after interval cover of the unique minimizer."
        ),
    }
