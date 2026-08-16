"""Cutoff-free degree-1 odd pivot assembly (WO-RH-10)."""
from __future__ import annotations

from typing import Any, Dict

import core
from archimedean import h_plus
from finite_weil import g0_odd_pivot, gp_odd_pivot
from interval_backend import require_flint, set_precision_bits, is_definitely_positive, backend_info


def ginf_odd_cutoff_free(L, arb, T_cap: int = 4000, n: int = 2000):
    """Approximate Ginf[q1,q1] by truncated Fourier integral (pending tail bound)."""
    from finite_weil import stable_products_even

    T_a = arb(T_cap)
    L_a = arb(L)
    acc = arb(0)
    for i in range(n + 1):
        t = T_a * i / n
        w = arb("0.5") if i in (0, n) else arb(1)
        hp = h_plus(t)
        if t.contains(0) or abs(t.mid()) <= max(t.rad(), 1e-30):
            hq1sq = arb(0)
        else:
            _, _, hbsq = stable_products_even(t, L_a, arb)
            hq1sq = hbsq * (t**2) / 4
        acc += w * hp * hq1sq
    return acc * (T_a / n) / arb.pi()


def assemble_O1(L, precision_bits: int = 256) -> Dict[str, Any]:
    _, arb, _, _ = require_flint()
    set_precision_bits(precision_bits)
    L_a = arb(L)
    g0 = g0_odd_pivot(L_a, arb)
    gp = gp_odd_pivot(L_a, arb)
    gi = ginf_odd_cutoff_free(L_a, arb)
    O1 = g0 - gp + gi
    return {
        "O1": O1,
        "G0": g0,
        "Gp": gp,
        "Ginf_trunc": gi,
        "normalization": core.NORMALIZATION,
        "evidence_class": "E3_PENDING_TAIL_BOUND",
        "rh_proof_claim": False,
        "backend": backend_info(precision_bits).to_dict(),
        "notebook_regression_target": "0.0142397928900392393805",
        "note": "Regenerated assembly; tail bound required before E1 promotion.",
    }


def certify_degree1_e1(precision_bits: int = 256) -> Dict[str, Any]:
    """Attempt cell-wide positivity; returns E1 only if definitely positive with bound."""
    import math

    # Sample endpoints + interior; full interval cover deferred.
    samples = [math.log(3), 1.1059498113, 1.20, math.log(4)]
    rows = []
    all_pos = True
    min_lower = None
    for L in samples:
        a = assemble_O1(L, precision_bits=precision_bits)
        o1 = a["O1"]
        low = float(o1.lower())
        rows.append({"L": L, "O1_lower": low, "O1_upper": float(o1.upper())})
        all_pos = all_pos and is_definitely_positive(o1)
        min_lower = low if min_lower is None else min(min_lower, low)
    # Do not claim E1 until tail bound + interval cover exist.
    return {
        "certificate_version": "0.2",
        "program": "RH/Weil degree-1 odd pivot",
        "work_order": "WO-RH-10",
        "evidence_class": "E3",
        "status": "ASSEMBLY_REGENERATED_PENDING_TAIL_AND_COVER",
        "hard_constraints_certified": False,
        "rh_proof_claim": False,
        "normalization": core.NORMALIZATION,
        "domain": {"L_interval": ["log(3)", "log(4)"]},
        "samples": rows,
        "min_sample_lower": min_lower,
        "all_samples_definitely_positive": all_pos,
        "backend": backend_info(precision_bits).to_dict(),
        "note": (
            "Pole + prime + truncated archimedean assembly implemented. "
            "E1 blocked on archimedean tail bound and full interval cover."
        ),
    }
