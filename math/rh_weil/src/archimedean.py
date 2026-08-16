"""Archimedean Mellin multiplier h_+(τ) via Arb digamma (E1-capable)."""
from __future__ import annotations

from typing import Any, Optional

from interval_backend import require_flint, set_precision_bits


def h_plus(tau: Any, precision_bits: Optional[int] = None):
    """h_+(τ) = Re ψ(1/4 + i τ/2) − log(π), Arb ball.

    Uses the current FLINT context precision unless ``precision_bits`` is
    explicitly supplied. Callers that set a working precision via
    ``set_precision_bits`` must not have it silently overwritten.
    """
    _, arb, acb, _ = require_flint()
    if precision_bits is not None:
        set_precision_bits(precision_bits)
    tau_a = arb(tau) if not hasattr(tau, "mid") else tau
    z = acb(arb("0.25"), tau_a / 2)
    return z.digamma().real - arb.pi().log()
