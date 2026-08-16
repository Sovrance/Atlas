"""Archimedean Mellin multiplier h_+(τ) via Arb digamma (E1-capable)."""
from __future__ import annotations

from typing import Any

from interval_backend import require_flint, set_precision_bits


def h_plus(tau: Any, precision_bits: int = 256):
    """h_+(τ) = Re ψ(1/4 + i τ/2) − log(π), Arb ball."""
    _, arb, acb, _ = require_flint()
    set_precision_bits(precision_bits)
    tau_a = arb(tau) if not hasattr(tau, "mid") else tau
    z = acb(arb("0.25"), tau_a / 2)
    return z.digamma().real - arb.pi().log()
