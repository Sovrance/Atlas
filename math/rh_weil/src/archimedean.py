"""Archimedean Mellin multiplier h_+(τ) via Arb digamma (E1-capable)."""
from __future__ import annotations

from typing import Any, Optional, Tuple

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


def h_plus_derivatives(tau: Any) -> Tuple[Any, Any, Any]:
    """Return (h_+, h_+', h_+'') with respect to τ.

    z = 1/4 + i τ/2,  dz/dτ = i/2.
    h'  = Re(ψ'(z)·i/2) = −Im(ψ'(z))/2
    h'' = Re(ψ''(z)·(i/2)²) = −Re(ψ''(z))/4
    """
    _, arb, acb, _ = require_flint()
    tau_a = arb(tau) if not hasattr(tau, "mid") else tau
    z = acb(arb("0.25"), tau_a / 2)
    h = z.digamma().real - arb.pi().log()
    p1 = z.polygamma(1)
    p2 = z.polygamma(2)
    hp = -p1.imag / 2
    hpp = -p2.real / 4
    return h, hp, hpp


def h_plus_log_majorant(T: float) -> float:
    """Pointwise majorant: h_+(t) ≤ log(t) + 1 for t ≥ T ≥ 2.

    Justification: Re ψ(z) = log|z| + O(1/|z|) with z=1/4+it/2, so
    h_+(t) = log(|t|/2) − log π + ε(t) with |ε|→0; the constant 1 covers
    the transient on [2, ∞) (validated by direct Arb sampling at certification).
    """
    import math

    if T < 2:
        raise ValueError("majorant stated for T ≥ 2")
    return math.log(T) + 1
