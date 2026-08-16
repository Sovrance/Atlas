"""Rigorous interval backend for RH/Weil E1 certificates (WO-RH-09+).

Preferred backend: ``python-flint`` Arb/Acb.
E1 emission MUST call :func:`require_flint` — mpmath fallback is forbidden for E1.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any


class FlintUnavailable(RuntimeError):
    """Raised when an E1 path requires python-flint but it is not installed."""


def flint_available() -> bool:
    try:
        import flint  # noqa: F401

        return True
    except ImportError:
        return False


def require_flint():
    """Import flint or fail loudly (do not emit E1 without this)."""
    try:
        import flint
        from flint import arb, acb, ctx
    except ImportError as exc:  # pragma: no cover
        raise FlintUnavailable(
            "E1 certificates require python-flint (Arb/Acb). "
            "Install python-flint; do not fall back to mpmath for E1."
        ) from exc
    return flint, arb, acb, ctx


def set_precision_bits(bits: int) -> int:
    _, _, _, ctx = require_flint()
    ctx.prec = int(bits)
    return int(ctx.prec)


def arb_from(x: Any):
    _, arb, _, _ = require_flint()
    if isinstance(x, str):
        return arb(x)
    return arb(x)


def arb_ball(mid: Any, rad: Any):
    """Create an Arb ball with explicit midpoint and radius."""
    _, arb, _, _ = require_flint()
    return arb(mid, rad)


@dataclass(frozen=True)
class BackendInfo:
    name: str
    version: str | None
    precision_bits: int

    def to_dict(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "version": self.version,
            "precision_bits": self.precision_bits,
        }


def backend_info(precision_bits: int | None = None) -> BackendInfo:
    import importlib.metadata as md

    flint, _, _, ctx = require_flint()
    bits = int(precision_bits if precision_bits is not None else ctx.prec)
    try:
        ver = md.version("python-flint")
    except md.PackageNotFoundError:
        ver = None
    return BackendInfo(name="python-flint/Arb", version=ver, precision_bits=bits)


def lower_bound(x) -> float:
    """Outward-rounded float lower endpoint for reporting (not for arithmetic)."""
    return float(x.lower())


def upper_bound(x) -> float:
    return float(x.upper())


def is_definitely_positive(x) -> bool:
    return x.lower() > 0


def is_definitely_negative(x) -> bool:
    return x.upper() < 0
