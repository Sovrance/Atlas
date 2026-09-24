"""Exact Fraction helpers for B15 (copies the b7_onsager/exact.py pattern)."""

from fractions import Fraction
from typing import Any

Q = Fraction


def q(x: Any) -> Fraction:
    """Coerce to Fraction (strings ``"p/q"`` accepted)."""
    return x if isinstance(x, Fraction) else Fraction(x)


def fmt(x: Fraction) -> str:
    """Lowest-terms ``"p/q"`` string with explicit denominator (pir/canonical
    convention) so certificate values never pass through floats."""
    x = q(x)
    return f"{x.numerator}/{x.denominator}"
