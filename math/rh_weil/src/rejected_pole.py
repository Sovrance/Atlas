"""ARCHIVAL — the REJECTED Candidate-B pole block. Never import from production.

WO-RH-17 rejected the block

    B_ij = (sqrt(3)/2) (E_i^+ E_j^+ + E_i^- E_j^-)

as the pole contribution to the finite Weil Gram. On the even sector it equals
the adopted Candidate A (``pole.pole_gram_entry``) multiplied by

    B/A = (sqrt(3)/2) cosh(L/2),

a factor that equals 1 at exactly one point, ``L = log 3``, and differs
elsewhere (+0.18% at L=1.1059498113, +2.66% at L=1.20, +8.25% at L=log 4).
A change of basis is a *constant* congruence, so an ``L``-dependent factor
cannot come from one: B is not a renormalised A, it is a multiplicative
calibration fitted at a single test point.

This module exists so the adjudication stays reproducible and so the rejected
expression can be exhibited in audit artifacts. It is **not** production code:

* ``tests/test_production_imports.py`` fails CI if any module under ``src/``
  other than this one imports it;
* nothing here may be used to build a certificate that is promoted.

Legitimate importers are the adjudication scripts
(``scripts/derive_normalization.py``, ``scripts/run_normalization_crosscheck.py``)
and the adjudication tests.

No RH proof claim is made by this module.
"""
from __future__ import annotations

from math import cosh, sqrt
from typing import Any, Tuple

import pole

LEGACY_POLE_SCALE = "sqrt(3)/2"
LEGACY_STATUS = "REJECTED_FITTED_CALIBRATION"
LEGACY_WORK_ORDER = "WO-RH-17"
CALIBRATION_FIXED_POINT = "L = log 3"
SUPERSEDED_BY = "pole.pole_gram_entry (Candidate A)"


def _sqrt3_over_2(carrier: Any) -> Any:
    three = 0 * carrier + 3
    return (three.sqrt() / 2) if hasattr(three, "sqrt") else (sqrt(3) / 2)


def legacy_pole_entry(i: str, j: str, L: Any) -> Any:
    """``(sqrt(3)/2)(E_i^+E_j^+ + E_i^-E_j^-)`` — REJECTED; audit use only."""
    s = sqrt(3) / 2 if isinstance(L, float) else _sqrt3_over_2(L)
    eip, eim = pole.laplace_plus(i, L), pole.laplace_minus(i, L)
    ejp, ejm = pole.laplace_plus(j, L), pole.laplace_minus(j, L)
    return s * (eip * ejp + eim * ejm)


def legacy_even_block(L: Any) -> Tuple[Any, Any, Any]:
    """The rejected even block ``(g00, g0b, gbb)`` — audit use only."""
    return (
        legacy_pole_entry("one", "one", L),
        legacy_pole_entry("one", "b", L),
        legacy_pole_entry("b", "b", L),
    )


def legacy_over_adopted_ratio(L: Any) -> Any:
    """Even-sector ratio ``B/A = (sqrt(3)/2) cosh(L/2)``.

    On the even sector ``h(L-x) = h(x)`` gives ``E^- = e^{-L/2}E^+``, hence
    ``A = 2 e^{-L/2} v+ v+^T`` and ``B = (sqrt(3)/2)(1 + e^{-L}) v+ v+^T``.
    The quotient is ``(sqrt(3)/4)(e^{L/2} + e^{-L/2}) = (sqrt(3)/2)cosh(L/2)``,
    which is 1 exactly at ``L = log 3`` and nowhere else.
    """
    if isinstance(L, float):
        return (sqrt(3) / 2) * cosh(L / 2)
    half = L / 2
    ch = half.cosh() if hasattr(half, "cosh") else cosh(float(half))
    return _sqrt3_over_2(L) * ch
