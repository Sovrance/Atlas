"""What the low moments do and do not pin down about a spectrum (§6, WO-RH-31).

The eigenvalues of a finite Hermitian ``G`` are the support of a measure::

    mu = sum_i delta_{lambda_i},    m_0 = n,    m_k = tr(G^k)

so asking "what spectra are consistent with these moments" *is* a truncated
moment problem, and Atlas already has an exact engine for that in
``b1_moment_solver``. §6 says reuse it, so the exact-rational path here calls
B1's Hankel/PSD/rank routines directly rather than growing a second solver.

Where the two engines divide
----------------------------
B1 is exact-rational by construction. The moments of a Weil block are interval
enclosures, not rationals, so a B1 verdict computed at one rational point inside
an enclosure would say nothing about the enclosure as a whole. The interval path
therefore asks the ENG-006 inertia engine for the definiteness of the same
Hankel and localizing matrices -- also an existing component, and one that fails
closed. Exact inputs go to B1; interval inputs go to the inertia engine; nothing
new is invented for either.

Which conclusions are actually available
----------------------------------------
The honest asymmetry, and the reason "insufficient information" is a first-class
result here:

* A localizing matrix that is **not** PSD is conclusive. If ``[m_{i+j+1}]`` fails
  to be PSD then no representing measure lives on ``[0, inf)``, so ``G`` has a
  negative eigenvalue. That is a proof.
* A localizing matrix that **is** PSD proves nothing on its own. With finitely
  many moments, PSD-ness of the truncated localizing matrix is necessary but not
  sufficient for support in ``[0, inf)`` -- sufficiency needs a flat extension
  (Curto-Fialkow), which four moments generally do not provide. So "the moments
  force PSD" is usually *not* a conclusion the data support, and this module
  says so rather than implying it.

Two-sided eigenvalue bounds come from ``m1`` and ``m2`` alone, via the
Wolkowicz-Styan inequalities: with ``mu = m1/n`` and ``s^2 = m2/n - mu^2``,

    mu - s*sqrt(n-1)  <=  lambda_min  <=  mu - s/sqrt(n-1)
    mu + s/sqrt(n-1)  <=  lambda_max  <=  mu + s*sqrt(n-1)

At ``n = 2`` these collapse to equalities, so for a 2x2 block the first two
moments determine the spectrum outright -- which is exactly the situation the
degree-3 odd block is in.

No RH proof claim is made by this module.
"""
from __future__ import annotations

import sys
from fractions import Fraction
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

_REPO = Path(__file__).resolve().parents[3]
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))

_PROGRAM = Path(__file__).resolve().parents[1]
if str(_PROGRAM) not in sys.path:
    sys.path.insert(0, str(_PROGRAM))

#: Result statuses. ``INSUFFICIENT`` is a certified answer, not a failure (§6).
CONCLUSIVE = "CONCLUSIVE"
INSUFFICIENT = "INSUFFICIENT_INFORMATION"


def b1_available() -> bool:
    try:
        import b1_moment_solver.exact  # noqa: F401

        return True
    except ImportError:
        return False


def hankel_matrix(moments: Sequence[Any], t: int) -> List[List[Any]]:
    """``M_t`` with ``(M_t)_{ij} = m_{i+j}``; ``moments`` is indexed from ``m_0``."""
    if len(moments) < 2 * t + 1:
        raise ValueError(f"need m_0..m_{2*t} for M_{t}, got {len(moments)}")
    return [[moments[i + j] for j in range(t + 1)] for i in range(t + 1)]


def localizing_matrix(moments: Sequence[Any], t: int, *, shift: Any = None,
                      flip: bool = False) -> List[List[Any]]:
    """Localizing matrix for ``x - shift`` (or ``shift - x`` when ``flip``).

    PSD-ness of this matrix is necessary for the representing measure to live on
    ``[shift, inf)`` (or ``(-inf, shift]``). Its *failure* is the conclusive
    direction.
    """
    need = 2 * t + 2
    if len(moments) < need:
        raise ValueError(f"need m_0..m_{need-1} for the localizing matrix, "
                         f"got {len(moments)}")
    out = []
    for i in range(t + 1):
        row = []
        for j in range(t + 1):
            v = moments[i + j + 1]
            if shift is not None:
                v = v - shift * moments[i + j]
            row.append(-v if flip else v)
        out.append(row)
    return out


def _is_exact(x) -> bool:
    return isinstance(x, (int, Fraction))


def psd_status(M: Sequence[Sequence[Any]]) -> Dict[str, Any]:
    """Definiteness of ``M``: B1 when the entries are exact, inertia otherwise."""
    exact = all(_is_exact(x) for row in M for x in row)
    if exact and b1_available():
        from b1_moment_solver.exact import psd_certificate

        status, pivots, rank = psd_certificate([[Fraction(x) for x in row] for row in M])
        return {
            "engine": "b1_moment_solver.exact.psd_certificate",
            "status": status,
            "is_psd": status in ("PSD_CERTIFIED", "PD_CERTIFIED"),
            "definitely_not_psd": status == "NOT_PSD_CERTIFIED",
            "rank": rank,
            "pivots": [str(p) for p in pivots],
        }
    from inertia.ldl import exact_inertia, interval_inertia

    res = exact_inertia(M) if exact else interval_inertia(M)
    is_psd = res.status == "PASS" and res.n_negative == 0
    not_psd = res.status == "PASS" and res.n_negative > 0
    return {
        "engine": "inertia.ldl." + ("exact_inertia" if exact else "interval_inertia"),
        "status": res.status,
        "is_psd": bool(is_psd),
        "definitely_not_psd": bool(not_psd),
        "inertia": list(res.signature) if res.signature else None,
        "blocker": res.blocker,
    }


def _arb_ops():
    from interval_backend import require_flint

    _, arb, _, _ = require_flint()
    return arb


def eigenvalue_bounds_from_two_moments(n: int, m1, m2) -> Dict[str, Any]:
    """Wolkowicz-Styan two-sided bounds on the extreme eigenvalues.

    Exact at ``n = 2``, where ``sqrt(n-1) = 1`` makes both inequalities
    equalities and the spectrum is fully determined by ``m1`` and ``m2``.
    """
    from interval_backend import interval_box

    arb = _arb_ops()
    if n < 1:
        raise ValueError("n must be >= 1")
    N = arb(int(n))
    a1 = m1 if hasattr(m1, "lower") else arb(str(m1))
    a2 = m2 if hasattr(m2, "lower") else arb(str(m2))
    mu = a1 / N
    var = a2 / N - mu * mu
    # The variance is a mean of squares about the mean, so it cannot be negative.
    # An enclosure can still dip below zero through rounding or through a wide
    # input; re-seating the lower end at 0 is outward-safe (the true value is in
    # the clamped interval too) and keeps the square root real.
    if var.lower() < 0:
        var = interval_box(0.0, max(0.0, float(var.upper())))
    s = var.sqrt()
    if n == 1:
        # A 1x1 matrix has one eigenvalue and it is the mean; there is no spread.
        lam_min_lo = lam_min_hi = lam_max_lo = lam_max_hi = mu
        exact = True
    else:
        k = (N - arb(1)).sqrt()
        lam_min_lo, lam_min_hi = mu - s * k, mu - s / k
        lam_max_lo, lam_max_hi = mu + s / k, mu + s * k
        exact = (n == 2)
    return {
        "n": int(n),
        "mean": [repr(float(mu.lower())), repr(float(mu.upper()))],
        "lambda_min": {"lo": repr(float(lam_min_lo.lower())),
                       "hi": repr(float(lam_min_hi.upper()))},
        "lambda_max": {"lo": repr(float(lam_max_lo.lower())),
                       "hi": repr(float(lam_max_hi.upper()))},
        "tight": exact,
        "method": "wolkowicz_styan_from_m1_m2",
        "note": ("at n = 2 the inequalities are equalities, so m1 and m2 determine "
                 "the spectrum exactly" if n == 2 else
                 "n = 1: the single eigenvalue is the mean" if n == 1 else
                 "two-sided bounds only; m1 and m2 do not determine the spectrum "
                 f"for n = {n}"),
    }


def spectrum_from_two_moments_n2(m1, m2) -> Dict[str, Any]:
    """The exact 2x2 spectrum: ``lambda = (m1 +- sqrt(2 m2 - m1^2)) / 2``."""
    arb = _arb_ops()
    a1 = m1 if hasattr(m1, "lower") else arb(str(m1))
    a2 = m2 if hasattr(m2, "lower") else arb(str(m2))
    disc = arb(2) * a2 - a1 * a1
    if disc.upper() < 0:
        return {"status": INSUFFICIENT,
                "blocker": "discriminant enclosure is negative; inputs are not the "
                           "moments of a real symmetric 2x2 matrix"}
    root = disc.sqrt()
    lo = (a1 - root) / arb(2)
    hi = (a1 + root) / arb(2)
    return {
        "status": CONCLUSIVE,
        "lambda_1": [repr(float(lo.lower())), repr(float(lo.upper()))],
        "lambda_2": [repr(float(hi.lower())), repr(float(hi.upper()))],
        "determined": True,
        "note": "m1 and m2 determine a 2x2 spectrum outright",
    }
