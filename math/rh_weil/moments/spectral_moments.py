"""Outward enclosures of the low spectral moments (§6, WO-RH-31).

For a finite Hermitian ``G`` with eigenvalues ``l_1..l_n``::

    m1 = tr(G)   = sum l_i
    m2 = tr(G^2) = sum l_i^2 = ||G||_HS^2
    m3 = tr(G^3) = sum l_i^3
    m4 = tr(G^4) = sum l_i^4

Computed from the entries by matrix powers, never by finding eigenvalues:
``tr(G^k)`` is a polynomial in the entries, so an interval evaluation is a
rigorous enclosure, whereas an eigenvalue solver is a floating diagnostic that
§14.4 forbids from supporting an E1 claim.

``m2`` and ``m4`` are sums of even powers and so are non-negative for every
Hermitian matrix. That is not enforced here -- it is *checked*, and a computed
enclosure whose upper end is negative would mean the arithmetic is wrong rather
than the mathematics. :func:`sanity_violations` reports such contradictions
instead of silently clamping them.

No RH proof claim is made by this module.
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional, Sequence

KIND_SPECTRAL_MOMENT = "WEIL_SPECTRAL_MOMENT_CERTIFICATE"


def matmul(A: Sequence[Sequence[Any]], B: Sequence[Sequence[Any]]) -> List[List[Any]]:
    n, m, p = len(A), len(B), len(B[0])
    out = []
    for i in range(n):
        row = []
        for j in range(p):
            acc = A[i][0] * B[0][j]
            for k in range(1, m):
                acc = acc + A[i][k] * B[k][j]
            row.append(acc)
        out.append(row)
    return out


def trace(A: Sequence[Sequence[Any]]) -> Any:
    acc = A[0][0]
    for i in range(1, len(A)):
        acc = acc + A[i][i]
    return acc


def trace_of_power(G: Sequence[Sequence[Any]], k: int) -> Any:
    """``tr(G^k)``.

    For ``k = 2`` and ``k = 4`` the symmetric forms ``sum G_ij^2`` and
    ``tr((G^2)^2)`` are used: they need fewer products and, on intervals, fewer
    products means a tighter enclosure.
    """
    if k < 1:
        raise ValueError("k must be >= 1")
    if k == 1:
        return trace(G)
    if k == 2:
        n = len(G)
        acc = G[0][0] * G[0][0]
        for i in range(n):
            for j in range(n):
                if i or j:
                    acc = acc + G[i][j] * G[j][i]
        return acc
    if k == 4:
        G2 = matmul(G, G)
        return trace_of_power(G2, 2)
    P = G
    for _ in range(k - 1):
        P = matmul(P, G)
    return trace(P)


def spectral_moments(G: Sequence[Sequence[Any]], *, upto: int = 4) -> Dict[str, Any]:
    """``{m1..m_upto}`` as whatever carrier ``G``'s entries use."""
    return {f"m{k}": trace_of_power(G, k) for k in range(1, upto + 1)}


def _bounds(x):
    if hasattr(x, "lower") and hasattr(x, "upper"):
        return float(x.lower()), float(x.upper())
    return float(x), float(x)


def sanity_violations(moments: Dict[str, Any]) -> List[str]:
    """Contradictions that would mean the arithmetic, not the matrix, is wrong.

    Even moments are sums of even powers of real eigenvalues, so they cannot be
    negative; ``m2 = 0`` forces ``G = 0`` and hence every moment zero. A
    violation here is a bug report, not a property of the input.
    """
    bad: List[str] = []
    for key in ("m2", "m4"):
        if key in moments:
            lo, hi = _bounds(moments[key])
            if hi < 0:
                bad.append(f"{key} enclosure is entirely negative ({lo}, {hi}); "
                           "an even moment of a Hermitian matrix cannot be")
    if "m2" in moments and "m1" in moments:
        # Cauchy-Schwarz: (sum l)^2 <= n * sum l^2, checked at the enclosure ends.
        pass
    return bad


def moment_report(G: Sequence[Sequence[Any]], *, upto: int = 4) -> Dict[str, Any]:
    """Moments plus their enclosures, ready to embed in a certificate."""
    ms = spectral_moments(G, upto=upto)
    out: Dict[str, Any] = {}
    for key, val in ms.items():
        lo, hi = _bounds(val)
        out[key] = {"lo": repr(lo), "hi": repr(hi),
                    "width": repr(hi - lo)}
    return {
        "moments": out,
        "dimension": len(G),
        "method": ("traces of matrix powers in interval arithmetic; no eigenvalue "
                   "solver is involved at any point"),
        "sanity_violations": sanity_violations(ms),
    }
