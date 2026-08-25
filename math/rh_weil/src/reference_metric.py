"""The reference Gram metric ``M`` for generalized spectral gaps (ENG-009 §WO-RH-57/58).

ENG-008 ended with a raw third minor around ``1e-14`` next to entries of order
``1e-1``, and §Mission is explicit about the trap: *raw determinant magnitude is
not a basis-invariant distance-to-failure observable*. The cure is to measure
``G`` against a mathematically defined positive reference form on the same
basis, because generalized eigenvalues of the pencil ``(G, M)`` are unchanged
under any simultaneous change of basis ``G -> S^T G S``, ``M -> S^T M S``.

The reference chosen here is the ordinary ``L^2`` Gram matrix of the basis on
the support interval:

    M_ij(L) = int_0^L h_i(x; L) h_j(x; L) dx.

Three facts make it the right choice, and this module derives all three
exactly rather than asserting them:

1. **It is exact.** Every basis element is a polynomial in ``x`` with
   ``L``-polynomial coefficients, so each entry is an exact polynomial in
   ``L`` with rational coefficients, derived from the same primitive table
   (``basis_algebra.BASIS_L_POLY``) the kernels come from.

2. **It is a monomial in ``L``.** Each basis element is homogeneous in
   ``(x, L)`` -- of degree 0, 1, 2, 3, 4 for one/q1/b/b3/b2 -- so
   ``M_ij(L) = m_ij * L^(d_i + d_j + 1)`` with ``m_ij`` a single rational.
   Equivalently ``M(L) = D(L)^T M(1) D(L)`` for the *real* diagonal
   ``D(L) = diag(L^(d_i + 1/2))``, which is invertible for every ``L > 0``.
   The homogeneity is checked at import time, not assumed.

3. **Its positivity is E0.** By that congruence, ``M(L)`` is positive definite
   for every ``L > 0`` iff the constant rational matrix ``M(1)`` is, and that
   is decided by exact rational Sylvester minors -- no interval arithmetic, no
   covers, no floating point anywhere. ``M(1)`` restricted to the even sector
   is the Hankel moment matrix of ``x(1-x)`` on ``[0, 1]``, which is why its
   minors are positive: it is the moment matrix of a measure with infinite
   support. The module still *computes* the minors rather than citing that.

Cross-parity entries vanish identically (an odd-about-midpoint integrand over
``[0, L]``), so ``M`` is parity-block-diagonal exactly as ``G`` is; the checks
verify that too.

No RH proof claim is made by this module. The metric is a *reference*, not a
result: it turns "the determinant is small" into the well-posed question "is
the pencil's smallest generalized eigenvalue small".
"""
from __future__ import annotations

from fractions import Fraction
from functools import lru_cache
from typing import Any, Dict, List, Sequence, Tuple

import basis_algebra
from basis_algebra import BASIS_L_POLY, BASIS_PARITY, LPoly, evaluate_l_poly

#: Homogeneity degree of each basis element in ``(x, L)`` jointly. Derived
#: below from the primitive table and verified, not free-standing data.
BASIS_DEGREE: Dict[str, int] = {}


def _derive_degrees() -> None:
    """Each element must be homogeneous: ``x``-power + ``L``-power constant."""
    for name, coeffs in BASIS_L_POLY.items():
        degrees = set()
        for xpow, lpoly in enumerate(coeffs):
            for lpow, c in lpoly.items():
                if c:
                    degrees.add(xpow + lpow)
        if len(degrees) != 1:
            raise AssertionError(
                f"basis element {name!r} is not homogeneous in (x, L): "
                f"degrees {sorted(degrees)}")
        BASIS_DEGREE[name] = degrees.pop()


_derive_degrees()


@lru_cache(maxsize=None)
def metric_l_poly(i: str, j: str) -> LPoly:
    """``M_ij(L) = int_0^L h_i h_j dx`` as an exact ``L``-polynomial."""
    for name in (i, j):
        if name not in BASIS_L_POLY:
            raise KeyError(f"unknown basis element {name!r}")
    pi, pj = BASIS_L_POLY[i], BASIS_L_POLY[j]
    out: LPoly = {}
    for xi, li in enumerate(pi):
        for xj, lj in enumerate(pj):
            xpow = xi + xj
            # int_0^L x^xpow dx = L^(xpow + 1) / (xpow + 1)
            for ki, ci in li.items():
                for kj, cj in lj.items():
                    c = ci * cj / (xpow + 1)
                    if not c:
                        continue
                    key = ki + kj + xpow + 1
                    got = out.get(key, Fraction(0)) + c
                    if got:
                        out[key] = got
                    else:
                        out.pop(key, None)
    return out


@lru_cache(maxsize=None)
def metric_monomial(i: str, j: str) -> Tuple[Fraction, int]:
    """``(m_ij, e_ij)`` with ``M_ij(L) = m_ij * L^e_ij`` -- exact.

    Homogeneity makes every entry a single monomial; cross-parity entries are
    exactly zero and come back as ``(0, e)`` with the degree the pairing would
    have had. A two-term polynomial here would mean the basis table changed
    incompatibly, and is an error rather than a fallback.
    """
    poly = metric_l_poly(i, j)
    expected = BASIS_DEGREE[i] + BASIS_DEGREE[j] + 1
    if not poly:
        return Fraction(0), expected
    if len(poly) != 1:
        raise AssertionError(
            f"M[{i},{j}] is not a monomial: powers {sorted(poly)}")
    (power, coeff), = poly.items()
    if power != expected:
        raise AssertionError(
            f"M[{i},{j}] has degree {power}, homogeneity predicts {expected}")
    return coeff, power


def metric_exact(i: str, j: str, L: Fraction) -> Fraction:
    """``M_ij`` at an exact rational ``L``."""
    coeff, power = metric_monomial(i, j)
    return coeff * Fraction(L) ** power


def metric_value(i: str, j: str, L: Any) -> Any:
    """``M_ij`` on the caller's carrier (float, mpmath, Arb ball)."""
    return evaluate_l_poly(metric_l_poly(i, j), L)


def metric_matrix_exact(basis: Sequence[str], L: Fraction) -> List[List[Fraction]]:
    return [[metric_exact(i, j, L) for j in basis] for i in basis]


def metric_matrix_over(basis: Sequence[str], L: Any) -> List[List[Any]]:
    """The metric on an interval (or any) carrier ``L``.

    Each entry is a monomial ``m L^e``, so the interval evaluation is a single
    power and a scalar multiple -- there is no dependency problem to manage.
    """
    return [[metric_value(i, j, L) for j in basis] for i in basis]


def unit_matrix(basis: Sequence[str]) -> List[List[Fraction]]:
    """``M(1)`` -- the constant rational matrix the congruence reduces to."""
    return [[metric_monomial(i, j)[0] for j in basis] for i in basis]


def _leading_minors_exact(m: List[List[Fraction]]) -> List[Fraction]:
    """Leading principal minors by exact fraction-free cofactor expansion."""
    out: List[Fraction] = []
    for k in range(1, len(m) + 1):
        out.append(_det_exact([row[:k] for row in m[:k]]))
    return out


def _det_exact(m: List[List[Fraction]]) -> Fraction:
    n = len(m)
    if n == 1:
        return m[0][0]
    total = Fraction(0)
    for col in range(n):
        minor = [row[:col] + row[col + 1:] for row in m[1:]]
        term = m[0][col] * _det_exact(minor)
        total += -term if col % 2 else term
    return total


def certify_positive_definite(basis: Sequence[str]) -> Dict[str, Any]:
    """E0 proof record that ``M(L)`` is PD on the block for every ``L > 0``.

    The argument, in the order the record states it:

    1. every entry is exactly ``m_ij L^(d_i + d_j + 1)`` (verified homogeneity);
    2. hence ``M(L) = D(L)^T M(1) D(L)`` with ``D(L) = diag(L^(d_i + 1/2))``,
       invertible for ``L > 0``;
    3. the leading principal minors of ``M(1)``, computed in exact rational
       arithmetic, are all positive, so ``M(1)`` is PD by Sylvester;
    4. positive definiteness transfers through the congruence
       (``AtlasRH.posDef_of_diagonal_congruence`` is the formal statement of
       this step for ``n = 3``).

    Nothing is floating point and nothing is an interval: the record's minors
    are exact fractions, serialized as strings.
    """
    m1 = unit_matrix(basis)
    minors = _leading_minors_exact(m1)
    if any(v <= 0 for v in minors):
        raise AssertionError(
            f"M(1) on {tuple(basis)} is not positive definite: minors {minors}")
    return {
        "reference_metric": "l2_gram_on_support",
        "definition": "M_ij(L) = int_0^L h_i(x; L) h_j(x; L) dx",
        "basis": list(basis),
        "homogeneity_degrees": {name: BASIS_DEGREE[name] for name in basis},
        "monomial_form": {
            f"{i}_{j}": {"coefficient": str(metric_monomial(i, j)[0]),
                         "l_power": metric_monomial(i, j)[1]}
            for i in basis for j in basis
        },
        "congruence": "M(L) = D(L)^T M(1) D(L), D(L) = diag(L^(d_i + 1/2))",
        "unit_matrix": [[str(v) for v in row] for row in m1],
        "unit_leading_minors": [str(v) for v in minors],
        "conclusion": "M(L) positive definite for every L > 0",
        "evidence_class": "E0",
        "arithmetic": "exact_rational",
    }


def rayleigh_quotient_exact(basis: Sequence[str],
                            g_entries: Dict[Tuple[str, str], Any],
                            v: Sequence[Fraction], L: Any) -> Any:
    """``v^T G v / v^T M v`` with ``v`` exact rational, on the caller's carrier.

    The workhorse of certified *upper* bounds on the generalized minimum: any
    vector gives one, and a rational vector keeps the numerator assembly exact
    in the coefficients.
    """
    num = None
    den = None
    for a, i in enumerate(basis):
        for b, j in enumerate(basis):
            c = Fraction(v[a]) * Fraction(v[b])
            if not c:
                continue
            key = (i, j) if (i, j) in g_entries else (j, i)
            gterm = g_entries[key] * c.numerator / c.denominator
            mterm = metric_value(i, j, L) * c.numerator / c.denominator
            num = gterm if num is None else num + gterm
            den = mterm if den is None else den + mterm
    return num / den
