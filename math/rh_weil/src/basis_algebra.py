"""Exact prime-overlap kernels, derived rather than tabulated (ENG-008 §WO-RH-48).

The prime block needs, for each pair of basis elements,

    K_ij(a; L) = int_0^{L-a} [ h_i(x) h_j(x+a) + h_j(x) h_i(x+a) ] dx

in three forms: as a closed form in ``a``, as a coefficient list in ``a`` (the
real-space archimedean route consumes that), and as the same list for
``d/dL K``. Until ENG-008 all three were hand-written tables, one entry per
pair, in three different modules -- so adding a basis element meant deriving six
polynomial identities by hand and pasting them into three places, and forgetting
one produced a ``KeyError`` at best.

That is the pattern §WO-RH-49 removed from the derivative machinery, for the
same reason: ENG-005 lost time to a hand-derived coefficient written as ``L/2``
instead of ``L``. So this module computes all three forms from the one thing
that is genuinely primitive -- the monomial coefficients of each basis element,
as exact polynomials in ``L`` -- using exact rational arithmetic at import time.

Every basis element here is a polynomial in ``x`` whose coefficients are
polynomials in ``L`` with rational coefficients, so ``K_ij`` is a polynomial in
``(a, L)`` with rational coefficients. That bivariate table is computed once,
exactly, and the three runtime forms are read off it. Nothing is rounded until
the caller's carrier is applied.

``tests/test_kernel_algebra.py`` checks the derived kernels against the six
hand-written closed forms in :mod:`core` that ENG-005 and ENG-006 verified
against SymPy, and against direct symbolic integration. Two independent routes,
neither trusting the other.

No RH proof claim is made by this module.
"""
from __future__ import annotations

from fractions import Fraction
from math import comb
from typing import Any, Dict, List, Sequence, Tuple

#: Each basis element as ``x``-coefficients, every one a polynomial in ``L``
#: given as ``{L_power: rational}``. This is the single primitive: everything
#: else in this module is derived from it.
#:
#:   one = 1
#:   q1  = x - L/2
#:   b   = L x - x^2
#:   b3  = -(L^2/2) x + (3L/2) x^2 - x^3
#:   b2  = L^2 x^2 - 2L x^3 + x^4
BASIS_L_POLY: Dict[str, Tuple[Dict[int, Fraction], ...]] = {
    "one": ({0: Fraction(1)},),
    "q1": ({1: Fraction(-1, 2)}, {0: Fraction(1)}),
    "b": ({}, {1: Fraction(1)}, {0: Fraction(-1)}),
    "b3": ({}, {2: Fraction(-1, 2)}, {1: Fraction(3, 2)}, {0: Fraction(-1)}),
    "b2": ({}, {}, {2: Fraction(1)}, {1: Fraction(-2)}, {0: Fraction(1)}),
}

BASIS_NAMES: Tuple[str, ...] = tuple(BASIS_L_POLY)

#: Parity about ``x = L/2``, checked exactly in the tests rather than asserted.
BASIS_PARITY: Dict[str, str] = {
    "one": "even", "q1": "odd", "b": "even", "b3": "odd", "b2": "even",
}


# --------------------------------------------------------------------------- #
# A minimal exact bivariate polynomial in (a, L)                               #
# --------------------------------------------------------------------------- #
BiPoly = Dict[Tuple[int, int], Fraction]  # (a_power, L_power) -> coefficient


def _bp_add(p: BiPoly, key: Tuple[int, int], c: Fraction) -> None:
    if not c:
        return
    got = p.get(key, Fraction(0)) + c
    if got:
        p[key] = got
    else:
        p.pop(key, None)


def _bp_mul(p: BiPoly, q: BiPoly) -> BiPoly:
    out: BiPoly = {}
    for (a1, l1), c1 in p.items():
        for (a2, l2), c2 in q.items():
            _bp_add(out, (a1 + a2, l1 + l2), c1 * c2)
    return out


def _bp_pow(p: BiPoly, n: int) -> BiPoly:
    out: BiPoly = {(0, 0): Fraction(1)}
    for _ in range(n):
        out = _bp_mul(out, p)
    return out


def _l_poly_to_bipoly(coeffs: Dict[int, Fraction]) -> BiPoly:
    return {(0, k): Fraction(v) for k, v in coeffs.items() if v}


# --------------------------------------------------------------------------- #
# The kernel, as an exact polynomial in (a, L)                                 #
# --------------------------------------------------------------------------- #
def _shifted_x_coeffs(name: str) -> List[BiPoly]:
    """``h(x + a)`` re-expanded in powers of ``x``.

    ``(x + a)^n = sum_r C(n, r) x^r a^{n-r}``, so the coefficient of ``x^r``
    collects an ``a^{n-r}`` from every higher monomial.
    """
    src = BASIS_L_POLY[name]
    out: List[BiPoly] = [dict() for _ in src]
    for n, lc in enumerate(src):
        base = _l_poly_to_bipoly(lc)
        for r in range(n + 1):
            factor = Fraction(comb(n, r))
            for (ap, lp), c in base.items():
                _bp_add(out[r], (ap + (n - r), lp), c * factor)
    return out


def _plain_x_coeffs(name: str) -> List[BiPoly]:
    return [_l_poly_to_bipoly(lc) for lc in BASIS_L_POLY[name]]


def _integrate_to_L_minus_a(prod: List[BiPoly]) -> BiPoly:
    """``int_0^{L-a} p(x) dx`` for ``p`` given as ``x``-coefficients."""
    dm = {(1, 0): Fraction(-1), (0, 1): Fraction(1)}  # (L - a)
    out: BiPoly = {}
    for k, ck in enumerate(prod):
        if not ck:
            continue
        power = _bp_pow(dm, k + 1)
        scaled = {key: c / (k + 1) for key, c in power.items()}
        for key, c in _bp_mul(ck, scaled).items():
            _bp_add(out, key, c)
    return out


def _multiply_x_series(p: List[BiPoly], q: List[BiPoly]) -> List[BiPoly]:
    out: List[BiPoly] = [dict() for _ in range(len(p) + len(q) - 1)]
    for m, pm in enumerate(p):
        if not pm:
            continue
        for n, qn in enumerate(q):
            if not qn:
                continue
            for key, c in _bp_mul(pm, qn).items():
                _bp_add(out[m + n], key, c)
    return out


_KERNEL_CACHE: Dict[Tuple[str, str], BiPoly] = {}


def kernel_bipoly(i: str, j: str) -> BiPoly:
    """``K_ij(a; L)`` as exact ``(a_power, L_power) -> Fraction``."""
    key = (i, j) if (i, j) in _KERNEL_CACHE else (j, i)
    if key in _KERNEL_CACHE:
        return _KERNEL_CACHE[key]
    for name in (i, j):
        if name not in BASIS_L_POLY:
            raise KeyError(f"unknown basis element {name!r}")
    hi, hj = _plain_x_coeffs(i), _plain_x_coeffs(j)
    si, sj = _shifted_x_coeffs(i), _shifted_x_coeffs(j)
    integrand = _multiply_x_series(hi, sj)
    for idx, term in enumerate(_multiply_x_series(hj, si)):
        for k, c in term.items():
            _bp_add(integrand[idx], k, c)
    out = _integrate_to_L_minus_a(integrand)
    _KERNEL_CACHE[(i, j)] = out
    return out


def kernel_degree_in_a(i: str, j: str) -> int:
    poly = kernel_bipoly(i, j)
    return max((ap for ap, _ in poly), default=0)


# --------------------------------------------------------------------------- #
# The three runtime forms                                                      #
# --------------------------------------------------------------------------- #
def _eval_L(terms: Sequence[Tuple[int, Fraction]], L: Any, zero: Any) -> Any:
    """``sum c * L^k`` on the caller's carrier, with exact rational scaling."""
    out = zero
    for k, c in terms:
        term = (zero + 1) if k == 0 else L ** k
        out = out + term * c.numerator / c.denominator
    return out


def kernel_coeffs_in_a(i: str, j: str, L: Any) -> List[Any]:
    """``K_ij(a; L)`` as coefficients in ``a``, ascending, on ``L``'s carrier."""
    poly = kernel_bipoly(i, j)
    deg = kernel_degree_in_a(i, j)
    zero = L * 0
    buckets: List[List[Tuple[int, Fraction]]] = [[] for _ in range(deg + 1)]
    for (ap, lp), c in poly.items():
        buckets[ap].append((lp, c))
    return [_eval_L(sorted(b), L, zero) for b in buckets]


def kernel_dL_coeffs_in_a(i: str, j: str, L: Any) -> List[Any]:
    """``d/dL K_ij(a; L)`` as coefficients in ``a``, ascending."""
    poly = kernel_bipoly(i, j)
    deg = kernel_degree_in_a(i, j)
    zero = L * 0
    buckets: List[List[Tuple[int, Fraction]]] = [[] for _ in range(deg + 1)]
    for (ap, lp), c in poly.items():
        if lp == 0:
            continue
        buckets[ap].append((lp - 1, c * lp))
    return [_eval_L(sorted(b), L, zero) for b in buckets]


def kernel_value(i: str, j: str, a: Any, L: Any) -> Any:
    """``K_ij(a; L)`` on the caller's carrier, by Horner in ``a``."""
    coeffs = kernel_coeffs_in_a(i, j, L)
    out = coeffs[-1]
    for c in reversed(coeffs[:-1]):
        out = out * a + c
    return out


def kernel_at_zero(i: str, j: str, L: Any) -> Any:
    """``K_ij(0; L)``, the value the real-space route subtracts against."""
    return kernel_coeffs_in_a(i, j, L)[0]


def kernel_exact(i: str, j: str, a: Fraction, L: Fraction) -> Fraction:
    """``K_ij(a; L)`` in exact rational arithmetic, for the E0 tests."""
    total = Fraction(0)
    for (ap, lp), c in kernel_bipoly(i, j).items():
        total += c * Fraction(a) ** ap * Fraction(L) ** lp
    return total
