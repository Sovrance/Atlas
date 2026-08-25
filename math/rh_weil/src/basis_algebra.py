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
from functools import lru_cache
from math import comb
from typing import Any, Dict, List, Optional, Sequence, Tuple

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
# Exact L-polynomials for the boundary quantities                              #
# --------------------------------------------------------------------------- #
# The derivative machinery needs `h(L; L)`, `(d_L h)(L; L)` and the total
# derivative `(d/dL)[h(L; L)]`. Each is a polynomial in `L`, and for most of the
# basis it is *identically zero*: `b`, `b3` and `b2` all vanish at `x = L`.
#
# Evaluating them on the carrier and letting the cancellation happen there is
# sound but catastrophically loose. On an `L`-ball of radius 1e-2, computing
# `b(L; L) = L*L - L*L` gives a ball of radius 2.2e-2 around zero instead of
# exact zero, and that width propagates through every derivative bound built on
# it -- measurably, when it was first tried: the degree-1 certified bound fell
# and the degree-2 one moved by 6%.
#
# So the cancellation happens here, in exact rational arithmetic, and the
# carrier only ever sees the simplified polynomial. An identically zero quantity
# comes back as an empty polynomial and evaluates to exact zero.
LPoly = Dict[int, Fraction]


def _lp_add(p: LPoly, k: int, c: Fraction) -> None:
    if not c:
        return
    got = p.get(k, Fraction(0)) + c
    if got:
        p[k] = got
    else:
        p.pop(k, None)


def differentiate_l(coeffs: Sequence[LPoly]) -> Tuple[LPoly, ...]:
    """``d/dL`` of each ``x``-coefficient, exactly."""
    out: List[LPoly] = []
    for lc in coeffs:
        d: LPoly = {}
        for k, c in lc.items():
            if k:
                _lp_add(d, k - 1, c * k)
        out.append(d)
    return tuple(out)


@lru_cache(maxsize=None)
def basis_dL_coeffs(name: str) -> Tuple[LPoly, ...]:
    """``d_L h`` as ``x``-coefficients, each an exact polynomial in ``L``."""
    return differentiate_l(BASIS_L_POLY[name])


@lru_cache(maxsize=None)
def basis_d2L_coeffs(name: str) -> Tuple[LPoly, ...]:
    """``d^2_L h``, likewise."""
    return differentiate_l(basis_dL_coeffs(name))


def _collapse_at_x_equals_L(coeffs: Sequence[LPoly]) -> LPoly:
    """``p(L; L)`` -- substitute ``x = L`` and collect, exactly."""
    out: LPoly = {}
    for n, lc in enumerate(coeffs):
        for k, c in lc.items():
            _lp_add(out, k + n, c)
    return out


@lru_cache(maxsize=None)
def endpoint_poly(name: str) -> LPoly:
    """``h(L; L)`` as an exact polynomial in ``L``. Empty means identically 0."""
    return _collapse_at_x_equals_L(BASIS_L_POLY[name])


@lru_cache(maxsize=None)
def endpoint_dL_inner_poly(name: str) -> LPoly:
    """``(d_L h)(L; L)`` -- the moving-integrand term, exactly."""
    return _collapse_at_x_equals_L(basis_dL_coeffs(name))


@lru_cache(maxsize=None)
def endpoint_total_dL_poly(name: str) -> LPoly:
    """``(d/dL)[h(L; L)] = (d_x h + d_L h)(L; L)``, exactly.

    Every element of the current basis makes this identically zero -- ``one`` is
    constant and the other four vanish at ``x = L`` for every ``L``, so their
    endpoint value is the zero polynomial and its derivative with it. The general
    form is computed anyway, because that is an accident of this basis rather
    than a fact about the construction.
    """
    src = BASIS_L_POLY[name]
    dx: List[LPoly] = []
    for n, lc in enumerate(src):
        if n == 0:
            continue
        dx.append({k: c * n for k, c in lc.items()})
    out = _collapse_at_x_equals_L(tuple(dx)) if dx else {}
    for k, c in endpoint_dL_inner_poly(name).items():
        _lp_add(out, k, c)
    return out


def evaluate_l_poly(poly: LPoly, L: Any) -> Any:
    """Evaluate an exact ``L``-polynomial on the caller's carrier.

    An empty polynomial is exactly zero, and returns the carrier's zero rather
    than a ball that merely happens to contain it.
    """
    zero = L * 0
    if not poly:
        return zero
    out = zero
    for k, c in sorted(poly.items()):
        term = (zero + 1) if k == 0 else L ** k
        out = out + term * c.numerator / c.denominator
    return out


def basis_coeffs_dL_on(name: str, L: Any) -> Tuple[Any, ...]:
    return tuple(evaluate_l_poly(lc, L) for lc in basis_dL_coeffs(name))


def basis_coeffs_d2L_on(name: str, L: Any) -> Tuple[Any, ...]:
    return tuple(evaluate_l_poly(lc, L) for lc in basis_d2L_coeffs(name))


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
# The (L - a)^m factorization                                                  #
# --------------------------------------------------------------------------- #
# Every kernel integrates over `[0, L - a]`, so `(L - a)` divides it -- the old
# hand-written closed forms all displayed that factor, e.g.
#
#     K_b3b3 = (L-a)^3 (L^4 + 3L^3 a - 15L^2 a^2 - 18L a^3 - 6a^4) / 420.
#
# On an exact carrier the factored and expanded forms are the same number. On an
# interval carrier they are not, and the difference is large: expanding puts `L`
# into every coefficient independently, so their widths stop cancelling. The
# first version of this module evaluated the expanded form and the prime block's
# enclosure widened by 3x for `K_q1q1` and 48x for `K_b3b3`, which cost the
# degree-3 determinant bound 26%.
#
# So the factorization is recovered here, automatically, by synthetic division in
# `a` at the root `a = L`, repeated while the remainder is exactly the zero
# polynomial. That is a statement about the exact rational polynomial, decided in
# exact arithmetic, and it costs nothing: `kernel_value` then evaluates
# `(L - a)^m * Q(a; L)` with `L - a` formed once.
def _synthetic_divide_by_a_minus_L(
    coeffs: Sequence[LPoly],
) -> Optional[List[LPoly]]:
    """Divide a polynomial in ``a`` by ``(a - L)``; ``None`` if it does not divide.

    Coefficients are ascending in ``a`` and are themselves exact polynomials in
    ``L``. Horner's scheme, run in exact rational arithmetic, so "does not
    divide" is a decision rather than a tolerance.
    """
    n = len(coeffs) - 1
    if n < 1:
        return None
    out: List[LPoly] = [dict() for _ in range(n)]
    carry: LPoly = dict(coeffs[n])
    for k in range(n - 1, -1, -1):
        out[k] = dict(carry)
        shifted: LPoly = {}
        for e, c in carry.items():
            _lp_add(shifted, e + 1, c)
        carry = dict(coeffs[k])
        for e, c in shifted.items():
            _lp_add(carry, e, c)
    return out if not carry else None


@lru_cache(maxsize=None)
def kernel_factored(i: str, j: str) -> Tuple[int, Tuple[LPoly, ...]]:
    """``K_ij = (L - a)^m * Q(a; L)``, as ``(m, Q-coefficients-in-a)``.

    ``m`` is the exact multiplicity of ``a = L`` as a root, found by repeated
    synthetic division. The sign is folded into ``Q``, since dividing by
    ``(a - L)`` differs from dividing by ``(L - a)`` by one factor of ``-1``.
    """
    poly = kernel_bipoly(i, j)
    deg = kernel_degree_in_a(i, j)
    coeffs: List[LPoly] = [dict() for _ in range(deg + 1)]
    for (ap, lp), c in poly.items():
        _lp_add(coeffs[ap], lp, c)
    m = 0
    while True:
        divided = _synthetic_divide_by_a_minus_L(coeffs)
        if divided is None:
            break
        coeffs = [{e: -c for e, c in lc.items()} for lc in divided]
        m += 1
    return m, tuple(coeffs)


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
    """``K_ij(a; L)``, evaluated in the factored form ``(L - a)^m Q(a; L)``.

    The factorization is not cosmetic on an interval carrier: see
    :func:`kernel_factored`. On an exact carrier it makes no difference, and
    ``tests/test_kernel_algebra.py`` checks the two agree exactly.
    """
    m, qcoeffs = kernel_factored(i, j)
    out = evaluate_l_poly(qcoeffs[-1], L)
    for lc in reversed(qcoeffs[:-1]):
        out = out * a + evaluate_l_poly(lc, L)
    if m:
        d = L - a
        out = out * d ** m
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
