"""An independent assembly of the 3x3 even block (ATLAS-RH-ENG-008 §WO-RH-48).

A cross-check is only worth running if it can fail. This module therefore
implements ``G = G0 - Gp + Ginf`` for ``{1, b, b^2}`` **without importing any of
the code it checks** -- not :mod:`basis_algebra`, not :mod:`pole`, not
:mod:`weil_entries`, not :mod:`archimedean_realspace`. A test asserts that, so
the independence cannot quietly lapse into a re-run of the same code.

Where the rigorous path uses exact bivariate integer arithmetic for the overlap
kernels and Arb for the quadrature, this one derives the kernels with SymPy
straight from the definition

    K_ij(a; L) = int_0^{L-a} [ h_i(x) h_j(x+a) + h_j(x) h_i(x+a) ] dx

and integrates in mpmath. Different algebra, different arithmetic, same object.

Evidence class is **E3**. mpmath never certifies anything in this program
(ENG-004 §5), and nothing here is promoted. Its job is regression: if the
rigorous assembly ever drifts, the two stop agreeing.

The archimedean term uses the ENG-005 real-space form

    Ginf_ij(L) = (K(0)/2) h_+(0) + int_0^L [K(0) - K(u)] w(u) du + K(0) S(L),
    w(u) = e^{-u/2} / (1 - e^{-2u}),   S(L) = sum_n e^{-(2n + 1/2) L} / (2n + 1/2)

which is the same *formula* as the rigorous route -- it has to be, or the two
would be computing different objects and agreement would mean nothing. What is
independent is every step of getting from that formula to a number.

No RH proof claim is made by this module.
"""
from __future__ import annotations

import math
from functools import lru_cache
from typing import Any, Dict, List, Sequence, Tuple

CLAIM_SCOPE = "finite_dimensional_weil_compression"
EVIDENCE_CLASS = "E3"

BASIS: Tuple[str, ...] = ("one", "b", "b2")

CELL: Tuple[float, float] = (math.log(3.0), math.log(4.0))


def require_mpmath():
    import mpmath

    return mpmath


def require_sympy():
    import sympy

    return sympy


# --------------------------------------------------------------------------- #
# Kernels, derived symbolically from the definition                            #
# --------------------------------------------------------------------------- #
def _basis_expr(name: str, sp, x, L):
    if name == "one":
        return sp.Integer(1)
    if name == "b":
        return x * (L - x)
    if name == "b2":
        return (x * (L - x)) ** 2
    if name == "bcube":  # ENG-010 §WO-RH-66: same definition route as the rest
        return (x * (L - x)) ** 3
    raise KeyError(f"unknown basis element {name!r}")


@lru_cache(maxsize=None)
def _kernel_lambdas(i: str, j: str):
    """``K_ij`` as coefficients in ``a``, each lambdified as a function of ``L``.

    Kept as compiled callables rather than round-tripped through strings: a
    symbol re-created by ``sympify`` carries no assumptions and therefore does
    not match the one the expression was built with, so substitution silently
    fails and leaves ``L`` in the result. Lambdifying against the original
    symbol cannot go wrong that way, and it is far faster besides.
    """
    sp = require_sympy()
    x, a, L = sp.symbols("x a L", positive=True)
    hi = _basis_expr(i, sp, x, L)
    hj = _basis_expr(j, sp, x, L)
    integrand = sp.expand(hi * hj.subs(x, x + a) + hj * hi.subs(x, x + a))
    K = sp.expand(sp.integrate(integrand, (x, 0, L - a)))
    coeffs = list(reversed(sp.Poly(K, a).all_coeffs()))
    return tuple(sp.lambdify(L, c, "mpmath") for c in coeffs)


@lru_cache(maxsize=None)
def _basis_lambda(name: str):
    sp = require_sympy()
    x, L = sp.symbols("x L", positive=True)
    return sp.lambdify((x, L), _basis_expr(name, sp, x, L), "mpmath")


def kernel_coeffs(i: str, j: str, L_val, mp) -> List[Any]:
    """``K_ij(a; L)`` as coefficients in ``a``, ascending."""
    L_m = mp.mpf(str(L_val))
    return [mp.mpf(f(L_m)) for f in _kernel_lambdas(i, j)]


def kernel(i: str, j: str, a_val, L_val, mp):
    """``K_ij(a; L)`` by Horner in ``a``."""
    coeffs = kernel_coeffs(i, j, L_val, mp)
    out = coeffs[-1]
    for c in reversed(coeffs[:-1]):
        out = out * a_val + c
    return out


def kernel_difference_over_u(i: str, j: str, u, L_val, mp):
    """``[K(0) - K(u)] / u``, as a polynomial -- no 0/0 at the origin."""
    coeffs = kernel_coeffs(i, j, L_val, mp)
    total = mp.mpf(0)
    power = mp.mpf(1)
    for c in coeffs[1:]:
        total += c * power
        power = power * u
    return -total


# --------------------------------------------------------------------------- #
# The three blocks                                                             #
# --------------------------------------------------------------------------- #
def prime_powers_below(L_value: float) -> List[Tuple[int, int]]:
    cap = int(math.floor(math.exp(L_value)))
    out = []
    for p in range(2, cap + 1):
        if any(p % d == 0 for d in range(2, int(p ** 0.5) + 1)):
            continue
        q = p
        while q <= cap and math.log(q) < L_value:
            out.append((q, p))
            q *= p
    return sorted(out)


def pole_entry(i: str, j: str, L_val, mp):
    """``G0_ij = E_i^+ E_j^- + E_i^- E_j^+`` by direct numerical quadrature.

    Quadrature rather than a closed form on purpose: the rigorous route
    evaluates a closed form for ``E^±``, so integrating numerically here checks
    that closed form instead of repeating it.
    """
    L_m = mp.mpf(str(L_val))

    def lap(name, sign):
        f = _basis_lambda(name)
        return mp.quad(lambda t: f(t, L_m) * mp.e ** (sign * t / 2), [0, L_m])

    return lap(i, 1) * lap(j, -1) + lap(i, -1) * lap(j, 1)


def prime_entry(i: str, j: str, L_val, mp, prime_powers=None):
    if prime_powers is None:
        prime_powers = prime_powers_below(float(L_val))
    total = mp.mpf(0)
    for q, p in prime_powers:
        total += (mp.log(p) / mp.sqrt(q)) * kernel(i, j, mp.log(q), L_val, mp)
    return total


def arch_entry(i: str, j: str, L_val, mp, *, series_terms: int = 200):
    """``Ginf_ij`` by the ENG-005 real-space form, in mpmath."""
    L_m = mp.mpf(str(L_val))
    K0 = kernel_coeffs(i, j, L_val, mp)[0]

    def w(u):
        return mp.e ** (-u / 2) / (1 - mp.e ** (-2 * u))

    def integrand(u):
        if u == 0:
            return mp.mpf(0)
        return kernel_difference_over_u(i, j, u, L_val, mp) * u * w(u)

    integral = mp.quad(integrand, [0, L_m / 4, L_m / 2, L_m])
    S = mp.nsum(lambda n: mp.e ** (-(2 * n + mp.mpf(1) / 2) * L_m)
                / (2 * n + mp.mpf(1) / 2), [0, series_terms])
    h_plus_0 = mp.re(mp.digamma(mp.mpf(1) / 4)) - mp.log(mp.pi)
    return K0 / 2 * h_plus_0 + integral + K0 * S


def gram_entry(i: str, j: str, L_val, mp, prime_powers=None):
    return (pole_entry(i, j, L_val, mp)
            - prime_entry(i, j, L_val, mp, prime_powers)
            + arch_entry(i, j, L_val, mp))


def gram_matrix(L_val, mp, *, dps: int = 40) -> List[List[Any]]:
    previous = mp.mp.dps
    mp.mp.dps = dps
    try:
        pp = prime_powers_below(float(L_val))
        n = len(BASIS)
        out = [[None] * n for _ in range(n)]
        for a in range(n):
            for b in range(a, n):
                v = gram_entry(BASIS[a], BASIS[b], L_val, mp, pp)
                out[a][b] = v
                out[b][a] = v
        return out
    finally:
        mp.mp.dps = previous


def leading_minors(M) -> List[Any]:
    a, b, c = M[0]
    _, d, e = M[1]
    _, _, f = M[2]
    return [a,
            a * d - b * b,
            a * (d * f - e * e) - b * (b * f - e * c) + c * (b * e - d * c)]
