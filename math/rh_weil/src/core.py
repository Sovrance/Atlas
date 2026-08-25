"""Exact/algebraic core for Atlas RH/Weil notebook integration.

No RH claim is made by this module. It contains finite-block identities only.
Stdlib-only by design so exact regression tests can run everywhere.
"""
from __future__ import annotations

from math import comb, sqrt

# Explicit assembly convention (normalization audit).
NORMALIZATION = "G = G0 - Gp + Ginf"
CLAIM_BOUNDARY = "finite_block_only_no_rh_proof"


def overlap_c(i: int, j: int, a, L):
    """C_ij(a;L) for numeric/symbolic-like scalar types supporting arithmetic."""
    d = L - a
    total = 0
    for r in range(j + 1):
        total += comb(j, r) * (a ** (j - r)) * (d ** (i + r + 1)) / (i + r + 1)
    return total


def kernel_ij(i: int, j: int, a, L):
    """K_ij = C_ij + C_ji."""
    return overlap_c(i, j, a, L) + overlap_c(j, i, a, L)


def kernel_00(a, L):
    return 2 * (L - a)


def kernel_01(a, L):
    return L * (L - a)


def kernel_11(a, L):
    d = L - a
    return d * d * (2 * L + a) / 3


def kernel_q1q1(a, L):
    d = L - a
    return d * (L * L - 2 * L * a - 2 * a * a) / 6


def kernel_q1b3(a, L):
    """``K_q1b3(a; L)``, ENG-006 §7 -- verified against SymPy in the test suite."""
    d = L - a
    return d * d * (L**3 + 2 * L**2 * a - 12 * L * a**2 - 6 * a**3) / 60


def kernel_b3b3(a, L):
    """``K_b3b3(a; L)``, ENG-006 §7 -- verified against SymPy in the test suite."""
    d = L - a
    return (d * d * d
            * (L**4 + 3 * L**3 * a - 15 * L**2 * a**2 - 18 * L * a**3 - 6 * a**4)
            / 420)


def kernel_0b(a, L):
    d = L - a
    return d * d * (L + 2 * a) / 3


def kernel_bb(a, L):
    d = L - a
    return d**3 * (L * L + 3 * L * a + a * a) / 15


def kernel_bubble_det(a, L):
    d = L - a
    return d**4 * (L * L - 2 * L * a - 14 * a * a) / 45


def scalar_curvature(L):
    """W00''(L) inside a prime-power cell, r = e^L."""
    from math import exp

    r = exp(L)
    return 2 * (r**3 - r - 1) / (sqrt(r) * (r * r - 1))


def scalar_curvature_r(r):
    """W00'' as a function of r = e^L (same algebraic expression)."""
    return 2 * (r**3 - r - 1) / (sqrt(r) * (r * r - 1))


def von_mangoldt_jump(q: int, p: int):
    """Downward derivative jump magnitude at q = p^k: 2 Λ(q)/√q with Λ=log p."""
    from math import log

    return 2 * log(p) / sqrt(q)


def q1_sign_threshold():
    return (sqrt(3) - 1) / 2


def bubble_det_threshold():
    return (sqrt(15) - 1) / 14


def degree2_raw_det(g00, odd_pivot, even_det, L):
    """D2 = E2 + L^2 G00 O1 (parity factorization)."""
    return even_det + L * L * g00 * odd_pivot


def degree2_full_det(odd_pivot, even_det):
    """det G[{1,x,x^2}] = O1 * E2."""
    return odd_pivot * even_det


def midpoint_reflection_q1(L, a):
    """q1 = x - L/2; K_q1q1 via monomial transform (FORMULAS / WO-RH-03)."""
    k00 = kernel_00(a, L)
    k01 = kernel_01(a, L)
    k11 = kernel_11(a, L)
    return k11 - L * k01 + (L * L / 4) * k00
