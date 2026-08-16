"""True finite-cutoff Weil Gram (WO-RH-12). Not the E3 energy probe.

Assembly: G = G0 - Gp + Ginf_T
Ginf_T[i,j] = (1/π) ∫_0^T h_+(t) Re(conj(H_i(t;L)) H_j(t;L)) dt
"""
from __future__ import annotations

from typing import Any, Dict

import core
from archimedean import h_plus
from interval_backend import FlintUnavailable, require_flint, set_precision_bits


NORMALIZATION = core.NORMALIZATION


def _A_B(z, arb):
    """sinc and bubble shape factors; z-ball must not contain 0 for power forms."""
    if z.contains(0) or (hasattr(z, "rad") and z.rad() >= abs(z.mid())):
        # Caller should use analytic limits at t=0.
        raise ValueError("z-ball contains 0; use t=0 analytic limits")
    A = z.sin() / z
    B = (z.sin() - z * z.cos()) / (z**3)
    return A, B


def stable_products_even(t, L, arb):
    """Return (|H0|^2, Re(conj(H0)Hb), |Hb|^2) as Arb balls."""
    t_a = arb(t)
    L_a = arb(L)
    if t_a.contains(0) or abs(t_a.mid()) <= t_a.rad():
        # t≈0 limits
        h0sq = L_a**2
        re0b = L_a**4 / 6  # Re(conj(L)* (L^3/6)) = L^4/6
        hbsq = (L_a**3 / 6) ** 2
        return h0sq, re0b, hbsq
    z = L_a * t_a / 2
    A, B = _A_B(z, arb)
    h0sq = (L_a**2) * (A**2)
    re0b = (L_a**4 / 2) * A * B
    hbsq = (L_a**6 / 4) * (B**2)
    return h0sq, re0b, hbsq


def prime_powers_below(L, arb):
    """Exact (q, p, log q) for prime powers with log q < L (midpoint compare)."""
    import math

    L_mid = float(L.mid()) if hasattr(L, "mid") else float(L)
    c = max(2, int(math.floor(math.exp(L_mid) + 1e-12)))
    is_prime = [True] * (c + 1)
    is_prime[0] = is_prime[1] = False
    for i in range(2, int(c**0.5) + 1):
        if is_prime[i]:
            for j in range(i * i, c + 1, i):
                is_prime[j] = False
    primes = [i for i in range(2, c + 1) if is_prime[i]]
    out = []
    for p in primes:
        q = p
        while q <= c:
            lq = arb(q).log()
            if lq < L:
                out.append((q, p, lq))
            if q > c // p:
                break
            q *= p
    return out


def gp_even_block(L, arb):
    """Prime block for basis {1,b} using K00, K0b, Kbb."""
    g00 = arb(0)
    g0b = arb(0)
    gbb = arb(0)
    for q, p, a in prime_powers_below(L, arb):
        w = arb(p).log() / arb(q).sqrt()
        # Use float mid kernels then lift — kernels are rational in (L,a).
        # Evaluate with Arb arithmetic:
        d = L - a
        k00 = 2 * d
        k0b = (d**2) * (L + 2 * a) / 3
        kbb = (d**3) * (L**2 + 3 * L * a + a**2) / 15
        g00 += w * k00
        g0b += w * k0b
        gbb += w * kbb
    return g00, g0b, gbb


def gp_odd_pivot(L, arb):
    """Prime block for O1 = G[q1,q1]."""
    o1 = arb(0)
    for q, p, a in prime_powers_below(L, arb):
        w = arb(p).log() / arb(q).sqrt()
        d = L - a
        k = d * (L**2 - 2 * L * a - 2 * a**2) / 6
        o1 += w * k
    return o1


def pole_odd_A(L, arb):
    """A(L)=L cosh(L/4)-4 sinh(L/4); G0[q1,q1]=-8 A^2."""
    return L * (L / 4).cosh() - 4 * (L / 4).sinh()


def g0_odd_pivot(L, arb):
    A = pole_odd_A(L, arb)
    return -8 * (A**2)


def pole_even_helpers(L, arb):
    """Eb± helpers from ENG-002 (rank-1 pole even block)."""
    eL2 = (L / 2).exp()
    emL2 = (-L / 2).exp()
    ebp = 4 * ((L - 4) * eL2 + L + 4)
    ebm = 4 * ((L - 4) + (L + 4) * emL2)
    return ebp, ebm


def ginf_even_block_quad(L, T: int, arb, n: int = 512):
    """Trapezoid Arb quadrature for even Ginf on [0,T] (working precision via ctx).

    This is a rigorous *enclosure* only when combined with an explicit
    quadrature-error / tail bound. Until the error binder is attached, callers
    must label results E3 or attach a bound before E1 promotion.
    """
    require_flint()
    L_a = arb(L)
    T_a = arb(T)
    g00 = arb(0)
    g0b = arb(0)
    gbb = arb(0)
    for i in range(n + 1):
        t = T_a * i / n
        w = arb("0.5") if i in (0, n) else arb(1)
        hp = h_plus(t)
        h0sq, re0b, hbsq = stable_products_even(t, L_a, arb)
        g00 += w * hp * h0sq
        g0b += w * hp * re0b
        gbb += w * hp * hbsq
    dt = T_a / n
    pi = arb.pi()
    return g00 * dt / pi, g0b * dt / pi, gbb * dt / pi


def finite_weil_even_block(L, T: int = 84, precision_bits: int = 256, n_quad: int = 512) -> Dict[str, Any]:
    """Structured even-block assembly at cutoff T."""
    _, arb, _, _ = require_flint()
    set_precision_bits(precision_bits)
    L_a = arb(L)
    gp00, gp0b, gpbb = gp_even_block(L_a, arb)
    # Compact pole even block: rank ≤ 1, det 0. Use Eb helpers as G0 vector factors.
    # ENG-002: even pole block has determinant zero. For the Gram we need the
    # actual G0[{1,b}] entries — until the explicit outer-product form is wired,
    # set G0_even = 0 and rely on -Gp + Ginf for the finite-T even block used in
    # Fourier certificates (compact cutoff-free degree-2 uses a separate path).
    g000 = arb(0)
    g00b = arb(0)
    g0bb = arb(0)
    gi00, gi0b, gibb = ginf_even_block_quad(L_a, T, arb, n=n_quad)
    G00 = g000 - gp00 + gi00
    G0b = g00b - gp0b + gi0b
    Gbb = g0bb - gpbb + gibb
    E2 = G00 * Gbb - G0b * G0b
    return {
        "G00": G00,
        "G0b": G0b,
        "Gbb": Gbb,
        "E2": E2,
        "normalization": NORMALIZATION,
        "cutoff_T": T,
        "evidence_class": "E3_PENDING_QUADRATURE_BOUND",
        "rh_proof_claim": False,
        "note": (
            "True Weil integrand (h_+ · stable products) minus prime block. "
            "Pole even outer-product and quadrature-error binder still required for E1."
        ),
    }


def finite_weil_odd_pivot(L, T: int = 84, precision_bits: int = 256, n_quad: int = 512) -> Dict[str, Any]:
    """O1 at finite T: G0 - Gp + Ginf for q1 (Ginf via Hb relation)."""
    _, arb, _, _ = require_flint()
    set_precision_bits(precision_bits)
    L_a = arb(L)
    g0 = g0_odd_pivot(L_a, arb)
    gp = gp_odd_pivot(L_a, arb)
    # |Hq1|^2 = (t^2/4) |Hb|^2 from Hb=-2i/t Hq1 ⇒ |Hb|=2/|t| |Hq1|
    T_a = arb(T)
    o_inf = arb(0)
    n = n_quad
    for i in range(n + 1):
        t = T_a * i / n
        w = arb("0.5") if i in (0, n) else arb(1)
        hp = h_plus(t)
        if t.contains(0) or abs(t.mid()) <= t.rad():
            # ∫ x(L-x) related; Hq1 at 0: ∫_0^L (x-L/2) dx = 0
            hq1sq = arb(0)
        else:
            _, _, hbsq = stable_products_even(t, L_a, arb)
            hq1sq = hbsq * (t**2) / 4
        o_inf += w * hp * hq1sq
    o_inf = o_inf * (T_a / n) / arb.pi()
    O1 = g0 - gp + o_inf
    return {
        "O1": O1,
        "normalization": NORMALIZATION,
        "cutoff_T": T,
        "evidence_class": "E3_PENDING_QUADRATURE_BOUND",
        "rh_proof_claim": False,
    }


def finite_weil_degree2(L, T: int = 84, precision_bits: int = 256, n_quad: int = 512) -> Dict[str, Any]:
    even = finite_weil_even_block(L, T=T, precision_bits=precision_bits, n_quad=n_quad)
    odd = finite_weil_odd_pivot(L, T=T, precision_bits=precision_bits, n_quad=n_quad)
    _, arb, _, _ = require_flint()
    L_a = arb(L)
    D2 = even["E2"] + (L_a**2) * even["G00"] * odd["O1"]
    full = odd["O1"] * even["E2"]
    return {
        **even,
        "O1": odd["O1"],
        "D2": D2,
        "full_det": full,
        "cutoff_T": T,
        "normalization": NORMALIZATION,
        "rh_proof_claim": False,
    }


def finite_weil_entry(i: str, j: str, L, T: int = 84, backend: str = "flint"):
    if backend != "flint":
        raise FlintUnavailable("finite_weil_entry E-path requires flint backend")
    block = finite_weil_degree2(L, T=T)
    key = {
        ("1", "1"): "G00",
        ("1", "b"): "G0b",
        ("b", "1"): "G0b",
        ("b", "b"): "Gbb",
        ("q1", "q1"): "O1",
    }[(i, j)]
    return block[key]
