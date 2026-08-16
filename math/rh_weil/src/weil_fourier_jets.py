"""Analytic L-jets of the finite-T even Weil Gram (WO-RH-15 support).

Differentiates G0, Gp, and the archimedean integrand under the integral sign.
Point Arb evaluations keep radii tiny; trapezoid remainders use M2 majorants
on the differentiated integrands.
"""
from __future__ import annotations

from typing import Any, Dict, Tuple

from archimedean import h_plus
from finite_weil import (
    g0_even_block,
    gp_even_block,
    pole_even_I0,
    pole_even_helpers,
)
from interval_backend import require_flint, set_precision_bits


def _dL_g0_even(L, arb) -> Tuple[Any, Any, Any]:
    """∂_L of even pole Gram (√3/2)(v₊v₊ᵀ+v₋v₋ᵀ)."""
    scale = arb(3).sqrt() / 2
    i0p, i0m = pole_even_I0(L, arb)
    ibp, ibm = pole_even_helpers(L, arb)
    # ∂_L I0+ = e^{L/2}, ∂_L I0- = e^{-L/2}
    d0p = (L / 2).exp()
    d0m = (-L / 2).exp()
    # Eb+ = 4[(L-4)e^{L/2}+L+4]
    # ∂_L Eb+ = 4[e^{L/2} + (L-4)e^{L/2}/2 + 1] = 4[e^{L/2}(1+(L-4)/2)+1]
    dbp = 4 * ((L / 2).exp() * (1 + (L - 4) / 2) + 1)
    # Eb- = 4[(L-4)+(L+4)e^{-L/2}]
    # ∂_L Eb- = 4[1 + e^{-L/2} - (L+4)e^{-L/2}/2] = 4[1 + e^{-L/2}(1-(L+4)/2)]
    dbm = 4 * (1 + (-L / 2).exp() * (1 - (L + 4) / 2))
    # ∂(v vᵀ) = v' vᵀ + v (v')ᵀ
    g00 = scale * (2 * i0p * d0p + 2 * i0m * d0m)
    g0b = scale * (d0p * ibp + i0p * dbp + d0m * ibm + i0m * dbm)
    gbb = scale * (2 * ibp * dbp + 2 * ibm * dbm)
    return g00, g0b, gbb


def _dL_gp_even(L, arb) -> Tuple[Any, Any, Any]:
    """∂_L of prime even block (kernels depend on L; active set via mid compare)."""
    from finite_weil import prime_powers_below

    g00 = arb(0)
    g0b = arb(0)
    gbb = arb(0)
    for q, p, a in prime_powers_below(L, arb):
        w = arb(p).log() / arb(q).sqrt()
        d = L - a
        # K00=2d ⇒ ∂_L=2
        # K0b=d²(L+2a)/3 ⇒ ∂_L = [2d(L+2a) + d²]/3 = d(2L+4a+d)/3 = d(2L+4a+L-a)/3 = d(3L+3a)/3 = d(L+a)
        # Kbb=d³(L²+3La+a²)/15
        # ∂_L Kbb = 3d²(L²+3La+a²)/15 + d³(2L+3a)/15
        dk00 = arb(2)
        dk0b = d * (L + a)
        dkbb = (
            3 * (d**2) * (L**2 + 3 * L * a + a**2) / 15
            + (d**3) * (2 * L + 3 * a) / 15
        )
        g00 += w * dk00
        g0b += w * dk0b
        gbb += w * dkbb
    return g00, g0b, gbb


def _dL_products_even(t, L, arb):
    """∂_L (|H0|², Re conj(H0)Hb, |Hb|²) at fixed t."""
    t_a = arb(t)
    L_a = arb(L)
    if abs(float(t_a.mid())) < 1e-14:
        # H0=L, Hb=L³/6; |H0|²=L² ⇒ ∂=2L
        # Re conj(H0)Hb = L⁴/6 ⇒ ∂=4 L³/6 = 2 L³/3
        # |Hb|²=L^6/36 ⇒ ∂=6 L^5/36 = L^5/6
        return 2 * L_a, (2 * L_a**3) / 3, (L_a**5) / 6
    # H0 = L e^{iz} A, z=Lt/2; ∂_L H0 = e^{i t L}
    # Use product rule on stable real forms:
    z = L_a * t_a / 2
    s = z.sin()
    c = z.cos()
    A = s / z
    B = (s - z * c) / (z**3)
    Ap = (z * c - s) / (z**2)
    Bp = (z * z * s - 3 * (s - z * c)) / (z**4)
    zt_L = t_a / 2  # ∂z/∂L
    # p00 = L² A²; ∂_L = 2L A² + L² * 2 A Ap * zt_L
    p00_L = 2 * L_a * (A**2) + 2 * (L_a**2) * A * Ap * zt_L
    # p0b = (L⁴/2) A B; ∂_L = (4 L³/2)AB + (L⁴/2)(Ap B + A Bp) zt_L
    p0b_L = 2 * (L_a**3) * A * B + (L_a**4 / 2) * (Ap * B + A * Bp) * zt_L
    # pbb = (L⁶/4) B²; ∂_L = (6 L⁵/4) B² + (L⁶/4)*2 B Bp * zt_L
    pbb_L = (3 * L_a**5 / 2) * (B**2) + (L_a**6 / 2) * B * Bp * zt_L
    return p00_L, p0b_L, pbb_L


def ginf_even_dL_quad(L, T: int, arb, n: int = 8192):
    """Trapezoid of (1/π)∫ h_+ ∂_L(products) on [0,T] (point nodes)."""
    L_a = arb(L)
    T_a = arb(T)
    g00 = arb(0)
    g0b = arb(0)
    gbb = arb(0)
    for i in range(n + 1):
        t = T_a * i / n
        w = arb("0.5") if i in (0, n) else arb(1)
        hp = h_plus(t)
        d00, d0b, dbb = _dL_products_even(t, L_a, arb)
        g00 += w * hp * d00
        g0b += w * hp * d0b
        gbb += w * hp * dbb
    dt = T_a / n
    pi = arb.pi()
    return g00 * dt / pi, g0b * dt / pi, gbb * dt / pi


def even_E2_and_derivative(
    L,
    T: int = 84,
    precision_bits: int = 192,
    n_quad: int = 16384,
) -> Dict[str, Any]:
    """Return E2 and a first L-derivative enclosure (quad remainder not yet tight)."""
    _, arb, _, _ = require_flint()
    set_precision_bits(precision_bits)
    L_a = arb(L)
    g0 = g0_even_block(L_a, arb)
    gp = gp_even_block(L_a, arb)
    from finite_weil import ginf_even_block_quad

    gi = ginf_even_block_quad(L_a, T, arb, n=n_quad)
    G00 = g0[0] - gp[0] + gi[0]
    G0b = g0[1] - gp[1] + gi[1]
    Gbb = g0[2] - gp[2] + gi[2]
    E2 = G00 * Gbb - G0b * G0b

    dg0 = _dL_g0_even(L_a, arb)
    dgp = _dL_gp_even(L_a, arb)
    dgi = ginf_even_dL_quad(L_a, T, arb, n=n_quad)
    dG00 = dg0[0] - dgp[0] + dgi[0]
    dG0b = dg0[1] - dgp[1] + dgi[1]
    dGbb = dg0[2] - dgp[2] + dgi[2]
    # E' = G00' Gbb + G00 Gbb' - 2 G0b G0b'
    E2p = dG00 * Gbb + G00 * dGbb - 2 * G0b * dG0b
    return {
        "G00": G00,
        "G0b": G0b,
        "Gbb": Gbb,
        "E2": E2,
        "E2_first": E2p,
        "dG00": dG00,
        "dG0b": dG0b,
        "dGbb": dGbb,
        "rh_proof_claim": False,
    }
