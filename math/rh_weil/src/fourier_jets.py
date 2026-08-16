"""Analytic support-length Taylor jets for H0/Hb (WO-RH-14).

No finite differences on the E1 path. Numerical FD is tests-only.
"""
from __future__ import annotations

from typing import Any, List


def require_mpmath():
    import mpmath as mp

    return mp


def H0_L_jets(t: Any, L: Any, order: int) -> List[Any]:
    """∂_L^n H0: n=0 → H0; n≥1 → (i t)^{n-1} e^{i t L}."""
    mp = require_mpmath()
    t_m = mp.mpf(t)
    L_m = mp.mpf(L)
    if t_m == 0:
        vals = [mp.mpc(L_m, 0)]
        if order >= 1:
            vals.append(mp.mpc(1, 0))
        for _ in range(2, order + 1):
            vals.append(mp.mpc(0, 0))
        return vals[: order + 1]
    e = mp.exp(1j * t_m * L_m)
    h0 = (e - 1) / (1j * t_m)
    vals = [h0]
    for n in range(1, order + 1):
        vals.append(((1j * t_m) ** (n - 1)) * e)
    return vals


def Hb_L_jets(t: Any, L: Any, order: int) -> List[Any]:
    """Analytic jets of Hb=∫_0^L x(L-x) e^{i t x} dx.

    ∂_L Hb = ∫_0^L x e^{itx} dx
    ∂_L² Hb = L e^{itL}
    ∂_L^n Hb = e^{itL} [ L (it)^{n-2} + (n-2)(it)^{n-3} ] for n≥3
    """
    mp = require_mpmath()
    t_m = mp.mpf(t)
    L_m = mp.mpf(L)
    # n=0
    if t_m == 0:
        hb0 = mp.mpc(L_m**3 / 6, 0)
    else:
        z = L_m * t_m / 2
        if abs(z) < mp.mpf("1e-8"):
            B = mp.mpf(1) / 3 - z * z / 30
        else:
            B = (mp.sin(z) - z * mp.cos(z)) / (z**3)
        hb0 = (L_m**3 / 2) * mp.exp(1j * z) * B
    vals = [hb0]
    if order < 1:
        return vals
    # n=1
    if t_m == 0:
        d1 = mp.mpc(L_m**2 / 2, 0)
    else:
        it = 1j * t_m
        d1 = mp.exp(it * L_m) * (it * L_m - 1) / (it**2) + 1 / (it**2)
    vals.append(d1)
    if order < 2:
        return vals
    # n=2
    if t_m == 0:
        d2 = mp.mpc(L_m, 0)
    else:
        d2 = L_m * mp.exp(1j * t_m * L_m)
    vals.append(d2)
    # n≥3
    for n in range(3, order + 1):
        if t_m == 0:
            vals.append(mp.mpc(0, 0))
        else:
            it = 1j * t_m
            e = mp.exp(it * L_m)
            vals.append(e * (L_m * (it ** (n - 2)) + (n - 2) * (it ** (n - 3))))
    return vals


def re_conj_product_jet(U_jets: List[Any], V_jets: List[Any], n: int):
    """∂_L^n Re(conj(U) V) via Leibniz."""
    mp = require_mpmath()
    total = mp.mpc(0)
    for k in range(n + 1):
        total += mp.binomial(n, k) * mp.conj(U_jets[k]) * V_jets[n - k]
    return mp.re(total)
