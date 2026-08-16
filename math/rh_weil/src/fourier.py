"""Stable Fourier forms and Taylor jets in support length L (WO-RH-05).

Requires mpmath. Emits floating/E3 scans unless interval coverage is closed.
No RH claim; finite cutoff T only.
"""
from __future__ import annotations

from typing import Any, Callable


def require_mpmath():
    try:
        import mpmath as mp
    except ImportError as exc:  # pragma: no cover
        raise ImportError("mpmath required for Fourier forms") from exc
    return mp


def H0(t: Any, L: Any):
    """Stable H0(t;L) = ∫_0^L e^{i t x} dx = L e^{i z} sinc(z), z=Lt/2."""
    mp = require_mpmath()
    t_m = mp.mpf(t)
    L_m = mp.mpf(L)
    if t_m == 0:
        return mp.mpc(L_m, 0)
    z = L_m * t_m / 2
    if abs(z) < mp.mpf("1e-18"):
        # Taylor: sinc(z) = 1 - z^2/6 + ...
        a = 1 - z * z / 6
    else:
        a = mp.sin(z) / z
    return L_m * mp.exp(1j * z) * a


def Hb(t: Any, L: Any):
    """Stable Hb for b=x(L-x): (L^3/2) e^{i z} B(z), B=(sin z - z cos z)/z^3."""
    mp = require_mpmath()
    t_m = mp.mpf(t)
    L_m = mp.mpf(L)
    if t_m == 0:
        return mp.mpc(L_m**3 / 6, 0)  # ∫_0^L x(L-x) dx = L^3/6
    z = L_m * t_m / 2
    if abs(z) < mp.mpf("1e-8"):
        # B(z) = 1/3 - z^2/30 + O(z^4)
        B = mp.mpf("1") / 3 - z * z / 30
    else:
        B = (mp.sin(z) - z * mp.cos(z)) / (z**3)
    return (L_m**3 / 2) * mp.exp(1j * z) * B


def H0_L_jet(t: Any, L: Any, order: int = 2):
    """Taylor jet of H0 in L at fixed t, by direct differentiation (not FD)."""
    mp = require_mpmath()
    t_m = mp.mpf(t)
    L_m = mp.mpf(L)
    # H0 = (e^{i t L} - 1)/(i t) for t≠0; ∂_L H0 = e^{i t L}
    # ∂_L² H0 = i t e^{i t L}
    if t_m == 0:
        # H0=L, H0'=1, H0''=0
        vals = [mp.mpc(L_m, 0), mp.mpc(1, 0)]
        for _ in range(2, order + 1):
            vals.append(mp.mpc(0, 0))
        return vals[: order + 1]
    e = mp.exp(1j * t_m * L_m)
    h = (e - 1) / (1j * t_m)
    d1 = e
    vals = [h, d1]
    # Higher: ∂_L^k H0 = (i t)^{k-1} e^{i t L} for k>=1
    for k in range(2, order + 1):
        vals.append(((1j * t_m) ** (k - 1)) * e)
    return vals


def Hb_L_jet(t: Any, L: Any, order: int = 1):
    """Low-order jet of Hb via stable form + mpmath diff (order≤1 analytic)."""
    mp = require_mpmath()
    # ∂_L Hb is available numerically from the closed form; for order 0 return Hb.
    vals = [Hb(t, L)]
    if order >= 1:
        # Differentiate: b=x(L-x) ⇒ ∂_L of integral uses Leibniz + parameter.
        # Use complex step / high-dps difference only as last resort — prefer
        # analytic: Hb = ∫_0^L x(L-x) e^{itx} dx.
        # ∂_L Hb = ∫_0^L x e^{itx} dx  (the upper-limit term vanishes: L(L-L)=0)
        t_m = mp.mpf(t)
        L_m = mp.mpf(L)
        if t_m == 0:
            d1 = mp.mpc(L_m**2 / 2, 0)
        else:
            # ∫_0^L x e^{itx} dx = e^{itL}(itL-1)/(it)^2 + 1/(it)^2
            it = 1j * t_m
            d1 = mp.exp(it * L_m) * (it * L_m - 1) / (it**2) + 1 / (it**2)
        vals.append(d1)
    return vals[: order + 1]


def fourier_even_gram_entries(L: Any, T: int, dps: int = 40):
    """Assemble floating G00, G0b, Gbb archimedean-style probes at cutoff T.

    This is a **diagnostic E3 scaffold** using |H|^2-type energy, not the full
    Weil G^0-G^p+G^∞ matrix. Full interval E2,84 coverage remains WO-RH-05 open
    work until the complete finite-Fourier Gram is regenerated.
    """
    mp = require_mpmath()
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        L_m = mp.mpf(L)
        # Trapezoid on a coarse t-grid in [0,T] — heuristic only.
        n = max(32, int(T))
        g00 = mp.mpf(0)
        g0b = mp.mpf(0)
        gbb = mp.mpf(0)
        for i in range(n + 1):
            t = mp.mpf(T) * i / n
            w = mp.mpf("0.5") if i in (0, n) else mp.mpf(1)
            h0 = H0(t, L_m)
            hb = Hb(t, L_m)
            # Use real parts of conjugated products (energy).
            g00 += w * (mp.re(h0) ** 2 + mp.im(h0) ** 2)
            g0b += w * (mp.re(h0) * mp.re(hb) + mp.im(h0) * mp.im(hb))
            gbb += w * (mp.re(hb) ** 2 + mp.im(hb) ** 2)
        dt = mp.mpf(T) / n
        g00 *= dt
        g0b *= dt
        gbb *= dt
        det = g00 * gbb - g0b * g0b
        return {
            "G00": g00,
            "G0b": g0b,
            "Gbb": gbb,
            "E2_probe": det,
            "evidence_class": "E3",
            "note": "Heuristic |H|^2 energy probe — not the certified Weil Gram.",
        }
    finally:
        mp.mp.dps = old


def scan_E2_probe(L_values, T: int = 84, dps: int = 25):
    """E3 scan of the even-block probe along L samples."""
    mp = require_mpmath()
    rows = []
    for L in L_values:
        ent = fourier_even_gram_entries(L, T, dps=dps)
        rows.append(
            {
                "L": str(L),
                "E2_probe": mp.nstr(ent["E2_probe"], 20),
                "sign": int(mp.sign(ent["E2_probe"])),
            }
        )
    return rows
