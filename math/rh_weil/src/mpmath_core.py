"""Atlas-native mpmath formulas for RH/Weil cross-checks.

Optional dependency: ``mpmath``. Exact stdlib tests must not import this module.
These implementations are independent of ``connes-cvs`` (and of its flint path).
"""
from __future__ import annotations

from typing import Any


def require_mpmath():
    try:
        import mpmath as mp
    except ImportError as exc:  # pragma: no cover
        raise ImportError(
            "mpmath is required for RH/Weil numeric cross-checks "
            "(optional research dependency)."
        ) from exc
    return mp


def atlas_h_plus(tau: Any, dps: int = 80):
    """Archimedean multiplier Re ψ(1/4 + i τ/2) − log(π), Atlas-native mpmath."""
    mp = require_mpmath()
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        tau_m = mp.mpf(tau)
        z = mp.mpc(mp.mpf("0.25"), tau_m / 2)
        return mp.re(mp.digamma(z)) - mp.log(mp.pi)
    finally:
        mp.mp.dps = old


def atlas_h0(t: Any, L: Any):
    """Stable H0(t;L) = ∫_0^L exp(i t x) dx (FORMULAS.md)."""
    mp = require_mpmath()
    t_m = mp.mpf(t)
    L_m = mp.mpf(L)
    if t_m == 0:
        return mp.mpc(L_m, 0)
    # Entire low-frequency form: L exp(i z) sinc(z), z = L t / 2
    z = L_m * t_m / 2
    if z == 0:
        return mp.mpc(L_m, 0)
    a = mp.sin(z) / z
    return L_m * mp.exp(1j * z) * a


def atlas_prime_powers_up_to(c: int, dps: int = 80):
    """Prime powers n≤c with exact n/base and high-precision (log n, Λ(n)/√n)."""
    if c < 2:
        raise ValueError("c must be >= 2")
    mp = require_mpmath()
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        is_prime = [True] * (c + 1)
        is_prime[0] = is_prime[1] = False
        for i in range(2, int(c**0.5) + 1):
            if is_prime[i]:
                step = i
                start = i * i
                for j in range(start, c + 1, step):
                    is_prime[j] = False
        primes = [i for i in range(2, c + 1) if is_prime[i]]
        base_by_n: dict[int, int] = {}
        for p in primes:
            pk = p
            while pk <= c:
                base_by_n[pk] = p
                if pk > c // p:
                    break
                pk *= p
        data = [
            (n, mp.log(n), mp.log(base_by_n[n]) / mp.sqrt(n))
            for n in sorted(base_by_n)
        ]
        return data, primes
    finally:
        mp.mp.dps = old


def atlas_scalar_arch_probe(L: Any, T: int, dps: int, h_plus_fn) -> Any:
    """Diagnostic scalar archimedean probe isolating ``h_plus``.

    I(L,T) = (1/π) ∫_0^T h_plus(t) |H0(t;L)|² dt

    This is an Atlas-side assembly convention for cross-validation only.
    It is not a positivity certificate and does not claim RH.
    """
    mp = require_mpmath()
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        L_m = mp.mpf(L)
        T_m = mp.mpf(T)

        def integrand(t):
            h0 = atlas_h0(t, L_m)
            return h_plus_fn(t, dps) * (mp.re(h0) ** 2 + mp.im(h0) ** 2)

        # Split at 0 already; mild singularity handling via open interval.
        total = mp.quad(integrand, [mp.mpf("1e-30"), T_m])
        return total / mp.pi
    finally:
        mp.mp.dps = old
