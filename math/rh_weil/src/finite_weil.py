"""True finite-cutoff Weil Gram (WO-RH-12). Not the E3 energy probe.

Assembly: G = G0 - Gp + Ginf_T
Ginf_T[i,j] = (1/π) ∫_0^T h_+(t) Re(conj(H_i(t;L)) H_j(t;L)) dt

Even pole block (ENG-004 §1): the ADOPTED Candidate A, delegated to
``pole.pole_gram_matrix``:

  G0_ij = E_i^+ E_j^- + E_i^- E_j^+,  E_i^± = int_0^L h_i(x) e^{±x/2} dx.

The former ``(sqrt(3)/2)(v+v+^T + v-v-^T)`` calibration was REJECTED by
WO-RH-17 and has been removed from this module; it survives only in the
archival ``rejected_pole`` module, which production may not import.
"""
from __future__ import annotations

from typing import Any, Dict, Tuple

import core
import pole
from archimedean import h_plus, h_plus_derivatives
from interval_backend import FlintUnavailable, require_flint, set_precision_bits


NORMALIZATION = core.NORMALIZATION
# ENG-004 §1: production carries the adopted candidate only. The rejected
# ``sqrt(3)/2`` scale is not defined here -- it lives in ``rejected_pole``.
POLE_CANDIDATE = pole.POLE_CANDIDATE          # "A"
POLE_STATUS = pole.POLE_STATUS                # "ADOPTED_WO_RH_17"
POLE_FORMULA = pole.POLE_FORMULA
POLE_EVEN_SCALE_STATUS = "REJECTED_WO_RH_17"  # verdict on the removed constant
POLE_EVEN_SCALE_SUPERSEDED_BY = "pole.pole_gram_entry (Candidate A)"


def _A_B(z, arb):
    """sinc and bubble shape factors; z-ball must not contain 0 for power forms."""
    if z.contains(0) or (hasattr(z, "rad") and z.rad() >= abs(z.mid())):
        raise ValueError("z-ball contains 0; use t=0 analytic limits")
    A = z.sin() / z
    B = (z.sin() - z * z.cos()) / (z**3)
    return A, B


def stable_products_even(t, L, arb):
    """Return (|H0|^2, Re(conj(H0)Hb), |Hb|^2) as Arb balls."""
    t_a = arb(t)
    L_a = arb(L)
    if t_a.contains(0) or abs(t_a.mid()) <= t_a.rad():
        h0sq = L_a**2
        re0b = L_a**4 / 6
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
    """Eb± helpers from ENG-002 (= I_b± exactly)."""
    eL2 = (L / 2).exp()
    emL2 = (-L / 2).exp()
    ebp = 4 * ((L - 4) * eL2 + L + 4)
    ebm = 4 * ((L - 4) + (L + 4) * emL2)
    return ebp, ebm


def pole_even_I0(L, arb):
    """I0± = ∫_0^L e^{±x/2} dx."""
    i0p = 2 * ((L / 2).exp() - 1)
    i0m = 2 * (1 - (-L / 2).exp())
    return i0p, i0m


def g0_even_block(L, arb):
    """Even pole Gram under Candidate A: ``G0_ij = E_i^+E_j^- + E_i^-E_j^+``.

    Delegates to ``pole.pole_gram_matrix`` so there is exactly one pole formula
    in the tree (ENG-004 §1). ``pole_even_I0``/``pole_even_helpers`` remain as
    the named ``E^±`` closed forms and are pinned against it by
    ``tests/test_pole_primitive.py``.
    """
    g = pole.pole_gram_matrix(("one", "b"), L)
    return g[0][0], g[0][1], g[1][1]


def _product_t_derivatives(t, L, arb, acb):
    """Return (p00,p0b,pbb) and first/second t-derivatives at a point-like t."""
    t_a = arb(t)
    L_a = arb(L)
    if abs(float(t_a.mid())) < 1e-14:
        p00 = L_a**2
        p0b = L_a**4 / 6
        pbb = (L_a**3 / 6) ** 2
        # At t=0: p00' = 0, p00'' = -L^4/6
        p00_t = arb(0)
        p0b_t = arb(0)
        pbb_t = arb(0)
        p00_tt = -(L_a**4) / 6
        p0b_tt = -(L_a**6) / 60
        pbb_tt = -(L_a**8) / 3780
        return (
            (p00, p0b, pbb),
            (p00_t, p0b_t, pbb_t),
            (p00_tt, p0b_tt, pbb_tt),
        )
    z = L_a * t_a / 2
    s = z.sin()
    c = z.cos()
    A = s / z
    B = (s - z * c) / (z**3)
    Ap = (z * c - s) / (z**2)
    # dB/dz: B = (sin - z cos)/z^3
    Bp = (z * s - 3 * (s - z * c) / z) / (z**3)
    # cleaner: Bp = (z^2 sin - 3 z (sin - z cos) wait)
    # d/dz[(sin-z cos)/z^3] = [(cos - cos + z sin)*z^3 - (sin-z cos)*3z^2]/z^6
    # = [z sin * z^3 - 3z^2 (sin - z cos)] / z^6 = [z^2 sin - 3(sin - z cos)] / z^4
    Bp = (z * z * s - 3 * (s - z * c)) / (z**4)
    App = (-z * s - 2 * c + 2 * s / z) / (z * z)
    # d²B/dz² via quotient; use numerical-free formula:
    # B' = (z² sin - 3 sin + 3 z cos) / z^4
    # Differentiate: num = z² s - 3 s + 3 z c, den = z^4
    num = z * z * s - 3 * s + 3 * z * c
    num_p = 2 * z * s + z * z * c - 3 * c + 3 * c - 3 * z * s
    # = 2z s + z² c - 3z s = z² c - z s
    num_p = z * z * c - z * s
    Bpp = (num_p * (z**4) - num * 4 * (z**3)) / (z**8)

    zt = L_a / 2
    p00 = (L_a**2) * (A**2)
    p0b = (L_a**4 / 2) * A * B
    pbb = (L_a**6 / 4) * (B**2)

    p00_z = 2 * (L_a**2) * A * Ap
    p0b_z = (L_a**4 / 2) * (Ap * B + A * Bp)
    pbb_z = (L_a**6 / 2) * B * Bp

    p00_zz = 2 * (L_a**2) * (Ap**2 + A * App)
    p0b_zz = (L_a**4 / 2) * (App * B + 2 * Ap * Bp + A * Bpp)
    pbb_zz = (L_a**6 / 2) * (Bp**2 + B * Bpp)

    return (
        (p00, p0b, pbb),
        (p00_z * zt, p0b_z * zt, pbb_z * zt),
        (p00_zz * zt * zt, p0b_zz * zt * zt, pbb_zz * zt * zt),
    )


def _integrand_second_derivatives(t, L, arb, acb):
    """(f00'', f0b'', fbb'') for f = h_+ * product / π."""
    h0, hp, hpp = h_plus_derivatives(t)
    (p00, p0b, pbb), (p00t, p0bt, pbbt), (p00tt, p0btt, pbbtt) = _product_t_derivatives(
        t, L, arb, acb
    )
    pi = arb.pi()

    def ftt(h, hp_, hpp_, p, pt, ptt):
        return (hpp_ * p + 2 * hp_ * pt + h * ptt) / pi

    return (
        ftt(h0, hp, hpp, p00, p00t, p00tt),
        ftt(h0, hp, hpp, p0b, p0bt, p0btt),
        ftt(h0, hp, hpp, pbb, pbbt, pbbtt),
    )


# Conservative majorants for |f'''| on t∈[0,84], L∈[log3,log4] (empirical ≪ these).
_M3_MAJORANT = (200.0, 50.0, 20.0)


def enclose_integrand_M2(L, T: float, n_sample: int = 4000) -> Tuple[Any, Any, Any]:
    """Rigorous upper bounds on max|f''| for the three even integrands on [0,T].

    Samples analytic f'' at a uniform grid and inflates by (h/2)·M3 majorant so the
    continuous max is enclosed.
    """
    _, arb, acb, _ = require_flint()
    L_a = arb(L)
    h = float(T) / n_sample
    m00 = 0.0
    m0b = 0.0
    mbb = 0.0
    for i in range(n_sample + 1):
        t = i * h
        f00tt, f0btt, fbbtt = _integrand_second_derivatives(t, L_a, arb, acb)
        m00 = max(m00, float(abs(f00tt).upper()))
        m0b = max(m0b, float(abs(f0btt).upper()))
        mbb = max(mbb, float(abs(fbbtt).upper()))
    return (
        arb(m00 + 0.5 * h * _M3_MAJORANT[0]),
        arb(m0b + 0.5 * h * _M3_MAJORANT[1]),
        arb(mbb + 0.5 * h * _M3_MAJORANT[2]),
    )


def ginf_even_block_quad(
    L,
    T: int,
    arb,
    n: int = 512,
    *,
    with_error_bound: bool = False,
    m2_bounds=None,
    n_m2_sample: int = 4000,
):
    """Trapezoid Arb quadrature for even Ginf on [0,T].

    When ``with_error_bound`` is True, adds a rigorous trapezoid remainder ball
    using composite error |(b-a) h²/12| · max|f''| with enclosed M2.
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
    g00 = g00 * dt / pi
    g0b = g0b * dt / pi
    gbb = gbb * dt / pi
    if with_error_bound:
        if m2_bounds is None:
            m2_bounds = enclose_integrand_M2(L_a, float(T), n_sample=n_m2_sample)
        h = dt
        # |E| ≤ (b-a) h²/12 M2 = T * (T/n)² / 12 * M2
        factor = T_a * (h**2) / 12
        g00 = g00 + arb(0, float((factor * m2_bounds[0]).upper()))
        g0b = g0b + arb(0, float((factor * m2_bounds[1]).upper()))
        gbb = gbb + arb(0, float((factor * m2_bounds[2]).upper()))
    return g00, g0b, gbb


def archimedean_tail_even(L, T: float, arb) -> Tuple[Any, Any, Any]:
    """Rigorous absolute-value majorant for ∫_T^∞ (h_+/π) · products dt.

    For t ≥ T ≥ 2/L: |A(z)| ≤ 1/|z|, |B(z)| ≤ 2/|z|² roughly.
    Uses |H0|² ≤ 4/t², |Re H0̄Hb| ≤ 2 L / t³ · L³? — see bounds below,
    and h_+(t) ≤ log(t) + 1 for t ≥ 2 (validated majorant on the critical strip).
    """
    L_a = arb(L)
    T_a = arb(T)
    if float(T_a.mid()) < 2.0:
        raise ValueError("tail bound requires T ≥ 2")
    # h_+(t) ≤ log(t) + 1 for t ≥ 2 (Re ψ(1/4+it/2) ≤ log|t/2| + C with C≤1+log2)
    # ∫_T^∞ (log t + 1) / t² dt = (log T + 1)/T + 1/T = (log T + 2)/T
    # |H0|² ≤ 4/t²  ⇒  ∫ h |H0|²/π ≤ (4/π) ∫ (log t+1)/t² = (4/π)(log T+2)/T
    logT = T_a.log()
    I_log_t2 = (logT + 2) / T_a
    tail00 = arb(4) / arb.pi() * I_log_t2
    # |Re conj(H0)Hb| ≤ (L^4/2) |A| |B| ≤ (L^4/2)(1/|z|)(2/|z|²) = L^4 / |z|³
    # z=Lt/2 ⇒ |z|=Lt/2 ⇒ L^4 / (L t/2)³ = L^4 * 8 / (L³ t³) = 8 L / t³
    # ∫_T^∞ (log t+1)/t³ dt ≤ ∫_T^∞ (log t+1)/t³ dt = (log T+1)/(2 T²) + 1/(4 T²)
    I_log_t3 = (logT + 1) / (2 * T_a**2) + 1 / (4 * T_a**2)
    tail0b = arb(8) * L_a / arb.pi() * I_log_t3
    # |Hb|² ≤ (L^6/4) (2/|z|²)² = (L^6/4)*4/|z|^4 = L^6 / (L t/2)^4 = L^6 * 16 / (L^4 t^4) = 16 L² / t^4
    # ∫ (log t+1)/t^4 dt from T: (log T+1)/(3 T³) + 1/(9 T³)
    I_log_t4 = (logT + 1) / (3 * T_a**3) + 1 / (9 * T_a**3)
    tailbb = arb(16) * (L_a**2) / arb.pi() * I_log_t4
    return tail00, tail0b, tailbb


def _require_backend(backend: str) -> None:
    if backend != "flint":
        raise FlintUnavailable(
            f"backend={backend!r} unsupported; E-path requires backend='flint'"
        )


def finite_weil_even_block(
    L,
    T: int = 84,
    precision_bits: int = 256,
    n_quad: int = 65536,
    backend: str = "flint",
    *,
    rigorous: bool = False,
    n_m2_sample: int = 4000,
) -> Dict[str, Any]:
    """Structured even-block assembly at cutoff T."""
    _require_backend(backend)
    _, arb, _, _ = require_flint()
    set_precision_bits(precision_bits)
    L_a = arb(L)
    gp00, gp0b, gpbb = gp_even_block(L_a, arb)
    g000, g00b, g0bb = g0_even_block(L_a, arb)
    gi00, gi0b, gibb = ginf_even_block_quad(
        L_a,
        T,
        arb,
        n=n_quad,
        with_error_bound=rigorous,
        n_m2_sample=n_m2_sample,
    )
    G00 = g000 - gp00 + gi00
    G0b = g00b - gp0b + gi0b
    Gbb = g0bb - gpbb + gibb
    E2 = G00 * Gbb - G0b * G0b
    # Rank-1 / det-0 regression for the pole block alone
    pole_det = g000 * g0bb - g00b * g00b
    evidence = "E1" if rigorous else "E3_PENDING_QUADRATURE_BOUND"
    return {
        "G00": G00,
        "G0b": G0b,
        "Gbb": Gbb,
        "E2": E2,
        "G0_00": g000,
        "G0_0b": g00b,
        "G0_bb": g0bb,
        "pole_det": pole_det,
        "pole_candidate": POLE_CANDIDATE,
        "pole_formula": POLE_FORMULA,
        "normalization": NORMALIZATION,
        "cutoff_T": T,
        "n_quad": n_quad,
        "evidence_class": evidence,
        "rh_proof_claim": False,
        "note": (
            "True Weil integrand with even pole outer-product (√3/2)(v₊v₊ᵀ+v₋v₋ᵀ). "
            + (
                "Trapezoid remainder enclosed via max|f''| majorant."
                if rigorous
                else "Set rigorous=True to attach quadrature-error binder for E1."
            )
        ),
    }


def finite_weil_odd_pivot(
    L,
    T: int = 84,
    precision_bits: int = 256,
    n_quad: int = 65536,
    backend: str = "flint",
    *,
    rigorous: bool = False,
) -> Dict[str, Any]:
    """O1 at finite T: G0 - Gp + Ginf for q1 (Ginf via Hb relation)."""
    _require_backend(backend)
    _, arb, _, _ = require_flint()
    set_precision_bits(precision_bits)
    L_a = arb(L)
    g0 = g0_odd_pivot(L_a, arb)
    gp = gp_odd_pivot(L_a, arb)
    T_a = arb(T)
    o_inf = arb(0)
    n = n_quad
    for i in range(n + 1):
        t = T_a * i / n
        w = arb("0.5") if i in (0, n) else arb(1)
        hp = h_plus(t)
        if t.contains(0) or abs(t.mid()) <= t.rad():
            hq1sq = arb(0)
        else:
            _, _, hbsq = stable_products_even(t, L_a, arb)
            hq1sq = hbsq * (t**2) / 4
        o_inf += w * hp * hq1sq
    o_inf = o_inf * (T_a / n) / arb.pi()
    if rigorous:
        # |Hq1|² ≤ |Hb|² t²/4 ≤ 4 L² / t² for large t; on [0,T] use
        # remainder via M2 of odd integrand ≤ M2_bb * (sup t²/4 factor handled loosely):
        # reuse bb M2 scaled by max(t²/4)≤ T²/4 is too crude; use M2_0b-style:
        m2 = enclose_integrand_M2(L_a, float(T))[2]  # share decay with Hb
        # |Hq1|² = (t²/4)|Hb|² ≤ (T²/4) |Hb|² so M2_odd ≤ (T²/4) M2_bb + cross terms;
        # conservative: M2_odd ≤ (T**2) * (M2_bb + 1)
        m2_odd = arb(float(T) ** 2) * (m2 + 1)
        factor = T_a * ((T_a / n) ** 2) / 12
        o_inf = o_inf + arb(0, float((factor * m2_odd).upper()))
    O1 = g0 - gp + o_inf
    return {
        "O1": O1,
        "normalization": NORMALIZATION,
        "cutoff_T": T,
        "n_quad": n_quad,
        "evidence_class": "E1" if rigorous else "E3_PENDING_QUADRATURE_BOUND",
        "rh_proof_claim": False,
    }


def finite_weil_degree2(
    L,
    T: int = 84,
    precision_bits: int = 256,
    n_quad: int = 65536,
    backend: str = "flint",
    *,
    rigorous: bool = False,
) -> Dict[str, Any]:
    _require_backend(backend)
    even = finite_weil_even_block(
        L,
        T=T,
        precision_bits=precision_bits,
        n_quad=n_quad,
        backend=backend,
        rigorous=rigorous,
    )
    odd = finite_weil_odd_pivot(
        L,
        T=T,
        precision_bits=precision_bits,
        n_quad=n_quad,
        backend=backend,
        rigorous=rigorous,
    )
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
    _require_backend(backend)
    block = finite_weil_degree2(L, T=T, backend=backend)
    key = {
        ("1", "1"): "G00",
        ("1", "b"): "G0b",
        ("b", "1"): "G0b",
        ("b", "b"): "Gbb",
        ("q1", "q1"): "O1",
    }[(i, j)]
    return block[key]
