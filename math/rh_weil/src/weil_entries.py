"""Finite Weil Gram entries for any basis pair (ATLAS-RH-ENG-005 §4/§5/§7).

Generalises the scalar canary's assembly to the whole degree-2 basis
``{one, q1, b}``:

    G_ij(L) = G0_ij(L) - Gp_ij(L) + Ginf_ij(L)

* ``G0`` — the adjudicated Candidate-A pole, via ``pole.pole_gram_entry``. There
  is no other pole implementation in the tree.
* ``Gp_ij(L) = sum_{q=p^k, log q < L} (log p / sqrt q) K_ij(log q; L)`` with the
  exact polynomial kernels.
* ``Ginf_ij(L) = (1/pi) int_0^T h_+(t) Re(conj(H_i) H_j) dt``, with ``T = inf``
  for the cutoff-free entries (§4, §5) and ``T = 84`` for the direct-Fourier
  entries (§7, §10).

Analytic continuation
---------------------
``acb.integral`` needs an integrand analytic in the integration variable, and
``Re`` is not. For a **real** basis and real ``t``, ``conj(H_i(t)) = H_i(-t)``, so

    Re(conj(H_i(t)) H_j(t)) = [H_i(-z) H_j(z) + H_i(z) H_j(-z)] / 2

continues analytically in ``z``. Likewise ``h_+`` continues as
``(psi(1/4+iz/2) + psi(1/4-iz/2))/2 - log pi``.

``H_i(z; L) = int_0^L h_i(x) e^{izx} dx`` is computed by ``pole.poly_exp_integral``
with ``a = iz``. That routine switches to a series with an explicit remainder ball
when ``|aL|`` is small, which is exactly the region ``acb.integral`` probes near
``t = 0`` where the endpoint closed form is a ``0/0`` ball. Its magnitude helper
returns an **upper** bound on ``|z|``, so the remainder is sized by the ball's
supremum rather than its midpoint.

No RH proof claim is made by this module.
"""
from __future__ import annotations

import math
from typing import Any, Dict, List, Optional, Sequence, Tuple

import core
import pole
from rigorous_integration import rigorous_panel_integral

BASIS: Tuple[str, ...] = ("one", "q1", "b")

#: Exact prime kernels ``K_ij(a; L)``, as named in ENG-005 §4/§5.
_KERNELS = {
    ("one", "one"): core.kernel_00,
    ("one", "b"): core.kernel_0b,
    ("b", "b"): core.kernel_bb,
    ("q1", "q1"): core.kernel_q1q1,
    # ENG-006 §7: the odd degree-3 partner b3(x) = x(L-x)(x-L/2).
    ("q1", "b3"): core.kernel_q1b3,
    ("b3", "b3"): core.kernel_b3b3,
}


# --------------------------------------------------------------------------- #
# Prime block                                                                  #
# --------------------------------------------------------------------------- #
def prime_powers_below(L_value: float) -> List[Tuple[int, int]]:
    """``(q, p)`` for prime powers ``q = p^k`` with ``log q < L``."""
    cap = int(math.floor(math.exp(L_value)))
    out = []
    for p in range(2, cap + 1):
        if any(p % d == 0 for d in range(2, int(p**0.5) + 1)):
            continue
        q = p
        while q <= cap and math.log(q) < L_value:
            out.append((q, p))
            q *= p
    return sorted(out)


def kernel(i: str, j: str, a: Any, L: Any) -> Any:
    """``K_ij(a; L)`` — exact polynomial overlap kernel."""
    key = (i, j) if (i, j) in _KERNELS else (j, i)
    if key not in _KERNELS:
        raise KeyError(f"no kernel for the pair {(i, j)!r}")
    return _KERNELS[key](a, L)


def prime_entry(i: str, j: str, L, arb,
                prime_powers: Optional[Sequence[Tuple[int, int]]] = None) -> Any:
    """``Gp_ij(L) = sum_q (log p / sqrt q) K_ij(log q; L)``."""
    if prime_powers is None:
        prime_powers = prime_powers_below(float(L))
    total = arb(0)
    for q, p in prime_powers:
        w = arb(p).log() / arb(q).sqrt()
        total += w * kernel(i, j, arb(q).log(), L)
    return total


# --------------------------------------------------------------------------- #
# Analytic continuation of the archimedean integrand                           #
# --------------------------------------------------------------------------- #
def h_plus_analytic(z, arb, acb, log_pi=None):
    """``(psi(1/4+iz/2) + psi(1/4-iz/2))/2 - log pi`` — analytic in ``z``."""
    log_pi = arb.pi().log() if log_pi is None else log_pi
    quarter = acb(arb("0.25"))
    i = acb(0, 1)
    return ((quarter + i * z / 2).digamma() + (quarter - i * z / 2).digamma()) / 2 - log_pi


def H_transform(name: str, z, L, acb):
    """``H_i(z; L) = int_0^L h_i(x) e^{izx} dx`` (complex, analytic in ``z``)."""
    coeffs = pole.basis_coeffs(name, L)
    return pole.poly_exp_integral([acb(c) for c in coeffs], acb(0, 1) * z, acb(L))


def spectral_product(i: str, j: str, z, L, acb):
    """``Re(conj(H_i) H_j)`` continued analytically: ``[H_i(-z)H_j(z) + H_i(z)H_j(-z)]/2``.

    Real and equal to the real part on the real axis, because the basis is real
    so ``conj(H_i(t)) = H_i(-t)`` there.
    """
    hi_p, hi_m = H_transform(i, z, L, acb), H_transform(i, -z, L, acb)
    hj_p, hj_m = H_transform(j, z, L, acb), H_transform(j, -z, L, acb)
    return (hi_m * hj_p + hi_p * hj_m) / 2


def arch_entry(i: str, j: str, L, T: float, arb, acb, *,
               panels=None, options=None):
    """``(1/pi) int_0^T h_+(t) Re(conj(H_i)H_j) dt`` as a rigorous enclosure.

    Returns ``(value, quadrature_record)``.
    """
    log_pi = arb.pi().log()
    L_a = arb(L) if not hasattr(L, "mid") else L

    def integrand(z, _analytic):
        return h_plus_analytic(z, arb, acb, log_pi) * spectral_product(i, j, z, L_a, acb)

    total, record = rigorous_panel_integral(integrand, T, acb,
                                            panels=panels, options=options)
    return total.real / arb.pi(), record


# --------------------------------------------------------------------------- #
# Assembly                                                                     #
# --------------------------------------------------------------------------- #
def gram_entry_truncated(i: str, j: str, L, T: float, arb, acb, *,
                         prime_powers=None, panels=None, options=None):
    """``G_ij`` with the archimedean integral truncated at ``T``.

    For the cutoff-free entries the caller adds the tail lemma's contribution;
    for the T=84 entries this *is* the object of study (§7, §10).
    """
    L_a = arb(L) if not hasattr(L, "mid") else L
    arch, record = arch_entry(i, j, L_a, T, arb, acb, panels=panels, options=options)
    g0 = pole.pole_gram_entry(i, j, L_a)
    gp = prime_entry(i, j, L_a, arb, prime_powers)
    return g0 - gp + arch, record


def block_truncated(L, T: float, arb, acb, *, prime_powers=None,
                    panels=None, options=None) -> Dict[str, Any]:
    """The four entries the program needs, plus the derived determinants.

    ``E2 = G00 Gbb - G0b^2`` is the compact even determinant; ``O1 = G[q1,q1]``
    is the odd pivot. ``det(G_deg<=2) = O1 * E2`` by the parity factorization —
    the pole and prime blocks are both parity block diagonal, so the full Gram
    splits and the determinant is the product of the two blocks.
    """
    if prime_powers is None:
        prime_powers = prime_powers_below(float(arb(L).mid()) if hasattr(L, "mid") else float(L))
    entries: Dict[str, Any] = {}
    records: Dict[str, Any] = {}
    for key, (i, j) in (("G00", ("one", "one")), ("G0b", ("one", "b")),
                        ("Gbb", ("b", "b")), ("O1", ("q1", "q1"))):
        val, rec = gram_entry_truncated(i, j, L, T, arb, acb,
                                        prime_powers=prime_powers,
                                        panels=panels, options=options)
        entries[key] = val
        records[key] = rec
    entries["E2"] = entries["G00"] * entries["Gbb"] - entries["G0b"] ** 2
    entries["det_deg2"] = entries["O1"] * entries["E2"]
    entries["_quadrature"] = records
    return entries
