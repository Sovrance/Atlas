"""Frozen normalization for the finite Weil quadratic form (WO-RH-17).

This module is the single mathematical source of truth for the finite Weil
Gram assembly. Every convention is stated explicitly; nothing here is fitted.

Adjudication summary (see ``docs/NORMALIZATION_ADJUDICATION_v0.1.md``):

* **Candidate A (ADOPTED)** — the pole block obtained from the explicit formula,

      G0_ij = Fhat_ij(i/2) + Fhat_ij(-i/2) = E_i^+ E_j^- + E_i^- E_j^+,
      E_i^± = int_0^L h_i(x) e^{±x/2} dx.

* **Candidate B (REJECTED)** — the legacy even block ``(sqrt(3)/2)(v+v+^T +
  v-v-^T)``. On the even sector it equals Candidate A multiplied by the
  *L-dependent* factor ``(sqrt(3)/2)cosh(L/2)``, which is 1 **only** at
  ``L = log 3``. It is therefore a multiplicative calibration fitted at a single
  test point, not a normalization convention, and is retained here only so the
  audit is reproducible.

No RH proof claim is made by this module.
"""
from __future__ import annotations

from math import cosh, exp, factorial, log, sqrt
from typing import Any, Callable, Dict, Sequence, Tuple

import pole

NORMALIZATION_VERSION = "v0.1"

# --------------------------------------------------------------------------- #
# 1. Conventions (all seven axes the work order requires to be pinned)         #
# --------------------------------------------------------------------------- #
FOURIER_CONVENTION = "Fhat(xi) = int_R F(x) e^{-i*xi*x} dx"
TILDE_CONVENTION = "tilde_h(x) = conj(h(-x)); real basis => tilde_h(x) = h(-x)"
CONVOLUTION_CONVENTION = "F_ij = h_i * tilde_h_j, supported on [-L, L]"
POLE_FORMULA = (
    "G0_ij = Fhat_ij(i/2) + Fhat_ij(-i/2) = E_i^+ E_j^- + E_i^- E_j^+, "
    "E_i^± = int_0^L h_i(x) e^{±x/2} dx"
)
PRIME_FORMULA = (
    "Gp_ij = sum_{q=p^k, log q < L} (log p / sqrt(q)) * K_ij(log q; L), "
    "K_ij(a;L) = int_0^{L-a} [h_i(x)h_j(x+a) + h_j(x)h_i(x+a)] dx"
)
ARCH_FORMULA = (
    "Ginf_ij(T) = (1/pi) int_0^T h_+(t) Re(conj(H_i(t;L)) H_j(t;L)) dt, "
    "one-sided frequency integral with the even/real basis"
)
ASSEMBLY = "G = G0 - Gp + Ginf"
FREQUENCY_INTEGRAL = "one-sided (0,T], real part taken; even basis"
BASIS_NORMALIZATION = "unnormalised polynomials on [0,L]: 1, q1 = x - L/2, b = x(L-x)"

CLAIM_BOUNDARY = "finite_block_only_no_rh_proof"
RH_PROOF_CLAIM = False

BASIS_NAMES: Tuple[str, ...] = ("one", "q1", "b")


# ENG-004 §1: basis + parity are owned by ``pole.py``.
basis_coeffs = pole.basis_coeffs
basis_parity = pole.basis_parity


# --------------------------------------------------------------------------- #
# 2. E^± integrals — exact closed form for any polynomial basis element        #
# --------------------------------------------------------------------------- #
# ENG-004 §1: the quadrature primitives live in ``pole.py``; re-exported here so
# existing importers keep working without a second implementation.
monomial_exp_integral = pole.monomial_exp_integral
poly_exp_integral = pole.poly_exp_integral


def E_pm(name: str, L: Any, sign: int) -> Any:
    """``E^±_i = int_0^L h_i(x) e^{±x/2} dx`` in closed form."""
    if sign == 1:
        return pole.laplace_plus(name, L)
    if sign == -1:
        return pole.laplace_minus(name, L)
    raise ValueError("sign must be +1 or -1")


def H_transform(name: str, t: Any, L: Any) -> Any:
    """``H_i(t;L) = int_0^L h_i(x) e^{i t x} dx`` in closed form (complex)."""
    return poly_exp_integral(basis_coeffs(name, L), 1j * t, L)


# --------------------------------------------------------------------------- #
# 3. Candidate A — the ADOPTED pole block (delegated to the canonical module)  #
# --------------------------------------------------------------------------- #
# ENG-004 §1: ``src/pole.py`` is the single implementation of the pole formula.
# This name is kept because the adjudication artifacts and cross-check reference
# it, but it must never grow a second copy of the formula.
def pole_entry(i: str, j: str, L: Any) -> Any:
    """``G0_ij = E_i^+ E_j^- + E_i^- E_j^+`` (Candidate A, adopted)."""
    return pole.pole_gram_entry(i, j, L)


# --------------------------------------------------------------------------- #
# 4. Candidate B is NOT here                                                   #
# --------------------------------------------------------------------------- #
# The rejected ``(sqrt(3)/2)`` block lives in ``src/rejected_pole.py``, which is
# archival and forbidden to production (tests/test_production_imports.py). Do not
# re-import it here: this module is production.
LEGACY_POLE_MODULE = "rejected_pole"
LEGACY_STATUS = "REJECTED_FITTED_CALIBRATION"
CALIBRATION_FIXED_POINT = "L = log 3"


# --------------------------------------------------------------------------- #
# 5. Prime block — K_ij kernel (shared real-space correlation)                 #
# --------------------------------------------------------------------------- #
def _poly_mul(p: Sequence[Any], q: Sequence[Any]) -> list:
    out = [0 * (p[0] + q[0])] * (len(p) + len(q) - 1)
    out = list(out)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            out[i + j] = out[i + j] + a * b
    return out


def _poly_shift(p: Sequence[Any], a: Any) -> list:
    """Coefficients of p(x + a)."""
    from math import comb

    out = [0 * (p[0] + a)] * len(p)
    out = list(out)
    for n, c in enumerate(p):
        for r in range(n + 1):
            out[r] = out[r] + c * comb(n, r) * (a ** (n - r))
    return out


def _poly_integral(p: Sequence[Any], upper: Any) -> Any:
    total = 0 * upper
    for n, c in enumerate(p):
        total = total + c * (upper ** (n + 1)) / (n + 1)
    return total


def kernel_K(i: str, j: str, a: Any, L: Any) -> Any:
    """``K_ij(a;L) = int_0^{L-a}[h_i(x)h_j(x+a) + h_j(x)h_i(x+a)]dx`` (exact)."""
    hi, hj = basis_coeffs(i, L), basis_coeffs(j, L)
    upper = L - a
    term1 = _poly_integral(_poly_mul(hi, _poly_shift(hj, a)), upper)
    term2 = _poly_integral(_poly_mul(hj, _poly_shift(hi, a)), upper)
    return term1 + term2


def prime_powers_below(L_value: float):
    """Exact ``(q, p, log q)`` for prime powers with ``log q < L``."""
    import math

    c = max(2, int(math.floor(math.exp(float(L_value)) + 1e-12)))
    sieve = [True] * (c + 1)
    if c >= 0:
        sieve[0] = False
    if c >= 1:
        sieve[1] = False
    for i in range(2, int(c**0.5) + 1):
        if sieve[i]:
            for m in range(i * i, c + 1, i):
                sieve[m] = False
    out = []
    for p in (i for i in range(2, c + 1) if sieve[i]):
        q = p
        while q <= c:
            lq = math.log(q)
            if lq < float(L_value):
                out.append((q, p, lq))
            if q > c // p:
                break
            q *= p
    return sorted(out)


def prime_entry(i: str, j: str, L: Any) -> Any:
    """``Gp_ij = sum_q (log p / sqrt q) K_ij(log q; L)``."""
    total = 0 * L
    for q, p, _ in prime_powers_below(float(L)):
        w = log(p) / sqrt(q)
        total = total + w * kernel_K(i, j, log(q), L)
    return total


# --------------------------------------------------------------------------- #
# 6. Identity used by the independent real-space pole route                    #
# --------------------------------------------------------------------------- #
POLE_REALSPACE_IDENTITY = (
    "G0_ij = int_0^L K_ij(a;L) * 2*cosh(a/2) da  "
    "(same K_ij as the prime block; equals E_i^+E_j^- + E_i^-E_j^+)"
)


def normalization_content() -> Dict[str, Any]:
    """The frozen, content-addressable normalization definition."""
    return {
        "normalization_version": NORMALIZATION_VERSION,
        "fourier_convention": FOURIER_CONVENTION,
        "tilde_convention": TILDE_CONVENTION,
        "convolution_convention": CONVOLUTION_CONVENTION,
        "pole_formula": POLE_FORMULA,
        "pole_realspace_identity": POLE_REALSPACE_IDENTITY,
        "prime_formula": PRIME_FORMULA,
        "archimedean_formula": ARCH_FORMULA,
        "assembly": ASSEMBLY,
        "frequency_integral": FREQUENCY_INTEGRAL,
        "basis_normalization": BASIS_NORMALIZATION,
        "claim_boundary": CLAIM_BOUNDARY,
        "rh_proof_claim": RH_PROOF_CLAIM,
    }


def normalization_id() -> str:
    """Content-addressed id of the frozen normalization."""
    import hashlib
    import json

    blob = json.dumps(normalization_content(), sort_keys=True, separators=(",", ":"))
    return "norm_sha256_" + hashlib.sha256(blob.encode("utf-8")).hexdigest()[:32]


# --------------------------------------------------------------------------- #
# 6. Quarantine registry (WO-RH-17 §3.4) — single source of truth              #
# --------------------------------------------------------------------------- #
# Certificates whose numbers were produced under the REJECTED even pole block.
# Kept here rather than in ``scripts/quarantine_normalization.py`` so that the
# certificate writer can enforce the quarantine at the point of write: a certify
# script that regenerates one of these files must not be able to silently drop
# the marker (the same hazard already fixed for ``work_order_status.json``).
QUARANTINE_STATE = "QUARANTINED_NORMALIZATION_ADJUDICATION"

QUARANTINED_CERTIFICATES: Tuple[str, ...] = (
    "e1_scalar_log3_log4.json",
    "e1_degree1_log3_log4.json",
    "e1_degree2_compact_log3_log4.json",
    "e1_fourier_T84_points.json",
    "e1_fourier_T84_uniform_degree2.json",
)

QUARANTINE_REASON = (
    "Depends on the even pole block under adjudication (WO-RH-17). The repository "
    "block (sqrt(3)/2)(v+v+^T+v-v-^T) equals the explicit-formula pole times "
    "(sqrt(3)/2)cosh(L/2), which is 1 only at L = log 3; values at any other L "
    "carry that factor. Regenerate under the frozen normalization (WO-RH-19/20) "
    "before any promotion."
)


#: Released from the WO-RH-17 quarantine by a later work order, after a rigorous
#: regeneration under the active normalization. A released file stays in
#: ``QUARANTINED_CERTIFICATES``: the legacy certifier that produced the disputed
#: version still exists, so an *unauthorised* write must still fail closed. Only a
#: body that carries ``quarantine_released`` and passes the promotion predicate is
#: left alone (see ``scripts/quarantine_normalization.py``).
RELEASED_CERTIFICATES: Dict[str, str] = {
    "e1_scalar_log3_log4.json": "ATLAS-RH-ENG-004 §4 scalar canary",
}


def is_quarantined_certificate(name: str) -> bool:
    """True if ``name`` must carry the WO-RH-17 quarantine marker when written."""
    return name in QUARANTINED_CERTIFICATES


def quarantine_block(body: Dict[str, Any]) -> Dict[str, Any]:
    """Build the ``quarantine`` sub-object, preserving ``body``'s claim as prior state."""
    return {
        "reason": QUARANTINE_REASON,
        "work_order": "WO-RH-17",
        "adjudication_certificate": "normalization_adjudication.json",
        "active_normalization_id": normalization_id(),
        "prior_state": {
            "hard_constraints_certified": body.get("hard_constraints_certified"),
            "status": body.get("status"),
        },
        "pre_quarantine_content_hash": body.get("content_hash"),
        "evidence_class_preserved": body.get("evidence_class"),
        "note": "not relabelled E3 on purpose; the historical claim is preserved as evidence",
        "numerically_unaffected_at": "L = log 3 (the calibration fixed point) — still "
                                     "not promotable until regenerated",
    }
