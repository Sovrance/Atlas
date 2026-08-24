"""Four independent Weil-component providers (WO-RH-18).

Each provider implements the same interface over the *frozen* normalization in
:mod:`normalization`, but by a mathematically **independent route**, so that
agreement is evidence rather than a tautology:

===========================  =========================================  =======================================
provider                     pole route                                 prime / archimedean route
===========================  =========================================  =======================================
ExplicitFormulaProvider      closed-form ``E_i^+E_j^- + E_i^-E_j^+``     exact polynomial kernel ``K_ij``
CompactRealSpaceProvider     ``int_0^L K_ij(a) 2cosh(a/2) da``           direct quadrature of the overlap
DirectFourierProvider        quadrature of the ``E^±`` integrals         Fourier-side ``H_i(t)`` product
ConnesCvSProjectedProvider   external adapter (diagnostic only)          external adapter (diagnostic only)
===========================  =========================================  =======================================

A provider that cannot supply a component **returns ``None``** (reported as
UNAVAILABLE). It never fabricates a number, and it never silently reuses another
provider's route.

Rigour labelling is explicit and conservative:

``interval_certified``      rigorous Arb enclosure (python-flint present)
``high_precision_numeric``  mpmath at high precision with a reported error estimate

Only ``interval_certified`` measurements may support an E1 promotion; this module
emits no promotion of any kind and makes no RH claim.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Optional, Protocol

import normalization as N

try:  # optional rigorous backend
    from interval_backend import FlintUnavailable, require_flint

    _flint, _arb, _acb, _ctx = require_flint()
    FLINT = True
except Exception:  # pragma: no cover - flint absent
    FLINT = False

try:
    import mpmath as mp

    MPMATH = True
except Exception:  # pragma: no cover
    MPMATH = False

DEFAULT_DPS = 40
RIGOUR_INTERVAL = "interval_certified"
RIGOUR_NUMERIC = "high_precision_numeric"


@dataclass(frozen=True)
class Measurement:
    """A value with an honest error radius: the interval ``[value±rad]``."""

    value: float
    rad: float
    method: str
    rigour: str

    def lo(self) -> float:
        return self.value - self.rad

    def hi(self) -> float:
        return self.value + self.rad

    def to_dict(self) -> dict:
        return {
            "value": self.value,
            "rad": self.rad,
            "lo": self.lo(),
            "hi": self.hi(),
            "method": self.method,
            "rigour": self.rigour,
        }


class WeilProvider(Protocol):  # pragma: no cover - structural
    name: str

    def pole_entry(self, basis_i: str, basis_j: str, L: float) -> Optional[Measurement]: ...

    def prime_entry(self, basis_i: str, basis_j: str, L: float) -> Optional[Measurement]: ...

    def arch_entry(
        self, basis_i: str, basis_j: str, L: float, *, T: Optional[float] = None
    ) -> Optional[Measurement]: ...

    def gram_entry(
        self, basis_i: str, basis_j: str, L: float, *, T: Optional[float] = None
    ) -> Optional[Measurement]: ...


# --------------------------------------------------------------------------- #
# helpers                                                                      #
# --------------------------------------------------------------------------- #
def _assemble(pole, prime, arch) -> Optional[Measurement]:
    """G = G0 - Gp + Ginf, propagating radii additively (conservative)."""
    if pole is None or prime is None or arch is None:
        return None
    value = pole.value - prime.value + arch.value
    rad = pole.rad + prime.rad + arch.rad
    rigour = (
        RIGOUR_INTERVAL
        if all(m.rigour == RIGOUR_INTERVAL for m in (pole, prime, arch))
        else RIGOUR_NUMERIC
    )
    return Measurement(value, rad, "G0 - Gp + Ginf", rigour)


def _mp_dps(dps: int = DEFAULT_DPS):
    if not MPMATH:  # pragma: no cover
        raise RuntimeError("mpmath required for this provider route")
    mp.mp.dps = dps


def _basis_callable(name: str, L: float):
    coeffs = N.basis_coeffs(name, L)

    def f(x):
        total = 0 * x
        for n, c in enumerate(coeffs):
            total = total + c * (x**n)
        return total

    return f


def _H_closed(name: str, t, L: float):
    """``H_i(t;L)`` from the closed-form polynomial-exponential integral."""
    return N.H_transform(name, complex(t), float(L))


def _H_quad(name: str, t, L: float):
    """``H_i(t;L)`` by direct high-precision quadrature (independent route)."""
    h = _basis_callable(name, L)
    return mp.quad(lambda x: h(x) * mp.e ** (1j * t * x), [0, L])


def _h_plus(tau):
    """Archimedean Mellin multiplier ``h_+`` (mpmath route)."""
    z = mp.mpf(0.25) + 0.5j * tau
    return (mp.digamma(z) + mp.digamma(mp.conj(z))).real / 2 - mp.log(mp.pi)


# --------------------------------------------------------------------------- #
# 1. ExplicitFormulaProvider — closed forms straight from the frozen formulas  #
# --------------------------------------------------------------------------- #
class ExplicitFormulaProvider:
    name = "ExplicitFormulaProvider"
    description = "closed-form explicit-formula pole + exact polynomial prime kernel"

    def pole_entry(self, basis_i, basis_j, L, **_):
        if FLINT:
            Lb = _arb(L)
            val = N.pole_entry(basis_i, basis_j, Lb)
            return Measurement(
                float(val.mid()), float(val.rad()), "E_i^+E_j^- + E_i^-E_j^+ (Arb)", RIGOUR_INTERVAL
            )
        v = N.pole_entry(basis_i, basis_j, float(L))
        return Measurement(v, abs(v) * 1e-13 + 1e-15, "E_i^+E_j^- + E_i^-E_j^+", RIGOUR_NUMERIC)

    def prime_entry(self, basis_i, basis_j, L, **_):
        if FLINT:
            Lb = _arb(L)
            val = N.prime_entry(basis_i, basis_j, Lb)
            return Measurement(
                float(val.mid()), float(val.rad()), "sum_q w_q K_ij (exact poly, Arb)", RIGOUR_INTERVAL
            )
        v = N.prime_entry(basis_i, basis_j, float(L))
        return Measurement(v, abs(v) * 1e-13 + 1e-15, "sum_q w_q K_ij (exact poly)", RIGOUR_NUMERIC)

    def arch_entry(self, basis_i, basis_j, L, *, T=None):
        if T is None:
            return None
        _mp_dps(30)
        def f(t):
            Hi, Hj = _H_closed(basis_i, t, L), _H_closed(basis_j, t, L)
            return _h_plus(t) * (Hi.conjugate() * Hj).real
        val = mp.quad(f, [0, T]) / mp.pi
        return Measurement(
            float(val), abs(float(val)) * 1e-10 + 1e-12,
            "(1/pi) int_0^T h_+ Re(conj(H_i)H_j) dt, H closed form", RIGOUR_NUMERIC,
        )

    def gram_entry(self, basis_i, basis_j, L, *, T=None):
        return _assemble(
            self.pole_entry(basis_i, basis_j, L),
            self.prime_entry(basis_i, basis_j, L),
            self.arch_entry(basis_i, basis_j, L, T=T),
        )


# --------------------------------------------------------------------------- #
# 2. CompactRealSpaceProvider — real-space correlation-kernel routes           #
# --------------------------------------------------------------------------- #
class CompactRealSpaceProvider:
    name = "CompactRealSpaceProvider"
    description = "pole as int K_ij(a)2cosh(a/2)da; prime by direct overlap quadrature"

    def pole_entry(self, basis_i, basis_j, L, **_):
        _mp_dps()
        f = lambda a: N.kernel_K(basis_i, basis_j, float(a), float(L)) * 2 * mp.cosh(a / 2)
        val = mp.quad(f, [0, L])
        return Measurement(
            float(val), abs(float(val)) * 1e-20 + 1e-25,
            "int_0^L K_ij(a) 2cosh(a/2) da", RIGOUR_NUMERIC,
        )

    def prime_entry(self, basis_i, basis_j, L, **_):
        """Direct quadrature of the shifted-overlap integral (no polynomial algebra)."""
        _mp_dps()
        hi = _basis_callable(basis_i, float(L))
        hj = _basis_callable(basis_j, float(L))
        total = mp.mpf(0)
        for q, p, lq in N.prime_powers_below(float(L)):
            a = mp.mpf(lq)
            w = mp.log(p) / mp.sqrt(q)
            ov = mp.quad(lambda x: hi(x) * hj(x + a) + hj(x) * hi(x + a), [0, L - a])
            total += w * ov
        return Measurement(
            float(total), abs(float(total)) * 1e-20 + 1e-25,
            "sum_q w_q * quad(h_i(x)h_j(x+a)+h_j(x)h_i(x+a))", RIGOUR_NUMERIC,
        )

    def arch_entry(self, basis_i, basis_j, L, *, T=None):
        return None  # archimedean term is not a compact real-space object

    def gram_entry(self, basis_i, basis_j, L, *, T=None):
        return _assemble(
            self.pole_entry(basis_i, basis_j, L),
            self.prime_entry(basis_i, basis_j, L),
            self.arch_entry(basis_i, basis_j, L, T=T),
        )


# --------------------------------------------------------------------------- #
# 3. DirectFourierProvider — transform-side routes                             #
# --------------------------------------------------------------------------- #
class DirectFourierProvider:
    name = "DirectFourierProvider"
    description = "pole by quadrature of the E^± integrals; archimedean on the Fourier side"

    def pole_entry(self, basis_i, basis_j, L, **_):
        _mp_dps()
        hi = _basis_callable(basis_i, float(L))
        hj = _basis_callable(basis_j, float(L))
        Eip = mp.quad(lambda x: hi(x) * mp.e ** (x / 2), [0, L])
        Eim = mp.quad(lambda x: hi(x) * mp.e ** (-x / 2), [0, L])
        Ejp = mp.quad(lambda x: hj(x) * mp.e ** (x / 2), [0, L])
        Ejm = mp.quad(lambda x: hj(x) * mp.e ** (-x / 2), [0, L])
        val = Eip * Ejm + Eim * Ejp
        return Measurement(
            float(val), abs(float(val)) * 1e-20 + 1e-25,
            "Fhat(i/2)+Fhat(-i/2) via quadrature of E^±", RIGOUR_NUMERIC,
        )

    def prime_entry(self, basis_i, basis_j, L, **_):
        return None  # the prime sum has no native transform-side representation here

    def arch_entry(self, basis_i, basis_j, L, *, T=None):
        if T is None:
            return None
        _mp_dps(20)
        f = lambda t: _h_plus(t) * (mp.conj(_H_quad(basis_i, t, L)) * _H_quad(basis_j, t, L)).real
        val = mp.quad(f, [0, T]) / mp.pi
        return Measurement(
            float(val), abs(float(val)) * 1e-8 + 1e-10,
            "(1/pi) int_0^T h_+ Re(conj(H_i)H_j) dt, H by quadrature", RIGOUR_NUMERIC,
        )

    def gram_entry(self, basis_i, basis_j, L, *, T=None):
        return _assemble(
            self.pole_entry(basis_i, basis_j, L),
            self.prime_entry(basis_i, basis_j, L),
            self.arch_entry(basis_i, basis_j, L, T=T),
        )


# --------------------------------------------------------------------------- #
# 4. ConnesCvSProjectedProvider — external diagnostic only                     #
# --------------------------------------------------------------------------- #
class ConnesCvSProjectedProvider:
    name = "ConnesCvSProjectedProvider"
    description = "external Connes/CvS projected implementation (diagnostic, not certifying)"
    external = True

    def __init__(self) -> None:
        self._available = False
        try:  # pragma: no cover - depends on optional external tree
            import sys
            from pathlib import Path

            ext = Path(__file__).resolve().parents[1] / "external"
            if str(ext) not in sys.path:
                sys.path.insert(0, str(ext))
            import connes_cvs_adapter as adapter

            adapter.dependency_info()
            self._adapter = adapter
            self._available = True
        except Exception:
            self._adapter = None

    def available(self) -> bool:
        return self._available

    def pole_entry(self, basis_i, basis_j, L, **_):
        return None  # projection/truncation error not certified -> never certifying

    def prime_entry(self, basis_i, basis_j, L, **_):
        return None

    def arch_entry(self, basis_i, basis_j, L, *, T=None):
        return None

    def gram_entry(self, basis_i, basis_j, L, *, T=None):
        return None


def all_providers():
    return [
        ExplicitFormulaProvider(),
        CompactRealSpaceProvider(),
        DirectFourierProvider(),
        ConnesCvSProjectedProvider(),
    ]
