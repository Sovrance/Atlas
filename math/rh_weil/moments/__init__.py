"""Spectral moments and the Atlas B1 adapter (ATLAS-RH-ENG-006 §6).

Low-order moments are a second certified channel alongside inertia: cheap to
enclose, and sometimes enough to recover the signature outright. When they are
not enough, that is reported as such rather than papered over.

No RH proof claim is made by this package.
"""
from .adapter import analyse  # noqa: F401
from .feasible_spectrum import CONCLUSIVE, INSUFFICIENT  # noqa: F401
from .spectral_moments import (  # noqa: F401
    KIND_SPECTRAL_MOMENT,
    moment_report,
    spectral_moments,
    trace_of_power,
)
