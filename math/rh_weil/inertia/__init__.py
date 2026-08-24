"""Certified inertia for finite Hermitian matrices (ATLAS-RH-ENG-006 §3/§4).

Positivity is one bit of a signature. This package certifies the whole thing --
``(n_positive, n_negative, n_zero)`` -- so a finite Weil block that fails to be
positive still yields a rigorous result instead of a dead end.

No RH proof claim is made by this package.
"""
from .certificate import (  # noqa: F401
    KIND_INERTIA,
    KIND_STRATIFICATION,
    build_inertia_certificate,
    build_stratification_certificate,
    satisfies_psd_requirement,
    validate_against_schema,
)
# NB: the bare name ``congruence`` is deliberately not re-exported here -- it
# would shadow the ``inertia.congruence`` submodule of the same name, so
# ``from inertia import congruence`` would hand back a function where callers
# reasonably expect the module. Reach it as ``inertia.congruence.congruence``.
from .congruence import charpoly_inertia, inertia_2x2  # noqa: F401
from .ldl import (  # noqa: F401
    INCONCLUSIVE,
    BallSignOracle,
    ExactSignOracle,
    InertiaResult,
    exact_inertia,
    interval_inertia,
    ldl_inertia,
)
from .stratify import (  # noqa: F401
    DEFAULT_POLICY,
    InertiaStratification,
    Stratum,
    TransitionRegion,
    certify_inertia_family,
)


def certify_inertia(matrix_ball, *, backend="arb", precision_bits=None):
    """Inertia of one Hermitian matrix (§3 primary API).

    ``backend="arb"`` treats the entries as balls and never claims exact zero;
    ``backend="exact"`` treats them as rationals and does.
    """
    if backend == "exact":
        return exact_inertia(matrix_ball)
    if backend != "arb":
        raise ValueError(f"unknown backend {backend!r}")
    if precision_bits is not None:
        from interval_backend import set_precision_bits

        set_precision_bits(precision_bits)
    return interval_inertia(matrix_ball)
