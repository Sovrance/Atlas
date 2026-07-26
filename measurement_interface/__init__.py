"""M7 -- the measurement-interface layer (SPEC §2 layer 3, executable).

Public API re-exported from :mod:`.engine`. See
``docs/M7-measurement-interface-scope.md``.
"""

from .engine import (
    FDT_THRESHOLD,
    PROMOTABLE_ROUTES,
    apparatus_limited_gate,
    apply_transfer,
    censoring_audit,
    dual_route_recovery,
    effective_rank,
    injection_recovery,
    invert_transfer,
    is_promotable_route,
    truncated_normal_mean_mle,
)

__all__ = [
    "FDT_THRESHOLD",
    "PROMOTABLE_ROUTES",
    "apparatus_limited_gate",
    "apply_transfer",
    "censoring_audit",
    "dual_route_recovery",
    "effective_rank",
    "injection_recovery",
    "invert_transfer",
    "is_promotable_route",
    "truncated_normal_mean_mle",
]
