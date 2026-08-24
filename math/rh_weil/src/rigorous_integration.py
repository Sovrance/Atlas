"""Canonical rigorous quadrature for the RH/Weil program (ATLAS-RH-ENG-005 §2).

Every E1 numerical integral in this program goes through
:func:`rigorous_panel_integral`. Two reasons, both learned the hard way:

**One-shot integration is not reliable.** ``acb.integral`` over the whole range
in a single call exhausts its evaluation budget on this integrand and returns a
*non-finite ball* — which is sound (an infinite enclosure encloses everything)
but useless. Splitting the range into panels converges every time. So panels are
the canonical path, not a workaround applied when the one-shot version happens to
fail.

**The schedule must be deterministic and recorded.** A certificate that says
"integrated with Arb" cannot be re-derived; one that records its exact panel
edges can. :func:`panel_schedule` is a pure function of ``T``, and every
certificate stores the schedule it used.

The old global composite-trapezoid path (``finite_weil.ginf_even_block_quad``
with ``with_error_bound=True``) must not emit E1 for this program. Its rigorous
remainder uses a single global ``M2`` over the whole range and returns a radius
of ~2e4 on this integrand — six orders of magnitude larger than the quantity
being bounded. :func:`assert_not_trapezoid_path` exists so callers can state that
constraint in code.

No RH proof claim is made by this module.
"""
from __future__ import annotations

from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

#: The fixed T=84 schedule named in ENG-005 §2. Dyadic from 1, so each panel
#: spans one octave of the oscillation and the adaptive integrator never has to
#: resolve a wildly varying scale inside a single panel.
PANELS_T84: Tuple[Tuple[float, float], ...] = (
    (0.0, 1.0), (1.0, 2.0), (2.0, 4.0), (4.0, 8.0),
    (8.0, 16.0), (16.0, 32.0), (32.0, 64.0), (64.0, 84.0),
)

TRAPEZOID_PATH = "finite_weil.ginf_even_block_quad"


class QuadratureFailure(RuntimeError):
    """A panel did not converge to a finite enclosure."""


def panel_schedule(T: float) -> List[Tuple[float, float]]:
    """Deterministic panel edges for ``[0, T]``.

    ``T = 84`` returns exactly :data:`PANELS_T84`. Otherwise: ``[0,1]``, then
    dyadic doubling to the last power of two below ``T``, then a final panel to
    ``T``. Pure function of ``T`` — the same ``T`` always yields the same edges,
    which is what makes a recorded schedule reproducible.
    """
    if T <= 0:
        raise ValueError(f"T must be positive, got {T!r}")
    if float(T) == 84.0:
        return list(PANELS_T84)
    edges: List[float] = [0.0, 1.0]
    while edges[-1] * 2 < T:
        edges.append(edges[-1] * 2)
    if edges[-1] < T:
        edges.append(float(T))
    return [(lo, hi) for lo, hi in zip(edges, edges[1:]) if hi > lo]


def rigorous_panel_integral(
    integrand: Callable[[Any, bool], Any],
    T: float,
    acb,
    *,
    panels: Optional[Sequence[Tuple[float, float]]] = None,
    options: Optional[Dict[str, Any]] = None,
) -> Tuple[Any, Dict[str, Any]]:
    """Integrate ``integrand`` over ``[0, T]`` panel by panel, rigorously.

    ``integrand`` takes ``(z, analytic)`` as ``acb.integral`` requires and must
    return a finite ball for every ball it is handed — including balls containing
    a removable singularity, which is why the callers here carry series branches.

    Returns ``(value, record)``. ``record`` is what the certificate stores: the
    exact schedule, each panel's enclosure radius, and the integrator options.

    Raises :class:`QuadratureFailure` on a non-finite panel rather than returning
    an infinite enclosure — a bound derived from one would be vacuously true.
    """
    schedule = list(panels) if panels is not None else panel_schedule(T)
    if not schedule:
        raise ValueError("empty panel schedule")
    if abs(schedule[0][0]) > 0 or abs(schedule[-1][1] - float(T)) > 1e-12:
        raise ValueError(f"schedule {schedule[0][0]}..{schedule[-1][1]} does not cover [0, {T}]")

    total = acb(0)
    rows: List[Dict[str, Any]] = []
    for lo, hi in schedule:
        piece = acb.integral(integrand, lo, hi, **(options or {}))
        if not piece.is_finite():
            raise QuadratureFailure(
                f"panel [{lo}, {hi}] did not converge to a finite enclosure; "
                "a bound derived from an infinite ball would be vacuous"
            )
        total += piece
        rows.append({"lo": lo, "hi": hi, "radius": float(piece.real.rad())})

    record = {
        "method": "arb_acb_integral_panelled",
        "T": float(T),
        "panel_schedule": [[lo, hi] for lo, hi in schedule],
        "n_panels": len(schedule),
        "panel_radii": [r["radius"] for r in rows],
        "max_panel_radius": max((r["radius"] for r in rows), default=0.0),
        "options": dict(options or {}),
        "trapezoid_path_used": False,
    }
    return total, record


def assert_not_trapezoid_path(record: Dict[str, Any]) -> None:
    """Refuse a quadrature record that came from the rejected trapezoid path."""
    if record.get("trapezoid_path_used") or record.get("method", "").startswith("trapezoid"):
        raise QuadratureFailure(
            f"{TRAPEZOID_PATH} may not emit E1 for this program (ENG-005 §2): its "
            "rigorous remainder uses one global M2 and is ~1e6 times the quantity "
            "being bounded"
        )
    if record.get("method") != "arb_acb_integral_panelled":
        raise QuadratureFailure(f"unrecognised quadrature method {record.get('method')!r}")
