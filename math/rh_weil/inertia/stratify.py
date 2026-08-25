"""Inertia of a parametrised Hermitian family over a real cell (§3, §9 Outcome B).

The uniform certificates in ENG-005 asked one yes/no question -- is this scalar
positive everywhere on the cell -- and answered it by refusing to accept a box
until it cleared. Inertia is a different shape of question: the answer is a
*label*, it is allowed to change along the cell, and a change is a result rather
than a failure. So the driver here does not try to force one answer over the
whole interval. It subdivides until each piece carries a determined signature,
then merges neighbours that agree, and reports the leftover pieces as transition
regions.

That leftover is the honest part. A signature can only change where the matrix
is singular, and an interval enclosure can never certify singularity (§14.3), so
the cells bracketing a crossing will refuse to resolve no matter how far they are
split. Their width is the resolution to which the transition has been located,
and they are reported as ``INCONCLUSIVE_TRANSITION_REGION`` -- not silently
absorbed into whichever neighbour happens to be adjacent.

No RH proof claim is made by this module.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable, Dict, List, Optional, Tuple

from .ldl import INCONCLUSIVE, InertiaResult, interval_inertia

#: Default subdivision policy for :func:`certify_inertia_family`.
DEFAULT_POLICY: Dict[str, Any] = {
    "initial_cells": 16,
    "max_depth": 24,
    #: Below this width a cell stops splitting and is reported as a transition
    #: region. Without it, a genuine crossing would recurse until the split
    #: points stop being distinguishable as floats.
    "min_width": 1e-12,
}


@dataclass
class Stratum:
    lo: float
    hi: float
    signature: Tuple[int, int, int]
    cells: int = 1

    def to_dict(self) -> Dict[str, Any]:
        return {
            "interval": [repr(self.lo), repr(self.hi)],
            "n_positive": self.signature[0],
            "n_negative": self.signature[1],
            "n_zero": self.signature[2],
            "cells_merged": self.cells,
        }


@dataclass
class TransitionRegion:
    lo: float
    hi: float
    blocker: str
    depth: int

    def to_dict(self) -> Dict[str, Any]:
        return {
            "interval": [repr(self.lo), repr(self.hi)],
            "width": repr(self.hi - self.lo),
            "status": "INCONCLUSIVE_TRANSITION_REGION",
            "reached_depth": self.depth,
            "blocker": self.blocker,
        }


@dataclass
class InertiaStratification:
    status: str
    cell: Tuple[float, float]
    strata: List[Stratum] = field(default_factory=list)
    transitions: List[TransitionRegion] = field(default_factory=list)
    boxes_examined: int = 0
    max_depth: int = 0
    policy: Dict[str, Any] = field(default_factory=dict)

    @property
    def is_constant(self) -> bool:
        return not self.transitions and len(self.strata) == 1

    def signature_if_constant(self) -> Optional[Tuple[int, int, int]]:
        return self.strata[0].signature if self.is_constant else None

    def to_dict(self) -> Dict[str, Any]:
        return {
            "status": self.status,
            "cell": [repr(self.cell[0]), repr(self.cell[1])],
            "constant_inertia": self.is_constant,
            "strata": [s.to_dict() for s in self.strata],
            "transition_regions": [t.to_dict() for t in self.transitions],
            "boxes_examined": self.boxes_examined,
            "max_subdivision_depth": self.max_depth,
            "subdivision_policy": dict(self.policy),
            "method": "interval_ldl_congruence_with_adaptive_subdivision",
            "coverage": (
                "the strata and transition regions together tile the cell exactly; "
                "every point of the cell lies in one of them"
            ),
        }


def certify_inertia_family(
    matrix_fn: Callable[[float, float], Any],
    L_interval: Tuple[float, float],
    *,
    subdivision_policy: Optional[Dict[str, Any]] = None,
    evaluate: Callable[[Any], InertiaResult] = interval_inertia,
) -> InertiaStratification:
    """Stratify ``L_interval`` by the inertia of ``matrix_fn(lo, hi)``.

    ``matrix_fn(lo, hi)`` must return a Hermitian matrix of balls enclosing the
    family over the whole sub-interval ``[lo, hi]``.
    """
    policy = dict(DEFAULT_POLICY)
    if subdivision_policy:
        policy.update(subdivision_policy)
    a, b = L_interval
    n0 = int(policy["initial_cells"])
    # Cell edges use the same ``a + (b-a)*k/n`` form the ENG-005 covers use, so
    # neighbouring cells share an endpoint bit-for-bit and the tiling below is
    # exact rather than approximate. Results are sorted before merging, so the
    # order they come off the stack does not matter.
    stack: List[Tuple[float, float, int]] = [
        (a + (b - a) * k / n0, a + (b - a) * (k + 1) / n0, 0)
        for k in range(n0)
    ]
    resolved: List[Tuple[float, float, Any, int, Optional[str]]] = []
    examined = 0
    deepest = 0

    while stack:
        lo, hi, depth = stack.pop()
        examined += 1
        deepest = max(deepest, depth)
        res = evaluate(matrix_fn(lo, hi))
        if res.status == "PASS":
            resolved.append((lo, hi, res.signature, depth, None))
            continue
        if depth >= int(policy["max_depth"]) or (hi - lo) <= float(policy["min_width"]):
            resolved.append((lo, hi, None, depth, res.blocker or INCONCLUSIVE))
            continue
        mid = (lo + hi) / 2
        stack.append((mid, hi, depth + 1))
        stack.append((lo, mid, depth + 1))

    resolved.sort(key=lambda r: r[0])
    strata: List[Stratum] = []
    transitions: List[TransitionRegion] = []
    for lo, hi, sig, depth, blocker in resolved:
        if sig is None:
            transitions.append(TransitionRegion(lo, hi, blocker, depth))
            continue
        if strata and strata[-1].signature == sig and strata[-1].hi == lo:
            strata[-1].hi = hi
            strata[-1].cells += 1
        else:
            strata.append(Stratum(lo, hi, sig))

    status = "PASS" if not transitions else "PASS_WITH_TRANSITION_REGIONS"
    return InertiaStratification(
        status=status,
        cell=(a, b),
        strata=strata,
        transitions=transitions,
        boxes_examined=examined,
        max_depth=deepest,
        policy=policy,
    )
