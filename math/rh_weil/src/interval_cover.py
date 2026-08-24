"""Generic adaptive interval cover for uniform bounds (ATLAS-RH-ENG-005 §8).

One branch-and-bound used by every uniform certificate in the program, so the
cutoff-free entries (§4/§5) and the direct-Fourier T=84 entries (§10) are
certified by the same code and report the same statistics.

The routine takes an evaluator ``box -> enclosure`` and nothing else. It does not
know or assume anything about the shape of the function: no convexity, no
monotonicity, no single minimum. That is deliberate — §8 forbids precommitting to
a topology, so what the certificate reports is the topology the cover *actually*
established: a uniform lower bound over an exhaustive set of boxes, with the box
count, the deepest split, and where the binding box sits.

A box whose enclosure clears the target is accepted. One that does not is split.
Failure to separate after ``max_depth`` splits is raised rather than papered
over — that is the ENG-005 §15 stop condition, and returning a bound that a
narrow region violates would be worse than reporting nothing.

No RH proof claim is made by this module.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, Callable, Dict, List, Optional, Tuple


class NotSeparated(ValueError):
    """The quantity could not be separated from the target on some sub-box."""


@dataclass
class CoverResult:
    quantity: str
    certified_lower_bound: float
    target: float
    boxes_examined: int
    boxes_accepted: int
    max_depth: int
    min_box: Tuple[float, float]
    min_box_enclosure: Tuple[float, float]
    max_width: float
    initial_boxes: int
    cell: Tuple[float, float]
    #: Best lower bound over every box not wholly contained in ``exclude``, so it
    #: bounds the quantity on ``cell \ exclude``. ``inf`` when no exclusion was
    #: requested, or when every box lies inside it.
    lower_bound_outside: float = math.inf
    exclude: Optional[Tuple[float, float]] = None
    extra: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        out = {
            "quantity": self.quantity,
            "certified_lower_bound": repr(self.certified_lower_bound),
            "target": repr(self.target),
            "cell": [repr(self.cell[0]), repr(self.cell[1])],
            "boxes_examined": self.boxes_examined,
            "boxes_accepted": self.boxes_accepted,
            "max_subdivision_depth": self.max_depth,
            "binding_box": [repr(self.min_box[0]), repr(self.min_box[1])],
            "binding_box_enclosure": [repr(self.min_box_enclosure[0]),
                                      repr(self.min_box_enclosure[1])],
            "max_box_enclosure_width": repr(self.max_width),
            "initial_boxes": self.initial_boxes,
            "method": "adaptive_interval_branch_and_bound",
            **({"exclude": [repr(self.exclude[0]), repr(self.exclude[1])],
                "certified_lower_bound_outside_exclude":
                    repr(self.lower_bound_outside)} if self.exclude else {}),
            "topology_proved": (
                "uniform lower bound over an exhaustive box cover of the closed "
                "cell; no convexity or monotonicity assumed"
            ),
        }
        out.update(self.extra)
        return out


def adaptive_cover(
    evaluate: Callable[[float, float], Tuple[float, float]],
    *,
    quantity: str,
    cell: Tuple[float, float],
    target: float = 0.0,
    initial_boxes: int = 24,
    max_depth: int = 22,
    exclude: Optional[Tuple[float, float]] = None,
) -> CoverResult:
    """Certify ``evaluate > target`` uniformly on ``cell``.

    ``evaluate(lo, hi)`` returns a rigorous ``(lower, upper)`` enclosure of the
    quantity over ``[lo, hi]``.

    ``exclude`` additionally reports the best bound over every box *not* wholly
    contained in that sub-interval — so the reported number bounds the quantity
    on all of ``cell`` minus ``exclude``, edge-straddling boxes included. Every box is still examined and still has to clear
    the target — the exclusion changes what is *reported*, never what is proved.
    It exists so a region governed by a separate, sharper argument can be taken
    out of the headline number without weakening the cover.
    """
    a, b = cell
    stack: List[Tuple[float, float, int]] = [
        (a + (b - a) * k / initial_boxes, a + (b - a) * (k + 1) / initial_boxes, 0)
        for k in range(initial_boxes)
    ]
    examined = accepted = 0
    deepest = 0
    best_lower = math.inf
    best_box = (a, b)
    best_enclosure = (0.0, 0.0)
    best_outside = math.inf
    widest = 0.0

    while stack:
        lo, hi, depth = stack.pop()
        examined += 1
        deepest = max(deepest, depth)
        vlo, vhi = evaluate(lo, hi)
        widest = max(widest, vhi - vlo)

        if vlo > target:
            accepted += 1
            if vlo < best_lower:
                best_lower, best_box, best_enclosure = vlo, (lo, hi), (vlo, vhi)
            # A box counts toward the outside bound unless it lies *wholly*
            # inside the excluded interval. A box straddling the edge contains
            # territory the other warrant does not govern, so its bound has to
            # be honoured here -- dropping it would leave that sliver unbounded.
            if exclude is not None and not (exclude[0] <= lo and hi <= exclude[1]):
                best_outside = min(best_outside, vlo)
            continue
        if depth >= max_depth:
            raise NotSeparated(
                f"{quantity} not separated from {target} on [{lo}, {hi}] after "
                f"{max_depth} splits; enclosure [{vlo}, {vhi}] — ENG-005 §15 stop "
                "condition"
            )
        mid = (lo + hi) / 2
        stack.append((lo, mid, depth + 1))
        stack.append((mid, hi, depth + 1))

    return CoverResult(
        quantity=quantity,
        certified_lower_bound=best_lower,
        target=target,
        boxes_examined=examined,
        boxes_accepted=accepted,
        max_depth=deepest,
        min_box=best_box,
        min_box_enclosure=best_enclosure,
        max_width=widest,
        initial_boxes=initial_boxes,
        cell=(a, b),
        lower_bound_outside=best_outside,
        exclude=exclude,
    )
