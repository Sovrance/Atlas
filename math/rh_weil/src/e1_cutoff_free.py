"""Cutoff-free degree-1 and compact degree-2 E1 (ATLAS-RH-ENG-005 §4/§5).

Certifies uniform positive lower bounds on the closed cell ``[log 3, log 4]`` for

* the odd pivot ``O1(L) = G[q1,q1]`` (degree 1, §4), and
* the compact even determinant ``E2(L) = G00 Gbb - G0b^2`` (degree 2, §5),

both assembled from the adjudicated Candidate-A pole and the **exact** real-space
archimedean term, so there is no frequency cutoff to bound away.

Method
------
Adaptive interval branch and bound. Each box uses the assembled mean-value form
from ``archimedean_realspace``: an exact point evaluation at the midpoint plus
the radius times an enclosure of the exact assembled derivative. Enclosure width
falls quadratically in the box radius (the derivative enclosure's own looseness
is what is left, and that is itself linear in the radius), so a box that fails
splits and its halves almost always succeed.

This is deliberately *not* the convexity argument ENG-004 used for the scalar
cell. Convexity there came from a closed form for ``G00''`` that has no analogue
for ``E2``, which is a product difference. Subdivision needs no such structure —
it only needs enclosures that shrink, which the real-space form provides and the
oscillatory frequency form does not.

The certificate records the topology actually established: a positive lower bound
by exhaustive cover, with the box count, the deepest split, and where the minimum
sits (ENG-005 §8's instruction not to precommit to a topology applies to the
T=84 uniform certificate, and the same honesty applies here).

No RH proof claim is made by this module.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

import archimedean_realspace as AR
import weil_entries as WE
from interval_backend import backend_info, interval_box, require_flint, set_precision_bits

CELL_LABEL = ("log(3)", "log(4)")
CLAIM_SCOPE = "finite_dimensional_weil_compression"

DEFAULT_PRECISION_BITS = 220
#: Start coarse; the search refines where it must.
DEFAULT_INITIAL_BOXES = 24
DEFAULT_MAX_DEPTH = 22


@dataclass
class CoverResult:
    quantity: str
    certified_lower_bound: float
    boxes_examined: int
    boxes_accepted: int
    max_depth: int
    min_box: Tuple[float, float]
    min_box_enclosure: Tuple[float, float]
    max_width: float
    initial_boxes: int
    precision_bits: int
    detail: List[Dict[str, Any]] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "quantity": self.quantity,
            "certified_lower_bound": repr(self.certified_lower_bound),
            "boxes_examined": self.boxes_examined,
            "boxes_accepted": self.boxes_accepted,
            "max_subdivision_depth": self.max_depth,
            "argmin_box": [repr(self.min_box[0]), repr(self.min_box[1])],
            "argmin_box_enclosure": [repr(self.min_box_enclosure[0]),
                                     repr(self.min_box_enclosure[1])],
            "max_box_enclosure_width": repr(self.max_width),
            "initial_boxes": self.initial_boxes,
            "precision_bits": self.precision_bits,
            "method": "adaptive_interval_branch_and_bound",
            "topology_proved": "uniform positive lower bound by exhaustive cover",
        }


def certify_positive(
    quantity: str,
    *,
    precision_bits: int = DEFAULT_PRECISION_BITS,
    initial_boxes: int = DEFAULT_INITIAL_BOXES,
    max_depth: int = DEFAULT_MAX_DEPTH,
    target: float = 0.0,
    cell: Optional[Tuple[float, float]] = None,
) -> CoverResult:
    """Certify ``quantity > target`` uniformly on the cell by exhaustive cover.

    ``quantity`` is a key of ``archimedean_realspace.block_centred``:
    ``"G00"``, ``"G0b"``, ``"Gbb"``, ``"O1"``, ``"E2"`` or ``"det_deg2"``.
    """
    _, arb, acb, _ = require_flint()
    set_precision_bits(precision_bits)

    a, b = cell if cell else (math.log(3.0), math.log(4.0))
    primes = WE.prime_powers_below((a + b) / 2)

    # (lo, hi, depth)
    stack: List[Tuple[float, float, int]] = [
        (a + (b - a) * k / initial_boxes, a + (b - a) * (k + 1) / initial_boxes, 0)
        for k in range(initial_boxes)
    ]
    examined = accepted = 0
    deepest = 0
    best_lower = math.inf
    best_box = (a, b)
    best_enclosure = (0.0, 0.0)
    widest = 0.0

    while stack:
        lo, hi, depth = stack.pop()
        examined += 1
        deepest = max(deepest, depth)
        mid, rad = (lo + hi) / 2, (hi - lo) / 2
        box = interval_box(lo, hi)
        blk = AR.block_centred(box, arb, acb, prime_powers=primes)
        val = blk[quantity]
        vlo, vhi = float(val.lower()), float(val.upper())
        widest = max(widest, vhi - vlo)

        if vlo > target:
            accepted += 1
            if vlo < best_lower:
                best_lower, best_box, best_enclosure = vlo, (lo, hi), (vlo, vhi)
            continue
        if depth >= max_depth:
            raise ValueError(
                f"{quantity} not separated from {target} on [{lo}, {hi}] after "
                f"{max_depth} splits; enclosure [{vlo}, {vhi}] — ENG-005 §15 stop "
                "condition"
            )
        stack.append((lo, mid, depth + 1))
        stack.append((mid, hi, depth + 1))

    return CoverResult(
        quantity=quantity,
        certified_lower_bound=best_lower,
        boxes_examined=examined,
        boxes_accepted=accepted,
        max_depth=deepest,
        min_box=best_box,
        min_box_enclosure=best_enclosure,
        max_width=widest,
        initial_boxes=initial_boxes,
        precision_bits=precision_bits,
    )


def parity_identities(L_values=None, precision_bits: int = DEFAULT_PRECISION_BITS):
    """Check ``D2 = E2 + L^2 G00 O1`` and ``det(G_deg<=2) = O1 * E2`` (§5).

    Both are consequences of the pole and prime blocks being parity block
    diagonal: the degree-2 Gram splits into the odd pivot and the even block, so
    its determinant is their product, and ``D2`` is the raw (unfactored) form.
    """
    import core

    _, arb, acb, _ = require_flint()
    set_precision_bits(precision_bits)
    rows = []
    for L in (L_values or (math.log(3.0), 1.2, 1.2828, math.log(4.0))):
        L_a = arb(repr(L))
        blk = AR.block_centred(L_a, arb, acb)
        g00, o1, e2 = blk["G00"], blk["O1"], blk["E2"]
        d2 = core.degree2_raw_det(g00, o1, e2, L_a)
        det = core.degree2_full_det(o1, e2)
        rows.append({
            "L": repr(L),
            "E2": repr(float(e2)),
            "O1": repr(float(o1)),
            "D2_raw": repr(float(d2)),
            "D2_matches_E2_plus_L2_G00_O1": bool((d2 - (e2 + L_a**2 * g00 * o1)).contains(0)),
            "det_matches_O1_times_E2": bool((det - blk["det_deg2"]).contains(0)),
        })
    return rows
