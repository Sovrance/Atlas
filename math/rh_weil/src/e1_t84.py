"""Direct-Fourier T=84 point and uniform E1 (ATLAS-RH-ENG-005 §7/§8/§10).

The T=84 objects are the truncated finite Weil matrix — a different object from
the cutoff-free entries of §4/§5, not an approximation of them. Their archimedean
term stops at ``T = 84`` by definition.

Point certificates (§7) are straight rigorous evaluations at chosen ``L``. The
uniform certificate (§10) covers the closed cell with the same adaptive
branch-and-bound every other uniform bound in the program uses, via the assembled
mean-value form: an exact point evaluation at each box midpoint plus the radius
times an enclosure of the **exact** first jet (``t84.entry_jet`` order 1, §9 —
no finite differences).

Choosing the topology (§8)
--------------------------
§8 says not to precommit to convexity or monotonicity, and to describe the
topology actually proved. :func:`describe_topology` reads the fresh scan and
classifies what it sees; the certificate then records both that classification
*and* the fact that the uniform bound itself was established by exhaustive cover,
which needs no topological assumption at all. If the scan and the cover ever
disagreed about where the minimum sits, the cover is the warrant and the scan is
a clue.

No RH proof claim is made by this module.
"""
from __future__ import annotations

import math
from typing import Any, Dict, List, Optional, Sequence, Tuple

import interval_cover as IC
import t84
import weil_entries as WE
from interval_backend import require_flint, set_precision_bits

CELL = (math.log(3.0), math.log(4.0))
CELL_LABEL = ("log(3)", "log(4)")
CLAIM_SCOPE = "finite_dimensional_weil_compression"
DEFAULT_PRECISION_BITS = t84.DEFAULT_PRECISION_BITS


# --------------------------------------------------------------------------- #
# §7 point certificates                                                        #
# --------------------------------------------------------------------------- #
def point_rows(points: Sequence[Tuple[str, float]], arb, acb, *,
               T: float = t84.T84, prime_powers=None) -> List[Dict[str, Any]]:
    """Rigorous Arb balls for the true finite Weil matrix at each point."""
    if prime_powers is None:
        prime_powers = WE.prime_powers_below(sum(CELL) / 2)
    rows = []
    for label, L in points:
        blk = t84.block_t84(arb(repr(L)), arb, acb, order=1, T=T,
                            prime_powers=prime_powers)
        rows.append({
            "label": label,
            "L": repr(L),
            "G00": [repr(float(blk["G00"].lower())), repr(float(blk["G00"].upper()))],
            "G0b": [repr(float(blk["G0b"].lower())), repr(float(blk["G0b"].upper()))],
            "Gbb": [repr(float(blk["Gbb"].lower())), repr(float(blk["Gbb"].upper()))],
            "O1": [repr(float(blk["O1"].lower())), repr(float(blk["O1"].upper()))],
            "E2": [repr(float(blk["E2"].lower())), repr(float(blk["E2"].upper()))],
            "E2_d1": repr(float(blk["E2_d1"])),
            "E2_definitely_positive": float(blk["E2"].lower()) > 0,
            "O1_definitely_positive": float(blk["O1"].lower()) > 0,
        })
    return rows


def selected_points(scan: Optional[Dict[str, Any]] = None) -> List[Tuple[str, float]]:
    """``log3``, the apparent bottleneck if the scan found one, ``1.20``, ``log4`` (§7)."""
    pts = [("log3", CELL[0])]
    if scan:
        bottleneck = scan.get("E2_min_on_grid", {})
        if bottleneck and not bottleneck.get("at_endpoint", True):
            pts.append(("apparent_bottleneck", float(bottleneck["L"])))
    pts.append(("1.20", 1.20))
    pts.append(("log4", CELL[1]))
    # keep them sorted and unique by value
    seen, out = set(), []
    for label, L in sorted(pts, key=lambda kv: kv[1]):
        key = round(L, 12)
        if key not in seen:
            seen.add(key)
            out.append((label, L))
    return out


# --------------------------------------------------------------------------- #
# §8 topology description                                                      #
# --------------------------------------------------------------------------- #
def describe_topology(scan: Dict[str, Any]) -> Dict[str, Any]:
    """Classify what the fresh Candidate-A scan shows, without assuming it."""
    stationary = scan.get("E2_stationary_points", [])
    curvature = scan.get("E2_curvature_changes", [])
    at_end = scan.get("E2_min_on_grid", {}).get("at_endpoint", None)
    if len(stationary) == 0:
        kind = "monotone_on_grid"
        why = ("E2' shows no sign change on the scan grid, so the grid minimum sits "
               "at an endpoint")
    elif len(stationary) == 1:
        kind = "single_interior_minimum_on_grid"
        why = "E2' changes sign exactly once on the scan grid"
    else:
        kind = "multiple_critical_points_on_grid"
        why = f"E2' changes sign {len(stationary)} times on the scan grid"
    return {
        "classification": kind,
        "reason": why,
        "stationary_points": stationary,
        "curvature_changes": curvature,
        "grid_min_at_endpoint": at_end,
        "warrant": (
            "The uniform bound is established by exhaustive interval cover, which "
            "assumes no topology. This classification describes what the fresh "
            "Candidate-A scan shows and is E3 evidence, not the warrant."
        ),
        "candidate_b_topology_reused": False,
    }


# --------------------------------------------------------------------------- #
# §10 uniform certificate                                                      #
# --------------------------------------------------------------------------- #
def _entry_evaluator(quantity: str, arb, acb, prime_powers, *, T: float,
                     options=None):
    """Box evaluator using the assembled mean-value form with exact jets."""

    def evaluate(lo: float, hi: float) -> Tuple[float, float]:
        mid, rad = (lo + hi) / 2, (hi - lo) / 2
        centre = t84.block_t84(arb(repr(mid)), arb, acb, order=0, T=T,
                               prime_powers=prime_powers, options=options,
                               keys=t84.EVEN_KEYS)[quantity]
        if rad == 0.0:
            return float(centre.lower()), float(centre.upper())
        box = arb(repr(mid), repr(rad))
        d1 = t84.block_t84(box, arb, acb, order=1, T=T, prime_powers=prime_powers,
                           options=options, keys=t84.EVEN_KEYS)[f"{quantity}_d1"]
        slope = max(abs(float(d1.lower())), abs(float(d1.upper())))
        enc = centre + arb(0, rad * slope)
        return float(enc.lower()), float(enc.upper())

    return evaluate


#: Working precision and integrator tolerance for the uniform cover.
#:
#: A box's enclosure width is set by the interval-L dependency, not by the
#: quadrature: at halfwidth 4e-4 the width is 4.175e-6 at rel_tol 1e-25 and
#: 4.174e-6 at rel_tol 1e-10 -- indistinguishable -- while the cost falls 6x
#: (29.5s to 4.9s per box). Spending 140 bits and 1e-25 here buys accuracy
#: nothing consumes; the residual quadrature contribution is ~1e-21 against a
#: 4e-6 width.
UNIFORM_PRECISION_BITS = 110
UNIFORM_INTEGRAL_OPTIONS = {"rel_tol": 1e-11}

#: E2_84 bottoms out at ~3.46e-6 at the left endpoint (fresh scan), and the box
#: width grows like ~48 r^2, so separation there needs r <~ 1.8e-4. Starting the
#: cover near that spacing means most boxes clear at depth 0 instead of being
#: split four or five times from a coarse start -- the same total work, minus the
#: discarded parents.
UNIFORM_INITIAL_BOXES = 448


def certify_uniform(
    quantity: str = "E2",
    *,
    T: float = t84.T84,
    precision_bits: int = UNIFORM_PRECISION_BITS,
    initial_boxes: int = UNIFORM_INITIAL_BOXES,
    max_depth: int = 8,
    target: float = 0.0,
    options=None,
) -> IC.CoverResult:
    """Uniform positive lower bound for a T=84 quantity over the closed cell."""
    if options is None:
        options = UNIFORM_INTEGRAL_OPTIONS
    _, arb, acb, _ = require_flint()
    set_precision_bits(precision_bits)
    primes = WE.prime_powers_below(sum(CELL) / 2)
    result = IC.adaptive_cover(
        _entry_evaluator(quantity, arb, acb, primes, T=T, options=options),
        quantity=f"{quantity}_T{int(T)}",
        cell=CELL,
        target=target,
        initial_boxes=initial_boxes,
        max_depth=max_depth,
    )
    result.extra.update({
        "T": T,
        "precision_bits": precision_bits,
        "box_form": "assembled mean-value: exact midpoint value + radius * exact-jet enclosure",
        "jets": "exact support-length jets (ENG-005 §9); no finite differences",
        "quadrature": "deterministic Arb panel schedule for T=84",
    })
    return result
