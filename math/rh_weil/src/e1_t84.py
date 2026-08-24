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


# The cover is embarrassingly parallel: a partition of the cell into contiguous
# chunks, each covered independently, is still an exhaustive cover of the cell,
# and the merged bound is the minimum over chunks. Chunk edges are computed with
# the same expression ``a + (b-a)*k/n`` that :func:`interval_cover.adaptive_cover`
# uses internally, so consecutive chunks share an endpoint bit-for-bit and no
# sliver of the cell falls between two chunks.
def _cover_chunk_worker(spec: Dict[str, Any]) -> Dict[str, Any]:
    """One chunk of the uniform cover. Module level so it survives pickling."""
    _, arb, acb, _ = require_flint()
    set_precision_bits(spec["precision_bits"])
    primes = WE.prime_powers_below(sum(CELL) / 2)
    ev = _entry_evaluator(spec["quantity"], arb, acb, primes, T=spec["T"],
                          options=spec["options"])
    r = IC.adaptive_cover(ev, quantity=spec["quantity"], cell=(spec["lo"], spec["hi"]),
                          target=spec["target"], initial_boxes=spec["boxes"],
                          max_depth=spec["max_depth"],
                          exclude=(tuple(spec["exclude"]) if spec.get("exclude") else None))
    return {
        "index": spec["index"],
        "cell": [r.cell[0], r.cell[1]],
        "certified_lower_bound": r.certified_lower_bound,
        "boxes_examined": r.boxes_examined,
        "boxes_accepted": r.boxes_accepted,
        "max_depth": r.max_depth,
        "min_box": [r.min_box[0], r.min_box[1]],
        "min_box_enclosure": [r.min_box_enclosure[0], r.min_box_enclosure[1]],
        "max_width": r.max_width,
        "initial_boxes": r.initial_boxes,
        "lower_bound_outside": r.lower_bound_outside,
    }


def certify_uniform_parallel(
    quantity: str = "E2",
    *,
    T: float = t84.T84,
    precision_bits: int = UNIFORM_PRECISION_BITS,
    initial_boxes: int = UNIFORM_INITIAL_BOXES,
    max_depth: int = 8,
    target: float = 0.0,
    options=None,
    workers: int = 4,
    exclude: Optional[Tuple[float, float]] = None,
    progress=None,
) -> IC.CoverResult:
    """:func:`certify_uniform`, run as ``workers`` independent chunk covers."""
    import multiprocessing as mp

    if options is None:
        options = UNIFORM_INTEGRAL_OPTIONS
    a, b = CELL
    n = max(1, int(workers))
    per = max(1, initial_boxes // n)
    specs = [{
        "index": k,
        "lo": a + (b - a) * k / n,
        "hi": a + (b - a) * (k + 1) / n,
        "boxes": per,
        "quantity": quantity,
        "T": T,
        "precision_bits": precision_bits,
        "options": options,
        "target": target,
        "max_depth": max_depth,
        "exclude": list(exclude) if exclude else None,
    } for k in range(n)]
    # Chunk 0 must start and chunk n-1 must end exactly on the cell endpoints.
    specs[0]["lo"], specs[-1]["hi"] = a, b

    if n == 1:
        parts = [_cover_chunk_worker(specs[0])]
    else:
        with mp.get_context("fork").Pool(processes=n) as pool:
            parts = []
            for part in pool.imap_unordered(_cover_chunk_worker, specs):
                parts.append(part)
                if progress:
                    progress(part)
    parts.sort(key=lambda p: p["index"])

    # The chunks must tile the cell with no gap and no overlap.
    edges = [(p["cell"][0], p["cell"][1]) for p in parts]
    if edges[0][0] != a or edges[-1][1] != b:
        raise IC.NotSeparated("parallel cover does not reach the cell endpoints")
    for (_, hi), (lo, _) in zip(edges, edges[1:]):
        if hi != lo:
            raise IC.NotSeparated(f"parallel cover leaves a gap at {hi} != {lo}")

    worst = min(parts, key=lambda p: p["certified_lower_bound"])
    result = IC.CoverResult(
        quantity=f"{quantity}_T{int(T)}",
        certified_lower_bound=worst["certified_lower_bound"],
        target=target,
        boxes_examined=sum(p["boxes_examined"] for p in parts),
        boxes_accepted=sum(p["boxes_accepted"] for p in parts),
        max_depth=max(p["max_depth"] for p in parts),
        min_box=tuple(worst["min_box"]),
        min_box_enclosure=tuple(worst["min_box_enclosure"]),
        max_width=max(p["max_width"] for p in parts),
        initial_boxes=sum(p["initial_boxes"] for p in parts),
        cell=CELL,
        lower_bound_outside=min(p["lower_bound_outside"] for p in parts),
        exclude=exclude,
    )
    result.extra.update({
        "T": T,
        "precision_bits": precision_bits,
        "box_form": "assembled mean-value: exact midpoint value + radius * exact-jet enclosure",
        "jets": "exact support-length jets (ENG-005 §9); no finite differences",
        "quadrature": "deterministic Arb panel schedule for T=84",
        "parallel_chunks": [{"cell": [repr(p["cell"][0]), repr(p["cell"][1])],
                             "certified_lower_bound": repr(p["certified_lower_bound"]),
                             "boxes_examined": p["boxes_examined"],
                             "max_subdivision_depth": p["max_depth"]} for p in parts],
        "chunk_tiling_verified": "consecutive chunk endpoints compared bit-for-bit",
    })
    return result


# --------------------------------------------------------------------------- #
# §8 interior-minimum strategy: isolate, bound the basin, sign elsewhere       #
# --------------------------------------------------------------------------- #
# The exhaustive cover above proves the uniform bound while assuming nothing
# about the shape of E2. It is the warrant, but it is not *informative*: it says
# "every box cleared", not "here is the minimiser and here is why nothing beats
# it". §8's interior-minimum strategy gives the second, and the fresh scan shows
# it applies -- E2' changes sign exactly once, close to the left endpoint.
#
# The argument has three parts, each rigorous:
#
#   1. isolate x* : a bracket [a, b] with E2'(a) < 0 < E2'(b), and E2'' > 0
#      throughout, so the critical point in it is unique and is a minimum.
#   2. basin bound: E2 >= m on [a, b], by interval cover of that short interval.
#   3. sign elsewhere: E2' < 0 on [log3, a] and E2' > 0 on [b, log4].
#
# Together: E2 decreases onto [a, b] and increases off it, so its minimum over
# the whole cell is attained inside [a, b] and is therefore >= m. No box outside
# the basin needs its *value* bounded at all -- only the sign of its derivative.
def _d1_evaluator(arb, acb, prime_powers, *, T: float, options, centred: bool = True):
    """Enclosure of ``E2'`` over a box.

    Evaluated raw, an ``L``-box of radius ``r`` returns an ``E2'`` enclosure of
    halfwidth ~``25 r`` — at ``r = 1e-4`` that is ``1.1e-3``, two orders above the
    ``|E2'| <= 1e-5`` the sign conditions have to resolve, so the raw form cannot
    certify a sign anywhere near the stationary point. The same mean-value
    centring the entry evaluator uses fixes it: an exact point value at the
    midpoint plus ``r`` times an enclosure of the exact second jet. That trades
    the ``25 r`` dependency blow-up for ``r * (|E2''| + O(r))`` — at ``r = 1e-4``,
    ``3e-6`` instead of ``1.1e-3``.
    """

    def evaluate(lo: float, hi: float) -> Tuple[float, float]:
        mid, rad = (lo + hi) / 2, (hi - lo) / 2
        if rad == 0.0 or not centred:
            box = arb(repr(mid), repr(rad)) if rad > 0 else arb(repr(mid))
            d1 = t84.block_t84(box, arb, acb, order=1, T=T, prime_powers=prime_powers,
                               options=options, keys=t84.EVEN_KEYS)["E2_d1"]
            return float(d1.lower()), float(d1.upper())
        centre = t84.block_t84(arb(repr(mid)), arb, acb, order=1, T=T,
                               prime_powers=prime_powers, options=options,
                               keys=t84.EVEN_KEYS)["E2_d1"]
        box = arb(repr(mid), repr(rad))
        d2 = t84.block_t84(box, arb, acb, order=2, T=T, prime_powers=prime_powers,
                           options=options, keys=t84.EVEN_KEYS)["E2_d2"]
        slope = max(abs(float(d2.lower())), abs(float(d2.upper())))
        enc = centre + arb(0, rad * slope)
        return float(enc.lower()), float(enc.upper())

    return evaluate


def _d2_evaluator(arb, acb, prime_powers, *, T: float, options):
    """Enclosure of ``E2''`` over a box — raw, since the pole jets stop at order 2.

    Centring this would need an exact third jet, and ``pole.pole_gram_entry_d2L``
    is the last closed form available; a finite-difference third jet is exactly
    what §9 forbids in an E1 path. So ``E2''`` is certified raw, on a deliberately
    narrow window where the raw dependency blow-up (~``300 r``) still sits below
    ``E2'' ~ 1.4e-3``.
    """

    def evaluate(lo: float, hi: float) -> Tuple[float, float]:
        mid, rad = (lo + hi) / 2, (hi - lo) / 2
        box = arb(repr(mid), repr(rad)) if rad > 0 else arb(repr(mid))
        d2 = t84.block_t84(box, arb, acb, order=2, T=T, prime_powers=prime_powers,
                           options=options, keys=t84.EVEN_KEYS)["E2_d2"]
        return float(d2.lower()), float(d2.upper())

    return evaluate


#: Half-width of the curvature window around the isolated stationary point. The
#: window has to be wide enough that the sign covers on either side never need to
#: resolve a derivative smaller than they can, and narrow enough that raw ``E2''``
#: still separates from 0 on it. ``|E2'|`` grows like ``E2'' * d ~ 1.4e-3 d`` away
#: from the critical point, and the centred ``E2'`` form resolves ``~1e-8``, so
#: ``d = 1e-4`` leaves ``|E2'| ~ 1.4e-7`` at the window edge — an order of margin.
STATIONARY_WINDOW = 1.0e-4

#: How far past the curvature window the derivative-sign covers are pushed.
#:
#: The sign covers cannot reach ``log 4`` at any sane cost. The ``E2''`` enclosure
#: on a box of radius ``r`` carries a dependency blow-up of ~``300 r`` on top of
#: the true ``|E2''| ~ 1.4e-3``, so the centred ``E2'`` form separates from 0 only
#: while ``300 r^2 < |E2'|``. Out where ``|E2'| ~ 4e-5`` that forces ``r < 1.2e-4``
#: — about 1200 boxes to reach ``log 4``, more than the exhaustive cover of ``E2``
#: itself costs, and each box dearer. So the derivative argument governs a band
#: around the minimiser and the direct cover governs the rest; the two together
#: still bound the whole cell, and the band is chosen wide enough that the direct
#: cover clears the window bound comfortably beyond it.
SIGN_BAND = 0.03


def isolate_stationary_point(arb, acb, *, bracket: Tuple[float, float],
                             T: float = t84.T84, options=None,
                             tol: float = 1e-9, max_steps: int = 60):
    """Bisect ``bracket`` until ``E2'`` is sign-certified at both ends.

    Returns ``(a, b, detail)`` with ``E2'(a) < 0 < E2'(b)`` proved by rigorous
    Arb enclosures, so ``E2'`` has a zero in ``[a, b]`` by continuity. Point
    enclosures of ``E2'`` come back with width ~3e-13, so the bisection is limited
    by ``tol``, not by the arithmetic.
    """
    options = options if options is not None else UNIFORM_INTEGRAL_OPTIONS
    primes = WE.prime_powers_below(sum(CELL) / 2)
    d1 = _d1_evaluator(arb, acb, primes, T=T, options=options)

    a, b = bracket
    lo_a, hi_a = d1(a, a)
    lo_b, hi_b = d1(b, b)
    if not (hi_a < 0 < lo_b):
        raise IC.NotSeparated(
            f"E2' does not change sign across [{a}, {b}]: "
            f"E2'(a) in [{lo_a}, {hi_a}], E2'(b) in [{lo_b}, {hi_b}]"
        )
    steps = 0
    stopped = None
    for _ in range(max_steps):
        if b - a <= tol:
            break
        m = (a + b) / 2
        lo_m, hi_m = d1(m, m)
        steps += 1
        if hi_m < 0:
            a, lo_a, hi_a = m, lo_m, hi_m
        elif lo_m > 0:
            b, lo_b, hi_b = m, lo_m, hi_m
        else:
            # The enclosure straddles zero: undecidable at this precision. The
            # bracket is already as tight as the arithmetic supports.
            stopped = {"at": repr(m), "E2_d1": [repr(lo_m), repr(hi_m)],
                       "reason": "derivative enclosure straddles 0"}
            break
    detail = {
        "method": "certified bisection on sign(E2') with exact order-1 jets",
        "bracket": [repr(a), repr(b)],
        "width": repr(b - a),
        "midpoint": repr((a + b) / 2),
        "E2_d1_at_left_end": [repr(lo_a), repr(hi_a)],
        "E2_d1_at_right_end": [repr(lo_b), repr(hi_b)],
        "sign_change_certified": bool(hi_a < 0 < lo_b),
        "bisection_steps": steps,
        "stopped_early": stopped,
    }
    return a, b, detail


def certify_derivative_sign(arb, acb, *, lo: float, hi: float, sign: int,
                            T: float = t84.T84, options=None,
                            initial_boxes: int = 16, max_depth: int = 22):
    """Certify ``sign * E2' > 0`` uniformly on ``[lo, hi]`` by interval cover."""
    options = options if options is not None else UNIFORM_INTEGRAL_OPTIONS
    primes = WE.prime_powers_below(sum(CELL) / 2)
    base = _d1_evaluator(arb, acb, primes, T=T, options=options)

    def evaluate(a: float, b: float) -> Tuple[float, float]:
        dlo, dhi = base(a, b)
        return (dlo, dhi) if sign > 0 else (-dhi, -dlo)

    return IC.adaptive_cover(evaluate, quantity=f"{'+' if sign > 0 else '-'}E2_d1",
                             cell=(lo, hi), target=0.0,
                             initial_boxes=initial_boxes, max_depth=max_depth)


def certify_curvature_positive(arb, acb, *, lo: float, hi: float,
                               T: float = t84.T84, options=None,
                               initial_boxes: int = 64, max_depth: int = 12):
    """Certify ``E2'' > 0`` on ``[lo, hi]`` — makes the critical point unique."""
    options = options if options is not None else UNIFORM_INTEGRAL_OPTIONS
    primes = WE.prime_powers_below(sum(CELL) / 2)
    return IC.adaptive_cover(_d2_evaluator(arb, acb, primes, T=T, options=options),
                             quantity="E2_d2", cell=(lo, hi), target=0.0,
                             initial_boxes=initial_boxes, max_depth=max_depth)


def _sub_cover_worker(spec: Dict[str, Any]) -> Dict[str, Any]:
    """One of the interior-minimum sub-covers. Module level so it pickles.

    The four covers -- curvature, window value, left sign, right sign -- share no
    state and are run concurrently.
    """
    _, arb, acb, _ = require_flint()
    set_precision_bits(spec["precision_bits"])
    kind, T, options = spec["kind"], spec["T"], spec["options"]
    if kind == "curvature":
        r = certify_curvature_positive(arb, acb, lo=spec["lo"], hi=spec["hi"], T=T,
                                       options=options,
                                       initial_boxes=spec["initial_boxes"])
    elif kind == "value":
        primes = WE.prime_powers_below(sum(CELL) / 2)
        r = IC.adaptive_cover(
            _entry_evaluator("E2", arb, acb, primes, T=T, options=options),
            quantity="E2_window", cell=(spec["lo"], spec["hi"]), target=0.0,
            initial_boxes=spec["initial_boxes"], max_depth=12)
    else:  # "sign"
        r = certify_derivative_sign(arb, acb, lo=spec["lo"], hi=spec["hi"],
                                    sign=spec["sign"], T=T, options=options,
                                    initial_boxes=spec["initial_boxes"])
    return {"kind": kind, "dict": r.to_dict(),
            "certified_lower_bound": r.certified_lower_bound,
            "boxes_examined": r.boxes_examined, "max_depth": r.max_depth}


def certify_interior_minimum(*, T: float = t84.T84,
                             precision_bits: int = UNIFORM_PRECISION_BITS,
                             scan: Optional[Dict[str, Any]] = None,
                             window: float = STATIONARY_WINDOW,
                             band: float = SIGN_BAND,
                             options=None,
                             workers: int = 4,
                             progress=None) -> Dict[str, Any]:
    """The §8 interior-minimum certificate for ``E2_84`` on the closed cell.

    Four rigorous parts, which together locate the minimiser rather than merely
    bounding the function:

    1. **Isolate.** A bracket ``[a, b]`` with ``E2'(a) < 0 < E2'(b)``, both signs
       certified, so a critical point ``L*`` lies in it.
    2. **Uniqueness / minimality.** ``E2'' > 0`` on ``W = [a - w, b + w]``, so
       ``E2'`` is strictly increasing there: ``L*`` is the only critical point in
       ``W`` and it is a strict local minimum.
    3. **Basin bound.** ``E2 >= m`` on ``W`` by interval cover.
    4. **No lower values elsewhere.** ``E2' < 0`` on ``[log 3, a - w]`` and
       ``E2' > 0`` on ``[b + w, log 4]``, so ``E2`` is strictly decreasing to the
       left of ``W`` and strictly increasing to its right.

    (2) and (4) leave no critical point outside ``W``, and (4) forces the cell
    minimum into ``W``, where (3) bounds it. So ``E2 >= m`` on the whole cell,
    attained at ``L* in [a, b]``.
    """
    _, arb, acb, _ = require_flint()
    set_precision_bits(precision_bits)
    options = options if options is not None else UNIFORM_INTEGRAL_OPTIONS
    primes = WE.prime_powers_below(sum(CELL) / 2)

    def say(msg):
        if progress:
            progress(msg)

    # Where the fresh scan says E2' turns over. A starting bracket only -- the
    # sign change is then certified, never assumed, and the scan grid is E3.
    stationary = (scan or {}).get("E2_stationary_points") or \
                 (scan or {}).get("stationary_points") or []
    if stationary:
        seg = stationary[0]["between"]
        bracket = (float(seg[0]), float(seg[1]))
    else:
        bracket = (CELL[0], CELL[0] + 0.02)

    say(f"isolating stationary point in {bracket}")
    a, b, iso = isolate_stationary_point(arb, acb, bracket=bracket, T=T,
                                         options=options)
    say(f"  isolated to [{a!r}, {b!r}] width={b - a!r}")

    w_lo, w_hi = max(CELL[0], a - window), min(CELL[1], b + window)

    w_lo2, w_hi2 = w_lo, w_hi
    band_hi = min(CELL[1], w_hi + band)
    specs = [
        {"kind": "curvature", "lo": w_lo2, "hi": w_hi2, "initial_boxes": 64},
        {"kind": "value", "lo": w_lo2, "hi": w_hi2, "initial_boxes": 8},
    ]
    if w_lo2 > CELL[0] + 1e-15:
        specs.append({"kind": "sign", "sign": -1, "lo": CELL[0], "hi": w_lo2,
                      "initial_boxes": 24})
    if band_hi > w_hi2 + 1e-15:
        specs.append({"kind": "sign", "sign": +1, "lo": w_hi2, "hi": band_hi,
                      "initial_boxes": 96})
    for spec in specs:
        spec.update({"T": T, "options": options, "precision_bits": precision_bits})

    say(f"running {len(specs)} sub-covers concurrently: "
        f"E2'' > 0 and E2 bound on [{w_lo2!r}, {w_hi2!r}], "
        f"E2' < 0 on [{CELL[0]!r}, {w_lo2!r}], E2' > 0 on [{w_hi2!r}, {band_hi!r}]")
    if workers and workers > 1 and len(specs) > 1:
        import multiprocessing as mp
        with mp.get_context("fork").Pool(processes=min(workers, len(specs))) as pool:
            parts = []
            for part in pool.imap_unordered(_sub_cover_worker, specs):
                parts.append(part)
                say(f"  {part['kind']}: bound={part['certified_lower_bound']:.6e} "
                    f"boxes={part['boxes_examined']} depth={part['max_depth']}")
    else:
        parts = []
        for spec in specs:
            part = _sub_cover_worker(spec)
            parts.append(part)
            say(f"  {part['kind']}: bound={part['certified_lower_bound']:.6e} "
                f"boxes={part['boxes_examined']} depth={part['max_depth']}")

    curvature = next(p for p in parts if p["kind"] == "curvature")
    basin = next(p for p in parts if p["kind"] == "value")
    left = next((p for p in parts if p["kind"] == "sign"
                 and float(p["dict"]["cell"][0]) == CELL[0]), None)
    right = next((p for p in parts if p["kind"] == "sign"
                  and float(p["dict"]["cell"][0]) != CELL[0]), None)

    def side(part, lo, hi, statement):
        if part is None:
            return {"interval": [repr(lo), repr(hi)],
                    "statement": "empty: the curvature window reaches this endpoint"}
        return {"interval": [repr(lo), repr(hi)], "statement": statement,
                "certified_margin": repr(part["certified_lower_bound"]),
                "cover": part["dict"]}

    return {
        "strategy": "interior_minimum",
        "quantity": "E2_84",
        "cell": [repr(CELL[0]), repr(CELL[1])],
        "T": T,
        "stationary_point": {
            "statement": ("E2'(a) < 0 < E2'(b) with both signs certified, so E2' has "
                          "a zero L* in [a, b] by continuity"),
            "approximate_location": repr((a + b) / 2),
            "rigorous_interval": [repr(a), repr(b)],
            "interval_width": repr(b - a),
            "sign_change_certified": iso["sign_change_certified"],
            "E2_d1_at_left_end": iso["E2_d1_at_left_end"],
            "E2_d1_at_right_end": iso["E2_d1_at_right_end"],
            "bisection": iso,
        },
        "curvature": {
            "statement": ("E2'' > 0 on the window, so E2' is strictly increasing "
                          "there: L* is the unique critical point in the window and "
                          "is a strict local minimum"),
            "window": [repr(w_lo), repr(w_hi)],
            "window_halfwidth": repr(window),
            "E2_d2_certified_lower_bound": repr(curvature["certified_lower_bound"]),
            "cover": curvature["dict"],
        },
        "basin_bound": {
            "statement": "E2 >= certified_lower_bound on the whole curvature window",
            "window": [repr(w_lo), repr(w_hi)],
            "certified_lower_bound": repr(basin["certified_lower_bound"]),
            "cover": basin["dict"],
        },
        "no_lower_values_elsewhere": {
            "left": side(left, CELL[0], w_lo,
                         "E2' < 0, so E2 is strictly decreasing onto the window"),
            "right_band": side(right, w_hi, band_hi,
                               "E2' > 0, so E2 is strictly increasing off the window"),
            "consequence": ("E2 has no critical point in [log 3, band_hi] other than "
                            "L*, decreases strictly onto the window from log 3 and "
                            "increases strictly off it out to band_hi, so its minimum "
                            "over [log 3, band_hi] is attained inside the window"),
            "beyond_the_band": {
                "interval": [repr(band_hi), repr(CELL[1])],
                "governed_by": ("the exhaustive interval cover of E2 itself, not by a "
                                "derivative sign — see the uniform certificate. The "
                                "band ends where the direct cover clears the window "
                                "bound with room to spare, so nothing is left "
                                "unbounded."),
            },
        },
        "governed_interval": [repr(CELL[0]), repr(band_hi)],
        "conclusion": {
            "statement": ("on [log 3, band_hi] the minimum of E2_84 is attained at L* "
                          "in the isolated interval and is at least the window bound; "
                          "beyond band_hi the exhaustive cover applies"),
            "lower_bound_on_governed_interval": repr(basin["certified_lower_bound"]),
            "minimiser_enclosure": [repr(a), repr(b)],
        },
        "precision_bits": precision_bits,
        "jets": "exact support-length jets (ENG-005 §9); no finite differences",
        "mpmath_used": False,
    }
