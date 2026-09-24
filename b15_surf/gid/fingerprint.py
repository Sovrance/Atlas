"""Fingerprint features from the label-stripped dataset only.

Features (all exact):
  pole_set          chords X_{i,j} in which the amplitude has a simple pole somewhere
                    on the registered probe line (single-variable model fit through
                    ``pir.symbolic.linear.solve``: polynomial model first, then
                    "alpha/(t+d) + poly" with the pole location d as a LINEAR unknown;
                    an overdetermined UNIQUE fit is the exactness certificate,
                    INCONSISTENT for both models -> NONIDENTIFIABLE(model order))
  pole_locations    {chord: location} — representation-specific (carries delta)
  zero_loci         registered loci on which every record value is exactly 0
  split_structure   exact rank of the m x m value matrix on each registered split
                    grid (RANK_TEST = b3_electroweak.rank.matrix_rank); rank 1 is the
                    factorization signature of eq. (1.7) of arXiv:2405.09608
  redundancy_metric (#triangulation terms) / (#independent fingerprint constraints) (P1)
Canonicalisation (R15 rule as code): every chord-labelled feature is reduced to its
lexicographically minimal image under the dihedral relabelling group D_n
(cyclic + reflection of the colour order) BEFORE hashing; ``canon.engine``'s kinds
(generator / choi / covariance / metric_pair, numpy-based) do not apply to
amplitude fingerprints, so the orbit-minimum canonicalisation is implemented
here exactly (see docs/sprint-B15-SURF-reconciliation.md).
"""

from __future__ import annotations

from fractions import Fraction
from math import comb
from typing import Dict, List, Optional, Tuple

from b3_electroweak.rank import matrix_rank
from pir.symbolic.linear import solve

from ..exact import fmt

POLY_DEGREE = 3          # registered maximal polynomial degree of the single-variable model


# --------------------------------------------------------------------------- #
# dihedral canonicalisation                                                    #
# --------------------------------------------------------------------------- #
def _norm(n: int, i: int, j: int) -> Tuple[int, int]:
    i = (i - 1) % n + 1
    j = (j - 1) % n + 1
    return (i, j) if i < j else (j, i)


def dihedral_maps(n: int):
    for k in range(n):
        yield lambda v, k=k: (v - 1 + k) % n + 1                 # rotation
        yield lambda v, k=k: (-(v - 1) + k) % n + 1              # reflection
    return


def canonical_chord_set(n: int, chords: List[Tuple[int, int]]) -> Tuple[Tuple[int, int], ...]:
    best = None
    for g in dihedral_maps(n):
        img = tuple(sorted(_norm(n, g(i), g(j)) for (i, j) in chords))
        if best is None or img < best:
            best = img
    return best


def canonical_labelled(n: int, items: Dict[Tuple[int, int], str]) -> Tuple[Tuple[int, int, str], ...]:
    best = None
    for g in dihedral_maps(n):
        img = tuple(sorted((*_norm(n, g(i), g(j)), v) for (i, j), v in items.items()))
        if best is None or img < best:
            best = img
    return best


def canonical_pair_sets(n: int, sets: List[List[Tuple[int, int]]]) -> Tuple[Tuple[Tuple[int, int], ...], ...]:
    best = None
    for g in dihedral_maps(n):
        img = tuple(sorted(tuple(sorted(_norm(n, g(i), g(j)) for (i, j) in s)) for s in sets))
        if best is None or img < best:
            best = img
    return best


# --------------------------------------------------------------------------- #
# single-variable pole model                                                   #
# --------------------------------------------------------------------------- #
def fit_single_variable(samples: List[Tuple[Fraction, Fraction]]) -> Dict:
    """samples = [(t, A(t))]. Returns {model: 'polynomial'|'pole', location: d or None}."""
    ts = [t for t, _ in samples]
    vs = [v for _, v in samples]
    # (i) polynomial model  A(t) = sum_{m<=D} g_m t^m
    A = [[t ** m for m in range(POLY_DEGREE + 1)] for t in ts]
    res = solve(A, vs)
    if res["status"] == "UNIQUE":
        return {"model": "polynomial", "location": None, "rank": res["rank"]}
    if res["status"] == "UNDERDETERMINED":
        return {"model": "underdetermined", "location": None, "rank": res["rank"]}
    # (ii) pole model  (t + d) A(t) = sum_{m<=D+1} g_m t^m  ->  A(t) d - sum g_m t^m = -t A(t)
    A2 = [[v] + [-(t ** m) for m in range(POLY_DEGREE + 2)] for t, v in samples]
    b2 = [-t * v for t, v in samples]
    res2 = solve(A2, b2)
    if res2["status"] == "UNIQUE":
        d = res2["solution"][0]
        return {"model": "pole", "location": -d, "rank": res2["rank"]}
    return {"model": "unresolved:" + res2["status"], "location": None, "rank": res2.get("rank")}


# --------------------------------------------------------------------------- #
# feature extraction                                                           #
# --------------------------------------------------------------------------- #
def _chord(s: str) -> Tuple[int, int]:
    i, j = s.split(",")
    return (int(i), int(j))


def features_for(ds: Dict, cls: str, n: int) -> Dict:
    recs = [r for r in ds["records"] if r["class"] == cls and r["n"] == n]
    if not recs:
        raise KeyError((cls, n))
    warnings: List[Dict[str, str]] = []

    # poles
    probes: Dict[Tuple[int, int], List[Tuple[Fraction, Fraction]]] = {}
    for r in recs:
        if r["kind"] == "pole_probe":
            probes.setdefault(_chord(r["probe"]["chord"]), []).append(
                (Fraction(r["probe"]["t"]), Fraction(r["value"])))
    pole_locations: Dict[Tuple[int, int], Fraction] = {}
    unresolved = []
    for ch, samples in sorted(probes.items()):
        fit = fit_single_variable(samples)
        if fit["model"] == "pole":
            pole_locations[ch] = fit["location"]
        elif fit["model"] != "polynomial":
            unresolved.append((ch, fit["model"]))
    if unresolved:
        warnings.append({"location": "b15_surf/gid/fingerprint.py::fit_single_variable",
                         "text": f"unresolved single-variable models at n={n}: {unresolved}"})
    pole_set = sorted(pole_locations)

    # zeros
    loci = ds["registered"]["loci"][str(n)]
    vanishing, nonvanishing_controls = [], []
    for lid, pairs in loci.items():
        vals = [Fraction(r["value"]) for r in recs if r["kind"] == "locus" and r["locus_id"] == lid]
        if vals and all(v == 0 for v in vals):
            vanishing.append([_chord(p) for p in pairs])
    ctrl_vals = [Fraction(r["value"]) for r in recs if r["kind"] == "control"]
    nonzero_rate = (Fraction(sum(1 for v in ctrl_vals if v != 0), len(ctrl_vals))
                    if ctrl_vals else None)

    # splits (exact rank of the grid matrix)
    split_ranks: Dict[str, int] = {}
    tri: Dict[str, Tuple[int, int, int]] = {}
    for r in recs:
        if r["kind"] in ("split_grid", "split_ctrl"):
            tri[r["split_id"]] = tuple(r["triangle"])
    for sid, t in tri.items():
        grid = [r for r in recs if r["kind"] in ("split_grid", "split_ctrl") and r["split_id"] == sid]
        m = max(g["grid"][0] for g in grid) + 1
        M = [[Fraction(0)] * m for _ in range(m)]
        for g in grid:
            M[g["grid"][0]][g["grid"][1]] = Fraction(g["value"])
        split_ranks[sid] = matrix_rank(M)
    split_structure = tuple(sorted((canonical_chord_set(n, [(t[0], t[1]), (t[1], t[2]), (t[0], t[2])]),
                                    split_ranks[sid])
                                   for sid, t in tri.items() if not sid.endswith("-ctrl")))
    split_control_ranks = {sid: split_ranks[sid] for sid in tri if sid.endswith("-ctrl")}

    # delta from pole locations (odd/odd and even/even chords)
    deltas = []
    for (i, j), loc in pole_locations.items():
        if i % 2 == 0 and j % 2 == 0:
            deltas.append(-loc)         # X_ee + delta = 0  ->  loc = -delta
        elif i % 2 == 1 and j % 2 == 1:
            deltas.append(loc)          # X_oo - delta = 0  ->  loc = +delta
    if not deltas:
        delta = {"verdict": "NONIDENTIFIABLE", "cause": "delta-underdetermined",
                 "detail": "no X_{e,e} / X_{o,o} pole present"}
    elif all(d == deltas[0] for d in deltas):
        delta = {"verdict": "FORCED", "value": fmt(deltas[0]), "n_constraints": len(deltas)}
    else:
        delta = {"verdict": "NONIDENTIFIABLE", "cause": "delta-inconsistent",
                 "detail": sorted({fmt(d) for d in deltas})}

    # delta-normalised location pattern (representation-specific feature)
    dval = Fraction(delta["value"]) if delta["verdict"] == "FORCED" else None
    pattern = {}
    for ch, loc in pole_locations.items():
        pattern[ch] = fmt(loc / dval) if (dval not in (None, 0)) else fmt(loc)

    n_atoms = len(pole_set) + len(vanishing) + len(split_structure)
    redundancy = Fraction(comb(2 * (n - 2), n - 2) // (n - 1), n_atoms) if n_atoms else None

    return {
        "n": n,
        "pole_set": canonical_chord_set(n, pole_set),
        "zero_loci": canonical_pair_sets(n, vanishing),
        "split_structure": split_structure,
        "split_control_ranks": split_control_ranks,
        "pole_location_pattern": canonical_labelled(n, pattern),
        "delta": delta,
        "control_nonzero_rate": fmt(nonzero_rate) if nonzero_rate is not None else None,
        "n_vanishing_loci": len(vanishing),
        "n_registered_loci": len(loci),
        "redundancy_metric": fmt(redundancy) if redundancy is not None else None,
        "raw_pole_set": [f"{i},{j}" for (i, j) in pole_set],
        "warnings": warnings,
    }


def delta_class(feat: Dict) -> str:
    d = feat.get("delta") or {}
    if d.get("verdict") == "FORCED":
        return "zero" if Fraction(d["value"]) == 0 else "nonzero"
    return "absent"


def atoms(feat: Dict, keys=("pole_set", "zero_loci", "split_structure", "delta_class")) -> set:
    out = set()
    if "delta_class" in keys:
        out.add(("delta_class", feat.get("delta_class") or delta_class(feat)))
    if "pole_set" in keys:
        out |= {("pole",) + ch for ch in feat["pole_set"]}
    if "zero_loci" in keys:
        out |= {("zero", s) for s in feat["zero_loci"]}
    if "split_structure" in keys:
        out |= {("split",) + s for s in feat["split_structure"]}
    return out
