"""Ground-truth amplitude dataset writer — data/b15/zero_dataset_v1.json.

Every value is an exact ``"p/q"`` string (pir/canonical convention); the file
is hash-stable canonical JSON. Record kinds (all registered in prereg-002 §ZERO):

  generic     seeded generic point                      (Route A == Route B checked in tests)
  locus       point on a registered hidden-zero locus   (value must be 0)
  control     point on a random linear constraint set   (negative control; generically nonzero)
  pole_probe  one chord varied over registered t-values (single-variable pole/regularity probe)
  split_grid  m x m grid of split kinematics            (rank-1 factorization probe, eq. 1.7)
  split_ctrl  m x m grid of NON-split kinematics        (negative control for split_grid)

Theories (classes): TrPhi3 (n = 4,5,6,8), NLSM (n = 4,6,8, via the delta -> inf limit,
eq. 2 of 2401.05483) and TrPhi3_shifted (finite registered delta, n = 4,6,8) — the
third class exists so that delta-recovery (WP3.4) has a positive instance; the
NLSM instance is registered to return NONIDENTIFIABLE(delta-underdetermined).
"""

from __future__ import annotations

import hashlib
import os
import random
import sys
from fractions import Fraction
from typing import Dict, List

from pir.canonical import canonical_json, sha256_hex

from ..exact import fmt
from .delta_shift import nlsm_from_shift, shifted_point
from .kinematics import Chord, Point, chords, random_point, random_rational
from .recursion import amplitude_B
from .splits import (_sub_chords, random_split_kinematics, registered_splits,
                     sub_point, surfaces)
from .zeros import (random_constraints, registered_loci, sample_on_constraints,
                    sample_on_locus)

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
DATASET_PATH = os.path.join(ROOT, "data", "b15", "zero_dataset_v1.json")

# ---- registered generator parameters (prereg-002 §ZERO) -------------------- #
SEED = 20260923
DELTA_SHIFTED = Fraction(7, 3)          # hidden from the GID pass; registered here
N_GENERIC = {"TrPhi3": {4: 60, 5: 200, 6: 200, 8: 200},
             "NLSM": {4: 60, 6: 200, 8: 200},
             "TrPhi3_shifted": {4: 60, 6: 200, 8: 200}}
N_LOCUS_POINTS = 3                       # per locus per theory
N_CONTROL_SETS = 3                       # random constraint systems per n
N_CONTROL_POINTS = 15                    # per control set
PROBE_T = [Fraction(v) for v in ("1/2", "1/3", "1/5", "1/7", "2/3", "3/2", "5/2", "4/3", "-1/2", "-3/4")]
SPLIT_GRID_M = 3
NONZERO_RATE_THRESHOLD = Fraction(95, 100)


def theory_value(theory: str, pt: Point, n: int) -> Fraction:
    if theory == "TrPhi3":
        return amplitude_B(pt, n)
    if theory == "NLSM":
        lim, low = nlsm_from_shift(pt, n)
        assert all(x == 0 for x in low)
        return lim
    if theory == "TrPhi3_shifted":
        return amplitude_B(shifted_point(pt, n, DELTA_SHIFTED), n)
    raise ValueError(theory)


def theory_ns(theory: str) -> List[int]:
    return sorted(N_GENERIC[theory])


def _pt_json(pt: Point) -> Dict[str, str]:
    return {f"{i},{j}": fmt(v) for (i, j), v in sorted(pt.items())}


def generator_sha256() -> str:
    h = hashlib.sha256()
    for name in sorted(os.listdir(HERE)):
        if name.endswith(".py"):
            with open(os.path.join(HERE, name), "rb") as f:
                h.update(name.encode() + b"\0" + f.read() + b"\0")
    return h.hexdigest()


def _pole_free(theory: str, pt: Point, n: int) -> bool:
    q = shifted_point(pt, n, DELTA_SHIFTED) if theory == "TrPhi3_shifted" else pt
    return all(v != 0 for v in q.values())


def _safe_probe_point(theory: str, n: int, rng: random.Random) -> Point:
    """Base point whose probe values never hit a (possibly shifted) pole."""
    while True:
        pt = random_point(n, rng)
        ok = True
        for ch in chords(n):
            for t in PROBE_T:
                q = dict(pt)
                q[ch] = t
                if theory == "TrPhi3_shifted":
                    q = shifted_point(q, n, DELTA_SHIFTED)
                if any(v == 0 for v in q.values()):
                    ok = False
        if ok:
            return pt


def build(seed: int = SEED) -> Dict:
    rng = random.Random(seed)
    records: List[Dict] = []
    loci_by_n = {n: registered_loci(n) for n in (4, 5, 6, 8)}
    splits_by_n = {n: registered_splits(n) for n in (5, 6, 8)}

    for theory in ("TrPhi3", "NLSM", "TrPhi3_shifted"):
        for n in theory_ns(theory):
            # generic points
            for _ in range(N_GENERIC[theory][n]):
                pt = random_point(n, rng)
                while not _pole_free(theory, pt, n):
                    pt = random_point(n, rng)
                records.append({"theory": theory, "n": n, "kind": "generic",
                                "point": _pt_json(pt), "value": fmt(theory_value(theory, pt, n)),
                                "on_zero_locus": False, "split_id": None, "locus_id": None})
            # registered zero loci
            for lid, pairs in loci_by_n[n].items():
                for _ in range(N_LOCUS_POINTS):
                    pt = sample_on_locus(n, pairs, rng)
                    while not _pole_free(theory, pt, n):
                        pt = sample_on_locus(n, pairs, rng)
                    records.append({"theory": theory, "n": n, "kind": "locus",
                                    "point": _pt_json(pt), "value": fmt(theory_value(theory, pt, n)),
                                    "on_zero_locus": True, "split_id": None, "locus_id": lid,
                                    "constraints": [f"c{a},{b}=0" for (a, b) in pairs]})
            # negative-control constraint sets (same count as a locus of the same n)
            m = max(len(p) for p in loci_by_n[n].values())
            for s in range(N_CONTROL_SETS):
                A = random_constraints(n, m, rng)
                for _ in range(N_CONTROL_POINTS):
                    pt = sample_on_constraints(n, A, rng)
                    while not _pole_free(theory, pt, n):
                        pt = sample_on_constraints(n, A, rng)
                    records.append({"theory": theory, "n": n, "kind": "control",
                                    "point": _pt_json(pt), "value": fmt(theory_value(theory, pt, n)),
                                    "on_zero_locus": False, "split_id": None,
                                    "locus_id": f"CTRL{n}-{s}"})
            # single-variable pole probes
            base = _safe_probe_point(theory, n, rng)
            for ch in chords(n):
                for t in PROBE_T:
                    q = dict(base)
                    q[ch] = t
                    records.append({"theory": theory, "n": n, "kind": "pole_probe",
                                    "point": _pt_json(q), "value": fmt(theory_value(theory, q, n)),
                                    "on_zero_locus": False, "split_id": None, "locus_id": None,
                                    "probe": {"chord": f"{ch[0]},{ch[1]}", "t": fmt(t)}})
            # split grids (+ non-split control grids)
            if n in splits_by_n:
                for sid, (i, j, k) in splits_by_n[n].items():
                    P, Q = surfaces(n, i, j, k)
                    from .splits import split_point
                    while True:      # resample until every grid point is pole-free
                        xs = [{c_: random_rational(rng) for c_ in _sub_chords(P)} for _ in range(SPLIT_GRID_M)]
                        ys = [{c_: random_rational(rng) for c_ in _sub_chords(Q)} for _ in range(SPLIT_GRID_M)]
                        if all(_pole_free(theory, split_point(n, i, j, k, xa, yb), n)
                               for xa in xs for yb in ys):
                            break
                    for a in range(SPLIT_GRID_M):
                        for b in range(SPLIT_GRID_M):
                            pt = split_point(n, i, j, k, xs[a], ys[b])
                            val = theory_value(theory, pt, n)
                            records.append({"theory": theory, "n": n, "kind": "split_grid",
                                            "point": _pt_json(pt), "value": fmt(val),
                                            "on_zero_locus": False, "split_id": sid, "locus_id": None,
                                            "grid": [a, b], "triangle": [i, j, k]})
                    # control: vary the P-only chords with a and the rest with b, no split relations
                    Pset = set(P)
                    while True:
                        base = random_point(n, rng)
                        varA = [dict((c_, random_rational(rng)) for c_ in chords(n) if c_[0] in Pset and c_[1] in Pset)
                                for _ in range(SPLIT_GRID_M)]
                        varB = [dict((c_, random_rational(rng)) for c_ in chords(n) if not (c_[0] in Pset and c_[1] in Pset))
                                for _ in range(SPLIT_GRID_M)]
                        if all(_pole_free(theory, {**base, **va, **vb}, n) for va in varA for vb in varB):
                            break
                    for a in range(SPLIT_GRID_M):
                        for b in range(SPLIT_GRID_M):
                            pt = dict(base)
                            pt.update(varA[a])
                            pt.update(varB[b])
                            records.append({"theory": theory, "n": n, "kind": "split_ctrl",
                                            "point": _pt_json(pt), "value": fmt(theory_value(theory, pt, n)),
                                            "on_zero_locus": False, "split_id": sid + "-ctrl", "locus_id": None,
                                            "grid": [a, b], "triangle": [i, j, k]})

    ds = {
        "schema": "b15-zero-dataset-v1",
        "seed": seed,
        "python_version": sys.version.split()[0],
        "generator_sha256": generator_sha256(),
        "theories": {t: {"n": theory_ns(t)} for t in N_GENERIC},
        "registered": {
            "delta_shifted": fmt(DELTA_SHIFTED),
            "probe_t": [fmt(t) for t in PROBE_T],
            "split_grid_m": SPLIT_GRID_M,
            "loci": {str(n): {lid: [f"{a},{b}" for (a, b) in pairs] for lid, pairs in L.items()}
                     for n, L in loci_by_n.items()},
            "splits": {str(n): {sid: list(t) for sid, t in S.items()} for n, S in splits_by_n.items()},
        },
        "records": records,
    }
    ds["dataset_sha256"] = sha256_hex(ds["records"])
    return ds


def write(path: str = DATASET_PATH, seed: int = SEED) -> Dict:
    ds = build(seed)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        f.write(canonical_json(ds))
        f.write("\n")
    return ds


def load(path: str = DATASET_PATH) -> Dict:
    import json
    with open(path) as f:
        return json.load(f)


def label_stripped(ds: Dict, class_map: Dict[str, str]) -> Dict:
    """The view the GID pass sees: theory names replaced by opaque class ids,
    registered delta removed."""
    out = {"schema": ds["schema"] + "-stripped", "records": []}
    for r in ds["records"]:
        r2 = dict(r)
        r2["class"] = class_map[r2.pop("theory")]
        out["records"].append(r2)
    reg = dict(ds["registered"])
    reg.pop("delta_shifted", None)
    out["registered"] = reg
    return out
