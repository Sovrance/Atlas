"""Blind identification: fingerprint -> KnownGrammarDB lookup -> candidate forest -> verdict.

* DB: ``pir.fingerprints.KnownGrammarDB`` populated from data/b15/known_grammars_v1.json
  (expected pole sets, zero loci, split ranks, deformation parameter and range,
  sources). Two-hash semantics: full = delta-insensitive invariants
  (``pole_set``, ``zero_loci``, ``split_structure`` — added to
  ``pir.fingerprints.INVARIANT_KEYS``); specific = also the delta-normalised
  ``pole_location_pattern``.
* Candidate forest: GVAR rules (``pir.candidates.apply_rules``) fired on the
  predicates read off the fingerprint; ``pir.candidates.evaluate`` gives the
  lattice verdict; ``lattice_fact`` lowers it to a PIR fact (re-homed to the
  ``analyst`` namespace / DOMAIN layer per work-order §6).
* Similarity / confidence (separate numbers, correlator named):
    correlator "b15.gid.jaccard_atoms_v1"  (primary; atoms = pole chords, vanishing
                loci, rank-1 splits, and the deformation-parameter CLASS
                zero / nonzero / absent read off the exact delta recovery)
      similarity  = |A_obs ∩ A_best| / |A_obs ∪ A_best|   (Jaccard over feature atoms)
      confidence  = fraction of the atoms that DISCRIMINATE best from runner-up on
                    which the observation agrees with best (1 when no runner-up)
    correlator "b15.gid.invariants_only_v1": delta-class atom removed (the
      delta-insensitive invariants alone — by construction cannot separate a base
      grammar from its deformation; registered to make that limit visible)
    correlator "b15.gid.zeros_only_v1": pole and delta atoms removed
      (registered secondary correlator: tests whether zeros + splits alone separate).
"""

from __future__ import annotations

import json
import os
from dataclasses import replace
from fractions import Fraction
from typing import Dict, List, Optional, Tuple

from pir.candidates import apply_rules, evaluate, lattice_fact
from pir.fingerprints import KnownGrammarDB, full_fingerprint, specific_fingerprint

from ..exact import fmt
from .fingerprint import (atoms, canonical_chord_set, canonical_labelled,
                          canonical_pair_sets)

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
KNOWN_GRAMMARS_PATH = os.path.join(ROOT, "data", "b15", "known_grammars_v1.json")

KNOWN_GRAMMARS = {
    "schema": "b15-known-grammars-v1",
    "families": {
        "TrPhi3": {
            "poles": "all_chords", "zeros": "all_registered_loci", "splits": "rank1_all_registered",
            "deformation": {"parameter": "delta", "value": "0/1", "range": "delta = 0"},
            "location_pattern": "zero",
            "sources": ["arXiv:2312.16282 §2-3", "arXiv:2405.09608 §2"],
        },
        "TrPhi3_shifted": {
            "poles": "all_chords", "zeros": "all_registered_loci", "splits": "rank1_all_registered",
            "deformation": {"parameter": "delta", "value": "nonzero-finite",
                            "range": "delta in Q \\ {0}; X_ee -> X_ee + delta, X_oo -> X_oo - delta"},
            "location_pattern": "parity_sign",
            "sources": ["arXiv:2312.16282 §7", "arXiv:2401.05483 eq. (1)"],
        },
        "NLSM": {
            "poles": "opposite_parity_chords", "zeros": "all_registered_loci",
            "splits": "rank1_all_registered",
            "deformation": {"parameter": "delta", "value": "infinity",
                            "range": "delta -> inf limit, eq. (2) of arXiv:2401.05483"},
            "location_pattern": "zero",
            "sources": ["arXiv:2401.05483 eq. (2)", "arXiv:2312.16282 §4, §7.4"],
        },
    },
    "gvar_rules": [
        {"rule_id": "GV-TRPHI3", "family": "TrPhi3",
         "requires_predicates": ["poles:oe", "poles:ee", "poles:oo", "zeros:all_registered",
                                 "splits:rank1", "delta:zero"]},
        {"rule_id": "GV-TRPHI3-SHIFTED", "family": "TrPhi3_shifted",
         "requires_predicates": ["poles:oe", "poles:ee", "poles:oo", "zeros:all_registered",
                                 "splits:rank1", "delta:nonzero"]},
        {"rule_id": "GV-NLSM", "family": "NLSM",
         "requires_predicates": ["poles:oe", "zeros:all_registered", "splits:rank1"],
         "forbids_predicates": ["poles:ee", "poles:oo"]},
    ],
}


def write_known_grammars(path: str = KNOWN_GRAMMARS_PATH) -> None:
    from pir.canonical import canonical_json
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        f.write(canonical_json(KNOWN_GRAMMARS) + "\n")


def load_known_grammars(path: str = KNOWN_GRAMMARS_PATH) -> Dict:
    with open(path) as f:
        return json.load(f)


# --------------------------------------------------------------------------- #
# expected features of a registered family at multiplicity n                  #
# --------------------------------------------------------------------------- #
def _chords(n: int) -> List[Tuple[int, int]]:
    return [(i, j) for i in range(1, n + 1) for j in range(i + 2, n + 1) if not (i == 1 and j == n)]


def expected_features(fam: Dict, n: int, registered: Dict) -> Dict:
    ch = _chords(n)
    if fam["poles"] == "all_chords":
        poles = ch
    elif fam["poles"] == "opposite_parity_chords":
        poles = [c for c in ch if (c[0] + c[1]) % 2 == 1]
    else:
        raise ValueError(fam["poles"])
    loci = [[tuple(int(x) for x in p.split(",")) for p in pairs]
            for pairs in registered["loci"][str(n)].values()]
    splits = registered["splits"].get(str(n), {})
    split_structure = tuple(sorted(
        (canonical_chord_set(n, [(t[0], t[1]), (t[1], t[2]), (t[0], t[2])]), 1)
        for t in splits.values()))
    dv = fam["deformation"]["value"]
    dclass = "zero" if dv == "0/1" else "absent" if dv == "infinity" else "nonzero"
    if fam["location_pattern"] == "zero":
        pattern = {c: "0/1" for c in poles}
    else:
        pattern = {c: ("-1/1" if (c[0] % 2 == 0 and c[1] % 2 == 0)
                       else "1/1" if (c[0] % 2 == 1 and c[1] % 2 == 1) else "0/1") for c in poles}
    return {"n": n, "pole_set": canonical_chord_set(n, poles),
            "zero_loci": canonical_pair_sets(n, loci),
            "split_structure": split_structure,
            "pole_location_pattern": canonical_labelled(n, pattern),
            "delta_class": dclass}


def _hashable(feat: Dict) -> Dict:
    return {k: feat[k] for k in ("n", "pole_set", "zero_loci", "split_structure",
                                 "pole_location_pattern")}


def build_db(known: Dict, n: int, registered: Dict) -> Tuple[KnownGrammarDB, Dict[str, Dict]]:
    db = KnownGrammarDB()
    feats = {}
    for name, fam in known["families"].items():
        if name == "NLSM" and n % 2:
            continue
        f = expected_features(fam, n, registered)
        feats[name] = f
        db.register(name, _hashable(f))
    return db, feats


# --------------------------------------------------------------------------- #
# predicates, similarity, verdict                                              #
# --------------------------------------------------------------------------- #
def predicates(feat: Dict, registered: Dict) -> set:
    n = feat["n"]
    pred = set()
    kinds = set()
    for (i, j) in feat["pole_set"]:
        kinds.add("ee" if (i % 2 == 0 and j % 2 == 0) else "oo" if (i % 2 == 1 and j % 2 == 1) else "oe")
    pred |= {f"poles:{k}" for k in kinds}
    if feat["n_vanishing_loci"] == feat["n_registered_loci"] and feat["n_registered_loci"] > 0:
        pred.add("zeros:all_registered")
    if feat["split_structure"] and all(r == 1 for _, r in feat["split_structure"]):
        pred.add("splits:rank1")
    elif not feat["split_structure"]:
        pred.add("splits:rank1")           # no registered split at this n: vacuous (n=4)
    d = feat["delta"]
    if d["verdict"] == "FORCED":
        pred.add("delta:zero" if Fraction(d["value"]) == 0 else "delta:nonzero")
    else:
        pred.add("delta:absent")
    return pred


def scores(feat: Dict, db_feats: Dict[str, Dict], keys) -> Dict:
    A = atoms(feat, keys)
    sims = {}
    for name, ef in db_feats.items():
        B = atoms(ef, keys)
        sims[name] = Fraction(len(A & B), len(A | B)) if (A | B) else Fraction(1)
    ranked = sorted(sims.items(), key=lambda kv: (-kv[1], kv[0]))
    best, s_best = ranked[0]
    ties = [nm for nm, s in ranked if s == s_best]
    if len(ranked) > 1 and len(ties) == 1:
        runner = ranked[1][0]
        disc = atoms(db_feats[best], keys) ^ atoms(db_feats[runner], keys)
        agree = sum(1 for a in disc if (a in A) == (a in atoms(db_feats[best], keys)))
        conf = Fraction(agree, len(disc)) if disc else Fraction(1)
    elif len(ties) > 1:
        conf = Fraction(0)
    else:
        conf = Fraction(1)
    return {"best": best, "ties": ties, "similarity": s_best, "confidence": conf,
            "all": {k: fmt(v) for k, v in sims.items()}}


def identify(feat: Dict, known: Dict, registered: Dict, threshold: Fraction = Fraction(9, 10)) -> Dict:
    n = feat["n"]
    db, db_feats = build_db(known, n, registered)
    obs = _hashable(feat)
    lookup = db.identify(obs)
    pred = predicates(feat, registered)
    cands = apply_rules(known["gvar_rules"], pred)
    lat = evaluate(cands, declared_interventions=[])
    prim = scores(feat, db_feats, ("pole_set", "zero_loci", "split_structure", "delta_class"))
    inv = scores(feat, db_feats, ("pole_set", "zero_loci", "split_structure"))
    sec = scores(feat, db_feats, ("zero_loci", "split_structure"))

    # --- pre-registered adjudication (prereg-002 §GID) ----------------------
    if len(lat.compatible) == 1 and prim["similarity"] >= threshold and prim["confidence"] >= threshold \
            and lookup.get("status") == "IDENTIFIED" and lookup["grammar_id"] == lat.compatible[0]:
        verdict, cause, cls = "FORCED", None, None
        selected = lat.compatible[0]
    elif len(lat.compatible) >= 2:
        verdict, cause, cls = "OBSERVATIONALLY_EQUIVALENT", None, sorted(lat.compatible)
        selected = None
    elif len(prim["ties"]) >= 2 and feat["delta"]["verdict"] == "FORCED":
        verdict, cause, cls = "AMBIGUOUS", None, None
        selected = None
    else:
        verdict = "NONIDENTIFIABLE"
        cause = ("no compatible grammar" if not lat.compatible else
                 f"threshold not met (similarity={fmt(prim['similarity'])}, confidence={fmt(prim['confidence'])})")
        cls = None
        selected = None

    fact = lattice_fact(lat, eqc_id=f"b15-gid-n{n}")
    fact = replace(fact, namespace="analyst", layer="DOMAIN",
                   measurement_interface=("label-stripped exact dataset data/b15/zero_dataset_v1.json",))
    return {
        "n": n, "full": full_fingerprint(obs), "specific": specific_fingerprint(obs),
        "db_lookup": lookup, "predicates": sorted(pred),
        "candidates": [{"family": c.family, "rule_id": c.rule_id} for c in cands],
        "lattice": {"verdict": lat.verdict, "compatible": lat.compatible, "detail": lat.detail},
        "primary": {"correlator": "b15.gid.jaccard_atoms_v1", "best": prim["best"],
                    "similarity": fmt(prim["similarity"]), "confidence": fmt(prim["confidence"]),
                    "all": prim["all"], "ties": prim["ties"]},
        "invariants_only": {"correlator": "b15.gid.invariants_only_v1", "best": inv["best"],
                            "similarity": fmt(inv["similarity"]), "confidence": fmt(inv["confidence"]),
                            "all": inv["all"], "ties": inv["ties"]},
        "secondary": {"correlator": "b15.gid.zeros_only_v1", "best": sec["best"],
                      "similarity": fmt(sec["similarity"]), "confidence": fmt(sec["confidence"]),
                      "all": sec["all"], "ties": sec["ties"]},
        "delta": feat["delta"],
        "verdict": verdict, "verdict_cause": cause, "verdict_class": cls,
        "selected": selected,
        "pir_fact": fact.to_dict(),
    }
