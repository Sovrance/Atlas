"""Two-hash grammar fingerprints (Stage 2 / P5b).

The BinDiff / Function-ID pattern: every grammar carries two canonical hashes.

* **full fingerprint** — *similarity-invariant*: computed only from
  representation-invariant features (the R15 canonicalization set: positivity-cone
  data, rank sequence, symmetry group, factorization poset, flow exponents). Two
  presentations of the same physics share it.
* **specific fingerprint** — finer: includes representation-specific detail, so
  it disambiguates grammars that happen to share invariants.

A :class:`KnownGrammarDB` indexes grammars by full fingerprint. B8 "blind grammar
identification" then becomes a fingerprint lookup. When a full-fingerprint bucket
holds grammars with *different* specific fingerprints, the match is a collision
and is filed as REPRESENTATION_DEPENDENT (invariants agree, representations do
not) rather than a confident identification.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional

from .canonical import sha256_hex

# The invariant feature keys that define the *full* (similarity-invariant) hash.
# B15-SURF (additive): amplitude fingerprints — pole set, hidden-zero loci and
# split (factorization) structure are basis-independent invariants (R37/P3).
# ``_project`` skips absent keys, so every pre-B15 full-fingerprint hash is
# unchanged (regression test: tests/test_b15_gid.py::T0).
INVARIANT_KEYS = ("psd_signature", "rank_sequence", "symmetry_group",
                  "factorization_poset", "flow_exponents", "spectral_class",
                  "pole_set", "zero_loci", "split_structure")


def _project(features: Dict, keys) -> Dict:
    return {k: features[k] for k in keys if k in features}


def full_fingerprint(features: Dict) -> str:
    """Similarity-invariant hash over the canonical-invariant subset only."""
    inv = _project(features, INVARIANT_KEYS)
    return "fpf_" + sha256_hex({"invariant": inv})[:16]


def specific_fingerprint(features: Dict) -> str:
    """Finer hash over the full feature dict (invariants + representation)."""
    return "fps_" + sha256_hex({"all": features})[:16]


@dataclass
class KnownGrammarDB:
    """Index grammars by full fingerprint; lookups return an identification."""

    _by_full: Dict[str, List[Dict]] = field(default_factory=dict)

    def register(self, grammar_id: str, features: Dict) -> Dict:
        entry = {"grammar_id": grammar_id,
                 "full": full_fingerprint(features),
                 "specific": specific_fingerprint(features)}
        self._by_full.setdefault(entry["full"], []).append(entry)
        return entry

    def identify(self, features: Dict) -> Dict:
        """Blind identification by full fingerprint.

        Returns {verdict, ...}: a single match -> ``AMBIGUOUS``? No -- a clean
        single match is an identification (``matched`` grammar id); multiple
        entries with the same full but different specific hashes ->
        ``REPRESENTATION_DEPENDENT`` (a collision, not an identification)."""
        full = full_fingerprint(features)
        spec = specific_fingerprint(features)
        bucket = self._by_full.get(full, [])
        if not bucket:
            return {"status": "NOT_DETECTED", "full": full, "specific": spec}
        specifics = {e["specific"] for e in bucket}
        exact = [e for e in bucket if e["specific"] == spec]
        if exact and len(specifics) == 1:
            return {"status": "IDENTIFIED", "grammar_id": exact[0]["grammar_id"],
                    "full": full, "specific": spec}
        if exact:
            return {"status": "IDENTIFIED", "grammar_id": exact[0]["grammar_id"],
                    "full": full, "specific": spec,
                    "note": "full-hash bucket also holds other specifics"}
        # Same invariants, different representation-specific hash: a collision.
        return {"status": "REPRESENTATION_DEPENDENT", "full": full, "specific": spec,
                "colliding_grammars": sorted(e["grammar_id"] for e in bucket),
                "reason": "full (invariant) fingerprint matches but specific "
                          "fingerprint differs — invariants agree, representation "
                          "does not"}
