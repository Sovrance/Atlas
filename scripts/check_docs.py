#!/usr/bin/env python3
"""ATLAS-RH-ENG-007 §3.3 / §14 (WO-RH-37, WO-RH-45) — documentation truth gate.

    python3 scripts/check_docs.py

Why this is a merge gate and not cleanup
----------------------------------------
Stale instructions in a live document are a correctness defect, not documentation debt.
A future agent reading "WO-RH-05 is open, do not start degree 3" will regenerate or
reinterpret mathematics that is already certified. ENG-007 §0 calls this *statement drift*
and treats it as the primary remaining risk, ahead of numerics.

The policy this enforces (DOCUMENTATION_TRUTH_PASS.md):

* every canonical doc declares `status: CURRENT | HISTORICAL | SUPERSEDED`;
* a historical file MAY contain stale statements, because history is evidence and
  WO-RH-17 forbids deleting contrary evidence -- what it may not do is present them
  as live instructions;
* every local link in a canonical doc resolves;
* live docs agree with `certificates/work_order_status.json`, which is machine-readable
  and generated, not prose;
* live RH docs carry the finite-dimensional / no-RH-proof boundary;
* documents that call a version "current" agree with each other and with the tree.

Scope is `docs/docs_status.json`. Checks are deliberately scoped to canonical/live docs
rather than grepped blindly across history (§DOCUMENTATION_TRUTH_PASS, "Do not simply grep
blindly across history").
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
REGISTRY = REPO / "docs" / "docs_status.json"

#: Claims that were true before ENG-005/ENG-006 and are false now. A live document
#: asserting any of these is actively misleading a future agent.
#: Each entry is (regex, human explanation of what is actually true).
STALE_PATTERNS: list[tuple[str, str]] = [
    (r"WO-RH-05[^.\n]{0,80}?\bopen\b",
     "WO-RH-05 was recovered by ENG-005 (cutoff-free uniform E1)"),
    (r"interval E1 \*\*open\*\*",
     "the interval E1 for WO-RH-05 is closed; ENG-005 released it"),
    (r"degree[ -]?3[^.\n]{0,80}?(do not start|must not start|not begin|blocked)",
     "degree 3 is implemented and E1-certified by ENG-006, inertia (2,0,0)"),
    # Both word orders. The notebook checkpoint phrases it as "Do not start degree 3",
    # with the instruction ahead of its subject, which a subject-first pattern misses.
    (r"(do not start|must not start|do not begin)[^.\n]{0,60}?degree[ -]?3",
     "degree 3 is implemented and E1-certified by ENG-006, inertia (2,0,0)"),
    (r"(Only after|only after) WO-RH-05 closes",
     "WO-RH-05 is closed"),
    (r"degree-?1[^.\n]{0,60}?remains quarantined",
     "degree-1 E1 was recovered and PROMOTED by ENG-005"),
    (r"T\s*=\s*84[^.\n]{0,60}?remains quarantined",
     "T=84 point and uniform E1 were recovered and PROMOTED by ENG-005"),
]

#: A live RH document must state the claim boundary somewhere. Any one of these suffices.
BOUNDARY_MARKERS = [
    r"[Nn]o RH proof",
    r"rh_proof_claim.{0,12}false",
    r"does not (prove|claim).{0,30}Riemann",
    r"finite-dimensional[^.\n]{0,60}only",
]

STATUS_RE = re.compile(r"^\s*(?:[-*]\s*)?\**status\**\s*:\s*\**(CURRENT|HISTORICAL|SUPERSEDED)\**",
                       re.IGNORECASE | re.MULTILINE)
LINK_RE = re.compile(r"\[[^\]]*\]\(([^)]+)\)")


class Report:
    def __init__(self) -> None:
        self.errors: list[str] = []
        self.notes: list[str] = []

    def fail(self, doc: str, msg: str) -> None:
        self.errors.append(f"{doc}: {msg}")

    def note(self, msg: str) -> None:
        self.notes.append(msg)


def is_historical(rel: str, roots: list[str], text: str = "") -> bool:
    """Historical by location OR by its own declaration.

    DOCUMENTATION_TRUTH_PASS: "A historical file may contain stale statements if it is
    clearly historical and no live entrypoint presents them as instructions." A document
    that declares `status: HISTORICAL` at the top has made itself clearly historical, so
    the stale-claim check is not applied to it -- deleting the contrary evidence it records
    would violate WO-RH-17.
    """
    if any(rel.startswith(root) for root in roots):
        return True
    m = STATUS_RE.search(text)
    return bool(m and m.group(1).upper() in {"HISTORICAL", "SUPERSEDED"})


def strip_code_fences(text: str) -> str:
    """Remove fenced code blocks.

    A command example or a quoted historical transcript is not an instruction, and matching
    inside one produces false positives that would train people to ignore this gate.
    """
    return re.sub(r"```.*?```", "", text, flags=re.DOTALL)


def check_links(rel: str, text: str, rep: Report) -> None:
    doc_dir = (REPO / rel).parent
    for target in LINK_RE.findall(text):
        if re.match(r"^(https?:|mailto:|#)", target):
            continue
        clean = target.split("#", 1)[0].split("?", 1)[0].strip()
        if not clean:
            continue
        resolved = (doc_dir / clean).resolve()
        if not resolved.exists():
            alt = (REPO / clean.lstrip("/")).resolve()
            if not alt.exists():
                rep.fail(rel, f"broken local link -> {target}")


def check_status_marker(rel: str, text: str, rep: Report) -> None:
    if not STATUS_RE.search(text):
        rep.fail(rel, "no `status: CURRENT|HISTORICAL|SUPERSEDED` declaration "
                      "(ENG-007 §14 requires one on every canonical doc)")


def check_stale(rel: str, text: str, rep: Report) -> None:
    body = strip_code_fences(text)
    for pattern, truth in STALE_PATTERNS:
        m = re.search(pattern, body, re.IGNORECASE)
        if m:
            line = body[: m.start()].count("\n") + 1
            rep.fail(rel, f"line ~{line}: stale claim {m.group(0)!r} -- {truth}")


def check_boundary(rel: str, text: str, rep: Report) -> None:
    if not any(re.search(p, text) for p in BOUNDARY_MARKERS):
        rep.fail(rel, "live RH doc does not state the finite-dimensional / no-RH-proof "
                      "claim boundary")


def check_work_order_table(rel: str, text: str, status: dict, rep: Report) -> None:
    """Every WO-RH-nn the RH README mentions with a status word must agree with the
    generated status file.

    The check is deliberately narrow: it verifies that the README does not contradict the
    machine-readable source, not that it transcribes every key. A README is allowed to
    summarise; it is not allowed to disagree.
    """
    orders = status.get("orders", {})
    body = strip_code_fences(text)
    for m in re.finditer(r"\b(WO-RH-\d\d)\b(?P<rest>[^\n]{0,120})", body):
        wo, rest = m.group(1), m.group("rest")
        actual = orders.get(wo)
        if actual is None:
            continue
        claims_open = re.search(r"\b(open|blocked|quarantined|pending|not started)\b",
                                rest, re.IGNORECASE)
        settled = actual.startswith(("done", "recovered"))
        if claims_open and settled:
            line = body[: m.start()].count("\n") + 1
            rep.fail(rel, f"line ~{line}: says {wo} is "
                          f"{claims_open.group(0)!r} but work_order_status.json has "
                          f"{actual!r}")


def check_atlas_version(rel: str, text: str, rep: Report) -> None:
    """The root README must link the newest constant-atlas as its primary document."""
    versions = sorted(
        (p.name for p in (REPO / "docs").glob("constant-atlas-v*.md")),
        key=lambda n: [int(x) for x in re.findall(r"\d+", n)],
    )
    if not versions:
        return
    newest = versions[-1]
    head = text.split("\n## ", 1)[0]
    linked = re.findall(r"docs/(constant-atlas-v[\d.]+\.md)", head)
    if linked and newest not in linked:
        rep.fail(rel, f"'Core documents' links {linked[0]} but the current atlas in the "
                      f"tree is {newest}")
    declared = re.findall(r"[Cc]urrent atlas:\s*\[`?docs/(constant-atlas-v[\d.]+\.md)",
                          text)
    if declared and declared[0] != newest:
        rep.fail(rel, f"declares current atlas {declared[0]}, tree has {newest}")
    if declared and linked and declared[0] != linked[0]:
        rep.fail(rel, f"self-contradictory: header links {linked[0]}, body declares "
                      f"{declared[0]} current")


def main() -> int:
    if not REGISTRY.exists():
        print(f"missing {REGISTRY.relative_to(REPO)}", file=sys.stderr)
        return 1
    reg = json.loads(REGISTRY.read_text(encoding="utf-8"))
    rep = Report()

    status_path = REPO / reg["work_order_status_source"]
    if not status_path.exists():
        rep.fail(reg["work_order_status_source"], "missing machine-readable status source")
        status = {}
    else:
        status = json.loads(status_path.read_text(encoding="utf-8"))

    if status.get("rh_proof_claim") is not False:
        rep.fail(reg["work_order_status_source"], "rh_proof_claim must be false")

    historical = reg.get("historical_roots", [])
    rh_live = set(reg.get("rh_live_docs", []))

    for rel in reg["canonical_docs"]:
        path = REPO / rel
        if not path.exists():
            rep.fail(rel, "canonical document does not exist")
            continue
        text = path.read_text(encoding="utf-8")
        check_status_marker(rel, text, rep)
        check_links(rel, text, rep)
        doc_is_historical = is_historical(rel, historical, text)
        if not doc_is_historical:
            check_stale(rel, text, rep)
        if rel in rh_live and not doc_is_historical:
            check_boundary(rel, text, rep)
        if rel == "math/rh_weil/README.md" and status:
            check_work_order_table(rel, text, status, rep)
        if rel == "README.md":
            check_atlas_version(rel, text, rep)

    # The registry's own claims must match the tree.
    latest = reg.get("latest_completed_work_order", "")
    if status and latest:
        spec = status.get("eng_spec", "")
        tag = latest.replace("ATLAS-RH-", "")
        if tag not in spec:
            rep.fail("docs/docs_status.json",
                     f"latest_completed_work_order {latest!r} but work_order_status.json "
                     f"eng_spec is {spec[:60]!r}")

    if rep.errors:
        print("DOCS TRUTH GATE: FAIL", file=sys.stderr)
        for e in rep.errors:
            print(f"  - {e}", file=sys.stderr)
        print(f"\n{len(rep.errors)} problem(s). Stale live documentation is a correctness "
              "defect: a future agent will act on it.", file=sys.stderr)
        return 1

    print(f"DOCS TRUTH GATE: OK ({len(reg['canonical_docs'])} canonical docs; links, "
          "status markers, work-order agreement and claim boundary all check out)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
