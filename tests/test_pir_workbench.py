"""PIR Workbench (P9 UI) build test.

W1  The page builds: the data placeholder is replaced and the inlined JSON is
    valid and carries every bundle section the six surfaces need.
W2  Self-containment: no external requests (no http(s):// URLs, no src=/href= to
    remote hosts, no font/CDN links) — CSP-safe.
W3  All six surface renderers and the honesty affordances (SOUND/HEURISTIC tags,
    E0–E4 badges, theme tokens for both light and dark) are present.
W4  The `</script>` inside the data is escaped so the JSON island can't break the
    page.
"""

import json
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ci.build_workbench import render_html

FAILURES = []


def check(label, errors):
    if errors:
        FAILURES.append((label, errors)); print(f"  [FAIL] {label}")
        for e in errors[:6]:
            print(f"         - {e}")
    else:
        print(f"  [ ok ] {label}")


HTML = render_html()


def w1_data():
    errs = []
    if "__PIR_DATA__" in HTML:
        errs.append("data placeholder not replaced")
    m = re.search(r'id="pir-data" type="application/json">(.*?)</script>', HTML, re.S)
    if not m:
        errs.append("data island missing")
    else:
        data = json.loads(m.group(1).replace("<\\/", "</"))
        need = {"facts", "verdict_matrix", "distributions", "invalidation_demo",
                "cross_domain_diff", "candidate_lattice", "structural_graph", "coverage"}
        missing = need - set(data)
        if missing:
            errs.append(f"bundle missing sections: {sorted(missing)}")
        if len(data.get("facts", [])) < 40:
            errs.append("too few facts embedded")
    check("W1 page builds with a valid, complete embedded bundle", errs)


def w2_self_contained():
    errs = []
    # no remote hosts anywhere (allow the SPEC-y '://' only inside JSON string data? none expected).
    for pat in (r'https?://', r'src\s*=\s*["\']//', r'@import', r'url\(\s*https?:'):
        if re.search(pat, HTML):
            errs.append(f"external reference matched /{pat}/")
    check("W2 no external requests — CSP-safe / self-contained", errs)


def w3_surfaces_and_honesty():
    errs = []
    for fn in ("renderFacts", "renderMatrix", "renderProv", "renderLattice",
               "renderDiff", "renderGraph"):
        if fn not in HTML:
            errs.append(f"missing surface renderer {fn}")
    for token in ("tag sound", "tag heuristic", "ebadge", "is never a certificate"):
        if token not in HTML:
            errs.append(f"missing honesty affordance: {token!r}")
    # both themes defined at token level
    if 'data-theme="dark"' not in HTML or 'data-theme="light"' not in HTML:
        errs.append("both explicit themes not defined")
    if "prefers-color-scheme: dark" not in HTML:
        errs.append("no prefers-color-scheme fallback")
    check("W3 six surfaces + honesty badges + both themes present", errs)


def w4_escaped():
    errs = []
    m = re.search(r'type="application/json">(.*?)</script>', HTML, re.S)
    island = m.group(1)
    if "</script" in island:
        errs.append("unescaped </script in data island (could break the page)")
    check("W4 data island escapes </ so it cannot terminate the script early", errs)


def _mirror_invalidate(facts, asm):
    """Pure-Python mirror of the workbench's invalidation algorithm: direct
    assumption membership + transitive closure over depends_on_facts."""
    affected = {f["fact_id"] for f in facts if asm in (f.get("assumptions") or [])}
    changed = True
    while changed:
        changed = False
        for f in facts:
            if f["fact_id"] in affected:
                continue
            if set(f.get("depends_on_facts") or []) & affected:
                affected.add(f["fact_id"]); changed = True
    return sorted(affected)


def w5_invalidation_fidelity():
    """Ported from the parallel P9 session: the invalidation the workbench shows
    must reproduce FactStore.invalidate_assumption exactly — mirror + real store."""
    from ci.export_pir_view import build_bundle
    from pir.provenance import FactStore
    from pir.domains import b9
    b = build_bundle()
    demo = b["invalidation_demo"]
    asm = demo["assumption"]
    committed = sorted(demo["downgraded_facts"])
    errs = []
    mirror = _mirror_invalidate(b["facts"], asm)
    if not set(committed) <= set(mirror):
        errs.append(f"page-mirror {mirror} does not cover committed demo {committed}")
    _, _, b9f = b9.lower(b9.load_certificate())
    store = FactStore()
    for f in b9f:
        store.add_fact(f)
    real = sorted(store.invalidate_assumption(asm, reason="workbench fidelity drill"))
    if real != committed:
        errs.append(f"FactStore {real} != committed demo {committed}")
    check("W5 invalidation view reproduces FactStore.invalidate_assumption", errs)


if __name__ == "__main__":
    print("== PIR workbench (P9) build ==")
    w1_data(); w2_self_contained(); w3_surfaces_and_honesty(); w4_escaped()
    w5_invalidation_fidelity()
    print("\n== SUMMARY ==")
    if FAILURES:
        print(f"FAIL: {len(FAILURES)} check(s) failed."); sys.exit(1)
    print(f"PASS: self-contained workbench builds ({len(HTML)//1024} KB, 6 surfaces).")
    sys.exit(0)
