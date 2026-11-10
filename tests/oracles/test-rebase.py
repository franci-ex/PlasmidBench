#!/usr/bin/env python3
"""Compare PlasmidBench's restriction enzymes and cut positions with Biopython.

Usage:  python3 -I test-rebase.py [rebase-mine.json]

Bio.Restriction is the oracle: its table and PlasmidBench's both come from
REBASE, but the two searches are written differently -- Biopython compiles a
regular expression for each enzyme, PlasmidBench walks the sequence against
IUPAC codes -- so agreement on where an enzyme cuts is real evidence and not
two copies of one mistake.

Three things are checked for every enzyme the two libraries share:

  the recognition sequence, letter for letter;
  the position of the cut within the site, on both strands;
  every cut position in five sequences, read as linear and as circular.

An enzyme only one library carries is counted and named, never passed over
in silence: a missing enzyme is a gap in the table, and a table nobody has
counted is not a validated table.
"""
import json
import os
import sys

from Bio.Seq import Seq
from Bio import Restriction

HERE = os.path.dirname(os.path.abspath(__file__))
dump = json.load(open(sys.argv[1] if len(sys.argv) > 1
                      else os.path.join(HERE, "rebase-mine.json")))

mine = {e["name"]: e for e in dump["enzymes"]}
theirs = {str(e): e for e in Restriction.AllEnzymes}

shared = sorted(set(mine) & set(theirs))
only_mine = sorted(set(mine) - set(theirs))
only_theirs = sorted(set(theirs) - set(mine))

ok = bad = 0
problems = []

# ---- 1. the recognition sequence and the cut inside it --------------------
for name in shared:
    m, t = mine[name], theirs[name]
    if m["site"].upper() != str(t.site).upper():
        bad += 1
        problems.append(f"  FAIL {name}: site {m['site']} vs REBASE {t.site}")
    else:
        ok += 1
    # Biopython's fst5 is the distance from the start of the site to the cut
    # in the top strand; fst3 is measured from the end of the site. Enzymes
    # that do not cut (or cut unpredictably) have fst5 None.
    if t.fst5 is None:
        continue
    n = len(str(t.site))
    # An enzyme that cuts within its site stores top/bot from the start of
    # the site; a type IIS enzyme stores o1/o2 from the end of it. Converting
    # here, on the oracle's side, keeps the dump a plain report.
    if m.get("top") is not None:
        got5, want5 = m["top"], t.fst5
        got3, want3 = m["bot"], n + t.fst3
        how = "within its site"
    else:
        got5, want5 = n + m["o1"], t.fst5
        got3, want3 = m["o2"], t.fst3
        how = "beyond its site"
    if got5 != want5:
        bad += 1
        problems.append(f"  FAIL {name} ({how}): top-strand cut at {got5} vs REBASE {want5}")
    else:
        ok += 1
    if got3 != want3:
        bad += 1
        problems.append(
            f"  FAIL {name} ({how}): bottom-strand cut at {got3} vs REBASE {want3}")
    else:
        ok += 1

# ---- 2. where each enzyme actually cuts -----------------------------------
checked = 0
for key, seq in dump["sequences"].items():
    s = Seq(seq)
    for name in shared:
        t = theirs[name]
        if t.fst5 is None:
            continue
        got = dump["cuts"].get(key, {}).get(name)
        got_lin = got["linear"] if got else []
        got_cir = got["circular"] if got else []
        # search() returns 1-based positions of the base after the break
        want_lin = sorted(p - 1 for p in t.search(s, linear=True))
        want_cir = sorted(p - 1 for p in t.search(s, linear=False))
        checked += 2
        if sorted(got_lin) != want_lin:
            bad += 1
            extra = sorted(set(got_lin) - set(want_lin))[:3]
            missing = sorted(set(want_lin) - set(got_lin))[:3]
            problems.append(
                f"  FAIL {name} on {key} linear: {len(got_lin)} cuts vs "
                f"{len(want_lin)}; extra {extra}, missing {missing}")
        else:
            ok += 1
        if sorted(got_cir) != want_cir:
            bad += 1
            extra = sorted(set(got_cir) - set(want_cir))[:3]
            missing = sorted(set(want_cir) - set(got_cir))[:3]
            problems.append(
                f"  FAIL {name} on {key} circular: {len(got_cir)} cuts vs "
                f"{len(want_cir)}; extra {extra}, missing {missing}")
        else:
            ok += 1

print(f"  {ok} agreed, {bad} differed over {len(shared)} shared enzymes "
      f"and {checked} cut lists")
print(f"    PlasmidBench carries {len(mine)}, Biopython's REBASE {len(theirs)}; "
      f"{len(only_mine)} only here, {len(only_theirs)} only there")
if only_mine:
    print("    only here:  " + ", ".join(only_mine[:12]) + ("..." if len(only_mine) > 12 else ""))
for p in problems[:25]:
    print(p)
sys.exit(1 if bad else 0)
