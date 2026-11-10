#!/usr/bin/env python3
"""Compare the commercial suppliers PlasmidBench lists for each enzyme with
REBASE as Biopython distributes it.

Usage:  python3 -I cmpsuppl.py <rebase-mine.json>

This matters at the bench rather than in the arithmetic: the enzyme panel
offers to filter by what a lab can actually buy, and a wrong supplier letter
sends someone looking for a tube from a company that does not sell it.

Biopython's Restriction_Dictionary carries REBASE's supplier table, which is
the same source PlasmidBench's letters came from but a separate copy of it.
Both the letter-to-company mapping and the per-enzyme letters are checked.
"""
import json
import sys

from Bio import Restriction
from Bio.Restriction.Restriction_Dictionary import suppliers

dump = json.load(open(sys.argv[1]))
mine_names = dump["supplierNames"]
ok = bad = 0
problems = []

# ---- the letter-to-company table -----------------------------------------
for code, name in sorted(mine_names.items()):
    if code not in suppliers:
        bad += 1
        problems.append(f"  FAIL supplier {code}: not a REBASE supplier code")
        continue
    want = suppliers[code][0]
    if name.strip() != want.strip():
        bad += 1
        problems.append(f"  FAIL supplier {code}: “{name}” vs REBASE “{want}”")
    else:
        ok += 1

# ---- the letters on each enzyme ------------------------------------------
missing_codes = set()
for e in dump["enzymes"]:
    name = e["name"]
    t = getattr(Restriction, name, None)
    if t is None:
        continue
    want = set(t.suppl)
    got = set(e["suppliers"] or "")
    # Only the companies PlasmidBench's table knows about can be compared;
    # a letter for a company it does not list is a gap in the names table,
    # counted separately rather than scored as a wrong supplier.
    unknown = want - set(mine_names)
    missing_codes |= unknown
    want &= set(mine_names)
    if got != want:
        bad += 1
        extra = "".join(sorted(got - want))
        lost = "".join(sorted(want - got))
        problems.append(
            f"  FAIL {name}: suppliers “{''.join(sorted(got))}” vs REBASE "
            f"“{''.join(sorted(want))}”"
            + (f"; lists {extra} that REBASE does not" if extra else "")
            + (f"; misses {lost}" if lost else ""))
    else:
        ok += 1

print(f"  {ok} agreed, {bad} differed over {len(mine_names)} suppliers "
      f"and {len(dump['enzymes'])} enzymes")
if missing_codes:
    named = ", ".join(f"{c} ({suppliers[c][0]})" for c in sorted(missing_codes))
    print(f"    REBASE codes PlasmidBench does not name: {named}")
for p in problems[:25]:
    print(p)
sys.exit(1 if bad else 0)
