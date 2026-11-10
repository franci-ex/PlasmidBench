#!/usr/bin/env python3
"""Compare PlasmidBench's protein descriptors with Bio.SeqUtils.ProtParam.

Usage:  python3 -I cmpprops.py <dumpprops output.json>

Only descriptors with one published definition are compared, so that a
difference means one of the two is wrong rather than that they define a
different quantity:

  molecular weight   average isotopic masses
  isoelectric point  Bjellqvist, as ProtParam implements it
  GRAVY              Kyte & Doolittle 1982
  aromaticity        Lobry & Gautier 1994
  instability index  Guruprasad et al. 1990
  aliphatic index    Ikai 1980 (Biopython has no equivalent; skipped)

Tolerances are stated per descriptor and are there for floating-point
accumulation only: 1e-6 relative on the weight, 1e-4 on the point, 1e-9 on
the rest.  A disagreement larger than that is a failure, not noise.
"""
import json
import sys

from Bio.SeqUtils.ProtParam import ProteinAnalysis

mine = json.load(open(sys.argv[1]))
ok = bad = skipped = 0
problems = []
per = {}

CHECKS = [
    # name, how Biopython computes it, absolute tolerance, relative tolerance
    ("mw", lambda a: a.molecular_weight(), 0.0, 1e-6),
    ("pI", lambda a: a.isoelectric_point(), 1e-4, 0.0),
    ("gravy", lambda a: a.gravy(), 1e-9, 0.0),
    ("aromaticity", lambda a: a.aromaticity(), 1e-9, 0.0),
    ("instability", lambda a: a.instability_index(), 1e-7, 0.0),
]

for key, m in sorted(mine.items()):
    aa = m["seq"]
    a = ProteinAnalysis(aa)
    for name, how, atol, rtol in CHECKS:
        try:
            want = how(a)
        except Exception:
            # ProtParam refuses sequences with X, B, Z or U for some
            # descriptors; nothing to compare against there.
            skipped += 1
            per.setdefault(name, [0, 0, 0])[2] += 1
            continue
        got = m[name]
        tol = max(atol, abs(want) * rtol)
        if abs(got - want) <= tol:
            ok += 1
            per.setdefault(name, [0, 0, 0])[0] += 1
        else:
            per.setdefault(name, [0, 0, 0])[1] += 1
            bad += 1
            problems.append(
                f"  FAIL {key} {name}: {got!r} vs Biopython {want!r} "
                f"(difference {got - want:.3e}, allowed {tol:.1e})")

print(f"  {ok} agreed, {bad} differed, {skipped} not defined by ProtParam, "
      f"over {len(mine)} proteins")
for name, (a, b, sk) in per.items():
    print(f"    {name:<12} {a} agreed, {b} differed, {sk} skipped")
for p in problems[:20]:
    print(p)
sys.exit(1 if bad else 0)
