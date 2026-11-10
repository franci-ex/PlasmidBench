#!/usr/bin/env python3
"""Compare PlasmidBench's SnapGene .dna and ABI .ab1 readers with Biopython's.

Usage:  python3 -I cmpreal.py <dumpreal output.json> <dir> [<dir> ...]

Both formats are binary and undocumented by their vendors; Biopython's
readers were written from the same reverse engineering but independently of
this project, which is what makes the comparison worth anything.  Biopython
is the oracle and this script never imports PlasmidBench.

For a .dna file: the sequence, the topology, and every feature's span.
For an .ab1 file: the called bases and the per-base quality values.
"""
import json
import os
import sys
import warnings

from Bio import SeqIO

warnings.simplefilter("ignore")

DUMP = sys.argv[1]
DIRS = sys.argv[2:]

mine = json.load(open(DUMP))
ok = bad = notes = 0
problems = []


def find(name):
    for d in DIRS:
        p = os.path.join(d, name)
        if os.path.exists(p):
            return p
    return None


for name, m in sorted(mine.items()):
    path = find(name)
    if path is None:
        continue
    fmt = "snapgene" if m.get("kind") == "dna" else "abi"
    try:
        theirs = list(SeqIO.parse(path, fmt))
    except Exception as e:
        # A file Biopython refuses says nothing about this reader; but if
        # PlasmidBench accepted it, that is worth seeing.
        notes += 1
        state = "also refused" if "error" in m else f"accepted it ({m.get('length')} bases)"
        problems.append(f"  note {name}: Biopython refuses it ({e}); PlasmidBench {state}")
        continue

    if "error" in m:
        bad += 1
        problems.append(f"  FAIL {name}: PlasmidBench raised {m['error']}, Biopython read it")
        continue
    if len(theirs) != 1:
        notes += 1
        problems.append(f"  note {name}: Biopython found {len(theirs)} records")
        continue
    t = theirs[0]

    tseq = str(t.seq).upper()
    if m["seq"] != tseq:
        n = min(len(m["seq"]), len(tseq))
        at = next((j for j in range(n) if m["seq"][j] != tseq[j]), n)
        bad += 1
        problems.append(
            f"  FAIL {name}: sequence differs at {at} "
            f"(mine {len(m['seq'])}, Biopython {len(tseq)})")
        continue
    ok += 1

    if m["kind"] == "dna":
        ttopo = "circular" if t.annotations.get("topology") == "circular" else "linear"
        if m["topology"] != ttopo:
            bad += 1
            problems.append(f"  FAIL {name}: topology {m['topology']}, Biopython {ttopo}")
        else:
            ok += 1
        tf = sorted((f.type, int(f.location.start), int(f.location.end))
                    for f in t.features if f.location is not None)
        mf = sorted((f["type"], f["start"], f["end"]) for f in m["features"])

        # SnapGene stores a feature that crosses the origin with its end
        # before its start. Biopython 1.88 flattens such a feature to
        # 0..length -- the whole plasmid -- which turns pFA-KanMX4's 19 bp
        # SP6 promoter into a 3,941 bp one. PlasmidBench keeps the wrap. The
        # test for which is right is the length: a promoter is not the size
        # of its plasmid. Each case is reported with that length so the
        # claim can be checked rather than taken on trust.
        L = len(tseq)
        flattened = []
        for f in m["features"]:
            if f["end"] <= f["start"] and (f["type"], 0, L) in tf:
                flattened.append((f["type"], f["label"], (L - f["start"]) + f["end"]))
        if flattened:
            tf = [x for x in tf if not any(x == (t0, 0, L) for t0, _, _ in flattened)]
            mf = [x for x in mf if not any(x[0] == t0 and x[2] <= x[1] for t0, _, _ in flattened)]
            for t0, lab, n in flattened:
                notes += 1
                problems.append(
                    f"  note {name}: {t0} \u201c{lab}\u201d crosses the origin and is "
                    f"{n} bp here; Biopython flattens it to the whole {L} bp molecule")

        if mf != tf:
            # SnapGene's own feature types are free text and the two readers
            # name a few of them differently; spans are the substance.
            if sorted(x[1:] for x in mf) == sorted(x[1:] for x in tf):
                notes += 1
                problems.append(f"  note {name}: same {len(mf)} spans, some types named differently")
            else:
                bad += 1
                problems.append(
                    f"  FAIL {name}: {len(mf)} feature spans vs Biopython {len(tf)}; "
                    f"first difference {next((a for a, b in zip(mf, tf) if a != b), None)}")
        else:
            ok += 1
    else:
        tq = t.letter_annotations.get("phred_quality")
        if tq is None:
            notes += 1
            problems.append(f"  note {name}: Biopython gives no quality values")
        elif list(m["qual"] or []) != list(tq):
            diff = sum(1 for a, b in zip(m["qual"] or [], tq) if a != b)
            bad += 1
            problems.append(
                f"  FAIL {name}: {diff} quality values differ "
                f"(mine {len(m['qual'] or [])}, Biopython {len(tq)})")
        else:
            ok += 1

print(f"  {ok} agreed, {bad} differed, {notes} noted across {len(mine)} binary files")
for p in problems[:25]:
    print(p)
sys.exit(1 if bad else 0)
