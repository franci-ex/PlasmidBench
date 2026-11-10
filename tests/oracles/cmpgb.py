#!/usr/bin/env python3
"""Compare PlasmidBench's GenBank reader against Biopython's, file by file.

Usage:  python3 -I cmpgb.py <corpus dir> <dumpgb output.json> [skip,list]

Biopython is the oracle here and this script is the only place the expected
values come from: it reads the same files with Bio.SeqIO and never imports or
executes any PlasmidBench code.  The JavaScript half (dumpgb.js) reports and
does not compare, so neither half can be tuned to agree with the other.

What is checked, per record: the sequence itself, the length, the topology,
how many features there are, and each feature's type and location -- the
location in the same 0-based end-exclusive form both libraries hold
internally, including the parts of a join() in order.

Known, deliberate differences are listed in DIVERGENCES with the reason; they
are reported as notes, not failures, so the count stays honest.
"""
import hashlib
import json
import os
import sys
import warnings

from Bio import BiopythonParserWarning, SeqIO

warnings.simplefilter("ignore", BiopythonParserWarning)
warnings.simplefilter("ignore", BiopythonParserWarning.__base__)

CORPUS, DUMP = sys.argv[1], sys.argv[2]
SKIP = set(x for x in (sys.argv[3].split(",") if len(sys.argv) > 3 else []) if x)

# Files where the two libraries legitimately disagree, and why.  Each entry
# has to name a reason; "it fails" is not one.
DIVERGENCES = {
    "addgene-plasmid-11664-sequence-180430.gbk": {
        "topology": (
            "Addgene and SnapGene write \"linear\" on the LOCUS line of plasmid "
            "exports and then annotate enhancer 8569..276, which runs off the "
            "end of an 8,596 bp molecule and back to base 1. Taken at its "
            "word the record is linear and that feature is impossible, so "
            "Biopython keeps it as two parts on a linear sequence while "
            "PlasmidBench concludes the molecule is circular and says so "
            "(topologyFromFeature). A backwards span cannot exist on a linear "
            "molecule, so it is the stronger evidence of the two."),
    },
}


def bioloc(feat):
    """The feature's location as start..end:strand, joins as parts in order."""
    if feat.location is None:       # a location Biopython could not parse
        return "?"
    parts = list(feat.location.parts)
    out = []
    for p in parts:
        strand = "-" if p.strand == -1 else "+"
        out.append(f"{int(p.start)}..{int(p.end)}:{strand}")
    return ",".join(out)


def sha(text):
    return hashlib.sha1(str(text).upper().encode()).hexdigest()[:16]


def topology_of(rec):
    """GenBank puts the topology on the LOCUS line; Biopython keeps it in
    annotations, absent meaning linear."""
    t = rec.annotations.get("topology", "linear")
    return str(t).lower()


mine = json.load(open(DUMP))
ok = bad = notes = 0
problems = []

for fname, recs in sorted(mine.items()):
    if fname in SKIP:
        continue
    path = os.path.join(CORPUS, fname)
    if not os.path.exists(path):
        continue
    try:
        theirs = list(SeqIO.parse(path, "genbank"))
    except Exception as e:                      # a file Biopython itself rejects
        notes += 1
        problems.append(f"  note {fname}: Biopython cannot read it ({e})")
        continue

    if isinstance(recs, dict) and "error" in recs:
        bad += 1
        problems.append(f"  FAIL {fname}: PlasmidBench raised {recs['error']}")
        continue

    if len(recs) != len(theirs):
        bad += 1
        problems.append(f"  FAIL {fname}: {len(recs)} records, Biopython says {len(theirs)}")
        continue

    for i, (m, t) in enumerate(zip(recs, theirs)):
        tag = fname if len(recs) == 1 else f"{fname}[{i}]"
        # A GenBank record may carry no ORIGIN at all (a CONTIG join, or a
        # header-only entry). Biopython reports that as an undefined
        # sequence; there is nothing to compare, so only the features are.
        try:
            tseq = str(t.seq).upper()
        except Exception:
            notes += 1
            problems.append(f"  note {tag}: no sequence in the file, features only")
            tseq = None

        if tseq is None:
            pass
        elif m["seq"] != tseq:
            # locate the first difference so the report is actionable
            n = min(len(m["seq"]), len(tseq))
            at = next((j for j in range(n) if m["seq"][j] != tseq[j]), n)
            bad += 1
            problems.append(
                f"  FAIL {tag}: sequence differs at {at} "
                f"(mine {len(m['seq'])} bp, Biopython {len(tseq)} bp)")
            continue
        else:
            ok += 1

        if m["topology"] != topology_of(t):
            why = DIVERGENCES.get(fname, {}).get("topology")
            if why:
                notes += 1
                problems.append(
                    f"  note {tag}: topology {m['topology']} vs Biopython "
                    f"{topology_of(t)} — known: {why}")
            else:
                bad += 1
                problems.append(
                    f"  FAIL {tag}: topology {m['topology']}, Biopython {topology_of(t)}")
        else:
            ok += 1

        # Biopython keeps every feature including source; so does the engine.
        # A feature PlasmidBench refused is listed in "skipped" and reported
        # as a note: refusing a location it cannot place on this sequence --
        # a part in another accession, a negative coordinate -- and saying so
        # is a different thing from losing it.
        tf = [(f.type, bioloc(f)) for f in t.features]
        mf = [(f["type"], f["loc"]) for f in m["features"]]
        if len(mf) != len(tf):
            skipped = m.get("skipped") or []
            if len(mf) + len(skipped) == len(tf):
                notes += 1
                problems.append(
                    f"  note {tag}: {len(skipped)} feature(s) declined and listed: "
                    + "; ".join(skipped[:3]))
            else:
                bad += 1
                problems.append(
                    f"  FAIL {tag}: {len(mf)} features + {len(skipped)} declined, "
                    f"Biopython {len(tf)}")
        else:
            ok += 1
            types_m = [a for a, _ in mf]
            types_t = [a for a, _ in tf]
            if types_m != types_t:
                bad += 1
                problems.append(f"  FAIL {tag}: feature types differ")
            elif tseq is not None:
                # the test that matters: the bases each feature covers
                wrong = []
                for fm, ft in zip(m["features"], t.features):
                    if ft.location is None:
                        continue
                    try:
                        want = sha(str(ft.extract(t.seq)))
                    except Exception:
                        continue
                    if fm["extract"] != want:
                        wrong.append((fm["type"], fm["loc"], bioloc(ft)))
                if wrong:
                    bad += 1
                    problems.append(
                        f"  FAIL {tag}: {len(wrong)} feature(s) extract different bases, "
                        f"first {wrong[0][0]} at {wrong[0][1]} (Biopython {wrong[0][2]})")
                else:
                    ok += 1

print(f"  {ok} agreed, {bad} differed, {notes} skipped by Biopython "
      f"across {len(mine)} files")
for p in problems[:25]:
    print(p)
sys.exit(1 if bad else 0)
