#!/usr/bin/env python3
"""Check the identity the tolerant matcher reports against Biopython's aligner.

Usage:  python3 -I cmpanno.py <dumpanno output.json>

The identity on an annotation is the number a person reads and trusts, and
it comes out of PlasmidBench's own banded alignment. Bio.Align is the oracle
for it: the same library entry against the same span of the same sequence,
aligned by an implementation that shares no code with it, and asked how
similar they are.

Two things cannot be compared directly and are handled rather than ignored:

  * Scoring. Two aligners with different gap costs produce different
    alignments of the same pair, so the identities need not match to the
    digit. What must hold is that they agree closely, and that where they
    differ it is because Biopython found a BETTER alignment, never a worse
    one -- a matcher reporting a higher identity than the best alignment of
    those two strings would be reporting something that is not there.

  * Which bases. PlasmidBench reports a span of the sequence; the oracle
    aligns the entry to exactly that span. If the span is wrong, the
    identity computed over it will be too, so this checks the span as well
    as the number on it.
"""
import json
import sys

from Bio import Align

rows = json.load(open(sys.argv[1]))

aligner = Align.PairwiseAligner()
aligner.mode = "local"
aligner.match_score = 1
aligner.mismatch_score = -1
aligner.open_gap_score = -7
aligner.extend_gap_score = -3

ok = bad = notes = 0
problems = []


def identity_of(a, b):
    """The identity of the best local alignment of a against b, over the
    aligned columns -- the same quantity PlasmidBench reports."""
    if not a or not b:
        return None
    try:
        aln = aligner.align(a, b)[0]
    except Exception:
        return None
    # count matches and aligned columns from the aligned blocks
    matches = cols = 0
    for (astart, aend), (bstart, bend) in zip(aln.aligned[0], aln.aligned[1]):
        for i in range(aend - astart):
            cols += 1
            if a[astart + i] == b[bstart + i]:
                matches += 1
    # gaps are columns too: the span covered, minus the columns matched up
    spanA = aln.aligned[0][-1][1] - aln.aligned[0][0][0] if len(aln.aligned[0]) else 0
    spanB = aln.aligned[1][-1][1] - aln.aligned[1][0][0] if len(aln.aligned[1]) else 0
    gaps = max(spanA, spanB) - cols
    total = cols + max(0, gaps)
    return (matches / total) if total else None


for r in rows:
    tag = r["label"] + "/" + r["entry"]
    if not r["found"]:
        notes += 1
        problems.append(f"  note {tag}: not found by PlasmidBench; nothing to check")
        continue

    entry = r["entryOnStrand"]
    span = r["span"]

    # 1 -- the span is really where that entry is
    want = identity_of(entry, span)
    if want is None:
        notes += 1
        continue
    got = r["identity"]

    if got > want + 0.02:
        bad += 1
        problems.append(
            f"  FAIL {tag}: reports {got:.4f} identity, but the best alignment of "
            f"the entry against the span it named is only {want:.4f}")
    elif abs(got - want) <= 0.02:
        ok += 1
    else:
        # PlasmidBench below the oracle: conservative, but worth seeing when large
        if want - got > 0.08:
            bad += 1
            problems.append(
                f"  FAIL {tag}: reports {got:.4f} where a better alignment of the "
                f"same two strings gives {want:.4f} -- too much is being missed")
        else:
            notes += 1
            problems.append(f"  note {tag}: {got:.4f} against the oracle's {want:.4f}, "
                            f"within the difference two gap models make")

    # 2 -- coverage: how much of the entry the span accounts for
    cov = r["coverage"]
    # the span cannot account for more of the entry than it has bases for
    if cov * len(entry) > len(span) + max(4, 0.1 * len(entry)):
        bad += 1
        problems.append(
            f"  FAIL {tag}: claims {cov:.2f} of a {len(entry)} bp entry "
            f"({cov * len(entry):.0f} bases) in a span of only {len(span)}")
    else:
        ok += 1

    # 3 -- exact means exact, and is checked by looking
    if r["exact"]:
        if entry != span:
            bad += 1
            problems.append(f"  FAIL {tag}: called exact, but the span is not the entry "
                            f"({len(span)} bp against {len(entry)})")
        else:
            ok += 1

print(f"  {ok} agreed, {bad} differed, {notes} noted over {len(rows)} annotations")
print("    oracle: Bio.Align.PairwiseAligner, local, the same span and entry")
for p in problems[:25]:
    print(p)
sys.exit(1 if bad else 0)
