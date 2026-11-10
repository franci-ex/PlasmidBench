"""Compare PlasmidBench's embedded BLOSUM62 against Biopython's, cell by cell.

The matrix is a table of 576 numbers typed into a source file. A single
wrong cell would quietly change every protein alignment and every
similarity column in a figure, and would be invisible on inspection.

  node dumpb62.js | python3 -I cmpb62.py
"""
import json
import sys

from Bio.Align import substitution_matrices


def main():
    mine = json.load(sys.stdin)
    ref = substitution_matrices.load("BLOSUM62")
    alpha = mine["alphabet"]
    bad = 0
    n = 0
    for a in alpha:
        for b in alpha:
            n += 1
            want = int(ref[a, b])
            got = mine["matrix"].get(a + b)
            if got != want:
                bad += 1
                if bad <= 5:
                    print("  DIFF %s%s: mine %s, Biopython %s" % (a, b, got, want))
    verdict = "PASS" if bad == 0 else "FAIL"
    print("  [blosum62]    %d cells vs Biopython, %d differ -> %s" % (n, bad, verdict))
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
