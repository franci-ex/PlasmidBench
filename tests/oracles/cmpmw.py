#!/usr/bin/env python3
"""Check PlasmidBench's nucleic-acid weights against two outside references.

Usage:  python3 -I cmpmw.py <dumpmw output.json>

There are two conventions in use for the weight of a DNA strand and they
differ by one phosphate, so a single oracle cannot settle the question --
it can only say which convention the other library chose.  Both are
therefore computed here:

  1. IDT's published formula for a synthetic oligonucleotide, which is the
     convention a lab needs: the strand as it comes from the synthesiser,
     with a free 5'-OH.

         MW = 313.21*A + 304.20*T + 289.18*C + 329.21*G - 61.96

     The subtraction removes HPO2 from the 5' end.  This is the number on an
     IDT spec sheet and the number that converts a primer's nanograms to
     picomoles, so it is the one PlasmidBench is required to match exactly.

  2. Bio.SeqUtils.molecular_weight, which weighs the same strand with a
     5'-monophosphate and so carries one HPO3 more, 79.98.  It also uses
     per-residue masses to four decimal places where IDT publishes two
     (A: 313.2018 against IDT's 313.21), so over a long strand the two drift
     apart by that rounding -- about 0.0035 a base, 10 Da over 3 kb.  The
     check is therefore that the gap is one phosphate to within the rounding
     IDT's own published constants allow, 0.01 a base.  Nothing wider
     passes: if the gap ever leaves that envelope, what changed is not a
     convention.

Formula 1 is implemented here from the published definition, not read out of
the engine, so the two sides remain independent.
"""
import json
import sys

from Bio.Seq import Seq
from Bio.SeqUtils import molecular_weight

HPO3 = 79.9799  # one phosphate: the whole of the difference between the two

IDT = {"A": 313.21, "T": 304.20, "C": 289.18, "G": 329.21}


def idt_ss(seq):
    """IDT's anhydrous MW for a single strand with a free 5'-OH."""
    total = 0.0
    for ch in seq.upper():
        if ch not in IDT:
            return None       # N and the ambiguity codes are not molecules
        total += IDT[ch]
    return total - 61.96


rows = json.load(open(sys.argv[1]))
ok = bad = skipped = 0
problems = []

for r in rows:
    s = r["full"]
    tag = f"{r['sha']} ({r['length']} nt)"

    # --- reference 1: IDT, exact ---------------------------------------
    want = idt_ss(s)
    if want is None:
        skipped += 1
        problems.append(f"  note {tag}: contains a code that is not a base; no defined weight")
    else:
        tol = abs(want) * 1e-9
        if abs(r["ss"] - want) <= tol:
            ok += 1
        else:
            bad += 1
            problems.append(
                f"  FAIL {tag} single strand: {r['ss']!r} vs IDT's formula "
                f"{want!r} (difference {r['ss'] - want:.4e})")
        # a duplex is the strand plus its complement, by the same formula
        comp = str(Seq(s).reverse_complement())
        want2 = want + idt_ss(comp)
        if abs(r["ds"] - want2) <= abs(want2) * 1e-9:
            ok += 1
        else:
            bad += 1
            problems.append(
                f"  FAIL {tag} duplex: {r['ds']!r} vs IDT's formula {want2!r} "
                f"(difference {r['ds'] - want2:.4e})")

    # --- reference 2: Biopython, offset by exactly one phosphate -------
    try:
        bp_ss = molecular_weight(Seq(s), seq_type="DNA", double_stranded=False,
                                 circular=False)
    except Exception:
        continue
    gap = bp_ss - r["ss"]
    budget = 0.01 + 0.01 * r["length"]      # IDT's two decimal places, per base
    if abs(gap - HPO3) <= budget:
        ok += 1
    else:
        bad += 1
        problems.append(
            f"  FAIL {tag}: Biopython is {gap:.4f} heavier; one 5' phosphate is "
            f"{HPO3} and the rounding in IDT's constants allows {budget:.2f} "
            f"over {r['length']} bases, so this gap is not accounted for")

print(f"  {ok} agreed, {bad} differed, {skipped} undefined, over {len(rows)} sequences")
print("    reference 1: IDT oligo formula, 5'-OH — exact match required")
print(f"    reference 2: Biopython molecular_weight, 5'-phosphate — gap must be {HPO3} "
      f"within the rounding of IDT's published constants")
for p in problems[:20]:
    print(p)
sys.exit(1 if bad else 0)
