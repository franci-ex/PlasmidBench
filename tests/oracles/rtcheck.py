#!/usr/bin/env python3
"""Read PlasmidBench's exported GenBank with Biopython and compare it with
the original file, also read with Biopython.

Usage:  python3 -I rtcheck.py <exported dir> <corpus dir> [<corpus dir> ...]

The exported file is named <original>.<record index>.gb, so each one can be
matched back to the record it came from. Both sides of the comparison are
read by Bio.SeqIO, which means PlasmidBench's reader takes no part in
judging PlasmidBench's writer: a writer that drops a feature and a reader
that cannot see it would otherwise agree with each other perfectly.

Checked: the sequence, the topology, and every feature's type and extracted
bases. A feature the original file declared and the export lost is a
failure, with the feature named.
"""
import json
import os
import sys
import warnings

from Bio import SeqIO

warnings.simplefilter("ignore")

EXPORTED, CORPUS = sys.argv[1], sys.argv[2:]
ok = bad = notes = 0
problems = []


def crossing(rec, L):
    """Features written as a join that runs to the end of the molecule and on
    from base 1 -- only possible on a circle."""
    out = []
    for f in rec.features:
        if f.location is None or len(f.location.parts) < 2:
            continue
        lo = min(int(p.start) for p in f.location.parts)
        hi = max(int(p.end) for p in f.location.parts)
        if lo == 0 and hi == L:
            out.append((f.type, len(f.extract(rec.seq))))
    return out


def original(name):
    base, idx = name.rsplit(".", 2)[0], int(name.rsplit(".", 2)[1])
    for d in CORPUS:
        p = os.path.join(d, base)
        if os.path.exists(p):
            fmt = ("snapgene" if base.endswith(".dna")
                   else "embl" if base.endswith(".embl")
                   else "fasta" if base.endswith((".fasta", ".fa", ".faa", ".fsa"))
                   else "genbank")
            try:
                return list(SeqIO.parse(p, fmt))[idx], base
            except Exception:
                return None, base
    return None, base


for name in sorted(os.listdir(EXPORTED)):
    if not name.endswith(".gb"):
        continue
    try:
        got = SeqIO.read(os.path.join(EXPORTED, name), "genbank")
    except Exception as e:
        bad += 1
        problems.append(f"  FAIL {name}: Biopython cannot read the export ({e})")
        continue
    want, base = original(name)
    if want is None:
        notes += 1
        continue

    try:
        wseq = str(want.seq).upper()
    except Exception:
        notes += 1
        continue
    if str(got.seq).upper() != wseq:
        bad += 1
        problems.append(
            f"  FAIL {base}: exported sequence is {len(got.seq)} bp, original {len(wseq)}")
        continue
    ok += 1

    wtopo = str(want.annotations.get("topology", "linear")).lower()
    gtopo = str(got.annotations.get("topology", "linear")).lower()
    if gtopo != wtopo:
        # Addgene and SnapGene write "linear" for plasmids that carry a
        # feature crossing the origin; PlasmidBench concludes circular and
        # writes that, which is a correction it makes deliberately.
        # Read as linear, such a feature is impossible and Biopython drops
        # it, so there is nothing in the original to look for -- the evidence
        # is in the export, which carries it as a join across the origin.
        wrapped = bool(crossing(got, len(wseq)))
        if gtopo == "circular" and wtopo == "linear" and wrapped:
            notes += 1
            problems.append(
                f"  note {base}: exported as circular; the original says linear but "
                f"annotates a feature across the origin")
        else:
            bad += 1
            problems.append(f"  FAIL {base}: exported topology {gtopo}, original {wtopo}")
    else:
        ok += 1

    # every feature the original declared, by the bases it covers
    def sig(rec):
        out = []
        for f in rec.features:
            if f.location is None:
                continue
            try:
                out.append((f.type, str(f.extract(rec.seq)).upper()))
            except Exception:
                pass
        return sorted(out)

    ws, gs = sig(want), sig(got)
    if ws != gs:
        lost = [t for t in ws if t not in gs]
        gained = [t for t in gs if t not in ws]
        # GenBank expects a source feature spanning the whole molecule, and
        # the formats PlasmidBench reads do not all carry one -- a SnapGene
        # .dna never does. Supplying it on export is correct, so one gained
        # source covering the whole sequence is not a lost or altered
        # feature. Anything else gained is.
        gained = [t for t in gained
                  if not (t[0] == "source" and len(t[1]) == len(wseq))]
        # The same Addgene correction seen in the topology check: a feature
        # that crosses the origin is unreadable in a linear record, so
        # Biopython has none to offer and the export appears to invent one.
        # It is reported with its length, which is what shows it is the real
        # feature and not an artefact -- a 304 bp enhancer, not a plasmid.
        if gained and str(want.annotations.get("topology", "linear")).lower() == "linear":
            cross = {t for t, _ in crossing(got, len(wseq))}
            rescued = [t for t in gained if t[0] in cross]
            if rescued:
                gained = [t for t in gained if t not in rescued]
                for t, b in rescued:
                    notes += 1
                    problems.append(
                        f"  note {base}: {t} of {len(b)} bp crosses the origin; the "
                        f"original calls itself linear so Biopython cannot hold it")
        if lost or gained:
            bad += 1
            detail = ""
            if lost:
                detail = f"; lost {lost[0][0]} of {len(lost[0][1])} bp"
            if gained:
                detail += f"; gained {gained[0][0]} of {len(gained[0][1])} bp"
            problems.append(
                f"  FAIL {base}: {len(ws)} features in, {len(gs)} out{detail}")
        else:
            ok += 1
    else:
        ok += 1

print(f"  {ok} agreed, {bad} differed, {notes} noted over "
      f"{len(os.listdir(EXPORTED))} exported records")
for p in problems[:25]:
    print(p)
sys.exit(1 if bad else 0)
