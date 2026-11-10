#!/usr/bin/env python3
"""Check PlasmidBench's SHA-256 against Python's hashlib.

Usage:  python3 -I cmpsha.py <dumpsha output.json>

hashlib is the oracle. The inputs include the padding boundaries -- 55, 56,
57, 63, 64 and 65 bytes -- because that is where a hand-written SHA-256 goes
wrong, and multi-byte UTF-8 including a surrogate pair, because that is
where the string-to-bytes step goes wrong.

The three seqChecksum rows are checked against the normalisation they claim:
cleaned, uppercased, U rewritten as T, prefixed by the topology.
"""
import hashlib
import json
import re
import sys

rows = json.load(open(sys.argv[1]))
ok = bad = 0
problems = []

for r in rows:
    if r.get("normalisedFrom"):
        # work out what should have been hashed, from the documented rule
        raw = r["text"]
        clean = re.sub(r"[^ACGTURYSWKMBDHVN]", "", raw.upper()).replace("U", "T")
        topo = "circular" if "circular" in r["label"] else "linear"
        want = "sha256:" + hashlib.sha256((topo + ":" + clean).encode()).hexdigest()
        got = r["sha"]
        label = r["label"]
    else:
        data = (r["repeat"]["ch"] * r["repeat"]["n"]) if r["repeat"] else r["text"]
        want = hashlib.sha256(data.encode("utf-8")).hexdigest()
        got = r["sha"]
        label = f"{len(data)} byte(s)"
    if got == want:
        ok += 1
    else:
        bad += 1
        problems.append(f"  FAIL {label}: {got} vs hashlib {want}")

print(f"  {ok} agreed, {bad} differed over {len(rows)} inputs "
      f"(including the 55/56/57/63/64/65-byte padding boundaries)")
for p in problems[:10]:
    print(p)
sys.exit(1 if bad else 0)
