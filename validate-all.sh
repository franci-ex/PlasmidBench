#!/bin/bash
# ============================================================================
# Every check this project has, in one run.
#
# Three kinds, in order of how much they prove:
#
#   1. unit suites        the engine against truth worked out in the test
#   2. outside oracles    the engine against Biopython 1.88, REBASE and IDT's
#                         published oligo formula, on real files
#   3. browser suites     the built page driven through Playwright
#
# The oracles are what make the rest worth anything: a test written beside
# the code it tests can only show the code agrees with itself. Each oracle is
# two halves that do not share any code -- a JavaScript half that reports
# what the engine says, and a Python half that computes the expected answer
# from an independent library and compares. Neither half can be tuned to
# agree with the other without someone noticing.
#
# Requirements:
#   pip install biopython        (1.88 is what the expected values come from)
#   Biopython's Tests/ corpus at $REAL  -- real GenBank, EMBL, SnapGene and
#   ABI files, including the malformed ones. Fetch with:
#     git clone --filter=blob:none --no-checkout --depth 1 \
#       https://github.com/biopython/biopython.git bp
#     cd bp && git sparse-checkout init --cone
#     git sparse-checkout set Tests/GenBank Tests/EMBL Tests/SnapGene \
#       Tests/Abi Tests/Fasta Tests/Nexus && git checkout
#
# Usage:  ./validate-all.sh [types|engine|oracles|browser]
# ============================================================================
set -u
cd "$(dirname "$0")"
REAL=${REAL:-/home/claude/realfiles/bp/Tests}
SP=${SP:-${TMPDIR:-/tmp}/pb-validate}
mkdir -p "$SP"
UT=${UT:-$SP/uitmp}; mkdir -p "$UT"
WHAT=${1:-all}
fails=0

note () { printf '\n== %s ==\n' "$1"; }
run  () { printf '  %-22s ' "$1"; shift
          if "$@" >"$SP/out.txt" 2>&1; then
            sed -n '/[^[:space:]]/{p;q;}' "$SP/out.txt" | sed 's/^ *//'
          else
            fails=$((fails+1)); echo FAILED
            sed 's/^/      /' "$SP/out.txt" | head -12
          fi; }

have_corpus () { [ -d "$REAL/GenBank" ]; }

# ----------------------------------------------------------------- 0. types --
if [ "$WHAT" = all ] || [ "$WHAT" = types ]; then
note "types (tsc --checkJs, no build step, nothing emitted)"
if [ -x ./node_modules/.bin/tsc ]; then
  printf '  %-22s ' "engine + store"
  n=$(./node_modules/.bin/tsc -p tsconfig.json 2>&1 | grep -c "error TS" || true)
  if [ "$n" -eq 0 ]; then echo "0 type errors"; else
    fails=$((fails+1)); echo "$n type error(s)"
    ./node_modules/.bin/tsc -p tsconfig.json 2>&1 | grep "error TS" | head -10 | sed 's/^/      /'
  fi
else
  echo "  skipped: npm install --no-save typescript @types/node"
fi
fi

# ---------------------------------------------------------------- 1. units --
if [ "$WHAT" = all ] || [ "$WHAT" = engine ]; then
note "engine unit suites"
for t in test.js test2.js test3.js test4.js test5.js test6.js test7.js test8.js \
         test9.js test10.js test13.js test14.js test15.js test16.js test17.js test18.js \
         test19.js; do
  printf '  %-22s ' "$t"
  r=$(node "tests/$t" 2>&1 | grep -oE '[0-9]+ passed, [0-9]+ failed' | tail -1)
  [ -z "$r" ] && { r="NO RESULT"; fails=$((fails+1)); }
  case "$r" in *", 0 failed") ;; *) fails=$((fails+1));; esac
  echo "$r"
done
fi

# -------------------------------------------------------------- 2. oracles --
if [ "$WHAT" = all ] || [ "$WHAT" = oracles ]; then
note "restriction enzymes and cut positions vs Bio.Restriction (REBASE)"
node tests/oracles/test-rebase-node.js >/dev/null && run "607 enzymes" python3 -I tests/oracles/test-rebase.py
run "supplier table" python3 -I tests/oracles/cmpsuppl.py tests/oracles/rebase-mine.json

note "protein descriptors vs Bio.SeqUtils.ProtParam"
if have_corpus; then
  node tests/oracles/dumpprops.js "$SP/props.json" $REAL/GenBank/*.gb $REAL/GenBank/*.gbk \
      $REAL/EMBL/*.embl $REAL/Fasta/* 2>/dev/null
else
  node tests/oracles/dumpprops.js "$SP/props.json"
fi
run "ProtParam" python3 -I tests/oracles/cmpprops.py "$SP/props.json"

note "oligo and duplex weight vs IDT's formula and Bio.SeqUtils"
node tests/oracles/dumpmw.js > "$SP/mw.json" && run "weights" python3 -I tests/oracles/cmpmw.py "$SP/mw.json"

if have_corpus; then
  note "GenBank reader vs Bio.SeqIO, on real and malformed files"
  F=$(ls $REAL/GenBank/*.gb $REAL/GenBank/*.gbk 2>/dev/null | grep -v -E "1MRR_A|output")
  node tests/oracles/dumpgb.js $F > "$SP/gb.json" && run "GenBank" python3 -I tests/oracles/cmpgb.py "$REAL/GenBank" "$SP/gb.json"

  note "EMBL reader vs Bio.SeqIO"
  node tests/oracles/dumpembl.js $REAL/EMBL/*.embl > "$SP/embl.json" &&
    run "EMBL" python3 -I tests/oracles/cmpembl.py "$REAL/EMBL" "$SP/embl.json"

  note "SnapGene .dna and ABI .ab1 vs Bio.SeqIO"
  node tests/oracles/dumpreal.js $REAL/SnapGene/*.dna $REAL/Abi/*.ab1 > "$SP/real.json" &&
    run "binary formats" python3 -I tests/oracles/cmpreal.py "$SP/real.json" "$REAL/SnapGene" "$REAL/Abi"

  note "GFF3, Benchling, openVectorEditor and INSDSeq XML"
  rm -rf "$SP/fmt"; python3 -I tests/oracles/mkfmt.py "$REAL/GenBank" "$SP/fmt" >/dev/null 2>&1
  run "other formats" node tests/oracles/cmpfmt.js "$SP/fmt"

  note "GenBank export, read back by Biopython"
  rm -rf "$SP/rt"
  node tests/oracles/exportgb.js "$SP/rt" $F $REAL/SnapGene/*.dna >/dev/null &&
    run "round trip" python3 -I tests/oracles/rtcheck.py "$SP/rt" "$REAL/GenBank" "$REAL/SnapGene"
else
  note "file-format oracles"
  echo "  skipped: no corpus at $REAL (see the header of this script)"
fi

note "SHA-256 vs Python's hashlib"
node tests/oracles/dumpsha.js "$SP/sha.json" && run "checksums" python3 -I tests/oracles/cmpsha.py "$SP/sha.json"

note "annotation identities vs Bio.Align.PairwiseAligner"
node tests/oracles/dumpanno.js "$SP/anno.json" 2>/dev/null &&
  run "tolerant matcher" python3 -I tests/oracles/cmpanno.py "$SP/anno.json"

note "BLOSUM62 vs Bio.Align.substitution_matrices"
node tests/oracles/dumpb62.js > "$SP/b62.json" && run "BLOSUM62" sh -c 'python3 -I tests/oracles/cmpb62.py < "$1"' _ "$SP/b62.json"
fi

# ------------------------------------------------------------- 3. browser ---
if [ "$WHAT" = all ] || [ "$WHAT" = browser ]; then
note "browser suites (Playwright against the built page)"
node build.js >/dev/null
for t in uitest.js uitest-dna.js uitest-zoom.js uitest-clone.js uitest-triage.js \
         uitest-edit.js uitest-export.js uitest-infusion.js uitest-sanger.js \
         uitest-tagswap.js uitest-express.js uitest-gel.js uitest-formats.js \
         uitest-save.js uitest-multisel.js uitest-figure.js uitest-verify.js \
         uitest-anatomy.js uitest-wrap.js uitest-primertm.js uitest-store.js uitest-escape.js uitest-elab.js uitest-profile.js uitest-annotate.js; do
  printf '  %-22s ' "$t"
  r=$(TMPDIR="$UT" timeout 300 node "tests/$t" 2>&1 | grep -oE '[0-9]+ passed, [0-9]+ failed' | tail -1)
  [ -z "$r" ] && { r="NO RESULT"; fails=$((fails+1)); }
  case "$r" in *", 0 failed") ;; *) fails=$((fails+1));; esac
  echo "$r"
done
fi

printf '\n'
if [ "$fails" -eq 0 ]; then echo "everything agreed."; else echo "$fails check(s) did not pass."; fi
exit $([ "$fails" -eq 0 ] && echo 0 || echo 1)
