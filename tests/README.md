# The checks

Four layers. `./validate-all.sh` runs them all in about 25 minutes; each can
be run on its own.

```sh
./validate-all.sh types     # tsc --checkJs over the engine. Must be 0.
./validate-all.sh engine    # 17 unit suites, ~8,900 assertions, ~1 minute
./validate-all.sh oracles   # the outside comparisons, ~15 minutes
./validate-all.sh browser   # 25 Playwright suites, ~8 minutes
```

## What you need

```sh
npm install --no-save playwright@1.48.2 typescript@5.6.3 @types/node@22
npx playwright install chromium
pip install biopython==1.88
```

Biopython is **pinned** because the expected values come from that version.
An unpinned oracle is not an oracle.

### The corpus

The file-format oracles read real — and deliberately malformed — GenBank,
EMBL, SnapGene and ABI files from Biopython's own test directory. It is not
vendored here: it is Biopython's, and fetching it keeps the provenance
obvious.

```sh
git clone --filter=blob:none --no-checkout --depth 1 \
  https://github.com/biopython/biopython.git bp
cd bp
git sparse-checkout init --cone
git sparse-checkout set Tests/GenBank Tests/EMBL Tests/SnapGene \
  Tests/Abi Tests/Fasta Tests/Nexus
git checkout
```

Then point the suite at it:

```sh
REAL=/path/to/bp/Tests ./validate-all.sh oracles
```

Without it the format oracles say they are skipping and the rest still runs.
**Fetch it before trusting a green run**: three real bugs in the readers were
found only by files in there, and two suites had been quietly passing on an
empty corpus for months.

## How the oracles are built, and why it matters

Each is **two halves that share no code**:

- a **JavaScript half** (`tests/oracles/dump*.js`) that reports what the
  engine says, and does not compare anything;
- a **Python half** (`tests/oracles/cmp*.py`) that computes the expected
  answer from an independent library and compares.

Neither half can be tuned to agree with the other without somebody noticing,
which is the whole point: a test written beside the code it tests can only
show that the code agrees with itself. That is true of any codebase and
especially of one written with AI assistance.

| oracle | against |
|---|---|
| `test-rebase.py` | `Bio.Restriction` — 607 enzymes, their sites, both cut offsets, and every cut position on five sequences read linearly and circularly |
| `cmpsuppl.py` | REBASE's supplier table |
| `cmpgb.py`, `cmpembl.py` | `Bio.SeqIO` on real and malformed files, comparing sequence, topology, feature count and **the bases each feature extracts** |
| `cmpreal.py` | `Bio.SeqIO` on SnapGene `.dna` and ABI `.ab1` — binary formats, where a wrong offset gives a sequence that still looks like DNA |
| `cmpprops.py` | `Bio.SeqUtils.ProtParam` — only descriptors with one published definition |
| `cmpmw.py` | IDT's published oligo formula **and** `Bio.SeqUtils.molecular_weight`, which uses the other end convention |
| `rtcheck.py` | the GenBank exporter, read back by `Bio.SeqIO` — so the engine's reader takes no part in judging the engine's writer |
| `cmpsha.py` | Python's `hashlib`, including the 55/56/57/63/64/65-byte padding boundaries |
| `cmpanno.py` | `Bio.Align.PairwiseAligner` — the identity of every annotation, against the exact span the matcher named |
| `cmpb62.py` | `Bio.Align.substitution_matrices` |
| `cmpfmt.js` | GFF3, Benchling, openVectorEditor and INSDSeq XML fixtures written by Biopython |

### Divergences are named, never silenced

Where the two libraries legitimately disagree, the comparison script holds
the case **with its reason**, and reports it as a note rather than a pass.
There are four. Two of them are places where this project is more correct:

- Biopython flattens a SnapGene feature that crosses the origin to the whole
  molecule, turning pFA-KanMX4's 19 bp SP6 promoter into a 3,941 bp one. The
  note reports the length, which is the evidence.
- Addgene and SnapGene write `linear` on the LOCUS line of plasmid exports
  and then annotate a feature that crosses the origin. This concludes the
  molecule is circular and records that it did; Biopython takes the LOCUS
  line at its word and drops the feature.

If you add a divergence, give it a reason. "It fails" is not one.

## The unit suites

`tests/test*.js`. Node, no framework, one function:

```js
ok('what should be true', condition, somethingUsefulWhenItIsNot);
```

The third argument matters more than it looks: a failure that prints the
actual value is a failure you can act on.

Several carry comments saying a previous expectation was wrong, what
replaced it, and what settled it. Those comments are evidence and should not
be tidied away.

## The browser suites

`tests/uitest*.js`. Playwright drives the built `plasmidbench.html` from a
`file://` URL, which is also how a lab machine opens it. They test what a
person sees, not internals:

- `uitest-escape.js` — a construct whose every field is markup, opened in
  every view and every tool. Nineteen assertions covering four hundred
  interpolation sites, and it keeps working as the interface changes.
- `uitest-store.js` — a construct survives a reload; a 400 kb construct does
  too; the old localStorage library is migrated; a browser that refuses to
  keep data says so.
- `uitest-elab.js` — the whole eLabFTW cycle against `tests/mock-elab.js`,
  in each of its CORS modes, including the assertion that no `PATCH`, no
  `DELETE` and no POST to an existing upload ever leaves the page.
- `uitest-profile.js` — a lab profile exports and imports, changes what the
  page calculates, and carries no credentials and no records.

### Running one

```sh
node tests/uitest-infusion.js
TMPDIR=/tmp/pb node tests/uitest-store.js   # if /tmp is small
```

Build first: the browser suites test the built file, not the source.

```sh
node build.js && node tests/uitest.js
```

## Measurement harnesses

Some decisions in this codebase were made by measuring and are written down
with their numbers — the seed length for tolerant annotation, the recall at
each divergence, the false-positive rate. `tests/test19.js` contains the
recall measurement so that it is re-run on every validation rather than
being a number in a comment that rots.

One thing to know before writing a test that needs random sequence: **the LCG
used in the older suites has a single cycle**, so two different seeds are two
starting points on one sequence. A "random" negative control drawn from it
can contain the library drawn from it — which is exactly what happened while
measuring false positives, and made 365 genuine matches look like a defect.
Use mulberry32, as `test19.js` does, when a sequence genuinely must not
contain something.
