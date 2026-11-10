# Changelog

## v0.4.0 — November 2026

The release that answers the external evaluation (`docs/evaluation.md`). The
theme of it is evidence rather than features: the project had 25 tools and
no way for anyone else to check that any of them was right.

### Validation — the part that matters

A validation harness in four layers, `./validate-all.sh`:

- **types** — `tsc --checkJs` over the engine from JSDoc and one `.d.ts`.
  **Zero errors.** No TypeScript to compile; `src/` stays plain JavaScript.
- **engine** — 17 unit suites, **8,892 assertions**.
- **oracles** — **11,000+ comparisons** against Biopython 1.88, REBASE, IDT's
  published oligo formula, Python's `hashlib` and `Bio.Align`. Each oracle is
  two halves that share no code, so neither can be tuned to agree with the
  other.
- **browser** — 25 Playwright suites, **850+ assertions**, driving the built
  file.

Every divergence from an oracle is named and justified in the comparison
script. Two of the four are places where this project is more correct than
Biopython.

### Defects it found

Each of these was live, and none would have been found without an outside
reference:

- **Protein records were silently mangled.** E, F, I, L, P and Q are not
  nucleotide codes, so DNA cleaning deleted two residues in five and left a
  string that still read as DNA. A 1,841-residue record became 1,244
  "bases".
- **An Addgene plasmid lost an annotation.** Addgene and SnapGene write
  `linear` on plasmid exports and then annotate features that cross the
  origin. Taken at its word the location is impossible and the feature was
  dropped. A backwards span cannot exist on a linear molecule, so it is now
  the stronger evidence; the topology is corrected and the correction is
  recorded.
- **Records with no `ORIGIN` lost every feature** — every location was
  range-checked against length zero.
- **`require('./engine7.js')` did not load the engine.** Modules 3 to 6 were
  not in the chain, so "importable into Node" held only if you required all
  files in order.
- **A ragged stock-list header sent you to the wrong freezer box.** A
  nine-column table with an eight-column header had the sequence land right
  by luck while the location read from the room-number column.
- **A long `/note` or `/label` was permanently truncated by reading the
  file** — the reader cut it to 80 characters and the shortened string went
  into exports.
- **A `/note` with no spaces came back altered from a round trip**, because
  the writer broke the line mid-word and a reader rejoins with a space.
- **`E.dimer` returned two different shapes** depending on whether it found
  complementarity.
- **The `.ab1` reader threw on a damaged file** — two tags come back as byte
  arrays and were treated as strings.
- **Two copies of `layoutLabels` had diverged**: the live one fans a crowded
  side into columns, the dead one did not.
- **The alignment view did not escape its row labels**, which are record
  names.

### Lab deployment (evaluation §4.1, §8.1)

- **The claude.ai runtime dependency is gone.** `src/store.js` makes
  **IndexedDB** the system of record. The 200,000-base ceiling and the
  5–10 MB localStorage budget are gone with it; a 400 kb construct now
  survives a reload.
- Settings migrate from the old `pb.*` keys on first load and still read
  synchronously, from an in-memory cache written through to IndexedDB.
- A host-side store is still read and mirrored to where one exists, so
  nothing saved that way is stranded. One setting turns it off.
- **The status line no longer claims an account it has not got**, and says
  plainly when the browser will not keep data at all.
- **Saving a file downloads it** when there is no host, instead of explaining
  that saving is unavailable.

### eLabFTW (§7)

- `src/elab.js` and the panel, under one rule: **read freely, write only by
  adding.** No `PATCH`, no `DELETE`, no replacing an attachment — checked
  against the wire in `tests/uitest-elab.js`.
- **The CORS question is a button.** Five steps, reported separately, because
  "blocked by CORS" and "your key was rejected" have different remedies and
  look identical from inside a browser. On a default install it names the two
  missing headers and points at the file route.
- `tests/mock-elab.js` reproduces every case — full CORS, default install, no
  CORS, read-only key, rejected key, corrupted download — so the integration
  is testable without an instance.
- **SHA-256 in the engine** (because `crypto.subtle` does not exist on a
  `file://` page), checked against `hashlib`. A download that does not match
  its checksum is refused, not opened.

### Tailoring for a lab

- `src/profile.js` and a panel: one exportable bundle of everything a lab has
  to set. **No credentials and no records**, enforced by pattern and checked
  by a test that plants both and searches the exported file.
- Suppliers now filter the enzyme list; the sequencing planner opens on your
  facility's real read length. Applying a profile takes effect immediately.

### Tolerant annotation (§10)

- `src/engine8.js`: k-mer index, two-hit rule, X-drop extension, banded
  Smith–Waterman, Aho–Corasick for the exact pass.
- Measured: **100% / 99% / 94% / 88–92%** recall at 0 / 2 / 5 / 10%
  divergence, 0.03 false positives per 5 kb, 7.6–199 ms per annotation
  against a 1,200-entry library.
- Library import from **ApE**, **FPbase** and its own format, every entry
  keeping its source and licence. No Addgene crawling.
- Reports identity, coverage and whether a match is partial — a fragment is a
  finding, not something to discard.
- Identities checked against `Bio.Align.PairwiseAligner`: 162 comparisons, no
  disagreements.

### Repository (§5, §9, §13)

`src/`, `tests/`, `tests/oracles/`, `docs/`; MIT licence; `NOTICE.md` naming
every data table's origin and terms; `AGENTS.md`; `CONTRIBUTING.md`; a user
manual, developer guide, architecture note and limitations list; GitHub
Actions running all four layers, and a release workflow that attaches the
built file.

### Measurements that changed a decision

Recorded because the obvious guess was wrong each time:

- The whole primer's melting temperature is **not** always higher than its
  annealing portion's — 0.4% of random primers invert it.
- The seed length for annotation was chosen from recall and false-positive
  rates across four values.
- A local aligner with ±1 scoring always over-extends; the mismatch penalty
  is now **derived from the identity asked for**.
- The extension and the alignment need **different** penalties — found by
  setting them equal and watching recall fall.
- **The LCG in the older test suites has a single cycle**, so two seeds are
  two starting points on one sequence. A "random" negative control contained
  365 of 1,200 library entries, which made 365 correct matches look like
  false positives.

### Also

- Reusable backbone primer pairs across a series of constructs, keeping their
  own names, with a warning when the placement and the chosen pair disagree.
- The two menu entries both called "Edit sequence" are now named for what
  they do.
- Labels are kept whole in the record and shortened only where drawn.

## Earlier

Developed from 2025 as a single-file workbench, reaching 25 tools, 10 views,
readers for a dozen formats and an `.xlsx` writer, with no tests. The
evaluation in `docs/evaluation.md` was written against that state.
