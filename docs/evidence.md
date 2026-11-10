# What is checked, and against what

A one-page answer to "is it right?", for someone deciding whether to trust a
number this tool gave them. The numbers below come from
`./validate-all.sh`; run it yourself rather than taking them from here.

## The layers

| layer | what it is | size |
|---|---|---|
| types | `tsc --checkJs` over the engine, from JSDoc and one `.d.ts` | **0 errors** |
| engine | 17 unit suites, Node | **8,892 assertions** |
| oracles | the engine against independent implementations | **11,000+ comparisons** |
| browser | 25 Playwright suites against the built file | **853 assertions** |

## Why the oracles are built the way they are

Each one is **two halves that share no code**:

- `tests/oracles/dump*.js` reports what the engine says, and compares
  nothing;
- `tests/oracles/cmp*.py` computes the expected answer from an independent
  library and compares.

Neither half can be tuned to agree with the other without somebody noticing.
A test written beside the code it tests can only show that the code agrees
with itself — true of any codebase, and the whole question for one written
with AI assistance.

## What each oracle covers

| what | against | result |
|---|---|---|
| 607 restriction enzymes: recognition site, both cut offsets, every cut position on five sequences, read linearly and circularly | `Bio.Restriction` / REBASE | 7,891 agreed, 0 differed |
| supplier letters and company names | REBASE's supplier table | 622 agreed, 0 differed |
| protein descriptors — MW, pI, GRAVY, aromaticity, instability — on 310 real proteins | `Bio.SeqUtils.ProtParam` | 892 agreed, 0 differed |
| oligo and duplex weight | IDT's published formula **and** `Bio.SeqUtils.molecular_weight` | 33 agreed, 0 differed |
| GenBank reader on 45 real and deliberately malformed files: sequence, topology, feature count, and the bases each feature extracts | `Bio.SeqIO` | 170 agreed, 0 differed |
| EMBL reader, 17 files | `Bio.SeqIO` | 168 agreed, 0 differed |
| SnapGene `.dna` and ABI `.ab1`, 16 binary files | `Bio.SeqIO` | 38 agreed, 0 differed |
| GFF3, Benchling, openVectorEditor and INSDSeq XML | fixtures written by Biopython | 1,646 checks, 0 differ |
| GenBank export, read back — **Biopython on both sides**, so the engine's reader takes no part in judging its writer | `Bio.SeqIO` | 128 agreed, 0 differed |
| SHA-256, including the 55/56/57/63/64/65-byte padding boundaries and multi-byte UTF-8 | Python `hashlib` | 20 agreed, 0 differed |
| the identity of every annotation, against the exact span the matcher named | `Bio.Align.PairwiseAligner` | 162 agreed, 0 differed |
| BLOSUM62, cell by cell | `Bio.Align.substitution_matrices` | 576 cells, 0 differ |

## Where it disagrees with the oracle, and why

Four cases. Each is held in the comparison script **with its reason** and
reported as a note, never silenced. Two are places where this project is
more correct:

1. **A SnapGene feature that crosses the origin.** Biopython 1.88 flattens
   it to `0..length` — the whole plasmid — which makes pFA-KanMX4's SP6
   promoter 3,941 bp instead of 19. PlasmidBench keeps the wrap. The note
   reports the length, which is the evidence: a promoter is not the size of
   its plasmid.
2. **Addgene's `linear` plasmids.** Addgene and SnapGene write `linear` on
   the LOCUS line and then annotate `enhancer 8569..276` on an 8,596 bp
   circle. Biopython takes the LOCUS line at its word and the feature is
   impossible. PlasmidBench concludes the molecule is circular and records
   that it did.
3. **Four protein patent records in `patents.embl`** whose FT block precedes
   the SQ line. Biopython returns no features; PlasmidBench reads the
   `source 1..358` that is plainly written there. Verified by reading the
   file by hand.
4. **Two conventions for oligo weight**, differing by one 5' phosphate. IDT's
   convention is exact here because it is the one that converts a primer's
   nanograms to picomoles; Biopython's is checked to be that one phosphate
   away, within the rounding IDT's own two-decimal constants allow.

## Measured, not assumed

Several decisions are recorded with their measurements, because the obvious
guess was wrong:

- **A tailed primer's annealing Tm against the whole primer's.** The whole
  primer is *nearly* always higher — measured at 99.6% of random primers,
  not 100%, so the code does not say "always".
- **Tolerant annotation**, over 200 trials per point on a 1,200-entry
  library: 100% / 99% / 94% / 88–92% recall at 0 / 2 / 5 / 10% divergence,
  every hit lying over the planted feature, **0.03 false positives per 5 kb**.
  7.6 ms per annotation at a 14-base seed, 199 ms at 12.
- **The cost of guessing an arm boundary** rather than knowing it: 0.7 °C,
  always low, against the 15.2 °C error of using the whole primer.

## The corpus

The format oracles read real and deliberately malformed files from
Biopython's own `Tests/` directory, fetched at validation time rather than
vendored. Three live defects in the readers were found only by files in
there, and two suites had been quietly passing on an empty corpus for
months — so **fetch it before trusting a green run.** One command, in
`tests/README.md`.

## What is not checked this way

Said plainly:

- **No oracle for the cloning simulators' products.** Expected products are
  worked out by hand in the tests. OpenCloning could serve as one and has
  not been used yet.
- **No oracle for the expression heuristics.** They are published heuristics
  combined by judgement, and there is no reference implementation to compare
  against. `docs/limitations.md` says what they are and are not.
- **No oracle for the gel mobility model.** It is a fit for recognising a
  band, not for measuring one, and is documented as such.
- **No ground-truth annotation benchmark.** This is the gap worth closing,
  and §12.1 of the evaluation argues it is also the publishable piece of
  work: a curated set of reference plasmids with verified feature maps,
  scored by precision and recall. It needs hand-checking against the
  literature, which is the part no amount of code produces.
