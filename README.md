# PlasmidBench

A molecular cloning workbench that runs from **one HTML file** in any
browser. No installation, no server, no account. Open it from a lab machine,
a memory stick or a static web server and it works, offline.

It designs primers and assemblies, simulates PCR and digests, reads and
writes the formats a lab actually has, interprets sequencing traces — and,
unusually, knows what *your* freezer holds: it searches your own oligo
inventory for tubes that can already read a region, and plans expression
against the strains and plasmids you have.

```sh
git clone https://github.com/<owner>/PlasmidBench.git
cd PlasmidBench
node build.js            # writes plasmidbench.html
open plasmidbench.html   # or just double-click it
```

Nothing else to install. `node` is needed only to join the source files;
the file it produces has no dependencies at all.

---

## What it does

| | |
|---|---|
| **Reads** | GenBank, EMBL, FASTA, SnapGene `.dna`, ABI `.ab1` traces, GFF3, four JSON dialects (Benchling, Open Vector Format, openVectorEditor/jbei, its own), NCBI INSDSeq and TinySeq XML |
| **Writes** | GenBank, EMBL, FASTA, GFF3, JSON, SVG and PNG figures, `.xlsx` |
| **Designs** | restriction cloning, seamless assembly (In-Fusion / Gibson / NEBuilder), Gateway BP and LR, TOPO/TA, site-directed mutagenesis, overlap extension, silent mutations, codon optimisation, back-translation |
| **Simulates** | PCR with tailed primers, digests with 607 REBASE enzymes, agarose gels, ligation and junction checking |
| **Interprets** | Sanger reads against a construct — what is vector, what is new, what differs, in which frame, with the ORFs |
| **Knows your lab** | your oligo inventory, plasmid collection and competent-cell stock, searched and planned against |
| **Annotates** | exact and **tolerant** feature search: variants and cloning fragments found and reported with their identity |
| **Connects** | eLabFTW: reads records and their files, adds results back, never replaces anything |

---

## Why one file

Most scientific web software needs installing, a backend, or a database.
This needs a browser. That buys four things worth having:

- **it runs offline**, which is what a tissue-culture room with no network is;
- **it is one file to hand to someone**, including a student on their first day;
- **it needs no administrator**, which is what a shared lab machine needs;
- **a release is a file**, so the version that produced a result can be kept.

The *source* is modular — 30 files joined by `build.js` — because a 1.4 MB
monolith would be impossible to maintain. Monolithic distribution is an
advantage; monolithic development would be a liability.

---

## Is it right?

That is the question that matters, and the answer is meant to be checkable
rather than promised.

```sh
./validate-all.sh
```

Four layers, about 25 minutes:

| layer | what it is | size |
|---|---|---|
| **types** | `tsc --checkJs` over the engine, from JSDoc | **0 errors** |
| **engine** | 17 unit suites | **8,892 assertions** |
| **oracles** | the engine against Biopython 1.88, REBASE, IDT's published formula, Python's `hashlib` | **11,000+ comparisons** |
| **browser** | 25 Playwright suites driving the built file | **850+ assertions** |

The oracles are the part that makes the rest worth anything. Each is **two
halves that share no code**: a JavaScript half that reports what the engine
says and does not compare, and a Python half that computes the expected
answer from an independent library and compares. Neither half can be tuned
to agree with the other.

| what | against | result |
|---|---|---|
| 607 restriction enzymes: site, both cut offsets, every cut position on five sequences, linear and circular | `Bio.Restriction` / REBASE | 7,891 agreed, 0 differed |
| supplier letters and company names | REBASE's supplier table | 622 agreed, 0 differed |
| protein descriptors (MW, pI, GRAVY, aromaticity, instability) on 310 real proteins | `Bio.SeqUtils.ProtParam` | 892 agreed, 0 differed |
| oligo and duplex weight | IDT's published formula **and** `Bio.SeqUtils` | 33 agreed, 0 differed |
| GenBank reader, 45 real and malformed files | `Bio.SeqIO` | 170 agreed, 0 differed |
| EMBL reader, 17 files | `Bio.SeqIO` | 168 agreed, 0 differed |
| SnapGene `.dna` and ABI `.ab1`, 16 binary files | `Bio.SeqIO` | 38 agreed, 0 differed |
| GFF3, Benchling, openVectorEditor, INSDSeq XML | Biopython-written fixtures | 1,646 checks, 0 differ |
| GenBank export, read back | `Bio.SeqIO` on both sides | 128 agreed, 0 differed |
| SHA-256 checksums | Python `hashlib` | 20 agreed, 0 differed |
| annotation identities | `Bio.Align.PairwiseAligner` | 162 agreed, 0 differed |
| BLOSUM62 | `Bio.Align.substitution_matrices` | 576 cells, 0 differ |

Every divergence from an oracle is **named and justified in the comparison
script**, never silenced. There are four, and two of them are places where
this project is more correct than Biopython — Biopython flattens a SnapGene
feature that crosses the origin to the whole plasmid, which turns
pFA-KanMX4's 19 bp SP6 promoter into a 3,941 bp one.

`tests/README.md` explains the one command that fetches the corpus.

---

## What it is not

Stated plainly, because a tool that overstates itself is dangerous in a lab.

- **Not a record system.** It holds no authoritative records and should not.
  Attribution, timestamps and an audit trail belong to eLabFTW, which this
  reads from and *adds to* — never replacing or patching anything. A
  construct does not become verified lab stock because a design tool
  produced it. See `docs/elabftw.md`.
- **Not a replacement for SnapGene.** SnapGene's editors are far more
  developed and it ships thousands of curated vectors. This ships one demo
  construct and 27 feature sequences.
- **Not a cloning-provenance system.** OpenCloning does that, with a Python
  backend and reproducible workflows. The four simulators here are enough;
  they need tests, not successors.
- **Not a database for 100,000 records.** It holds everything in the
  browser. Past a few thousand records that is the wrong architecture, and
  eLabFTW or OpenCloningDB is the right one.
- **Not BLAST.** The tolerant annotator reports identity, coverage and the
  chance of its seed being coincidence. It does not report an E-value,
  because an E-value needs a database model it does not have.

---

## Making it yours

A lab's defaults are not the author's. **Tools → Lab profile** holds the
lab name, the bench conditions every melting temperature is calculated at,
which suppliers you can actually buy from, what your sequencing facility
really gives you, your spreadsheet layouts, and your eLabFTW address — and
exports the lot as one file another lab can import.

It carries **no credentials and no records**. That is enforced by pattern in
`src/profile.js` and checked by a test that plants an API key and an oligo
inventory and then searches the exported file for them.

These are not decorative. Setting your suppliers filters the enzyme list to
what you can order; setting your facility's read length changes what the
sequencing planner draws.

---

## Documentation

| | |
|---|---|
| `docs/manual.md` | using it, with the things that are easy to get wrong |
| `docs/developer-guide.md` | how to change it, how to run the checks, how to add a format |
| `docs/architecture.md` | the 30 source files, the build, and why it is shaped this way |
| `docs/elabftw.md` | the integration, and the CORS question answered |
| `docs/limitations.md` | what it does not do, and what is known to be approximate |
| `AGENTS.md` | the conventions. Read before changing anything |
| `NOTICE.md` | where every data table came from and what its terms are |
| `docs/evaluation.md` | the external evaluation this version was built to answer |

---

## Licence and credit

MIT (`LICENSE`). The data tables came from REBASE, NCBI, Kazusa, IDT and
published papers, and have their own terms: **read `NOTICE.md` before
redistributing the built file or publishing with it.**

Parts of this codebase were written with the help of Anthropic's Claude. The
scientific choices and the validation strategy are the author's. Every
algorithm in the engine is explained beside it and every numerical claim has
an outside oracle — which is what makes AI-assisted scientific code
reviewable at all.

If you use PlasmidBench in work you publish, please cite REBASE and
Biopython as well; this would not be checkable without them.
