# PlasmidBench: Project Evaluation and Development Roadmap

Measurements taken from `PlasmidBench.v26.html` (22,366 lines, 1.22 MB, 16
script blocks), built from 17 JavaScript source files — seven engine modules, a
spreadsheet writer, eight interface modules, and a read-interpretation module —
plus an HTML/CSS shell, joined by `build.js`.

---

## 1. Overall Assessment

PlasmidBench is a capable and well-conceived piece of scientific software. It
contains 25 tools across six groups, 10 main views, and readers for GenBank,
EMBL, FASTA, SnapGene `.dna`, `.ab1` traces, GFF3, four JSON dialects
(Benchling, Open Vector Format, jbei/openVectorEditor, and its own) and NCBI
INSDSeq/TinySeq XML. It bundles 612 REBASE restriction enzymes with 15 supplier
codes and includes a from-scratch `.xlsx` writer.

The architecture is sound. The 9,421-line engine contains **zero** references
to `document`, `innerHTML` or `addEventListener`, and each engine module exports
under CommonJS. Scientific logic is cleanly separated from the interface and is
already importable into Node. The interface reaches the engine only through its
public namespace (297 call sites).

The breadth of functionality and the quality of the separation are both
significant achievements.

The priority now is **evidence** rather than structure: the project has no
tests, no public development history, no independent validation, and one
runtime dependency that prevents laboratory deployment. Addressing these four
items is the shortest path to software the laboratory can rely on and others
can maintain.

---

## 2. The Single-File Distribution Is a Genuine Strength

Most scientific web applications require installation, a backend server, or a
database. PlasmidBench runs from one HTML file in any browser.

1. **Portability** — one file to download or share, no installation.
2. **Offline operation** — core analysis and design work without a network.
3. **Accessibility** — no administrator privileges or server infrastructure.
4. **Low maintenance for users** — distribute as a versioned HTML release.
5. **Educational value** — low barrier for students exploring the calculations.

The combination of modular source and single-file output is the right
arrangement and should be preserved. Monolithic *distribution* is an advantage;
monolithic *development* would be a liability, and the project avoids it.

Modular source also has a practical benefit for AI-assisted work: editing a
1.22 MB built file costs roughly fifty times the tokens of editing a 300-line
module, on every change.

---

## 3. Position Relative to SnapGene and OpenCloning

### SnapGene

SnapGene is a mature commercial application with extensively developed sequence
and map editors, comprehensive cloning tools, and broad format support.
PlasmidBench implements a meaningful subset.

Two gaps account for most of SnapGene's practical advantage:

| | PlasmidBench | SnapGene |
|---|---|---|
| Built-in annotated vectors | 1 demo construct | thousands, curated |
| Features available for auto-annotation | 27 (`FEATURE_LIB`) | thousands, curated |
| Writes SnapGene `.dna` | no (reads only) | native |

Section 7 addresses the annotation gap, and the recommendation is not to build
a feature library.

SnapGene is a useful benchmark for functionality and usability, but it should
not be the target, and the laboratory holds no licence, so it cannot serve as a
validation reference (section 6).

### OpenCloning

OpenCloning is an open-source cloning platform with reproducible workflows,
construct provenance, a Python backend, and an established eLabFTW integration.

Its cloning simulators overlap PlasmidBench's, and it already has the eLabFTW
integration the laboratory wants. **Further development effort should not go
into PlasmidBench's cloning simulation tools.** The four that exist are
sufficient; they need tests, not successors.

### PlasmidBench's distinctive territory

Three capabilities have no equivalent in either package:

1. **Primer stock search** — searching the laboratory's accumulated oligo
   inventory for tubes that can already read a chosen region on either strand.
2. **Plasmid and competent-cell stock** — planning against what the freezer
   holds.
3. **Expression strategy** — assessing likely solubility and selecting
   construct and host from available strains, including co-expression.

These lab-aware features are the defensible direction. SnapGene designs in the
abstract and has no model of a freezer.

### Recommended positioning

**A lightweight, portable, offline-capable DNA workbench that is aware of what
the laboratory physically holds, complementing SnapGene and OpenCloning and
integrating with an eLabFTW-centred record system.**

The advantage is portability, offline operation and freezer awareness — not the
number of cloning methods.

---

## 4. Architecture: Four Items Needing Attention

### 4.1 The claude.ai runtime dependency

The code calls:

```js
window.claude.use('db')          // shared plasmid store
window.claude.use('user')        // identity and write permission
window.claude.use('downloads')   // file delivery
```

This is the Claude Artifacts runtime API, available only when the file is
served as an artifact on claude.ai. Opened locally or served from a laboratory
machine, `window.claude` is undefined, persistence falls back to
`localStorage`, and the interface reports "Saved in this browser only."

Three consequences:

- The implied multi-device synchronisation cannot be deployed on laboratory
  infrastructure.
- Records written through that path would sit in a third-party store outside
  laboratory control, which is incompatible with GLP.
- A future maintainer would not expect the dependency.

It is roughly twenty lines. Replacing it with IndexedDB (section 8) and, when
connected, eLabFTW (section 6) should happen in the first development cycle.

### 4.2 Shared mutable state and global re-rendering

| Measurement | Count |
|---|---|
| `state.<field> = …` assignments | 167 |
| `renderAll()` calls | 35 |

Tools mutate one shared state object and then redraw the whole interface.
Nothing records what changed, so no view can subscribe to a narrow slice and
redraw cost is global. This limits both interface testing and performance on
larger data sets.

### 4.3 Interface built by string concatenation

| Measurement | Count |
|---|---|
| HTML string literals | 1,033 |
| `innerHTML` assignments | 86 |
| `addEventListener` calls | 279 |
| Distinct `data-act` values | 16 |

Correct escaping depends on calling `esc()` at every interpolation. For a local
file the risk is low; once PlasmidBench serves laboratory data from a server, a
missed call becomes an injection site. The low `data-act` count indicates event
delegation was begun but not completed — finishing it would remove most
listener re-binding.

### 4.4 Duplicated helpers

These are each defined two or three times across modules: `esc`, `f1`, `f2`,
`niceStep`, `layoutLabels`, `place`, `planText`, `exportPlan`, `run`. A fix
applied to one copy leaves the others unchanged. Consolidate into a shared
utility module.

Minor: the tools menu lists two different entries both labelled "Edit
sequence" — one replaces the whole sequence, the other edits in place. Rename
one.

### 4.5 Types: JSDoc before TypeScript

A TypeScript migration is not yet warranted. **JSDoc type annotations with
`checkJs` enabled** deliver most of the benefit: types live in comments, files
stay plain `.js`, no migration is needed, and editors and CI check the same
way. The annotations convert to TypeScript readily if the project outgrows them.

The existing JavaScript is substantial and working, the project has one
developer, and build complexity carries real cost when the maintainer may
change. Performance is not a language question — the only demanding operation
is searching 100,000 records, which belongs on a server.

---

## 5. Development Practice and the Publication Clock

### 5.1 Publish the repository now

The Journal of Open Source Software revised its policy in January 2026, moving
from effort proxies such as lines of code toward research impact, intellectual
contribution and open-source practice. It also introduced an eligibility rule:

> Privately developed projects are not eligible until there is a public record
> of open development, with **at least six months of public history** before
> submission, including evidence of releases and public issues or pull
> requests.

The norm is spreading — rOpenSci is testing comparable policies. The clock
costs nothing to start and cannot be back-dated. **Publishing the repository
should precede any further code change.**

### 5.2 Repository contents

- The 17 modular source files and `build.js`
- Automated tests and fixtures
- Documentation: user manual, developer guide, architecture note
- Sample plasmid and primer files
- Issue tracking, tagged releases, built HTML attached to each release
- A licence (section 9)

### 5.3 Workflow and versioning

Feature branches merged through pull requests, one module per request so
changes stay reviewable. Semantic versioning:

- `v0.1.0` — repository public, build reproducible, documented
- `v0.2.0` — characterization tests; claude.ai dependency removed
- `v0.3.0` — validated engine; IndexedDB persistence
- `v0.4.0` — eLabFTW file-exchange workflow
- `v1.0.0` — validated stable release

Each release carries notes, known limitations, test results and examples.

### 5.4 AI-assisted development

AI assistance should remain part of the workflow. Reviewers and employers now
assume it, which changes what is worth demonstrating: generating code is no
longer evidence of skill; validating it is. Five practices:

1. **Keep an `AGENTS.md` in the repository** recording conventions, the test
   oracle, and the rule that engine modules stay free of DOM access. AI
   assistants read it each session, so discipline persists without supervision.
2. **Never generate code and its tests in one pass.** They will agree with each
   other rather than with reality. Tests derive from an independent source.
3. **One module per pull request**, small enough to review properly.
4. **Be able to explain every algorithm in the engine.** An unexplainable
   nearest-neighbour calculation cannot be validated by anyone. This, not
   architecture, is the main risk in AI-assisted scientific software.
5. **Disclose AI involvement** in the repository and any publication.

---

## 6. Scientific Validation and Testing

No tests exist. The engine is pure, side-effect-free and already exports under
CommonJS, so it is testable immediately with no refactoring. This is the
cheapest high-value work available and should precede further restructuring.

### 6.1 Biopython is the right oracle

The laboratory holds no SnapGene licence, so SnapGene cannot serve as a
reference. Biopython is a natural fit: PlasmidBench's enzyme table is REBASE as
distributed with Biopython 1.88. A short Python script can generate expected cut
sites, fragment sizes, translations and GenBank round-trips as fixtures — free,
scriptable, reproducible. Published values serve as the secondary reference for
melting temperature: SantaLucia 1998 for nearest-neighbour thermodynamics,
Owczarzy for salt correction, NCBI genetic code tables.

### 6.2 Characterization tests first

Before engine code changes, capture current behaviour as golden files: for a
fixed set of inputs, record exported GenBank, computed cut lists and primer
melting temperatures. These form the regression net. Without it, neither a
developer nor an AI assistant can distinguish a successful refactor from a
silent regression, and assistants fail confidently.

Testing therefore precedes restructuring in the roadmap.

### 6.3 Priority fixtures

The read-interpretation module is the most urgent target. It maps genomic
positions to codons and must handle reverse-strand CDSs, spliced
`join(a..b,c..d)` features, origin-wrapping genes, and all three together. The
module's own documentation names pGro7's CmR gene as the case that combines
them. That is an excellent first fixture: the logic is hard to verify by
inspection and the edge cases are already enumerated.

### 6.4 Test categories

| Category | Validation approach |
|---|---|
| Sequence operations | Independently verified examples |
| Restriction digestion | Biopython on known plasmids; expected fragment sizes |
| Primer melting temperature | Published SantaLucia/Owczarzy values |
| PCR simulation | Known primer-template combinations |
| Cloning operations | Expected products, checked by hand |
| Read interpretation | Reverse-strand, spliced and origin-wrapping CDSs |
| Alignment | Mismatches, insertions, deletions |
| Circular DNA | Origin-spanning features and cut sites |
| Format import | Sequence, topology, features, annotations preserved |
| Round trips | Export then re-import without loss |
| Inventory | Duplicates, missing data, storage failures |

A set of representative plasmids should be maintained as regression fixtures.
Section 10 describes how that set can serve three purposes at once.

**Scientific validation is a development requirement, not a final step.**

---

## 7. Integration with eLabFTW

eLabFTW is the source of truth for official laboratory records. PlasmidBench
and OpenCloning read from it and propose additions to it.

### 7.1 Division of responsibilities

| System | Responsibility |
|---|---|
| eLabFTW | Laboratory database: plasmid and primer records, samples, inventory, experiments, audit trail |
| OpenCloning | Reproducible cloning design and construct provenance |
| PlasmidBench | Portable analysis, visualization, design, offline work, lab-stock awareness |

### 7.2 GLP requirements and the rule they imply

Good Laboratory Practice expects records that are attributable, legible,
contemporaneous, original and accurate: every change traceable to a named person
with a timestamp, no silent overwriting, retained history. eLabFTW provides
this; PlasmidBench does not and should not attempt to.

The resulting design rule is stronger than "eLabFTW is authoritative":

> **Every write passes through eLabFTW's API so its audit trail records it.** A
> tool holding authoritative records in browser storage and synchronising later
> breaks both attribution and contemporaneity.

PlasmidBench should be read-mostly: fetch, analyse, propose. A construct does
not become verified laboratory stock because a design tool produced it.

### 7.3 Integration approach

A dedicated connector service adds a third system to install, secure and
maintain. Two cheaper options come first:

**Option A — File exchange, no new infrastructure.** Download GenBank from
eLabFTW, open it in PlasmidBench, upload the result as an attachment. This works
today with no code and is the appropriate first proof of concept.

**Option B — Browser calls eLabFTW's API directly.** Each user supplies their
own personal API token, stored in their own browser. No new service, and
permissions are per-user, which GLP favours.

**The deciding question is whether eLabFTW's CORS configuration permits
browser-origin requests.** That is a half-day test and should run before any
connector is designed, because the answer determines the architecture. Only if
both options fail is a server-side component justified, and then as the smallest
possible proxy rather than a synchronisation service.

Credentials must never be embedded in the distributed HTML. Offline edits are
submitted as proposed revisions, never automatic overwrites.

### 7.4 Data standards

- GenBank for annotated sequences; FASTA for plain sequence
- OpenCloning JSON for cloning history
- Stable eLabFTW identifiers for registered constructs and primers
- Separate identifiers for physical samples and aliquots
- Sequence checksums and explicit version numbers
- Explicit links between constructs, experiments and sequencing results

### 7.5 Git-backed YAML

Git provides versioning, diffs and commit-level attribution with no server to
maintain. It does not provide indexed search, concurrent editing, or sane merges
across tens of thousands of records, so it is not the operational store. Two
roles suit it well:

1. **Archival export from eLabFTW** — a periodic YAML dump committed to a
   repository is a permanent, diffable, vendor-independent record, valuable for
   GLP and for any future migration.
2. **Curated reference data** — section 10.

---

## 8. Deployment and Persistence

### 8.1 Standalone mode

Distribute `PlasmidBench.html` as a release artifact: no installation, offline,
local editing and analysis, format import and export, independent of eLabFTW.

**Replace `localStorage` with IndexedDB.** Current persistence uses 18 `pb.*`
`localStorage` keys. `localStorage` is a synchronous string store of roughly
5–10 MB — both a capacity and a responsiveness limit. IndexedDB is asynchronous,
far larger and indexable. The claude.ai store should be removed at the same
time; its limits illustrate why that path was never viable for laboratory use:
500 documents per load and a 200,000-base ceiling per record.

### 8.2 Laboratory-connected mode

```text
Laboratory Server
├── eLabFTW         authoritative records, inventory, audit trail
└── OpenCloning     cloning design and provenance

Researcher's browser
└── PlasmidBench.html
    ├── offline standalone mode
    └── optional direct API access to eLabFTW (per-user token)
```

The same HTML can be served from any static web server. Docker is not needed to
run it, and is justified only if the CORS test fails and a minimal proxy becomes
necessary.

### 8.3 Scaling to 10,000–100,000 records

- eLabFTW remains authoritative; benchmark its search and inventory performance
  at realistic record counts before committing.
- Never load a whole collection into the browser. The current design holds all
  records in memory and must change.
- Paginate and use indexed search; retrieve subsets on demand.
- Evaluate **OpenCloningDB** before building redundant database infrastructure.
- Add a dedicated sequence-search service only if measurements show the need.

The standalone HTML should not be the primary database for 100,000 records.

---

## 9. Licensing and Attribution

Required before the repository goes public:

- **Choose a licence.** MIT or BSD-2 keeps the project maximally reusable. Note
  that pLannotate is GPL-3.0, so incorporating its *code* would impose
  copyleft; adopting its approach does not.
- **Check REBASE redistribution terms** for the bundled enzyme table (612
  enzymes, 15 supplier codes, obtained via Biopython). Academic use is
  straightforward; publication and redistribution warrant a check.
- **Attribute data sources**: REBASE, Biopython, NCBI genetic code tables,
  Kazusa codon usage, and any imported feature libraries.
- **Record AI involvement.**

---

## 10. Auto-Annotation: Use Existing Libraries

`FEATURE_LIB` holds 27 entries of the form `{ name, type, seq }`. The matcher
is stronger than the library it serves: exact substring search on both strands,
origin-spanning aware for circular sequences, with a `chanceHits` function
computing how often a motif of a given length occurs by coincidence, and a
context requirement for short motifs.

Crawling NCBI and Addgene to harvest feature-to-sequence correspondences is not
the right approach, for three reasons.

### 10.1 Curated open libraries already exist

- **pLannotate** (McGuffie & Barrick, *NAR* 2021; GPL-3.0) — the state of the
  art open annotator, drawing on a SnapGene-derived feature set, FPbase,
  Swiss-Prot and Rfam, and additionally detecting partial feature fragments left
  from cloning history.
- **GenoLIB** (Adames et al., *NAR* 2015) — ~2,000 widely used plasmids reduced
  to 1,967 non-redundant parts in 112 categories, in SBOL. Paper is CC BY 4.0.
- **ApE** — ships a plain-text, user-editable default feature library and a
  fluorescent-protein library, in a format close to PlasmidBench's own. Likely
  the easiest direct import.

One caveat: pLannotate's primary feature set traces, via GenoLIB, to SnapGene's
curation. Academic use is very likely fine; a published bundled library warrants
a review of the chain of terms. Cleanly licensed sources for parts of the space:
**FPbase** (open API), **Rfam**, **Swiss-Prot**.

### 10.2 Addgene access goes through its Developers Portal

Addgene requires login to view sequences on its website. Bulk and programmatic
access requires requesting a scope, accepting a data-access licence, and about
five business days for approval; data arrives as JSON refreshed daily. Antibody
and affinity-reagent sequences carry terms forbidding redistribution outside the
laboratory.

Crawling would circumvent an access control. The legitimate route is free and
takes a week.

### 10.3 The gap is the algorithm, not the data

pLannotate uses BLAST, DIAMOND and Infernal — homology search. PlasmidBench uses
exact substring matching, so a CMV promoter variant one base from the stored
sequence is not found. BLAST cannot run in a browser offline, so pLannotate's
databases are not directly usable by PlasmidBench's matcher; converting them to
exact-match entries discards the tolerance that makes them effective.

| Approach | Effort | Offline | Recall |
|---|---|---|---|
| Curate a few hundred canonical exact entries | low | yes | moderate |
| Approximate matching in JavaScript (k-mer seed-and-extend, then banded alignment) | moderate | yes | good |
| Call pLannotate's service when online | very low | no | best |

Recommended: the first now, the second as the substantive project, the third as
an optional "deep annotate" action when a network is available.

**No offline, browser-based, tolerant plasmid annotator exists.** pLannotate
needs a server, SnapGene a licence, ApE an install. The gap is real, well
bounded, and suited to PlasmidBench's distinctive strength.

Practical sequence: request Addgene Developers Portal access now; meanwhile
import ApE's feature file and FPbase to reach a few hundred entries cheaply;
then build the approximate matcher. Store the library as a YAML data file with
per-entry provenance, version-controlled and reviewable in pull requests, rather
than as a literal in the source. Beyond roughly a thousand entries, replace the
per-feature linear scan with a single-pass multi-pattern search (Aho–Corasick).

---

## 11. Succession

The student will graduate. If PlasmidBench becomes laboratory infrastructure,
the laboratory inherits it. This supports several recommendations already made:
plain JavaScript with JSDoc types rather than a framework, documentation written
for the next maintainer rather than the author, and tests as the mechanism by
which someone who did not write the code can change it safely.

It also argues for a defined minimum outcome:

> **Minimum durable outcome:** a public repository containing the modular
> source, a working `build.js`, characterization tests for the engine, the
> claude.ai dependency removed, a developer guide, and a tagged release with the
> built HTML attached.

This is achievable within the next cycle and leaves something another person can
take up. Everything beyond it is additive.

---

## 12. Publication Prospects

Implementation effort is no longer the contribution; producing a plasmid
workbench is now inexpensive. What remains scarce is validation, curation,
domain judgement and demonstrated use. Four options, strongest first:

**1. A benchmark of plasmid feature annotation.** No head-to-head accuracy
comparison appears to exist. pLannotate's own paper compares against SnapGene
and PlasMapper on three plasmids, qualitatively, without metrics. Existing
plasmid benchmarks — PlasEval and the 2025 detection and reconstruction
comparisons — address a different problem: detecting plasmids in short-read
data, not annotating engineered constructs.

Curating a ground-truth set of reference plasmids with verified feature maps and
scoring tools by precision, recall and F1 is publishable, and an AI assistant
cannot produce it: the work is domain judgement — construction history traceable
from the literature, sequence verified, features checked by hand.

It serves three purposes simultaneously: the publication, PlasmidBench's
validation fixtures, and the laboratory's verified vector records for eLabFTW.

**2. The offline browser annotator** (section 10.3). Modest but real novelty.
ApE was published in *Frontiers in Bioinformatics* in 2022, long after release,
largely for being free, used and well described.

**3. The lab-aware features.** Searching an oligo inventory for sequencing
primers, and selecting expression hosts from laboratory stock, appears
unpublished. Suitable for a short application note.

**4. A general description of the workbench.** Not publishable on its own.

Breadth now reads differently than it once did: twenty-five tools without tests
reads as a risk rather than as productivity. The signals of credibility are
tests, external validation, provenance, documentation and visible development
history.

Publication should not drive the technical roadmap — the laboratory needs a GLP
system, and that is the real deliverable — but the overlap is large enough to
pursue deliberately.

---

## 13. Roadmap

Phases 1–4 and the annotation work are a realistic year. Phases 5–7 depend on
the CORS test result and available time.

| Phase | Focus | Deliverable |
|---|---|---|
| **1** | Publish repository, documentation, reproducible build, licence | Public repository with tagged release; six-month clock started |
| **2** | Characterization tests against Biopython; remove claude.ai dependency | Regression net; no third-party runtime dependency |
| **3** | Engine consolidation: de-duplicate helpers, JSDoc types, module interfaces | Validated, type-checked engine |
| **4** | IndexedDB persistence; complete event delegation; escaping audit | Stable standalone release |
| **5** | eLabFTW file exchange, CORS test, direct API if permitted | Working laboratory workflow, no new services |
| **6** | Feature library import and approximate matcher | Offline tolerant annotation |
| **7** | Ground-truth plasmid set and benchmark | Validation fixtures and a publishable result |

### Immediate deliverables

**1 — Public repository.** Modular source, `build.js`, README, build
instructions, licence, issue tracking, `AGENTS.md`, and an initial tagged
release with the built HTML attached. First, because the publication clock
cannot be back-dated.

**2 — Characterization tests.** Golden files for GenBank export, restriction cut
lists and primer melting temperatures, with expected values from Biopython.
Priority: the read-interpretation module (section 6.3) and the format parsers,
where malformed real-world files cause most failures.

**3 — Remove the claude.ai dependency** and move browser persistence to
IndexedDB.

**4 — Documentation.** User manual, developer guide, an architecture note
describing the 17 source files and the build, and an explicit list of supported
formats and known limitations.

**5 — eLabFTW proof of concept.** The file-exchange workflow: a plasmid
exported from eLabFTW, analysed in PlasmidBench, returned as a proposed new
construct without overwriting the original. Run the CORS test and report the
result.

---

## 14. Long-Term Direction

PlasmidBench need not replace SnapGene or duplicate OpenCloning. Its value rests
on three strengths:

1. **Standalone portability** — a complete workbench from one HTML file, no
   installation, no server.
2. **Scientific reliability** — independently validated molecular biology
   functions, continuously tested.
3. **Laboratory awareness and interoperability** — working with eLabFTW,
   OpenCloning and standard formats, and knowing what the laboratory holds,
   without creating a competing record system.

The intended ecosystem:

- **eLabFTW** — authoritative records, inventory, experiments, audit trail
- **OpenCloning** — reproducible cloning workflows and provenance
- **SnapGene** — commercial reference for comprehensive design and analysis
- **PlasmidBench** — portable, offline, laboratory-aware analysis and design

---

## Summary

The architecture is right: modular source, a clean engine free of interface
code, already importable for testing, and a single-file distribution that is a
genuine advantage. The recent read-interpretation work — mapping genomic
positions to codons correctly across reverse-strand, spliced and
origin-wrapping features — is careful and addresses a question that matters at
the bench.

What is missing is evidence rather than structure: tests, a public development
history, independent validation, and removal of one runtime dependency.

Further features are now unlikely to add as much value as demonstrating that
the existing ones are correct. The next phase should emphasise publishing the
repository, independent scientific validation, and the thinnest workable
eLabFTW integration.

**The next significant achievement should be a public, version-controlled,
independently tested repository, and a standalone HTML release that someone
other than its author can maintain.**
