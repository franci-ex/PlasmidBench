# How PlasmidBench is put together

For whoever maintains this next. It says what each file is, why the shape is
what it is, and which parts are load-bearing.

## The shape in one paragraph

Thirty source files are joined by `build.js` into one HTML file with no
dependencies. Eight of them are the engine: the science, with no interface
code in it at all, which is why it can be required into Node and tested.
Twelve are the interface, which reaches the engine only through its public
namespace. The rest are the shell, the stylesheet, a spreadsheet writer, and
the three modules that talk to the outside: persistence, eLabFTW, and the
lab profile.

```
                         build.js
                            │
   ┌────────────────────────┼────────────────────────┐
   │                        │                        │
part1-head.html      src/engine*.js ×8        src/app-*.js ×12
(shell + CSS)        the science               the interface
part2-body.html      no DOM, no window         reaches the engine
                     240 exports               only through Bio.*
                            │
                     src/store.js      IndexedDB
                     src/elab.js       eLabFTW
                     src/profile.js    the lab profile
                     src/xlsx.js       .xlsx writer
                            │
                     plasmidbench.html  (1.4 MB, one file)
```

## The engine

Eight modules, 10,500 lines, **zero** references to `document`,
`addEventListener` or `innerHTML`. Each requires its predecessor, so

```js
const Bio = require('./src/engine8.js');   // all 240 exports
```

loads the whole thing. (That chain was broken until November 2026 — modules
3 to 6 were not in it — so "the engine is importable" held only if you
required all eight in order. If you add a module, add it to the chain and to
`tsconfig.json`.)

| file | what is in it |
|---|---|
| `engine.js` | sequence operations, the genetic codes, GenBank/EMBL/FASTA reading and writing, locations and features, 607 enzymes, the exact feature matcher |
| `engine2.js` | primers and melting temperatures, PCR, assembly design, inverse PCR, Sanger read interpretation, ABI traces, reaction setup |
| `engine3.js` | the lab's stock: oligo inventory, plasmid and strain lists, the sequencing planner |
| `engine4.js` | protein descriptors, expression strategy, host selection, codon usage |
| `engine5.js` | the JSON and XML dialects, SnapGene `.dna`, GFF3 |
| `engine6.js` | alignment and the figures drawn from it |
| `engine7.js` | BLAST result parsing, checksums |
| `engine8.js` | tolerant annotation: the k-mer index, banded Smith-Waterman, Aho-Corasick, library import |

### Conventions it holds to

**Positions are 0-based, ranges end-exclusive, everywhere.** GenBank's
1-based inclusive numbers exist only in the text a reader consumes and a
writer produces. Every long-lived bug in this project has been an off-by-one
at that boundary.

**A feature can cross the origin, written two ways.** `end <= start`, and
`end > length`. Both occur in real files. `E.featurePieces` resolves either,
and nothing should compare a `start` to an `end` without going through it.
This was the cause of a whole gene silently vanishing from assembly
products: the length printed correctly and the location printed correctly,
and only the extracted bases were short.

**Some functions return an array that also carries properties.** `E.pcr`
returns the products with `found`, `sites` and `truncated` on the array.
That is written out as an intersection type in `types.d.ts` rather than left
to be discovered.

**Types live in comments.** `types.d.ts` holds the shapes that cross module
boundaries; the rest is JSDoc. `npx tsc -p tsconfig.json` checks them and
must report zero errors. There is no TypeScript to compile and `src/` stays
plain JavaScript that Node can require as it is.

## The interface

Twelve modules. They share one `state` object and redraw through
`PB.renderAll()`. Markup is built by joining strings, which means **every
interpolation of a value from outside goes through `esc()`** — four hundred
of them, which is why the coverage is behavioural: `tests/uitest-escape.js`
loads a construct whose name, description, every feature label, every
qualifier and every primer name is markup, then opens every view and every
tool and asks whether anything ran.

| file | what it draws |
|---|---|
| `app-core.js` | state, persistence, undo, the library, the enzyme rail, file saving |
| `app-views1.js` | the sequence and feature views |
| `app-map.js` | the plasmid map, including the label layout that fans a crowded side into columns |
| `app-views2.js` | primers, gels, sequencing verification, alignment |
| `app-tools.js` | the tools menu and most of the dialogs |
| `app-seq.js` | the sequencing planner |
| `app-express.js` | the expression panel |
| `app-annotate.js` | the annotate panel |
| `app-elab.js` | the eLabFTW panel |
| `app-profile.js` | the lab profile panel |
| `app-init.js` | start-up, theme, keyboard |
| `xlsx.js` | a spreadsheet writer, written from scratch, because a dependency would defeat the single file |

### Known rough edges

Honest about these, since the next person will meet them:

- **One shared mutable state object, and global redraws.** 167 assignments to
  `state.*` and 35 calls to `renderAll()`. Nothing records what changed, so
  no view can subscribe to a slice and redraw cost is global. It is fine at
  this scale and it limits both interface testing and performance on large
  data.
- **Event delegation is begun, not finished.** 16 `data-act` values against
  279 `addEventListener` calls. Finishing it would remove most listener
  re-binding.
- **Markup by string concatenation.** Correct because `esc()` is called
  everywhere, and that correctness rests on a test rather than on the
  structure of the code.

## Persistence

`store.js`. **IndexedDB is the system of record**, in the browser the file
is opened in. It works with no network, no account and no host, which is
what a page served from a lab machine needs.

localStorage held the library until November 2026 and could not: it is a
synchronous string store of five to ten megabytes for the whole origin,
settings and sequences together, and a write that does not fit throws where
it stands. A single 100 kb plasmid is 100 kb of JSON.

Two things worth knowing about how it is arranged:

- **Settings still read synchronously.** Around a hundred call sites read a
  setting in the middle of drawing something. They are held in memory, seeded
  from localStorage at start-up (which is already synchronous and already
  holds the old values), and written through to IndexedDB in the background.
- **A host-side store is mirrored to when one exists** — the claude.ai
  artifact runtime is the one that does — so that constructs saved that way
  are not stranded and someone who has been opening the page on two devices
  does not silently lose the second. `pb.mirrorHost = false` turns it off for
  a lab install, and nothing depends on it.

When IndexedDB is refused — a private window, blocked site data — the page
says so as a warning and tells you to export, rather than appearing to save
and losing the lot on close.

## The build

`build.js` concatenates the parts in order, stripping the `require` lines
that only exist for Node. It takes about a second. The output is
`plasmidbench.html`, 1.4 MB, which is what a release attaches.

The only external resource the built file asks for is a Google Fonts
stylesheet, and it falls back cleanly when that is unreachable — which is
what offline means. (Test suites filter that one console error.)

## Validation

`validate-all.sh`, four layers, about 25 minutes. See
`docs/developer-guide.md` for how to run parts of it, and `tests/README.md`
for the corpus.

The architecturally important part: **the oracles are two halves that share
no code.** A JavaScript half reports what the engine says and does not
compare; a Python half computes the expected answer from Biopython and
compares. That is the only arrangement in which the tests are evidence
rather than a restatement of the code — which matters especially for a
codebase written with AI assistance.
