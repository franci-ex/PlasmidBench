# What PlasmidBench does not do, and what is approximate

Written down because a tool that overstates itself is dangerous at a bench,
and because the next person deserves to know where the edges are without
having to find them.

## Not a record system

It holds no authoritative records and should not. Attribution, timestamps,
no silent overwriting, retained history — eLabFTW provides those and this
does not. It reads from eLabFTW and **adds** to it, never replacing or
patching. A construct does not become verified lab stock because a design
tool produced it.

If you are using this as your record of what is in the freezer, stop. Export
to eLabFTW.

## Not a database

Everything is held in the browser. IndexedDB will take a great deal more
than localStorage would, but the design loads whole collections into memory,
and past a few thousand records that is the wrong architecture. For 10,000 to
100,000 records, eLabFTW or OpenCloningDB is the answer and this should be
reading from them.

## The feature library is small

Twenty-seven built-in entries. SnapGene ships thousands, curated. A plasmid
annotated here with only the built-in library will come back nearly empty,
and the panel says so rather than letting you conclude there is nothing
there.

Import one — ApE's feature file reads directly, FPbase has an open API — and
the imported entries carry their source and licence with them.

## The tolerant annotator is not BLAST

It reports identity, coverage, and the chance of its seed being a
coincidence. It does **not** report an E-value, because an E-value needs a
database model it does not have, and a number that looks like a statistic
but is not one is worse than no number.

Its measured behaviour, on a 1,200-entry library and 5 kb sequences:

| divergence | found |
|---|---|
| exact | 100% |
| 1 base in 50 | 99% |
| 1 in 20 | 94% |
| 1 in 10 | 88–92% |

0.03 false positives per 5 kb. Past roughly 15% divergence a 12-base seed
has little chance of landing cleanly, and the thing will be missed. If you
need homology search proper, pLannotate's service does it better and this
page offers no substitute for it when you are online.

Two specific behaviours worth knowing:

- **An insertion larger than a few bases splits a feature into fragments**
  rather than being reported as one feature with a hole. That is deliberate:
  a terminator with thirty bases of unrelated sequence in the middle is two
  pieces of one, which is the more useful finding.
- **The band width is the largest indel an alignment can contain.** A hit
  that reached the edge of it says so.

## Melting temperature

Nearest-neighbour (SantaLucia 1998) with Owczarzy salt corrections, checked
against published values. The usual caveats apply and are not this project's
to fix:

- It is a **prediction**, and predictions from different tools disagree by a
  degree or two. A primer set that works is the only real evidence.
- For a **primer with a 5' tail**, two temperatures matter: the annealing
  portion governs the early cycles and the whole primer the later ones, once
  the tail is in the product. Both are shown, and the annealing one is the
  one to set the block by. This was wrong until November 2026 — the library
  showed only the whole primer, which can be fifteen degrees out.
- Nothing models **secondary structure in the template**, which is a common
  reason a predicted-fine primer fails.

## Expression strategy

The solubility assessment is a set of published heuristics — hydropathy,
low-complexity content, disorder propensity, cysteine count — not a
prediction from a trained model. It is a way of ordering constructs by how
much trouble they are likely to be, and it will be wrong about individual
proteins. Treat it as a reason to try two constructs rather than as a reason
not to try one.

## Sanger interpretation

- The read-against-construct matching uses **k-mer containment at k=18**
  with a contiguity requirement. It is reliable on a good read and will give
  up on a bad one, which is the right failure.
- **The confident read length and dead zone are a setting**, because
  facilities differ a great deal. The defaults are textbook; set yours in
  the lab profile or the planner will draw a plan your facility cannot
  deliver.
- Base calling is taken from the trace file as the instrument called it.
  Nothing re-calls peaks.

## Gels

The mobility model is a fit to log(size) over an agarose gel's resolving
window. It is for *recognising* a band you are looking at, not for measuring
one: do not read a fragment size off the simulated gel and treat it as data.

## Formats

- **Writes GenBank, EMBL, FASTA, GFF3 and JSON. Does not write SnapGene
  `.dna`** — it reads it. SnapGene reads GenBank, so GenBank is the route
  back.
- A **protein** GenBank or FASTA record is read and recognised as protein,
  and then refused for DNA operations, with Back-translate offered. It is not
  silently treated as DNA. (It was, until November 2026: two residues in five
  were being deleted by DNA cleaning.)
- Features **referring to other accessions** — `join(AL121804.2:41..610,...)`
  — are declined and listed in `skippedFeatures`, because they cannot be
  placed on this sequence. Declining them and saying so is the intent; they
  are not lost silently.
- **Fuzzy and malformed locations** (`one-of(1888,1901)..2479`,
  `complement(-2..492)`) are declined the same way.

## The interface

- **One shared mutable state object, redrawn globally.** Fine at this scale;
  it limits performance on large data and makes narrow interface tests
  harder.
- **Markup is built by joining strings**, so correctness rests on `esc()`
  being called at every interpolation. It is, and
  `tests/uitest-escape.js` is what keeps it so, but the structure of the code
  does not enforce it.
- **No mobile layout.** It expects a desktop-sized window.

## Browser storage

- IndexedDB is per-browser and per-origin. Opening the file from a different
  folder is a different origin, and your library will not be there. **Export
  anything you intend to keep.**
- A private window, or blocked site data, means nothing is kept. The page
  says so as a warning rather than appearing to save.
- Nothing synchronises between machines on its own. That is what eLabFTW is
  for.

## Known rough edges in the code

Listed so they are not discovered as surprises:

- Event delegation is begun and not finished: 16 `data-act` values against
  279 `addEventListener` calls.
- `src/app-tools.js` is large and holds most of the dialogs.
- The `preview.html` and screenshot scripts in `tests/` are development aids,
  not part of the build.
