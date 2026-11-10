# Using PlasmidBench

Open `plasmidbench.html`. That is the whole installation.

Everything is kept in that browser, on that machine. **Export anything you
intend to keep** — and if the lab runs eLabFTW, put it there, because that is
where records belong.

---

## Getting a sequence in

**Import**, or drop a file on the window. It reads GenBank, EMBL, FASTA,
SnapGene `.dna`, ABI `.ab1` traces, GFF3, Benchling and openVectorEditor
JSON, Open Vector Format, and NCBI INSDSeq/TinySeq XML. Paste works too:
**Tools → Replace the whole sequence**.

A few things it will tell you that are worth reading:

- **"The LOCUS line says linear but a feature crosses the origin, so this is
  being read as circular."** Addgene and SnapGene both write `linear` on
  plasmid exports. The correction is right, and it is announced rather than
  made quietly.
- **"N features were declined"** with a list. A location referring to another
  accession, or a negative coordinate, cannot be placed on this sequence.
  They are named so that you know what is missing.
- **A protein file** is recognised as protein and refused for DNA work, with
  **Back-translate** offered.

---

## Designing an assembly

**Tools → Seamless assembly** (In-Fusion, Gibson, NEBuilder — the same
design, with the arm length set by the kit).

Pick where the insert goes: a selection on the map, a feature to replace, or
positions typed in. The panel then shows the primers, the junctions, and the
plasmid you will get.

### Three things in here that are easy to get wrong

**The melting temperature of a tailed primer.** A primer with a 5' homology
arm anneals only along its 3' part in the first cycles, and *that* part sets
the temperature the block must run at. The whole primer's is usually fifteen
degrees higher, and a block set to it anneals above what holds, which gives
nothing. The library shows **both**, annealing first. Set the block by the
annealing one.

**Reusing a backbone pair across a series.** If you are putting a dozen
different inserts into the same place in the same vector, you want the same
two backbone primers every time — one tube, one order, one linearisation
that you know works. The **"Use primers you already have"** dropdown in the
assembly panel offers any saved pair that opens the loaded vector, and
choosing one keeps its own name rather than relabelling it after today's
insert. Saving the design again will not add a second copy of it.

When a pair is chosen, the placement above it is no longer deciding anything
— and if the two disagree, the panel says so and says which won. Read that
warning; it is the one that would otherwise be found on a gel.

**Replacing an annotated feature.** Selecting a CDS and asking for it to be
replaced removes it. If the feature crosses the origin, that still holds —
there was a bug here until November 2026 where a wrapping selection removed
nothing and the primers kept the old sequence.

---

## Primers

**Primers** in the view bar. Each row shows the annealing temperature, the
whole primer's underneath, and what the split is based on:

| what it says | what it means |
|---|---|
| `58.6` / `full 72.1` / `15 nt arm, found` | the arm was not recorded, so it was found by searching the vector for it |
| `…, set by you` | you told it where the arm ends |
| one number only | no overhang: the whole primer anneals |
| `No overhang is known` on hover | click the temperature to say where the arm ends |

A primer more than 3 °C from your target is flagged, and a pair more than
4 °C apart is flagged, because that is the pair that gives one product and
one smear.

An arm found by searching errs **low by about a degree**, never high. A
temperature that is slightly too low costs nothing; one that is too high
costs the reaction.

---

## Checking sequencing

**Verify** in the view bar. Drop in `.ab1` traces or `.seq` files.

It tells you, for each read:

- **what matched the construct**, and where;
- **what is new** — present in the read, not in the construct, with the
  clean length of it;
- **what differs**, base by base;
- **what is unreadable**, and whether it sits inside an otherwise good run;
- **the ORFs**, in all three frames, with the amino acids **and** the
  nucleotides, each copyable.

"What of this read is actually in the plasmid" is its own section, because
that is the question you came with.

### Blast

**Run BLAST** opens NCBI with the read in it. Paste the result back
(**Paste BLAST results**) and the hits are read in and shown against the
construct. The page does not call NCBI itself — it holds no network
credentials and makes no outbound requests of its own.

---

## Annotating

**Tools → Annotate.** Two searches, and the panel says which found what:

- **exact** — in the library letter for letter;
- **approximate** — something very like it: a variant, or a fragment a
  cloning step left behind.

A row reading `95.8% identical`, `62% of it` is telling you something true
and specific. A row reading `exact`, `all of it` is telling you something
different. The distinction is the point.

**The 27 built-in entries are not enough** for a real plasmid, and the panel
says so. Import a library — ApE's feature file reads directly, and FPbase
has an open API for fluorescent proteins. Every imported entry keeps its
source and licence, shown in a *From* column.

Nothing is written into your construct until you tick it. What the match was
stays on the feature afterwards, so a call made at 91% identity still says so
next year.

---

## What the lab holds

Three imports, once each, and then every plan knows about them:

- **Tools → Primer stock** — your oligo inventory. **Tools → Sanger primers
  from stock** then searches all of it for tubes that can already read a
  region, on both strands, and tells you which box.
- **Tools → Plasmid stock** — the collection.
- **Tools → Competent cells** — the strains, with what each has already spent.

The importer reads a spreadsheet as it is, and if the header does not line up
with the rows it says so rather than guessing: a PDF-extracted header with
run-together labels used to shift the location column and send you to the
wrong box.

---

## Making it yours

**Tools → Lab profile.** The lab name, the bench conditions **every melting
temperature on the page is calculated at**, which suppliers you can buy
from, what your sequencing facility really gives you, your spreadsheet
layouts, your eLabFTW address.

Two of those change results immediately:

- **Suppliers** filters the enzyme list to what you can order. Without it you
  are offered all 607, which is how a design comes to rest on a tube nobody
  can buy.
- **Read length** changes what the sequencing planner draws. The default is
  700 confident bases; if your facility gives 500, a plan drawn at 700 leaves
  a gap you find after paying for the reads.

**Export this profile** gives one file another lab can import. It carries no
API keys and no records — only configuration.

Importing shows you every setting it would change, with the old and new
value, before changing anything.

---

## eLabFTW

**Tools → eLabFTW.** Put in your lab's address and **your own** API key
(eLabFTW → Settings → API keys), choose whether to remember it, and **test
the connection**.

The test tells you which of five things is true, and that matters because
"the server blocks browsers" and "your key was rejected" look identical from
inside a browser if you only try once. If it says **blocked by CORS**, the
panel names the two headers an administrator needs to add, and describes the
file route that works meanwhile with nothing installed.

Once connected: search for a record, open the GenBank attached to it, work on
it, and **attach the result back**. The result is *added* — a new attachment
plus a comment saying what made it, when, and the sequence checksum. Nothing
already there is changed or replaced, and a download whose bytes do not match
the checksum eLabFTW holds is refused rather than opened.

`docs/elabftw.md` has the full account.

---

## Where your work lives

Bottom of the window:

| what it says | what it means |
|---|---|
| **Saved in this browser.** | IndexedDB, on this machine, in this browser. Export to move it. |
| **… Also mirrored to your account** | you are in the hosted artifact, and it reaches your other devices too |
| **This browser will not keep data** | a private window, or blocked site data. **Export before closing.** |

IndexedDB is per-origin: opening the file from a different folder is a
different origin and your library will not be there. Nothing synchronises
between machines on its own — eLabFTW is for that.

---

## Things worth knowing

- **Undo** works on sequence edits, annotation and the structural tools.
  Ctrl/Cmd-Z.
- **Saving a file** just downloads it when the page is opened normally. Inside
  an embedded viewer a click can be dropped silently, so there you get a link
  to click yourself — which is honest rather than claiming a save that did
  not happen.
- **The simulated gel is for recognising a band, not measuring one.** Do not
  read a fragment size off it and treat it as data.
- **The expression assessment is a set of published heuristics**, not a
  prediction from a trained model. It is a reason to try two constructs, not
  a reason not to try one.
- **It works offline.** The only thing it asks the network for is a font.

`docs/limitations.md` is the full list of what it does not do, and what is
approximate. It is worth ten minutes before you rely on any of this.
