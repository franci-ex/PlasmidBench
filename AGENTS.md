# Working on PlasmidBench

Read this before changing anything. It is for whoever works on this next,
and it is read by AI assistants at the start of every session, which is most
of why it exists: discipline that lives in a file persists, and discipline
that lives in someone's head does not survive them graduating.

The project's own evaluation (`docs/evaluation.md`) argues that in
AI-assisted scientific software the risk is not architecture — it is code
nobody can explain, validated by tests that agree with it rather than with
reality. Everything below is aimed at that.

---

## The one rule

**A test and the code it tests must not come from the same place.**

Generating code and its tests in one pass produces two things that agree
with each other and need not agree with anything else. So:

- Expected values come from an **outside source**: Biopython, REBASE, a
  published formula, a vendor's spec sheet, or arithmetic worked out in the
  test from first principles and written down there.
- The oracles in `tests/oracles/` are **two halves that share no code**. A
  JavaScript half reports what the engine says and does not compare. A
  Python half computes the expected answer from an independent library and
  compares. Neither half can be tuned to agree with the other without
  someone noticing.
- **Never change a test expectation to make it pass.** Find out which side
  is wrong first. If it is the test, say in a comment what settled it — the
  source, the arithmetic, the date. There are a dozen such comments in
  `tests/` and each one is a small piece of evidence.

If you find yourself about to write "the old expectation was wrong", write
*why*, and name what you checked it against.

---

## The shape of it

```
src/
  engine.js … engine8.js   the science. No DOM, no window, no interface.
  store.js                 IndexedDB persistence
  elab.js                  eLabFTW API client
  profile.js               the lab profile
  xlsx.js                  a spreadsheet writer, from scratch
  app-*.js                 the interface
  part1-head.html          the shell and the stylesheet
  part2-body.html
build.js                   joins them into one HTML file
tests/                     unit suites and browser suites
tests/oracles/             the outside comparisons
docs/                      manual, developer guide, architecture note
```

**The engine must stay free of interface code.** No `document`, no `window`,
no `addEventListener`, in any `src/engine*.js`. It is checked:

```sh
grep -nE '\b(document|window|addEventListener|innerHTML)\b' src/engine*.js
```

That must print nothing. It is why the engine can be required into Node and
tested at all, and it is the single most valuable property this codebase has.

`require('./src/engine8.js')` loads the whole engine: each module requires
its predecessor. (That chain was broken until November 2026 — modules 3 to 6
were not in it — so "the engine is importable" was true only if you required
all eight files in the right order. If you add a module, add it to the chain.)

---

## Conventions that are not negotiable

**Positions are 0-based and ranges are end-exclusive. Everywhere.** GenBank's
1-based inclusive numbers exist only in the text a reader consumes and a
writer produces. Every bug in the first year of this project that took more
than an hour to find was an off-by-one at that boundary.

**A feature can cross the origin, and there are two ways files write it.**
`end <= start`, and `end > length`. Both occur in real files; both are legal;
`E.featurePieces` resolves either. Never compare a feature's `start` to its
`end` without going through it.

**A value that came from a file is text until proven otherwise.** Every
interpolation into markup goes through `esc()`. `tests/uitest-escape.js`
loads a construct whose every field is markup and opens every view and every
tool, which is the only way to cover four hundred interpolation sites.

**Say what is not known.** A number with no basis, a feature called with no
evidence, a file "saved" that was not — those are worse than an absence. The
interface has several places where the honest answer is "this could not be
determined", and they are deliberate.

---

## Before you commit

```sh
./validate-all.sh                 # everything: ~25 minutes
./validate-all.sh types           # tsc --checkJs, must be 0 errors
./validate-all.sh engine          # the unit suites, ~1 minute
./validate-all.sh oracles         # the outside comparisons
./validate-all.sh browser         # Playwright against the built file
```

All four layers green, or say in the pull request which is not and why.

`npx tsc -p tsconfig.json` must report **zero** errors. The types live in
JSDoc and in `types.d.ts`; there is no TypeScript to compile and `src/` stays
plain JavaScript.

---

## Pull requests

**One module per pull request**, small enough that someone can actually read
it. The evaluation is right that this is the part AI assistance makes harder,
not easier: it is easy to produce a thousand lines that look fine.

A pull request should say:

1. what changed and why,
2. what it was validated against — name the oracle or the source,
3. anything that was measured, with the numbers.

"Measured, not assumed" is the house style. Several decisions in this
codebase are written down with their measurements because the obvious guess
was wrong:

- the whole primer's melting temperature is *not* always higher than its
  annealing portion's — measured at 0.4% of random primers;
- the seed length for tolerant annotation, chosen from recall and false-
  positive rates across four values rather than picked;
- the false-positive rate of that matcher, which looked like 15.6 per 5 kb
  until the random control turned out to share a cycle with the library.

That last one is worth knowing about generally: **the LCG in the older test
suites has a single cycle**, so two different seeds are two starting points
on one sequence. If a test needs sequence that genuinely does not contain
something, use mulberry32 as `tests/test19.js` does.

---

## Things deliberately not done

Say so rather than quietly doing them:

- **No crawling Addgene.** Their Developers Portal is free and takes about a
  week; crawling would circumvent an access control.
- **No hand-built feature library.** Curated open ones exist (ApE, FPbase,
  GenoLIB) and are imported, with provenance per entry.
- **No TypeScript migration.** JSDoc plus `checkJs` gives most of the benefit
  with no build step, and converts later if the project outgrows it.
- **No new cloning simulators.** The four that exist are enough; OpenCloning
  covers that ground. They need tests, not successors.
- **No authoritative records in this tool.** eLabFTW holds the records and
  the audit trail. This reads from it and *adds* to it, never replaces or
  patches. A construct does not become verified lab stock because a design
  tool produced it.
- **No credentials in the built file.** There is nowhere in the source for an
  API key to live. The lab profile cannot export one; that is enforced by
  pattern in `src/profile.js`, not by remembering.

---

## If you are an AI assistant

Everything above applies to you, and these in particular:

1. **Read the comment before changing the code.** Most non-obvious lines in
   this codebase have a note saying what went wrong before and what the
   measurement was. Several look like mistakes and are not.
2. **Do not generate a test and the code it tests in one pass.** Write the
   code, then validate it against something that existed before you did.
3. **If a test fails, find out which side is wrong.** Do not adjust the
   expectation to match the output. Say what settled it.
4. **Measure before claiming.** "This should be faster" and "this will find
   more" are hypotheses. There are measurement harnesses in
   `tests/oracles/`; add to them.
5. **Be able to explain every line you add.** An unexplainable
   nearest-neighbour calculation cannot be validated by anyone, which makes
   it worthless in a lab regardless of whether it is right.
6. **Record that you helped**, in the pull request. `NOTICE.md` says AI
   assistance is part of this project's history; it is not a thing to hide,
   and it is also not a thing to leave implicit.
