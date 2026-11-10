# Changing PlasmidBench

Read `AGENTS.md` first — it has the conventions, and the one rule about
tests. This is the practical companion: how to build it, how to run the
checks, and how to do the things you are most likely to need to do.

## Setting up

```sh
git clone https://github.com/<owner>/PlasmidBench.git
cd PlasmidBench
node build.js            # writes plasmidbench.html
```

That is all that is needed to build. For the checks:

```sh
npm install --no-save playwright@1.48.2 typescript@5.6.3 @types/node@22
npx playwright install chromium
pip install biopython==1.88
```

and the corpus, which `tests/README.md` describes in one command.

## The loop

```sh
# change a file in src/
node build.js
node tests/test.js              # the suite covering what you touched
./validate-all.sh engine        # all seventeen, about a minute
./validate-all.sh               # everything, about 25 minutes, before a PR
```

Build before running a browser suite: they test the built file, not the
source. Forgetting this produces a confusing green run against yesterday's
code.

## Where things are

`docs/architecture.md` has the full map. The short version:

- **Science** goes in `src/engine*.js`. No DOM. If you need the interface,
  you are in the wrong file.
- **Interface** goes in `src/app-*.js`, and reaches the engine only through
  `Bio.*`.
- **A new tool** is a function in `src/app-tools.js` (or its own
  `src/app-*.js` if it is large) plus a row in the `TOOLS` table.
- **A new module** must be added to `build.js`, to the `require` chain of the
  next engine module, and to `tsconfig.json`.

## How to do the common things

### Add a sequence format

1. A reader in `src/engine5.js` (or `engine.js` for a text format), returning
   `SeqRecord[]`.
2. Register the extension where `PB.handleFiles` dispatches, in
   `src/app-tools.js`.
3. **Write the oracle before you write the reader.** If Biopython reads the
   format, `tests/oracles/cmpgb.py` is the pattern: a JavaScript half that
   reports and a Python half that compares. If nothing else reads it, say so
   in the test and compare against files you can check by hand.
4. Add a round trip if you also wrote an exporter, and make the reading side
   of the comparison Biopython — never your own reader, which would hide a
   writer that drops a feature behind a reader that cannot see it.

A format where no independent implementation exists is a format you cannot
validate. Say that in the test rather than leaving it implied.

### Add an enzyme

Don't — the table is REBASE and is compared against Biopython's copy of it.
Users add their own through the enzyme panel, which stores them in
`customEnzymes`. If REBASE itself has changed, regenerate the table and let
`tests/oracles/test-rebase.py` tell you what moved.

### Change a calculation

1. Find the published source and put it in the comment, with the paper.
2. Add or extend the oracle. If there is no outside source for the quantity,
   that is worth knowing before you ship it.
3. If a result changes, **say by how much** in the pull request. Someone's
   plasmid was designed with the old number.

### Change the interface

The browser suites test what a person sees. Add to the one that covers the
area; if a behaviour cannot be expressed as something a person would notice,
that is a sign it should not be a behaviour.

Every value from a file goes through `esc()`. `tests/uitest-escape.js` will
catch a missed one, but only for a view it opens — if you add a view, add it
to the list there.

### Add a setting to the lab profile

1. A row in `FIELDS` in `src/profile.js`, with a line saying what it is **for**
   rather than what it is called. That line is what the panel shows.
2. Make it take effect: either read it at the point of use, or add it to
   `PB.reloadSettings()` in `src/app-core.js` if it is cached at start-up. A
   setting that needs a reload is a setting people stop trusting.
3. If it is a credential or a record list, it must not be exportable — check
   that `NEVER_EXPORT` in `src/profile.js` already catches it, and add a
   pattern if not.

## Debugging

The engine is plain Node, so the fastest way in is a one-liner:

```sh
node -e "
  const B = require('./src/engine8.js');
  const fs = require('fs');
  const rec = B.parseGenBank(fs.readFileSync('tests/fixtures/pGro7.gb','utf8'))[0];
  console.log(rec.name, rec.seq.length, rec.topology, rec.features.length);
  const cat = rec.features.find(f => /CmR|cat/i.test(f.label||''));
  console.log(JSON.stringify(B.featurePieces(cat, rec.seq.length)));
"
```

For the interface, open `plasmidbench.html` and use `window.PB` — the whole
state and every render function are on it, which is deliberate and is what
the browser suites drive.

```js
PB.state.records           // what is loaded
PB.active()                // the current construct
PB.renderAll()             // redraw
PB.store.allRecords()      // what IndexedDB holds
```

A Playwright suite that fails is easiest to read with the browser open:

```js
// at the top of the suite
const browser = await chromium.launch({ headless: false, slowMo: 300 });
```

## Things that will catch you

- **Two ways a feature crosses the origin**: `end <= start`, and
  `end > length`. Both are legal and both occur. Go through
  `E.featurePieces`.
- **0-based, end-exclusive, everywhere.** The 1-based numbers belong to
  GenBank's text, not to the engine.
- **`E.pcr` returns an array** that also carries `found`, `sites` and
  `truncated`. There is no `.products`.
- **Settings read synchronously** from an in-memory cache that `store.js`
  fills. `PB.lsSet` is fire-and-forget; it does not await IndexedDB.
- **The built file is what ships.** A change that works in `src/` and not in
  the build is usually a `require` line the build strips, or a module not
  added to `build.js`.
- **`crypto.subtle` does not exist on a `file://` page.** That is why SHA-256
  is in the engine. Anything else needing crypto has the same problem.
- **Google Fonts is the only external request**, and it fails offline on
  purpose. Test suites filter that one console error; do not add more.

## Releasing

```sh
git tag -a v0.4.0 -m "what changed"
git push origin v0.4.0
```

The release workflow builds and attaches `PlasmidBench-v0.4.0.html`. The
point of attaching the built file is that a result can be traced to the exact
file that produced it — so do not retag, and do not replace a release asset.

Version numbers follow the roadmap in `docs/evaluation.md` §5.3.
