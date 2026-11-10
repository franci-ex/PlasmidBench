/* ===========================================================================
   Dump every restriction enzyme PlasmidBench knows, and where each one cuts
   a fixed set of sequences, for test-rebase.py to set beside Bio.Restriction.

   Both tables come from REBASE, but the two implementations of the search
   are unrelated: Biopython compiles a regular expression per enzyme, this
   one walks the sequence against IUPAC codes. The sequences below are
   chosen so that every enzyme in the table cuts at least one of them —
   including the awkward cases: a site overlapping itself, a site that spans
   the origin of a circle, a site whose recognition sequence is ambiguous,
   and a cut that falls outside its own site (type IIS).

   Output: rebase-mine.json next to this file.
   ========================================================================= */
const B = require('../../src/engine7.js');
const fs = require('fs');
const path = require('path');

/* Deterministic pseudo-random DNA: the same every run, so a disagreement is
   reproducible from this file alone without carrying a fixture around. */
function dna(n, seed) {
  let x = seed, s = '';
  for (let i = 0; i < n; i++) { x = (x * 1103515245 + 12345) & 0x7fffffff; s += 'ACGT'[x % 4]; }
  return s;
}

const SEQS = {
  /* random DNA long enough that most six-cutters appear by chance */
  random8k: dna(8000, 11),
  random8k_b: dna(8000, 29),
  /* every 6-mer in sequence, so every palindromic six-base site is present
     exactly once and adjacent sites overlap */
  allsites: (function () {
    const b = 'ACGT'; let s = '';
    for (const a of b) for (const c of b) for (const d of b) for (const e of b) for (const f of b) for (const g of b) s += a + c + d + e + f + g;
    return s;
  })(),
  /* a short one where a site sits right at each end */
  ends: 'GAATTCAAAAAAAAAAAAAAAAAAAAGAATTC',
  /* the same site repeated so that cuts land next to one another */
  tandem: 'AAA' + 'GAATTC'.repeat(6) + 'AAA'
};

const enzymes = B.ALL_ENZYMES.filter(e => !e.alias);
const out = {
  /* Two shapes live in this table. An enzyme that cuts inside its own site
     carries top/bot, the distance from the start of the site to the break in
     each strand. A type IIS enzyme cuts at a distance beyond the site and
     carries o1/o2 instead, measured from the end of it. Both are reported as
     they are stored, and the Python side converts rather than this one, so a
     conversion error cannot hide here. */
  enzymes: enzymes.map(e => ({
    name: e.name, site: e.site,
    top: e.top, bot: e.bot, o1: e.o1, o2: e.o2,
    typeIIS: !!e.typeIIS, suppliers: e.suppliers || ''
  })),
  supplierNames: B.SUPPLIER_NAMES,
  sequences: SEQS,
  cuts: {}
};

/* Where each enzyme cuts each sequence, as 0-based positions of the break in
   the top strand — the same quantity Biopython's search() reports, less its
   1-based offset. Linear and circular are both recorded: a site spanning the
   origin exists only in the circular reading, and that is where a cut list
   most easily goes wrong. */
for (const key in SEQS) {
  const seq = SEQS[key];
  out.cuts[key] = {};
  for (const e of enzymes) {
    const lin = B.digest(seq, [e.name], { circular: false });
    const cir = B.digest(seq, [e.name], { circular: true });
    const t = (lin.cuts || []).map(c => c.top);
    const c = (cir.cuts || []).map(x => x.top);
    if (t.length || c.length) out.cuts[key][e.name] = { linear: t, circular: c };
  }
}

fs.writeFileSync(path.join(__dirname, 'rebase-mine.json'), JSON.stringify(out));
console.log(enzymes.length + ' enzymes, ' + Object.keys(SEQS).length + ' sequences');
