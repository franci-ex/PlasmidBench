/* ===========================================================================
   Report what the tolerant matcher found, with the bases it aligned, for
   cmpanno.py to check against Biopython's own aligner.

   The identity a matcher reports is the one number a person reads off the
   screen and trusts, and it comes out of this project's own alignment. So
   it needs an oracle: Bio.Align.PairwiseAligner, aligning the same library
   entry against the same span of the same sequence, and asking whether it
   agrees about how similar they are.

   This half only reports. The expected values are computed in the Python
   half from the sequences, and nothing here tells it what to expect.
   ========================================================================= */
const B = require('../../src/engine8.js');
const fs = require('fs');

function stream(seed) {
  let a = seed >>> 0;
  return () => {
    a = (a + 0x6D2B79F5) >>> 0;
    let t = a;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}
const dnaFrom = r => n => { let s = ''; for (let i = 0; i < n; i++) s += 'ACGT'[Math.floor(r() * 4)]; return s; };

/* A library, and sequences built from it with known damage: substitutions,
   insertions, deletions, truncations, and combinations. */
const libR = stream(606), libDna = dnaFrom(libR);
const lib = [];
for (let i = 0; i < 60; i++) {
  lib.push({ name: 'entry' + i, type: 'misc_feature',
             seq: libDna([30, 48, 80, 150, 400][i % 5]) });
}
const idx = B.buildFeatureIndex(lib);

const cases = [];
const qR = stream(909), qDna = dnaFrom(qR);
const damage = {
  exact: s => s,
  subs2: s => s.split('').map((c, i) => i % 50 === 7 ? 'ACGT'[(('ACGT'.indexOf(c)) + 1) % 4] : c).join(''),
  subs5: s => s.split('').map((c, i) => i % 20 === 3 ? 'ACGT'[(('ACGT'.indexOf(c)) + 2) % 4] : c).join(''),
  del3: s => s.slice(0, Math.floor(s.length / 2)) + s.slice(Math.floor(s.length / 2) + 3),
  ins4: s => s.slice(0, Math.floor(s.length / 2)) + 'GTCA' + s.slice(Math.floor(s.length / 2)),
  trimFront: s => s.slice(Math.max(14, Math.floor(s.length * 0.2))),
  trimBack: s => s.slice(0, -Math.max(14, Math.floor(s.length * 0.2))),
  revcomp: s => B.revcomp(s)
};

for (const [label, f] of Object.entries(damage)) {
  for (let i = 0; i < 60; i += 7) {
    const entry = lib[i];
    const planted = f(entry.seq);
    const at = 600;
    const seq = qDna(at) + planted + qDna(700);
    const res = B.annotateApprox(seq, idx, {});
    const h = res.hits.find(x => x.name === entry.name);
    if (!h) { cases.push({ label: label, entry: entry.name, found: false,
                           entrySeq: entry.seq, planted: planted }); continue; }
    /* the bases the matcher says the feature occupies */
    const span = h.wrapped ? seq.slice(h.start) + seq.slice(0, h.end) : seq.slice(h.start, h.end);
    cases.push({
      label: label, entry: entry.name, found: true,
      strand: h.strand, identity: h.identity, coverage: h.coverage, gaps: h.gaps,
      exact: h.exact, partial: h.partial, how: h.how,
      entrySeq: entry.seq,
      /* on the strand the matcher reported, so the Python side aligns the
         same two strings and not a reverse complement of one of them */
      entryOnStrand: h.strand === 1 ? entry.seq : B.revcomp(entry.seq),
      span: span,
      start: h.start, end: h.end
    });
  }
}

fs.writeFileSync(process.argv[2] || '/dev/stdout', JSON.stringify(cases));
console.error(cases.length + ' cases, ' + cases.filter(c => c.found).length + ' found');
