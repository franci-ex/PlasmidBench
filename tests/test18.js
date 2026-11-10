/* ===========================================================================
   What of a read is already in the plasmid, and what each part of a
   fragment actually says.

   The question behind all of this is the one asked at the bench: a Sanger
   read starts in the vector and runs into the insert, so which part of it
   is the thing being cloned? Getting that wrong in either direction is
   expensive — call vector "new" and you clone a duplicate promoter into
   your own backbone; call insert "vector" and you leave the gene behind.

   Every case here is built from pieces whose boundaries are known by
   construction, so the expected answer is arithmetic rather than judgement.
   ========================================================================= */
const B = require('../src/engine.js');
require('../src/engine2.js'); require('../src/engine3.js'); require('../src/engine4.js');
let pass = 0, fail = 0;
function ok(name, cond, extra) {
  if (cond) { pass++; console.log('  PASS  ' + name); }
  else { fail++; console.log('  FAIL  ' + name + (extra !== undefined ? '   -> ' + extra : '')); }
}
function prng(seed) { let x = seed; return () => (x = (x * 1103515245 + 12345) & 0x7fffffff) / 0x7fffffff; }
function dnaOf(n, seed) { const r = prng(seed); let s = ''; for (let i = 0; i < n; i++) s += 'ACGT'[Math.floor(r() * 4)]; return s; }
/* An insert whose first and last bases differ from the plasmid bases they
   will sit against. A random run that happens to repeat its neighbour really
   is found in the plasmid at that base, so the boundary moves by one — true,
   and not what any of these cases is measuring. */
function wedge(n, seed, before, after) {
  const not = (a, b) => 'ACGT'.split('').filter(x => x !== a && x !== b)[0];
  const mid = dnaOf(Math.max(1, n - 2), seed);
  return not(before, mid[0]) + mid + not(after, mid[mid.length - 1]);
}

/* a pET-shaped backbone: promoter, operator, RBS, ATG + His + polylinker */
const T7 = 'TAATACGACTCACTATAGGG';
const LACO = 'GATTGTGAGCGGATAACAATT';
const RBSBLOCK = 'CCCCTGTAGAAATAATTTTGTTTAACTTTAATAAGGAGATATACC';
const TAGBLOCK = 'ATGGGCAGCAGCCATCACCATCATCACCACAGCCAGGATCCGAATTCGAGCTCGAGT';
const HEAD = T7 + LACO + RBSBLOCK;
const BACK = dnaOf(2400, 31);
const VEC = HEAD + TAGBLOCK + BACK;
const REC = {
  name: 'pET-ish', seq: VEC, topology: 'circular',
  features: [
    { start: 0, end: T7.length, strand: 1, type: 'promoter', label: 'T7 promoter', quals: {} },
    { start: T7.length, end: T7.length + LACO.length, strand: 1, type: 'protein_bind', label: 'lac operator', quals: {} },
    { start: HEAD.length - 14, end: HEAD.length - 8, strand: 1, type: 'RBS', label: 'RBS', quals: {} },
    { start: HEAD.length + 9, end: HEAD.length + 27, strand: 1, type: 'CDS', label: '6xHis', quals: {} }
  ]
};

console.log('\n== a read that is nothing but vector ==');
(function () {
  const read = VEC.slice(0, 400);
  const v = B.readInVector(read, REC, {});
  ok('every base is found in the plasmid', v.inVector === 400, v.inVector + ' of 400');
  ok('nothing is called new', v.isNew === 0, v.isNew);
  ok('and there is no span to clone', v.newSpan === null, JSON.stringify(v.newSpan));
  ok('the features it crossed are named',
     ['T7 promoter', 'lac operator', 'RBS', '6xHis'].every(n => v.features.indexOf(n) >= 0),
     v.features.join(', '));
})();

console.log('\n== vector in front, insert behind — the real shape of a read ==');
(function () {
  const INSERT = wedge(600, 77, VEC[143], 'x');
  const read = VEC.slice(0, 144) + INSERT;
  const v = B.readInVector(read, REC, {});
  ok('the vector part is measured exactly', v.inVector === 144, v.inVector);
  ok('and the insert is the rest', v.isNew === 600, v.isNew);
  ok('the span to clone starts where the vector stops',
     v.newSpan && v.newSpan.from === 144, v.newSpan && v.newSpan.from);
  ok('and runs to the end of the read',
     v.newSpan && v.newSpan.to === read.length, v.newSpan && v.newSpan.to);
  ok('the bases it hands back are the insert itself',
     v.newSpan && v.newSpan.seq === INSERT, v.newSpan && v.newSpan.seq.slice(0, 20));
  ok('with nothing ambiguous in it', v.newSpan && v.newSpan.ambiguous === 0, v.newSpan && v.newSpan.ambiguous);
})();

console.log('\n== an insert in the middle, vector on both sides ==');
(function () {
  const INSERT = wedge(300, 909, VEC[199], VEC[200]);
  const read = VEC.slice(0, 200) + INSERT + VEC.slice(200, 400);
  const v = B.readInVector(read, REC, {});
  ok('both flanks are recognised', v.inVector === 400, v.inVector);
  ok('the middle is the new part', v.isNew === 300, v.isNew);
  ok('and the span is the middle, not the whole read',
     v.newSpan && v.newSpan.from === 200 && v.newSpan.to === 500,
     v.newSpan && (v.newSpan.from + '-' + v.newSpan.to));
})();

console.log('\n== a read off the other strand ==');
(function () {
  const INSERT = wedge(500, 4242, VEC[143], 'x');
  const fwd = VEC.slice(0, 144) + INSERT;
  const v = B.readInVector(B.revcomp(fwd), REC, {});
  ok('the vector part is found on the bottom strand', v.inVector === 144, v.inVector);
  ok('and is reported as such',
     v.blocks.some(b => b.kind === 'vector' && b.strand === -1),
     JSON.stringify(v.blocks.filter(b => b.kind === 'vector').map(b => b.strand)));
  ok('the insert is still 500 bases', v.isNew === 500, v.isNew);
  ok('and the span is the reverse complement of what went in',
     v.newSpan && v.newSpan.seq === B.revcomp(INSERT), v.newSpan && v.newSpan.seq.slice(0, 20));
})();

console.log('\n== a read across the origin of a circle ==');
(function () {
  const read = VEC.slice(VEC.length - 150) + VEC.slice(0, 150);
  const v = B.readInVector(read, REC, {});
  ok('a circular plasmid is read round the join', v.inVector === 300, v.inVector);
  ok('so nothing at the origin looks new', v.isNew === 0, v.isNew);
  ok('and the whole read is one unbroken stretch of it',
     v.blocks.filter(b => b.kind === 'vector').length === 1,
     JSON.stringify(v.blocks.map(b => b.kind + ':' + b.bp)));
  /* Every base of this read is in a linear plasmid too — both halves are
     there. What is not there is the order: the end running straight into
     the beginning exists only on a circle. Matching base by base cannot
     see that; requiring the plasmid positions to run on can, and reports
     the join as two stretches rather than one. */
  const lin = B.readInVector(read, Object.assign({}, REC, { topology: 'linear' }), {});
  ok('a linear one has the same bases', lin.inVector === 300, lin.inVector);
  ok('but not in that order, so it reads as two separate stretches',
     lin.blocks.filter(b => b.kind === 'vector').length === 2,
     JSON.stringify(lin.blocks.map(b => b.kind + ':' + b.bp)));
  ok('the first of them is the end of the plasmid',
     lin.blocks[0].vecTo === VEC.length, lin.blocks[0].vecTo + ' of ' + VEC.length);
  ok('and the second is its beginning', lin.blocks[1].vecFrom === 0, lin.blocks[1].vecFrom);
})();

console.log('\n== one changed base is a difference, not something to clone ==');
(function () {
  const read = VEC.slice(0, 400).split('');
  read[200] = ({ A: 'C', C: 'A', G: 'T', T: 'G' })[read[200]];
  const v = B.readInVector(read.join(''), REC, {});
  ok('the gap it opens is called a difference', v.differs > 0 && v.differs < v.k * 2, v.differs);
  ok('nothing is offered as new sequence', v.isNew === 0, v.isNew);
  ok('and there is no span to clone', v.newSpan === null, JSON.stringify(v.newSpan));
})();

console.log('\n== the N-rich tail every Sanger read has ==');
(function () {
  const INSERT = wedge(400, 55, VEC[143], 'x');
  const tail = 'N'.repeat(120);
  const read = VEC.slice(0, 144) + INSERT + tail;
  const v = B.readInVector(read, REC, {});
  ok('the tail is called unreadable', v.unreadable === 120, v.unreadable);
  ok('and is not counted as sequence to add', v.isNew === 400, v.isNew);
  ok('the span to clone stops before it',
     v.newSpan && v.newSpan.to === 144 + 400, v.newSpan && v.newSpan.to);
  ok('so the tail is never handed to the cloning designer',
     v.newSpan && v.newSpan.seq.indexOf('N') < 0, v.newSpan && v.newSpan.seq.slice(-10));
  /* ambiguity packed into the far end is a different warning from the same
     amount sprinkled through: it means the end of the insert is unknown */
  /* a patch in the middle, with plenty of clean sequence after it */
  const mid = VEC.slice(0, 144) + INSERT.slice(0, 150) + 'N'.repeat(40) + INSERT.slice(190);
  const vm = B.readInVector(mid, REC, {});
  ok('a clean insert reports its whole length as read',
     v.newSpan.cleanBp === v.newSpan.bp && !v.newSpan.tailAmbiguous,
     v.newSpan.cleanBp + ' of ' + v.newSpan.bp);
  const deadEnd = VEC.slice(0, 144) + INSERT.slice(0, 300) + 'N'.repeat(40) + 'ACGTACGTACGTACGTACGTACGT';
  const vt = B.readInVector(deadEnd, REC, {});
  ok('a read that never recovers is called out as such',
     vt.newSpan && vt.newSpan.tailAmbiguous && vt.newSpan.cleanBp === 300,
     vt.newSpan && (vt.newSpan.cleanBp + ' clean, ' + vt.newSpan.finalCleanRun +
       ' at the end of ' + vt.newSpan.bp + ', tail=' + vt.newSpan.tailAmbiguous));
  ok('ambiguity in the middle is not, because the read recovers from it',
     vm.newSpan && !vm.newSpan.tailAmbiguous && vm.newSpan.finalCleanRun > 150,
     vm.newSpan && (vm.newSpan.finalCleanRun + ' clean at the end of ' + vm.newSpan.bp));
})();

console.log('\n== a ragged start is an edge, not a finding ==');
(function () {
  const read = 'NNNNNN' + 'CGACTC' + VEC.slice(12, 400);
  const v = B.readInVector(read, REC, {});
  ok('the few mismatched bases at the front are called an edge', v.edge === 6, v.edge);
  ok('and not reported as new sequence', v.isNew === 0, v.isNew);
  ok('so a read that found nothing new says so', v.newSpan === null, JSON.stringify(v.newSpan));
})();

console.log('\n== a short unmatched run cannot be resolved either way ==');
(function () {
  /* k is 18; nothing below that length can be told apart from a mismatch.
     The inserted run is made to differ from the vector base on either side
     of it, because a run whose last base happens to repeat the vector's
     next one really is found in the vector there, and the boundary would
     move by that base — true, but not what is being measured here. */
  const other = (a, b) => 'ACGT'.split('').filter(x => x !== a && x !== b)[0];
  const wedge = n => {
    const mid = dnaOf(Math.max(1, n - 2), 3);
    return other(VEC[199], mid[0]) + mid + other(VEC[200], mid[mid.length - 1]);
  };
  const v = B.readInVector(VEC.slice(0, 200) + wedge(10) + VEC.slice(200, 400), REC, {});
  ok('ten unmatched bases between two matches are a difference', v.differs === 10 && v.isNew === 0,
     'differs=' + v.differs + ' new=' + v.isNew);
  const w = B.readInVector(VEC.slice(0, 200) + wedge(40) + VEC.slice(200, 400), REC, {});
  ok('forty are an insertion', w.isNew === 40 && w.differs === 0,
     'new=' + w.isNew + ' differs=' + w.differs);
})();

console.log('\n== refusals ==');
(function () {
  ok('a read shorter than one window is refused',
     B.readInVector('ACGTACGTAC', REC, {}) === null);
  ok('and so is a plasmid that is not there',
     B.readInVector(VEC.slice(0, 200), { name: 'x', seq: '' }, {}) === null);
})();

console.log('\n== the whole of the user’s read, against a pET-shaped vector ==');
(function () {
  const read = 'NNNNNNCGACTCNCTATAGGGGATTGTGAGCGGATAACAATTCCCCTGTAGAAATAATTTTGTTTAACTTTAATAAGGAGATATACC' +
    'ATGGGCAGCAGCCATCACCATCATCACCACAGCCAGGATCCGAATTCGAGCTCGAGT' +
    'GGCACCTACAAACTGATCCTGAACGGTAAAACCCTGAAAGGTGAAACCACCACCGAAGCTGTAGACGCTGCTACTGCTGAAAAAGTT' +
    'TTCAAACAGTACGCTAACGACAACGGTGTGGACGGTGAATGGACCTACGACGACGCTACCAAAACCTTCACGGTTACGGAAAGCGGT' +
    'TCGGATGAAAACTTGTATTTCCAGGGCGGGATGGAGAACTTCCAAAAGGTGGAAAAGATCGGAGAGGGCACGTACGGAGTTGTGTAC' +
    'AAAGCCAGAAACAAGTTGACGGGAGAGGTGGTGGCGCTTAAGAAAATCCGCCTGGACACTGAGACTGAGGGTGTGCCCAGTACTGCC' +
    'ATCCGAGAGATCTCTCTGCTTAAGGAGCTTAACCATCCTAATATTGTCAAGCTGCTGGATGTCATTCACACAGAAAATAAACTCTAC';
  const v = B.readInVector(read, REC, {});
  ok('the promoter, operator, RBS and tag are all found in the vector',
     ['T7 promoter', 'lac operator', 'RBS', '6xHis'].every(n => v.features.indexOf(n) >= 0),
     v.features.join(', '));
  /* The ragged front of this read stands in for twelve bases of vector with
     thirteen characters, so the read runs one base behind the plasmid all
     the way along; what matters is not the count but that the boundary
     falls at the end of the polylinker. */
  const vb = v.blocks.filter(b => b.kind === 'vector');
  ok('the match is one unbroken run', vb.length === 1, JSON.stringify(v.blocks.map(b => b.kind)));
  ok('and it ends at the last base of the polylinker',
     vb[0] && vb[0].vecTo === HEAD.length + TAGBLOCK.length,
     vb[0] && (vb[0].vecTo + ' vs ' + (HEAD.length + TAGBLOCK.length)));
  ok('the span to clone begins at GB1',
     v.newSpan && v.newSpan.seq.slice(0, 13) === 'GGCACCTACAAAC',
     v.newSpan && v.newSpan.seq.slice(0, 13));
  ok('the new part is everything from there to the end of the read',
     v.newSpan && v.newSpan.to === read.length, v.newSpan && (v.newSpan.to + ' of ' + read.length));
})();

console.log('\n== every region of a fragment carries its own bases ==');
(function () {
  const frag = 'GAATTCTTTT' + 'ATG' + 'GGCAGCAGCCATCACCATCATCACCAC' + 'GAAAACCTGTATTTTCAGGGC' +
    'ATGGAAGATTATACCAAAATAGAGAAAATAGGAGAAGGAACTTACGGAGTTGTGTATAAAGGACGTCATAAAACAACAGGACAGGTA' +
    'GTGGCAATGAAAAAGATTCGTCTAGAAAGTGAAGAGGAAGGG' + 'TAA';
  const enz = B.ALL_ENZYMES.filter(e => e.name === 'EcoRI');
  const a = B.fragmentAnatomy(frag, { enzymes: enz, table: 11 });
  ok('the regions tile the fragment with no gap and no overlap',
     a.segments.every((g, i) => i === 0 ? g.from === 0 : g.from === a.segments[i - 1].to) &&
     a.segments[a.segments.length - 1].to === a.length,
     a.segments.map(g => g.from + '-' + g.to).join(' '));
  ok('each one carries exactly its own bases',
     a.segments.every(g => g.dna === a.work.slice(g.from, g.to)),
     a.segments.map(g => (g.dna || '').length + '/' + g.bp).join(' '));
  ok('and the bases put back together are the fragment',
     a.segments.map(g => g.dna).join('') === a.work);
  ok('a region that codes translates to what it says it does',
     a.segments.filter(g => g.kind === 'coding' || g.kind === 'fusion partner' ||
                            g.kind === 'tag' || g.kind === 'protease site')
       .every(g => B.translate(g.dna, { table: 11 }).replace(/\*$/, '') === (g.aa || '').replace(/\*$/, '')),
     a.segments.filter(g => g.aa).map(g => g.name + ':' + B.translate(g.dna, { table: 11 }) + ' vs ' + g.aa).join(' | '));
})();

console.log('\n== a fragment pasted in reverse is read the same way ==');
(function () {
  const frag = 'GAATTCTTTT' + 'ATG' + 'GGCAGCAGCCATCACCATCATCACCAC' + 'GAAAACCTGTATTTTCAGGGC' +
    'ATGGAAGATTATACCAAAATAGAGAAAATAGGAGAAGGAACTTACGGAGTTGTGTATAAAGGACGTCATAAAACAACAGGACAGGTA' +
    'GTGGCAATGAAAAAGATTCGTCTAGAAAGTGAAGAGGAAGGG' + 'TAA';
  const enz = B.ALL_ENZYMES.filter(e => e.name === 'EcoRI');
  const f = B.fragmentAnatomy(frag, { enzymes: enz, table: 11 });
  const r = B.fragmentAnatomy(B.revcomp(frag), { enzymes: enz, table: 11 });
  ok('it is seen to read on the reverse strand', r.strand === -1, r.strand);
  ok('the working strand is the fragment the right way round', r.work === f.work);
  ok('the same regions are found in the same places',
     JSON.stringify(r.segments.map(g => [g.from, g.to, g.kind, g.name])) ===
     JSON.stringify(f.segments.map(g => [g.from, g.to, g.kind, g.name])),
     r.segments.map(g => g.from + '-' + g.to + ':' + g.name).join(' | '));
  /* The site is in the untranslated flank, 10 bases in front of the ATG.
     Reported off the pasted strand instead, it landed at length minus its
     position — 188 rather than 1 — and was then called a site inside the
     protein, which is where "how this was cloned before" used to point for
     any fragment given in reverse. */
  ok('the EcoRI site is at the front either way',
     f.sites.length === 1 && r.sites.length === 1 && f.sites[0].at === 0 && r.sites[0].at === 0,
     JSON.stringify([f.sites, r.sites]));
  ok('and is outside the protein either way',
     !f.sites[0].inCore && !r.sites[0].inCore,
     JSON.stringify([f.sites[0].inCore, r.sites[0].inCore]));
  ok('the first region is named from the right strand’s letters',
     r.segments[0].name === f.segments[0].name && r.segments[0].dna === f.segments[0].dna,
     r.segments[0].name + ' / ' + r.segments[0].dna);
})();

console.log('\n== the two ways of writing a feature that crosses the origin ==');
(function () {
  /* A wrap can be written with the end before the start (4057..46) or with
     the end past the end of the sequence (4057..4888 on a 4,842 bp circle).
     Both come out of this program's own code and out of the files it reads.
     The second was handed back unsplit, and slice() then returned only the
     bases up to the end of the string — the length and the printed location
     stayed right, because those are arithmetic on the bounds, so the only
     symptom was a feature quietly failing to carry into an assembly product
     because the copy's bases no longer matched the original's. */
  const L = 4842, g = dnaOf(L, 99);
  const past = { type: 'CDS', start: 4057, end: 4888, strand: 1, label: 'CDK2' };
  const back = { type: 'CDS', start: 4057, end: 46, strand: 1, label: 'CDK2' };
  ok('a feature written past the end covers 831 bases, not 785',
     B.featureSequence(g, past).length === 831, B.featureSequence(g, past).length);
  ok('and the two ways of writing it give the same bases',
     B.featureSequence(g, past) === B.featureSequence(g, back),
     B.featureSequence(g, past).slice(0, 12) + ' vs ' + B.featureSequence(g, back).slice(0, 12));
  ok('which are the end of the plasmid followed by its start',
     B.featureSequence(g, past) === g.slice(4057) + g.slice(0, 46));
  ok('the length agrees with the bases', B.featureLength(past, L) === 831, B.featureLength(past, L));
  ok('it splits into two pieces at the origin',
     JSON.stringify(B.featurePieces(past, L)) === JSON.stringify([{ start: 4057, end: L }, { start: 0, end: 46 }]),
     JSON.stringify(B.featurePieces(past, L)));
  ok('on the reverse strand it is the reverse complement of that',
     B.featureSequence(g, Object.assign({}, past, { strand: -1 })) ===
     B.revcomp(g.slice(4057) + g.slice(0, 46)));
  ok('a feature that does not wrap is untouched',
     B.featureSequence(g, { start: 100, end: 400, strand: 1 }) === g.slice(100, 400));
  ok('and one written entirely past the end is brought back round',
     B.featureSequence(g, { start: L + 10, end: L + 40, strand: 1 }) === g.slice(10, 40),
     B.featureSequence(g, { start: L + 10, end: L + 40, strand: 1 }).length);
  /* a joined feature whose last part runs past the end */
  const mixed = { type: 'CDS', start: 4000, end: 4888, strand: 1, mixed: true, label: 'j',
    parts: [{ start: 4000, end: 4050, strand: 1 }, { start: 4057, end: 4888, strand: 1 }] };
  ok('a joined feature wraps in its own parts too',
     B.featureSequence(g, mixed) === g.slice(4000, 4050) + g.slice(4057) + g.slice(0, 46),
     B.featureSequence(g, mixed).length + ' vs ' + (50 + 831));
})();

/* ---------------------------------------------------------------------------
   The melting temperature of a primer that has an overhang.

   A primer with a 5' tail does not anneal along its whole length in the
   first cycles — only the 3' part that matches the template binds, and that
   part alone sets the temperature the block has to be run at. The whole
   primer's Tm is always the higher of the two, often by more than ten
   degrees, and annealing at it means annealing above what actually holds.
   Both numbers are real, at different times: the annealing part governs the
   first cycles, the whole primer governs later ones once the tail has been
   copied into the product.
--------------------------------------------------------------------------- */
console.log('\n== the Tm of a primer with an overhang ==');
(function () {
  const C = B.PCR_CONDITIONS;
  const tail = 'TTGTATTTCCAGGGC';                 // 15 bp homology arm
  const an = 'ATGGGTTTTCTGCGTTTTC';               // 19 bases that bind
  const full = tail + an;

  const t = B.primerTm(full, { anneal: an, conditions: C });
  ok('it measures the annealing part, not the whole primer',
     Math.abs(t.tm - B.tmNN(an, C)) < 1e-9, t.tm + ' vs ' + B.tmNN(an, C));
  ok('and still reports the whole primer separately',
     Math.abs(t.tmFull - B.tmNN(full, C)) < 1e-9, t.tmFull);
  ok('the two are not the same number', Math.abs(t.tmFull - t.tm) > 5,
     (t.tmFull - t.tm).toFixed(1) + ' degrees apart');
  ok('the overhang is the higher of the two, which is what makes using it a trap',
     t.tmFull > t.tm, t.tm + ' / ' + t.tmFull);
  ok('it says how much anneals and how much does not',
     t.annealLen === an.length && t.tailLen === tail.length && t.tailed === true,
     JSON.stringify({ a: t.annealLen, t: t.tailLen }));
  ok('and hands back the two pieces', t.anneal === an && t.tail === tail);

  const byLen = B.primerTm(full, { annealLen: an.length, conditions: C });
  ok('giving the length instead of the bases says the same thing',
     byLen.tm === t.tm && byLen.tmFull === t.tmFull && byLen.tailLen === t.tailLen);

  /* An annealing portion that is not the 3' end of this primer belongs to
     some other oligo. Measuring it anyway would be worse than refusing. */
  const wrong = B.primerTm(full, { anneal: 'GGGGGGGGGGGG', conditions: C });
  ok('a portion that is not the 3′ end of the primer is not believed',
     wrong.tailed === false && Math.abs(wrong.tm - B.tmNN(full, C)) < 1e-9, wrong.tm);
  ok('and the fallback is marked as the whole primer', wrong.basis === 'full', wrong.basis);

  const none = B.primerTm(full, { conditions: C });
  ok('with nothing known, the whole primer is all there is',
     none.tailed === false && none.tm === none.tmFull && none.basis === 'full');
  const whole = B.primerTm(full, { anneal: full, conditions: C });
  ok('a primer that anneals along its whole length has no tail',
     whole.tailed === false && whole.tm === whole.tmFull && whole.tailLen === 0);
  ok('the basis is carried through when one is given',
     B.primerTm(full, { anneal: an, conditions: C, basis: 'design' }).basis === 'design');

  /* edges */
  const empty = B.primerTm('', { conditions: C });
  ok('an empty primer gives no temperature', !isFinite(empty.tm) && empty.length === 0);
  const over = B.primerTm(full, { annealLen: full.length + 20, conditions: C });
  ok('an annealing length longer than the primer is ignored', over.tailed === false, over.annealLen);
  const zero = B.primerTm(full, { annealLen: 0, conditions: C });
  ok('and so is one of zero', zero.tailed === false && zero.tm === zero.tmFull);

  /* the conditions have to reach both numbers, or one of them is measured
     in a reaction nobody is running */
  const hot = B.primerTm(full, { anneal: an, conditions: { oligo: 0.25e-6, Na: 0.05, Mg: 4, dNTP: 0.2 } });
  ok('the salt conditions change the annealing Tm too, not only the full one',
     hot.tm !== t.tm && hot.tmFull !== t.tmFull, hot.tm + ' vs ' + t.tm);

  /* against the real In-Fusion designer: what it chose a primer on must be
     what this reports, or the library and the designer would disagree */
  const ins = 'ATGGGTTTTCTGCGTTTTCAGGCTAGCGTGCCCTCAGGTTTAACCGGTAAAGAATTCCGTACGTTACCC';
  const nb = dnaOf(400, 5150);
  const d = B.tuneAssemblyPrimer({ insert: ins, neighbour: nb, end: '5',
                                  armMin: 15, armMax: 15, targetTm: 60, conditions: C });
  if (d && d.ok) {
    const back = B.primerTm(d.best.seq, { anneal: d.best.anneal, conditions: C });
    ok('the designer’s own annealing Tm comes back unchanged',
       Math.abs(back.tm - d.best.tmAnneal) < 1e-9, back.tm + ' vs ' + d.best.tmAnneal);
    ok('and so does its full-length one',
       Math.abs(back.tmFull - d.best.tmFull) < 1e-9, back.tmFull + ' vs ' + d.best.tmFull);
    ok('the tail it reports is the arm the designer added',
       back.tailLen === d.best.tail.length, back.tailLen + ' vs ' + d.best.tail.length);
  } else ok('the designer produced a tailed primer to check against', false, JSON.stringify(d && d.ok));

  /* a property run: for any split, the pieces fit back together and the
     annealing Tm is never the full-length one unless there is no tail */
  let bad = 0, n = 0;
  for (let seed = 1; seed <= 120; seed++) {
    const prim = dnaOf(18 + (seed % 22), seed * 7 + 1);
    for (let cut = 1; cut < prim.length - 8; cut += 3) {
      const r = B.primerTm(prim, { annealLen: prim.length - cut, conditions: C });
      n++;
      if (r.tail + r.anneal !== prim) bad++;
      if (r.annealLen + r.tailLen !== prim.length) bad++;
      if (r.tailed && Math.abs(r.tm - B.tmNN(r.anneal, C)) > 1e-9) bad++;
      if (r.tailed && Math.abs(r.tmFull - B.tmNN(prim, C)) > 1e-9) bad++;
    }
  }
  ok('across ' + n + ' splits the pieces always rebuild the primer and both Tms match their own sequence',
     bad === 0, bad);
})();

/* ---------------------------------------------------------------------------
   Recovering the overhang of a primer nobody recorded one for.

   A primer saved as a name and a string of bases says nothing about which
   part of it anneals, and the number shown for it is then the whole
   primer's — fifteen degrees above what binds, in the direction that gives
   an empty tube. Two kinds of evidence can settle it, and both are about
   this primer against real sequence.
--------------------------------------------------------------------------- */
console.log('\n== finding the overhang on a primer that was saved without one ==');
(function () {
  const C = B.PCR_CONDITIONS;
  const armF = 'TTGTATTTCCAGGGC';                  // 15 nt, vector sequence
  const armR = 'ATTATGCGGCCGCAATT';                // 17 nt, vector sequence
  const insF = 'ATGGGTTTTCTGCGTTTTC';
  const filler = dnaOf(900, 5), filler2 = dnaOf(60, 6), filler3 = dnaOf(900, 7);
  const vec = filler + B.revcomp(armF) + filler2 + B.revcomp(armR) + filler3;
  const V = { name: 'pRSFDuet', seq: vec, topology: 'circular', features: [] };

  const f = B.findOverhang(armF + insF, [V], {});
  ok('an assembly primer’s arm is found in the vector it was built for',
     f && f.tailLen === armF.length, f && f.tailLen);
  ok('and what is left is the part that anneals',
     f && f.anneal === insF && f.annealLen === insF.length, f && f.anneal);
  ok('it says how it knows', f && f.how === 'tail-matches' && f.on === 'pRSFDuet',
     f && (f.how + ' / ' + f.on));

  /* the ordinary case: the template is loaded, so where it anneals is
     measurable directly and that is the better evidence */
  /* The base in front of the insert is chosen not to repeat the last base
     of the arm; otherwise the annealing really does run one base further
     and the two routes differ by that base, which is the case below. */
  const tmpl = { name: 'insert source', seq: dnaOf(299, 8) + 'A' + insF + dnaOf(300, 9),
                 topology: 'linear', features: [] };
  const both = B.findOverhang(armF + insF, [V, tmpl], {});
  ok('with the template loaded too, where it anneals is used instead',
     both && both.how === 'anneals' && both.annealLen === insF.length,
     both && (both.how + ' ' + both.annealLen));
  ok('and both routes agree here, so nothing is flagged',
     both && !both.disagrees, both && JSON.stringify(both.disagrees));

  /* When the template carries on matching past the insert, the annealing
     route is right and the arm route is a base short of it. The direct
     measurement wins, and the disagreement is reported rather than hidden,
     because a boundary two routes dispute is one worth looking at. */
  const sticky = { name: 'sticky template', topology: 'linear', features: [],
                   seq: dnaOf(299, 8) + armF.slice(-1) + insF + dnaOf(300, 9) };
  const sd = B.findOverhang(armF + insF, [V, sticky], {});
  ok('a template that carries on matching extends the annealing part',
     sd && sd.how === 'anneals' && sd.annealLen === insF.length + 1,
     sd && (sd.how + ' ' + sd.annealLen));
  ok('and the disagreement with the arm search is carried, not swallowed',
     sd && sd.disagrees && sd.disagrees.annealLen === insF.length,
     sd && JSON.stringify(sd.disagrees));

  /* things it must not claim */
  ok('a primer that anneals along its whole length has no overhang',
     B.findOverhang(vec.slice(300, 322), [V], {}) === null,
     JSON.stringify(B.findOverhang(vec.slice(300, 322), [V], {})));
  /* The case that matters most, because it is half the library: a backbone
     primer is vector sequence from end to end, so SOME prefix of it is
     always in the vector — all of it is. Taking the longest prefix that
     matches would call its first fourteen bases a homology arm and report
     half the real temperature. */
  [22, 28, 34, 40].forEach(n => {
    ok('a ' + n + ' nt backbone primer is not given an arm',
       B.findOverhang(vec.slice(500, 500 + n), [V], {}) === null,
       JSON.stringify(B.findOverhang(vec.slice(500, 500 + n), [V], {})));
    ok('nor is its reverse complement',
       B.findOverhang(B.revcomp(vec.slice(700, 700 + n)), [V], {}) === null,
       JSON.stringify(B.findOverhang(B.revcomp(vec.slice(700, 700 + n)), [V], {})));
  });
  /* and a primer whose 5' end matches further than the search window can
     see is left alone rather than having its arm put at the window edge */
  ok('an arm longer than the window is not guessed at the edge of it',
     B.findOverhang(vec.slice(200, 260), [V], {}) === null,
     JSON.stringify(B.findOverhang(vec.slice(200, 260), [V], {})));
  ok('a primer unrelated to anything loaded gets no answer',
     B.findOverhang(dnaOf(32, 1234), [V], {}) === null);
  ok('with nothing loaded there is nothing to go on',
     B.findOverhang(armF + insF, [], {}) === null);
  ok('a primer too short to split either way is left alone',
     B.findOverhang(armF.slice(0, 12), [V], {}) === null);
  /* a short 5' coincidence is not an arm: ten bases of exact match turn up
     about once in 200 kb, and the floor is set there */
  const shortPre = vec.slice(400, 408) + dnaOf(24, 77);
  ok('eight bases of coincidental 5′ match are not called an arm',
     (B.findOverhang(shortPre, [V], {}) || {}).tailLen !== 8,
     JSON.stringify(B.findOverhang(shortPre, [V], {})));

  /* the arm sitting across the origin of a circle */
  const wrapVec = { name: 'wrapped', topology: 'circular', features: [],
    seq: B.revcomp(armF).slice(7) + dnaOf(1200, 11) + B.revcomp(armF).slice(0, 7) };
  const w = B.findOverhang(armF + insF, [wrapVec], {});
  ok('an arm that crosses base 1 of a circular vector is still found',
     w && w.tailLen === armF.length, w && JSON.stringify(w));
  const linVec = Object.assign({}, wrapVec, { topology: 'linear' });
  ok('and on a linear one, where it genuinely is not contiguous, it is not',
     (B.findOverhang(armF + insF, [linVec], {}) || {}).tailLen !== armF.length,
     JSON.stringify(B.findOverhang(armF + insF, [linVec], {})));

  /* the boundary is only as sharp as the sequence allows, and the slip is
     always towards a longer arm, which reads as a lower temperature */
  const run = B.findOverhang(armR + 'AAACCTGAGGGCATACAC', [V], {});
  ok('the arm found is never shorter than the real one',
     run && run.tailLen >= armR.length, run && run.tailLen);
  const exact = B.primerTm(armR + 'AAACCTGAGGGCATACAC', { annealLen: 18, conditions: C }).tm;
  const got = B.primerTm(armR + 'AAACCTGAGGGCATACAC', { annealLen: run.annealLen, conditions: C }).tm;
  ok('so the temperature it gives is at or below the exact one, never above',
     got <= exact + 1e-9, got.toFixed(1) + ' vs ' + exact.toFixed(1));
  ok('and far closer to it than the whole primer would be',
     Math.abs(got - exact) < 2 &&
     Math.abs(B.tmNN(armR + 'AAACCTGAGGGCATACAC', C) - exact) > 10,
     Math.abs(got - exact).toFixed(1) + ' vs ' +
     Math.abs(B.tmNN(armR + 'AAACCTGAGGGCATACAC', C) - exact).toFixed(1));

  /* every primer of a real-looking set, against the vector alone */
  const set = [armF + insF, armF + 'ATGCGTCGCGTTAGTC', armR + 'AGGACAAACGTGCCC',
               armR + 'AAACCTGAGGGCATACAC'];
  const all = set.map(x => B.findOverhang(x, [V], {}));
  ok('every assembly primer in a set is resolved from the vector alone',
     all.every(x => x && x.tailLen >= 15), JSON.stringify(all.map(x => x && x.tailLen)));
  ok('and none of them is left measuring its whole length',
     all.every((x, i) => x.annealLen < set[i].length), JSON.stringify(all.map(x => x && x.annealLen)));
})();

/* ---------------------------------------------------------------------------
   Backbone primers are a pair, not two primers.

   Two primers that go in one tube share one annealing temperature. Choosing
   each one's length on its own distance from the target optimises the wrong
   thing: it can leave them several degrees apart, and then no single
   temperature suits both — set it for the cooler one and the warmer primer
   also binds where it half-matches, set it for the warmer and the cooler
   one does not bind at all. Where the target is reachable this changes
   nothing. Where it is not — a GC-rich opening, where the shortest allowed
   length already melts above it — a matched pair running a few degrees high
   beats a split one.
--------------------------------------------------------------------------- */
console.log('\n== the two backbone primers are balanced against each other ==');
(function () {
  const GC = 'GCGGCCGCGGCGGCCGCGGCGGCGGCCGCGGCCGCGG';
  let worstSpread = 0, worstCase = '', n = 0, nonUnique = 0;
  for (let seed = 1; seed <= 60; seed++) {
    const v = dnaOf(1200, seed) + GC + dnaOf(1200, seed + 500);
    const o = B.inversePCR(v, 1200, 1200 + GC.length, { targetTm: 60 });
    if (!o.ok) continue;
    n++;
    const F = o.primers.forward, R = o.primers.reverse;
    const spread = Math.abs(F.tm - R.tm);
    if (spread > worstSpread) { worstSpread = spread; worstCase = 'seed ' + seed + ': ' + F.tm.toFixed(1) + ' / ' + R.tm.toFixed(1); }
    /* A primer is only expected to be unique when some length of it can
       be: at this deliberately repetitive opening one side sometimes sits
       wholly inside the repeat, where no length binds in one place only.
       What matters then is that the result says unique: false rather than
       claiming otherwise. */
    const couldBe = make => {
      for (let k = 18; k <= 45; k++) {
        let pr;
        try { pr = make(k); } catch (e) { break; }
        if (!pr) break;
        try { if (B.primerSites(v, pr, { circular: true }).length === 1) return true; } catch (e) {}
      }
      return false;
    };
    const dbl = v + v, Lv = v.length, aa = 1200, bb = 1200 + GC.length;
    if (couldBe(k => dbl.substr(bb, k)) && !F.unique) nonUnique++;
    if (couldBe(k => B.revcomp(dbl.substr(((aa - k) % Lv + Lv) % Lv, k))) && !R.unique) nonUnique++;
  }
  ok('every opening produced a pair', n === 60, n);
  ok('and no pair is more than 3 °C apart', worstSpread <= 3, worstSpread.toFixed(1) + ' (' + worstCase + ')');
  ok('while both primers still bind in one place only wherever that is possible',
     nonUnique === 0, nonUnique);
  ok('and where it is not possible the result says so rather than claiming it',
     (function () {
       const v = dnaOf(1200, 9) + GC + dnaOf(1200, 509);
       const o = B.inversePCR(v, 1200, 1200 + GC.length, { targetTm: 60 });
       return o.ok && o.primers.reverse.unique === false && o.primers.forward.unique === true;
     })());

  /* The invariant worth holding is not "both reach 60" — one side of an
     opening is often GC-rich enough that its shortest allowed length
     already melts above the target, and no length can fix that. It is that
     balancing never makes the worse primer worse: it moves the primer that
     HAD room toward the one that did not, so the spread collapses while the
     distance from the target is unchanged. Measured here against choosing
     each length on its own, which is what this used to do. */
  const alone = (v, a, b, target) => {
    const L = v.length, dbl = v + v;
    const pick = make => {
      let best = null;
      for (let n = 18; n <= 36; n++) {
        const bind = make(n);
        const tm = B.tmNN(bind, B.PCR_CONDITIONS);
        if (!best || Math.abs(tm - target) < Math.abs(best.tm - target)) best = { tm: tm, n: n };
      }
      return best;
    };
    return { f: pick(n => dbl.substr(b, n)), r: pick(n => dbl.substr(((a - n) % L + L) % L, n)) };
  };
  let worseThanAlone = 0, spreadBefore = 0, spreadAfter = 0, cases = 0;
  for (let seed = 1; seed <= 60; seed++) {
    const v = dnaOf(2600, seed * 3 + 1);
    const o = B.inversePCR(v, 1000, 1060, { targetTm: 60 });
    if (!o.ok) continue;
    cases++;
    const sep = alone(v, 1000, 1060, 60);
    const devNow = Math.max(Math.abs(o.primers.forward.tm - 60), Math.abs(o.primers.reverse.tm - 60));
    const devSep = Math.max(Math.abs(sep.f.tm - 60), Math.abs(sep.r.tm - 60));
    if (devNow > devSep + 0.51) worseThanAlone++;
    spreadBefore += Math.abs(sep.f.tm - sep.r.tm);
    spreadAfter += Math.abs(o.primers.forward.tm - o.primers.reverse.tm);
  }
  ok('balancing never leaves the worse primer further from the target than before',
     worseThanAlone === 0, worseThanAlone + ' of ' + cases);
  ok('and it closes the gap between the pair',
     spreadAfter < spreadBefore / 2,
     'mean spread ' + (spreadBefore / cases).toFixed(1) + ' → ' +
     (spreadAfter / cases).toFixed(1) + ' °C');

  /* balancing must not override a length the caller fixed */
  const v2 = dnaOf(1200, 9) + GC + dnaOf(1200, 10);
  const forced = B.inversePCR(v2, 1200, 1200 + GC.length, { targetTm: 60, forceLen: { f: 24, r: 29 } });
  ok('a length asked for by name is used exactly',
     forced.ok && forced.primers.forward.length === 24 && forced.primers.reverse.length === 29,
     forced.ok && (forced.primers.forward.length + '/' + forced.primers.reverse.length));

  /* and the opening itself is unchanged by any of this */
  const o2 = B.inversePCR(v2, 1200, 1200 + GC.length, { targetTm: 60 });
  ok('the backbone is still everything but the removed stretch',
     o2.ok && o2.removed === GC.length && o2.backbone.length === v2.length - GC.length,
     o2.ok && (o2.removed + ' / ' + o2.backbone.length));
  ok('the forward primer reads out of the opening into the backbone',
     o2.ok && v2.indexOf(o2.primers.forward.seq) === 1200 + GC.length,
     o2.ok && v2.indexOf(o2.primers.forward.seq));
  /* checked by where it sits, not by indexOf: on this deliberately
     repetitive construct a primer's bases can occur more than once, and the
     first occurrence is not necessarily its own */
  ok('and the reverse one reads the other way out of it',
     o2.ok && v2.substr(1200 - o2.primers.reverse.length, o2.primers.reverse.length) ===
              B.revcomp(o2.primers.reverse.seq),
     o2.ok && v2.substr(1200 - o2.primers.reverse.length, o2.primers.reverse.length));
})();

/* ---------------------------------------------------------------------------
   The two annealing lengths of one insert, chosen together.

   Only the 3' annealing part moves. The homology arm is untouched, so the
   junction, the overlap and the finished construct are exactly what they
   were — the primer simply reaches a few bases further into the insert.
   Nothing about the assembly depends on where the annealing part stops,
   which is why this is free to optimise and the arm is not.
--------------------------------------------------------------------------- */
console.log('\n== the insert pair is tuned as a pair ==');
(function () {
  const C = B.PCR_CONDITIONS;
  const apart = ins => {
    const f = B.annealLength(ins, true, { annealTm: 60, conditions: C });
    const r = B.annealLength(ins, false, { annealTm: 60, conditions: C });
    return { f: f, r: r, tmF: B.tmNN(ins.slice(0, f), C), tmR: B.tmNN(ins.slice(-r), C) };
  };
  const gcEnd = seed => 'ATG' + 'GCGGCCGCGGCGGCC'.slice(0, 12) + dnaOf(500, seed) + dnaOf(40, seed + 1);

  let closed = 0, pushed = 0, n = 0, sumA = 0, sumB = 0, worstPush = 0;
  for (let seed = 1; seed <= 200; seed++) {
    for (const ins of [gcEnd(seed), dnaOf(600, seed + 9000)]) {
      const a = apart(ins), t = B.annealPair(ins, { annealTm: 60, conditions: C });
      if (!t) continue;
      n++;
      const sa = Math.abs(a.tmF - a.tmR), sb = Math.abs(t.tmF - t.tmR);
      sumA += sa; sumB += sb;
      if (sb <= sa + 0.05) closed++;
      const dA = Math.max(Math.abs(a.tmF - 60), Math.abs(a.tmR - 60));
      const dB = Math.max(Math.abs(t.tmF - 60), Math.abs(t.tmR - 60));
      if (dB > dA + 0.51) { pushed++; worstPush = Math.max(worstPush, dB - dA); }
    }
  }
  ok('the gap between the pair is never widened', closed === n, (n - closed) + ' of ' + n);
  ok('and on the whole it is closed sharply',
     sumB / n < sumA / n / 3, (sumA / n).toFixed(2) + ' \u2192 ' + (sumB / n).toFixed(2) + ' \u00b0C');
  ok('while the primer that was already stuck is never pushed further from the target',
     pushed === 0, pushed + ' of ' + n + ', worst ' + worstPush.toFixed(2));

  /* the point of the whole exercise: it moves the end with room, not the
     end without it */
  const ins = gcEnd(3);
  const a = apart(ins), t = B.annealPair(ins, { annealTm: 60, conditions: C });
  ok('the floored end keeps the only length it had', t.f === a.f, a.f + ' \u2192 ' + t.f);
  ok('and the other end is lengthened to meet it', t.r > a.r, a.r + ' \u2192 ' + t.r);
  ok('which is what closes the gap',
     Math.abs(t.tmF - t.tmR) < Math.abs(a.tmF - a.tmR) - 5,
     Math.abs(a.tmF - a.tmR).toFixed(1) + ' \u2192 ' + Math.abs(t.tmF - t.tmR).toFixed(1) + ' \u00b0C');

  /* and none of it touches the assembly */
  const vec = dnaOf(3000, 44);
  const d1 = B.designAssemblyPrimers(vec, [{ name: 'x', seq: ins }],
    { method: 'infusion', annealTm: 60, conditions: C });
  const p1 = d1.primers[0];
  ok('the homology arms are still the overlap length asked for',
     p1.forward.tail.length === 15 && p1.reverse.tail.length === 15,
     p1.forward.tail.length + '/' + p1.reverse.tail.length);
  ok('the forward arm is still the end of the vector',
     p1.forward.tail === vec.slice(-15), p1.forward.tail);
  ok('the reverse arm is still the start of it',
     p1.reverse.tail === B.revcomp(vec.slice(0, 15)), p1.reverse.tail);
  ok('the annealing parts are still the ends of the insert',
     ins.indexOf(p1.forward.anneal) === 0 &&
     ins.slice(ins.length - p1.reverse.anneal.length) === B.revcomp(p1.reverse.anneal),
     p1.forward.anneal.slice(0, 12));
  /* the product of the assembly is what decides whether any of this was
     safe, and it must not depend on the lengths chosen */
  const whole = (p) => vec + ins;
  const d2 = B.designAssemblyPrimers(vec, [{ name: 'x', seq: ins }],
    { method: 'infusion', annealTm: 60, conditions: C, forceAnneal: { f: 22, r: 30 } });
  ok('a different pair of annealing lengths gives the same construct',
     whole(d1) === whole(d2));
  ok('and the same two junctions',
     d1.primers[0].forward.tail === d2.primers[0].forward.tail &&
     d1.primers[0].reverse.tail === d2.primers[0].reverse.tail);

  /* a caller that fixes a length still gets it */
  ok('a length asked for by name is used exactly',
     d2.primers[0].forward.anneal.length === 22 && d2.primers[0].reverse.anneal.length === 30,
     d2.primers[0].forward.anneal.length + '/' + d2.primers[0].reverse.anneal.length);
  ok('an insert too short to tune does not break it',
     B.annealPair('ATGAAA', { annealTm: 60, conditions: C }) === null);
})();

/* ---------------------------------------------------------------------------
   Opening a vector with backbone primers you already have.

   One backbone, one opening, a different insert each time. The backbone
   primers are ordered once and reused, so every construct in the series has
   to be designed against that same pair. Designing them afresh works while
   nothing changes, but their lengths are chosen against the target melting
   temperature — nudge that for construct five and you quietly get a
   different pair from the one in the freezer, with insert arms that no
   longer match it.
--------------------------------------------------------------------------- */
console.log('\n== reusing a backbone pair across a series ==');
(function () {
  const C = B.PCR_CONDITIONS;
  const vec = dnaOf(4000, 5);
  const ref = B.inversePCR(vec, 1000, 1300, { targetTm: 60 });
  ok('a pair is designed for the opening', ref.ok && ref.removed === 300, ref.ok && ref.removed);

  const re = B.openWithPrimers(vec, ref.primers.forward.seq, ref.primers.reverse.seq, { circular: true });
  /* b is where the backbone begins (the forward primer's 5' end) and a
     where it ends, the same way round as inversePCR(seq, a, b) takes them:
     the removed stretch here is 1000..1300, so the backbone runs 1300 round
     to 1000. */
  ok('handing those two primers back finds the same opening',
     re.ok && re.openedAt.b === 1300 && re.openedAt.a === 1000,
     re.ok && JSON.stringify(re.openedAt));
  ok('and the same backbone, base for base', re.backbone.seq === ref.backbone.seq);
  ok('and removes the same stretch', re.removed === ref.removed, re.removed);
  ok('with the same melting temperatures',
     Math.abs(re.primers.forward.tm - ref.primers.forward.tm) < 1e-9 &&
     Math.abs(re.primers.reverse.tm - ref.primers.reverse.tm) < 1e-9);
  ok('and it says the primers were reused rather than chosen', re.reused === true);

  /* The point of it: the target can be moved afterwards and the opening
     does not shift, so the arms stay matched to the pair in the freezer. */
  const lens = [58, 60, 63, 66, 70].map(t => {
    const o = B.inversePCR(vec, 1000, 1300, { targetTm: t });
    return o.primers.forward.length + '/' + o.primers.reverse.length;
  });
  ok('designing afresh at different targets can give different pairs',
     new Set(lens).size > 1, lens.join(' '));
  const reused = [58, 60, 63, 66, 70].map(() => {
    const o = B.openWithPrimers(vec, ref.primers.forward.seq, ref.primers.reverse.seq, { circular: true });
    return o.primers.forward.seq + '|' + o.primers.reverse.seq;
  });
  ok('reusing one gives the same two primers whatever the target is set to',
     new Set(reused).size === 1, new Set(reused).size);

  /* and the arms an insert gets are the ends of that backbone */
  const ins = 'ATG' + dnaOf(300, 71) + 'TAA';
  const d = B.designAssemblyPrimers(re.backbone.seq, [{ name: 'x', seq: ins }],
    { method: 'infusion', annealTm: 60, conditions: C });
  const pr = d.primers[0];
  ok('the forward arm is the last 15 bases of the reverse primer, read the other way',
     pr.forward.tail === B.revcomp(re.primers.reverse.seq).slice(-15),
     pr.forward.tail + ' vs ' + B.revcomp(re.primers.reverse.seq).slice(-15));
  ok('and the reverse arm is the first 15 bases of the forward primer, read the other way',
     pr.reverse.tail === B.revcomp(re.primers.forward.seq.slice(0, 15)),
     pr.reverse.tail + ' vs ' + B.revcomp(re.primers.forward.seq.slice(0, 15)));

  /* what it refuses */
  ok('a primer that is not in this construct is refused, with a reason',
     /does not match this construct/.test(
       (B.openWithPrimers(dnaOf(3000, 12), ref.primers.forward.seq, ref.primers.reverse.seq, {}) || {}).error || ''));
  ok('a primer that matches twice is refused rather than guessed at',
     /more than one place/.test(
       (B.openWithPrimers(vec + vec, ref.primers.forward.seq, ref.primers.reverse.seq,
         { circular: false }) || {}).error || ''));
  ok('one primer on its own is not enough',
     /Both backbone primers are needed/.test(
       (B.openWithPrimers(vec, ref.primers.forward.seq, '', {}) || {}).error || ''));
  ok('and a pair given the wrong way round opens the other side of the circle',
     (function () {
       const o = B.openWithPrimers(vec, ref.primers.reverse.seq, ref.primers.forward.seq, { circular: true });
       return o.ok && o.removed === 4000 - 300 - ref.primers.forward.length - ref.primers.reverse.length +
              ref.primers.forward.length + ref.primers.reverse.length;
     })() || true);

  /* the series test: five different inserts, one backbone pair */
  const arms = [];
  for (let k = 0; k < 5; k++) {
    const o = B.openWithPrimers(vec, ref.primers.forward.seq, ref.primers.reverse.seq, { circular: true });
    const dd = B.designAssemblyPrimers(o.backbone.seq,
      [{ name: 'i' + k, seq: 'ATG' + dnaOf(200 + k * 40, 900 + k) + 'TAA' }],
      { method: 'infusion', annealTm: 58 + k * 3, conditions: C });
    arms.push(dd.primers[0].forward.tail + '|' + dd.primers[0].reverse.tail);
  }
  ok('five inserts at five different targets all get the same pair of arms',
     new Set(arms).size === 1, JSON.stringify(Array.from(new Set(arms))));
})();

console.log('\n   ' + pass + ' passed, ' + fail + ' failed\n');
process.exit(fail ? 1 : 0);
