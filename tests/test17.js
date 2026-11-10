/* ===========================================================================
   What a read says about the protein, and BLAST results pasted back.

   The hard part is the map from a genomic coordinate to a codon, because a
   coding feature can be reverse-stranded, spliced, wrapped across the
   origin, or all three. The map is checked here against featureSequence(),
   which the validation suite already compares base for base against
   Biopython over the whole GenBank corpus — so if the map is wrong, it is
   wrong against an independently verified answer, not against itself.

   The consequence calls are checked against codons worked out by hand.
   ========================================================================= */
const B = require('../src/engine7.js');
let pass = 0, fail = 0;
function ok(name, cond, extra) {
  if (cond) { pass++; }
  else { fail++; console.log('  FAIL  ' + name + (extra !== undefined ? '   -> ' + extra : '')); }
}
function section(s) { console.log('\n== ' + s + ' =='); }

/* a deterministic pseudo-random sequence, so a failure is reproducible */
let seed = 20261106;
function rnd() { seed = (seed * 1103515245 + 12345) & 0x7fffffff; return seed / 0x7fffffff; }
function randSeq(n) { let s = ''; for (let i = 0; i < n; i++) s += 'ACGT'[Math.floor(rnd() * 4)]; return s; }

/* ------------------------------------------------------- the coordinate map */
section('the genomic position to codon map');
{
  const rec = { seq: randSeq(600), features: [] };
  const L = rec.seq.length;
  const cases = [
    ['plain forward', { type: 'CDS', label: 'f1', start: 30, end: 120, strand: 1 }],
    ['plain reverse', { type: 'CDS', label: 'f2', start: 30, end: 120, strand: -1 }],
    ['spliced forward', { type: 'CDS', label: 'f3', start: 30, end: 180, strand: 1,
      parts: [{ start: 30, end: 75 }, { start: 120, end: 180 }] }],
    ['spliced reverse', { type: 'CDS', label: 'f4', start: 30, end: 180, strand: -1,
      parts: [{ start: 30, end: 75 }, { start: 120, end: 180 }] }],
    ['wrapped forward', { type: 'CDS', label: 'f5', start: 570, end: 30, strand: 1,
      parts: [{ start: 570, end: 600 }, { start: 0, end: 30 }] }],
    ['wrapped reverse', { type: 'CDS', label: 'f6', start: 570, end: 30, strand: -1,
      parts: [{ start: 570, end: 600 }, { start: 0, end: 30 }] }]
  ];
  cases.forEach(([name, f]) => {
    const idx = B.cdsIndex(rec, f);
    let built = '';
    for (let i = 0; i < idx.length; i++) built += idx.base(i);
    ok(name + ': the map rebuilds exactly what featureSequence gives',
      built === B.featureSequence(rec.seq, f), name + ' got ' + built.slice(0, 24) + '… want ' +
      B.featureSequence(rec.seq, f).slice(0, 24) + '…');
    ok(name + ': its length is the feature’s length',
      idx.length === B.featureLength(f, L), idx.length + ' vs ' + B.featureLength(f, L));
    ok(name + ': every genomic position maps back to the base it holds', (() => {
      for (let i = 0; i < idx.length; i++) {
        const back = idx.codonOf(idx.pos[i]);
        if (!back || back.cdsPos !== i) return false;
      }
      return true;
    })());
    ok(name + ': a position outside the feature maps to nothing', (() => {
      for (let g = 0; g < L; g++) {
        if (idx.at[g] !== undefined) continue;
        if (idx.codonOf(g) !== null) return false;
      }
      return true;
    })());
  });
  /* the whole translated protein agrees with translating featureSequence */
  cases.forEach(([name, f]) => {
    const idx = B.cdsIndex(rec, f);
    let fromMap = '';
    for (let k = 0; k * 3 + 2 < idx.length; k++) fromMap += B.translate(B.cdsCodon(idx, k));
    const direct = B.translate(B.featureSequence(rec.seq, f));
    ok(name + ': codon by codon it translates to the same protein',
      fromMap === direct.slice(0, fromMap.length), name);
  });
}

/* ------------------------------------------------------ consequence calling */
section('what a difference does to the protein');
{
  /* ATG AAA TTT TGG TAA — Met Lys Phe Trp stop, at position 10 */
  const seq = 'CCCCCCCCCC' + 'ATGAAATTTTGGTAA' + 'CCCCCCCCCC';
  const rec = { seq: seq, features: [] };
  const f = { type: 'CDS', label: 'tiny', start: 10, end: 25, strand: 1 };
  const eff = (refPos, read, type) => B.diffEffect(rec, f, { type: type || 'mismatch', refPos: refPos, ref: seq[refPos], read: read });

  ok('the reference protein is MKFW*', B.translate(B.featureSequence(seq, f)) === 'MKFW*',
    B.translate(B.featureSequence(seq, f)));
  /* AAA -> AAG is still Lys */
  let e = eff(15, 'G');
  ok('a third-base change that keeps the residue is silent', e && e.kind === 'silent', e && e.text);
  ok('and it names the right residue', e && e.residue === 2, e && e.residue);
  /* AAA -> GAA is Glu */
  e = eff(13, 'G');
  ok('a first-base change that swaps the residue is missense', e && e.kind === 'missense', e && e.text);
  ok('and it reads K2E', e && e.from === 'K' && e.to === 'E' && e.residue === 2, e && e.text);
  /* TGG -> TGA is a stop */
  e = eff(21, 'A');
  ok('a change that makes a stop codon is nonsense', e && e.kind === 'nonsense', e && e.text);
  ok('and it is at residue 4', e && e.residue === 4, e && e.residue);
  /* TAA -> TAC loses the stop */
  e = eff(23, 'C');
  ok('a change that loses the stop is called as such', e && e.kind === 'readthrough', e && e.text);
  e = eff(14, 'G', 'deletion');
  ok('a single-base deletion in a CDS is a frameshift', e && e.kind === 'frameshift', e && e.text);
  ok('a difference outside the CDS has no consequence for it', eff(3, 'G') === null);

  /* the same gene on the other strand must give the same residue numbers */
  const rc = B.revcomp(seq);
  const rrec = { seq: rc, features: [] };
  const rf = { type: 'CDS', label: 'tiny', start: rc.length - 25, end: rc.length - 10, strand: -1 };
  ok('the reverse-strand copy translates identically',
    B.translate(B.featureSequence(rc, rf)) === 'MKFW*', B.translate(B.featureSequence(rc, rf)));
  /* position of codon 2 base 1 on the reverse strand */
  const gpos = rc.length - 1 - 13;
  const re = B.diffEffect(rrec, rf, { type: 'mismatch', refPos: gpos, ref: rc[gpos], read: B.revcomp('G') });
  ok('a reverse-strand change is complemented before the codon is read',
    re && re.kind === 'missense' && re.from === 'K' && re.to === 'E' && re.residue === 2,
    re && re.text);
}

/* a CDS that wraps the origin, which slicing start..end cannot even read */
section('a coding sequence across the origin');
{
  const tail = 'ATGAAA', head = 'TTTTGGTAA';
  const seq = tail + 'GGGGGGGGGGGGGGGGGGGG' + head;   /* CDS = last 6 + first 9 */
  const rec = { seq: seq, features: [] };
  const f = { type: 'CDS', label: 'wrap', start: seq.length - 9, end: 6, strand: 1,
    parts: [{ start: seq.length - 9, end: seq.length }, { start: 0, end: 6 }] };
  ok('naive slicing would read nothing', seq.slice(f.start, f.end) === '');
  /* the coding sequence is the last nine bases then the first six:
     TTT TGG TAA ATG AAA, so the stop sits in the middle of it */
  ok('the map reads it correctly', B.translate(B.featureSequence(seq, f)) === 'FW*MK',
    B.translate(B.featureSequence(seq, f)));
  const idx = B.cdsIndex(rec, f);
  ok('and the map agrees', (() => { let s = ''; for (let i = 0; i < idx.length; i++) s += idx.base(i); return s; })()
    === B.featureSequence(seq, f));
  /* base 1 of the plasmid is coding position 10, so residue 4 — the ATG
     that a reader would wrongly call residue 1 if they sliced naively */
  const e = B.diffEffect(rec, f, { type: 'mismatch', refPos: 0, ref: 'A', read: 'G' });
  ok('a change at base 1 of the plasmid is placed in the right codon',
    e && e.residue === 4 && e.text === 'missense \u00b7 M4V', e && JSON.stringify(e));
}

/* ----------------------------------------------------------- read overlap */
section('what one read covers');
{
  const seq = randSeq(1000);
  const rec = { seq: seq, features: [
    { type: 'CDS', label: 'geneA', start: 100, end: 400, strand: 1 },
    { type: 'promoter', label: 'pro', start: 140, end: 180, strand: 1 },
    { type: 'CDS', label: 'geneB', start: 600, end: 900, strand: -1 }
  ] };
  const aln = { refStart: 150, refEnd: 500, differences: [], alignment: null };
  const ov = B.readOverlap(rec, aln);
  ok('only the features it touches come back', ov.length === 2, ov.map(o => o.label).join(','));
  const a = ov.filter(o => o.label === 'geneA')[0];
  ok('the covered bases are counted, not the whole feature', a.bases === 250, a.bases);
  ok('a feature the read only clips is counted for the part it clipped',
    ov.filter(o => o.label === 'pro')[0].bases === 30, ov.filter(o => o.label === 'pro')[0].bases);
  ok('and reported as a fraction of it', Math.abs(a.fraction - 250 / 300) < 1e-9, a.fraction);
  ok('a coding feature reports which residues the read reached',
    a.residues && a.residues.from === 17 && a.residues.to === 100, JSON.stringify(a.residues));
  ok('a non-coding feature reports none', ov.filter(o => o.label === 'pro')[0].residues === null);
  ok('the list is ordered by how much was covered', ov[0].bases >= ov[1].bases);
  ok('geneB, which the read never reached, is absent', !ov.some(o => o.label === 'geneB'));
}

/* ---------------------------------------------------------- read protein */
section('the protein the read itself supports');
{
  const cds = 'ATGAAATTTTGGCCCGGGAAATTTTGA';          /* M K F W P G K F stop */
  const seq = 'TTTTTTTTTT' + cds + 'AAAAAAAAAA';
  const rec = { seq: seq, features: [] };
  const f = { type: 'CDS', label: 'p', start: 10, end: 10 + cds.length, strand: 1 };
  /* a read covering only codons 2-5, with one mismatch making K2E */
  const refStart = 13, refEnd = 25;
  const alnA = [], alnB = [];
  for (let g = refStart; g < refEnd; g++) { alnA.push(seq[g]); alnB.push(g === 13 ? 'G' : seq[g]); }
  const aln = { refStart: refStart, refEnd: refEnd, differences: [{ type: 'mismatch', refPos: 13, ref: 'A', read: 'G' }],
                alignment: { a: alnA.join(''), b: alnB.join('') } };
  const p = B.readProtein(rec, aln, f);
  ok('the reference protein is the whole gene', p.reference === 'MKFWPGKF*', p.reference);
  ok('positions the read never reached are not filled in from the reference',
    p.observed[0] === '.' && p.observed[p.observed.length - 1] === '.', p.observed);
  ok('the covered residues are called from the read', p.observed.indexOf('E') !== -1, p.observed);
  ok('the change is reported with its residue number',
    p.changes.length === 1 && p.changes[0].residue === 2 && p.changes[0].from === 'K' && p.changes[0].to === 'E',
    JSON.stringify(p.changes));
  ok('coverage is counted in residues actually read', p.covered === 4, p.covered);
  ok('and as a percentage of the protein', Math.abs(p.percent - 4 / 9 * 100) < 1e-9, p.percent);
  ok('it is not called identical when it is not', p.identical === false);
  /* a clean read over the same stretch */
  const clean = { refStart: refStart, refEnd: refEnd, differences: [],
                  alignment: { a: alnA.join(''), b: alnA.join('') } };
  ok('a clean read over part of a gene reports no changes', B.readProtein(rec, clean, f).changes.length === 0);
}

/* ------------------------------------------------------------ BLAST input */
section('BLAST results pasted back');
{
  const hitTable = [
    '# BLASTN 2.15.0+',
    '# Query: my_read',
    'my_read\tNM_001301717.2\t99.286\t560\t4\t0\t1\t560\t201\t760\t0.0\t1023',
    'my_read\tNM_001301717.2\t97.000\t100\t3\t0\t600\t699\t900\t999\t1e-40\t180',
    'my_read\tXM_011542432.1\t94.500\t400\t20\t2\t10\t409\t5\t404\t3e-170\t615'
  ].join('\n');
  const r1 = B.parseBlast(hitTable);
  ok('a hit table is recognised', r1.format === 'Hit table', r1.format);
  ok('its rows are grouped into one entry per subject', r1.hits.length === 2, r1.hits.length);
  ok('best bit score first', r1.hits[0].accession === 'NM_001301717.2', r1.hits[0].accession);
  ok('the two HSPs of that subject are both kept', r1.hits[0].hsps.length === 2, r1.hits[0].hsps.length);
  ok('and the best one sets the headline numbers',
    r1.hits[0].bits === 1023 && r1.hits[0].evalue === 0 && r1.hits[0].identity === 99.286,
    JSON.stringify([r1.hits[0].bits, r1.hits[0].evalue, r1.hits[0].identity]));

  const descCsv = 'Description,Scientific Name,Max Score,Total Score,Query Cover,E value,Per. Ident,Accession\n' +
    '"Homo sapiens kinase, mRNA",Homo sapiens,1023,1203,98%,0.0,99.29,NM_001301717.2\n' +
    '"Pan troglodytes kinase, transcript variant 1",Pan troglodytes,615,615,72%,3e-170,94.50,XM_011542432.1';
  const r2 = B.parseBlast(descCsv);
  ok('a Descriptions CSV is recognised', r2.format === 'Descriptions table', r2.format);
  ok('a description containing a comma is not split across columns',
    r2.hits[0].title === 'Homo sapiens kinase, mRNA', r2.hits[0].title);
  ok('the accession comes through', r2.hits[0].accession === 'NM_001301717.2', r2.hits[0].accession);
  ok('query cover is read as a number', r2.hits[0].cover === 98, r2.hits[0].cover);
  ok('e value is read in scientific notation', r2.hits[1].evalue === 3e-170, r2.hits[1].evalue);

  const json2 = JSON.stringify({ BlastOutput2: [{ report: { results: { search: { hits: [
    { num: 1, len: 2400, description: [{ accession: 'NM_001301717.2', title: 'Homo sapiens kinase' }],
      hsps: [{ bit_score: 1023.5, evalue: 0, identity: 556, align_len: 560, query_from: 1, query_to: 560,
               hit_from: 201, hit_to: 760, qseq: 'ACGT', hseq: 'ACGT' }] },
    { num: 2, len: 2100, description: [{ accession: 'XM_011542432.1', title: 'Pan troglodytes kinase' }],
      hsps: [{ bit_score: 615, evalue: 3e-170, identity: 378, align_len: 400 }] }
  ] } } } }] });
  const r3 = B.parseBlast(json2);
  ok('single-file JSON is recognised', r3.format === 'Single-file JSON', r3.format);
  ok('its hits come through in score order',
    r3.hits.length === 2 && r3.hits[0].accession === 'NM_001301717.2', JSON.stringify(r3.hits.map(h => h.accession)));
  ok('identity is worked out from the counts, not copied',
    Math.abs(r3.hits[0].identity - 556 / 560 * 100) < 1e-9, r3.hits[0].identity);
  ok('the subject length is kept', r3.hits[0].subjectLength === 2400, r3.hits[0].subjectLength);
  ok('the aligned sequences are kept when the format carries them', r3.hits[0].hitSeq === 'ACGT');

  const pairwise = [
    'Sequences producing significant alignments:',
    '',
    '>NM_001301717.2 Homo sapiens kinase, mRNA',
    'Length=2400',
    '',
    ' Score = 1023 bits (554),  Expect = 0.0',
    ' Identities = 556/560 (99%), Gaps = 0/560 (0%)',
    '',
    '>XM_011542432.1 Pan troglodytes kinase',
    'Length=2100',
    '',
    ' Score = 615 bits (333),  Expect = 3e-170',
    ' Identities = 378/400 (95%), Gaps = 2/400 (0%)'
  ].join('\n');
  const r4 = B.parseBlast(pairwise);
  ok('pairwise text off the screen is recognised', r4.format === 'Pairwise text', r4.format);
  ok('both subjects are found', r4.hits.length === 2, r4.hits.length);
  ok('the accession is split from the title',
    r4.hits[0].accession === 'NM_001301717.2' && /Homo sapiens kinase/.test(r4.hits[0].title),
    JSON.stringify([r4.hits[0].accession, r4.hits[0].title]));
  ok('the score and e value are read', r4.hits[0].bits === 1023 && r4.hits[0].evalue === 0,
    JSON.stringify([r4.hits[0].bits, r4.hits[0].evalue]));
  ok('identity is read as a percentage', r4.hits[0].identity === 99, r4.hits[0].identity);
  ok('subject length is read', r4.hits[0].subjectLength === 2400, r4.hits[0].subjectLength);

  let threw = '';
  try { B.parseBlast('the quick brown fox jumps'); } catch (e) { threw = e.code; }
  ok('something that is not a BLAST result is refused, not half-read', threw === 'unknown', threw);
  threw = '';
  try { B.parseBlast('{"nothing":1}'); } catch (e) { threw = e.code; }
  ok('JSON with no hits in it says so', threw === 'nohits', threw);
  ok('an empty paste gives nothing rather than throwing', B.parseBlast('   ').hits.length === 0);
}

/* ---------------------------------------------------------------------------
   The bases behind each translated frame.

   A peptide alone is half of what a verification gives you; the other half
   is the DNA it came from, which is what goes into a primer or an order.
   The invariant that makes it safe to hand over is that it translates back:
   a minus frame is read off the reverse complement, and bases taken from the
   wrong strand, or three off, look perfectly ordinary and mean something
   else entirely.
--------------------------------------------------------------------------- */
console.log('\n== the bases behind each frame ==');
(function () {
  const check = (seq, table, note) => {
    const frames = B.frameTranslate(seq, { table: table });
    let bad = [];
    frames.forEach(f => {
      if (B.translate(f.longestNt, { table: table }) !== f.longest) bad.push(f.label + ' longest');
      if (B.translate(f.cleanNt, { table: table }) !== f.clean) bad.push(f.label + ' clean');
      if (f.longestNt.length !== f.longest.length * 3) bad.push(f.label + ' longest length');
      if (f.cleanNt.length !== f.clean.length * 3) bad.push(f.label + ' clean length');
      if (f.fromMet) {
        if (B.translate(f.fromMetNt, { table: table }) !== f.fromMet) bad.push(f.label + ' fromMet');
        if (!/^ATG/.test(f.fromMetNt)) bad.push(f.label + ' fromMet start');
      } else if (f.fromMetNt) bad.push(f.label + ' fromMet should be empty');
    });
    ok(note, bad.length === 0, bad.join(', '));
    return frames;
  };

  /* a plain coding read, forward */
  const orf = 'ATGGGCAGCAGCCATCACCATCATCACCACGAAAACCTGTATTTTCAGGGCATGGAAGATTATACC' +
              'AAAATAGAGAAAATAGGAGAAGGAACTTACGGAGTTGTGTATAAAGGACGTCATAAAACAACAGGA' +
              'CAGGTAGTGGCAATGAAAAAGATTCGTCTAGAAAGTGAAGAGGAAGGGTAA';
  const fw = check('CCTTGACG' + orf + 'TTAGCCAT', 1, 'every frame of a forward read gives back its own bases');
  /* the same read written backwards: the frame that wins is a minus one,
     and its bases are the reverse complement, not what is in the file */
  const rv = check(B.revcomp('CCTTGACG' + orf + 'TTAGCCAT'), 1, 'and so does every frame of the same read reversed');
  const best = rv.slice().sort((a, b) => b.score - a.score)[0];
  ok('the called frame of a reversed read is a minus frame', best.dir === -1, best.label);
  ok('its bases are the reverse complement of the read, not a slice of it',
     B.translate(best.clean ? best.cleanNt : best.longestNt, {}) === (best.clean || best.longest));
  const fwBest = fw.slice().sort((a, b) => b.score - a.score)[0];
  ok('and read either way round, the same bases come back',
     (fwBest.clean ? fwBest.cleanNt : fwBest.longestNt) ===
     (best.clean ? best.cleanNt : best.longestNt),
     (fwBest.clean ? fwBest.cleanNt : fwBest.longestNt).slice(0, 24) + ' vs ' +
     (best.clean ? best.cleanNt : best.longestNt).slice(0, 24));
  ok('the bases stop where the peptide does, with no stop codon on the end',
     !/(TAA|TAG|TGA)$/.test(fwBest.clean ? fwBest.cleanNt : fwBest.longestNt) ||
     B.translate(fwBest.clean ? fwBest.cleanNt : fwBest.longestNt, {}).indexOf('*') < 0,
     (fwBest.clean ? fwBest.cleanNt : fwBest.longestNt).slice(-9));

  /* ambiguous bases, a bacterial table, and the ragged lengths that leave a
     trailing base or two outside the last whole codon */
  check('NNNN' + orf + 'NN', 11, 'a read with ambiguous ends behaves the same under table 11');
  for (let pad = 0; pad < 3; pad++)
    check('A'.repeat(pad) + orf + 'C'.repeat(pad), 1,
      'and so does one whose length leaves ' + pad + ' spare base(s) past the last codon');

  /* property check across random reads: no frame may ever hand back bases
     that translate to something other than the peptide beside them */
  let n = 0, bad = 0;
  let x = 20260107;
  const rnd = () => (x = (x * 1103515245 + 12345) & 0x7fffffff) / 0x7fffffff;
  for (let t = 0; t < 250; t++) {
    let s2 = '';
    const len = 120 + Math.floor(rnd() * 700);
    for (let i = 0; i < len; i++) s2 += (rnd() < 0.04 ? 'N' : 'ACGT'[Math.floor(rnd() * 4)]);
    B.frameTranslate(s2, { table: 1 }).forEach(f => {
      n++;
      if (B.translate(f.longestNt, {}) !== f.longest) bad++;
      if (B.translate(f.cleanNt, {}) !== f.clean) bad++;
      if (f.fromMet && B.translate(f.fromMetNt, {}) !== f.fromMet) bad++;
    });
  }
  ok('across ' + n + ' frames of random reads, the bases always translate back', bad === 0, bad);
})();

console.log('\n   ' + pass + ' passed, ' + fail + ' failed\n');
process.exit(fail ? 1 : 0);
