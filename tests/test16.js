/* ===========================================================================
   Alignment figures.

   A figure is read, not computed with, so a wrong one is not caught by the
   reader the way a wrong number is. These tests pin the things a wrong
   figure would get wrong quietly: which columns count as conserved under
   each published rule, what the consensus line says, where a structure
   element lands once the gaps of its row are opened up, and that the
   monochrome export really has no hue left in it.
   ========================================================================= */
const B = require('../src/engine6.js');
let pass = 0, fail = 0;
function ok(name, cond, extra) {
  if (cond) { pass++; }
  else { fail++; console.log('  FAIL  ' + name + (extra !== undefined ? '   -> ' + extra : '')); }
}
function section(s) { console.log('\n== ' + s + ' =='); }

/* ----------------------------------------------------------- BLOSUM62 */
section('BLOSUM62');
ok('the matrix is symmetric in all 576 cells', (() => {
  const al = 'ARNDCQEGHILKMFPSTWYVBZX*';
  for (const a of al) for (const b of al) if (B.blosum62(a, b) !== B.blosum62(b, a)) return false;
  return true;
})());
/* anchors a transcription error would move */
[['W', 'W', 11], ['C', 'C', 9], ['A', 'A', 4], ['A', 'W', -3], ['I', 'V', 3], ['L', 'M', 2],
 ['F', 'Y', 3], ['D', 'E', 2], ['K', 'R', 2], ['G', 'G', 6], ['P', 'P', 7], ['*', '*', 1],
 ['A', '*', -4]].forEach(([a, b, v]) =>
  ok('BLOSUM62 ' + a + b + ' = ' + v, B.blosum62(a, b) === v, B.blosum62(a, b)));
ok('an unknown letter is not scored as a match', B.blosum62('J', 'A') === null);

/* ------------------------------------------------------ column scores */
section('column scores, by the published rules');
const col = (s) => s.split('');
ok('an identical column scores 1 in every mode',
  ['strict', 'multalin', 'physchem', 'blosum62', 'purine']
    .every(m => B.alnColumnScore(col('WWWW'), m) === 1));
ok('a column of one residue and three gaps scores 0 (nothing to conserve)',
  B.alnColumnScore(col('W---'), 'strict') === 0, B.alnColumnScore(col('W---'), 'strict'));
ok('three of four identical is 0.75, not 1',
  B.alnColumnScore(col('WWWY'), 'strict') === 0.75, B.alnColumnScore(col('WWWY'), 'strict'));
ok('a gap dilutes the column rather than being ignored',
  B.alnColumnScore(col('WWW-'), 'strict') === 0.75, B.alnColumnScore(col('WWW-'), 'strict'));

ok('MultAlin: I and V are one group', B.alnColumnScore(col('IVIV'), 'multalin') === 1);
ok('MultAlin: L and M are one group', B.alnColumnScore(col('LMLM'), 'multalin') === 1);
ok('MultAlin: F and Y are one group', B.alnColumnScore(col('FYFY'), 'multalin') === 1);
ok('MultAlin: N D Q E B Z are one group', B.alnColumnScore(col('NDQE'), 'multalin') === 1);
ok('MultAlin: I and L are NOT one group', B.alnColumnScore(col('ILIL'), 'multalin') === 0.5,
  B.alnColumnScore(col('ILIL'), 'multalin'));

ok('physico-chemical: H K R', B.alnColumnScore(col('HKRH'), 'physchem') === 1);
ok('physico-chemical: A V L I M', B.alnColumnScore(col('AVLI'), 'physchem') === 1);
ok('physico-chemical: F Y W', B.alnColumnScore(col('FYWF'), 'physchem') === 1);
ok('physico-chemical: S T N Q', B.alnColumnScore(col('STNQ'), 'physchem') === 1);
ok('physico-chemical: P G C', B.alnColumnScore(col('PGCP'), 'physchem') === 1);
ok('physico-chemical: D and K are not equivalent',
  B.alnColumnScore(col('DKDK'), 'physchem') === 0.5, B.alnColumnScore(col('DKDK'), 'physchem'));

ok('BLOSUM62 mode: every pair positive gives 1', B.alnColumnScore(col('IVIV'), 'blosum62') === 1);
/* The score is the fraction of residue PAIRS in the column that score above
   zero in BLOSUM62 — the same thing BLAST marks with a '+'. A column of
   W P W P therefore is not zero: W/W (+11) and P/P (+7) are two of its six
   pairs, so it scores 2/6. Zero needs a column with no positive pair at all. */
ok('BLOSUM62 mode: unrelated residues give 0',
  B.alnColumnScore(col('WPGC'), 'blosum62') === 0, B.alnColumnScore(col('WPGC'), 'blosum62'));
ok('BLOSUM62 mode: a residue against itself is always positive, so a repeated ' +
  'pair of unrelated residues is not zero',
  Math.abs(B.alnColumnScore(col('WPWP'), 'blosum62') - 2 / 6) < 1e-12,
  B.alnColumnScore(col('WPWP'), 'blosum62'));
ok('BLOSUM62 mode: a half-conserved column lands between',
  (() => { const v = B.alnColumnScore(col('IVWP'), 'blosum62'); return v > 0 && v < 1; })(),
  B.alnColumnScore(col('IVWP'), 'blosum62'));

ok('purine mode: a transition A/G counts as conserved', B.alnColumnScore(col('AGAG'), 'purine') === 1);
ok('purine mode: a transition C/T counts as conserved', B.alnColumnScore(col('CTCT'), 'purine') === 1);
ok('purine mode: a transversion A/C does not',
  B.alnColumnScore(col('ACAC'), 'purine') === 0.5, B.alnColumnScore(col('ACAC'), 'purine'));
ok('purine mode: U sits with the pyrimidines', B.alnColumnScore(col('CUCU'), 'purine') === 1);

/* -------------------------------------------------------- consensus */
section('the consensus line');
const consOf = rows => {
  const L = rows[0].length;
  return B.alnConsensus(Array.from({ length: L }, (_, i) => rows.map(r => r[i])));
};
ok('strict identity is upper case', consOf(['AW', 'AW']) === 'AW', consOf(['AW', 'AW']));
ok('I/V is !', consOf(['I', 'V']) === '!', consOf(['I', 'V']));
ok('L/M is $', consOf(['L', 'M']) === '$', consOf(['L', 'M']));
ok('F/Y is %', consOf(['F', 'Y']) === '%', consOf(['F', 'Y']));
ok('N/D/Q/E is #', consOf(['N', 'D', 'Q', 'E']) === '#', consOf(['N', 'D', 'Q', 'E']));
ok('a majority residue is lower case',
  consOf(['A', 'A', 'A', 'W']) === 'a', consOf(['A', 'A', 'A', 'W']));
ok('no majority gives a blank', consOf(['A', 'W', 'P', 'G']) === ' ', JSON.stringify(consOf(['A', 'W', 'P', 'G'])));
ok('a column of gaps gives a blank', consOf(['-', '-']) === ' ', JSON.stringify(consOf(['-', '-'])));

/* -------------------------------------------------- secondary structure */
section('the structure track');
const ss1 = B.parseSecondaryStructure('HHHEEE', 6, 'ACDEFG');
ok('DSSP letters become runs', ss1.runs.length === 2 && ss1.runs[0].kind === 'helix' && ss1.runs[1].kind === 'strand',
  JSON.stringify(ss1.runs));
ok('a run carries its half-open span', ss1.runs[1].start === 3 && ss1.runs[1].end === 6, JSON.stringify(ss1.runs[1]));
/* the string is given per residue; the gaps of its row have to be opened up */
const ss2 = B.parseSecondaryStructure('HHHEEE', 8, 'AC--DEFG');
ok('a per-residue string is spread across the gaps of its row',
  ss2.kinds[0] === 'helix' && ss2.kinds[1] === 'helix' &&
  ss2.kinds[2] === null && ss2.kinds[3] === null && ss2.kinds[4] === 'helix' && ss2.kinds[5] === 'strand',
  JSON.stringify(ss2.kinds));
ok('3-10 and pi helices are told apart from an alpha helix', (() => {
  const s = B.parseSecondaryStructure('HHGGII', 6, 'ACDEFG');
  return s.runs.map(r => r.kind).join(',') === 'helix,helix310,helixPi';
})(), JSON.stringify(B.parseSecondaryStructure('HHGGII', 6, 'ACDEFG').runs.map(r => r.kind)));
const db = B.parseSecondaryStructure('((((....))))', 12, 'ACGUACGUACGU');
ok('dot-bracket is recognised and read as stems', db.dotBracket === true &&
  db.runs.length === 2 && db.runs[0].kind === 'stem5' && db.runs[1].kind === 'stem3',
  JSON.stringify(db.runs));
ok('no structure given means no track', B.parseSecondaryStructure('   ', 3, 'ACG') === null);

/* ------------------------------------------------------------ analysis */
section('the analysis a figure is drawn from');
const prot = [
  { name: 'a', seq: 'MKVLAAGIVG--LNL' },
  { name: 'b', seq: 'MKVLAAGIVG--LNL' },
  { name: 'c', seq: 'MKILSAGVVGKKMNI' }
];
const A = B.alnAnalyse(prot, { mode: 'physchem', threshold: 0.7 });
ok('rows are padded to one width', A.rows.every(r => r.seq.length === A.width));
/* Counted by hand down the three rows: columns 0 1 3 5 6 8 9 13 are identical
   in all three (M K _ L _ A G _ V G _ _ _ N _), which is 8 of 15. Columns 2 4
   7 12 14 differ and 10 11 are gaps. The 9 that used to stand here was a
   miscount of this very alignment. */
ok('strictly conserved columns are found',
  A.strict.filter(Boolean).length === 8, A.strict.filter(Boolean).length);
ok('a column with a gap is never strictly conserved',
  !A.strict[10] && !A.strict[11], JSON.stringify([A.strict[10], A.strict[11]]));
ok('frames are contiguous runs of conserved columns, not single cells',
  A.frames.length >= 1 && A.frames.every(f => f.end > f.start), JSON.stringify(A.frames));
ok('identity is reported as a percentage of columns',
  Math.abs(A.identity - 8 / 15 * 100) < 1e-9, A.identity);
ok('a lower threshold cannot shade fewer columns', (() => {
  const lo = B.alnAnalyse(prot, { mode: 'physchem', threshold: 0.3 });
  const hi = B.alnAnalyse(prot, { mode: 'physchem', threshold: 0.95 });
  return lo.over.filter(Boolean).length >= hi.over.filter(Boolean).length;
})());

/* ------------------------------------------------------------ molecule */
section('protein or nucleotide');
ok('a peptide is read as protein', B.looksProteinSeq('MKVLAAGIVGLNLGGKLYEQ'));
ok('a gene is not', !B.looksProteinSeq('ATGGCTAGCAAAGGAGAAGAACTTTTCACT'));
ok('a short ACGT-only peptide is still read as nucleotide (the common case)',
  !B.looksProteinSeq('ACGTACGTACGT'));
ok('RNA is not read as protein', !B.looksProteinSeq('AUGGCUAGCAAAGGAGAAGAA'));
ok('the mode lists keep the groupings apart',
  B.ALN_MODES.protein.indexOf('purine') === -1 && B.ALN_MODES.nucleotide.indexOf('physchem') === -1,
  JSON.stringify(B.ALN_MODES));

/* -------------------------------------------------------------- the SVG */
section('the drawing');
const svg = B.alnFigureSVG(prot, { mode: 'physchem', ss: 'HHHHHEEEEEEEE', title: 'test' });
ok('it is a standalone SVG document', /^<\?xml/.test(svg) && /<svg[^>]+xmlns=/.test(svg));
ok('it carries its own size', /width="\d+"/.test(svg) && /height="\d+"/.test(svg));
ok('a strictly conserved column is a red box', svg.indexOf('#e02020') !== -1);
ok('conserved runs are framed in blue', svg.indexOf('#2b4fd8') !== -1);
ok('every sequence name is in it', prot.every(p => svg.indexOf('>' + p.name + '<') !== -1));
ok('the title is in it', svg.indexOf('>test<') !== -1);
ok('the letters are real text, not paths', (svg.match(/<text/g) || []).length > 30);
const monoSvg = B.alnFigureSVG(prot, { mode: 'physchem', mono: true });
ok('monochrome leaves no colour behind', (() => {
  const cols = monoSvg.match(/#[0-9a-f]{6}/gi) || [];
  return cols.every(c => {
    const v = [1, 3, 5].map(i => parseInt(c.substr(i, 2), 16));
    return Math.max.apply(null, v) - Math.min.apply(null, v) <= 8;
  });
})(), (monoSvg.match(/#[0-9a-f]{6}/gi) || []).filter(c => {
  const v = [1, 3, 5].map(i => parseInt(c.substr(i, 2), 16));
  return Math.max.apply(null, v) - Math.min.apply(null, v) > 8;
}).join(' '));
ok('monochrome is still the same alignment', prot.every(p => monoSvg.indexOf('>' + p.name + '<') !== -1));

let threw = '';
try { B.alnFigureSVG([{ name: 'only', seq: 'ACGT' }], {}); } catch (e) { threw = e.message; }
ok('one sequence is not an alignment', /at least two/.test(threw), threw);

section('the nucleotide figure');
const dna = [
  { name: 'x', seq: 'ATGGCTAGCAAAGG' },
  { name: 'y', seq: 'ATGGCTAGCAAAGG' },
  { name: 'z', seq: 'ATGGCCAGTAAAGG' }
];
const dsvg = B.alnFigureSVG(dna, { nucleotideColour: true });
ok('a nucleotide alignment draws', dsvg.indexOf('>x<') !== -1 && dsvg.indexOf('>A<') !== -1);
ok('RNA shows uracil in place of thymine', (() => {
  const r = B.alnFigureSVG(dna, { molecule: 'rna' });
  return r.indexOf('>U<') !== -1 && r.indexOf('>T<') === -1;
})());
ok('DNA keeps thymine', dsvg.indexOf('>T<') !== -1);
ok('a protein grouping asked for on DNA falls back rather than shading wrongly', (() => {
  const a1 = B.alnFigureSVG(dna, { mode: 'physchem' });
  const a2 = B.alnFigureSVG(dna, { mode: 'strict' });
  return a1 === a2;
})());
ok('purine colouring puts purines and pyrimidines in different colours',
  dsvg.indexOf('#c02626') !== -1 || dsvg.indexOf('#2b4fd8') !== -1);

console.log('\n   ' + pass + ' passed, ' + fail + ' failed\n');
process.exit(fail ? 1 : 0);
