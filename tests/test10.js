/* ===========================================================================
   PlasmidBench — searching, across one construct and across the library.

   Pattern hits are checked against a plain scan written out separately, and
   every hit is checked by reading the bases back out of the sequence: a hit
   must contain what was asked for, on the strand it claims.
   =========================================================================== */
const E = require('../src/engine2.js');
let pass = 0, fail = 0;
function ok(cond, what, extra) {
  if (cond) pass++; else { fail++; console.log('FAIL: ' + what + (extra !== undefined ? '  [' + extra + ']' : '')); }
}
function prng(seed) { let x = seed; return () => (x = (x * 1103515245 + 12345) & 0x7fffffff) / 0x7fffffff; }
const r = prng(2024), rnd = n => Math.floor(r() * n);
const dna = n => { let s = ''; for (let i = 0; i < n; i++) s += 'ACGT'[rnd(4)]; return s; };
const take = (s, h) => h.end > h.start ? s.slice(h.start, h.end) : s.slice(h.start) + s.slice(0, h.end);

/* ---------------------------------------------------- 1. plain patterns */
(function () {
  const s = 'AAGGATCCTTGGATCCAA';
  const h = E.findPattern(s, 'GGATCC');
  ok(h.length === 2, 'both copies of a palindrome are found', h.length);
  ok(h.every(x => x.strand === 1), 'a palindrome is reported once per place, not once per strand',
     JSON.stringify(h.map(x => x.strand)));
  ok(h[0].start === 2 && h[0].end === 8, 'at the right place', JSON.stringify(h[0]));

  const t = 'TTTTGAATTCTTTT';
  ok(E.findPattern(t, 'GAATTC').length === 1, 'EcoRI site found once');
  /* AAACCGGTTTAAA was a badly chosen example: CCGGT sits at 3..8 on the top
     strand, but its own reverse complement ACCGG sits at 2..7 on the top
     strand as well, so a search for either query honestly finds two places.
     ACGGCGTTTAAA has no such accident. */
  const asym = 'AAACGGCATTTAAA';
  const f = E.findPattern(asym, 'CGGCA');
  ok(f.length === 1 && f[0].strand === 1, 'a non-palindrome is found on the strand it reads on',
     JSON.stringify(f));
  const rv = E.findPattern(asym, E.revcomp('CGGCA'));
  ok(rv.length === 1 && rv[0].strand === -1 && rv[0].start === f[0].start,
     'and its reverse complement is found at the same place on the other strand', JSON.stringify(rv[0]));
  const both = E.findPattern('AAACCGGTTTAAA', 'CCGGT');
  ok(both.length === 2, 'a query whose own reverse complement also reads on the top strand ' +
     'is reported at both places', JSON.stringify(both));

  ok(E.findPattern('AAGGATCCTT', 'GGWTCC').length === 2, 'an IUPAC code matches the bases it stands for');
  ok(E.findPattern('AAGGCTCCTT', 'GGWTCC').length === 0, 'and only those bases');
  ok(E.findPattern('ACGT', 'NNNN').length >= 1, 'N matches anything');
  ok(E.findPattern('ACGT', 'ACGTA').length === 0, 'a query longer than the sequence finds nothing');
  ok(E.findPattern('ACGT', '').length === 0, 'an empty query finds nothing');
  ok(E.findPattern('acgtacgt', 'ACGT').length === 2, 'the sequence is read case-insensitively');
})();

/* --------------------------------------- 2. against an independent scan */
(function () {
  const IU = { A: 'A', C: 'C', G: 'G', T: 'T', R: 'AG', Y: 'CT', S: 'CG', W: 'AT', K: 'GT',
               M: 'AC', B: 'CGT', D: 'AGT', H: 'ACT', V: 'ACG', N: 'ACGT' };
  function scan(s, p, circ) {
    const L = s.length, n = p.length, out = [];
    [1, -1].forEach(strand => {
      const q = strand === 1 ? p : E.revcomp(p);
      if (strand === -1 && E.revcomp(p) === p) return;
      const last = circ ? L - 1 : L - n;
      for (let i = 0; i <= last; i++) {
        let okHit = true;
        for (let k = 0; k < n; k++) {
          const b = s[circ ? (i + k) % L : i + k];
          if (IU[q[k]].indexOf(b) < 0) { okHit = false; break; }
        }
        if (okHit) out.push(strand + ':' + i);
      }
    });
    return out.sort();
  }
  let cases = 0, found = 0;
  for (let t = 0; t < 400; t++) {
    const circ = t % 2 === 0;
    const L = 60 + rnd(400), s = dna(L);
    /* half the queries are taken from the sequence so hits are common */
    let p;
    if (r() < 0.6) { const a = rnd(L - 8); p = (s + s).slice(a, a + 4 + rnd(8)); if (r() < 0.5) p = E.revcomp(p); }
    else p = dna(4 + rnd(6));
    if (r() < 0.3) { const k = rnd(p.length); p = p.slice(0, k) + 'RYSWKMN'[rnd(7)] + p.slice(k + 1); }
    const got = E.findPattern(s, p, { circular: circ }).map(h => h.strand + ':' + h.start).sort();
    const want = scan(s, p, circ);
    ok(JSON.stringify(got) === JSON.stringify(want), 'scan case ' + t + (circ ? ' circular' : ' linear'),
       JSON.stringify(got) + ' vs ' + JSON.stringify(want));
    /* and every hit really contains the query */
    E.findPattern(s, p, { circular: circ }).forEach(h => {
      const sub = take(s, h);
      const q = h.strand === 1 ? p : E.revcomp(p);
      let okSub = sub.length === p.length;
      for (let k = 0; okSub && k < q.length; k++) if (IU[q[k]].indexOf(sub[k]) < 0) okSub = false;
      ok(okSub, 'scan case ' + t + ': the hit really contains the query', sub + ' vs ' + q);
      found++;
    });
    cases++;
  }
  console.log('  ' + cases + ' sequences scanned, ' + found + ' hits read back');
})();

/* ------------------------------------------------- 3. peptide searching */
(function () {
  let checked = 0;
  for (let t = 0; t < 200; t++) {
    const circ = t % 2 === 1;
    const L = 300 + rnd(600), s = dna(L);
    const dir = r() < 0.5 ? 1 : -1;
    const base = dir === 1 ? s : E.revcomp(s);
    const f = rnd(3), at = f + 3 * rnd(20);
    const pep = E.translate(base.slice(at, at + 24));
    if (/[*X]/.test(pep) || pep.length < 8) continue;
    const hits = E.findPeptide(s, pep, { circular: circ });
    ok(hits.length >= 1, 'peptide case ' + t + ': the peptide is found');
    hits.forEach(h => {
      const nt = take(s, h);
      const aa = E.translate(h.strand === 1 ? nt : E.revcomp(nt));
      ok(aa === pep, 'peptide case ' + t + ': the hit translates back to the peptide', aa + ' vs ' + pep);
    });
    checked++;
  }
  ok(checked > 120, 'enough peptide cases ran', checked);
  console.log('  ' + checked + ' peptides found and translated back');

  /* a peptide spanning base 1 of a circle */
  const c = dna(300);
  const rot = c.slice(290) + c.slice(0, 290);
  const pep = E.translate(rot.slice(0, 21));
  if (!/[*X]/.test(pep)) {
    const h = E.findPeptide(c, pep, { circular: true });
    ok(h.some(x => x.wraps), 'a peptide across base 1 is found there', JSON.stringify(h));
    h.filter(x => x.wraps).forEach(x => {
      const nt = take(c, x);
      ok(E.translate(x.strand === 1 ? nt : E.revcomp(nt)) === pep, 'and it translates back');
    });
  }
  ok(E.findPeptide('ACGTACGTACGT', 'M').length === 0, 'a one-residue query is refused');
})();

/* --------------------------------------------- 4. what a query is taken for */
(function () {
  const enz = [{ name: 'BamHI', site: 'GGATCC' }, { name: 'EcoRI', site: 'GAATTC' }];
  const k = q => E.classifyQuery(q, { enzymes: enz }).kinds;
  ok(k('GGATCC').indexOf('dna') >= 0 && k('GGATCC').indexOf('text') >= 0, 'bases are read as bases and as text');
  ok(k('BamHI').indexOf('enzyme') >= 0, 'an enzyme name is recognised');
  ok(k('bamhi').indexOf('enzyme') >= 0, 'whatever its case');
  ok(k('MEDYTK').indexOf('protein') >= 0, 'a peptide with E and M is read as protein');
  ok(k('MEDYTK').indexOf('dna') < 0, 'and not as bases');
  ok(k('ACGT').indexOf('protein') < 0, 'ACGT is not read as protein');
  ok(k('kanR').indexOf('text') >= 0, 'a name is read as text');
  ok(k('a').length === 0, 'a single character is not a query');
  ok(k('').length === 0, 'nothing is not a query');
})();

/* ------------------------------------------------ 5. across the library */
(function () {
  const enz = [{ name: 'BamHI', site: 'GGATCC' }];
  const A = { id: 'a', name: 'pVectorOne', seq: 'AAAA' + 'GGATCC' + dna(200) + 'GAATTC' + dna(100), topology: 'circular',
              features: [{ start: 4, end: 10, strand: 1, type: 'misc_feature', label: 'BamHI site', quals: {} },
                         { start: 50, end: 90, strand: -1, type: 'CDS', label: 'kanR', quals: { product: 'aminoglycoside kinase' } }] };
  const B2 = { id: 'b', name: 'pOther', seq: dna(150) + 'GGATCC' + dna(150), topology: 'linear',
               features: [{ start: 10, end: 40, strand: 1, type: 'promoter', label: 'T7 promoter', quals: {} }] };
  const primers = [{ name: 'kanR_check_F', seq: 'GGATCCAAGGTTCCAAGGTT' }, { name: 'other', seq: dna(20) }];

  let res = E.searchLibrary([A, B2], 'GGATCC', { enzymes: enz, primers: primers });
  const inA = res.hits.filter(h => h.recId === 'a'), inB = res.hits.filter(h => h.recId === 'b');
  ok(inA.some(h => h.kind === 'dna' && h.start === 4), 'the bases are found in the first construct');
  ok(inB.some(h => h.kind === 'dna'), 'and in the second');
  ok(inA.some(h => h.kind === 'enzyme'), 'the same query also finds the enzyme sites');
  /* A base query finds bases and enzyme sites. It does not find the feature
     labelled "BamHI site", because that label does not contain the string
     GGATCC — asking it to would mean quietly translating the query into an
     enzyme name and then matching on that, which is a guess. The label match
     is checked below with a query that really is in the label. */
  ok(!inA.some(h => h.kind === 'feature'), 'a base query does not match a label it is not in',
     JSON.stringify(inA.map(h => h.kind)));
  const byName = E.searchLibrary([A, B2], 'BamHI', { enzymes: enz, primers: primers }).hits;
  ok(byName.some(h => h.kind === 'feature' && h.recId === 'a'),
     'and a feature whose label contains the query is found', JSON.stringify(byName.map(h => h.kind)));
  const byProduct = E.searchLibrary([A, B2], 'kinase', { enzymes: enz, primers: primers }).hits;
  ok(byProduct.some(h => h.kind === 'feature' && h.recId === 'a'),
     'as is one whose product qualifier contains it', JSON.stringify(byProduct.map(h => h.kind)));
  ok(res.primers.some(p => p.name === 'kanR_check_F'), 'and a primer containing it');
  ok(res.hits.every(h => h.recName), 'every hit says which construct it is in');

  res = E.searchLibrary([A, B2], 'kanR', { enzymes: enz, primers: primers });
  ok(res.hits.some(h => h.kind === 'feature' && h.recId === 'a'), 'a name finds the feature');
  ok(res.primers.some(p => p.name === 'kanR_check_F'), 'and the primer named after it');
  ok(!res.hits.some(h => h.kind === 'dna'), 'a word that is not bases is not searched as bases');

  res = E.searchLibrary([A, B2], 'aminoglycoside', {});
  ok(res.hits.some(h => h.kind === 'feature'), 'a qualifier is searched too');

  res = E.searchLibrary([A, B2], 'pOther', {});
  ok(res.hits.some(h => h.kind === 'name' && h.recId === 'b'), 'a construct name is found');

  res = E.searchLibrary([A, B2], 'zzzznothing', { enzymes: enz, primers: primers });
  ok(res.total === 0, 'a query that matches nothing finds nothing');

  /* every hit's range can be read back out of the construct it names */
  res = E.searchLibrary([A, B2], 'GGATCC', { enzymes: enz, primers: primers });
  res.hits.filter(h => h.kind === 'dna').forEach(h => {
    const rec = h.recId === 'a' ? A : B2;
    ok(take(rec.seq, h) === 'GGATCC', 'a dna hit reads back as the query', take(rec.seq, h));
  });
})();

console.log('\ntest10 (search): ' + pass + ' passed, ' + fail + ' failed');
process.exit(fail ? 1 : 0);
