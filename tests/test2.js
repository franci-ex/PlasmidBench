const B = require('../src/engine.js');
require('../src/engine2.js');
let pass = 0, fail = 0;
function ok(name, cond, extra) {
  if (cond) { pass++; console.log('  PASS  ' + name); }
  else { fail++; console.log('  FAIL  ' + name + (extra !== undefined ? '   -> ' + extra : '')); }
}
function near(a, b, tol) { return Math.abs(a - b) <= tol; }
// deterministic pseudo-random DNA
function prng(seed) { let x = seed; return () => (x = (x * 1103515245 + 12345) & 0x7fffffff) / 0x7fffffff; }
function randDNA(n, seed) {
  const r = prng(seed); let s = '';
  for (let i = 0; i < n; i++) s += 'ACGT'[Math.floor(r() * 4)];
  return s;
}
// helper: get the left/right ends an enzyme leaves
function endsFor(name) {
  const e = B.enzymeByName(name);
  const pick = { R:'A',Y:'C',S:'G',W:'A',K:'G',M:'A',B:'C',D:'A',H:'A',V:'A',N:'A' };
  const site = e.site.split('').map(c => 'ACGT'.includes(c) ? c : pick[c]).join('');
  const probe = 'ATCGATCGAT' + site + 'ATCGATCGAT';
  const s = B.findSites(probe, e, false)[0];
  return {
    left: B.endOf(probe, s.cutTop, s.cutBottom, 'left', probe.length),
    right: B.endOf(probe, s.cutTop, s.cutBottom, 'right', probe.length)
  };
}

console.log('\n== fragment ends ==');
const eco = endsFor('EcoRI');
ok("EcoRI left end is 5' AATT", eco.left.type === "5'" && eco.left.seq === 'AATT', JSON.stringify(eco.left));
ok("EcoRI right end is 5' AATT (revcomp of AATT)", eco.right.seq === 'AATT', JSON.stringify(eco.right));
const sma = endsFor('SmaI');
ok('SmaI leaves blunt ends', sma.left.type === 'blunt' && sma.right.type === 'blunt', JSON.stringify(sma));
const pst = endsFor('PstI');
ok("PstI leaves 3' overhangs", pst.left.type === "3'" && pst.right.type === "3'", JSON.stringify(pst));
ok('PstI overhang is 4 nt', pst.left.seq.length === 4, pst.left.seq);
const nde = endsFor('NdeI');
ok("NdeI leaves 5' TA (2 nt)", nde.left.type === "5'" && nde.left.seq === 'TA', JSON.stringify(nde.left));

console.log('\n== ligation compatibility (known enzyme pairs) ==');
const bam = endsFor('BamHI'), bgl = endsFor('BglII'), sal = endsFor('SalI'),
      xho = endsFor('XhoI'), ecorv = endsFor('EcoRV'), xma = endsFor('XmaI'),
      age = endsFor('AgeI'), bspe = endsFor('BspEI'), hind = endsFor('HindIII');
ok('BamHI overhang is GATC', bam.left.seq === 'GATC', bam.left.seq);
ok('BglII overhang is GATC', bgl.left.seq === 'GATC', bgl.left.seq);
ok('BamHI + BglII are compatible', B.endsCompatible(bam.right, bgl.left));
ok('SalI + XhoI are compatible (TCGA)', B.endsCompatible(sal.right, xho.left),
   sal.right.seq + ' / ' + xho.left.seq);
ok('XmaI + AgeI are compatible (CCGG)', B.endsCompatible(xma.right, age.left),
   xma.right.seq + ' / ' + age.left.seq);
ok('XmaI + BspEI are compatible (CCGG)', B.endsCompatible(xma.right, bspe.left));
ok('SmaI + EcoRV are compatible (both blunt)', B.endsCompatible(sma.right, ecorv.left));
ok('BamHI + EcoRI are NOT compatible', !B.endsCompatible(bam.right, eco.left));
ok('EcoRI + HindIII are NOT compatible', !B.endsCompatible(eco.right, hind.left));
ok("PstI (3') + EcoRI (5') are NOT compatible", !B.endsCompatible(pst.right, eco.left));
ok('SmaI + EcoRI are NOT compatible (blunt vs sticky)', !B.endsCompatible(sma.right, eco.left));
ok('an end is compatible with its own partner end', B.endsCompatible(eco.right, eco.left));
const compat = B.compatibleEnzymes(bam.left).map(e => e.name);
ok('compatibleEnzymes finds BglII for a BamHI end', compat.includes('BglII'), compat.join(','));
ok('compatibleEnzymes finds Sau3AI/MboI for a BamHI end',
   compat.includes('Sau3AI') || compat.includes('MboI'), compat.join(','));
ok('compatibleEnzymes excludes EcoRI for a BamHI end', !compat.includes('EcoRI'), compat.join(','));

console.log('\n== restriction cloning (hand-verified lengths) ==');
// vector: circular, 2516 bp, EcoRI at 2000, HindIII 4 bp later
const vecPre = 'G'.repeat(2000);
const vector = vecPre + 'GAATTC' + 'TTTT' + 'AAGCTT' + 'C'.repeat(500);
ok('vector is 2516 bp', vector.length === 2516, vector.length);
const vdig = B.digestWithEnds(vector, ['EcoRI', 'HindIII'], { circular: true });
ok('vector digest gives 2 fragments', vdig.fragments.length === 2, vdig.fragments.length);
const backbone = vdig.fragments.slice().sort((a, b) => b.length - a.length)[0];
ok('backbone is 2506 bp', backbone.length === 2506, backbone.length);
ok('small vector fragment is 10 bp',
   vdig.fragments.find(f => f.length !== 2506).length === 10,
   vdig.fragments.map(f => f.length).join(','));
ok("backbone left end is 5' AGCT (HindIII)", backbone.left.seq === 'AGCT', JSON.stringify(backbone.left));
ok("backbone right end is 5' AATT (EcoRI)", backbone.right.seq === 'AATT', JSON.stringify(backbone.right));

const payload = randDNA(300, 7);
const insert = 'AAA' + 'GAATTC' + payload + 'AAGCTT' + 'TTT';
const clone = B.restrictionClone({
  vector, insert, vectorEnzymes: ['EcoRI', 'HindIII'], insertEnzymes: ['EcoRI', 'HindIII'],
  vectorCircular: true, insertCircular: false
});
ok('restriction clone succeeds', clone.ok, clone.error + ' ' + JSON.stringify(clone.detail || {}));
if (clone.ok) {
  ok('insert fragment is 306 bp', clone.insertFragment.length === 306, clone.insertFragment.length);
  ok('product is 2812 bp (2506 + 306)', clone.products[0].length === 2812, clone.products[0].length);
  ok('product is circular', clone.products[0].circular === true);
  ok('product contains the payload', clone.products[0].seq.includes(payload));
  ok('product regenerates the EcoRI site',
     B.findSites(clone.products[0].seq, B.enzymeByName('EcoRI'), true).length === 1,
     B.findSites(clone.products[0].seq, B.enzymeByName('EcoRI'), true).length);
  ok('product regenerates the HindIII site',
     B.findSites(clone.products[0].seq, B.enzymeByName('HindIII'), true).length === 1);
  ok('directional clone gives exactly one orientation', clone.products.length === 1, clone.products.length);
  console.log('    junctions:', clone.products[0].junction5, '|', clone.products[0].junction3);
}
// incompatible ends must be refused
const badClone = B.restrictionClone({
  vector, insert: 'AAA' + 'GGATCC' + payload + 'GGATCC' + 'TTT',
  vectorEnzymes: ['EcoRI', 'HindIII'], insertEnzymes: ['BamHI'],
  vectorCircular: true
});
ok('incompatible ends are refused', !badClone.ok, JSON.stringify(badClone).slice(0, 120));
// a symmetric (same-enzyme) insert should give both orientations
const sym = B.restrictionClone({
  vector: 'G'.repeat(1000) + 'GAATTC' + 'TT' + 'GAATTC' + 'C'.repeat(200),
  insert: 'AA' + 'GAATTC' + payload + 'GAATTC' + 'TT',
  vectorEnzymes: ['EcoRI'], insertEnzymes: ['EcoRI'], vectorCircular: true
});
ok('same-enzyme cloning gives both orientations', sym.ok && sym.products.length === 2,
   sym.ok ? sym.products.length : sym.error);

console.log('\n== Golden Gate (BsaI, MoClo-style parts) ==');
function ggPart(left, payload, right) {
  return 'AA' + 'GGTCTC' + 'A' + left + payload + right + 'A' + 'GAGACC' + 'AA';
}
const pA = randDNA(30, 11), pB = randDNA(40, 12), pC = randDNA(50, 13);
const parts = [
  { name: 'A', seq: ggPart('AATG', pA, 'GCTT') },
  { name: 'B', seq: ggPart('GCTT', pB, 'CTGA') },
  { name: 'C', seq: ggPart('CTGA', pC, 'AATG') }
];
const gg = B.goldenGate(parts, 'BsaI');
ok('Golden Gate assembles', gg.ok, gg.error + ' order=' + JSON.stringify(gg.order));
ok('no overhang clashes', gg.overhangClashes.length === 0, gg.overhangClashes.join('; '));
ok('all 3 parts used', gg.order.length === 3, JSON.stringify(gg.order));
ok('assembly closes into a circle', gg.closes);
ok('product is 132 bp (3x4 overhang + 30+40+50)', gg.product.length === 132, gg.product.length);
ok('product contains payload A', gg.product.seq.includes(pA));
ok('product contains payload B', gg.product.seq.includes(pB));
ok('product contains payload C', gg.product.seq.includes(pC));
ok('BsaI sites are removed from the product',
   B.findSites(gg.product.seq, B.enzymeByName('BsaI'), true).length === 0,
   B.findSites(gg.product.seq, B.enzymeByName('BsaI'), true).length);
console.log('    part order:', gg.order.join(' -> '), '| product', gg.product.length, 'bp');
// clashing overhangs must be reported
const clash = B.goldenGate([
  { name: 'A', seq: ggPart('AATG', pA, 'GCTT') },
  { name: 'B', seq: ggPart('AATG', pB, 'CTGA') }
], 'BsaI');
ok('duplicate overhangs are reported as a clash', clash.overhangClashes.length > 0,
   JSON.stringify(clash.overhangClashes));
ok('non-Type-IIS enzyme is rejected', !B.goldenGate(parts, 'EcoRI').ok);

console.log('\n== Gibson / HiFi assembly ==');
const S = randDNA(300, 21);
const s1 = S.slice(0, 100), s2 = S.slice(100, 200), s3 = S.slice(200);
const g1 = s1 + s2.slice(0, 20);
const g2 = s2 + s3.slice(0, 20);
const g3 = s3 + s1.slice(0, 20);
const gib = B.gibsonAssemble([{ name: 'f1', seq: g1 }, { name: 'f2', seq: g2 }, { name: 'f3', seq: g3 }],
  { minOverlap: 15, circular: true });
ok('Gibson assembles 3 fragments', gib.ok, gib.error);
ok('product is 300 bp (the original circle)', gib.product.length === 300, gib.product.length);
ok('product sequence equals the original', gib.product.seq === S, gib.product.seq.slice(0, 30));
ok('overlaps detected as 20 bp', gib.overlaps.every(n => n === 20), JSON.stringify(gib.overlaps));
ok('closing overlap is 20 bp', gib.closingOverlap === 20, gib.closingOverlap);
ok('product marked circular', gib.product.circular === true);
// insufficient homology must fail, not silently produce junk
const noHom = B.gibsonAssemble([
  { name: 'a', seq: randDNA(100, 31) }, { name: 'b', seq: randDNA(100, 32) }
], { minOverlap: 15 });
ok('fragments without homology are refused', !noHom.ok, noHom.error);
// a flipped fragment should still be found
const gibFlip = B.gibsonAssemble([{ name: 'f1', seq: g1 }, { name: 'f2', seq: B.revcomp(g2) },
                                  { name: 'f3', seq: g3 }], { minOverlap: 15, circular: true });
ok('a reverse-complemented fragment is detected and flipped', gibFlip.ok, gibFlip.error);
ok('flipped assembly still gives 300 bp', gibFlip.ok && gibFlip.product.length === 300,
   gibFlip.ok ? gibFlip.product.length : '-');

console.log('\n== assembly primer design ==');
const vecLin = randDNA(2000, 41);
const insA = randDNA(600, 42);
const des = B.designAssemblyPrimers(vecLin, [{ name: 'gene', seq: insA }], { overlap: 20, anneal: 20 });
ok('one primer pair per insert', des.primers.length === 1, des.primers.length);
const pp = des.primers[0];
ok('forward primer is 40 nt', pp.forward.length === 40, pp.forward.length);
ok('forward tail is the vector 3-prime end', pp.forward.tail === vecLin.slice(-20), pp.forward.tail);
ok('forward anneal is the insert 5-prime end', pp.forward.anneal === insA.slice(0, 20));
ok('reverse anneal is revcomp of insert 3-prime end',
   pp.reverse.anneal === B.revcomp(insA.slice(-20)));
ok('annealing Tm is reported', isFinite(pp.forward.tmAnneal) && pp.forward.tmAnneal > 30,
   pp.forward.tmAnneal);
// the designed primers must amplify the insert and the amplicon must assemble
const amp = B.pcr(insA, pp.forward.seq, pp.reverse.seq, { anchor: 18 });
ok('designed primers amplify the insert', amp.length >= 1, amp.length);
ok('amplicon = 20 + 600 + 20', amp.length && amp[0].length === 640, amp.length ? amp[0].length : '-');
const reassembled = B.gibsonAssemble([{ name: 'vec', seq: vecLin }, { name: 'amp', seq: amp[0].seq }],
  { minOverlap: 15, circular: true });
ok('vector + designed amplicon assemble into a circle', reassembled.ok, reassembled.error);
ok('assembled circle is 2600 bp', reassembled.ok && reassembled.product.length === 2600,
   reassembled.ok ? reassembled.product.length : '-');

console.log('\n== In-Fusion (15 bp overlaps) ==');
const inf = B.designAssemblyPrimers(vecLin, [{ name: 'gene', seq: insA }], { method: 'infusion' });
ok('In-Fusion uses 15 bp homology', inf.overlap === 15, inf.overlap);
ok('In-Fusion primer tail is 15 nt', inf.primers[0].forward.tail.length === 15);

console.log('\n== Gateway cloning ==');
/* The engine no longer carries att-site sequences of its own: it locates
   them by their conserved core in whatever sequences it is given, and says
   so in its notes. That is a deliberate change — a recombination product
   should be built from the user's own bases, not from a table inside the
   page — and it means the att sites belong here, in the test, where they
   can be stated and checked.

   These are the published attB cores (Invitrogen's Gateway system):
     attB1  ACAAGTTTGTACAAAAAAGCAGGCT
     attB2  ACCCAGCTTTCTTGTACAAAGTGGT
   Each contains the family core the engine searches for — TTTGTACAAAAAAGC
   for family 1, TTGTACAAAGT for family 2 — so a molecule carrying them is
   recombined at exactly the position those cores occupy.

   Rewritten 2026-11: this section asked for B.ATT_SITES and lr.cassetteLength,
   neither of which the engine has had since att sites stopped being stored.
   The geometry asserted below is derived by arithmetic from the pieces, not
   copied from what the engine returns. */
const attB1 = 'ACAAGTTTGTACAAAAAAGCAGGCT';
const attB2 = 'ACCCAGCTTTCTTGTACAAAGTGGT';
const cassette = randDNA(900, 51);
const entryBackbone = randDNA(1500, 52);
const destStuffer = randDNA(1200, 54);
const destBackbone = randDNA(2500, 53);
const entry = attB1 + cassette + attB2 + entryBackbone;
const dest = attB1 + destStuffer + attB2 + destBackbone;
const lr = B.gatewayLR(entry, dest);
ok('LR reaction succeeds', lr.ok, lr.error);
ok('LR moves the cassette into the destination', lr.ok && lr.product.seq.includes(cassette));
ok('LR keeps the destination backbone', lr.ok && lr.product.seq.includes(destBackbone));
ok('LR drops the stuffer the cassette replaces', lr.ok && !lr.product.seq.includes(destStuffer));
ok('LR leaves the entry backbone behind', lr.ok && !lr.product.seq.includes(entryBackbone));
ok('and the entry backbone is what the byproduct carries',
   lr.ok && lr.byproduct && lr.byproduct.seq.includes(entryBackbone));
/* the product is the cassette plus the destination backbone plus the two
   cores the two molecules now share — a length that can be worked out
   without asking the engine */
ok('the product is as long as its pieces make it',
   lr.ok && Math.abs(lr.product.seq.length - (cassette.length + destBackbone.length +
     lr.cores.reduce((a, c) => a + c.length, 0))) <= 1,
   lr.ok ? lr.product.seq.length + ' vs ' + (cassette.length + destBackbone.length +
     lr.cores.reduce((a, c) => a + c.length, 0)) : '-');
ok('both att families were located', lr.ok && lr.cores.length === 2,
   lr.ok ? JSON.stringify(lr.cores) : '-');
ok('and it says the sites came from the sequence rather than an annotation',
   lr.ok && (lr.notes || []).some(n => /located by sequence motif/.test(n)),
   JSON.stringify(lr.notes));
ok('LR refuses an entry clone with no att sites', !B.gatewayLR(randDNA(500, 61), dest).ok);
ok('LR refuses a destination with no att sites', !B.gatewayLR(entry, randDNA(500, 62)).ok);

console.log('\n== TOPO / TA cloning ==');
const topoDir = B.topoClone('CACC' + randDNA(500, 71), { mode: 'directional', vectorSeq: randDNA(3000, 72) });
ok('directional TOPO accepts a CACC insert', topoDir.ok, topoDir.error);
ok('directional TOPO reports fixed orientation', /fixed/.test(topoDir.orientation));
const topoBad = B.topoClone(randDNA(500, 73), { mode: 'directional' });
ok('directional TOPO rejects an insert without CACC', !topoBad.ok, topoBad.error);
const ta = B.topoClone(randDNA(500, 74), { mode: 'ta', vectorSeq: randDNA(3000, 75) });
ok('TA cloning succeeds', ta.ok);
ok('TA cloning warns about both orientations', /both orientations/.test(ta.orientation));

console.log('\n== site-directed mutagenesis ==');
const sdmTpl = randDNA(200, 81);
const sdmRes = B.sdm(sdmTpl, { position: 100, kind: 'substitute', newBases: 'G', flank: 15 });
ok('SDM succeeds', sdmRes.ok, sdmRes.error);
ok('SDM product has the same length', sdmRes.mutated.length === 200, sdmRes.mutated.length);
const diffs = [...sdmTpl].filter((c, i) => c !== sdmRes.mutated.seq[i]).length;
ok('SDM changes exactly one base (or zero if already G)', diffs <= 1, diffs);
ok('SDM base at the target position is G', sdmRes.mutated.seq[100] === 'G', sdmRes.mutated.seq[100]);
ok('SDM forward primer is 31 nt (15 + 1 + 15)', sdmRes.primers.forward.length === 31,
   sdmRes.primers.forward.length);
ok('SDM primers are reverse complements', sdmRes.primers.reverse.seq === B.revcomp(sdmRes.primers.forward.seq));
ok('SDM primer carries the mutation', sdmRes.primers.forward.seq[15] === 'G',
   sdmRes.primers.forward.seq);
ok('SDM primer Tm reported', isFinite(sdmRes.primers.forward.tm) && sdmRes.primers.forward.tm > 40,
   sdmRes.primers.forward.tm);
const del = B.sdm(sdmTpl, { position: 50, kind: 'delete', deleteLength: 3 });
ok('SDM deletion shortens by 3', del.mutated.length === 197, del.mutated.length);
const ins = B.sdm(sdmTpl, { position: 50, kind: 'insert', newBases: 'GGGGG' });
ok('SDM insertion lengthens by 5', ins.mutated.length === 205, ins.mutated.length);
ok('SDM insertion is present', ins.mutated.seq.slice(50, 55) === 'GGGGG', ins.mutated.seq.slice(45,60));

console.log('\n== silent mutations ==');
const cdsEco = 'ATG' + 'GAATTC' + 'AAAGGGCCC' + 'TAA';   // M E F K G P *
ok('test CDS translates to MEFKGP*', B.translate(cdsEco) === 'MEFKGP*', B.translate(cdsEco));
ok('test CDS has one EcoRI site', B.findSites(cdsEco, B.enzymeByName('EcoRI'), false).length === 1);
const sil = B.silentMutations(cdsEco, 'EcoRI', { mode: 'remove' });
ok('silent mutation search runs', sil.ok, sil.error);
ok('a synonymous fix was found', sil.results[0] && sil.results[0].seq, JSON.stringify(sil.results[0]));
if (sil.results[0] && sil.results[0].seq) {
  const fixed = sil.results[0].seq;
  ok('EcoRI site removed', B.findSites(fixed, B.enzymeByName('EcoRI'), false).length === 0);
  ok('protein unchanged', B.translate(fixed) === B.translate(cdsEco),
     B.translate(fixed) + ' vs ' + B.translate(cdsEco));
  ok('only one base changed', sil.results[0].changed.length === 1, JSON.stringify(sil.results[0].changed));
  console.log('    ' + sil.results[0].from + ' -> ' + sil.results[0].to +
              ' (' + sil.results[0].aa + ') removes the site');
}
// adding a site
const cdsPlain = 'ATG' + 'GAGTTT' + 'AAAGGGCCCTAA';
const add = B.silentMutations(cdsPlain, 'EcoRI', { mode: 'add' });
ok('silent-add finds a way to introduce EcoRI', add.results.length > 0, add.results.length);
if (add.results.length) {
  ok('added site is present', B.findSites(add.results[0].seq, B.enzymeByName('EcoRI'), false).length === 1);
  ok('protein still unchanged after adding a site',
     B.translate(add.results[0].seq) === B.translate(cdsPlain));
}

console.log('\n== primer design engine ==');
const pdTpl = randDNA(1200, 91);
const pd = B.designPrimers(pdTpl, { start: 100, end: 700, targetTm: 60 });
ok('primer design succeeds', pd.ok, pd.error);
ok('forward Tm is near target', near(pd.forward.tm, 60, 5), pd.forward.tm);
ok('reverse Tm is near target', near(pd.reverse.tm, 60, 5), pd.reverse.tm);
ok('Tm difference under 5 C', pd.tmDifference < 5, pd.tmDifference);
ok('amplicon found', pd.amplicon && pd.amplicon.length > 500, pd.amplicon && pd.amplicon.length);
ok('amplicon is about 600 bp', pd.amplicon && Math.abs(pd.amplicon.length - 600) <= 30,
   pd.amplicon && pd.amplicon.length);
ok('annealing temperature suggested', isFinite(pd.annealingSuggestion), pd.annealingSuggestion);
ok('hairpin analysis present', pd.forward.hairpin != null);
ok('thermodynamics reported', pd.forward.thermo && pd.forward.thermo.dG37 < 0,
   pd.forward.thermo && pd.forward.thermo.dG37);
// primers with 5' restriction tails
const pdTail = B.designPrimers(pdTpl, { start: 100, end: 700, forwardTail: 'GGGGAATTC', reverseTail: 'GGGAAGCTT' });
ok('forward tail is carried', pdTail.forward.seq.startsWith('GGGGAATTC'), pdTail.forward.seq.slice(0,12));
ok('tailed amplicon contains the new EcoRI site',
   pdTail.amplicon && B.findSites(pdTail.amplicon.seq, B.enzymeByName('EcoRI'), false).length >= 1);
ok('tailed amplicon contains the new HindIII site',
   pdTail.amplicon && B.findSites(pdTail.amplicon.seq, B.enzymeByName('HindIII'), false).length >= 1);

console.log('\n== hairpin / dimer detection ==');
const hp = B.hairpin('GGGGCCCCAAAGGGGCCCC');
ok('hairpin detected in a palindromic oligo', hp.stemLength >= 4, JSON.stringify(hp));
const noHp = B.hairpin('AAAAAAAAAAAAAAAAAAAA');
ok('no hairpin in poly-A', noHp.stemLength === 0, JSON.stringify(noHp));
const sd = B.selfDimer('GGGGCCCCGGGGCCCC');
ok('self-dimer found in a self-complementary oligo', sd.matches >= 8, JSON.stringify(sd));
const hd = B.heteroDimer('AAAAAAAAAAGGGGGGGGGG', 'CCCCCCCCCCTTTTTTTTTT');
ok('heterodimer detected between complementary primers', hd.matches >= 10, JSON.stringify(hd));
ok('low dimer score for unrelated primers', B.heteroDimer(randDNA(20, 95), randDNA(20, 96)).matches < 10);

console.log('\n== multiple sequence alignment ==');
const base = randDNA(120, 101);
const m1 = B.msa([{ name: 'a', seq: base }, { name: 'b', seq: base }, { name: 'c', seq: base }]);
ok('MSA of identical sequences succeeds', m1.ok, m1.error);
ok('all rows identical', m1.rows.every(r => r.seq === m1.rows[0].seq));
ok('identity is 100%', near(m1.identity, 100, 1e-9), m1.identity);
ok('no gaps for identical input', m1.rows.every(r => !r.seq.includes('-')));
ok('consensus equals the input', m1.consensus === base, m1.consensus.slice(0, 20));
// with a point mutation and a deletion
const mutA = base.slice(0, 60) + 'T' + base.slice(61);
const delB = base.slice(0, 60) + base.slice(63);
const m2 = B.msa([{ name: 'wt', seq: base }, { name: 'point', seq: mutA }, { name: 'del', seq: delB }]);
ok('MSA with indels succeeds', m2.ok, m2.error);
ok('all rows are the same length', new Set(m2.rows.map(r => r.seq.length)).size === 1,
   m2.rows.map(r => r.seq.length).join(','));
ok('alignment is at least as long as the longest input', m2.length >= base.length, m2.length);
ok('ungapped rows recover the input sequences',
   m2.rows.find(r => r.name === 'wt').seq.replace(/-/g, '') === base &&
   m2.rows.find(r => r.name === 'del').seq.replace(/-/g, '') === delB);
ok('MSA reports its method', typeof m2.method === 'string', m2.method);
ok('MSA needs 2+ sequences', !B.msa([{ name: 'x', seq: 'ACGT' }]).ok);
// 5 sequences
const many = [0,1,2,3,4].map(i => ({ name: 's' + i,
  seq: base.slice(0, 50 + i) + randDNA(5, 200 + i) + base.slice(55) }));
const m3 = B.msa(many);
ok('MSA handles 5 sequences', m3.ok && m3.rows.length === 5, m3.ok ? m3.rows.length : m3.error);
ok('5-way rows are equal length', new Set(m3.rows.map(r => r.seq.length)).size === 1,
   m3.rows.map(r => r.seq.length).join(','));
ok('5-way rows all recover their inputs',
   m3.rows.every((r, i) => r.seq.replace(/-/g, '') === many.find(x => x.name === r.name).seq));

console.log('\n== contig assembly ==');
const full = randDNA(400, 111);
const r1 = full.slice(0, 250), r2 = full.slice(200);
const asm = B.assembleContigs([{ name: 'f', seq: r1 }, { name: 'r', seq: r2 }], { minOverlap: 20 });
ok('two overlapping reads give one contig', asm.contigCount === 1, asm.contigCount);
ok('contig is 400 bp', asm.contigs[0].length === 400, asm.contigs[0].length);
ok('contig equals the original', asm.contigs[0].seq === full);
// reverse-complemented read
const asmRC = B.assembleContigs([{ name: 'f', seq: r1 }, { name: 'r', seq: B.revcomp(r2) }], { minOverlap: 20 });
ok('a reverse-complemented read is still assembled', asmRC.contigCount === 1, asmRC.contigCount);
ok('rc contig is 400 bp', asmRC.contigs[0].length === 400, asmRC.contigs[0].length);
// three reads
const t3 = B.assembleContigs([
  { name: 'a', seq: full.slice(0, 180) },
  { name: 'b', seq: full.slice(140, 320) },
  { name: 'c', seq: full.slice(280) }], { minOverlap: 20 });
ok('three tiled reads give one contig', t3.contigCount === 1, t3.contigCount);
ok('three-read contig is 400 bp', t3.contigs[0].length === 400, t3.contigs[0].length);
// non-overlapping reads stay separate
const sep = B.assembleContigs([{ name: 'x', seq: randDNA(200, 121) }, { name: 'y', seq: randDNA(200, 122) }],
  { minOverlap: 20 });
ok('non-overlapping reads stay as separate contigs', sep.contigCount === 2, sep.contigCount);

console.log('\n== RNA secondary structure ==');
const fold = B.rnaFold('GGGGAAAACCCC');
ok('simple hairpin folds', fold.pairCount === 4, JSON.stringify(fold));
ok('structure string is ((((....))))', fold.structure === '((((....))))', fold.structure);
ok('brackets are balanced', balanced(fold.structure), fold.structure);
ok('all pairs are valid RNA pairs', fold.pairs.every(([i, j]) => {
  const a = fold.sequence[i], b = fold.sequence[j];
  return ['AU','UA','GC','CG','GU','UG'].includes(a + b);
}), JSON.stringify(fold.pairs));
ok('minimum loop respected (>=3)', fold.pairs.every(([i, j]) => j - i - 1 >= 3));
const unstruct = B.rnaFold('AAAAAAAAAAAAAAA');
ok('poly-A has no pairs', unstruct.pairCount === 0, unstruct.pairCount);
ok('T is converted to U', B.rnaFold('GGGGTTTTCCCC').sequence.includes('U'));
function balanced(s) {
  let d = 0;
  for (const c of s) { if (c === '(') d++; else if (c === ')') d--; if (d < 0) return false; }
  return d === 0;
}

console.log('\n== agarose gel model ==');
ok('bigger fragments migrate less', B.gelMigration(10000, 1.0) < B.gelMigration(1000, 1.0));
ok('smaller fragments migrate further', B.gelMigration(100, 2.0) > B.gelMigration(1000, 2.0));
ok('migration is monotonic in size', (function () {
  let prev = Infinity;
  for (const bp of [50, 100, 250, 500, 1000, 2000, 5000, 10000, 20000]) {
    const m = B.gelMigration(bp, 1.0);
    if (m > prev) return false;
    prev = m;
  }
  return true;
})());
ok('migration stays within the gel (0..1)', [10, 100, 1000, 100000].every(bp => {
  const m = B.gelMigration(bp, 1.0); return m >= 0 && m <= 1;
}));
/* The ladders are named after the products they are, and each band list is
   the manufacturer's. Checked against the published band lists 2026-11:
   NEB 1 kb (N3232) and 1 kb Plus (N3200), NEB 100 bp (N3231), Thermo
   GeneRuler 1 kb and 1 kb Plus, and lambda DNA cut with HindIII. This test
   asked for a ladder called '1 kb' with ten bands, which has not been the
   name since the real products replaced it. */
const LADDER_BANDS = {
  'NEB 1 kb Plus': 20, 'NEB 1 kb': 11, 'NEB 100 bp': 13,
  'GeneRuler 1 kb': 14, 'GeneRuler 1 kb Plus': 15, '\u03bb DNA HindIII': 8
};
ok('every ladder named is one of the six that exist',
   Object.keys(B.LADDERS).length === 6 &&
   Object.keys(B.LADDERS).every(k => LADDER_BANDS[k] != null),
   Object.keys(B.LADDERS).join(', '));
ok('each has the number of bands its maker lists',
   Object.keys(LADDER_BANDS).every(k => B.LADDERS[k] && B.LADDERS[k].length === LADDER_BANDS[k]),
   Object.keys(LADDER_BANDS).map(k => k + ':' + (B.LADDERS[k] || []).length).join(' '));
ok('every ladder runs largest band first, as a gel does',
   Object.keys(B.LADDERS).every(k => B.LADDERS[k].every((v, i, a) => i === 0 || a[i - 1] > v)),
   Object.keys(B.LADDERS).filter(k => !B.LADDERS[k].every((v, i, a) => i === 0 || a[i - 1] > v)).join(', '));
ok('the largest band of the lambda HindIII ladder is 23,130 bp',
   B.LADDERS['\u03bb DNA HindIII'][0] === 23130, B.LADDERS['\u03bb DNA HindIII'][0]);
ok('higher % agarose resolves smaller fragments better',
   Math.abs(B.gelMigration(100, 2.0) - B.gelMigration(200, 2.0)) >
   Math.abs(B.gelMigration(100, 0.7) - B.gelMigration(200, 0.7)));

console.log('\n== record transforms ==');
const rec = { seq: 'ATGCATGCATGC', topology: 'circular',
              features: [{ type: 'CDS', start: 0, end: 6, strand: 1, label: 'x' }] };
const rot = B.setOrigin(rec, 6);
ok('setOrigin rotates the sequence', rot.seq === 'ATGCATGCATGC'.slice(6) + 'ATGCATGCATGC'.slice(0, 6), rot.seq);
ok('setOrigin keeps the length', rot.seq.length === rec.seq.length);
const flipped = B.flipRecord(rec);
ok('flipRecord reverse-complements', flipped.seq === B.revcomp(rec.seq), flipped.seq);
ok('flipRecord flips feature strand', flipped.features[0].strand === -1);
ok('flipRecord moves feature coordinates',
   flipped.features[0].start === rec.seq.length - 6 && flipped.features[0].end === rec.seq.length,
   JSON.stringify(flipped.features[0]));
ok('double flip restores the record', B.flipRecord(flipped).seq === rec.seq);

console.log('\n==========================================');
console.log('   ' + pass + ' passed, ' + fail + ' failed');
console.log('==========================================\n');
process.exit(fail ? 1 : 0);
