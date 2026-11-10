const B = require('../src/engine.js');
let pass = 0, fail = 0;
function ok(name, cond, extra) {
  if (cond) { pass++; console.log('  PASS  ' + name); }
  else { fail++; console.log('  FAIL  ' + name + (extra !== undefined ? '   -> ' + extra : '')); }
}
function near(a, b, tol) { return Math.abs(a - b) <= tol; }

console.log('\n== genetic code tables ==');
for (const t of B.geneticCodes()) {
  // each table's aa + starts strings must be 64 chars
}
const codes = { 1:'Standard', 11:'Bacterial' };
ok('ATG -> M (table 1)', B.translate('ATG') === 'M', B.translate('ATG'));
ok('TGA -> * (table 1)', B.translate('TGA') === '*', B.translate('TGA'));
ok('TGG -> W', B.translate('TGG') === 'W', B.translate('TGG'));
ok('AAA -> K', B.translate('AAA') === 'K', B.translate('AAA'));
ok('TTT -> F', B.translate('TTT') === 'F', B.translate('TTT'));
ok('GGG -> G', B.translate('GGG') === 'G', B.translate('GGG'));
ok('CGA -> R', B.translate('CGA') === 'R', B.translate('CGA'));
ok('AGT -> S', B.translate('AGT') === 'S', B.translate('AGT'));
ok('GTG -> V', B.translate('GTG') === 'V', B.translate('GTG'));
ok('CCC -> P', B.translate('CCC') === 'P', B.translate('CCC'));
ok('CAT -> H', B.translate('CAT') === 'H', B.translate('CAT'));
ok('GAT -> D', B.translate('GAT') === 'D', B.translate('GAT'));
ok('GAA -> E', B.translate('GAA') === 'E', B.translate('GAA'));
ok('TGT -> C', B.translate('TGT') === 'C', B.translate('TGT'));
ok('ATT -> I', B.translate('ATT') === 'I', B.translate('ATT'));
ok('ACG -> T', B.translate('ACG') === 'T', B.translate('ACG'));
ok('AAT -> N', B.translate('AAT') === 'N', B.translate('AAT'));
ok('CAA -> Q', B.translate('CAA') === 'Q', B.translate('CAA'));
ok('TAT -> Y', B.translate('TAT') === 'Y', B.translate('TAT'));
ok('GCA -> A', B.translate('GCA') === 'A', B.translate('GCA'));
ok('CTG -> L', B.translate('CTG') === 'L', B.translate('CTG'));
// mito table 2: TGA = W, AGA/AGG = stop, ATA = M
ok('table2 TGA -> W', B.translate('TGA', {table:2}) === 'W', B.translate('TGA',{table:2}));
ok('table2 AGA -> *', B.translate('AGA', {table:2}) === '*', B.translate('AGA',{table:2}));
ok('table2 ATA -> M', B.translate('ATA', {table:2}) === 'M', B.translate('ATA',{table:2}));
// table 6 ciliate: TAA and TAG = Q
ok('table6 TAA -> Q', B.translate('TAA', {table:6}) === 'Q', B.translate('TAA',{table:6}));
ok('table6 TAG -> Q', B.translate('TAG', {table:6}) === 'Q', B.translate('TAG',{table:6}));
// table 5 invertebrate mito: AGA = S
ok('table5 AGA -> S', B.translate('AGA', {table:5}) === 'S', B.translate('AGA',{table:5}));
// EGFP N-terminus: ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGG -> MVSKGEELFTG
ok('EGFP N-term translation',
   B.translate('ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGG') === 'MVSKGEELFTG',
   B.translate('ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGG'));
// His6
ok('His6 tag translates to HHHHHH',
   B.translate('CATCATCATCATCATCAT') === 'HHHHHH', B.translate('CATCATCATCATCATCAT'));
// FLAG DYKDDDDK
ok('FLAG tag translates to DYKDDDDK',
   B.translate('GATTACAAGGATGACGACGATAAG') === 'DYKDDDDK', B.translate('GATTACAAGGATGACGACGATAAG'));
// HA tag YPYDVPDYA
ok('HA tag translates to YPYDVPDYA',
   B.translate('TACCCATACGATGTTCCAGATTACGCT') === 'YPYDVPDYA', B.translate('TACCCATACGATGTTCCAGATTACGCT'));
// Myc EQKLISEEDL
ok('Myc tag translates to EQKLISEEDL',
   B.translate('GAACAAAAACTCATCTCAGAAGAGGATCTG') === 'EQKLISEEDL', B.translate('GAACAAAAACTCATCTCAGAAGAGGATCTG'));
// ambiguity: GCN -> A (all alanine)
ok('GCN -> A (unambiguous product)', B.translate('GCN') === 'A', B.translate('GCN'));
ok('NNN -> X', B.translate('NNN') === 'X', B.translate('NNN'));
// toStop
ok('toStop truncates', B.translate('ATGAAATGACCC', {toStop:true}) === 'MK', B.translate('ATGAAATGACCC',{toStop:true}));

console.log('\n== reverse complement ==');
ok('revcomp ACGT', B.revcomp('ACGT') === 'ACGT', B.revcomp('ACGT'));
ok('revcomp AAAA -> TTTT', B.revcomp('AAAA') === 'TTTT', B.revcomp('AAAA'));
ok('revcomp GAATTC palindrome', B.revcomp('GAATTC') === 'GAATTC');
ok('revcomp IUPAC W->W', B.revcomp('W') === 'W', B.revcomp('W'));
ok('revcomp IUPAC R->Y', B.revcomp('R') === 'Y', B.revcomp('R'));
ok('revcomp IUPAC B->V', B.revcomp('B') === 'V', B.revcomp('B'));
ok('double revcomp is identity',
   B.revcomp(B.revcomp('ATGCGTACGTTAGCNRYKMWS')) === 'ATGCGTACGTTAGCNRYKMWS');

console.log('\n== GC content ==');
ok('GC of GCGC = 100', B.gcContent('GCGC') === 100);
ok('GC of ATAT = 0', B.gcContent('ATAT') === 0);
ok('GC of ATGC = 50', B.gcContent('ATGC') === 50);
ok('GC of 60% case', near(B.gcContent('GGGCCATTAT'), 50, 0.001), B.gcContent('GGGCCATTAT'));

console.log('\n== Tm: internal invariants ==');
const p1 = 'GTAAAACGACGGCCAGT';          // M13 -20 forward
const p2 = 'CAGGAAACAGCTATGAC';          // M13 reverse
const p3 = 'TAATACGACTCACTATAGGG';       // T7 promoter primer
const p4 = 'ATCGATCGATCGATCGATCG';
ok('Tm(seq) == Tm(revcomp) [NN symmetry]',
   near(B.tmNN(p1), B.tmNN(B.revcomp(p1)), 1e-9),
   B.tmNN(p1) + ' vs ' + B.tmNN(B.revcomp(p1)));
ok('Tm rises with salt',
   B.tmNN(p1, {Na:0.05}) < B.tmNN(p1, {Na:1.0}),
   B.tmNN(p1,{Na:0.05}) + ' < ' + B.tmNN(p1,{Na:1.0}));
ok('Tm rises with oligo conc',
   B.tmNN(p1, {oligo:0.1e-6}) < B.tmNN(p1, {oligo:1e-6}));
ok('longer primer, higher Tm',
   B.tmNN('ATCGATCGAT') < B.tmNN('ATCGATCGATCGATCGATCG'));
ok('GC-rich > AT-rich at equal length',
   B.tmNN('GCGCGCGCGCGCGCGCGCGC') > B.tmNN('ATATATATATATATATATAT'));
console.log('    M13-F  (17nt, 58.8% GC) Tm =', B.tmNN(p1).toFixed(1));
console.log('    M13-R  (17nt, 47.1% GC) Tm =', B.tmNN(p2).toFixed(1));
console.log('    T7     (20nt, 40.0% GC) Tm =', B.tmNN(p3).toFixed(1));
console.log('    50%GC  (20nt)           Tm =', B.tmNN(p4).toFixed(1));
// published/accepted ranges (NN, 50 mM Na+, 0.25 uM): M13F ~ 50-56, T7 ~ 48-54
ok('M13-F Tm in accepted range (48-58)', B.tmNN(p1) > 48 && B.tmNN(p1) < 58, B.tmNN(p1));
ok('M13-R Tm in accepted range (44-54)', B.tmNN(p2) > 44 && B.tmNN(p2) < 54, B.tmNN(p2));
ok('T7 Tm in accepted range (44-54)', B.tmNN(p3) > 44 && B.tmNN(p3) < 54, B.tmNN(p3));

console.log('\n== Tm: Mg2+ / PCR conditions (Owczarzy 2008) ==');
const pc = B.PCR_CONDITIONS;
console.log('    M13-F @ PCR conditions Tm =', B.tmNN(p1, pc).toFixed(1));
console.log('    T7    @ PCR conditions Tm =', B.tmNN(p3, pc).toFixed(1));
ok('Mg2+ raises Tm above monovalent-only', B.tmNN(p3, pc) > B.tmNN(p3),
   B.tmNN(p3, pc) + ' vs ' + B.tmNN(p3));
ok('M13-F at PCR conditions is 52-62', B.tmNN(p1, pc) > 52 && B.tmNN(p1, pc) < 62, B.tmNN(p1, pc));
ok('T7 at PCR conditions is 48-58', B.tmNN(p3, pc) > 48 && B.tmNN(p3, pc) < 58, B.tmNN(p3, pc));
// Mg2+ must never push Tm above the 1 M Na+ ceiling at these salt levels
ok('Mg-corrected Tm stays below the 1 M Na+ value',
   B.tmNN(p3, pc) < B.nnThermo(p3).tm1M,
   B.tmNN(p3, pc) + ' vs ' + B.nnThermo(p3).tm1M);
ok('Tm is monotonic in Mg across 1-10 mM', (function () {
  let prev = -Infinity;
  for (const m of [1, 1.5, 2, 3, 4, 6, 8, 10]) {
    const t = B.tmNN(p3, { Mg: m });
    if (t <= prev) return false;
    prev = t;
  }
  return true;
})());
ok('dNTP chelation lowers Tm vs free Mg',
   B.tmNN(p3, {Mg:1.5, dNTP:1.0}) < B.tmNN(p3, {Mg:1.5, dNTP:0}),
   B.tmNN(p3,{Mg:1.5,dNTP:1.0}) + ' vs ' + B.tmNN(p3,{Mg:1.5,dNTP:0}));
ok('Mg correction preserves revcomp symmetry',
   near(B.tmNN(p1, pc), B.tmNN(B.revcomp(p1), pc), 1e-9));
ok('more Mg -> higher Tm', B.tmNN(p3,{Mg:1.0}) < B.tmNN(p3,{Mg:4.0}));
// thermodynamics exposed and self-consistent: dG37 = dH - T*dS
const th = B.nnThermo(p3);
console.log('    T7 thermo: dH =', th.dH.toFixed(1), 'kcal/mol, dS =', th.dS.toFixed(1),
            'cal/mol/K, dG37 =', th.dG37.toFixed(2), 'kcal/mol');
ok('T7 dH matches hand calculation (-148.8)', near(th.dH, -148.8, 0.05), th.dH);
ok('T7 dS matches hand calculation (-411.2)', near(th.dS, -411.2, 0.05), th.dS);
ok('dG37 = dH - 310.15*dS/1000', near(th.dG37, th.dH - 310.15*th.dS/1000, 1e-9));
ok('dH is negative (duplex formation is exothermic)', th.dH < 0);
ok('Wallace rule 20mer 50% GC = 60', B.tmWallace(p4) === 60, B.tmWallace(p4));
ok('ambiguous base -> NaN Tm', Number.isNaN(B.tmNN('ATCGNATCG')));

console.log('\n== restriction site finding ==');
const t1 = 'AAAAGAATTCAAAA';   // EcoRI at index 4
let s = B.findSites(t1, B.enzymeByName('EcoRI'), false);
ok('EcoRI found once', s.length === 1, JSON.stringify(s));
ok('EcoRI recognition at index 4', s[0].pos === 4, s[0].pos);
ok('EcoRI top cut at 5 (G^AATTC)', s[0].cutTop === 5, s[0].cutTop);
ok('EcoRI bottom cut at 9', s[0].cutBottom === 9, s[0].cutBottom);
// none
ok('EcoRI absent in poly-A', B.findSites('AAAAAAAAAA', B.enzymeByName('EcoRI'), false).length === 0);
// SmaI blunt
let sm = B.findSites('TTCCCGGGTT', B.enzymeByName('SmaI'), false);
ok('SmaI blunt: top==bottom cut', sm[0].cutTop === sm[0].cutBottom, JSON.stringify(sm));
// PstI 3' overhang
let ps = B.findSites('AACTGCAGAA', B.enzymeByName('PstI'), false);
ok('PstI gives 3-prime overhang (bottom<top)', ps[0].cutBottom < ps[0].cutTop, JSON.stringify(ps));
// NdeI CA^TATG -> 2nt 5' overhang TA
let nd = B.findSites('GGCATATGCC', B.enzymeByName('NdeI'), false);
ok('NdeI 2nt 5-prime overhang', nd[0].cutBottom - nd[0].cutTop === 2, JSON.stringify(nd));
// BsaI type IIS forward: GGTCTCN^NNNN
let bs = B.findSites('AAGGTCTCATCGAAA', B.enzymeByName('BsaI'), false);
ok('BsaI site found', bs.length >= 1, JSON.stringify(bs));
ok('BsaI leaves 4nt 5-prime overhang', bs[0].cutBottom - bs[0].cutTop === 4, JSON.stringify(bs));
ok('BsaI cuts downstream of site', bs[0].cutTop > bs[0].pos + 6 - 1, JSON.stringify(bs));
// BsaI reverse orientation (GAGACC on top)
let bsr = B.findSites('AAAATCGAGAGACCAA', B.enzymeByName('BsaI'), false);
ok('BsaI reverse site found', bsr.length >= 1, JSON.stringify(bsr));
ok('BsaI reverse still 4nt 5-prime overhang',
   bsr[0].cutBottom - bsr[0].cutTop === 4, JSON.stringify(bsr));
ok('BsaI reverse cuts upstream of site', bsr[0].cutTop < bsr[0].pos, JSON.stringify(bsr));
// circular wrap-around site
const wrap = 'ATTCAAAAAAAAAAAAAAAAGA'; // GA + ATTC spans origin -> GAATTC
ok('circular site spanning origin is found',
   B.findSites(wrap, B.enzymeByName('EcoRI'), true).length === 1,
   JSON.stringify(B.findSites(wrap, B.enzymeByName('EcoRI'), true)));
ok('same site NOT found when linear',
   B.findSites(wrap, B.enzymeByName('EcoRI'), false).length === 0);

console.log('\n== digest ==');
// linear, 1 cut -> 2 fragments summing to length
const lin = 'A'.repeat(100) + 'GAATTC' + 'T'.repeat(100);
let d = B.digest(lin, ['EcoRI'], {circular:false});
ok('linear 1 cut -> 2 fragments', d.fragments.length === 2, d.fragments.length);
ok('linear fragment sizes sum to template',
   d.fragments.reduce((a,f)=>a+f.length,0) === lin.length,
   d.fragments.map(f=>f.length).join(','));
console.log('    linear EcoRI fragments:', d.fragments.map(f=>f.length).join(' + '));
// circular, 1 cut -> 1 linear fragment of full length
let dc = B.digest(lin, ['EcoRI'], {circular:true});
ok('circular 1 cut -> 1 fragment', dc.fragments.length === 1, dc.fragments.length);
ok('circular 1 cut fragment = full length', dc.fragments[0].length === lin.length, dc.fragments[0].length);
// circular, 2 cuts -> 2 fragments summing to length
const circ2 = 'GAATTC' + 'A'.repeat(300) + 'AAGCTT' + 'C'.repeat(500);
let dc2 = B.digest(circ2, ['EcoRI','HindIII'], {circular:true});
ok('circular 2 cuts -> 2 fragments', dc2.fragments.length === 2, dc2.fragments.length);
ok('circular 2-cut sizes sum to plasmid length',
   dc2.fragments.reduce((a,f)=>a+f.length,0) === circ2.length,
   dc2.fragments.map(f=>f.length).join(','));
console.log('    circular EcoRI+HindIII fragments:', dc2.fragments.map(f=>f.length).join(' + '), '(total', circ2.length + ')');
// uncut
let du = B.digest('A'.repeat(50), ['EcoRI'], {circular:true});
ok('no sites -> single uncut fragment', du.fragments.length === 1 && du.fragments[0].note === 'uncut');
// overhang description sanity
ok('EcoRI fragment end reported as 5-prime AATT',
   dc2.fragments.some(f => /5' AATT/.test(f.leftEnd) || /5' AATT/.test(f.rightEnd)),
   JSON.stringify(dc2.fragments.map(f=>[f.leftEnd,f.rightEnd])));

console.log('\n== pUC19 reference checks ==');
// pUC19 is 2686 bp. Known single sites incl. EcoRI, SacI, KpnI, BamHI, XbaI,
// SalI, PstI, SphI, HindIII (the MCS), and AatII, NdeI, etc.
// Build a synthetic construct instead and assert geometry, then test real MCS order.
const mcs = 'GAATTCGAGCTCGGTACCCGGGGATCCTCTAGAGTCGACCTGCAGGCATGCAAGCTT'; // pUC19 MCS (EcoRI..HindIII)
const order = ['EcoRI','SacI','KpnI','XmaI','SmaI','BamHI','XbaI','SalI','PstI','SphI','HindIII'];
let positions = order.map(n => {
  const si = B.findSites(mcs, B.enzymeByName(n), false);
  return { n, pos: si.length ? si[0].pos : -1, count: si.length };
});
for (const p of positions) ok('pUC19 MCS contains ' + p.n, p.count >= 1, JSON.stringify(p));
const monotonic = positions.filter(p=>p.pos>=0).every((p,i,arr) => i===0 || arr[i-1].pos <= p.pos);
ok('pUC19 MCS sites appear in canonical 5->3 order', monotonic,
   positions.map(p=>p.n+'@'+p.pos).join(' '));
console.log('    MCS site order:', positions.map(p=>p.n+'@'+p.pos).join(' '));

console.log('\n== ORF finder ==');
const orfSeq = 'TTT' + 'ATG' + 'AAA'.repeat(40) + 'TAA' + 'GGG';
let orfs = B.findORFs(orfSeq, {minAA: 10});
ok('ORF found', orfs.length >= 1, orfs.length);
const fwdOrf = orfs.find(o => o.strand === 1);
ok('ORF starts at index 3', fwdOrf && fwdOrf.start === 3, fwdOrf && fwdOrf.start);
ok('ORF protein starts with M', fwdOrf && fwdOrf.protein[0] === 'M', fwdOrf && fwdOrf.protein.slice(0,5));
ok('ORF length = 41 aa (M + 40K)', fwdOrf && fwdOrf.length === 41, fwdOrf && fwdOrf.length);
ok('ORF protein has no internal stop', fwdOrf && !fwdOrf.protein.includes('*'), fwdOrf && fwdOrf.protein);
ok('reverse-strand ORFs detected on revcomp input',
   B.findORFs(B.revcomp(orfSeq), {minAA:10}).some(o => o.strand === -1));

console.log('\n== pairwise alignment ==');
let a1 = B.align('ACGTACGT', 'ACGTACGT');
ok('identical seqs -> 100% identity', near(a1.identity, 100, 1e-9), a1.identity);
ok('identical seqs -> no gaps', !a1.a.includes('-') && !a1.b.includes('-'));
let a2 = B.align('ACGTACGT', 'ACGTTCGT');
ok('one mismatch -> 87.5% identity', near(a2.identity, 87.5, 0.01), a2.identity);
let a3 = B.align('ACGTACGT', 'ACGTCGT');  // deletion
ok('deletion produces a gap', a3.a.includes('-') || a3.b.includes('-'), a3.a + ' / ' + a3.b);
ok('aligned strings are equal length', a3.a.length === a3.b.length);
let a4 = B.align('TTTTACGTACGTTTTT', 'ACGTACGT', {mode:'local'});
ok('local alignment finds the embedded match', a4.identity === 100, a4.identity + ' ' + a4.a);
// known Smith-Waterman: classic example
let a5 = B.align('GGTTGACTA', 'TGTTACGG', {mode:'local', match:3, mismatch:-3, gap:-2});
ok('SW classic example score = 13', a5.score === 13, a5.score);
console.log('    SW classic:', a5.a, '/', a5.b, 'score', a5.score);

console.log('\n== PCR ==');
const tpl = 'CCCCC' + 'ATGGCGAATTCGGCGTAAACCGGTTACG' + 'TTTTTTTTTT' + 'GGCATTACGGCATTAAGCTT' + 'AAAAA';
const fw = 'ATGGCGAATTCGGCGTAAA';
const rv = B.revcomp('GGCATTACGGCATTAAGCTT');
let prods = B.pcr(tpl, fw, rv, {anchor:15});
ok('PCR yields a product', prods.length >= 1, prods.length);
if (prods.length) {
  const expected = tpl.indexOf('GGCATTACGGCATTAAGCTT') + 20 - tpl.indexOf(fw);
  ok('PCR product length matches template span', prods[0].length === expected,
     prods[0].length + ' vs ' + expected);
  ok('PCR product starts with fwd primer', prods[0].seq.startsWith(fw));
  ok('PCR product ends with revcomp of rev primer', prods[0].seq.endsWith(B.revcomp(rv)));
}
// primer with a 5' tail (e.g. adding a restriction site) must appear in product
let prods2 = B.pcr(tpl, 'GGGGGCTCGAG' + fw.slice(0,18), rv, {anchor:15});
ok('5-prime tail is carried into the amplicon',
   prods2.length >= 1 && prods2[0].seq.startsWith('GGGGGCTCGAG'),
   prods2.length ? prods2[0].seq.slice(0,20) : 'none');
// no product when primers don't match
ok('no product for non-matching primers',
   B.pcr(tpl, 'TTTTTTTTTTTTTTTTTTTT', rv, {anchor:15}).length === 0);

console.log('\n== Gibson primer design ==');
const vecUp = 'A'.repeat(40) + 'CCCGGGTTTAAACCTGCAGG';
const vecDown = 'TTAATTAAGCGGCCGCTTTT' + 'G'.repeat(40);
const ins = 'ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCTAA';
let gp = B.gibsonPrimers(ins, vecUp, vecDown, {overlap:20, anneal:20});
ok('fwd primer = 20nt vector tail + 20nt insert anneal', gp.insertF.length === 40, gp.insertF.length);
ok('fwd primer tail matches vector 3-prime end',
   gp.insertF.seq.startsWith(vecUp.slice(-20)), gp.insertF.seq.slice(0,20));
ok('fwd primer anneals to insert start',
   ins.startsWith(gp.insertF.anneal), gp.insertF.anneal);
ok('rev primer anneal is revcomp of insert end',
   gp.insertR.anneal === B.revcomp(ins.slice(-20)), gp.insertR.anneal);
ok('rev primer tail is revcomp of vector 5-prime end',
   gp.insertR.seq.startsWith(B.revcomp(vecDown.slice(0,20))));
// the designed primers must actually amplify the insert in the PCR simulator
let gibProd = B.pcr(ins, gp.insertF.seq, gp.insertR.seq, {anchor:18});
ok('designed Gibson primers amplify the insert in silico', gibProd.length >= 1, gibProd.length);
if (gibProd.length) {
  const expLen = 20 + ins.length + 20;
  ok('Gibson amplicon = tail + insert + tail', gibProd[0].length === expLen,
     gibProd[0].length + ' vs ' + expLen);
  // and the junctions must reconstruct seamlessly
  const assembled = vecUp + ins + vecDown;
  ok('amplicon overlaps reconstruct the seamless junction',
     assembled.includes(vecUp.slice(-20) + ins.slice(0,20)) &&
     assembled.includes(ins.slice(-20) + vecDown.slice(0,20)));
}

console.log('\n== CRISPR guides ==');
const target = 'GGGCCCAAATTTGGGCCCAAATTTGGGCCCAAATTTGGGCCCAAATTTGGG';
let guides = B.findGuides(target, {pam:'NGG', length:20});
ok('guides found', guides.length > 0, guides.length);
ok('every guide is 20nt', guides.every(g => g.protospacer.length === 20));
ok('every PAM matches NGG', guides.every(g => /.(GG)$/.test(g.pam)), guides.slice(0,3).map(g=>g.pam).join(','));
ok('guides found on both strands',
   guides.some(g=>g.strand===1) && guides.some(g=>g.strand===-1));
// a guide's protospacer+PAM must exist in the sequence (or its revcomp)
const allOk = guides.every(g => {
  const unit = g.protospacer + g.pam;
  return target.includes(unit) || B.revcomp(target).includes(unit);
});
ok('every protospacer+PAM is present in the template', allOk);

console.log('\n== silent mutation / codon tools ==');
ok('Leu has 6 synonymous codons', B.synonymousCodons('L',1).length === 6, B.synonymousCodons('L',1).join(','));
ok('Met has 1 codon', B.synonymousCodons('M',1).length === 1);
ok('Trp has 1 codon', B.synonymousCodons('W',1).length === 1);
ok('Ser has 6 codons', B.synonymousCodons('S',1).length === 6);
ok('Stop has 3 codons (table 1)', B.synonymousCodons('*',1).length === 3);
// codon optimization must preserve the protein exactly
const cds = 'ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACTAA';
const opt = B.codonOptimize(cds, B.ECOLI_USAGE, 1);
ok('codon-optimized CDS has same length', opt.length === cds.length, opt.length);
ok('codon optimization preserves the protein',
   B.translate(opt) === B.translate(cds),
   B.translate(opt) + ' vs ' + B.translate(cds));
ok('codon optimization actually changed codons', opt !== cds);
const optH = B.codonOptimize(cds, B.HUMAN_USAGE, 1);
ok('human optimization also preserves protein', B.translate(optH) === B.translate(cds));

console.log('\n== intron splicing ==');
const genomic = 'ATGGCC' + 'GTAAGTTTTTTTTTTTTTTTTTTTTAG' + 'GATTAA';
let mrna = B.spliceIntrons(genomic, [{start:6, end:6+27}]);
ok('intron removed', mrna === 'ATGGCCGATTAA', mrna);
ok('spliced product translates cleanly', B.translate(mrna) === 'MAD*', B.translate(mrna));

console.log('\n== molecular weight ==');
const mw = B.ssDNAWeight('ACGT');
ok('ssDNA MW of ACGT is plausible (1150-1300)', mw > 1150 && mw < 1300, mw);
ok('dsDNA MW ~ 2x ssDNA for palindrome',
   near(B.dsDNAWeight('ACGT'), 2*B.ssDNAWeight('ACGT'), 1), B.dsDNAWeight('ACGT'));
ok('MW scales with length', B.ssDNAWeight('A'.repeat(100)) > B.ssDNAWeight('A'.repeat(50)));

console.log('\n== parsers ==');
const fasta = '>pTest some plasmid\nATGCATGCATGC\nGGGCCC\n>second\nAAAATTTT\n';
let fr = B.parseFASTA(fasta);
ok('FASTA: 2 records', fr.length === 2, fr.length);
ok('FASTA: name parsed', fr[0].name === 'pTest', fr[0].name);
ok('FASTA: sequence concatenated across lines', fr[0].seq === 'ATGCATGCATGCGGGCCC', fr[0].seq);
ok('FASTA: second record', fr[1].seq === 'AAAATTTT', fr[1].seq);

const gb = [
'LOCUS       pTestPlasmid            60 bp    ds-DNA     circular SYN 06-OCT-2026',
'DEFINITION  A small test plasmid.',
'FEATURES             Location/Qualifiers',
'     promoter        1..18',
'                     /label="T7 promoter"',
'     CDS             complement(25..45)',
'                     /label="test CDS"',
'ORIGIN',
'        1 taatacgact cactataggg aattcgatat ccatcatcat catcatcatg gatcctaaaa',
'//'
].join('\n');
let gr = B.parseGenBank(gb);
ok('GenBank: 1 record', gr.length === 1, gr.length);
ok('GenBank: name', gr[0].name === 'pTestPlasmid', gr[0].name);
ok('GenBank: circular topology', gr[0].topology === 'circular', gr[0].topology);
ok('GenBank: 60 bp sequence recovered', gr[0].seq.length === 60, gr[0].seq.length);
ok('GenBank: sequence content', gr[0].seq.startsWith('TAATACGACTCACTATAGGG'), gr[0].seq.slice(0,20));
ok('GenBank: 2 features', gr[0].features.length === 2, JSON.stringify(gr[0].features));
ok('GenBank: feature label', gr[0].features[0].label === 'T7 promoter', gr[0].features[0].label);
ok('GenBank: feature coords 0-based start', gr[0].features[0].start === 0 && gr[0].features[0].end === 18,
   JSON.stringify(gr[0].features[0]));
ok('GenBank: complement strand parsed', gr[0].features[1].strand === -1, gr[0].features[1].strand);
// round-trip
const rt = B.toGenBank({name:'pRT', seq:gr[0].seq, topology:'circular', features:gr[0].features, description:'x'});
const rt2 = B.parseGenBank(rt);
ok('GenBank round-trip preserves sequence', rt2[0].seq === gr[0].seq, rt2[0].seq.length);
/* The exporter writes the mandatory source feature that GenBank requires and
   that NCBI will reject a file without, so a re-parse comes back with one
   feature more than went in. The application layer drops it again when it
   builds a record (it is boilerplate, not annotation), which is why the map
   does not draw a band round the whole plasmid. Compare the real features. */
const rtReal = rt2[0].features.filter(f => f.type !== 'source');
ok('GenBank round-trip preserves feature count', rtReal.length === gr[0].features.length,
   rt2[0].features.map(f => f.type).join(','));
ok('GenBank round-trip writes the mandatory source feature',
   rt2[0].features.some(f => f.type === 'source' && f.start === 0 && f.end === gr[0].seq.length),
   JSON.stringify(rt2[0].features.find(f => f.type === 'source')));
ok('GenBank round-trip preserves every feature exactly',
   JSON.stringify(rtReal.map(f => [f.type, f.label, f.start, f.end, f.strand])) ===
   JSON.stringify(gr[0].features.map(f => [f.type, f.label, f.start, f.end, f.strand])),
   JSON.stringify(rtReal.map(f => [f.type, f.label, f.start, f.end, f.strand])));
ok('GenBank round-trip preserves topology', rt2[0].topology === 'circular', rt2[0].topology);
const fa = B.toFASTA('x', 'A'.repeat(150));
ok('FASTA export wraps lines', fa.split('\n')[1].length === 70, fa.split('\n')[1].length);
ok('FASTA export round-trips', B.parseFASTA(fa)[0].seq === 'A'.repeat(150));

console.log('\n== feature auto-detection ==');
const plas = 'TAATACGACTCACTATAG' + 'GGGG' + 'CATCATCATCATCATCAT' + 'TTTT' + B.revcomp('GATTACAAGGATGACGACGATAAG');
let det = B.detectFeatures(plas);
ok('T7 promoter detected', det.some(f=>f.name==='T7 promoter'), det.map(f=>f.name).join(','));
ok('His6 detected', det.some(f=>f.name==='His6 tag'));
ok('FLAG detected on reverse strand',
   det.some(f=>f.name==='FLAG tag' && f.strand===-1),
   JSON.stringify(det.filter(f=>f.name==='FLAG tag')));
ok('detected coords are correct for T7',
   det.find(f=>f.name==='T7 promoter').start === 0);

console.log('\n== enzyme survey ==');
let sv = B.surveyEnzymes(mcs, {circular:false});
const unique = sv.filter(e => e.count === 1).map(e=>e.name);
ok('survey finds unique cutters in MCS', unique.includes('EcoRI') && unique.includes('HindIII'),
   unique.slice(0,12).join(','));
ok('survey reports zero for absent enzyme',
   sv.find(e=>e.name==='NotI').count === 0, sv.find(e=>e.name==='NotI').count);

console.log('\n== enzyme database integrity ==');
let dbBad = [];
for (const e of B.ALL_ENZYMES) {
  if (!/^[ACGTURYSWKMBDHVN]+$/.test(e.site)) dbBad.push(e.name + ' bad site');
  if (!e.typeIIS) {
    if (e.top < 0 || e.top > e.site.length) dbBad.push(e.name + ' top cut out of range');
    if (e.bot < 0 || e.bot > e.site.length) dbBad.push(e.name + ' bottom cut out of range');
  }
}
ok('all enzyme sites are valid IUPAC and cuts in range', dbBad.length === 0, dbBad.join('; '));
// each enzyme must actually cut its own recognition site
let noncut = [];
for (const e of B.ALL_ENZYMES) {
  // flanks long enough for the farthest-reaching Type IIS enzyme (MmeI ~26 bp)
  const flank = 'ATCGATCGATCGATCGATCGATCGATCGATCG';
  const probe = flank + concreteSite(e.site) + flank;
  if (B.findSites(probe, e, false).length === 0) noncut.push(e.name);
}
function concreteSite(site) {
  const pick = {R:'A',Y:'C',S:'G',W:'A',K:'G',M:'A',B:'C',D:'A',H:'A',V:'A',N:'A'};
  return site.split('').map(c => 'ACGT'.includes(c) ? c : (pick[c]||'A')).join('');
}
ok('every enzyme finds its own site', noncut.length === 0, noncut.join(','));
console.log('    enzyme DB size:', B.ALL_ENZYMES.length);

console.log('\n==========================================');
console.log('   ' + pass + ' passed, ' + fail + ' failed');
console.log('==========================================\n');
process.exit(fail ? 1 : 0);
