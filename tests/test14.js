/* ===========================================================================
   Expression strategy: tests.

   Two different kinds of claim are made by this part of the tool, and they
   are tested differently.

   The descriptors are arithmetic, and they are checked against Biopython
   separately (cmpprops.py) to the last decimal. What is tested here is
   that they are read correctly: that a membrane protein is called one,
   that a protein with thirty cysteines is sent to an oxidising strain, and
   that the advice follows the numbers rather than being attached to them.

   The compatibility rules are not predictions at all. Two plasmids of one
   origin family will not be co-maintained and one drug will not hold two
   plasmids; those are facts, and the tests state them as such.
   ========================================================================= */
const E = require('../src/engine4.js');

let pass = 0, fail = 0;
function ok(name, cond, extra) {
  if (cond) { pass++; }
  else { fail++; console.log('  FAIL  ' + name + (extra !== undefined ? '   -> ' + extra : '')); }
}
function section(s) { console.log('\n== ' + s + ' =='); }
let seed = 31337;
const rnd = () => { seed = (seed * 1103515245 + 12345) & 0x7fffffff; return seed / 0x7fffffff; };
const AAs = 'ACDEFGHIKLMNPQRSTVWY';
const randProt = n => { let p = 'M'; for (let i = 1; i < n; i++) p += AAs[Math.floor(rnd() * 20)]; return p; };

/* ======================================================================= */
section('protein descriptors');
{
  /* a hand-checkable case: one of each residue */
  const all = AAs;
  const p = E.proteinProps(all);
  ok('every residue is counted', p.length === 20, p.length);
  ok('the weight is the monomers less one water per bond',
     Math.abs(p.mw - (AAs.split('').reduce((t, a) => t + E.AA_MONOMER[a], 0) - 19 * E.WATER)) < 1e-9, p.mw);
  ok('GRAVY is the mean of the hydropathy scale',
     Math.abs(p.gravy - AAs.split('').reduce((t, a) => t + E.KD[a], 0) / 20) < 1e-12, p.gravy);
  ok('the cysteine count is right', p.cys === 1, p.cys);

  const poly = E.proteinProps('D'.repeat(50));
  /* Biopython's Bio.SeqUtils.IsoelectricPoint gives 4.05 for D*50 and for
     E*50, 11.6902 for K*50 and 6.7846 for one of each residue; this engine
     returns the same four figures. The old "< 3.2" here was a guess at what
     a poly-acid should titrate to, not a measured value. */
  ok('an all-aspartate chain titrates where Biopython puts it',
    Math.abs(poly.pI - 4.05) < 0.01, poly.pI);
  const polyK = E.proteinProps('K'.repeat(50));
  ok('an all-lysine chain has a very high pI', polyK.pI > 10.5, polyK.pI);
  ok('and they bracket a mixed sequence', (() => {
    const m = E.proteinProps(randProt(200));
    return m.pI > poly.pI && m.pI < polyK.pI;
  })());

  const hydrophobic = E.proteinProps('V'.repeat(100));
  ok('an all-valine chain is strongly hydrophobic', hydrophobic.gravy > 4, hydrophobic.gravy);
  const charged = E.proteinProps('DEKR'.repeat(25));
  ok('a charged chain is strongly hydrophilic', charged.gravy < -3, charged.gravy);

  ok('an empty input gives nothing back', E.solubilityOutlook('') === null);
}

/* ======================================================================= */
section('reading what the numbers mean');
{
  /* a protein that is mostly membrane helix */
  const tmProt = 'M' + ('VLIAFVLIAFVLIAFVLIAFVLIAF' + 'KEKEKEKEKE').repeat(4);
  const tmOut = E.solubilityOutlook(tmProt);
  ok('several long hydrophobic stretches are found', tmOut.tm.length >= 3, tmOut.tm.length);
  ok('and the verdict says this is the wrong target for soluble work',
     tmOut.band === 'wrong target', tmOut.band);
  ok('the summary says a fusion will not help', /membrane protein/.test(tmOut.summary), tmOut.summary);
  const strat = E.strategyFor(tmOut, []);
  ok('no fusion is offered for it', strat.options.length === 1 && strat.options[0].tag === null,
     JSON.stringify(strat.options.map(o => o.tag)));

  /* a soluble, well-behaved domain */
  const gb1 = 'MTYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE';
  const good = E.solubilityOutlook(gb1);
  ok('a small soluble domain reads as straightforward', good.band === 'straightforward', good.band);
  ok('no membrane helix is claimed in it', good.tm.length === 0, good.tm.length);
  const gstrat = E.strategyFor(good, []);
  ok('a plain tag is the first suggestion for it', gstrat.options[0].tag === 'His',
     gstrat.options[0].tag);

  /* cysteine-rich: should be sent towards an oxidising cytoplasm */
  const cysProt = 'M' + 'CDEKCGSTCNQKRCDEHGSC'.repeat(10);
  const cysOut = E.solubilityOutlook(cysProt);
  ok('a cysteine-rich protein is noticed', cysOut.props.cys >= 6, cysOut.props.cys);
  const cstrat = E.strategyFor(cysOut, []);
  ok('and an oxidising strain is suggested',
     cstrat.conditions.some(c => /oxidising/i.test(c.what)), JSON.stringify(cstrat.conditions.map(c => c.what)));
  ok('with thioredoxin among the fusions',
     cstrat.options.some(o => o.tag === 'Trx'), JSON.stringify(cstrat.options.map(o => o.tag)));

  /* very large: splitting should be suggested */
  const big = E.solubilityOutlook(randProt(1200));
  const bstrat = E.strategyFor(big, []);
  ok('a very large protein is told to be split',
     bstrat.conditions.some(c => /Split it/.test(c.what)), JSON.stringify(bstrat.conditions.map(c => c.what)));

  /* pI sitting on the buffer pH */
  const nearPI = E.solubilityOutlook(gb1, { bufferPH: 4.6 });
  ok('a pI next to the working pH is flagged',
     nearPI.factors.some(f => /Isoelectric/.test(f.name) && f.verdict === 'watch'),
     JSON.stringify(nearPI.factors.filter(f => /Isoelectric/.test(f.name))));
  const farPI = E.solubilityOutlook(gb1, { bufferPH: 8.0 });
  ok('and one well away from it is not',
     farPI.factors.some(f => /Isoelectric/.test(f.name) && f.verdict === 'good'));

  /* DNA in, protein out, with the codon reading that only DNA allows */
  const dna = 'ATG' + 'AGGAGGCGACGGATAGGA'.repeat(12);
  const dOut = E.solubilityOutlook(dna);
  ok('DNA is translated before it is read', dOut.protein[0] === 'M' && dOut.dna, dOut.protein.slice(0, 6));
  ok('and the rare-codon load is reported for it',
     dOut.factors.some(f => /Rare codons/.test(f.name)), JSON.stringify(dOut.factors.map(f => f.name)));
  ok('a protein given as protein gets no codon reading',
     !E.solubilityOutlook(gb1).factors.some(f => /Rare codons/.test(f.name)));
  const rare = E.rareCodonLoad('AGGAGGAGGATGATGATG');
  ok('rare codons are counted and their run measured', rare.count === 3 && rare.worstRun === 3,
     JSON.stringify(rare));
}

/* ======================================================================= */
section('knowing what a vector is');
{
  const cases = [
    ['pET28a-LIC', 'ColE1', 'Kan'], ['pET28-MHL', 'ColE1', 'Kan'],
    ['pET21a', 'ColE1', 'Amp'], ['pET15b', 'ColE1', 'Amp'],
    ['pNIC-CH', 'ColE1', 'Kan'], ['pNIC28-Bsa4', 'ColE1', 'Kan'],
    ['pGEX-6p-1', 'ColE1', 'Amp'], ['pGEX-4T-1', 'ColE1', 'Amp'],
    ['pACYCDuet-1', 'p15A', 'Cm'], ['pRSF-Duet', 'RSF1030', 'Kan'],
    ['pCDFDuet-1', 'CloDF13', 'Spec'], ['pETDuet-1', 'ColE1', 'Amp'],
    ['pGro7 (pACYC)', 'p15A', 'Cm'], ['pRARE2', 'p15A', 'Cm'],
    ['p6HISMAL-MHL', 'ColE1', 'Amp'], ['pET28-MKH8SUMO', 'ColE1', 'Kan'],
    ['pET15NUSA-MHL', 'ColE1', 'Amp'], ['pET15TRX-MHL', 'ColE1', 'Amp'],
    ['pGB1', 'ColE1', 'Kan'], ['pET28NGB-MHL', 'ColE1', 'Kan'],
    ['pQE30', 'ColE1', 'Amp'], ['pTXB1', 'ColE1', 'Amp'],
    ['pFastBac1', 'baculo', 'Amp'], ['pcDNA3.1', 'mammal', 'Amp'],
    ['pDONR223', 'ColE1', 'Kan'], ['pUC57', 'pUChigh', 'Amp']
  ];
  let wrong = 0;
  cases.forEach(([name, ori, marker]) => {
    const k = E.classifyVector(name);
    if (!k) { wrong++; console.log('  not recognised: ' + name); return; }
    if (k.ori !== ori || k.marker !== marker) {
      wrong++; console.log('  ' + name + ' -> ' + k.ori + '/' + k.marker + ', expected ' + ori + '/' + marker);
    }
  });
  ok(cases.length + ' vectors are placed on the right origin and marker', wrong === 0, wrong + ' wrong');
  ok('a name nobody recognises comes back as unknown rather than as a guess',
     E.classifyVector('pSomethingNobodyHasHeardOf123') === null);
  ok('the periplasmic pMAL is told apart from the cytoplasmic one',
     E.classifyVector('pMAL-p2E').periplasm === true && !E.classifyVector('pMAL-c2X').periplasm);
  ok('a fusion vector reports the tag it carries',
     E.classifyVector('p6HISMAL-MHL').tags.indexOf('MBP') >= 0,
     JSON.stringify(E.classifyVector('p6HISMAL-MHL').tags));
  ok('helper plasmids are marked as helpers, not as expression vectors',
     E.classifyVector('pGro7').helper === 'chaperones' && E.classifyVector('pRARE2').helper === 'tRNA');
}

/* ======================================================================= */
section('what cannot share a cell');
{
  const V = n => { const k = E.classifyVector(n); return { name: n, ori: k.ori, marker: k.marker, group: k.group, needsT7: k.needsT7, copy: k.copy }; };

  const duet = E.compatibility([V('pETDuet-1'), V('pACYCDuet-1'), V('pRSF-Duet'), V('pCDFDuet-1')],
    { host: 'BL21(DE3)' });
  ok('the four Duet plasmids are compatible, which is what they are for',
     duet.ok, JSON.stringify(duet.problems.map(p => p.what)));

  const twoColE1 = E.compatibility([V('pET28a-LIC'), V('pGEX-6p-1')], { host: 'BL21(DE3)' });
  ok('two ColE1 backbones are refused even with different markers',
     !twoColE1.ok && twoColE1.problems.some(p => p.kind === 'origin'),
     JSON.stringify(twoColE1.problems.map(p => p.kind)));
  ok('and the reason given is the origin family, not the antibiotic',
     twoColE1.problems.filter(p => p.kind === 'origin').length === 1 &&
     !twoColE1.problems.some(p => p.kind === 'marker'),
     JSON.stringify(twoColE1.problems.map(p => p.kind)));

  const sameDrug = E.compatibility([V('pACYCDuet-1'), V('pRARE2')], { host: 'BL21(DE3)' });
  ok('two chloramphenicol plasmids are refused', !sameDrug.ok &&
     sameDrug.problems.some(p => p.kind === 'marker'), JSON.stringify(sameDrug.problems.map(p => p.kind)));

  /* the trap: a Rosetta host already holds a p15A chloramphenicol plasmid */
  const rosetta = E.compatibility([V('pACYCDuet-1')], { host: 'Rosetta2(DE3)' });
  ok('pACYCDuet in a Rosetta host is caught, though only one plasmid was named',
     !rosetta.ok && rosetta.problems.length >= 2, JSON.stringify(rosetta.problems.map(p => p.kind)));
  ok('because the host carries one of its own',
     rosetta.problems.some(p => /rare-codon/.test(p.what)), JSON.stringify(rosetta.problems.map(p => p.what)));
  const plain = E.compatibility([V('pACYCDuet-1')], { host: 'BL21(DE3)' });
  ok('and in a plain BL21 the same plasmid is fine', plain.ok, JSON.stringify(plain.problems));

  ok('ampicillin is called out for long cold inductions',
     E.compatibility([V('pGEX-6p-1')], { host: 'BL21(DE3)' }).notes.some(n => /beta-lactamase|Carbenicillin/i.test(n)));
  ok('two T7 plasmids are noted as competing for the same polymerase',
     E.compatibility([V('pETDuet-1'), V('pACYCDuet-1')], { host: 'BL21(DE3)' })
       .notes.some(n => /same polymerase/.test(n)));
}

/* ======================================================================= */
section('planning a co-expression');
{
  const inv = ['pET28a-LIC', 'pACYCDuet-1', 'pRSF-Duet', 'pGEX-6p-1'].map(v => ({
    vector: v, kb: E.classifyVector(v), loc: 'L3 / BOX1 / A1'
  }));
  const two = E.coexpressionPlan([{ name: 'A' }, { name: 'B' }], inv, { host: 'BL21(DE3)' });
  ok('two partners are placed', two.assigned.length === 2 && !two.unplaced.length);
  ok('on different origins', two.assigned[0].slot.ori !== two.assigned[1].slot.ori);
  ok('and different markers', two.assigned[0].slot.marker !== two.assigned[1].slot.marker);
  ok('with no conflicts left', two.compatibility.ok, JSON.stringify(two.compatibility.problems));

  const five = E.coexpressionPlan([1, 2, 3, 4, 5].map(i => ({ name: 'p' + i })), inv, { host: 'BL21(DE3)' });
  ok('a fifth partner cannot be placed, because there is no fifth slot',
     five.unplaced.length === 1, five.unplaced.length);
  ok('and the way out is named', five.advice.some(a => /two cassettes|pST39/.test(a.text)),
     JSON.stringify(five.advice.map(a => a.text)));

  const inRosetta = E.coexpressionPlan([{ name: 'A' }, { name: 'B' }, { name: 'C' }, { name: 'D' }],
    inv, { host: 'Rosetta2(DE3)' });
  ok('a host carrying its own plasmid leaves one slot fewer',
     inRosetta.slotsAvailable === 3 && inRosetta.unplaced.length === 1,
     JSON.stringify({ slots: inRosetta.slotsAvailable, unplaced: inRosetta.unplaced.length }));

  const chapInRosetta = E.coexpressionPlan([{ name: 'A' }], inv, { host: 'Rosetta2(DE3)', chaperones: true });
  ok('chaperones in a Rosetta host are refused outright',
     !chapInRosetta.compatibility.ok &&
     chapInRosetta.compatibility.problems.some(p => p.kind === 'helper'),
     JSON.stringify(chapInRosetta.compatibility.problems.map(p => p.kind)));
  const chapInBL21 = E.coexpressionPlan([{ name: 'A' }], inv, { host: 'BL21(DE3)', chaperones: true });
  ok('but in a plain host they are fine', chapInBL21.compatibility.ok,
     JSON.stringify(chapInBL21.compatibility.problems));
  ok('with the separate inducer explained',
     chapInBL21.advice.some(a => /arabinose/.test(a.text)));

  /* the one that needs the most protein goes on the highest-copy backbone */
  const ranked = E.coexpressionPlan([{ name: 'minor', want: 1 }, { name: 'major', want: 5 }],
    inv, { host: 'BL21(DE3)' });
  const major = ranked.assigned.filter(a => a.partner.name === 'major')[0];
  const minor = ranked.assigned.filter(a => a.partner.name === 'minor')[0];
  ok('the partner you need more of gets the higher copy number',
     major.slot.copy > minor.slot.copy, major.slot.copy + ' vs ' + minor.slot.copy);
}

/* ======================================================================= */
section('reading a plasmid inventory');
{
  const tsv = [
    'Plasmid ID\tConstruct ID\tVector\tResistance\tTag\tShelf\tBox\tPosition\tExpression? (Y/N)',
    'YTC-0001\tUbV RFWD3\tpET53\tAmp\tN-6his\tL3\tWYTC-001\tA3\tHigh, soluble, LOBSTR',
    'YTC-0005\tRFWD3\tpET28-MHL\tKan\tN-6his-TEV\tL3\tWYTC-001\tA7\tLow/no, LOBSTR',
    'YTC-0009\tsomething\tpACYCDuet-1\tChl\t\tL3\tWYTC-001\tB1\t',
    'YTC-0010\tno vector here\t\tKan\t\tL3\tWYTC-001\tB2\t'
  ].join('\n');
  const r = E.parsePlasmidStock(tsv);
  ok('the rows with a vector are read', r.plasmids.length === 3, r.plasmids.length);
  ok('and the one without is reported, not dropped silently',
     r.skipped.length === 1 && /no vector/.test(r.skipped[0].why), JSON.stringify(r.skipped));
  ok('each is placed on its backbone', r.plasmids[1].kb.marker === 'Kan', JSON.stringify(r.plasmids[1].kb));
  ok('the tube location is assembled', /WYTC-001/.test(r.plasmids[0].loc) && /A3/.test(r.plasmids[0].loc),
     r.plasmids[0].loc);

  const o1 = r.plasmids[0].outcome, o2 = r.plasmids[1].outcome;
  ok('"High, soluble, LOBSTR" is read as soluble, high, in LOBSTR',
     o1.soluble === true && o1.level === 'high' && o1.host === 'LOBSTR', JSON.stringify(o1));
  ok('"Low/no" is read as no expression', o2.level === 'none', JSON.stringify(o2));
  ok('a blank note gives no outcome at all', r.plasmids[2].outcome === null);
  ok('an uncertain note is left as unknown rather than guessed',
     E.expressionOutcome('Clear? insoluble?, BL21').soluble === null,
     JSON.stringify(E.expressionOutcome('Clear? insoluble?, BL21')));
  ok('"18C" is picked up as a temperature',
     E.expressionOutcome('High. 18C, LOBSTR RIL').temperature === 18);

  const sum = E.vectorSummary(r.plasmids);
  ok('the inventory is summarised by vector family', sum.length === 3, sum.length);
  const kan = sum.filter(v => v.name === 'pET (kanamycin)')[0];
  ok('and the lab’s own outcomes are counted against each',
     kan && kan.noExpression === 1 && kan.soluble === 0, JSON.stringify(kan && { s: kan.soluble, n: kan.noExpression }));

  ok('a paste with no header is refused rather than misread',
     E.parsePlasmidStock('just\tsome\trandom\ttext').plasmids.length === 0);
}

/* ======================================================================= */
section('the advice follows the reading');
{
  const inv = ['pET28a-LIC', 'p6HISMAL-MHL', 'pET28-MKH8SUMO', 'pGB1', 'pET15NUSA-MHL'].map(v => ({
    vector: v, kb: E.classifyVector(v), loc: 'L3'
  }));
  /* something plainly hard: hydrophobic and large */
  let hard = 'M';
  for (let i = 0; i < 600; i++) hard += 'VILFAWMVILFAGVILFCY'[Math.floor(rnd() * 19)];
  const ho = E.solubilityOutlook(hard);
  const hs = E.strategyFor(ho, inv);
  ok('a hard protein is offered MBP first', hs.options[0].tag === 'MBP', hs.options[0].tag);
  ok('and MBP is found on the shelf', hs.options[0].have.length > 0, hs.options[0].have.length);
  ok('with a cold induction among the conditions',
     hs.conditions.some(c => /cold/i.test(c.what)), JSON.stringify(hs.conditions.map(c => c.what)));
  ok('and a lighter induction too',
     hs.conditions.some(c => /lightly/i.test(c.what)), JSON.stringify(hs.conditions.map(c => c.what)));

  const gb1 = 'MTYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE';
  const es = E.strategyFor(E.solubilityOutlook(gb1), inv);
  ok('an easy protein is not sent straight to the biggest fusion',
     es.options[0].tag === 'His', es.options[0].tag);
  ok('every option says why it is there',
     es.options.every(o => o.why && o.why.length > 10));
  ok('an empty shelf still gives advice, just with nothing to point at',
     E.strategyFor(E.solubilityOutlook(gb1), []).options.every(o => o.have.length === 0));
}

console.log('\n   ' + pass + ' passed, ' + fail + ' failed\n');
process.exit(fail ? 1 : 0);
