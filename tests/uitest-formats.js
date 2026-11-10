/* ===========================================================================
   Every format, through the real import path, in a real browser.

   test15.js checks the readers; this checks what happens when a file of each
   kind is dropped on the page: that it lands in the library, that the map
   draws it without throwing, that the wrong kinds are refused with something
   worth reading, and that the monochrome export really has no colour left.

   The pGro7 file is the one that started this: a SnapGene-exported GenBank
   plasmid with a mandatory `source` feature spanning all 5,467 bases and a
   CmR gene written as join(5249..5467,1..438), which wraps the origin.
   ========================================================================= */
const { chromium } = require('playwright');
const path = require('path'), fs = require('fs');
let pass = 0, fail = 0;
function ok(name, cond, extra) {
  if (cond) { pass++; console.log('  PASS  ' + name); }
  else { fail++; console.log('  FAIL  ' + name + (extra !== undefined ? '   -> ' + extra : '')); }
}
/* pGro7 as SnapGene exports it: a mandatory source feature over all 5,467
   bases, and CmR written as join(5249..5467,1..438) across the origin */
const GB = path.resolve(__dirname, 'fixtures/pGro7.gb');
const SEQ = 'ATGGCTAGCAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGGCACAAATTTTCTGTC';

(async () => {
  const browser = await chromium.launch();
  const page = await browser.newPage({ viewport: { width: 1500, height: 1050 } });
  const errors = [];
  page.on('pageerror', e => errors.push('pageerror: ' + String(e)));
  page.on('console', m => { if (m.type() === 'error' && !/Failed to load resource/.test(m.text())) errors.push('console: ' + m.text()); });
  await page.goto('file://' + path.resolve('plasmidbench.html'));
  await page.waitForTimeout(700);

  const imp = (text, name) => page.evaluate(([t, n]) => window.PB.importText(t, n), [text, name]);
  const active = () => page.evaluate(() => {
    const r = window.PB.active();
    return { name: r.name, len: r.seq.length, topo: r.topology, nf: r.features.length,
      types: r.features.map(f => f.type),
      feats: r.features.map(f => ({ l: f.label, s: f.start, e: f.end, st: f.strand, c: f.color || null })) };
  });
  const toastText = () => page.evaluate(() => (document.querySelector('#toast') || {}).textContent || '');

  /* ------------------------------------------------- the file that failed */
  console.log('\n== the GenBank plasmid that started this ==');
  const gbText = fs.readFileSync(GB, 'utf8');
  ok('it imports', (await imp(gbText, 'pGro7.gb')) === 1);
  let a = await active();
  ok('with its full length', a.len === 5467, a.len);
  ok('as a circle', a.topo === 'circular', a.topo);
  ok('the mandatory `source` feature is not drawn as a band round the whole plasmid',
    a.types.indexOf('source') === -1, JSON.stringify(a.types));
  ok('so the count is of real features', a.nf === 5, a.nf);
  const cmr = a.feats.filter(f => f.l === 'CmR')[0];
  ok('the gene written as join(5249..5467,1..438) keeps its wrap',
    cmr && cmr.s === 5248 && cmr.e === 438, JSON.stringify(cmr));
  ok('and its bases are the two pieces in order, 657 of them', await page.evaluate(() => {
    const r = window.PB.active(), f = r.features.filter(x => x.label === 'CmR')[0];
    const s = window.Bio.featureSequence(r.seq, f);
    return s.length === 657 && s === r.seq.slice(5248) + r.seq.slice(0, 438);
  }));
  ok('the toast names the format it was read as', /GenBank/.test(await toastText()), await toastText());
  ok('the map draws it', await page.evaluate(async () => {
    window.PB.state.view = 'map'; window.PB.renderAll();
    await new Promise(r => setTimeout(r, 400));
    return document.querySelectorAll('#view-map svg.map .feat').length === 5;
  }), await page.evaluate(() => document.querySelectorAll('#view-map svg.map .feat').length));
  ok('the feature key has no "Other" row for it', await page.evaluate(() =>
    !/\bOther\b/.test((document.querySelector('#inspector') || {}).textContent || '')));

  /* ------------------------------------------------------- every format */
  console.log('\n== every format on the list ==');
  const cases = [
    ['GenBank .gbk',
      'LOCUS       pGbk                      90 bp    DNA     circular SYN 01-JAN-2024\n' +
      'FEATURES             Location/Qualifiers\n     source          1..90\n                     /organism="synthetic"\n' +
      '     CDS             4..30\n                     /label=gfp\n' +
      'ORIGIN\n        1 ' + SEQ.toLowerCase() + '\n//\n', 'p.gbk', { len: 90, nf: 1, topo: 'circular' }],
    ['ApE .ape (GenBank inside)',
      'LOCUS       pApe        90 bp ds-DNA     circular     01-JAN-2024\n' +
      'FEATURES             Location/Qualifiers\n     misc_feature    4..30\n' +
      '                     /label=ape\n                     /ApEinfo_fwdcolor=#ff6600\n' +
      'ORIGIN\n        1 ' + SEQ.toLowerCase() + '\n//\n', 'p.ape', { len: 90, nf: 1, color: '#ff6600' }],
    ['FASTA .fna', '>pFna a plasmid\n' + SEQ + '\n', 'p.fna', { len: 90 }],
    ['EMBL .embl',
      'ID   pEmbl; SV 1; circular; DNA; STD; SYN; 90 BP.\n' +
      'FH   Key             Location/Qualifiers\n' +
      'FT   CDS             4..30\n' +
      'FT                   /label="emb"\n' +
      'SQ   Sequence 90 BP;\n     ' + SEQ.toLowerCase().slice(0, 60) + ' 60\n     ' +
      SEQ.toLowerCase().slice(60) + ' 90\n//\n', 'p.embl', { len: 90, nf: 1, topo: 'circular' }],
    ['GFF3 with bases .gff3',
      '##gff-version 3\n##sequence-region pGff 1 90\n' +
      'pGff\tt\tCDS\t4\t30\t.\t+\t0\tID=a;Name=gff-cds\n##FASTA\n>pGff\n' + SEQ + '\n', 'p.gff3', { len: 90, nf: 1 }],
    ['Benchling JSON',
      JSON.stringify({ name: 'pBench', bases: SEQ, isCircular: true,
        annotations: [{ name: 'bench', type: 'CDS', start: 3, end: 30, strand: 1 }] }), 'p.json',
      { len: 90, nf: 1, topo: 'circular' }],
    ['openVectorEditor JSON',
      JSON.stringify({ name: 'pOve', sequence: SEQ, circular: false,
        features: [{ name: 'ove', type: 'CDS', start: 3, end: 29, forward: true }] }), 'p.json',
      { len: 90, nf: 1, topo: 'linear' }],
    ['INSDSeq XML',
      '<?xml version="1.0"?><INSDSet><INSDSeq><INSDSeq_locus>pXml</INSDSeq_locus>' +
      '<INSDSeq_topology>circular</INSDSeq_topology><INSDSeq_feature-table><INSDFeature>' +
      '<INSDFeature_key>CDS</INSDFeature_key><INSDFeature_location>4..30</INSDFeature_location>' +
      '<INSDFeature_quals><INSDQualifier><INSDQualifier_name>label</INSDQualifier_name>' +
      '<INSDQualifier_value>xml</INSDQualifier_value></INSDQualifier></INSDFeature_quals>' +
      '</INSDFeature></INSDSeq_feature-table><INSDSeq_sequence>' + SEQ.toLowerCase() +
      '</INSDSeq_sequence></INSDSeq></INSDSet>', 'p.xml', { len: 90, nf: 1, topo: 'circular' }],
    ['TinySeq XML',
      '<?xml version="1.0"?><TSeqSet><TSeq><TSeq_accver>NM_1.1</TSeq_accver>' +
      '<TSeq_sequence>' + SEQ + '</TSeq_sequence></TSeq></TSeqSet>', 'p.xml', { len: 90 }],
    ['bare bases pasted', SEQ, 'pasted', { len: 90 }]
  ];
  for (const [label, text, name, want] of cases) {
    const n = await imp(text, name);
    const r = await active();
    const good = n === 1 && r.len === want.len &&
      (want.nf === undefined || r.nf === want.nf) &&
      (want.topo === undefined || r.topo === want.topo) &&
      (want.color === undefined || r.feats[0].c === want.color);
    ok(label + ' imports and lands right', good, JSON.stringify({ n: n, got: r }));
    /* every one of them has to draw without throwing */
    await page.evaluate(() => { window.PB.state.view = 'map'; window.PB.renderAll(); });
    await page.waitForTimeout(120);
  }
  ok('and the bases are identical whichever door they came through', await page.evaluate(SEQ => {
    const seqs = window.PB.state.records.filter(r => r.seq.length === 90).map(r => r.seq);
    return seqs.length >= 8 && seqs.every(s => s === SEQ);
  }, SEQ), await page.evaluate(() => window.PB.state.records.filter(r => r.seq.length === 90).length));
  ok('and so is the feature they carry, to the base', await page.evaluate(SEQ => {
    const want = SEQ.slice(3, 30);
    const got = window.PB.state.records.filter(r => r.seq.length === 90 && r.features.length === 1)
      .map(r => window.Bio.featureSequence(r.seq, r.features[0]));
    return got.length >= 6 && got.every(s => s === want);
  }, SEQ), await page.evaluate(SEQ => JSON.stringify(window.PB.state.records
    .filter(r => r.seq.length === 90 && r.features.length === 1)
    .map(r => r.name + ':' + window.Bio.featureSequence(r.seq, r.features[0]).length)), SEQ));

  /* -------------------------------------------------- GFF without bases */
  console.log('\n== a GFF that is annotations only ==');
  await page.evaluate(SEQ => {
    window.PB.addRecord(window.PB.makeRecord({ name: 'host', seq: SEQ, topology: 'circular', features: [] }));
  }, SEQ);
  const added = await imp('##gff-version 3\nhost\tt\tCDS\t4\t30\t.\t+\t0\tID=a;Name=from-gff\n', 'ann.gff3');
  ok('its annotations are applied to the construct already open', added === 1, added);
  ok('and land on the right bases', await page.evaluate(SEQ => {
    const r = window.PB.active(), f = r.features.filter(x => x.label === 'from-gff')[0];
    return !!f && window.Bio.featureSequence(r.seq, f) === SEQ.slice(3, 30);
  }, SEQ));
  const refused = await imp('##gff-version 3\nhost\tt\tCDS\t400\t900\t.\t+\t0\tID=a;Name=elsewhere\n', 'wrong.gff3');
  ok('annotations that run past the end of it are refused, not clamped to fit', refused === 0, refused);
  ok('and the construct is left alone', await page.evaluate(() =>
    !window.PB.active().features.some(f => f.label === 'elsewhere')));
  ok('with a message that says why', /past the end/.test(await toastText()), await toastText());

  /* ----------------------------------------------------- refusals */
  console.log('\n== what cannot be read says so ==');
  const bad = [
    ['{"colour":"blue"}', 'x.json', /no sequence|No sequence/],
    ['{oops', 'x.json', /not valid JSON|could not be read/],
    ['<?xml version="1.0"?><nothing/>', 'x.xml', /No sequence|INSDSeq/]
  ];
  for (const [text, name, re] of bad) {
    await imp(text, name);
    const t = await toastText();
    ok(name + ' with ' + text.slice(0, 14) + '… is refused with a reason', re.test(t), t);
  }
  ok('a CLC file is named as a private binary format with a way out', await page.evaluate(() =>
    /Export.*GenBank/i.test(window.Bio.UNREADABLE.clc)), await page.evaluate(() => window.Bio.UNREADABLE.clc));
  ok('and an SVG is called a picture rather than a sequence', await page.evaluate(() =>
    /picture/.test(window.Bio.UNREADABLE.svg)));

  /* ----------------------------------------------------- monochrome */
  console.log('\n== the monochrome image export ==');
  await page.evaluate(() => {
    const r = window.PB.state.records.filter(x => x.seq.length === 5467)[0];
    window.PB.state.activeId = r.id; window.PB.state.view = 'map'; window.PB.renderAll();
  });
  await page.waitForTimeout(400);
  const svgs = await page.evaluate(() => ({
    light: window.PB.mapSVGString({ whole: true, theme: 'light', background: 'solid' }),
    mono: window.PB.mapSVGString({ whole: true, theme: 'mono', background: 'solid' })
  }));
  const colours = s => {
    const out = new Set();
    (s.match(/rgba?\([^)]*\)|#[0-9a-f]{3,8}/gi) || []).forEach(c => out.add(c.toLowerCase()));
    return out;
  };
  const isGrey = c => {
    const m = /^rgba?\(\s*(\d+)[,\s]+(\d+)[,\s]+(\d+)/.exec(c);
    if (m) return Math.max(+m[1], +m[2], +m[3]) - Math.min(+m[1], +m[2], +m[3]) <= 8;
    const h = /^#([0-9a-f]{6})$/.exec(c);
    if (h) { const v = [0, 2, 4].map(i => parseInt(h[1].substr(i, 2), 16)); return Math.max(...v) - Math.min(...v) <= 8; }
    const s3 = /^#([0-9a-f]{3})$/.exec(c);
    if (s3) { const v = s3[1].split('').map(x => parseInt(x + x, 16)); return Math.max(...v) - Math.min(...v) <= 8; }
    return true;
  };
  const lightCols = Array.from(colours(svgs.light)), monoCols = Array.from(colours(svgs.mono));
  ok('the light export does have colour in it', lightCols.some(c => !isGrey(c)), lightCols.filter(c => !isGrey(c)).slice(0, 4).join(' '));
  ok('the monochrome export has none left anywhere', monoCols.every(isGrey), monoCols.filter(c => !isGrey(c)).slice(0, 6).join(' '));
  ok('it is still the same map, not an empty one', /Chaperonin GroEL/.test(svgs.mono) && /CmR/.test(svgs.mono));
  ok('its background is white', /fill="#ffffff"/.test(svgs.mono));
  ok('every feature tone is dark enough for the white label on it to read', await page.evaluate(() => {
    const s = window.PB.mapSVGString({ whole: true, theme: 'mono', background: 'solid' });
    const tones = (s.match(/stroke:rgb\((\d+), ?(\d+), ?(\d+)\)/g) || [])
      .map(m => +/rgb\((\d+)/.exec(m)[1])
      .filter(v => v > 0 && v < 250);
    return tones.length > 0 && tones.every(v => v <= 160);
  }));
  ok('an imported SnapGene or ApE colour is converted too, not left in', await page.evaluate(() => {
    const rec = window.PB.state.records.filter(r => r.features.some(f => f.color === '#ff6600'))[0];
    if (!rec) return false;
    window.PB.state.activeId = rec.id; window.PB.renderAll();
    const s = window.PB.mapSVGString({ whole: true, theme: 'mono', background: 'solid' });
    return s.indexOf('255, 102, 0') === -1 && s.indexOf('#ff6600') === -1;
  }));
  ok('the page is left in the theme it was in', await page.evaluate(() =>
    document.documentElement.getAttribute('data-theme') !== 'mono'));
  ok('the dialog offers it', await page.evaluate(async () => {
    window.PB.exportMapDialog();
    await new Promise(r => setTimeout(r, 300));
    const sel = document.querySelector('#meTheme');
    const has = !!sel && Array.from(sel.options).some(o => o.value === 'mono');
    window.PB.closeModal();
    return has;
  }));

  ok('nothing threw anywhere in the run', errors.length === 0, errors.slice(0, 3).join(' | '));
  console.log('\n   ' + pass + ' passed, ' + fail + ' failed\n');
  await browser.close();
  process.exit(fail ? 1 : 0);
})().catch(e => { console.error('HARNESS ERROR', e); process.exit(2); });
