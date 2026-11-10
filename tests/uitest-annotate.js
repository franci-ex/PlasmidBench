/* ===========================================================================
   The annotate panel.

   The engine behind it is tested in test19.js and against Biopython's
   aligner in cmpanno.py. What is tested here is what a person sees and what
   lands in their construct, and the thing that matters most is honesty:

     - a variant is shown as a variant, with its identity, not as a match
     - a fragment is shown as a fragment, with how much of the entry is there
     - nothing is written into the construct until it is chosen
     - what the match was is kept on the feature afterwards, so a call made
       at 91% identity still says so a year later
     - a library of eight-base "features" is called out as such, because
       those are found in everything

   The warning about the built-in library being small is also checked. It is
   27 entries, it will miss most of what is in a real plasmid, and a panel
   that did not say so would be quietly misleading.
   ========================================================================= */
const { chromium } = require('playwright');
const path = require('path');
let pass = 0, fail = 0;
function ok(name, cond, extra) {
  if (cond) { pass++; console.log('  PASS  ' + name); }
  else { fail++; console.log('  FAIL  ' + name + (extra !== undefined ? '   -> ' + extra : '')); }
}

const PAGE = 'file://' + path.resolve(__dirname, '..', 'plasmidbench.html');

(async () => {
  const browser = await chromium.launch();
  const page = await browser.newPage({ viewport: { width: 1500, height: 1150 } });
  const errors = [];
  page.on('pageerror', e => errors.push('pageerror: ' + String(e)));
  page.on('console', m => {
    if (m.type() === 'error' && !/Failed to load resource|ERR_(TUNNEL|NAME|INTERNET|CONNECTION)/.test(m.text()))
      errors.push(m.text());
  });
  await page.goto(PAGE);
  await page.waitForTimeout(1100);

  /* a construct holding one exact feature, one variant and one fragment */
  await page.evaluate(() => {
    const P = window.PB, B = window.Bio;
    const T7 = B.FEATURE_LIB.find(f => f.name === 'T7 promoter').seq;
    const TERM = B.FEATURE_LIB.find(f => f.name === 'T7 terminator').seq;
    const LAC = B.FEATURE_LIB.find(f => f.name === 'lac promoter').seq;
    /* mulberry32, so the filler cannot accidentally contain library entries */
    let a = 24680;
    const r = () => {
      a = (a + 0x6D2B79F5) >>> 0;
      let t = a;
      t = Math.imul(t ^ (t >>> 15), t | 1);
      t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
    const dna = n => { let s = ''; for (let i = 0; i < n; i++) s += 'ACGT'[Math.floor(r() * 4)]; return s; };
    const variant = TERM.split('');
    variant[9] = variant[9] === 'A' ? 'C' : 'A';
    variant[31] = variant[31] === 'G' ? 'T' : 'G';
    const rec = P.makeRecord({
      name: 'annotation test',
      seq: dna(400) + T7 + dna(300) + variant.join('') + dna(300) + B.revcomp(LAC).slice(0, 22) + dna(400),
      topology: 'circular', features: []
    });
    P.state.records = [rec]; P.state.activeId = rec.id;
    P.lsSet('featureLib', null);
    P.renderAll();
  });
  await page.waitForTimeout(400);

  console.log('\n== the panel says what library it is working from ==');
  await page.evaluate(() => window.PB.toolAnnotate());
  await page.waitForTimeout(400);
  const opened = await page.evaluate(() => {
    const m = document.querySelector('#modalRoot');
    return { text: m ? m.textContent : '', hasK: !!document.querySelector('#anK') };
  });
  ok('it opens', opened.hasK);
  ok('and warns that the built-in library is small',
     /built-in entries are loaded/.test(opened.text) &&
     /miss most of what is in a real plasmid/.test(opened.text), '');
  ok('naming the libraries worth importing instead',
     /ApE/.test(opened.text) && /FPbase/.test(opened.text), '');
  ok('and it explains what approximate search adds over exact',
     /variants and fragments/.test(opened.text), '');

  console.log('\n== what it finds, and how it describes each one ==');
  await page.evaluate(() => document.querySelector('#anGo').click());
  await page.waitForTimeout(1200);
  const rows = await page.evaluate(() =>
    Array.from(document.querySelectorAll('#anOut tbody tr')).map(tr => ({
      name: tr.children[1].textContent.trim(),
      where: tr.children[3].textContent.trim(),
      strand: tr.children[4].textContent.trim(),
      match: tr.children[6].textContent.replace(/\s+/g, ' ').trim(),
      from: tr.children[7].textContent.trim()
    })));
  ok('three features are listed', rows.length === 3, JSON.stringify(rows.map(r => r.name)));
  const t7 = rows.find(r => /T7 promoter/.test(r.name));
  const term = rows.find(r => /terminator/.test(r.name));
  const lac = rows.find(r => /lac/.test(r.name));
  ok('the exact one is labelled exact', t7 && /exact/.test(t7.match), t7 && t7.match);
  ok('and as all of the entry', t7 && /all of it/.test(t7.match), t7 && t7.match);
  ok('the variant is not labelled a match but given its identity',
     term && /identical/.test(term.match) && !/exact/.test(term.match), term && term.match);
  ok('to a tenth of a per cent, which is 95.8 for two changes in 48 bases',
     term && /95\.8% identical/.test(term.match), term && term.match);
  ok('the fragment says how much of the entry is there',
     lac && /% of it/.test(lac.match) && !/all of it/.test(lac.match), lac && lac.match);
  ok('and the reverse-strand one says so', lac && /−/.test(lac.strand), lac && lac.strand);
  ok('every row says where the entry came from',
     rows.every(r => r.from.length > 2), JSON.stringify(rows.map(r => r.from)));
  const summary = await page.evaluate(() => document.querySelector('#anOut .note').textContent);
  ok('the summary counts exact against approximate',
     /1 exact, 2 approximate/.test(summary), summary);
  ok('and how long it took', /in \d+ ms/.test(summary), summary);

  console.log('\n== nothing is written into the construct until it is chosen ==');
  const before = await page.evaluate(() => window.PB.active().features.length);
  ok('searching alone adds no features', before === 0, before);
  await page.evaluate(() => {
    /* untick the fragment, keep the other two */
    const boxes = Array.from(document.querySelectorAll('[data-anno]'));
    const rows2 = Array.from(document.querySelectorAll('#anOut tbody tr'));
    rows2.forEach((tr, i) => {
      if (/lac/.test(tr.children[1].textContent)) boxes[i].checked = false;
    });
    document.querySelector('#anAdd').click();
  });
  await page.waitForTimeout(500);
  const added = await page.evaluate(() => {
    const f = window.PB.active().features;
    return f.map(x => ({ label: x.label, note: (x.quals || {}).note, kept: x.annotated,
                         start: x.start, end: x.end, strand: x.strand }));
  });
  ok('only the ticked ones are added', added.length === 2,
     JSON.stringify(added.map(a => a.label)));
  ok('and the unticked one is not', !added.some(a => /lac/.test(a.label)), '');

  console.log('\n== the evidence stays on the feature ==');
  const v = added.find(a => /terminator/.test(a.label));
  const e = added.find(a => /promoter/.test(a.label));
  ok('a variant feature records the identity it was called at',
     v && /95\.8% identical/.test(v.note), v && v.note);
  ok('and over how much of the library entry', v && /100% of/.test(v.note), v && v.note);
  ok('in machine-readable form as well, for anything that reads the record later',
     v && v.kept && Math.abs(v.kept.identity - 46 / 48) < 1e-9 && v.kept.exact === false,
     v && JSON.stringify(v.kept));
  ok('an exact feature says it was exact', e && /exact match/.test(e.note), e && e.note);
  ok('the features are placed where they were found',
     e && e.start === 400 && e.end === 418, e && e.start + '-' + e.end);
  ok('and the step can be undone like any other edit',
     await page.evaluate(() => {
       window.PB.undo();
       return window.PB.active().features.length;
     }) === 0);

  console.log('\n== importing a library, with its provenance ==');
  const imported = await page.evaluate(async () => {
    const P = window.PB, B = window.Bio;
    /* an ApE-shaped file: name, sequence, type */
    const text = [
      'CMV promoter\tCGCAAATGGGCGGTAGGCGTGTACGGTGGGAGGTCTATATAAGCAGAGCT\tpromoter',
      'SV40 polyA\tAACTTGTTTATTGCAGCTTATAATGGTTACAAATAAAGCAATAGCATCAC\tpolyA_signal',
      'tiny thing\tGGATCCAA\tmisc_feature'
    ].join('\n');
    const r = B.parseFeatureLibrary(text, { source: 'ApE default features', licence: 'see ApE' });
    P.lsSet('featureLib', { source: 'ApE default features', licence: 'see ApE', entries: r.entries });
    return { count: r.count, coincidental: r.coincidental, shortest: r.shortest,
             source: r.entries[0].source, licence: r.entries[0].licence };
  });
  ok('the file reads', imported.count === 3, imported.count);
  ok('each entry keeps the source it came from',
     imported.source === 'ApE default features', imported.source);
  ok('and its licence, which differs between sources', imported.licence === 'see ApE',
     imported.licence);
  ok('an eight-base entry is counted as one that turns up by chance',
     imported.coincidental >= 1 && imported.shortest === 8,
     JSON.stringify([imported.coincidental, imported.shortest]));

  /* and the imported entries are then findable */
  const foundImported = await page.evaluate(async () => {
    const P = window.PB;
    const rec = P.active();
    const cmv = 'CGCAAATGGGCGGTAGGCGTGTACGGTGGGAGGTCTATATAAGCAGAGCT';
    rec.seq = rec.seq.slice(0, 200) + cmv + rec.seq.slice(200);
    P.touch(rec);
    P.closeModal();
    P.toolAnnotate();
    await new Promise(r => setTimeout(r, 300));
    document.querySelector('#anGo').click();
    await new Promise(r => setTimeout(r, 1200));
    const names = Array.from(document.querySelectorAll('#anOut tbody tr'))
      .map(tr => tr.children[1].textContent.trim());
    const froms = Array.from(document.querySelectorAll('#anOut tbody tr'))
      .map(tr => tr.children[7].textContent.trim());
    const warn = (document.querySelector('#modalRoot') || {}).textContent || '';
    return { names: names, froms: froms,
             stillWarns: /built-in entries are loaded/.test(warn) };
  });
  ok('an imported feature is found', foundImported.names.some(n => /CMV/.test(n)),
     JSON.stringify(foundImported.names));
  ok('and the table says which library it came from',
     foundImported.froms.some(f => /ApE/.test(f)), JSON.stringify(foundImported.froms));
  ok('the warning about the small built-in library is gone, since one is loaded',
     !foundImported.stillWarns, '');

  console.log('\n== the search settings mean what they say ==');
  const settings = await page.evaluate(async () => {
    const set = (id, v) => { document.querySelector(id).value = String(v); };
    set('#anId', 0.99);
    document.querySelector('#anGo').click();
    await new Promise(r => setTimeout(r, 1000));
    const strict = Array.from(document.querySelectorAll('#anOut tbody tr'))
      .map(tr => tr.children[1].textContent.trim());
    set('#anId', 0.85);
    document.querySelector('#anGo').click();
    await new Promise(r => setTimeout(r, 1000));
    const loose = Array.from(document.querySelectorAll('#anOut tbody tr'))
      .map(tr => tr.children[1].textContent.trim());
    return { strict: strict, loose: loose };
  });
  ok('demanding 99% identity drops the 95.8% variant',
     !settings.strict.some(n => /terminator/.test(n)) &&
     settings.loose.some(n => /terminator/.test(n)),
     JSON.stringify(settings));
  ok('while the exact ones stay', settings.strict.some(n => /T7 promoter/.test(n)),
     JSON.stringify(settings.strict));

  console.log('\n== console ==');
  ok('no page errors anywhere in that', errors.length === 0, errors.slice(0, 3).join(' | '));

  console.log('\n   ' + pass + ' passed, ' + fail + ' failed\n');
  await browser.close();
  process.exit(fail ? 1 : 0);
})();
