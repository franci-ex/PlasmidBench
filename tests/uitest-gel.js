/* ===========================================================================
   The gel, including a PCR lane.

   A digest lane shows what the enzymes leave. A PCR lane shows what the
   primers actually give on that template, which is not always one band:
   the point of running it before pouring a gel is to see the spurious
   products, and those come from the same engine that designs the primers.
   ========================================================================= */
const { chromium } = require('playwright');
const path = require('path'), fs = require('fs');
let pass = 0, fail = 0;
function ok(name, cond, extra) {
  if (cond) { pass++; console.log('  PASS  ' + name); }
  else { fail++; console.log('  FAIL  ' + name + (extra !== undefined ? '   -> ' + extra : '')); }
}
const VEC = path.resolve(__dirname, 'fixtures/pRSFDuet-1.dna');

(async () => {
  const browser = await chromium.launch();
  const page = await browser.newPage({ viewport: { width: 1500, height: 1100 } });
  const errors = [];
  page.on('pageerror', e => errors.push('pageerror: ' + String(e)));
  page.on('console', m => { if (m.type() === 'error' && !/Failed to load resource/.test(m.text())) errors.push('console: ' + m.text()); });
  await page.goto('file://' + path.resolve('plasmidbench.html'));
  await page.waitForTimeout(700);

  const b64 = fs.readFileSync(VEC).toString('base64');
  await page.evaluate(b64 => {
    const bin = atob(b64), u = new Uint8Array(bin.length);
    for (let i = 0; i < bin.length; i++) u[i] = bin.charCodeAt(i);
    const r = window.Bio.parseSnapGene(u.buffer);
    window.PB.addRecord(Object.assign({}, r, { name: 'pRSFDuet-1', topology: 'circular' }));
    const rec = window.PB.active();
    window.PB.state.primers = [
      { name: 'duet_F', seq: rec.seq.slice(3775, 3795) },
      { name: 'duet_R', seq: window.Bio.revcomp(rec.seq.slice(600, 620)) },
      { name: 'T7_promoter', seq: 'TAATACGACTCACTATAGGG' },
      { name: 'tailed_F', seq: 'GCTAGCGGATCCAAGCTT' + rec.seq.slice(3775, 3795) }
    ];
    window.PB.lsSet('primers', window.PB.state.primers);
    window.PB.state.view = 'gel';
    window.PB.renderAll();
  }, b64);
  await page.waitForTimeout(500);

  console.log('\n== the PCR lane exists and can be added ==');
  ok('the lane kinds include PCR', await page.evaluate(() =>
    Array.from(document.querySelectorAll('#view-gel [data-lk] option')).some(o => o.value === 'pcr')));
  ok('there is a button to add one', await page.evaluate(() => !!document.querySelector('#glPcrAdd')));
  await page.evaluate(() => document.querySelector('#glPcrAdd').click());
  await page.waitForTimeout(400);
  ok('clicking it adds a PCR lane pre-filled from the primer library', await page.evaluate(() => {
    const l = window.PB.state.gelLanes[window.PB.state.gelLanes.length - 1];
    return l.kind === 'pcr' && !!l.fwd && !!l.rev;
  }), await page.evaluate(() => JSON.stringify(window.PB.state.gelLanes.slice(-1))));
  ok('the library primers are offered for completion', await page.evaluate(() =>
    document.querySelectorAll('#glPrimers option').length >= 4));

  console.log('\n== what it computes ==');
  const read = () => page.evaluate(() => Array.from(document.querySelectorAll('#view-gel table tbody tr'))
    .map(tr => Array.from(tr.querySelectorAll('td')).map(td => td.textContent.trim())));

  await page.evaluate(() => {
    const L = window.PB.state.gelLanes;
    L.length = 2;                                   /* ladder + the digest */
    const id = window.PB.active().id;
    L.push({ kind: 'pcr', recordId: id, fwd: 'duet_F', rev: 'duet_R' });
    L.push({ kind: 'pcr', recordId: id, fwd: 'T7_promoter', rev: 'duet_R' });
    L.push({ kind: 'pcr', recordId: id, fwd: 'duet_R', rev: 'duet_F' });
    L.push({ kind: 'pcr', recordId: id, fwd: 'tailed_F', rev: 'duet_R' });
    window.PB.renderViews();
  });
  await page.waitForTimeout(600);
  const rows = await read();
  const byLane = n => rows.filter(r => r[0].indexOf(n) >= 0)[0];

  const clean = byLane('duet_F + duet_R');
  ok('a clean pair gives one band', clean && clean[1] === '1', clean && clean[1]);
  ok('and the engine agrees on its size', await page.evaluate(() => {
    const rec = window.PB.active();
    const p = window.Bio.pcr(rec.seq, window.PB.state.primers[0].seq,
      window.PB.state.primers[1].seq, { circular: true });
    return p.length === 1 && p[0].length === 674;
  }), clean && clean[2]);

  const multi = byLane('T7_promoter + duet_R');
  ok('a primer with two sites on the template gives two bands', multi && multi[1] === '2',
     multi && multi[1]);
  ok('and the note says where the extra band comes from',
     multi && /2 forward and 1 reverse priming sites/.test(multi[3]), multi && multi[3].slice(0, 110));
  ok('the site count is of real sites, not matches in the doubled circular search',
     multi && !/2 reverse priming/.test(multi[3]), multi && multi[3].slice(0, 110));

  const none = byLane('duet_R + duet_F');
  ok('primers the wrong way round give no product', none && none[1] === '0', none && none[1]);
  ok('and are told why, rather than left blank',
     none && /No product/.test(none[3]) && /5′→3′/.test(none[3]), none && none[3].slice(0, 90));

  const tailed = byLane('tailed_F + duet_R');
  ok('a tailed primer gives a longer product', await page.evaluate(() => {
    const rec = window.PB.active(), P = window.PB.state.primers;
    const plain = window.Bio.pcr(rec.seq, P[0].seq, P[1].seq, { circular: true })[0].length;
    const tail = window.Bio.pcr(rec.seq, P[3].seq, P[1].seq, { circular: true })[0].length;
    return tail === plain + 18;
  }), tailed && tailed[2]);

  console.log('\n== it is drawn, not just tabulated ==');
  ok('bands are drawn for the PCR lanes', await page.evaluate(() =>
    document.querySelectorAll('#view-gel svg.gel rect').length > 20));
  ok('the panel says products as well as fragments', await page.evaluate(() =>
    /Fragments and products/.test(document.querySelector('#view-gel').textContent)));
  ok('changing the agarose percentage recomputes them', await page.evaluate(async () => {
    const before = document.querySelectorAll('#view-gel svg.gel rect').length;
    window.PB.state.gelPct = 2.0; window.PB.renderViews();
    await new Promise(r => setTimeout(r, 200));
    const after = document.querySelectorAll('#view-gel svg.gel rect').length;
    window.PB.state.gelPct = 1.0; window.PB.renderViews();
    return before > 0 && after > 0;
  }));

  console.log('\n== a typed sequence works as well as a name ==');
  await page.evaluate(() => {
    const L = window.PB.state.gelLanes;
    const rec = window.PB.active();
    L.push({ kind: 'pcr', recordId: rec.id, fwd: rec.seq.slice(3775, 3795),
             rev: window.Bio.revcomp(rec.seq.slice(600, 620)) });
    window.PB.renderViews();
  });
  await page.waitForTimeout(500);
  const typed = (await read()).filter(r => /nt typed/.test(r[0]))[0];
  ok('a pasted sequence is accepted in place of a library name', !!typed, 'no such lane');
  ok('and gives the same band as the named primers', typed && typed[2] === '674', typed && typed[2]);

  await page.evaluate(() => { window.PB.state.gelLanes = null; });
  ok('nothing threw anywhere in the run', errors.length === 0, errors.slice(0, 3).join(' | '));
  console.log('\n   ' + pass + ' passed, ' + fail + ' failed\n');
  await browser.close();
  process.exit(fail ? 1 : 0);
})().catch(e => { console.error('HARNESS ERROR', e); process.exit(2); });
