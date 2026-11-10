const { chromium } = require('playwright');
const path = require('path');
const fs = require('fs');
const B = require('../src/engine.js'); require('../src/engine2.js');

// Wrap the artifact file in the same skeleton the platform applies at publish time,
// so what we test is what the viewer renders.
function wrap(body) {
  return '<!doctype html><html><head><meta charset="utf8">' +
    '<meta name="viewport" content="width=device-width,initial-scale=1,viewport-fit=cover">' +
    '<style>:root{color-scheme:light;padding-top:env(safe-area-inset-top,0px);padding-bottom:env(safe-area-inset-bottom,0px)}' +
    'body{margin:0;font:14px system-ui;background:#fafaf9}img{max-width:100%}[hidden]{display:none!important}</style>' +
    '</head><body>' + body + '</body></html>';
}

let pass = 0, fail = 0;
function ok(name, cond, extra) {
  if (cond) { pass++; console.log('  PASS  ' + name); }
  else { fail++; console.log('  FAIL  ' + name + (extra !== undefined ? '   -> ' + extra : '')); }
}

(async () => {
  const src = fs.readFileSync(path.join(__dirname, '..', 'plasmidbench.html'), 'utf8');
  fs.writeFileSync(path.join(__dirname, 'preview.html'), wrap(src));

  const browser = await chromium.launch();
  const page = await browser.newPage({ viewport: { width: 1440, height: 900 } });
  const errors = [], warnings = [];
  /* The page asks Google Fonts for IBM Plex and falls back to a system font
     when it cannot have it. In this sandbox it never can, and that refused
     request is not a fault in the page. */
  page.on('console', m => {
    if (m.type() === 'error' && !/Failed to load resource|ERR_(TUNNEL|NAME|INTERNET|CONNECTION)/.test(m.text())) errors.push(m.text());
    if (m.type() === 'warning') warnings.push(m.text());
  });
  page.on('pageerror', e => errors.push('PAGEERROR: ' + e.message));

  await page.goto('file://' + path.join(__dirname, 'preview.html'));
  await page.waitForTimeout(900);

  console.log('\n== boot ==');
  ok('no console errors on load', errors.length === 0, errors.slice(0, 4).join(' | '));
  ok('app frame rendered', await page.locator('.app').count() === 1);
  ok('title is PlasmidBench', (await page.title()) === 'PlasmidBench', await page.title());
  const name = await page.locator('#recName').inputValue();
  ok('demo construct loaded', /demo/i.test(name), name);
  const meta = await page.locator('#recMeta').innerText();
  ok('header shows bp / topology / GC', /bp/.test(meta) && /circular/.test(meta) && /% GC/.test(meta), meta.replace(/\n/g, ' '));
  const libCount = await page.locator('.lib-item').count();
  ok('library has the demo entry', libCount >= 1, libCount);
  ok('engine exposed on the page', await page.evaluate(() => typeof window.Bio.tmNN === 'function'));

  // sequence length sanity + verify the engine agrees in-browser
  const demoLen = await page.evaluate(() => window.PB.active().seq.length);
  ok('demo construct has a real sequence', demoLen > 2000, demoLen);
  /* Compare the browser against Node rather than against a number typed in
     here: the point is that the bundled engine and the module under test are
     the same code, and a hardcoded figure only goes stale when the default
     reaction conditions are retuned. */
  const inBrowserTm = await page.evaluate(() => window.Bio.tmNN('GTAAAACGACGGCCAGT'));
  const inNodeTm = B.tmNN('GTAAAACGACGGCCAGT');
  ok('Tm computes the same in the browser as in Node',
     Math.abs(inBrowserTm - inNodeTm) < 1e-9, inBrowserTm + ' vs ' + inNodeTm);

  console.log('\n== map view ==');
  ok('circular map drawn', await page.locator('#view-map svg.map').count() === 1);
  const featPaths = await page.locator('#view-map svg.map g.feat').count();
  ok('features drawn on the map', featPaths >= 5, featPaths);
  /* innerText is an HTML notion; the centre label is an SVG <text>, so read
     textContent. */
  ok('centre label present',
     (await page.evaluate(() => { const e = document.querySelector('#view-map .ctrname'); return e ? e.textContent : ''; })).length > 3);
  const enzTicks = await page.locator('#view-map .enztick').count();
  ok('enzyme ticks drawn', enzTicks > 0, enzTicks);
  // no NaN in any drawn geometry
  const svgHtml = await page.locator('#view-map svg.map').innerHTML();
  ok('no NaN in map geometry', !/NaN/.test(svgHtml), (svgHtml.match(/.{0,40}NaN.{0,40}/) || [''])[0]);
  ok('no undefined in map geometry', !/undefined/.test(svgHtml));

  await page.click('[data-mm="linear"]');
  await page.waitForTimeout(250);
  ok('linear map drawn', await page.locator('#view-map svg.map').count() === 1);
  const linHtml = await page.locator('#view-map svg.map').innerHTML();
  ok('no NaN in linear map', !/NaN/.test(linHtml));
  await page.click('[data-mm="circular"]');
  await page.waitForTimeout(200);

  console.log('\n== sequence view ==');
  await page.click('[data-v="seq"]');
  await page.waitForTimeout(400);
  const lines = await page.locator('#seqInner .seqline').count();
  ok('sequence lines rendered', lines > 10, lines);
  ok('line numbers present', (await page.locator('#seqInner .sl-num').first().innerText()).trim() === '1');
  ok('complement strand shown', await page.locator('#seqInner .sl-strand.bottom').count() > 0);
  ok('translation rows shown', await page.locator('#seqInner .sl-aa').count() > 0);
  ok('feature bars shown', await page.locator('#seqInner .sl-featbar').count() > 0);
  // the rendered bases must match the record exactly
  const firstLine = await page.evaluate(() => {
    const el = document.querySelector('#seqInner .sl-strand');
    return { text: el.textContent, start: +el.getAttribute('data-ls'), n: +el.getAttribute('data-ln') };
  });
  const expectFirst = await page.evaluate(n => window.PB.active().seq.substr(0, n), firstLine.n);
  ok('rendered bases match the record', firstLine.text === expectFirst, firstLine.text.slice(0, 20) + ' vs ' + expectFirst.slice(0, 20));

  // drag a selection across the first line
  const box = await page.locator('#seqInner .sl-strand').first().boundingBox();
  await page.mouse.move(box.x + 2, box.y + box.height / 2);
  await page.mouse.down();
  await page.mouse.move(box.x + box.width * 0.5, box.y + box.height / 2, { steps: 6 });
  await page.mouse.up();
  await page.waitForTimeout(250);
  const sel = await page.evaluate(() => window.PB.state.selection);
  ok('dragging creates a selection', sel && sel.end > sel.start, JSON.stringify(sel));
  ok('selection is within the sequence', sel && sel.end <= demoLen, JSON.stringify(sel));
  const status = await page.locator('#seqStatus').innerText();
  ok('status bar reports the selection', /selected/.test(status), status.replace(/\n/g, ' ').slice(0, 90));
  const insp = await page.locator('#inspector').innerText();
  ok('inspector shows selection stats', /Selection/.test(insp) && /GC/.test(insp));

  // find
  await page.fill('#sqFind', 'GAATTC');
  await page.click('#sqFindBtn');
  await page.waitForTimeout(350);
  const sel2 = await page.evaluate(() => window.PB.state.selection);
  const foundSeq = await page.evaluate(s => window.PB.active().seq.slice(s.start, s.end), sel2);
  ok('find locates an EcoRI site', foundSeq === 'GAATTC', foundSeq);
  // find by enzyme name
  await page.fill('#sqFind', 'HindIII');
  await page.click('#sqFindBtn');
  await page.waitForTimeout(350);
  const sel3 = await page.evaluate(() => window.PB.state.selection);
  const foundSeq3 = await page.evaluate(s => window.PB.active().seq.slice(s.start, s.end), sel3);
  ok('find by enzyme name works', foundSeq3 === 'AAGCTT', foundSeq3);

  console.log('\n== other views ==');
  for (const [v, probe] of [
    ['features', 'table.tbl'], ['enzymes', 'table.tbl'], ['orfs', '.panel'],
    ['primers', '.panel'], ['gel', 'svg.gel'], ['verify', '.panel'],
    ['align', '.panel'], ['history', '.panel']
  ]) {
    await page.click('[data-v="' + v + '"]');
    await page.waitForTimeout(320);
    const n = await page.locator('#view-' + v + ' ' + probe).count();
    const errTxt = await page.locator('#view-' + v + ' .note.bad').count();
    ok('view "' + v + '" renders', n > 0 && errTxt === 0, 'probe=' + n + ' errorNotes=' + errTxt);
  }

  console.log('\n== enzymes table ==');
  await page.click('[data-v="enzymes"]');
  await page.waitForTimeout(400);
  const enzRows = await page.locator('#view-enzymes tbody tr').count();
  ok('enzyme table populated', enzRows > 5, enzRows);
  const enzTxt = await page.locator('#view-enzymes').innerText();
  ok('unique cutters counted', /unique cutters/.test(enzTxt));
  ok('fragment sizes shown', /\d{2,}/.test(enzTxt));

  console.log('\n== gel ==');
  await page.click('[data-v="gel"]');
  await page.waitForTimeout(400);
  const bands = await page.locator('#view-gel svg.gel rect').count();
  ok('gel bands drawn', bands > 10, bands);
  const gelHtml = await page.locator('#view-gel svg.gel').innerHTML();
  ok('no NaN in gel geometry', !/NaN/.test(gelHtml));
  const laneCfg = await page.locator('#view-gel .lane-cfg').count();
  ok('lane configuration rows present', laneCfg >= 2, laneCfg);

  console.log('\n== ORFs ==');
  await page.click('[data-v="orfs"]');
  await page.waitForTimeout(400);
  const orfTxt = await page.locator('#view-orfs').innerText();
  ok('ORF finder ran', /reading frames|No ORFs/.test(orfTxt));
  const orfRows = await page.locator('#view-orfs tbody tr').count();
  ok('ORFs found in the demo construct', orfRows > 0, orfRows);

  console.log('\n== primer design ==');
  await page.click('[data-v="primers"]');
  await page.waitForTimeout(350);
  await page.fill('#pdStart', '100');
  await page.fill('#pdEnd', '700');
  await page.click('#pdGo');
  await page.waitForTimeout(500);
  const pTxt = await page.locator('#view-primers').innerText();
  /* innerText gives back what is rendered, and the section headings are
     upper-cased in CSS, so match without regard to case. */
  ok('primer pair designed', /forward primer/i.test(pTxt) && /Tm \(annealing\)/i.test(pTxt),
     pTxt.replace(/\s+/g, ' ').slice(0, 120));
  ok('suggested annealing temperature shown', /Suggested Ta/.test(pTxt));
  ok('product size reported', /Product/.test(pTxt));
  const tmVals = (pTxt.match(/(\d{2}\.\d) °C/g) || []).map(s => parseFloat(s));
  ok('designed Tm values are plausible', tmVals.length > 0 && tmVals.every(t => t > 45 && t < 80), tmVals.join(','));

  console.log('\n== restriction cloning end to end ==');
  // build a vector and insert in the library, then clone them through the UI
  await page.evaluate(() => {
    const B = window.Bio, PB = window.PB;
    function rnd(n, seed) { let x = seed, s = ''; for (let i = 0; i < n; i++) { x = (x * 1103515245 + 12345) & 0x7fffffff; s += 'ACGT'[Math.floor(x / 0x7fffffff * 4)]; } return s; }
    const vec = PB.makeRecord({ name: 'TestVector', seq: 'G'.repeat(2000) + 'GAATTC' + 'TTTT' + 'AAGCTT' + 'C'.repeat(500), topology: 'circular' });
    const ins = PB.makeRecord({ name: 'TestInsert', seq: 'AAA' + 'GAATTC' + rnd(300, 7) + 'AAGCTT' + 'TTT', topology: 'linear' });
    PB.state.records.push(vec, ins);
    PB.renderAll();
  });
  await page.waitForTimeout(300);
  await page.evaluate(() => window.PB.openTools());
  await page.waitForTimeout(250);
  const toolCards = await page.locator('.tool-card').count();
  ok('tools menu lists operations', toolCards >= 20, toolCards);
  // click "Restriction cloning"
  await page.locator('.tool-card', { hasText: 'Restriction cloning' }).first().click();
  await page.waitForTimeout(300);
  ok('restriction cloning dialog opened', await page.locator('#rcGo').count() === 1);
  /* Playwright's selectOption takes a literal label, not a pattern, so pick
     the option by matching its text here and passing its value. */
  const pickByText = async (sel, re) => {
    const v = await page.evaluate(([sel, src]) => {
      const o = Array.from(document.querySelector(sel).options).find(x => new RegExp(src).test(x.textContent));
      return o ? o.value : null;
    }, [sel, re.source]);
    if (v === null) throw new Error('no option matching ' + re + ' in ' + sel);
    await page.selectOption(sel, v);
  };
  await pickByText('#rcVec_src', /TestVector/);
  await pickByText('#rcIns_src', /TestInsert/);
  await page.fill('#rcVecEnz', 'EcoRI HindIII');
  await page.fill('#rcInsEnz', 'EcoRI HindIII');
  await page.click('#rcGo');
  await page.waitForTimeout(500);
  const rcTxt = await page.locator('#toolResult').innerText();
  ok('ligation reported as working', /Ligation works/.test(rcTxt), rcTxt.slice(0, 160).replace(/\n/g, ' '));
  ok('product size is 2812 bp', /2,812 bp/.test(rcTxt), (rcTxt.match(/[\d,]+ bp circular/) || [''])[0]);
  ok('junctions reported', /AATT/.test(rcTxt) && /AGCT/.test(rcTxt));
  await page.locator('[data-rcsave="0"]').first().click();
  await page.waitForTimeout(500);
  const newLen = await page.evaluate(() => window.PB.active().seq.length);
  ok('saved construct is 2812 bp', newLen === 2812, newLen);
  const savedTopo = await page.evaluate(() => window.PB.active().topology);
  ok('saved construct is circular', savedTopo === 'circular', savedTopo);
  const histLen = await page.evaluate(() => window.PB.active().history.length);
  ok('history records the cloning step', histLen >= 1, histLen);

  console.log('\n== history view ==');
  await page.click('[data-v="history"]');
  await page.waitForTimeout(300);
  const hTxt = await page.locator('#view-history').innerText();
  ok('history shows the cloning operation', /Restriction cloning/.test(hTxt), hTxt.slice(0, 120).replace(/\n/g, ' '));

  console.log('\n== verify / .ab1 path ==');
  // feed a read through the public API and confirm the view reports it
  await page.evaluate(() => {
    const rec = window.PB.active();
    const good = rec.seq.slice(100, 700);
    const bad = good.slice(0, 300) + (good[300] === 'T' ? 'C' : 'T') + good.slice(301);
    window.PB.addRead({ name: 'clone7_fwd', seq: bad });
  });
  await page.waitForTimeout(500);
  const vTxt = await page.locator('#view-verify').innerText();
  ok('read listed in verify table', /clone7_fwd/.test(vTxt));
  ok('exactly one difference found', /1 difference/.test(vTxt), (vTxt.match(/\d+ differences?/) || [''])[0]);
  const diffRows = await page.locator('#view-verify tbody tr').count();
  ok('difference table rendered', diffRows >= 2, diffRows);
  ok('alignment block rendered', await page.locator('#view-verify .aln-box').count() >= 1);
  ok('mismatch highlighted in the alignment', await page.locator('#view-verify .aln-box .mm').count() >= 1);

  console.log('\n== align view ==');
  await page.click('[data-v="align"]');
  await page.waitForTimeout(300);
  await page.evaluate(() => {
    document.querySelector('[data-src="a"]').value = '';
    document.querySelector('[data-src="b"]').value = '';
    document.querySelector('[data-txt="a"]').value = 'ACGTACGTACGTACGT';
    document.querySelector('[data-txt="b"]').value = 'ACGTACGTTCGTACGT';
    window.PB.state.alignState.aText = 'ACGTACGTACGTACGT';
    window.PB.state.alignState.bText = 'ACGTACGTTCGTACGT';
    window.PB.state.alignState.aSrc = ''; window.PB.state.alignState.bSrc = '';
  });
  await page.click('#alGo');
  await page.waitForTimeout(400);
  const aTxt = await page.locator('#view-align').innerText();
  ok('pairwise alignment ran', /Identity/.test(aTxt), aTxt.slice(0, 120).replace(/\n/g, ' '));
  ok('identity is 93.8% for one mismatch in 16', /93\.8/.test(aTxt), (aTxt.match(/9\d\.\d%/) || [''])[0]);

  console.log('\n== theme and responsive ==');
  // dark theme: body background must change and text must stay readable
  const lightBg = await page.evaluate(() => getComputedStyle(document.body).backgroundColor);
  await page.evaluate(() => document.documentElement.setAttribute('data-theme', 'dark'));
  await page.waitForTimeout(220);
  const darkBg = await page.evaluate(() => getComputedStyle(document.body).backgroundColor);
  ok('dark theme changes the background', lightBg !== darkBg, lightBg + ' -> ' + darkBg);
  const darkFg = await page.evaluate(() => getComputedStyle(document.body).color);
  ok('body sets an explicit background', darkBg !== 'rgba(0, 0, 0, 0)' && darkBg !== 'transparent', darkBg);
  function lum(rgb) {
    const m = rgb.match(/\d+/g).map(Number);
    return (0.2126 * m[0] + 0.7152 * m[1] + 0.0722 * m[2]) / 255;
  }
  ok('dark theme has light text on a dark ground', lum(darkFg) > lum(darkBg) + 0.4,
     'fg=' + lum(darkFg).toFixed(2) + ' bg=' + lum(darkBg).toFixed(2));
  await page.evaluate(() => document.documentElement.setAttribute('data-theme', 'light'));
  await page.waitForTimeout(200);
  const lFg = await page.evaluate(() => getComputedStyle(document.body).color);
  const lBg = await page.evaluate(() => getComputedStyle(document.body).backgroundColor);
  ok('light theme has dark text on a light ground', lum(lBg) > lum(lFg) + 0.4,
     'fg=' + lum(lFg).toFixed(2) + ' bg=' + lum(lBg).toFixed(2));
  await page.evaluate(() => document.documentElement.removeAttribute('data-theme'));

  // phone width: no horizontal page scroll
  await page.setViewportSize({ width: 390, height: 844 });
  await page.waitForTimeout(500);
  for (const v of ['map', 'seq', 'features', 'enzymes', 'gel', 'primers', 'align']) {
    await page.click('[data-v="' + v + '"]');
    await page.waitForTimeout(260);
    const over = await page.evaluate(() =>
      document.documentElement.scrollWidth - document.documentElement.clientWidth);
    ok('no horizontal page scroll at 390px in "' + v + '"', over <= 1, 'overflow ' + over + 'px');
  }
  const railHidden = await page.evaluate(() => {
    const r = document.querySelector('#rail');
    return r.getBoundingClientRect().left < -50;
  });
  ok('library rail is off-canvas on a phone', railHidden);
  await page.click('#btnRail');
  await page.waitForTimeout(350);
  const railOpen = await page.evaluate(() => document.querySelector('#rail').getBoundingClientRect().left >= -1);
  ok('rail opens from the menu button', railOpen);

  await page.setViewportSize({ width: 1440, height: 900 });
  await page.waitForTimeout(300);

  console.log('\n== screenshots ==');
  await page.click('[data-v="map"]');
  await page.waitForTimeout(500);
  await page.screenshot({ path: path.join(__dirname, 'shot-map.png') });
  await page.click('[data-v="seq"]');
  await page.waitForTimeout(450);
  await page.screenshot({ path: path.join(__dirname, 'shot-seq.png') });
  await page.evaluate(() => document.documentElement.setAttribute('data-theme', 'dark'));
  await page.click('[data-v="gel"]');
  await page.waitForTimeout(500);
  await page.screenshot({ path: path.join(__dirname, 'shot-gel-dark.png') });
  await page.evaluate(() => document.documentElement.removeAttribute('data-theme'));
  await page.click('[data-v="verify"]');
  await page.waitForTimeout(450);
  await page.screenshot({ path: path.join(__dirname, 'shot-verify.png') });
  console.log('  wrote 4 screenshots');

  console.log('\n== console cleanliness ==');
  ok('no console errors across the whole run', errors.length === 0, errors.slice(0, 5).join(' | '));

  await browser.close();
  console.log('\n==========================================');
  console.log('   ' + pass + ' passed, ' + fail + ' failed');
  console.log('==========================================\n');
  process.exit(fail ? 1 : 0);
})().catch(e => { console.error('HARNESS ERROR', e); process.exit(2); });
