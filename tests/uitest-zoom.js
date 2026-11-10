const { chromium } = require('playwright');
const path = require('path'), fs = require('fs');
let pass = 0, fail = 0;
function ok(name, cond, extra) {
  if (cond) { pass++; console.log('  PASS  ' + name); }
  else { fail++; console.log('  FAIL  ' + name + (extra !== undefined ? '   -> ' + extra : '')); }
}
function wrap(body) {
  return '<!doctype html><html><head><meta charset="utf8">' +
    '<meta name="viewport" content="width=device-width,initial-scale=1,viewport-fit=cover">' +
    '<style>:root{color-scheme:light;padding-top:env(safe-area-inset-top,0px);padding-bottom:env(safe-area-inset-bottom,0px)}' +
    'body{margin:0;font:14px system-ui;background:#fafaf9}img{max-width:100%}[hidden]{display:none!important}</style>' +
    '</head><body>' + body + '</body></html>';
}
const V = () => 'window.PB.mapView(window.PB.active())';

(async () => {
  fs.writeFileSync(path.join(__dirname, 'preview.html'),
    wrap(fs.readFileSync(path.join(__dirname, '..', 'plasmidbench.html'), 'utf8')));
  const browser = await chromium.launch();
  const page = await browser.newPage({ viewport: { width: 1440, height: 940 } });
  const errors = [];
  const ignorable = t => /ERR_TUNNEL_CONNECTION_FAILED|ERR_NAME_NOT_RESOLVED|fonts\.(googleapis|gstatic)\.com|Failed to load resource/.test(t);
  page.on('console', m => { if (m.type() === 'error' && !ignorable(m.text())) errors.push(m.text()); });
  page.on('pageerror', e => errors.push('PAGEERROR: ' + e.message));
  await page.goto('file://' + path.join(__dirname, 'preview.html'));
  await page.waitForTimeout(800);
  await page.evaluate(() => {
    const d = window.PB.state.records.find(r => r.builtin);
    window.PB.state.activeId = d.id;
    window.PB.state.view = 'map';
    window.PB.renderAll();
  });
  await page.waitForTimeout(600);

  const demoLen = await page.evaluate(() => window.PB.active().seq.length);

  console.log('\n== linear view: semantic zoom ==');
  await page.click('[data-mm="linear"]');
  await page.waitForTimeout(500);
  let v = await page.evaluate(() => JSON.parse(JSON.stringify(window.PB.mapView(window.PB.active()))));
  ok('opens showing the whole construct', Math.round(v.span) === demoLen && v.start === 0,
     JSON.stringify(v));
  ok('zoom reads 100% when fully zoomed out',
     (await page.locator('#mapZoomPct').textContent()).trim() === '100%',
     await page.locator('#mapZoomPct').textContent());

  const box = await page.locator('#mapWrap svg.map').boundingBox();
  // measure glyph size before zooming, to prove it does not scale
  const fontBefore = await page.evaluate(() => {
    const t = document.querySelector('#mapWrap svg.map text.ticklbl');
    return t ? getComputedStyle(t).fontSize : null;
  });

  await page.mouse.move(box.x + box.width / 2, box.y + box.height * 0.7);
  for (let i = 0; i < 6; i++) { await page.mouse.wheel(0, -400); await page.waitForTimeout(90); }
  await page.waitForTimeout(400);
  v = await page.evaluate(() => JSON.parse(JSON.stringify(window.PB.mapView(window.PB.active()))));
  ok('scrolling shows fewer bases', v.span < demoLen * 0.5, v.span);
  ok('the window stays inside the construct', v.start >= 0 && v.start + v.span <= demoLen + 1,
     JSON.stringify(v));
  const fontAfter = await page.evaluate(() => {
    const t = document.querySelector('#mapWrap svg.map text.ticklbl');
    return t ? getComputedStyle(t).fontSize : null;
  });
  ok('ruler text keeps a constant size while zooming (semantic, not geometric)',
     fontBefore && fontAfter && fontBefore === fontAfter, fontBefore + ' -> ' + fontAfter);
  ok('no scale() transform is used in the linear view',
     !(await page.locator('#mapWrap svg.map').innerHTML()).includes('scale('));

  // keep zooming until the bases appear
  for (let i = 0; i < 26 && v.span > 120; i++) {
    await page.mouse.wheel(0, -400);
    await page.waitForTimeout(60);
    v = await page.evaluate(() => JSON.parse(JSON.stringify(window.PB.mapView(window.PB.active()))));
  }
  await page.waitForTimeout(400);
  const seqRows = await page.locator('#mapWrap svg.map text.seqrow').count();
  ok('the sequence appears at high zoom', seqRows >= 1, 'span=' + Math.round(v.span) + ' rows=' + seqRows);
  const bases = await page.locator('#mapWrap svg.map text.seqrow tspan').count();
  ok('individual bases are drawn', bases > 10, bases);
  // the drawn bases must be the real ones at that coordinate
  const check = await page.evaluate(() => {
    const row = document.querySelector('#mapWrap svg.map text.seqrow');
    const letters = Array.from(row.querySelectorAll('tspan')).map(t => t.textContent).join('');
    const v = window.PB.mapView(window.PB.active());
    const rec = window.PB.active();
    const from = Math.max(0, Math.floor(v.start) - 1);
    return { letters: letters.slice(0, 40), expect: rec.seq.substr(from, letters.length).slice(0, 40) };
  });
  ok('drawn bases match the construct at that position', check.letters === check.expect,
     check.letters + ' vs ' + check.expect);
  const comp = await page.locator('#mapWrap svg.map text.seqrow.dim').count();
  ok('the complement strand is drawn too', comp >= 1, comp);

  console.log('\n== linear view: pan ==');
  const beforeStart = (await page.evaluate(() => window.PB.mapView(window.PB.active()).start));
  const box2 = await page.locator('#mapWrap svg.map').boundingBox();
  await page.mouse.move(box2.x + box2.width / 2, box2.y + box2.height * 0.75);
  await page.mouse.down();
  await page.mouse.move(box2.x + box2.width / 2 - 160, box2.y + box2.height * 0.75, { steps: 8 });
  await page.mouse.up();
  await page.waitForTimeout(350);
  const afterStart = await page.evaluate(() => window.PB.mapView(window.PB.active()).start);
  ok('dragging left moves further along the sequence', afterStart > beforeStart,
     beforeStart + ' -> ' + afterStart);
  const spanAfterPan = await page.evaluate(() => window.PB.mapView(window.PB.active()).span);
  ok('panning does not change the zoom level', Math.abs(spanAfterPan - v.span) < 1e-6);

  console.log('\n== overview strip ==');
  ok('overview strip drawn', await page.locator('#mapWrap svg.map [data-mini]').count() === 1);
  const miniBox = await page.locator('#mapWrap svg.map [data-mini]').boundingBox();
  await page.mouse.click(miniBox.x + miniBox.width * 0.8, miniBox.y + miniBox.height / 2);
  await page.waitForTimeout(350);
  const jumped = await page.evaluate(() => window.PB.mapView(window.PB.active()).start);
  ok('clicking the overview jumps along the construct', jumped > demoLen * 0.5,
     jumped + ' of ' + demoLen);

  console.log('\n== fit selection ==');
  await page.evaluate(() => {
    window.PB.state.selection = { start: 300, end: 420 };
    window.PB.renderHeader(); window.PB.renderMap();
  });
  await page.waitForTimeout(400);
  await page.click('[data-zm="sel"]');
  await page.waitForTimeout(400);
  const fitV = await page.evaluate(() => JSON.parse(JSON.stringify(window.PB.mapView(window.PB.active()))));
  ok('fit-selection frames the selection', fitV.start < 300 && fitV.start + fitV.span > 420,
     JSON.stringify(fitV));
  ok('fit-selection zooms in tightly', fitV.span < 400, fitV.span);

  console.log('\n== zoom limits ==');
  await page.click('[data-zm="fit"]');
  await page.waitForTimeout(250);
  for (let i = 0; i < 10; i++) { await page.click('[data-zm="out"]'); }
  await page.waitForTimeout(350);
  const out = await page.evaluate(() => JSON.parse(JSON.stringify(window.PB.mapView(window.PB.active()))));
  ok('cannot zoom out past the whole construct', Math.round(out.span) <= demoLen && out.start === 0,
     JSON.stringify(out));
  for (let i = 0; i < 30; i++) { await page.click('[data-zm="in"]'); }
  await page.waitForTimeout(400);
  const deep = await page.evaluate(() => JSON.parse(JSON.stringify(window.PB.mapView(window.PB.active()))));
  ok('zooming in stops at the base level, never inverts', deep.span >= 20 && deep.span > 0, deep.span);
  ok('window still inside the construct at max zoom',
     deep.start >= 0 && deep.start + deep.span <= demoLen + 1, JSON.stringify(deep));

  console.log('\n== circular view: labels hold their size ==');
  await page.click('[data-mm="circular"]');
  await page.waitForTimeout(500);
  // effective on-screen font size = declared size x the element's screen CTM.
  // (a bounding box is no good here: a textPath label's box spans its whole arc)
  const effFont = () => page.evaluate(() => {
    const t = document.querySelector('#mapWrap svg.map text.enzlbl');
    if (!t) return 0;
    const fs = parseFloat(getComputedStyle(t).fontSize);
    const m = t.getScreenCTM();
    return m ? fs * m.a : 0;
  });
  const cFontBefore = await effFont();
  const cBox = await page.locator('#mapWrap svg.map').boundingBox();
  await page.mouse.move(cBox.x + cBox.width / 2, cBox.y + cBox.height / 2);
  /* On the circle a plain wheel turns the plasmid; it is a pinch — a wheel
     event with ctrlKey set, which is how a trackpad pinch arrives — that
     zooms. Check both, since a plain wheel that zoomed would make the
     plasmid impossible to turn. */
  const rotBefore = await page.evaluate(() => window.PB.mapView(window.PB.active()).rot || 0);
  for (let i = 0; i < 3; i++) { await page.mouse.wheel(0, -400); await page.waitForTimeout(80); }
  await page.waitForTimeout(350);
  const afterPlain = await page.evaluate(() => {
    const v = window.PB.mapView(window.PB.active());
    return { k: v.k, rot: v.rot || 0 };
  });
  ok('a plain wheel turns the circle rather than zooming it',
     Math.abs(afterPlain.rot - rotBefore) > 1e-6 && Math.abs(afterPlain.k - 1) < 1e-9,
     JSON.stringify(afterPlain));
  /* Each pinch step is capped at 3x, so five steps of -400 land at the 120x
     ceiling, where no enzyme label is drawn at all and there is nothing left
     to measure. Two gentle steps give a few-fold zoom, which is what the
     label-size check is about. */
  await page.keyboard.down('Control');
  for (let i = 0; i < 2; i++) { await page.mouse.wheel(0, -100); await page.waitForTimeout(90); }
  await page.keyboard.up('Control');
  await page.waitForTimeout(450);
  const ck = await page.evaluate(() => window.PB.mapView(window.PB.active()).k);
  ok('circular view zooms in on a pinch', ck > 1.5, ck);
  const cFontAfter = await effFont();
  ok('circular labels keep their on-screen size at ' + ck.toFixed(1) + '× zoom',
     cFontBefore > 0 && cFontAfter > 0 && Math.abs(cFontAfter - cFontBefore) <= 0.6,
     cFontBefore.toFixed(2) + 'px -> ' + cFontAfter.toFixed(2) + 'px');
  const svgHtml = await page.locator('#mapWrap svg.map').innerHTML();
  ok('no NaN in the zoomed circular geometry', !/NaN/.test(svgHtml));

  console.log('\n== state is per construct ==');
  await page.evaluate(() => {
    const r = window.PB.makeRecord({ name: 'ZoomOther', seq: 'ACGT'.repeat(500), topology: 'circular' });
    window.PB.state.records.push(r);
    window.PB.state.activeId = r.id;
    window.PB.renderAll();
  });
  await page.waitForTimeout(500);
  const other = await page.evaluate(() => JSON.parse(JSON.stringify(window.PB.mapView(window.PB.active()))));
  ok('a different construct starts at its own default zoom', other.k === 1 || other.span,
     JSON.stringify(other));
  await page.evaluate(() => {
    const d = window.PB.state.records.find(r => r.builtin);
    window.PB.state.activeId = d.id;
    window.PB.renderAll();
  });
  await page.waitForTimeout(500);
  const back = await page.evaluate(() => window.PB.mapView(window.PB.active()).k);
  ok('returning to a construct restores its zoom', back > 1.5, back);

  console.log('\n== phone width ==');
  await page.evaluate(() => { window.PB.state.selection = { start: 10, end: 2000 }; window.PB.renderHeader(); });
  await page.setViewportSize({ width: 390, height: 844 });
  await page.waitForTimeout(600);
  for (const m of ['circular', 'linear']) {
    await page.click('[data-mm="' + m + '"]');
    await page.waitForTimeout(450);
    const over = await page.evaluate(() =>
      document.documentElement.scrollWidth - document.documentElement.clientWidth);
    ok('no horizontal page scroll in the ' + m + ' map at 390px', over <= 1, over);
  }
  await page.setViewportSize({ width: 1440, height: 940 });
  await page.waitForTimeout(400);

  console.log('\n== screenshots ==');
  await page.click('[data-mm="linear"]');
  await page.waitForTimeout(400);
  await page.evaluate(() => {
    const v = window.PB.mapView(window.PB.active());
    v.start = 120; v.span = 150; window.PB.scheduleMapRender();
  });
  await page.waitForTimeout(500);
  await page.screenshot({ path: path.join(__dirname, 'zoom-bases.png'), clip: { x: 237, y: 86, width: 892, height: 430 } });
  await page.evaluate(() => {
    const v = window.PB.mapView(window.PB.active());
    v.start = 0; v.span = 2312; window.PB.scheduleMapRender();
  });
  await page.waitForTimeout(450);
  await page.screenshot({ path: path.join(__dirname, 'zoom-linear.png'), clip: { x: 237, y: 86, width: 892, height: 380 } });
  await page.click('[data-mm="circular"]');
  await page.waitForTimeout(500);
  await page.evaluate(() => { window.PB.resetMapZoom(); window.PB.scheduleMapRender(); });
  await page.waitForTimeout(400);
  await page.screenshot({ path: path.join(__dirname, 'zoom-circ.png') });

  ok('no console errors throughout', errors.length === 0, errors.slice(0, 4).join(' | '));
  await browser.close();
  console.log('\n==========================================');
  console.log('   ' + pass + ' passed, ' + fail + ' failed');
  console.log('==========================================\n');
  process.exit(fail ? 1 : 0);
})().catch(e => { console.error('HARNESS ERROR', e); process.exit(2); });
