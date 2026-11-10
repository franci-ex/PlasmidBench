/* ===========================================================================
   Selecting more than one feature at a time.

   Ctrl, or cmd on a Mac, adds a feature to the selection instead of
   replacing it, and takes it back out if it is already in. The thing to
   guard is not the clicking but the consequences: that every selected
   feature is painted on the map AND in the sequence, that the hundred
   places which clear `selectedFeature` still clear the whole set, and
   that a feature deleted out from under a selection cannot leave a stale
   index behind for something else to read.
   ========================================================================= */
const { chromium } = require('playwright');
const path = require('path');
let pass = 0, fail = 0;
function ok(name, cond, extra) {
  if (cond) { pass++; console.log('  PASS  ' + name); }
  else { fail++; console.log('  FAIL  ' + name + (extra !== undefined ? '   -> ' + extra : '')); }
}

(async () => {
  const browser = await chromium.launch();
  const page = await browser.newPage({ viewport: { width: 1500, height: 1050 } });
  const errors = [];
  page.on('pageerror', e => errors.push('pageerror: ' + String(e)));
  page.on('console', m => { if (m.type() === 'error' && !/Failed to load resource/.test(m.text())) errors.push('console: ' + m.text()); });
  await page.goto('file://' + path.resolve('plasmidbench.html'));
  await page.waitForTimeout(800);

  const sel = () => page.evaluate(() => window.PB.featureSelection(window.PB.active()));
  const clickFeat = (i, mods) => page.evaluate(([i, mods]) => {
    const g = document.querySelector('#view-map [data-fi="' + i + '"]');
    g.dispatchEvent(new MouseEvent('click', Object.assign({ bubbles: true }, mods)));
  }, [i, mods || {}]);

  await page.evaluate(() => { window.PB.state.view = 'map'; window.PB.renderAll(); });
  await page.waitForTimeout(500);

  console.log('\n== on the map ==');
  await clickFeat(1); await page.waitForTimeout(200);
  ok('a plain click selects one', JSON.stringify(await sel()) === '[1]', JSON.stringify(await sel()));
  await clickFeat(2, { ctrlKey: true }); await page.waitForTimeout(200);
  await clickFeat(3, { ctrlKey: true }); await page.waitForTimeout(200);
  ok('ctrl-click adds to the selection', JSON.stringify(await sel()) === '[1,2,3]', JSON.stringify(await sel()));
  await clickFeat(4, { metaKey: true }); await page.waitForTimeout(200);
  ok('cmd does the same as ctrl, for a Mac', JSON.stringify(await sel()) === '[1,2,3,4]', JSON.stringify(await sel()));
  await clickFeat(2, { ctrlKey: true }); await page.waitForTimeout(200);
  ok('ctrl-clicking one already in takes it back out',
    JSON.stringify(await sel()) === '[1,3,4]', JSON.stringify(await sel()));
  ok('and the one still shown in the inspector is one of those left',
    await page.evaluate(() => window.PB.featureSelection(window.PB.active()).indexOf(window.PB.state.selectedFeature) !== -1));

  console.log('\n== what gets painted ==');
  ok('every selected feature is marked on the map',
    await page.evaluate(() => document.querySelectorAll('#view-map g.feat.sel').length) === 3,
    await page.evaluate(() => document.querySelectorAll('#view-map g.feat.sel').length));
  ok('and each one gets its own arc on the ring',
    await page.evaluate(() => document.querySelectorAll('#view-map .selarc').length) === 3,
    await page.evaluate(() => document.querySelectorAll('#view-map .selarc').length));
  ok('the header counts them and names them', await page.evaluate(() =>
    /3 features selected/.test(document.querySelector('.topbar, header, #header').textContent)),
    await page.evaluate(() => (document.querySelector('.topbar, header, #header').textContent.match(/\d+ features selected[^·]*·[^|]{0,60}/) || [''])[0]));
  ok('the bases add up to the three features, not to one span', await page.evaluate(() => {
    const rec = window.PB.active();
    const want = window.PB.featureSelection(rec)
      .reduce((t, i) => t + window.Bio.featureLength(rec.features[i], rec.seq.length), 0);
    const got = window.PB.featureSelectionPieces(rec).reduce((t, p) => t + (p.end - p.start), 0);
    return want === got && want > 0;
  }));

  await page.evaluate(() => { window.PB.state.view = 'seq'; window.PB.renderAll(); });
  await page.waitForTimeout(500);
  ok('the sequence view highlights all three, not just the last',
    await page.evaluate(() => document.querySelectorAll('#view-seq [style*="var(--sel)"]').length) >= 3,
    await page.evaluate(() => document.querySelectorAll('#view-seq [style*="var(--sel)"]').length));

  await page.evaluate(() => { window.PB.state.mapMode = 'linear'; window.PB.state.view = 'map'; window.PB.renderAll(); });
  await page.waitForTimeout(500);
  ok('the linear map shades all three too', await page.evaluate(() =>
    Array.from(document.querySelectorAll('#view-map rect'))
      .filter(r => (r.getAttribute('fill') || '').indexOf('sel-strong') >= 0).length) === 3,
    await page.evaluate(() => Array.from(document.querySelectorAll('#view-map rect'))
      .filter(r => (r.getAttribute('fill') || '').indexOf('sel-strong') >= 0).length));
  await page.evaluate(() => { window.PB.state.mapMode = 'circular'; window.PB.renderAll(); });
  await page.waitForTimeout(400);

  console.log('\n== from the Features tab ==');
  await page.evaluate(() => { window.PB.state.selectedFeature = null; window.PB.state.selectedFeatures = []; window.PB.state.view = 'features'; window.PB.renderAll(); });
  await page.waitForTimeout(400);
  const rowClick = (i, mods) => page.evaluate(([i, mods]) => {
    document.querySelector('#view-features [data-sel="' + i + '"]')
      .dispatchEvent(new MouseEvent('click', Object.assign({ bubbles: true }, mods)));
  }, [i, mods || {}]);
  await rowClick(0, { ctrlKey: true }); await page.waitForTimeout(300);
  ok('ctrl-clicking a row does not throw you into the sequence tab',
    await page.evaluate(() => window.PB.state.view) === 'features',
    await page.evaluate(() => window.PB.state.view));
  await rowClick(2, { ctrlKey: true }); await page.waitForTimeout(300);
  ok('and builds up a selection down the list', JSON.stringify(await sel()) === '[0,2]', JSON.stringify(await sel()));
  await rowClick(4); await page.waitForTimeout(400);
  ok('a plain click still picks one and opens it', JSON.stringify(await sel()) === '[4]', JSON.stringify(await sel()));
  ok('in the sequence view, as it always did', await page.evaluate(() => window.PB.state.view) === 'seq');

  console.log('\n== the set cannot go stale ==');
  await page.evaluate(() => {
    window.PB.state.view = 'map';
    const rec = window.PB.active();
    window.PB.state.selectedFeatures = [0, 1, 2];
    window.PB.state.selectedFeature = 2;
    /* a feature removed out from under the selection */
    rec.features.splice(1, 1);
    window.PB.renderAll();
  });
  await page.waitForTimeout(400);
  ok('an index past the end of a shortened feature list is dropped', await page.evaluate(() => {
    const rec = window.PB.active();
    return window.PB.featureSelection(rec).every(i => !!rec.features[i]);
  }));
  await page.evaluate(() => { window.PB.state.selectedFeature = null; window.PB.renderAll(); });
  await page.waitForTimeout(300);
  ok('clearing the single selection clears the whole set',
    JSON.stringify(await sel()) === '[]', JSON.stringify(await sel()));
/* Clearing the feature selection is not the same as clearing the selection.
     A plain click on a feature earlier in this run also selected its bases,
     and those bases are still selected — the arc that remains is the base
     range, not a stale feature. Dropping it here would mean a user who
     ctrl-clicked their way out of a feature set silently lost the span they
     had marked. Clearing the range as well leaves no arcs at all. */
  ok('the only arc left belongs to the base range, not to a feature',
    await page.evaluate(() => document.querySelectorAll('#view-map .selarc').length) === 1 &&
    await page.evaluate(() => !!window.PB.state.selection),
    JSON.stringify(await page.evaluate(() => ({
      arcs: document.querySelectorAll('#view-map .selarc').length,
      range: window.PB.state.selection }))));
  await page.evaluate(() => { window.PB.state.selection = null; window.PB.renderAll(); });
  await page.waitForTimeout(350);
  ok('and with the range cleared too the map goes back to no arcs',
    await page.evaluate(() => document.querySelectorAll('#view-map .selarc').length) === 0,
    JSON.stringify(await page.evaluate(() => ({
      arcs: document.querySelectorAll('#view-map .selarc').length,
      sel: window.PB.featureSelection(window.PB.active()),
      range: window.PB.state.selection }))));

  console.log('\n== the image export is unaffected ==');
  ok('a map image drawn with no selection has no arcs in it', await page.evaluate(() => {
    const rec = window.PB.active();
    window.PB.state.selectedFeatures = [0, 1];
    window.PB.state.selectedFeature = 1;
    const svg = window.PB.mapSVGString({ whole: true, theme: 'light', background: 'solid' });
    return svg.indexOf('selarc') === -1;
  }));
  ok('and the selection is still there afterwards',
    JSON.stringify(await sel()) === '[0,1]', JSON.stringify(await sel()));

  ok('nothing threw anywhere in the run', errors.length === 0, errors.slice(0, 3).join(' | '));
  console.log('\n   ' + pass + ' passed, ' + fail + ' failed\n');
  await browser.close();
  process.exit(fail ? 1 : 0);
})().catch(e => { console.error('HARNESS ERROR', e); process.exit(2); });
