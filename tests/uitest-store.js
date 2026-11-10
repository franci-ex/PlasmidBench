/* ===========================================================================
   Where the work is kept.

   Until now the library lived in localStorage: one synchronous string store,
   five to ten megabytes for the whole origin, settings and sequences
   together, and a write that does not fit throws on the spot. For a page
   that is meant to run on a lab machine and hold a few dozen plasmids and a
   few thousand oligos, that ceiling is reachable in an afternoon, and the
   failure is silent in the worst way: the save appears to happen.

   So it is IndexedDB now, and these are the tests that say so. They are
   written against what a person would notice rather than against the API:
   a construct survives a reload; a construct too big for the old store
   survives too; the old store's contents are brought across rather than
   abandoned; and when the browser refuses to keep anything at all — a
   private window, blocked site data — the page says that plainly instead of
   claiming the work is saved.

   The page is loaded from a file:// URL, which is its own origin, so the
   database is this test's alone.
   ========================================================================= */
const { chromium } = require('playwright');
const path = require('path');
let pass = 0, fail = 0;
function ok(name, cond, extra) {
  if (cond) { pass++; console.log('  PASS  ' + name); }
  else { fail++; console.log('  FAIL  ' + name + (extra !== undefined ? '   -> ' + extra : '')); }
}

const PAGE = 'file://' + path.resolve(__dirname, '..', 'plasmidbench.html');

/* A sequence larger than the host store's 200,000-base ceiling and larger
   than a comfortable share of localStorage: this is the case the change is
   for. Built in the page to keep it out of the test's own memory. */
const BIG_BP = 400000;

(async () => {
  const browser = await chromium.launch();
  const ctx = await browser.newContext({ viewport: { width: 1400, height: 1000 } });
  const page = await ctx.newPage();
  const errors = [];
  page.on('pageerror', e => errors.push('pageerror: ' + String(e)));
  page.on('console', m => {
    if (m.type() === 'error' && !/Failed to load resource|ERR_(TUNNEL|NAME|INTERNET|CONNECTION)/.test(m.text()))
      errors.push(m.text());
  });

  console.log('\n== the store opens, and says what it did ==');
  await page.goto(PAGE);
  await page.waitForTimeout(1200);
  const opened = await page.evaluate(() => ({
    present: !!window.PB.store,
    available: window.PB.store && window.PB.store.isAvailable(),
    report: window.PB.state.storeReport || null,
    status: (document.querySelector('#storeState') || {}).textContent || ''
  }));
  ok('IndexedDB is open', opened.present && opened.available === true, JSON.stringify(opened.report));
  ok('and the status line does not claim an account it has not got',
     !/your account/i.test(opened.status), opened.status);
  ok('it says where the work is instead', /saved in this browser/i.test(opened.status), opened.status);

  console.log('\n== a construct survives a reload ==');
  await page.evaluate(() => {
    const P = window.PB;
    const rec = P.makeRecord({ name: 'KeptAcrossReload', seq: 'ATGAAACGCATTAGCACCACC'.repeat(40),
                               topology: 'circular', features: [] });
    P.addRecord(rec);
    P.persist(rec);
  });
  await page.waitForTimeout(600);
  await page.reload();
  await page.waitForTimeout(1400);
  const back = await page.evaluate(() => {
    const r = window.PB.state.records.find(x => x.name === 'KeptAcrossReload');
    return r ? { len: r.seq.length, topo: r.topology } : null;
  });
  ok('it is there after the reload', !!back, JSON.stringify(back));
  ok('with its sequence intact', back && back.len === 840, back && back.len);
  ok('and its topology', back && back.topo === 'circular', back && back.topo);

  console.log('\n== a construct the old store could never have held ==');
  const big = await page.evaluate(async (bp) => {
    const P = window.PB;
    let s = ''; while (s.length < bp) s += 'ACGTTGCAAGGCCTTAACGG';
    s = s.slice(0, bp);
    const rec = P.makeRecord({ name: 'TooBigForLocalStorage', seq: s, topology: 'circular', features: [] });
    P.addRecord(rec); P.persist(rec);
    await new Promise(r => setTimeout(r, 900));
    /* what localStorage holds now, for comparison: the old single-key array
       either refused this or dropped it */
    let lsLen = 0;
    try { lsLen = (localStorage.getItem('pb.library.v1') || '').length; } catch (e) { lsLen = -1; }
    const rows = await P.store.allRecords();
    const row = rows.find(x => x.name === 'TooBigForLocalStorage');
    return { stored: !!row, storedLen: row ? row.seq.length : 0, lsLen: lsLen, asked: bp };
  }, BIG_BP);
  ok('IndexedDB keeps all ' + BIG_BP.toLocaleString() + ' bases of it',
     big.stored && big.storedLen === big.asked, JSON.stringify(big));
  await page.reload();
  await page.waitForTimeout(1600);
  const bigBack = await page.evaluate(() => {
    const r = window.PB.state.records.find(x => x.name === 'TooBigForLocalStorage');
    return r ? r.seq.length : 0;
  });
  ok('and it is still whole after a reload', bigBack === BIG_BP, bigBack);

  console.log('\n== what an earlier version left behind is brought across ==');
  /* A fresh origin holding only the old localStorage shapes: the library
     under one key, and a couple of settings. */
  const ctx2 = await browser.newContext({ viewport: { width: 1400, height: 1000 } });
  const page2 = await ctx2.newPage();
  page2.on('pageerror', e => errors.push('pageerror(2): ' + String(e)));
  await page2.goto(PAGE);
  await page2.waitForTimeout(400);
  await page2.evaluate(async () => {
    await window.PB.store._reset();
    localStorage.clear();
    localStorage.setItem('pb.library.v1', JSON.stringify([
      { id: 'legacy-1', name: 'FromTheOldStore', seq: 'ATGCGTCGCGTTAGTCAGGCTAGC', topology: 'linear',
        features: [], updatedAt: 1 }
    ]));
    localStorage.setItem('pb.conditions', JSON.stringify({ naPlus: 50, mgPlus: 2.5 }));
  });
  await page2.reload();
  await page2.waitForTimeout(1600);
  const moved = await page2.evaluate(async () => {
    const rows = await window.PB.store.allRecords();
    return {
      report: window.PB.state.storeReport,
      inStore: rows.some(r => r.id === 'legacy-1'),
      onScreen: window.PB.state.records.some(r => r.name === 'FromTheOldStore'),
      setting: window.PB.lsGet('conditions', null),
      status: (document.querySelector('#storeState') || {}).textContent || ''
    };
  });
  ok('the old library is now a row in IndexedDB', moved.inStore, JSON.stringify(moved.report));
  ok('and it is on screen', moved.onScreen);
  ok('the settings came across too', moved.setting && moved.setting.naPlus === 50,
     JSON.stringify(moved.setting));
  ok('and the move is reported rather than done silently',
     /moved out of the old store/.test(moved.status), moved.status);
  /* the old key is left alone: an older copy of the page opened on the same
     browser must still find its library */
  const kept = await page2.evaluate(() => !!localStorage.getItem('pb.library.v1'));
  ok('the old key is left where it was, for an older copy of the page', kept);
  await ctx2.close();

  console.log('\n== a browser that will not keep anything says so ==');
  const ctx3 = await browser.newContext({ viewport: { width: 1400, height: 1000 } });
  const page3 = await ctx3.newPage();
  /* indexedDB.open made to fail the way a private window makes it fail */
  await page3.addInitScript(() => {
    try {
      Object.defineProperty(window, 'indexedDB', {
        configurable: true,
        get: () => ({ open: () => { throw new Error('access denied'); },
                      deleteDatabase: () => ({}) })
      });
    } catch (e) { /* nothing to do */ }
  });
  await page3.goto(PAGE);
  await page3.waitForTimeout(1200);
  const refused = await page3.evaluate(() => ({
    available: window.PB.store.isAvailable(),
    report: window.PB.state.storeReport,
    status: (document.querySelector('#storeState') || {}).textContent || '',
    cls: (document.querySelector('#storeState') || {}).className || ''
  }));
  ok('the page still loads', refused.report != null, JSON.stringify(refused.report));
  ok('it knows the store is unavailable', refused.available === false, refused.available);
  ok('it says the browser will not keep data', /will not keep data/i.test(refused.status), refused.status);
  ok('and tells you to export before closing', /export/i.test(refused.status), refused.status);
  ok('marked as a warning, not as ordinary status', /warn/.test(refused.cls), refused.cls);
  /* and settings still work, from localStorage alone, so the page is usable */
  const stillWorks = await page3.evaluate(() => {
    window.PB.lsSet('probe', { a: 1 });
    return window.PB.lsGet('probe', null);
  });
  ok('settings still read and write without it', stillWorks && stillWorks.a === 1,
     JSON.stringify(stillWorks));
  await ctx3.close();

  console.log('\n== saving a file with no host at all ==');
  /* On a lab machine there is no window.claude. The page must use the
     browser's own download rather than telling the person saving is
     unavailable. */
  const dl = await page.evaluate(() => ({ standalone: window.PB.standalone(),
                                          hasNative: typeof window.PB.nativeSave === 'function' }));
  ok('the page knows it is standalone', dl.standalone === true);
  ok('and has a download path of its own', dl.hasNative);
  const download = await Promise.all([
    page.waitForEvent('download', { timeout: 8000 }).catch(() => null),
    page.evaluate(() => window.PB.offerFile('probe.gb', 'LOCUS probe 10 bp\n//\n', 'probe'))
  ]);
  ok('a save produces a real download, with no dialog',
     download[0] !== null && download[1] === true,
     download[0] ? download[0].suggestedFilename() : 'no download event');
  const noModal = await page.evaluate(() => {
    const m = document.querySelector('#modalRoot');
    return !m || m.hidden;
  });
  ok('and nothing claims saving is unavailable here', noModal);

  console.log('\n== console ==');
  ok('no page errors anywhere in that', errors.length === 0, errors.slice(0, 3).join(' | '));

  console.log('\n   ' + pass + ' passed, ' + fail + ' failed\n');
  await browser.close();
  process.exit(fail ? 1 : 0);
})();
