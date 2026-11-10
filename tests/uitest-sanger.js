/* ===========================================================================
   Sequencing from stock, through the interface: import an inventory,
   plan a region on a real vector, and check that what comes out is the
   same as what the engine says, that the picture is drawn, that the
   spreadsheet is a real workbook, and that a gap can be filled.
   ========================================================================= */
const { chromium } = require('playwright');
const path = require('path'), fs = require('fs');
let pass = 0, fail = 0;
function ok(name, cond, extra) {
  if (cond) { pass++; console.log('  PASS  ' + name); }
  else { fail++; console.log('  FAIL  ' + name + (extra !== undefined ? '   -> ' + extra : '')); }
}
const VEC = path.resolve(__dirname, 'fixtures/pRSFDuet-1.dna');
const STOCK = path.resolve(__dirname, 'fixtures/Primer_stock.tsv');

(async () => {
  const browser = await chromium.launch();
  const page = await browser.newPage({ viewport: { width: 1500, height: 1000 } });
  const errors = [];
  page.on('pageerror', e => errors.push('pageerror: ' + String(e)));
  page.on('console', m => { if (m.type() === 'error' && !/Failed to load resource/.test(m.text())) errors.push('console: ' + m.text()); });
  await page.goto('file://' + path.resolve('plasmidbench.html'));
  await page.waitForTimeout(700);

  /* a real vector */
  const b64 = fs.readFileSync(VEC).toString('base64');
  await page.evaluate(b64 => {
    const bin = atob(b64), buf = new Uint8Array(bin.length);
    for (let i = 0; i < bin.length; i++) buf[i] = bin.charCodeAt(i);
    const r = window.Bio.parseSnapGene(buf.buffer);
    window.PB.addRecord(Object.assign({}, r, { name: 'pRSFDuet-1', topology: 'circular' }));
    window.PB.renderAll();
  }, b64);
  await page.waitForTimeout(300);

  console.log('\n== importing an inventory through the dialog ==');
  const stockText = fs.readFileSync(STOCK, 'utf8');
  await page.evaluate(() => window.PB.stockImportDialog());
  await page.waitForTimeout(250);
  ok('the stock dialog opens', await page.evaluate(() => !!document.querySelector('#skText')));
  await page.evaluate(t => { document.querySelector('#skText').value = t; }, stockText);
  await page.evaluate(() => {
    const b = Array.from(document.querySelectorAll('#modalRoot button')).find(x => x.textContent.trim() === 'Check it');
    b.click();
  });
  await page.waitForTimeout(600);
  const prev = await page.evaluate(() => {
    const el = document.querySelector('#skPrev');
    return { text: el ? el.textContent.slice(0, 160) : '', rows: el ? el.querySelectorAll('tbody tr').length : 0 };
  });
  /* The inventory is a generated stand-in of the same shape: 874 oligos,
     use. The engine reads it to the same counts outside the browser, so a
     disagreement here is the dialog's doing and not the parser's. */
  ok('it previews what it read', /Read 874 oligos/.test(prev.text), prev.text.slice(0, 90));
  ok('and shows a sample of them', prev.rows === 12, prev.rows);
  /* Nothing is skipped in this inventory once the header is read properly:
     "Sequence 5'->3'" is the sequence column, "Name" is the name, and
     "Where it is" is the location. Before the header reader allowed a unit
     or an arrow after the word, none of those were recognised, the file fell
     through to shape-guessing, and every oligo came in named after its
     freezer position. */
  ok('and nothing in this inventory is skipped', !/skipped \d+ row/.test(prev.text) ||
     /skipped 0 row/.test(prev.text), prev.text.slice(0, 160));
  ok('the oligos keep their own names, not their freezer positions',
     await page.evaluate(() => {
       const ps = window.PB.stockPrimers ? window.PB.stockPrimers() : [];
       return ps.length ? !/^L\d+ \//.test(ps[0].name) : null;
     }) !== false);
  await page.evaluate(() => {
    const b = Array.from(document.querySelectorAll('#modalRoot button')).find(x => x.textContent.trim() === 'Import');
    b.click();
  });
  await page.waitForTimeout(500);
  const held = await page.evaluate(() => window.PB.stockCount());
  ok('the stock is held after importing', held === 874, held);
  /* An oligo inventory of 874 rows is far too big for the localStorage
     mirror, which is the reason the store moved to IndexedDB: this test
     used to read pb.stock out of localStorage and so was testing where the
     list was kept rather than that it was kept. Reloading the page and
     counting again is the thing a person would actually notice, and it
     passes wherever the list ends up living. */
  const persisted = await page.evaluate(() => new Promise(resolve => {
    /* Read it straight out of IndexedDB, in a connection of this test's
       own, so that what is checked is the stored copy and not the page's
       memory. A reload would prove the same thing but would throw away the
       record and selection the rest of this suite goes on to use. */
    const req = indexedDB.open('plasmidbench');
    req.onerror = () => resolve(-1);
    req.onsuccess = () => {
      const db = req.result;
      if (!db.objectStoreNames.contains('settings')) { resolve(-2); return; }
      const g = db.transaction('settings', 'readonly').objectStore('settings').get('stock');
      g.onsuccess = () => resolve(g.result && Array.isArray(g.result.value) ? g.result.value.length : 0);
      g.onerror = () => resolve(-3);
    };
  }));
  ok('and kept for next time, in IndexedDB rather than localStorage', persisted === 874, persisted);
  ok('the inventory is not dumped into the working primer library',
     await page.evaluate(() => (window.PB.state.primers || []).length) < 5);

  console.log('\n== planning a region ==');
  await page.evaluate(() => {
    const P = window.PB;
    P.state.selection = { start: 100, end: 1200 };
    P.closeModal(); P.toolSanger();
  });
  await page.waitForTimeout(350);
  ok('the sequencing tool opens on the selection',
     await page.evaluate(() => { const s = document.querySelector('#sqWhat'); return !!s && s.value === 'sel'; }));
  await page.evaluate(() => {
    const b = Array.from(document.querySelectorAll('#modalRoot button')).find(x => x.textContent.trim() === 'Find primers');
    b.click();
  });
  await page.waitForTimeout(1800);
  const res = await page.evaluate(() => {
    const out = document.querySelector('#sqOut');
    const last = window.PB.sangerLast();
    return {
      text: out ? out.textContent : '',
      chips: Array.from(out.querySelectorAll('.chip')).map(c => c.textContent),
      svg: !!out.querySelector('svg.covbar'),
      lanes: out.querySelectorAll('svg.covbar rect').length,
      tables: out.querySelectorAll('table').length,
      fwd: last.plan.plan.forward.length, rev: last.plan.plan.reverse.length,
      cov: last.plan.coverage, cands: last.plan.candidates.length,
      refused: last.plan.unusable.length
    };
  });
  ok('it finds usable primers in the inventory', res.cands >= 5, res.cands);
  ok('it plans forward reads', res.fwd >= 1, res.fwd);
  ok('and names the coverage on both strands', res.chips.some(c => /both strands/.test(c)),
     JSON.stringify(res.chips));
  ok('the coverage picture is drawn', res.svg && res.lanes >= 2, JSON.stringify({ svg: res.svg, lanes: res.lanes }));
  ok('tables of primers are shown', res.tables >= 1, res.tables);
  ok('every tube says where it is in the freezer',
     await page.evaluate(() => {
       const last = window.PB.sangerLast();
       return last.plan.plan.forward.concat(last.plan.plan.reverse).every(c => !!c.loc);
     }));
  ok('a primer with more than one site is refused, not planned', res.refused >= 1, res.refused);
  ok('and nothing with a repeat site is in the plan',
     await page.evaluate(() => {
       const last = window.PB.sangerLast();
       return last.plan.plan.forward.concat(last.plan.plan.reverse).every(c => c.sites === 1);
     }));
  ok('the T7 promoter primer is among the refusals, because a Duet vector has two of them',
     await page.evaluate(() => window.PB.sangerLast().plan.unusable.some(c => /T7 promoter/i.test(c.name || ''))),
     await page.evaluate(() => JSON.stringify(window.PB.sangerLast().plan.unusable.map(c => c.name))));

  /* what is on screen must be what the engine computed */
  const agree = await page.evaluate(() => {
    const last = window.PB.sangerLast();
    const txt = document.querySelector('#sqOut').textContent;
    const cov = last.plan.coverage;
    const pc = v => Math.round((v / cov.length) * 100) + '%';
    return { both: txt.indexOf('both strands ' + pc(cov.both)) >= 0,
             fwd: txt.indexOf('forward ' + pc(cov.forward)) >= 0,
             /* the screen writes long numbers with thousands separators, so
                compare against the formatted figure as well as the raw one */
             len: txt.indexOf(cov.length + ' bp region') >= 0 ||
                  txt.indexOf(cov.length.toLocaleString('en-US') + ' bp region') >= 0,
             shown: (txt.match(/[\d,]+ bp region/) || [''])[0], want: cov.length };
  });
  ok('the numbers on screen are the engine’s own', agree.both && agree.fwd && agree.len, JSON.stringify(agree));

  console.log('\n== clicking a primer shows where it anneals ==');
  await page.evaluate(() => {
    const b = document.querySelector('#sqOut [data-show]');
    if (b) b.click();
  });
  await page.waitForTimeout(300);
  const selNow = await page.evaluate(() => {
    const last = window.PB.sangerLast(), s = window.PB.state.selection;
    const c = last.plan.candidates[+document.querySelector('#sqOut [data-show]').getAttribute('data-show')];
    return { ok: s && c && s.start === c.site.start && s.end === c.site.end };
  });
  ok('it selects the annealing footprint, not the whole primer', selNow.ok, JSON.stringify(selNow));

  console.log('\n== exporting the plan ==');
  await page.evaluate(() => {
    window.__caught = null;
    window.PB.offerBlob = async (fn, blob) => {
      const b = new Uint8Array(await blob.arrayBuffer());
      let txt = ''; for (let i = 0; i < b.length; i++) txt += String.fromCharCode(b[i]);
      window.__caught = { fn: fn, n: b.length, head: [b[0], b[1]], txt: txt };
      return true;
    };
  });
  await page.evaluate(() => document.querySelector('#sqXls').click());
  await page.waitForTimeout(900);
  const wb = await page.evaluate(() => window.__caught && { fn: window.__caught.fn, n: window.__caught.n, head: window.__caught.head, hasPlan: /sheet1\.xml/.test(window.__caught.txt), hasAbout: /sheet2\.xml/.test(window.__caught.txt) });
  ok('the plan exports as a workbook', !!wb && /\.xlsx$/.test(wb.fn) && wb.head[0] === 0x50 && wb.head[1] === 0x4b,
     wb ? JSON.stringify(wb).slice(0, 120) : 'nothing offered');
  ok('with a plan sheet and an about sheet', wb && wb.hasPlan && wb.hasAbout, JSON.stringify(wb && { p: wb.hasPlan, a: wb.hasAbout }));

  console.log('\n== a region the freezer cannot reach ==');
  await page.evaluate(() => {
    const P = window.PB;
    P.closeModal();
    P.state.selection = { start: 2400, end: 3400 };   /* out in the backbone */
    P.toolSanger();
  });
  await page.waitForTimeout(300);
  await page.evaluate(() => {
    const b = Array.from(document.querySelectorAll('#modalRoot button')).find(x => x.textContent.trim() === 'Find primers');
    b.click();
  });
  await page.waitForTimeout(1800);
  const gapState = await page.evaluate(() => {
    const out = document.querySelector('#sqOut');
    const last = window.PB.sangerLast();
    return { hasGapBtn: !!out.querySelector('[data-fill]'),
             gapsF: last.plan.gaps.forward.length, gapsR: last.plan.gaps.reverse.length,
             text: out.textContent.slice(0, 200) };
  });
  ok('a stretch nothing reads is reported as a gap', gapState.gapsF + gapState.gapsR >= 1,
     JSON.stringify({ f: gapState.gapsF, r: gapState.gapsR }));
  if (gapState.hasGapBtn) {
    await page.evaluate(() => document.querySelector('#sqOut [data-fill]').click());
    await page.waitForTimeout(1200);
    const designed = await page.evaluate(() => {
      const t = document.querySelector('#modalRoot');
      const rows = Array.from(t.querySelectorAll('tbody tr')).map(r =>
        Array.from(r.querySelectorAll('td')).map(c => c.textContent.trim()));
      return { title: (t.querySelector('h3') || {}).textContent || '', rows: rows };
    });
    ok('it designs a primer to order for that gap', designed.rows.length >= 1 &&
       /^[ACGT]{18,26}$/.test(designed.rows[0][1] || ''), JSON.stringify(designed.rows[0] || designed.title));
    const uniq = await page.evaluate(() => {
      const rec = window.PB.active();
      const rows = Array.from(document.querySelectorAll('#modalRoot tbody tr'));
      return rows.map(r => {
        const seq = r.querySelectorAll('td')[1].textContent.trim();
        return window.Bio.primerSites(rec.seq, seq, { circular: true }).length;
      });
    });
    ok('and each designed primer anneals in exactly one place', uniq.every(n => n === 1), JSON.stringify(uniq));
    await page.evaluate(() => {
      const b = Array.from(document.querySelectorAll('#modalRoot button')).find(x => /Add to my primers/.test(x.textContent));
      if (b) b.click();
    });
    await page.waitForTimeout(500);
    ok('it can be added to the working primers', await page.evaluate(() => (window.PB.state.primers || []).length) >= 1);
  }

  console.log('\n== settings are respected ==');
  await page.evaluate(() => { window.PB.closeModal(); window.PB.state.selection = { start: 100, end: 1200 }; window.PB.toolSanger(); });
  await page.waitForTimeout(300);
  await page.fill('#sqGood', '300');
  await page.evaluate(() => {
    const b = Array.from(document.querySelectorAll('#modalRoot button')).find(x => x.textContent.trim() === 'Find primers');
    b.click();
  });
  await page.waitForTimeout(1800);
  const tight = await page.evaluate(() => {
    const l = window.PB.sangerLast();
    return { good: l.plan.settings.good, both: l.plan.coverage.both, fwd: l.plan.coverage.forward };
  });
  ok('a shorter confident read is used when asked', tight.good === 300, tight.good);

  await page.evaluate(() => window.PB.closeModal());
  ok('nothing threw anywhere in the run', errors.length === 0, errors.slice(0, 3).join(' | '));
  console.log('\n   ' + pass + ' passed, ' + fail + ' failed\n');
  await browser.close();
  process.exit(fail ? 1 : 0);
})().catch(e => { console.error('HARNESS ERROR', e); process.exit(2); });
