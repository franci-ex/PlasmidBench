/* ===========================================================================
   Replacing a region that crosses base 1.

   A circular plasmid has no end, but its coordinates do, and a gene put in
   by a previous assembly very often straddles the join — which makes the
   gene somebody wants to swap out exactly the one most likely to wrap.
   Given only a From and a To, such a region arrives with To below From.
   Read as an ordinary range that is empty, and the vector is opened with
   nothing taken out: the primers amplify the whole plasmid, the insert goes
   in beside the old gene, and the product carries both. It looks right —
   there are primers, there is a map, the insert is really in it — and the
   only visible sign is that the product is longer than it should be.

   Every case here is built so the arithmetic is known in advance: a
   4,842 bp plasmid carrying an 831 bp CDS at 4,058..4,842 + 1..46.
   ========================================================================= */
const { chromium } = require('playwright');
const path = require('path');
let pass = 0, fail = 0;
function ok(name, cond, extra) {
  if (cond) { pass++; console.log('  PASS  ' + name); }
  else { fail++; console.log('  FAIL  ' + name + (extra !== undefined ? '   -> ' + extra : '')); }
}

const L = 4842, GENE = 831, INS = 'ATGGGTTTTCTGCGTTTTCAGGCTAGCGTGCCCTCAGGTTTAACCGGTAAA';

(async () => {
  const browser = await chromium.launch();
  const page = await browser.newPage({ viewport: { width: 1500, height: 1150 } });
  const errors = [];
  page.on('pageerror', e => errors.push('pageerror: ' + String(e)));
  page.on('console', m => { if (m.type() === 'error' && !/Failed to load resource|ERR_(TUNNEL|NAME|INTERNET|CONNECTION)/.test(m.text())) errors.push(m.text()); });
  await page.goto('file://' + path.resolve(__dirname, '..', 'plasmidbench.html'));
  await page.waitForTimeout(800);

  /* the construct, and where the wrapping gene sits in it */
  const made = await page.evaluate(([L, GENE]) => {
    const P = window.PB;
    let x = 7; const r = () => (x = (x * 1103515245 + 12345) & 0x7fffffff) / 0x7fffffff;
    const dna = n => { let s = ''; for (let i = 0; i < n; i++) s += 'ACGT'[Math.floor(r() * 4)]; return s; };
    const cds = 'ATG' + dna(GENE - 6).replace(/TAA|TAG|TGA/g, 'CCC') + 'TAA';
    const seq = cds.slice(GENE - 46) + dna(L - GENE) + cds.slice(0, GENE - 46);
    window.__cds = cds;
    const rec = P.makeRecord({ name: 'YTC-0528 pRSFDuet', seq: seq, topology: 'circular', features: [
      { start: 4057, end: 4057 + GENE, strand: 1, type: 'CDS', label: 'CDK2', quals: { codon_start: '1' } },
      { start: 1200, end: 2000, strand: 1, type: 'CDS', label: 'lacI', quals: {} },
      { start: 2400, end: 3000, strand: 1, type: 'CDS', label: 'KanR', quals: {} } ] });
    P.state.records = [rec]; P.state.activeId = rec.id; P.renderAll();
    return { len: seq.length };
  }, [L, GENE]);
  ok('the test construct is the size it should be', made.len === L, made.len);

  const select = (start, end) => page.evaluate(([s, e]) => {
    window.PB.state.selection = { start: s, end: e };
    window.PB.renderAll();
    return window.PB.normSel(window.PB.active(), window.PB.state.selection);
  }, [start, end]);

  async function design(opts) {
    opts = opts || {};
    await page.evaluate(() => {
      const m = document.querySelector('#modalRoot');
      if (m && !m.hidden) window.PB.closeModal();
      /* back to the vector: saving a product opens it */
      const v = window.PB.state.records.find(r => r.name === 'YTC-0528 pRSFDuet');
      if (v) { window.PB.state.activeId = v.id; window.PB.renderAll(); }
    });
    await page.waitForTimeout(250);
    if (opts.sel) await select(opts.sel[0], opts.sel[1]);
    await page.evaluate(() => window.PB.openTools());
    await page.waitForTimeout(250);
    await page.locator('.tool-card', { hasText: 'Seamless assembly' }).first().click();
    await page.waitForTimeout(550);
    await page.selectOption('#gbIns_src', '');
    await page.fill('#gbIns_txt', opts.insert === undefined ? INS : opts.insert);
    if (opts.insert !== '') await page.fill('#gbInsName', 'SPDYC');
    await page.waitForTimeout(650);
    const got = await page.evaluate(src => {
      const sel = document.querySelector('#gbWhere');
      const o = Array.from(sel.options).find(x => new RegExp(src).test(x.textContent));
      if (!o) return null;
      sel.value = o.value; sel.dispatchEvent(new Event('change', { bubbles: true }));
      return o.textContent;
    }, opts.where.source);
    if (!got) throw new Error('no placement matching ' + opts.where + ' in that list');
    await page.waitForTimeout(700);
    if (opts.from != null) { await page.fill('#gbFrom', String(opts.from)); await page.waitForTimeout(250); }
    if (opts.to != null) { await page.fill('#gbTo', String(opts.to)); await page.waitForTimeout(650); }
    return { label: got, text: await page.evaluate(() => document.querySelector('#gbLive').textContent.replace(/\s+/g, ' ')) };
  }
  const saveAndRead = async () => {
    await page.click('#gbSave'); await page.waitForTimeout(800);
    return page.evaluate(ins => {
      const rec = window.PB.active(), d = rec.seq + rec.seq;
      return { name: rec.name, len: rec.seq.length,
               hasInsert: d.indexOf(ins) >= 0,
               oldGeneStillThere: d.indexOf(window.__cds.slice(100, 400)) >= 0,
               feats: rec.features.map(f => f.label) };
    }, INS);
  };

  console.log('\n== the gene that crosses base 1 ==');
  const sel = await select(4057, 4057 + GENE);
  ok('the selection is seen to wrap', sel.wraps === true && sel.length === GENE, JSON.stringify(sel));
  const wrapRun = await design({ where: /replace what is selected on the map/ });
  ok('the placement names both halves of it',
     /4,058–4,842 and 1–46 \(across base 1\), 831 bp/.test(wrapRun.label), wrapRun.label);
  ok('it removes the region rather than opening beside it',
     /removing 831 bp/.test(wrapRun.text), wrapRun.text.slice(0, 160));
  ok('and does not claim to have removed nothing',
     !/nothing removed/.test(wrapRun.text), wrapRun.text.slice(0, 160));
  const wrapLen = L - GENE + INS.length;
  ok('the product is vector minus the gene plus the insert',
     wrapRun.text.indexOf(wrapLen.toLocaleString('en-US') + ' bp circular') >= 0,
     (wrapRun.text.match(/[\d,]+ bp circular/) || [''])[0] + '   want ' + wrapLen);

  const prod = await saveAndRead();
  ok('the saved construct is that length', prod.len === wrapLen, prod.len + ' vs ' + wrapLen);
  ok('the insert really is in it', prod.hasInsert);
  ok('the gene that was replaced is gone from the sequence', !prod.oldGeneStillThere);
  ok('and gone from the feature list too',
     prod.feats.indexOf('CDK2') < 0, prod.feats.join(', '));
  ok('while the rest of the plasmid is untouched',
     prod.feats.indexOf('lacI') >= 0 && prod.feats.indexOf('KanR') >= 0, prod.feats.join(', '));

  console.log('\n== an ordinary selection still behaves ==');
  const plainRun = await design({ sel: [1200, 2000], where: /replace what is selected on the map/ });
  ok('it removes exactly what was selected', /removing 800 bp/.test(plainRun.text), plainRun.text.slice(0, 140));
  const plainLen = L - 800 + INS.length;
  ok('and the product is the right size',
     plainRun.text.indexOf(plainLen.toLocaleString('en-US') + ' bp circular') >= 0,
     (plainRun.text.match(/[\d,]+ bp circular/) || [''])[0] + '   want ' + plainLen);
  const prod2 = await saveAndRead();
  ok('the gene it replaced is the one that went',
     prod2.feats.indexOf('lacI') < 0 && prod2.feats.indexOf('CDK2') >= 0, prod2.feats.join(', '));

  console.log('\n== a range typed the same way round ==');
  const typed = await design({ sel: [0, 0], where: /a range I type/, from: 4058, to: 46 });
  ok('a From above a To runs the range across base 1, as it does everywhere else',
     /removing 831 bp/.test(typed.text) && /across base 1/.test(typed.text),
     typed.text.slice(0, 180));
  ok('and not as an opening that removes nothing',
     !/nothing removed/.test(typed.text), typed.text.slice(0, 140));

  console.log('\n== selecting the whole plasmid is refused, not silently ignored ==');
  const whole = await design({ sel: [0, L], where: /replace what is selected on the map/ });
  ok('it says there would be no backbone left',
     /no backbone left to amplify/.test(whole.text), whole.text.slice(0, 200));
  ok('and offers nothing to save',
     await page.evaluate(() => !document.querySelector('#gbSave')));

  console.log('\n== leaving the wrapping region out, with no insert ==');
  const vecOnly = await design({ sel: [4057, 4057 + GENE], insert: '',
                                 where: /everything except what is selected on the map/ });
  ok('the backbone alone is the plasmid minus that region',
     vecOnly.text.indexOf((L - GENE).toLocaleString('en-US') + ' bp') >= 0,
     (vecOnly.text.match(/[\d,]+ bp/) || [''])[0] + '   want ' + (L - GENE));
  ok('and it is still described as removing 831 bp',
     /831 bp/.test(vecOnly.text), vecOnly.text.slice(0, 160));

  console.log('\n== a linear construct has no base 1 to cross ==');
  await page.evaluate(() => {
    const m = document.querySelector('#modalRoot'); if (m && !m.hidden) window.PB.closeModal();
    const v = window.PB.state.records.find(r => r.name === 'YTC-0528 pRSFDuet');
    v.topology = 'linear'; window.PB.state.activeId = v.id; window.PB.state.selection = null;
    window.PB.renderAll();
  });
  await page.waitForTimeout(300);
  await page.evaluate(() => window.PB.openTools());
  await page.waitForTimeout(250);
  await page.locator('.tool-card', { hasText: 'Seamless assembly' }).first().click();
  await page.waitForTimeout(600);
  const linOpts = await page.$$eval('#gbWhere option', os => os.map(o => o.textContent));
  /* There is nothing to open on a linear backbone — it has two ends
     already — so no range placement is offered and the wrap question
     cannot arise. */
  ok('a linear backbone is used as it is, with no region to remove',
     linOpts.length === 1 && /already linear/.test(linOpts[0]), JSON.stringify(linOpts));

  console.log('\n== console ==');
  ok('no page errors anywhere in that', errors.length === 0, errors.slice(0, 3).join(' | '));

  console.log('\n   ' + pass + ' passed, ' + fail + ' failed\n');
  await browser.close();
  process.exit(fail ? 1 : 0);
})();
