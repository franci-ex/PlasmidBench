/* ===========================================================================
   Expression strategy, through the interface, against the real collection.
   ========================================================================= */
const { chromium } = require('playwright');
const path = require('path'), fs = require('fs');
let pass = 0, fail = 0;
function ok(name, cond, extra) {
  if (cond) { pass++; console.log('  PASS  ' + name); }
  else { fail++; console.log('  FAIL  ' + name + (extra !== undefined ? '   -> ' + extra : '')); }
}
/* The lab's own plasmid collection. It is what makes this suite worth
   running: whether a tag is "on the shelf" and whether a chaperone plasmid
   clashes with the host both depend on what is actually in the freezer. */
const STOCK = path.resolve(__dirname, 'fixtures/Plasmid_stock.tsv');
const GB1 = 'MTYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE';

(async () => {
  const browser = await chromium.launch();
  const page = await browser.newPage({ viewport: { width: 1500, height: 1050 } });
  const errors = [];
  page.on('pageerror', e => errors.push('pageerror: ' + String(e)));
  page.on('console', m => { if (m.type() === 'error' && !/Failed to load resource/.test(m.text())) errors.push('console: ' + m.text()); });
  await page.goto('file://' + path.resolve('plasmidbench.html'));
  await page.waitForTimeout(700);

  console.log('\n== importing the plasmid collection ==');
  const tsv = fs.existsSync(STOCK) ? fs.readFileSync(STOCK, 'utf8') : null;
  if (!tsv) { console.log('  (no stock file; skipping)'); }
  else {
    await page.evaluate(() => window.PB.plasmidImportDialog());
    await page.waitForTimeout(250);
    ok('the plasmid stock dialog opens', await page.evaluate(() => !!document.querySelector('#psText')));
    await page.evaluate(t => { document.querySelector('#psText').value = t; }, tsv);
    await page.evaluate(() => Array.from(document.querySelectorAll('#modalRoot button'))
      .find(b => b.textContent.trim() === 'Check it').click());
    await page.waitForTimeout(900);
    const prev = await page.evaluate(() => {
      const el = document.querySelector('#psPrev');
      return { text: el ? el.textContent.slice(0, 220) : '', rows: el ? el.querySelectorAll('tbody tr').length : 0 };
    });
    ok('it reports what it read', /Read \d+ plasmids/.test(prev.text), prev.text.slice(0, 80));
    ok('and recognises most of the backbones', /9[0-9]%\) are backbones/.test(prev.text), prev.text.slice(0, 160));
    ok('summarising them by vector family', prev.rows > 5, prev.rows);
    await page.evaluate(() => Array.from(document.querySelectorAll('#modalRoot button'))
      .find(b => b.textContent.trim() === 'Import').click());
    await page.waitForTimeout(700);
    ok('the collection is held', await page.evaluate(() => window.PB.plasmidCount()) > 500,
       await page.evaluate(() => window.PB.plasmidCount()));
  }

  console.log('\n== reading a sequence ==');
  await page.evaluate(() => { window.PB.closeModal(); window.PB.toolExpression(); });
  await page.waitForTimeout(350);
  ok('the tool opens', await page.evaluate(() => !!document.querySelector('#exSrc')));
  await page.selectOption('#exSrc', 'paste');
  await page.waitForTimeout(150);
  await page.evaluate(g => { document.querySelector('#exPaste').value = g; }, GB1);
  await page.evaluate(() => Array.from(document.querySelectorAll('#modalRoot button'))
    .find(b => b.textContent.trim() === 'Work it out').click());
  await page.waitForTimeout(1200);
  const easy = await page.evaluate(() => {
    const out = document.querySelector('#exOut');
    const l = window.PB.expressionLast();
    return { text: out.textContent, band: l.outlook.band, tables: out.querySelectorAll('table').length,
             first: l.strategy.options[0].tag,
             mw: l.outlook.props.mwKDa, pI: l.outlook.props.pI };
  });
  ok('a small soluble domain reads as straightforward', easy.band === 'straightforward', easy.band);
  ok('its numbers are reported', Math.abs(easy.mw - 6.2) < 0.2 && Math.abs(easy.pI - 4.55) < 0.1,
     JSON.stringify({ mw: easy.mw, pI: easy.pI }));
  ok('a plain tag is the first suggestion', easy.first === 'His', easy.first);
  ok('the reading, the constructs and the conditions are all shown', easy.tables >= 3, easy.tables);
  ok('and the limits of the method are stated', /composition cannot tell you/i.test(easy.text));

  console.log('\n== something that will fight back ==');
  await page.evaluate(() => {
    let s = 7, r = () => (s = (s * 1103515245 + 12345) & 0x7fffffff) / 0x7fffffff;
    let p = 'M'; for (let i = 0; i < 700; i++) p += 'VILFAWMVILFAGVILFCY'[Math.floor(r() * 19)];
    document.querySelector('#exPaste').value = p;
  });
  await page.evaluate(() => Array.from(document.querySelectorAll('#modalRoot button'))
    .find(b => b.textContent.trim() === 'Work it out').click());
  await page.waitForTimeout(1400);
  const hard = await page.evaluate(() => {
    const l = window.PB.expressionLast();
    const out = document.querySelector('#exOut');
    return { band: l.outlook.band, first: l.strategy.options[0].tag,
             shelf: l.strategy.options[0].have.length,
             conditions: l.strategy.conditions.map(c => c.what),
             text: out.textContent };
  });
  ok('a hydrophobic 700-mer reads as hard', hard.band === 'hard', hard.band);
  ok('MBP is the first thing suggested for it', hard.first === 'MBP', hard.first);
  ok('and it is found in the collection', hard.shelf > 0, hard.shelf);
  ok('with a cold, light induction', hard.conditions.some(c => /cold/i.test(c)) &&
     hard.conditions.some(c => /lightly/i.test(c)), JSON.stringify(hard.conditions));

  console.log('\n== co-expression, and the conflicts that matter ==');
  await page.evaluate(() => {
    document.querySelector('#exPartners').value = 'cyclin A\nCDK2';
    const h = document.querySelector('#exHost');
    h.value = 'BL21(DE3)';
  });
  await page.evaluate(() => Array.from(document.querySelectorAll('#modalRoot button'))
    .find(b => b.textContent.trim() === 'Work it out').click());
  await page.waitForTimeout(1400);
  const co = await page.evaluate(() => {
    const l = window.PB.expressionLast();
    return { n: l.coex.assigned.length,
             oris: l.coex.assigned.map(a => a.slot.ori),
             markers: l.coex.assigned.map(a => a.slot.marker),
             problems: l.coex.compatibility.problems.length,
             text: document.querySelector('#exOut').textContent };
  });
  ok('three partners are each given a backbone', co.n === 3, co.n);
  ok('all on different origins', new Set(co.oris).size === 3, JSON.stringify(co.oris));
  ok('and different drugs', new Set(co.markers).size === 3, JSON.stringify(co.markers));
  ok('with nothing left conflicting', co.problems === 0, co.problems);

  /* the trap: chaperones in a host that already carries a Cm plasmid */
  /* Pick the host from the list rather than typing its name into the select:
     the option values are 'n:<strain>' for the built-in strains and
     'i:<index>' for one off the lab's own cell sheet, so assigning a bare
     name leaves the select empty and the plan ends up with no host at all —
     and then of course nothing clashes. */
  const hostPick = await page.evaluate(() => {
    const sel = document.querySelector('#exHost');
    const o = Array.from(sel.options).find(x => /Rosetta/i.test(x.textContent));
    if (!o) return null;
    sel.value = o.value;
    sel.dispatchEvent(new Event('change', { bubbles: true }));
    document.querySelector('#exChap').checked = true;
    document.querySelector('#exPartners').value = '';
    return { value: sel.value, label: o.textContent };
  });
  ok('a Rosetta host can be chosen from the strain list', !!hostPick,
     JSON.stringify(hostPick));
  await page.evaluate(() => Array.from(document.querySelectorAll('#modalRoot button'))
    .find(b => b.textContent.trim() === 'Work it out').click());
  await page.waitForTimeout(1400);
  const clash = await page.evaluate(() => {
    const l = window.PB.expressionLast();
    return { problems: l.coex.compatibility.problems.map(p => p.kind),
             text: document.querySelector('#exOut').textContent };
  });
  ok('chaperones in a Rosetta host are refused', clash.problems.indexOf('helper') >= 0,
     JSON.stringify(clash.problems));
  /* The wording now names the host's own plasmid and what it carries rather
     than only the origin and the drug, which is the useful form: it tells you
     which of the two to drop. */
  ok('and the clash is spelled out on screen',
     /both want p15A and chloramphenicol/.test(clash.text) &&
     /rare-codon tRNAs/.test(clash.text),
     (clash.text.match(/.{0,40}p15A.{0,90}/) || [''])[0]);
  ok('and it says what to do about it',
     /BL21\(DE3\) or LOBSTR/.test(clash.text) && /drop the chaperones/.test(clash.text),
     (clash.text.match(/.{0,30}LOBSTR.{0,80}/) || [''])[0]);

  console.log('\n== exporting it ==');
  await page.evaluate(() => {
    window.__caught = null;
    window.PB.offerBlob = async (fn, blob) => {
      const b = new Uint8Array(await blob.arrayBuffer());
      let t = ''; for (let i = 0; i < b.length; i++) t += String.fromCharCode(b[i]);
      window.__caught = { fn: fn, n: b.length, head: [b[0], b[1]], txt: t };
      return true;
    };
  });
  await page.evaluate(() => document.querySelector('#exXls').click());
  await page.waitForTimeout(900);
  const wb = await page.evaluate(() => window.__caught && { fn: window.__caught.fn, n: window.__caught.n,
    head: window.__caught.head, sheets: (window.__caught.txt.match(/sheet\d+\.xml/g) || []).length });
  ok('the plan exports as a workbook', !!wb && /\.xlsx$/.test(wb.fn) && wb.head[0] === 0x50,
     wb ? JSON.stringify({ fn: wb.fn, n: wb.n }) : 'nothing offered');
  ok('with a sheet for each part of it', wb && wb.sheets >= 4, wb && wb.sheets);

  await page.evaluate(() => window.PB.closeModal());
  ok('nothing threw anywhere in the run', errors.length === 0, errors.slice(0, 3).join(' | '));
  console.log('\n   ' + pass + ' passed, ' + fail + ' failed\n');
  await browser.close();
  process.exit(fail ? 1 : 0);
})().catch(e => { console.error('HARNESS ERROR', e); process.exit(2); });
