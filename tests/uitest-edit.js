// In-place editing and undo/redo, through the real dialogs, keys and buttons.
const { chromium } = require('playwright');
const path = require('path');
const B = require('../src/engine.js'); require('../src/engine2.js');
let pass = 0, fail = 0;
function ok(name, cond, extra) {
  if (cond) { pass++; console.log('  PASS  ' + name); }
  else { fail++; console.log('  FAIL  ' + name + (extra !== undefined ? '   -> ' + extra : '')); }
}
function prng(seed) { let x = seed; return () => (x = (x * 1103515245 + 12345) & 0x7fffffff) / 0x7fffffff; }
function randDNA(n, seed) { const r = prng(seed); let s = ''; for (let i = 0; i < n; i++) s += 'ACGT'[Math.floor(r() * 4)]; return s; }

const SEQ = randDNA(1200, 5);
/* a feature that spans base 1, one in the middle, one near the end */
const FEATS = [
  { start: 1150, end: 40, strand: 1, type: 'CDS', label: 'spanner' },
  { start: 300, end: 420, strand: -1, type: 'misc_feature', label: 'middle' },
  { start: 1000, end: 1100, strand: 1, type: 'promoter', label: 'late' }
];

(async () => {
  const browser = await chromium.launch();
  const page = await browser.newPage({ viewport: { width: 1400, height: 1000 } });
  const errors = [];
  page.on('pageerror', e => errors.push(String(e)));
  page.on('console', m => {
    if (m.type() === 'error' && !/Failed to load resource|ERR_(TUNNEL|NAME|INTERNET|CONNECTION)/.test(m.text())) errors.push(m.text());
  });
  await page.goto('file://' + path.resolve('plasmidbench.html'));
  await page.waitForTimeout(600);

  const snap = () => page.evaluate(() => {
    const r = window.PB.active();
    return { seq: r.seq, name: r.name, topo: r.topology,
             feats: r.features.map(f => ({ l: f.label, s: f.start, e: f.end, st: f.strand, t: f.type })) };
  });
  const featSeqs = () => page.evaluate(() => {
    const r = window.PB.active(), B = window.Bio, o = {};
    r.features.forEach(f => { o[f.label] = B.featureSequence(r.seq, f); });
    return o;
  });
  const undoBtn = () => page.evaluate(() => ({
    u: !document.querySelector('#btnUndo').disabled, r: !document.querySelector('#btnRedo').disabled,
    ut: document.querySelector('#btnUndo').title, rt: document.querySelector('#btnRedo').title }));
  const clickModal = label => page.evaluate(l => {
    const b = Array.from(document.querySelectorAll('#modalRoot button')).find(x => x.textContent.trim() === l);
    if (!b) throw new Error('no button ' + l); b.click();
  }, label);
  const toView = v => page.evaluate(v => { const P = window.PB; P.state.view = v; P.renderTabs(); P.renderViews(); }, v);

  await page.evaluate(([seq, feats]) => {
    window.PB.addRecord(window.PB.makeRecord({ name: 'pTest', seq, topology: 'circular', features: feats }));
  }, [SEQ, FEATS]);
  await page.waitForTimeout(300);
  const s0 = await snap(), fs0 = await featSeqs();

  console.log('\n== nothing to undo at first ==');
  let ub = await undoBtn();
  ok('undo and redo start disabled', !ub.u && !ub.r, JSON.stringify(ub));

  console.log('\n== insert before base 1, inside a feature that spans it ==');
  await page.evaluate(() => { window.PB.state.selection = { start: 0, end: 1 }; window.PB.toolInsertDelete(); });
  await page.waitForTimeout(300);
  await page.selectOption('#idKind', 'insert');
  await page.fill('#idFrom', '1');
  await page.fill('#idSeq', 'GGGGGCCCCC');
  await page.waitForTimeout(150);
  const preview = await page.evaluate(() => document.querySelector('#idPreview').textContent);
  ok('the dialog previews the result', /10 bp go in before base 1/.test(preview) && /1,210 bp/.test(preview), preview);
  await clickModal('Apply'); await page.waitForTimeout(400);
  const s1 = await snap(), fs1 = await featSeqs();
  ok('the construct is edited in place, not copied', s1.name === 'pTest' && s1.seq.length === 1210,
     s1.name + ' ' + s1.seq.length);
  ok('the new bases are base 1', s1.seq.slice(0, 10) === 'GGGGGCCCCC');
  ok('the feature spanning base 1 takes the new bases',
     fs1.spanner === SEQ.slice(1150) + 'GGGGGCCCCC' + SEQ.slice(0, 40), fs1.spanner.length);
  ok('a feature clear of the edit keeps its exact bases', fs1.middle === fs0.middle && fs1.late === fs0.late);
  ub = await undoBtn();
  ok('undo is now offered, naming the edit', ub.u && /Inserted 10 bp before base 1/.test(ub.ut), ub.ut);

  console.log('\n== undo and redo ==');
  await page.click('#btnUndo'); await page.waitForTimeout(300);
  let s = await snap();
  ok('undo restores the sequence exactly', s.seq === s0.seq);
  ok('undo restores every feature exactly', JSON.stringify(s.feats) === JSON.stringify(s0.feats));
  ub = await undoBtn();
  ok('redo is offered after an undo', ub.r && /Inserted 10 bp/.test(ub.rt), ub.rt);
  await page.click('#btnRedo'); await page.waitForTimeout(300);
  s = await snap();
  ok('redo puts the edit back exactly', s.seq === s1.seq && JSON.stringify(s.feats) === JSON.stringify(s1.feats));

  console.log('\n== several kinds of change, undone in reverse ==');
  await page.fill('#recName', 'pTest renamed');
  await page.evaluate(() => document.querySelector('#recName').dispatchEvent(new Event('change')));
  await page.waitForTimeout(250);
  const s2 = await snap();
  await page.evaluate(() => { window.PB.state.selection = null; window.PB.setOriginDialog(); });
  await page.waitForTimeout(250);
  await page.fill('#soPos', '601');
  await clickModal('Rotate'); await page.waitForTimeout(350);
  const s3 = await snap();
  ok('the rotation happened', s3.seq !== s2.seq && s3.seq === s2.seq.slice(600) + s2.seq.slice(0, 600));
  await toView('seq'); await page.waitForTimeout(250);
  await page.evaluate(() => { window.PB.state.selection = { start: 100, end: 150 }; });
  await page.keyboard.press('Delete'); await page.waitForTimeout(350);
  const s4 = await snap();
  ok('Delete removes the selected bases in place', s4.seq.length === s3.seq.length - 50 &&
     s4.seq === s3.seq.slice(0, 100) + s3.seq.slice(150), s4.seq.length);
  /* a change made by any tool that saves through touch() is captured too */
  await page.evaluate(() => {
    const P = window.PB, r = P.active();
    r.features.push({ start: 10, end: 30, strand: 1, type: 'misc_feature', label: 'added', quals: {}, color: null, note: '' });
    P.addHistory(r, 'Feature added', 'added'); P.touch(r); P.renderAll();
  });
  const s5 = await snap();
  ok('five steps are on the stack', (await page.evaluate(() => window.PB.undoState().steps)) === 5,
     await page.evaluate(() => window.PB.undoState().steps));

  const back = [s4, s3, s2, s1, s0];
  for (let i = 0; i < back.length; i++) {
    await page.keyboard.press('Control+z'); await page.waitForTimeout(250);
    const now = await snap();
    ok('Ctrl+Z step ' + (i + 1) + ' returns exactly to the state before that change',
       now.seq === back[i].seq && now.name === back[i].name && JSON.stringify(now.feats) === JSON.stringify(back[i].feats),
       now.name + ' ' + now.seq.length);
  }
  ub = await undoBtn();
  ok('after undoing everything, undo is disabled again', !ub.u && ub.r);
  const fwd = [s1, s2, s3, s4, s5];
  for (let i = 0; i < fwd.length; i++) {
    await page.keyboard.press('Control+Shift+z'); await page.waitForTimeout(250);
    const now = await snap();
    ok('redo step ' + (i + 1) + ' reapplies that change exactly',
       now.seq === fwd[i].seq && now.name === fwd[i].name && JSON.stringify(now.feats) === JSON.stringify(fwd[i].feats));
  }

  console.log('\n== a new change clears the redo steps ==');
  await page.keyboard.press('Control+z'); await page.keyboard.press('Control+z'); await page.waitForTimeout(400);
  await page.evaluate(() => { window.PB.state.selection = { start: 5, end: 6 }; });
  await page.keyboard.press('Backspace'); await page.waitForTimeout(300);
  ub = await undoBtn();
  ok('redo is gone once a new edit is made', !ub.r && ub.u);

  console.log('\n== the history log keeps a record of undo and redo ==');
  const hist = await page.evaluate(() => window.PB.active().history.map(h => h.op));
  ok('undo and redo appear in the construct history', hist.indexOf('Undone') >= 0 && hist.indexOf('Redone') >= 0);

  console.log('\n== each construct has its own steps ==');
  const other = await page.evaluate(() => {
    const P = window.PB;
    const a = P.active().id;
    P.addRecord(P.makeRecord({ name: 'pOther', seq: 'ACGT'.repeat(100), topology: 'linear', features: [] }));
    return { a: a, b: P.active().id };
  });
  await page.waitForTimeout(300);
  ub = await undoBtn();
  ok('a different construct has nothing to undo', !ub.u, ub.ut);
  await page.keyboard.press('Control+z'); await page.waitForTimeout(250);
  const stillOther = await page.evaluate(() => window.PB.active().name);
  ok('undo there does not reach into the first construct', stillOther === 'pOther');
  await page.evaluate(id => { window.PB.state.activeId = id; window.PB.renderAll(); }, other.a);
  await page.waitForTimeout(250);
  ub = await undoBtn();
  ok('switching back, the first construct still has its steps', ub.u);

  console.log('\n== an edit across base 1 from the dialog ==');
  const before = await snap(), fsB = await featSeqs();
  const L = before.seq.length;
  await page.evaluate(() => { window.PB.state.selection = null; window.PB.toolInsertDelete(); });
  await page.waitForTimeout(250);
  await page.selectOption('#idKind', 'delete');
  await page.fill('#idFrom', String(L - 9));
  await page.fill('#idTo', '10');
  await page.waitForTimeout(150);
  const pv = await page.evaluate(() => document.querySelector('#idPreview').textContent);
  ok('the preview says the range crosses base 1', /across base 1/.test(pv) && /Removes 20 bp/.test(pv), pv);
  await clickModal('Apply'); await page.waitForTimeout(350);
  const after = await snap();
  ok('twenty bases are removed across the origin', after.seq === before.seq.slice(10, L - 10), after.seq.length);
  await page.keyboard.press('Control+z'); await page.waitForTimeout(300);
  const rest = await snap();
  ok('and undo restores them', rest.seq === before.seq && JSON.stringify(rest.feats) === JSON.stringify(before.feats));

  console.log('\n== pasting DNA opens the editor rather than changing anything ==');
  await toView('seq');
  await page.evaluate(() => { window.PB.state.selection = { start: 200, end: 230 }; });
  const lenBefore = (await snap()).seq.length;
  await page.evaluate(() => {
    const dt = new DataTransfer(); dt.setData('text/plain', '>frag\nAAAACCCCGGGGTTTT\n');
    document.body.dispatchEvent(new ClipboardEvent('paste', { clipboardData: dt, bubbles: true }));
  });
  await page.waitForTimeout(300);
  const dlg = await page.evaluate(() => ({ kind: document.querySelector('#idKind') && document.querySelector('#idKind').value,
    from: document.querySelector('#idFrom') && document.querySelector('#idFrom').value,
    to: document.querySelector('#idTo') && document.querySelector('#idTo').value,
    seq: document.querySelector('#idSeq') && document.querySelector('#idSeq').value }));
  ok('the editor opens on Replace with the selection and the pasted bases',
     dlg.kind === 'replace' && dlg.from === '201' && dlg.to === '230' && dlg.seq === 'AAAACCCCGGGGTTTT', JSON.stringify(dlg));
  ok('nothing changed before confirming', (await snap()).seq.length === lenBefore);
  await clickModal('Cancel'); await page.waitForTimeout(200);

  console.log('\n== Ctrl+Z in a text box is the browser\'s, not the app\'s ==');
  const stepsBefore = await page.evaluate(() => window.PB.undoState().steps);
  await page.focus('#recName');
  await page.keyboard.press('Control+z'); await page.waitForTimeout(200);
  ok('the construct is untouched by Ctrl+Z typed in a field',
     (await page.evaluate(() => window.PB.undoState().steps)) === stepsBefore);
  await page.evaluate(() => document.activeElement.blur());

  console.log('\n== the demo construct never changes ==');
  await page.evaluate(() => { const P = window.PB; P.state.activeId = P.state.records.find(r => r.builtin).id; P.renderAll(); });
  await page.waitForTimeout(250);
  const demo0 = await snap();
  await toView('seq');
  await page.evaluate(() => { window.PB.state.selection = { start: 10, end: 20 }; });
  await page.keyboard.press('Delete'); await page.waitForTimeout(250);
  ok('Delete on the demo leaves it as it was', (await snap()).seq === demo0.seq);
  await page.evaluate(() => window.PB.toolInsertDelete());
  await page.waitForTimeout(200);
  const forced = await page.evaluate(() => ({ c: document.querySelector('#idCopy').checked, d: document.querySelector('#idCopy').disabled }));
  ok('the editor forces a new construct for the demo', forced.c && forced.d);
  await clickModal('Cancel');

  ok('no errors anywhere in the run', errors.length === 0, errors.slice(0, 3).join(' | '));
  console.log('\n   ' + pass + ' passed, ' + fail + ' failed\n');
  await browser.close();
  process.exit(fail ? 1 : 0);
})();
