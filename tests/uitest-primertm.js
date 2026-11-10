/* ===========================================================================
   The melting temperature shown for a stored primer.

   A primer with a 5' overhang anneals only along its 3' part in the first
   cycles, and that part alone sets the temperature the block has to run at.
   The whole primer's Tm is always higher — on a real assembly primer by
   fifteen degrees — and a block set to it anneals above what holds, which
   gives nothing at all. So the number the library shows has to be the
   annealing one wherever it can be known.

   The hard case, and the one this suite is mostly about, is a primer saved
   before the program recorded any of this: a name and a string of bases.
   Nothing in the string says where the arm ends. What settles it is that an
   assembly arm is vector sequence, so it is found in the vector the primer
   was built for — the construct most likely to be open.
   ========================================================================= */
const { chromium } = require('playwright');
const path = require('path');
let pass = 0, fail = 0;
function ok(name, cond, extra) {
  if (cond) { pass++; console.log('  PASS  ' + name); }
  else { fail++; console.log('  FAIL  ' + name + (extra !== undefined ? '   -> ' + extra : '')); }
}

const ARM_F = 'TTGTATTTCCAGGGC';          // 15 nt
const ARM_R = 'ATTATGCGGCCGCAATT';        // 17 nt
const INS_F = 'ATGGGTTTTCTGCGTTTTC';

(async () => {
  const browser = await chromium.launch();
  const page = await browser.newPage({ viewport: { width: 1500, height: 1150 } });
  const errors = [];
  page.on('pageerror', e => errors.push('pageerror: ' + String(e)));
  page.on('console', m => { if (m.type() === 'error' && !/Failed to load resource|ERR_(TUNNEL|NAME|INTERNET|CONNECTION)/.test(m.text())) errors.push(m.text()); });
  await page.goto('file://' + path.resolve(__dirname, '..', 'plasmidbench.html'));
  await page.waitForTimeout(800);

  /* a vector carrying both arms, and primers stored the old way: name and
     sequence, nothing else */
  await page.evaluate(([aF, aR, iF]) => {
    const P = window.PB, B = window.Bio;
    let x = 99; const r = () => (x = (x * 1103515245 + 12345) & 0x7fffffff) / 0x7fffffff;
    const dna = n => { let s = ''; for (let i = 0; i < n; i++) s += 'ACGT'[Math.floor(r() * 4)]; return s; };
    const vec = dna(900) + B.revcomp(aF) + dna(60) + B.revcomp(aR) + dna(900);
    const rec = P.makeRecord({ name: 'FIX-0528 pRSFDuet', seq: vec, topology: 'circular', features: [] });
    P.state.records = [rec]; P.state.activeId = rec.id;
    window.__vec = vec;
    P.state.primers = [
      { name: 'SPDYC_F', seq: aF + iF },
      { name: 'SPDYC_R', seq: aR + 'AAACCTGAGGGCATACAC' },
      { name: 'SPDYC_bb_F', seq: vec.slice(300, 320) },
      { name: 'Ordered_elsewhere', seq: 'GCTTACGGTACCTATTGGCAGTTAC' }
    ];
    P.state.overhangCache = null;
    P.state.view = 'primers'; P.renderAll();
  }, [ARM_F, ARM_R, INS_F]);
  await page.waitForTimeout(700);

  const readLib = () => page.evaluate(() =>
    Array.from(document.querySelectorAll('#view-primers table.tbl tbody tr')).map(tr => {
      const c = Array.from(tr.children), btn = c[3].querySelector('button');
      const lines = btn ? btn.innerText.split('\n').map(x => x.trim()).filter(Boolean) : [c[3].textContent.trim()];
      return { name: (tr.querySelector('input') || {}).value, lines: lines,
               title: c[3].getAttribute('title') || '' };
    }));

  console.log('\n== a primer saved before any of this was recorded ==');
  const lib = await readLib();
  const F = lib.find(r => r.name === 'SPDYC_F');
  ok('the library reports a temperature for it', F && F.lines.length, JSON.stringify(F));
  ok('and it is the annealing one, not the whole primer’s',
     F && Math.abs(+F.lines[0] - 58.6) < 0.2, F && F.lines[0]);
  ok('the whole primer’s is still shown, underneath',
     F && /^full 7\d\.\d$/.test(F.lines[1] || ''), F && F.lines[1]);
  ok('the two are far enough apart to matter',
     F && Math.abs(+F.lines[0] - +(F.lines[1] || '').replace('full ', '')) > 10,
     F && F.lines.join(' / '));
  ok('it says the arm was found rather than presenting it as given',
     F && /nt arm, found/.test(F.lines[2] || ''), F && F.lines[2]);
  ok('and the 15 bases it found are the arm',
     F && /^15 nt arm/.test(F.lines[2] || ''), F && F.lines[2]);
  ok('hovering explains which number is which and where the split came from',
     F && /3′ 19 bases are what binds/.test(F.title) &&
          /FIX-0528 pRSFDuet/.test(F.title), F && F.title.slice(0, 120));
  ok('and warns that the boundary errs low, never high',
     F && /lower temperature, never a higher one/.test(F.title), F && F.title.slice(-140));

  console.log('\n== primers that have no overhang are left alone ==');
  const bb = lib.find(r => r.name === 'SPDYC_bb_F');
  ok('a primer that is simply vector sequence gets one number',
     bb && bb.lines.length === 1, bb && JSON.stringify(bb.lines));
  const el = lib.find(r => r.name === 'Ordered_elsewhere');
  ok('and so does one that matches nothing loaded',
     el && el.lines.length === 1, el && JSON.stringify(el.lines));
  ok('whose tooltip says no overhang is known and offers to be told',
     el && /No overhang is known/.test(el.title) && /click to set one/i.test(el.title),
     el && el.title.slice(0, 110));

  console.log('\n== setting the overhang by hand ==');
  await page.evaluate(() => {
    const rows = Array.from(document.querySelectorAll('#view-primers table.tbl tbody tr'));
    const i = rows.findIndex(tr => (tr.querySelector('input') || {}).value === 'Ordered_elsewhere');
    rows[i].querySelectorAll('button')[0].click();
  });
  await page.waitForTimeout(400);
  const dlg = await page.evaluate(() => {
    const m = document.querySelector('#modalRoot');
    return { open: !!m && !m.hidden, title: (m.querySelector('h3') || {}).textContent || '',
             hasField: !!m.querySelector('#poTail'), out: (m.querySelector('#poOut') || {}).textContent || '' };
  });
  ok('the Tm opens a dialog for it', dlg.open && dlg.hasField, JSON.stringify(dlg));
  ok('named after the primer', /Ordered_elsewhere/.test(dlg.title), dlg.title);
  ok('which starts by reporting the whole primer', /25 nt anneals/.test(dlg.out), dlg.out);
  const live = await page.evaluate(() => {
    const t = document.querySelector('#poTail');
    t.value = '7'; t.dispatchEvent(new Event('input', { bubbles: true }));
    return { out: document.querySelector('#poOut').textContent,
             coloured: !!document.querySelector('#poSeq span') };
  });
  ok('typing an overhang updates what it leaves, before saving',
     /18 nt anneals/.test(live.out), live.out);
  ok('and marks the overhang in the sequence', live.coloured);
  await page.evaluate(() => {
    const b = Array.from(document.querySelectorAll('#modalRoot button')).find(x => x.textContent.trim() === 'Save');
    b.click();
  });
  await page.waitForTimeout(500);
  const after = await readLib();
  const el2 = after.find(r => r.name === 'Ordered_elsewhere');
  ok('saving it changes the number the library shows',
     el2 && el2.lines.length >= 2 && el2.lines[0] !== el.lines[0],
     el2 && JSON.stringify(el2.lines));
  ok('and says the split was set by you', el2 && /set by you/.test(el2.lines[2] || ''),
     el2 && JSON.stringify(el2.lines));
  ok('it survives a re-render', await page.evaluate(() => {
    window.PB.renderAll();
    return window.PB.state.primers.find(p => p.name === 'Ordered_elsewhere').annealLen === 18;
  }));

  console.log('\n== and it can be told there is none ==');
  await page.waitForTimeout(300);
  await page.evaluate(() => {
    const rows = Array.from(document.querySelectorAll('#view-primers table.tbl tbody tr'));
    const i = rows.findIndex(tr => (tr.querySelector('input') || {}).value === 'SPDYC_F');
    rows[i].querySelectorAll('button')[0].click();
  });
  await page.waitForTimeout(350);
  await page.evaluate(() => {
    const t = document.querySelector('#poTail');
    t.value = '0'; t.dispatchEvent(new Event('input', { bubbles: true }));
    Array.from(document.querySelectorAll('#modalRoot button')).find(x => x.textContent.trim() === 'Save').click();
  });
  await page.waitForTimeout(500);
  const after2 = await readLib();
  const F2 = after2.find(r => r.name === 'SPDYC_F');
  ok('setting it to zero reports the whole primer again',
     F2 && F2.lines.length === 1 && Math.abs(+F2.lines[0] - 72.1) < 0.3, F2 && JSON.stringify(F2.lines));
  ok('and the arm search does not quietly put it back',
     await page.evaluate(() => {
       window.PB.state.overhangCache = null; window.PB.renderAll();
       const p = window.PB.state.primers.find(x => x.name === 'SPDYC_F');
       return p.noOverhang === true && !p.annealLen;
     }));

  console.log('\n== console ==');
  ok('no page errors anywhere in that', errors.length === 0, errors.slice(0, 3).join(' | '));

  console.log('\n   ' + pass + ' passed, ' + fail + ' failed\n');
  await browser.close();
  process.exit(fail ? 1 : 0);
})();
