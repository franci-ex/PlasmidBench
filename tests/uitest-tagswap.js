/* ===========================================================================
   Swapping a solubility tag, through the interface, and then swapping the
   protein behind it on the product — the two-stage job.

   The thing being tested is not that primers come out, but that the tool
   refuses to let a frameshift through quietly: the map looks identical
   either way, and only the translation tells them apart.
   ========================================================================= */
const { chromium } = require('playwright');
const path = require('path');
let pass = 0, fail = 0;
function ok(name, cond, extra) {
  if (cond) { pass++; console.log('  PASS  ' + name); }
  else { fail++; console.log('  FAIL  ' + name + (extra !== undefined ? '   -> ' + extra : '')); }
}

/* a construct shaped like the real thing: tag - TEV - protein */
const BUILD = `
  const bt = aa => { const c = { A:'GCG',R:'CGT',N:'AAC',D:'GAT',C:'TGC',Q:'CAG',E:'GAA',
    G:'GGC',H:'CAT',I:'ATT',L:'CTG',K:'AAA',M:'ATG',F:'TTT',P:'CCG',S:'AGC',T:'ACC',
    W:'TGG',Y:'TAT',V:'GTG' }; return aa.split('').map(x => c[x]).join(''); };
  let s = 7777; const rnd = () => (s = (s * 1103515245 + 12345) & 0x7fffffff) / 0x7fffffff;
  const rs = n => { let o = ''; for (let i = 0; i < n; i++) o += 'ACGT'[Math.floor(rnd() * 4)]; return o; };
  const prot = n => { const a = 'ARNDCQEGHILKMFPSTWYV'; let p = 'M';
    for (let i = 1; i < n; i++) p += a[Math.floor(rnd() * 20)]; return bt(p); };
  const GB1 = bt('MTYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE');
  const TEV = bt('ENLYFQG');
  const POI = prot(300).slice(3);
  const CDS = GB1 + TEV + POI + 'TAA';
  const up = rs(300), down = rs(1400);
  const seq = up + CDS + down;
  const MBP = 'ATG' + prot(370).slice(3);
  const donor = rs(900) + MBP + rs(1100);
  window.__fx = { seq, up, GB1, TEV, POI, CDS, MBP, donor,
                  gb1a: up.length, gb1b: up.length + GB1.length,
                  cdsEnd: up.length + CDS.length };
`;

(async () => {
  const browser = await chromium.launch();
  const page = await browser.newPage({ viewport: { width: 1500, height: 1050 } });
  const errors = [];
  page.on('pageerror', e => errors.push('pageerror: ' + String(e)));
  page.on('console', m => { if (m.type() === 'error' && !/Failed to load resource/.test(m.text())) errors.push('console: ' + m.text()); });
  await page.goto('file://' + path.resolve('plasmidbench.html'));
  await page.waitForTimeout(700);

  await page.evaluate(BUILD);
  await page.evaluate(() => {
    const P = window.PB, f = window.__fx;
    P.addRecord(P.makeRecord({
      name: 'pGB1-TEV-POI', seq: f.seq, topology: 'circular',
      features: [{ start: f.up.length, end: f.cdsEnd, strand: 1, type: 'CDS',
                   label: 'GB1-TEV-POI', quals: { codon_start: '1' } },
                 { start: f.gb1a, end: f.gb1b, strand: 1, type: 'CDS', label: 'GB1 tag' }]
    }));
    P.addRecord(P.makeRecord({
      name: 'pDonor-MBP', seq: f.donor, topology: 'circular',
      features: [{ start: 900, end: 900 + f.MBP.length, strand: 1, type: 'CDS', label: 'MBP' }]
    }));
    P.renderAll();
  });
  await page.waitForTimeout(400);

  console.log('\n== stage 1: amplify the backbone without the tag, put a new tag in ==');
  /* the vector is the construct, opened by a typed range = the old tag;
     the insert is the new tag, pasted as it would be from the donor */
  const setup = async (from, to) => {
    await page.evaluate(() => { const P = window.PB; P.closeModal(); P.toolGibson ? 0 : 0; });
    await page.evaluate(() => {
      const P = window.PB;
      P.state.activeId = P.state.records.find(r => r.name === 'pGB1-TEV-POI').id;
      P.renderAll();
    });
    await page.evaluate(() => {
      const b = Array.from(document.querySelectorAll('#modalRoot [data-tool]'));
      window.PB.openTools();
    });
    await page.evaluate(() => {
      const btns = Array.from(document.querySelectorAll('#modalRoot [data-tool]'));
      const b = btns.find(x => /In-Fusion/.test(x.textContent));
      b.click();
    });
    await page.waitForTimeout(400);
    /* paste the new tag as the insert */
    await page.evaluate(() => {
      const t = document.querySelector('#gbIns_txt');
      t.value = window.__fx.MBP;
      t.dispatchEvent(new Event('input', { bubbles: true }));
    });
    await page.waitForTimeout(300);
    /* choose "a range I type" and give it the old tag's bounds */
    await page.evaluate(() => {
      const sel = document.querySelector('#gbWhere');
      const i = Array.from(sel.options).findIndex(o => /a range I type/.test(o.textContent));
      sel.value = String(Array.from(sel.options)[i].value);
      sel.dispatchEvent(new Event('change', { bubbles: true }));
    });
    await page.waitForTimeout(300);
    await page.evaluate(([a, b]) => {
      const f = document.querySelector('#gbFrom'), t = document.querySelector('#gbTo');
      f.value = String(a); f.dispatchEvent(new Event('input', { bubbles: true }));
      t.value = String(b); t.dispatchEvent(new Event('input', { bubbles: true }));
    }, [from, to]);
    await page.waitForTimeout(900);
  };

  const fx = await page.evaluate(() => window.__fx);
  await setup(fx.gb1a + 1, fx.gb1b);        // 1-based, inclusive: exactly the tag

  const good = await page.evaluate(() => {
    const live = document.querySelector('#gbLive') || document.querySelector('#gbDesign');
    const t = live ? live.textContent : '';
    return {
      inFrame: /In frame\./.test(t),
      broken: /This breaks the protein/.test(t),
      residues: (t.match(/from ([\d,]+) to ([\d,]+) residues/) || []).slice(1, 3),
      tev: /TEV site behind it survives/.test(t),
      join: /read as protein/.test(t)
    };
  });
  ok('a clean tag swap is reported as in frame', good.inFrame && !good.broken, JSON.stringify(good));
  ok('it says how the protein changes size', good.residues.join('->') === '363->677', JSON.stringify(good.residues));
  ok('and that the TEV site behind the tag survives', good.tev);
  ok('the junction is shown as protein', good.join);

  console.log('\n== the same swap, one base out ==');
  await setup(fx.gb1a + 1, fx.gb1b + 1);
  const bad = await page.evaluate(() => {
    const live = document.querySelector('#gbLive') || document.querySelector('#gbDesign');
    const t = live ? live.textContent : '';
    return { broken: /This breaks the protein/.test(t),
             frame: /wrong frame/.test(t),
             stop: /stop codon now falls/.test(t),
             tev: /TEV site \(ENLYFQ\) that sat behind the tag is gone/.test(t) };
  });
  ok('one base out is called out as breaking the protein', bad.broken, JSON.stringify(bad));
  ok('the frame is named as the reason', bad.frame, bad.frame);
  ok('the premature stop is reported', bad.stop, bad.stop);
  ok('and the lost TEV site is reported', bad.tev, bad.tev);
  ok('the two cases are told apart only by the translation, not the primers',
     good.inFrame && bad.broken);

  console.log('\n== stage 2: the product can be used to swap the protein ==');
  await setup(fx.gb1a + 1, fx.gb1b);
  const saved = await page.evaluate(() => {
    const b = document.querySelector('#gbSave');
    if (!b) return { saved: false };
    b.click();
    return { saved: true };
  });
  await page.waitForTimeout(800);
  const after = await page.evaluate(() => {
    const P = window.PB, f = window.__fx;
    const rec = P.state.records[P.state.records.length - 1];
    /* The construct is circular, and the swap can leave the new tag sitting
       across base 1 — perfectly legal, and how a plasmid really is. Read the
       sequence doubled, which is how any lookup on a circle has to be done;
       the old single-copy indexOf found the tag at the very end and then ran
       off the string looking for what follows it. */
    const d = rec.seq + rec.seq;
    const idx = d.indexOf(f.MBP);
    const tev = d.indexOf(f.TEV);
    return { n: P.state.records.length, name: rec.name, len: rec.seq.length,
             hasTag: idx >= 0, tagAt: idx,
             tevRightAfter: idx >= 0 && d.slice(idx + f.MBP.length, idx + f.MBP.length + f.TEV.length) === f.TEV,
             afterTag: idx >= 0 ? d.slice(idx + f.MBP.length, idx + f.MBP.length + 48) : '',
             tevAt: tev, gap: (tev >= 0 && idx >= 0) ? tev - (idx + f.MBP.length) : null,
             poiKept: d.indexOf(f.POI) >= 0 };
  });
  ok('the new construct is saved', saved.saved && after.n >= 3, JSON.stringify(after));
  ok('it carries the new tag', after.hasTag, after.tagAt);
  ok('with the TEV site immediately behind it', after.tevRightAfter,
     'gap=' + after.gap + ' after the tag: ' + after.afterTag);
  ok('and the protein still there to swap next', after.poiKept);

  /* now run the second swap on that product: replace the protein */
  const second = await page.evaluate(() => {
    const P = window.PB, f = window.__fx;
    const rec = P.state.records[P.state.records.length - 1];
    const start = rec.seq.indexOf(f.POI);
    const check = window.Bio.swapFrameCheck(rec.seq, start, start + f.POI.length,
      'ATG'.repeat(1) + f.POI.slice(3, 600), { circular: true, features: rec.features || [] });
    return { applies: check.applies, ok: check.ok, why: check.why,
             before: check.residuesBefore, after: check.residuesAfter };
  });
  ok('the product can be checked for a second swap on top', second.applies,
     second.why || JSON.stringify(second));
  ok('and that second swap is judged on its own frame', typeof second.ok === 'boolean',
     JSON.stringify(second));

  await page.evaluate(() => window.PB.closeModal());
  ok('nothing threw anywhere in the run', errors.length === 0, errors.slice(0, 3).join(' | '));
  console.log('\n   ' + pass + ' passed, ' + fail + ' failed\n');
  await browser.close();
  process.exit(fail ? 1 : 0);
})().catch(e => { console.error('HARNESS ERROR', e); process.exit(2); });
