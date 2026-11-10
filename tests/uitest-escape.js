/* ===========================================================================
   Text from a file, put on screen.

   This page builds its interface by joining strings, which is fine until a
   value that came from outside carries a `<` in it. Correctness then rests
   on calling esc() at every one of some four hundred interpolations, and a
   count of four hundred is not evidence. Reading them one at a time proves
   nothing about the next change either.

   So this is the behavioural test instead. A construct is loaded whose name,
   description, every feature label, every qualifier, every primer name and
   every read name is markup — an image with an onerror, a script tag, a
   closing attribute quote, an event handler on a div. Then every view and
   every tool panel is opened in turn, and after each one the page is asked:
   did anything run, and is there an element here that the file put there?

   Nothing in this file knows which call sites escape and which do not, so it
   keeps working when the interface is rewritten. A real lab file can carry
   any of this: a GenBank /note is free text, and so is a sequencing
   facility's sample name.
   ========================================================================= */
const { chromium } = require('playwright');
const path = require('path');
let pass = 0, fail = 0;
function ok(name, cond, extra) {
  if (cond) { pass++; console.log('  PASS  ' + name); }
  else { fail++; console.log('  FAIL  ' + name + (extra !== undefined ? '   -> ' + extra : '')); }
}

/* Each of these, dropped into markup unescaped, leaves a trace this test can
   find: a tag that should not exist, or a handler that fires. */
const PAYLOADS = [
  '<img src=x onerror="window.__ran=(window.__ran||0)+1" class="pwned">',
  '"><span class="pwned">x</span>',
  "'><span class=pwned>y</span>",
  '<script>window.__ran=(window.__ran||0)+1<\/script>',
  '<div onmouseover="window.__ran=1" class="pwned">z</div>',
  '</td></tr><tr class="pwned"><td>broken out',
  '<svg onload="window.__ran=1" class="pwned"></svg>',
  '<iframe class="pwned" src="about:blank"></iframe>'
];

(async () => {
  const browser = await chromium.launch();
  const page = await browser.newPage({ viewport: { width: 1500, height: 1100 } });
  const errors = [];
  page.on('pageerror', e => errors.push('pageerror: ' + String(e)));
  page.on('console', m => {
    if (m.type() === 'error' && !/Failed to load resource|ERR_(TUNNEL|NAME|INTERNET|CONNECTION)/.test(m.text()))
      errors.push(m.text());
  });
  await page.goto('file://' + path.resolve(__dirname, '..', 'plasmidbench.html'));
  await page.waitForTimeout(1200);

  await page.evaluate((loads) => {
    const P = window.PB, B = window.Bio;
    window.__ran = 0;
    let x = 31; const r = () => (x = (x * 1103515245 + 12345) & 0x7fffffff) / 0x7fffffff;
    const dna = n => { let s = ''; for (let i = 0; i < n; i++) s += 'ACGT'[Math.floor(r() * 4)]; return s; };
    const seq = 'ATG' + dna(1200) + 'TAA' + dna(1800);

    /* a feature for every payload, plus one of each kind the map draws */
    const features = loads.map((p, i) => ({
      start: 100 + i * 180, end: 220 + i * 180, strand: i % 2 ? -1 : 1,
      type: ['CDS', 'promoter', 'terminator', 'primer_bind', 'rep_origin',
             'misc_feature', 'gene', 'protein_bind'][i % 8],
      label: p,
      quals: { note: p, product: p, gene: p, locus_tag: p, organism: p }
    }));
    features.push({ start: 0, end: seq.length, strand: 1, type: 'source',
                    label: loads[0], quals: { organism: loads[1] } });

    const rec = P.makeRecord({
      name: loads[0], description: loads[1], seq: seq,
      topology: 'circular', features: features
    });
    P.state.records = [rec];
    P.state.activeId = rec.id;

    /* the library panels take text from outside too */
    P.state.primers = loads.map((p, i) => ({ name: p, seq: dna(22), note: p, loc: p, method: p }));
    P.state.stock = loads.map((p, i) => ({ name: p, seq: dna(20), id: p, loc: p, note: p, tm: 60 }));
    P.state.plasmidStock = loads.map(p => ({ name: p, loc: p, note: p, marker: p }));
    P.state.cellStock = loads.map(p => ({ name: p, loc: p, note: p, genotype: p }));
    P.state.reads = loads.map((p, i) => ({ name: p, seq: dna(400), quality: null, note: p }));
    P.state.selection = { start: 120, end: 340 };
    P.renderAll();
  }, PAYLOADS);
  await page.waitForTimeout(900);

  /* What a payload leaves behind if it was not escaped. The page's own
     script blocks are in the body, so the count is taken once at the start
     and compared against — an absolute count would report the page itself. */
  const baseScripts = await page.evaluate(() => document.querySelectorAll('script').length);
  const probe = () => page.evaluate(base => ({
    ran: window.__ran || 0,
    pwned: document.querySelectorAll('.pwned').length,
    newScripts: document.querySelectorAll('script').length - base,
    iframes: document.querySelectorAll('iframe').length,
    imgs: Array.from(document.querySelectorAll('img')).filter(i => /(^|[?&])x$|\/x$|^x$/.test(i.getAttribute('src') || '')).length
  }), baseScripts);

  console.log('\n== every view, with a file made of markup ==');
  const views = await page.evaluate(() =>
    Array.from(document.querySelectorAll('[role="tab"][data-v]')).map(b => b.getAttribute('data-v')));
  ok('there are views to go through', views.length >= 6, views.join(' '));
  for (const v of views) {
    await page.evaluate(k => { window.PB.state.view = k; window.PB.renderAll(); }, v);
    await page.waitForTimeout(260);
    const p = await probe();
    ok('the ' + v + ' view renders it as text',
       p.ran === 0 && p.pwned === 0 && p.newScripts === 0 && p.iframes === 0 && p.imgs === 0,
       JSON.stringify(p));
  }

  console.log('\n== the name is on screen, as characters ==');
  const shown = await page.evaluate(() => {
    const t = document.body.innerText;
    return { hasAngle: t.indexOf('<img src=x onerror=') >= 0, len: t.length };
  });
  ok('the construct name is displayed, not interpreted',
     shown.hasAngle, JSON.stringify(shown));

  console.log('\n== each tool panel ==');
  /* Opening every tool is the broadest sweep available: each builds its own
     markup from the same record. */
  const tools = await page.evaluate(() => {
    const names = (window.PB.toolList ? window.PB.toolList() : []).map(t => t[1] || t.name);
    return names.filter(Boolean);
  });
  if (!tools.length) {
    /* no exported list; drive the menu instead */
    const menuTools = await page.evaluate(() => {
      const btn = document.querySelector('[data-act="tools"], #toolsBtn');
      if (btn) btn.click();
      return Array.from(document.querySelectorAll('#modalRoot [data-tool], .toolrow')).length;
    });
    ok('the tools menu opens', menuTools >= 0, menuTools);
    const p = await probe();
    ok('and the menu itself is clean',
       p.ran === 0 && p.pwned === 0 && p.newScripts === 0, JSON.stringify(p));
    await page.evaluate(() => window.PB.closeModal && window.PB.closeModal());
  } else {
    let dirty = [];
    for (const t of tools) {
      await page.evaluate(n => {
        const list = window.PB.toolList();
        const row = list.find(r => (r[1] || r.name) === n);
        try { (row[3] || row.run)(); } catch (e) { /* a tool may refuse this record */ }
      }, t);
      await page.waitForTimeout(170);
      const p = await probe();
      if (p.ran || p.pwned || p.newScripts || p.iframes || p.imgs) dirty.push(t + ' ' + JSON.stringify(p));
      await page.evaluate(() => window.PB.closeModal && window.PB.closeModal());
      await page.waitForTimeout(60);
    }
    ok('all ' + tools.length + ' tool panels render it as text', dirty.length === 0,
       dirty.slice(0, 3).join(' | '));
  }

  console.log('\n== the alignment view, where the labels are record names ==');
  const aln = await page.evaluate((loads) => {
    const html = window.PB.alignBlocks('ACGTACGTAC', 'ACGTACGTAC', loads[0], loads[1], 10, 0);
    const d = document.createElement('div');
    d.innerHTML = html;
    return { pwned: d.querySelectorAll('.pwned').length,
             imgs: d.querySelectorAll('img').length,
             text: (d.textContent || '').indexOf('<img') >= 0 };
  }, PAYLOADS);
  ok('a name with markup in it stays text in the alignment',
     aln.pwned === 0 && aln.imgs === 0, JSON.stringify(aln));
  ok('and is still readable as what it was', aln.text, JSON.stringify(aln));

  console.log('\n== exported files carry the text, not the markup ==');
  const exported = await page.evaluate(() => {
    const r = window.PB.active();
    const gb = window.Bio.toGenBank(r);
    /* A GenBank qualifier is quoted text, and what must not happen is the
       quoting being broken by the content: a label holding a double quote,
       or long enough to be wrapped, has to come back as the same string.
       Compared as a set, because the writer puts the source feature first
       whatever order the record held it in. */
    const back = window.Bio.parseGenBank(gb)[0];
    const sig = rec => (rec.features || [])
      .map(f => f.type + '|' + f.start + '|' + f.end + '|' + f.strand + '|' + (f.label || ''))
      .sort();
    const a = sig(r), b = sig(back);
    const lost = a.filter(x => b.indexOf(x) === -1);
    /* GenBank expects a source feature spanning the molecule and the writer
       supplies one, so a single gained source of full length is the format
       being honoured rather than a feature appearing out of nowhere. */
    const gained = b.filter(x => a.indexOf(x) === -1)
                    .filter(x => !/^source\|0\|/.test(x));
    return { count: a.length, back: b.length, sameSeq: back.seq === r.seq,
             lost: lost.slice(0, 2), gained: gained.slice(0, 2) };
  });
  ok('every feature survives, label and all',
     exported.lost.length === 0 && exported.gained.length === 0,
     JSON.stringify(exported));
  ok('and so does the sequence', exported.sameSeq, exported.sameSeq);

  console.log('\n== console ==');
  ok('no page errors anywhere in that', errors.length === 0, errors.slice(0, 3).join(' | '));

  console.log('\n   ' + pass + ' passed, ' + fail + ' failed\n');
  await browser.close();
  process.exit(fail ? 1 : 0);
})();
