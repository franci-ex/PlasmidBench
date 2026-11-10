/* ===========================================================================
   The cloning-site picker, judged by what it builds.

   The sites offered come from the vector's own annotation, and each one has
   to land where its label promises. This suite used to read the position out
   of the dialog's own From/To fields, which only proved the dialog agreed
   with itself. The placement is checked through the finished construct
   instead: pick "after 6xHis, in frame" and the tag and the insert must
   translate as one open reading frame; pick "from the ATG, untagged" and the
   tag must be gone while the ribosome binding site and its spacer survive
   untouched. Those are the two things a cloner is actually trusting.

   (The dialog's controls were renamed when the placement list replaced the
   old linearise-then-assemble flow: #gbLinSite became #gbWhere, #gbLinFrom
   and #gbLinTo became #gbFrom and #gbTo, and #gbLinGo and #gbGo went away
   because nothing waits for a button any more.)
   ========================================================================= */
const { chromium } = require('playwright');
const path = require('path'), fs = require('fs');
const B = require('../src/engine.js'); require('../src/engine2.js');
let pass = 0, fail = 0;
function ok(name, cond, extra) {
  if (cond) { pass++; console.log('  PASS  ' + name); }
  else { fail++; console.log('  FAIL  ' + name + (extra !== undefined ? '   -> ' + extra : '')); }
}
function prng(seed) { let x = seed; return () => (x = (x * 1103515245 + 12345) & 0x7fffffff) / 0x7fffffff; }
function randDNA(n, seed) { const r = prng(seed); let s = ''; for (let i = 0; i < n; i++) s += 'ACGT'[Math.floor(r() * 4)]; return s; }
function wrap(body) {
  return '<!doctype html><html><head><meta charset="utf8"><meta name="viewport" content="width=device-width,initial-scale=1">' +
    '<style>:root{color-scheme:light}body{margin:0;font:14px system-ui;background:#fafaf9}[hidden]{display:none!important}</style>' +
    '</head><body>' + body + '</body></html>';
}

/* A miniature expression vector, annotated the way a real one is:
   promoter, RBS, then an MCS holding ATG + a 6xHis tag + a polylinker. */
const PRO = 'TAATACGACTCACTATAGG';
const RBS = 'AAGGAG';
const LEADER = 'AACTTTAAT' + RBS + 'ATATACC';          // RBS with a 7 nt spacer
const TAG = 'ATGGGCAGCAGC' + 'CATCACCATCATCACCAC' + 'AGCCAGGATCCG';  // ATG GSS His6 SQDP
const POLY = 'AATTCGAGCTCGGCGCGCCTGCAGGTCGACAAGCTTGCGGCCGCATAA';
const BACK = randDNA(1800, 31);
const SEQ = PRO + LEADER + TAG + POLY + BACK;
// 1-based coordinates of the annotated parts
const p = { pro: 1, proEnd: PRO.length };
p.rbs = PRO.length + 10; p.rbsEnd = p.rbs + 5;
p.mcs = PRO.length + LEADER.length + 1;
p.atg = p.mcs; p.atgEnd = p.atg + 2;
p.his = p.mcs + 12; p.hisEnd = p.his + 17;
p.mcsEnd = p.mcs + TAG.length + POLY.length - 1;
const GB = (() => {
  let s = 'LOCUS       pTest  ' + SEQ.length + ' bp    DNA     circular SYN 01-JAN-2026\n' +
    'FEATURES             Location/Qualifiers\n' +
    '     promoter        ' + p.pro + '..' + p.proEnd + '\n                     /label="T7 promoter"\n' +
    '     RBS             ' + p.rbs + '..' + p.rbsEnd + '\n                     /label="RBS"\n' +
    '     misc_feature    ' + p.mcs + '..' + p.mcsEnd + '\n                     /label="MCS-1"\n' +
    '     CDS             ' + p.atg + '..' + p.atgEnd + '\n                     /label="ATG"\n                     /codon_start=1\n' +
    '     CDS             ' + p.his + '..' + p.hisEnd + '\n                     /label="6xHis"\n                     /codon_start=1\n' +
    'ORIGIN\n';
  for (let i = 0; i < SEQ.length; i += 60) s += String(i + 1).padStart(9) + ' ' + SEQ.slice(i, i + 60).match(/.{1,10}/g).join(' ') + '\n';
  return s + '//\n';
})();
// a clean ORF to clone in
const INS = 'ATG' + randDNA(297, 77).replace(/TAA|TAG|TGA/g, 'CCC') + 'TAA';

(async () => {
  fs.writeFileSync(path.join(__dirname, 'preview.html'),
    wrap(fs.readFileSync(path.join(__dirname, '..', 'plasmidbench.html'), 'utf8')));
  const browser = await chromium.launch();
  const page = await browser.newPage({ viewport: { width: 1480, height: 1100 } });
  const errors = [];
  const ignorable = t => /ERR_TUNNEL|ERR_NAME|fonts\.(googleapis|gstatic)|Failed to load resource/.test(t);
  page.on('console', m => { if (m.type() === 'error' && !ignorable(m.text())) errors.push(m.text()); });
  page.on('pageerror', e => errors.push('PAGEERROR: ' + e.message));
  await page.goto('file://' + path.join(__dirname, 'preview.html'));
  await page.waitForTimeout(700);
  ok('the test vector imports', await page.evaluate(t => window.PB.importText(t, 'test'), GB) === 1);

  /* The option VALUES are indices into the placement list; their positions in
     the list are not, because the options are grouped into optgroups in a
     reading order. Always pick by value. */
  const pickWhere = async (re) => {
    const v = await page.evaluate(src => {
      const sel = document.querySelector('#gbWhere');
      const o = Array.from(sel.options).find(x => new RegExp(src).test(x.textContent));
      if (!o) return null;
      sel.value = o.value;
      sel.dispatchEvent(new Event('change', { bubbles: true }));
      return { value: o.value, label: o.textContent };
    }, re.source);
    if (!v) throw new Error('no placement matching ' + re);
    await page.waitForTimeout(500);
    return v;
  };

  /* Saving a construct opens it, which is right — you want to look at what
     you just built. Designing a second clone from the same vector therefore
     means going back to the vector first, exactly as a person would. */
  async function openTool() {
    await page.evaluate(() => { const m = document.querySelector('#modalRoot'); if (m && !m.hidden) window.PB.closeModal(); });
    await page.evaluate(() => {
      const P = window.PB;
      const vec = P.state.records.find(r => /pTest/.test(r.name) && !/MyPOI/.test(r.name));
      if (vec) { P.state.activeId = vec.id; P.renderAll(); }
    });
    await page.waitForTimeout(300);
    await page.evaluate(() => window.PB.openTools()); await page.waitForTimeout(250);
    await page.locator('.tool-card', { hasText: 'Seamless assembly' }).first().click();
    await page.waitForTimeout(500);
    await page.selectOption('#gbIns_src', '');
    await page.fill('#gbIns_txt', INS);
    await page.fill('#gbInsName', 'MyPOI');
    await page.waitForTimeout(600);
  }

  console.log('\n== the sites offered come from the vector, not from guesswork ==');
  await openTool();
  const opts = await page.$$eval('#gbWhere option', os => os.map(o => o.textContent));
  ok('a site is offered for the annotated MCS', opts.some(o => /MCS-1/.test(o)), opts.join(' | '));
  ok('one option fuses in frame after the tag', opts.some(o => /in frame/.test(o)),
     opts.filter(o => /frame/.test(o)).join(' | '));
  ok('one option replaces from the ATG', opts.some(o => /from the ATG/.test(o)),
     opts.filter(o => /ATG/.test(o)).join(' | '));
  ok('one option replaces the whole site', opts.some(o => /whole site/.test(o)),
     opts.filter(o => /whole site/.test(o)).join(' | '));
  ok('a position of your own is still offered',
     opts.some(o => /a single base I type/.test(o)) && opts.some(o => /a range I type/.test(o)),
     opts.filter(o => /I type/.test(o)).join(' | '));
  ok('and the tag the fusion is behind is named in the label',
     opts.some(o => /after 6xHis, in frame/.test(o)), opts.filter(o => /in frame/.test(o)).join(' | '));

  console.log('\n== a position you type goes where you put it ==');
  /* This is the one placement whose coordinates the dialog does show, so the
     From field is the right thing to read for it. */
  await pickWhere(/a single base I type/);
  ok('the From field appears for it',
     await page.evaluate(() => !document.querySelector('#gbFromWrap').hidden));
  await page.fill('#gbFrom', String(p.hisEnd + 1));
  await page.waitForTimeout(500);
  ok('and nothing is removed by a single-base opening',
     await page.evaluate(() => document.querySelector('#gbToWrap').hidden));

  console.log('\n== the in-frame fusion really is a fusion ==');
  await openTool();
  await pickWhere(/after 6xHis, in frame/);
  ok('it assembles with no further clicks', await page.locator('#gbSave').count() === 1);
  await page.click('#gbSave'); await page.waitForTimeout(800);
  const r1 = await page.evaluate(() => {
    const x = window.PB.active();
    return { name: x.name, seq: x.seq, len: x.seq.length,
             poi: x.features.filter(f => f.label === 'MyPOI').map(f => ({ t: f.type, s: f.start, e: f.end })),
             his: x.features.some(f => f.label === '6xHis') };
  });
  ok('the construct is named after the insert', /MyPOI/.test(r1.name), r1.name);
  ok('the insert is annotated once, as a CDS', r1.poi.length === 1 && r1.poi[0].t === 'CDS', JSON.stringify(r1.poi));
  ok('the annotation covers exactly the insert',
     r1.poi.length === 1 && B.featureSequence(r1.seq, { start: r1.poi[0].s, end: r1.poi[0].e, strand: 1 }) === INS,
     r1.poi.length === 1 ? B.featureSequence(r1.seq, { start: r1.poi[0].s, end: r1.poi[0].e, strand: 1 }).slice(0, 24) : '');
  ok('the 6xHis tag is still annotated', r1.his);
  ok('the length is vector + insert', r1.len === SEQ.length + INS.length, r1.len + ' vs ' + (SEQ.length + INS.length));
  ok('the insert lands immediately behind the tag',
     (r1.seq + r1.seq).indexOf('CATCACCATCATCACCAC') >= 0 &&
     (r1.seq + r1.seq).indexOf(INS) > (r1.seq + r1.seq).indexOf('CATCACCATCATCACCAC'),
     (r1.seq + r1.seq).indexOf(INS) + ' vs ' + (r1.seq + r1.seq).indexOf('CATCACCATCATCACCAC'));
  /* The real test of "in frame": read from the vector's own ATG straight
     through into the insert and the protein must be the tag's residues
     followed by the insert's, with no stop in between and no frameshift. */
  const fusion = (() => {
    /* Start from the vector's own ATG. The whole TAG string is no longer
       contiguous in the product and never should be — the insert goes in
       between the His6 run and the polylinker that followed it, which is
       exactly what "after 6xHis" means. */
    const d = r1.seq + r1.seq + r1.seq, i = d.indexOf('ATGGGCAGCAGCCATCACCATCATCACCAC');
    if (i < 0) return '(vector ATG and tag not found)';
    const pr = B.translate(d.slice(i, i + 3000), { table: 1 });
    const stop = pr.indexOf('*');
    return stop < 0 ? pr : pr.slice(0, stop);
  })();
  const insAA = B.translate(INS, { table: 1 }).replace(/\*$/, '');
  ok('the tag and the insert translate as one open reading frame with no stop between them',
     /^MGSSHHHHHH/.test(fusion) && fusion.endsWith(insAA),
     fusion.slice(0, 32) + ' … ends with insert: ' + fusion.endsWith(insAA));
  ok('and the insert is translated whole, not clipped',
     fusion.length >= insAA.length + 6, fusion.length + ' vs insert ' + insAA.length);

  console.log('\n== the untagged option keeps the ribosome binding site intact ==');
  await openTool();
  const whereNow = await page.$$eval('#gbWhere option', os => os.map(o => o.textContent));
  ok('the untagged placement is still offered the second time round',
     whereNow.some(o => /from the ATG/.test(o)), whereNow.join(' | ').slice(0, 300));
  await pickWhere(/from the ATG/);
  await page.click('#gbSave'); await page.waitForTimeout(800);
  const r2 = await page.evaluate(() => { const x = window.PB.active(); return { seq: x.seq, len: x.seq.length }; });
  const dd = r2.seq + r2.seq;
  const j = dd.indexOf(INS);
  ok('the insert is in the saved construct', j >= 0, j);
  ok('the insert ATG sits straight after the native spacer',
     j >= 0 && dd.slice(j - LEADER.length, j) === LEADER,
     j >= 0 ? dd.slice(j - LEADER.length, j) + ' vs ' + LEADER : '');
  ok('so the RBS and its spacing are untouched',
     j >= 0 && dd.slice(j - LEADER.length, j).indexOf(RBS) >= 0,
     j >= 0 ? dd.slice(j - LEADER.length, j) : '');
  ok('the tag is gone', dd.indexOf('CATCACCATCATCACCAC') === -1);
  ok('no errors anywhere in the run', errors.length === 0, errors.slice(0, 2).join(' | '));

  console.log('\n==========================================');
  console.log('   ' + pass + ' passed, ' + fail + ' failed');
  console.log('==========================================\n');
  await browser.close();
  process.exit(fail ? 1 : 0);
})();
