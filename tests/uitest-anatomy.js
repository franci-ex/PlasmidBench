/* ===========================================================================
   Reading a fragment, through the interface.

   The thing being checked is the question a person actually arrives with:
   "I pasted my sequencing read — which part of this do I put into my
   vector?" A read off a sequencing primer starts inside the backbone, so
   the honest answer needs two sequences compared, not one described. The
   fixture is a pET-shaped vector and a read of a GB1–TEV–CDK1 construct
   made from it, so every boundary below is known by construction.
   ========================================================================= */
const { chromium } = require('playwright');
const path = require('path');
let pass = 0, fail = 0;
function ok(name, cond, extra) {
  if (cond) { pass++; console.log('  PASS  ' + name); }
  else { fail++; console.log('  FAIL  ' + name + (extra !== undefined ? '   -> ' + extra : '')); }
}
function prng(seed) { let x = seed; return () => (x = (x * 1103515245 + 12345) & 0x7fffffff) / 0x7fffffff; }
function dnaOf(n, seed) { const r = prng(seed); let s = ''; for (let i = 0; i < n; i++) s += 'ACGT'[Math.floor(r() * 4)]; return s; }

const HEAD = 'TAATACGACTCACTATAGGG' + 'GATTGTGAGCGGATAACAATT' +
             'CCCCTGTAGAAATAATTTTGTTTAACTTTAATAAGGAGATATACC';
const TAG = 'ATGGGCAGCAGCCATCACCATCATCACCACAGCCAGGATCCGAATTCGAGCTCGAGT';
const VEC = HEAD + TAG + dnaOf(2400, 31);
const GB1 = 'GGCACCTACAAACTGATCCTGAACGGTAAAACCCTGAAAGGTGAAACCACCACCGAAGCTGTAGACGCTGCTACTGCTGAAAAAGTT' +
            'TTCAAACAGTACGCTAACGACAACGGTGTGGACGGTGAATGGACCTACGACGACGCTACCAAAACCTTCACGGTTACGGAA';
const LINK = 'AGCGGTTCGGAT';
const TEV = 'GAAAACTTGTATTTCCAGGGC';
const CDK = 'GGGATGGAGAACTTCCAAAAGGTGGAAAAGATCGGAGAGGGCACGTACGGAGTTGTGTACAAAGCCAGAAACAAGTTGACGGGAGAG' +
            'GTGGTGGCGCTTAAGAAAATCCGCCTGGACACTGAGACTGAGGGTGTGCCCAGTACTGCCATCCGAGAGATCTCTCTGCTTAAGGAG' +
            'CTTAACCATCCTAATATTGTCAAGCTGCTGGATGTCATTCACACAGAAAATAAACTCTACCTGGTTTTTGAATTTCTGCACCAAGAT';
const INSERT = GB1 + LINK + TEV + CDK;
/* the read as it comes back: a ragged N-rich primer start, the vector, the
   insert, and an N-rich tail where the trace died */
const READ = 'NNNNNNCGACTCN' + VEC.slice(12, HEAD.length + TAG.length) + INSERT + 'NNNNNNNNNNNNNNNNNNNN';

(async () => {
  const browser = await chromium.launch();
  const page = await browser.newPage({ viewport: { width: 1450, height: 1150 } });
  const errors = [];
  page.on('pageerror', e => errors.push('pageerror: ' + String(e)));
  page.on('console', m => { if (m.type() === 'error' && !/Failed to load resource|ERR_(TUNNEL|NAME|INTERNET|CONNECTION)/.test(m.text())) errors.push(m.text()); });
  await page.goto('file://' + path.resolve(__dirname, '..', 'plasmidbench.html'));
  await page.waitForTimeout(800);

  const open = async (read) => {
    await page.evaluate(() => { const m = document.querySelector('#modalRoot'); if (m && !m.hidden) window.PB.closeModal(); });
    await page.evaluate(r => window.PB.toolAnatomy(r), read);
    await page.waitForTimeout(1200);
  };

  console.log('\n== with no construct loaded it describes the fragment only ==');
  await page.evaluate(() => { window.PB.state.records = []; window.PB.state.activeId = null; window.PB.renderAll(); });
  await page.waitForTimeout(300);
  await open(READ);
  const bare = await page.evaluate(() => document.querySelector('#modalRoot').textContent);
  ok('the regions are still found', /GB1 \(protein G B1 domain\)/.test(bare));
  ok('and nothing is claimed about a plasmid that is not there',
     !/already in/.test(bare) && !/are not in/.test(bare),
     (bare.match(/.{0,40}(already in|are not in).{0,40}/) || [''])[0]);

  console.log('\n== with the backbone loaded it says which part is new ==');
  await page.evaluate(([vec, head, tag]) => {
    const rec = window.PB.makeRecord({ name: 'pET-28a(+) backbone', seq: vec, topology: 'circular', features: [
      { start: 0, end: 20, strand: 1, type: 'promoter', label: 'T7 promoter', quals: {} },
      { start: 20, end: 41, strand: 1, type: 'protein_bind', label: 'lac operator', quals: {} },
      { start: 72, end: 78, strand: 1, type: 'RBS', label: 'RBS', quals: {} },
      { start: head + 9, end: head + 27, strand: 1, type: 'CDS', label: '6xHis', quals: {} },
      { start: head + 39, end: head + tag, strand: 1, type: 'misc_feature', label: 'MCS-1', quals: {} }] });
    window.PB.state.records = [rec]; window.PB.state.activeId = rec.id; window.PB.renderAll();
  }, [VEC, HEAD.length, TAG.length]);
  await page.waitForTimeout(400);
  await open(READ);

  const r = await page.evaluate(() => {
    const m = document.querySelector('#modalRoot');
    const rows = Array.from(m.querySelectorAll('#anTable tbody tr:not(.anseq)'))
      .map(tr => Array.from(tr.children).map(td => td.textContent.trim()));
    return { text: m.textContent, rows: rows,
      cols: Array.from(m.querySelectorAll('#anTable thead th')).map(x => x.textContent.trim()),
      hasClone: !!m.querySelector('#anNewClone'), hasCopy: !!m.querySelector('#anNewCopy'),
      strip: m.querySelectorAll('[title*="matches pET"]').length };
  });
  ok('the comparison has its own column in the table',
     r.cols.some(c => /In pET-28a/.test(c)), JSON.stringify(r.cols));
  ok('the vector’s His tag is marked as already there',
     /already there/.test((r.rows.find(x => /His tag/.test(x[2])) || [])[4] || ''),
     JSON.stringify(r.rows.find(x => /His tag/.test(x[2]))));
  ok('and the polylinker in front of the insert too',
     /already there/.test((r.rows.find(x => x[0] === '118–144') || [])[4] || ''),
     JSON.stringify(r.rows.find(x => x[0] === '118–144')));
  ['GB1 (protein G B1 domain)', 'TEV protease site'].forEach(nm => {
    ok(nm + ' is marked as not in it',
       /new — not in it/.test((r.rows.find(x => x[2].indexOf(nm) === 0) || [])[4] || ''),
       JSON.stringify(r.rows.find(x => x[2].indexOf(nm) === 0)));
  });
  const big = r.rows.filter(x => x[2] === 'protein').sort((a, b) =>
    +b[1].replace(/,/g, '') - +a[1].replace(/,/g, ''))[0];
  ok('the big coding region is marked as not in it',
     big && /new — not in it/.test(big[4] || ''), JSON.stringify(big));

  const want = { from: 13 + (HEAD.length + TAG.length - 12), bp: INSERT.length };
  ok('it states which bases are not in the plasmid, and they are the insert',
     r.text.indexOf('Bases ' + (want.from + 1).toLocaleString('en-US') + '–' +
       (want.from + want.bp).toLocaleString('en-US') + ' (' + want.bp.toLocaleString('en-US') +
       ' bp) are not in pET-28a(+) backbone') >= 0,
     (r.text.match(/Bases [\d,]+–[\d,]+ \([\d,]+ bp\) are not in [^.]*/) || [''])[0] +
       '   want from ' + (want.from + 1) + ' bp ' + want.bp);
  ok('it names the plasmid features the rest of the read sat on',
     /T7 promoter/.test(r.text) && /lac operator/.test(r.text) && /6xHis/.test(r.text),
     (r.text.match(/plasmid’s own [^—]{0,90}/) || [''])[0]);
  ok('it says why the front of the read is not something to clone',
     /read starts inside the vector/.test(r.text));
  ok('and offers to clone only the part the plasmid lacks', r.hasClone && r.hasCopy);
  ok('the strip says where each matched stretch sits in the plasmid', r.strip >= 1, r.strip);

  console.log('\n== every region hands over its own sequence ==');
  const gb1Row = await page.evaluate(() => {
    const rows = Array.from(document.querySelectorAll('#anTable tbody tr:not(.anseq)'));
    const i = rows.findIndex(tr => /GB1 \(protein G B1 domain\)/.test(tr.children[2].textContent));
    rows[i].querySelector('[data-anshow]').click();
    return i;
  });
  await page.waitForTimeout(350);
  const shown = await page.evaluate(gi => {
    const row = document.querySelector('[data-anseq="' + gi + '"]');
    const blocks = Array.from(row.querySelectorAll('.seqblock')).map(b => b.textContent.trim());
    return { open: !row.hidden, dna: blocks[0], aa: blocks[1],
      label: row.querySelector('label').textContent,
      btns: Array.from(row.querySelectorAll('button')).map(b => b.textContent.trim()) };
  }, gb1Row);
  ok('the region opens to show its bases', shown.open);
  ok('and the bases are the whole region, not a preview',
     shown.dna === GB1, shown.dna ? shown.dna.length + ' bp, starts ' + shown.dna.slice(0, 16) : shown.dna);
  ok('its translation is given whole as well',
     shown.aa && shown.aa.length === GB1.length / 3 && /^GTYKLILNGKTL/.test(shown.aa) && /TFTVTE$/.test(shown.aa),
     shown.aa);
  ok('the heading says which bases they are',
     /bases 145–312 \(168 bp\)/i.test(shown.label), shown.label);
  ok('and each can be copied on its own',
     shown.btns.indexOf('Copy the bases') >= 0 && shown.btns.indexOf('Copy the protein') >= 0,
     shown.btns.join(' | '));
  const closed = await page.evaluate(gi => {
    const rows = Array.from(document.querySelectorAll('#anTable tbody tr:not(.anseq)'));
    rows[gi].querySelector('[data-anshow]').click();
    return document.querySelector('[data-anseq="' + gi + '"]').hidden;
  }, gb1Row);
  ok('and it closes again', closed);

  console.log('\n== the N-rich tail is not offered as sequence to clone ==');
  ok('the tail is never part of the span',
     !/are not in pET-28a\(\+\) backbone/.test('') &&
     await page.evaluate(() => {
       const t = document.querySelector('#modalRoot').textContent;
       const m = t.match(/Bases [\d,]+–([\d,]+) \(/);
       return m ? +m[1].replace(/,/g, '') : 0;
     }) === want.from + want.bp,
     'span must end at ' + (want.from + want.bp));

  console.log('\n== a read of the plasmid itself has nothing to add ==');
  await open(VEC.slice(0, 700));
  const same = await page.evaluate(() => {
    const m = document.querySelector('#modalRoot');
    return { text: m.textContent, clone: !!m.querySelector('#anNewClone') };
  });
  ok('it says so plainly', /All of this is already in pET-28a\(\+\) backbone/.test(same.text),
     (same.text.match(/.{0,20}already in.{0,80}/) || [''])[0]);
  ok('and offers nothing to clone', !same.clone);

  console.log('\n== console ==');
  ok('no page errors anywhere in that', errors.length === 0, errors.slice(0, 3).join(' | '));

  console.log('\n   ' + pass + ' passed, ' + fail + ' failed\n');
  await browser.close();
  process.exit(fail ? 1 : 0);
})();
