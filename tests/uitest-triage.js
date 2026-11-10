// Read triage in the verify view: the verdict a read gets, the warnings it
// carries, and what leaves the tool in the exported FASTA.
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

const BT = { A: 'GCT', R: 'CGT', N: 'AAT', D: 'GAT', C: 'TGT', Q: 'CAA', E: 'GAA',
             G: 'GGT', H: 'CAT', I: 'ATT', L: 'CTG', K: 'AAA', M: 'ATG', F: 'TTT',
             P: 'CCT', S: 'TCT', T: 'ACT', W: 'TGG', Y: 'TAT', V: 'GTT' };
const back = p => p.split('').map(c => BT[c] || 'NNN').join('');
const CDK1 = 'MEDYTKIEKIGEGTYGVVYKGRHKTTGQVVAMKKIRLESEEEGVPSTAIREISLLKELRHPNIVSLQDVLMQDSRLYLIFEFLSMDLKKYLDSIPPGQYMDSSLVKSYLYQILQGIVFCHSRRVLHRDLKPQNLLIDDKGTIKLADFGLARAFGIPIRVYTHEVVTLWYRSPEVLLGSARYSTPVDIWSIGTIFAELATKKPLFHGDSEIDQLFRIFRALGTPNNEVWPEVESLQDYKNTFPKWKPGSLASHVKNLDENGLDLLSKMLIYDPAKRISGKMALNHPYFNDLDNQIKKM';
const LEADER = 'ATGGGCAGCAGCCATCATCATCATCATCACAGCAGCGGCGTGGATCTGGGTACCGAGAACCTGTACTTCCAGAGC';
const TERM = 'CTAGCATAACCCCTTGGGGCCTCTAAACGGGTCTTGAGGGGTTTTTTGCTGAAAGGAGG';
const VECUP = 'GATCCGGCGTAGAGGATCGAGATCTCGATCCCGCGAAATTAATACGACTCACTATAGGG';
const GOOD = VECUP + LEADER + back(CDK1) + 'TAA' + TERM;
const NONCODING = randDNA(700, 991);

/* the reads we hand it: one clean tagged read, a dead trace, three copies of
   one template, a non-coding backbone read, and a read with a single-base
   deletion partway through */
const DELAT = VECUP.length + LEADER.length + 300;
const READS = [
  { name: 'clone1_T7', seq: GOOD },
  { name: 'clone2_T7', seq: (function () { let s = ''; for (let i = 0; i < GOOD.length; i++) s += (i % 7 === 0 ? 'N' : GOOD[i]); return s; })() },
  { name: 'clone3_T7', seq: NONCODING },
  { name: 'clone4_T7', seq: NONCODING.slice(12, 690) },
  { name: 'clone5_T7', seq: B.revcomp(NONCODING.slice(4)) },
  { name: 'clone6_T7', seq: GOOD.slice(0, DELAT) + GOOD.slice(DELAT + 1) }
];

(async () => {
  const browser = await chromium.launch();
  const page = await browser.newPage({ viewport: { width: 1500, height: 1250 } });
  const errors = [];
  page.on('pageerror', e => errors.push(String(e)));
  /* The page asks Google Fonts for IBM Plex and falls back to a system font
     when it cannot have it; in this sandbox it never can, and that failed
     request is not a fault in the page. */
  page.on('console', m => { if (m.type() === 'error' && !/Failed to load resource/.test(m.text())) errors.push(m.text()); });
  await page.goto('file://' + path.resolve('plasmidbench.html'));
  await page.waitForTimeout(700);

  console.log('\n== triage of a mixed set of reads ==');
  await page.evaluate(([seq, reads]) => {
    window.PB.addRecord({ name: 'CDK1 construct', seq: seq, circular: true, features: [] });
    reads.forEach(r => window.PB.addRead({ name: r.name, seq: r.seq, source: 'test' }));
  }, [GOOD, READS]);
  await page.waitForTimeout(400);
  await page.click('[data-v="verify"]');
  await page.waitForTimeout(1200);

  const rows = await page.evaluate(() => {
    const cache = window.PB.state.triageCache;
    return cache.rows.map(r => ({
      read: r.read, v: r.t.verdict, frame: r.t.chosen.label, clean: r.t.chosen.cleanLen,
      stops: r.t.chosen.stops, amb: +(r.t.ambiguousFrac * 100).toFixed(1),
      pinned: r.t.chosen.pinned, dup: r.dupGroup || 0, flags: r.t.flags.map(f => f.level),
      shift: r.t.shift ? r.t.shift.at : null
    }));
  });
  console.log(rows.map(r => '    ' + r.read.padEnd(11) + r.v.padEnd(17) + r.frame +
    ' clean=' + String(r.clean).padStart(4) + ' stops=' + String(r.stops).padStart(3) +
    ' N=' + String(r.amb).padStart(5) + '%' + (r.pinned ? ' pinned' : '') +
    (r.dup ? ' dup' + r.dup : '') + (r.shift !== null ? ' shift@' + r.shift : '')).join('\n'));

  const byName = {};
  rows.forEach(r => { byName[r.read] = r; });
  ok('six reads are triaged', rows.length === 6, rows.length);
  ok('the clean tagged read is called clean', byName.clone1_T7.v === 'clean', byName.clone1_T7.v);
  ok('and its frame is pinned by the tag', byName.clone1_T7.pinned);
  ok('the dead trace is condemned', byName.clone2_T7.v === 'low-quality', byName.clone2_T7.v);
  ok('and it is condemned despite having few stops',
     byName.clone2_T7.stops < 25, 'stops=' + byName.clone2_T7.stops);
  ok('the dead trace carries a hard warning',
     byName.clone2_T7.flags.indexOf('bad') >= 0, byName.clone2_T7.flags.join(','));
  ok('the non-coding read is not presented as a protein',
     ['open-only', 'frame-ambiguous', 'no-protein'].indexOf(byName.clone3_T7.v) >= 0, byName.clone3_T7.v);
  ok('and it carries a warning rather than nothing',
     byName.clone3_T7.flags.some(f => f === 'warn' || f === 'bad'), byName.clone3_T7.flags.join(','));
  ok('the one-base deletion is caught', byName.clone6_T7.shift !== null, byName.clone6_T7.shift);
  ok('and the deletion is located near the right base',
     Math.abs(byName.clone6_T7.shift - DELAT) <= 60, byName.clone6_T7.shift + ' vs ' + DELAT);

  console.log('\n== the same DNA under three names ==');
  const dupNames = await page.evaluate(() => (window.PB.state.triageCache.dups || []).map(g => g.map(x => (x && x.read) || x)));
  console.log('    DUPGROUPS ' + JSON.stringify(dupNames));
  const dups = await page.evaluate(() => (window.PB.state.triageCache.dups || []).map(g => g.length));
  /* Two groups, not one. clone3/4/5 are three views of the same non-coding
     stretch, which is the group this section was written for. clone1 and
     clone6 are the same 1,275 base construct differing by a single deleted
     base, and a detector that did not group those would miss two reads off
     one template whenever a polymerase slipped — so the group is right. What
     matters is that the panel does not then invite the reader to ignore
     clone6, which is the read that found the frameshift. */
  ok('both duplicate groups are found', dups.length === 2 &&
     dups.slice().sort().join(',') === '2,3', JSON.stringify(dups));
  ok('the three views of one stretch are grouped together',
     dupNames.some(g => g.length === 3), JSON.stringify(dupNames));
  ok('and the construct read twice is grouped as well',
     dupNames.some(g => g.length === 2), JSON.stringify(dupNames));
  const split = await page.evaluate(() => {
    const el = document.querySelector('[data-dupsplit]');
    return el ? el.textContent : '';
  });
  ok('a group whose reads disagree is called out separately',
     /did not get the same verdict/.test(split), split.slice(0, 80));
  ok('and it names both verdicts rather than picking one',
     /clean/.test(split) && /frameshift/.test(split), split.slice(0, 160));
  ok('and it says not to treat either read as confirming the other',
     /confirming the other/.test(split), split.slice(0, 60));
  ok('all three copies are tagged in the table',
     ['clone3_T7', 'clone4_T7', 'clone5_T7'].every(n => byName[n].dup === byName.clone3_T7.dup &&
       byName[n].dup > 0), JSON.stringify(['clone3_T7', 'clone4_T7', 'clone5_T7'].map(n => byName[n].dup)));
  ok('the reverse-complement copy groups with the others', byName.clone5_T7.dup > 0);
  const banner = await page.evaluate(() => {
    const el = document.querySelector('#view-verify');
    return el ? el.textContent : '';
  });
  /* The banner was reworded to say what to do about it rather than stating a
     rule, so match the meaning: reads off one template are not independent
     evidence about different clones. */
  ok('the view warns that independent clones cannot give the same sequence',
     /not independent evidence/.test(banner) && /different clones/.test(banner),
     banner.slice(banner.indexOf('Some reads are the same DNA'), banner.indexOf('Some reads are the same DNA') + 200));
  ok('the duplicate banner names the reads', /clone3_T7, clone4_T7, clone5_T7/.test(banner));

  console.log('\n== what the user sees and exports ==');
  const shown = await page.evaluate(() => {
    const t = document.querySelector('#view-verify').textContent;
    return { clean: /Clean coding frame/.test(t), dead: /Trace failed/.test(t),
             nTrans: /N translates to X/.test(t), worth: /worth sending to BLAST/.test(t),
             counts: (t.match(/\d+ (clean coding frame|trace failed|open stretch[^<]*)/g) || []).length };
  });
  ok('the clean verdict is on screen', shown.clean);
  ok('the failed trace is on screen', shown.dead);
  ok('the N-versus-stop trap is explained on screen', shown.nTrans);
  ok('the panel says how many reads are worth BLASTing', shown.worth);

  const fasta = await page.evaluate(() => {
    const items = window.PB.state.trItems;
    return items.map(x => {
      const pep = x.f.clean || x.f.longest;
      const pins = (x.f.landmarks || []).filter(l => l.pins).map(l => l.name.replace(/\s+/g, '_'));
      return '>' + x.read.replace(/\s+/g, '_') + '_frame' + x.f.label.replace('−', 'minus') +
        ' ' + pep.length + 'aa stops=' + x.f.stops +
        ' ambiguous=' + (x.t.ambiguousFrac * 100).toFixed(1) + '%' +
        (x.f.fit === null ? '' : ' composition=' + x.f.fit.toFixed(2)) +
        ' verdict=' + x.t.verdict + (pins.length ? ' pinned_by=' + pins.join(',') : '');
    });
  });
  console.log(fasta.map(h => '    ' + h).join('\n'));
  ok('every FASTA header carries a verdict', fasta.every(h => /verdict=/.test(h)));
  ok('every FASTA header carries the ambiguous-base fraction', fasta.every(h => /ambiguous=/.test(h)));
  ok('the pinned read says what pinned it', fasta.some(h => /pinned_by=/.test(h)));
  ok('no header contains a unicode minus that would break a parser',
     fasta.every(h => h.indexOf('−') < 0));

  console.log('\n== the filter ==');
  await page.click('#trOnlyOk'); await page.waitForTimeout(500);
  /* Scoped to the triage table itself. The verify view has several tables in
     it now — frames, proteins, coverage, BLAST hits — and a selector that
     swept them all up was reading frame labels as if they were read names. */
  const filtered = await page.evaluate(() => {
    const rows = Array.from(document.querySelectorAll('#trTable tbody tr'))
      .filter(tr => tr.children.length > 3);
    return rows.map(tr => tr.children[0].textContent.trim()).filter(Boolean);
  });
  ok('the filter hides the reads not worth BLASTing',
     filtered.indexOf('clone2_T7') < 0 && filtered.some(n => /clone1_T7/.test(n)),
     filtered.join(','));
  await page.click('#trOnlyOk'); await page.waitForTimeout(400);

  console.log('\n== all six frames ==');
  await page.click('[data-trm="all"]'); await page.waitForTimeout(700);
  const six = await page.evaluate(() => window.PB.state.trItems.length);
  ok('six frames per read are listed', six === 36, six);
  await page.click('[data-trm="best"]'); await page.waitForTimeout(500);

  await page.screenshot({ path: 'shot-triage.png', fullPage: false });
  ok('no errors anywhere in the run', errors.length === 0, errors.slice(0, 3).join(' | '));

  console.log('\n   ' + pass + ' passed, ' + fail + ' failed\n');
  await browser.close();
  process.exit(fail ? 1 : 0);
})();
