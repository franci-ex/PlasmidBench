/* ===========================================================================
   Sequence verification, through the interface.

   Two things are being checked here, and they are the two halves of what a
   verification actually has to answer. First: which protein did this read
   land on, and what does the read say about it — not "98.7% identity" over
   the whole plasmid, which tells you nothing about whether your insert is
   right, but residue numbers in named features. Second: the BLAST round
   trip. The page has no network access, so the trip is one paste, and the
   paste has to survive all four shapes NCBI's Download menu produces.

   The read is seeded, not random: two mismatches at known offsets inside
   EGFP, placed so the expected consequences can be worked out by hand and
   written into the assertions below rather than read off the screen.
   ========================================================================= */
const { chromium } = require('playwright');
const path = require('path');
let pass = 0, fail = 0;
function ok(name, cond, extra) {
  if (cond) { pass++; console.log('  PASS  ' + name); }
  else { fail++; console.log('  FAIL  ' + name + (extra !== undefined ? '   -> ' + extra : '')); }
}

/* Hit table (text) — what Download -> Hit Table (text) gives you. */
const HITTABLE = `# blastn
# Iteration: 0
# Query: Query_1
# Database: nt
# Fields: query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
# 3 hits found
Query_1\tMN908947.3\t99.833\t600\t1\t0\t1\t600\t21563\t22162\t0.0\t1103
Query_1\tU55761.1\t98.500\t600\t9\t0\t1\t600\t1\t600\t0.0\t1060
Query_1\tAY741333.1\t87.167\t600\t77\t0\t1\t600\t102\t701\t2e-180\t636
`;

/* Descriptions (CSV) — the other common download. */
const DESCCSV = `Description,Scientific Name,Max Score,Total Score,Query Cover,E value,Per. Ident,Acc. Len,Accession
"Aequorea victoria green fluorescent protein (GFP) mRNA, complete cds",Aequorea victoria,1103,1103,100%,0.0,99.83,730,U55761.1
"Synthetic construct enhanced green fluorescent protein (EGFP) gene, complete cds",synthetic construct,1060,1060,100%,0.0,98.50,717,AY741333.1
`;

/* Single-file JSON — the shape BLAST's own API returns. */
const JSON2 = JSON.stringify({
  BlastOutput2: [{
    report: {
      results: {
        search: {
          query_title: 'EGFP_F.ab1', query_len: 600,
          hits: [{
            num: 1,
            description: [{ accession: 'U55761.1', title: 'green fluorescent protein mRNA', sciname: 'Aequorea victoria' }],
            len: 730,
            hsps: [{ num: 1, bit_score: 1103.2, evalue: 0, identity: 599, align_len: 600,
                     query_from: 1, query_to: 600, hit_from: 1, hit_to: 600, gaps: 0 }]
          }]
        }
      }
    }
  }]
});

/* Pairwise text — what people copy straight off the results page. */
const PAIRWISE = `
Sequences producing significant alignments:
                                                                   Score     E
Sequences producing significant alignments:                       (Bits)  Value

> U55761.1 Aequorea victoria green fluorescent protein (GFP) mRNA, complete cds
Length=730

 Score = 1103 bits (597),  Expect = 0.0
 Identities = 599/600 (99%), Gaps = 0/600 (0%)
 Strand=Plus/Plus

Query  1     ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGAC  60
             ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Sbjct  1     ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGAC  60
`;

(async () => {
  const browser = await chromium.launch();
  const page = await browser.newPage({ viewport: { width: 1500, height: 1100 } });
  const errors = [];
  page.on('pageerror', e => errors.push('pageerror: ' + String(e)));
  page.on('console', m => { if (m.type() === 'error' && !/Failed to load resource/.test(m.text())) errors.push('console: ' + m.text()); });
  await page.goto('file://' + path.resolve('plasmidbench.html'));
  await page.waitForTimeout(800);

  /* ---------- a read on a known construct ---------- */
  console.log('\n== which protein did this read land on ==');
  const seeded = await page.evaluate(() => {
    const rec = window.PB.active();
    if (!rec) return { error: 'no demo record' };
    const egfp = rec.features.find(f => /EGFP|GFP/i.test(f.label || '') && window.Bio.isCodingFeature(f));
    if (!egfp) return { error: 'no EGFP in the demo' };
    const idx = window.Bio.cdsIndex(rec, egfp);
    /* Read 600 bp of the reference starting 180 bp before the CDS, then
       break two codons inside it. Residue 54 and residue 105 of EGFP, in
       the read's own coordinates, are worked out from the CDS index so the
       test does not depend on where the demo happens to put the gene. */
    const g54 = idx.pos[(54 - 1) * 3 + 1];      // middle base of codon 54
    const g105 = idx.pos[(105 - 1) * 3 + 1];    // middle base of codon 105
    const start = Math.max(0, Math.min(g54, g105) - 120);
    const seq = rec.seq.slice(start, start + 600).split('');
    const put = (g, b) => { const o = g - start; if (o >= 0 && o < seq.length) seq[o] = b; };
    const flip = c => ({ A: 'C', C: 'A', G: 'T', T: 'G' }[c] || 'A');
    put(g54, flip(rec.seq[g54]));
    put(g105, flip(rec.seq[g105]));
    window.PB.addRead({ name: 'EGFP_F.ab1', seq: seq.join(''), source: 'test' });
    return { gene: egfp.label, start: start, len: 600 };
  });
  ok('the demo has a coding EGFP to read against', !seeded.error, seeded.error);
  await page.waitForTimeout(900);

  const detail = await page.evaluate(() => {
    window.PB.state.verifyDetail = 0;
    window.PB.renderViews();
    return null;
  });
  await page.waitForTimeout(500);

  const cov = await page.evaluate(() => {
    const host = document.querySelector('#view-verify') || document.body;
    const svg = host.querySelector('svg.covmap');
    const tbls = Array.from(host.querySelectorAll('table.tbl'));
    const covTbl = tbls.find(t => /Residues reached/.test(t.textContent));
    const rows = covTbl ? Array.from(covTbl.querySelectorAll('tbody tr')).map(tr =>
      Array.from(tr.children).map(td => td.textContent.trim())) : [];
    return {
      hasSvg: !!svg,
      ticks: svg ? svg.querySelectorAll('.covdiff').length : -1,
      readSpans: svg ? svg.querySelectorAll('.covread').length : -1,
      rows: rows,
      text: host.textContent
    };
  });
  ok('it draws where the read landed', cov.hasSvg);
  ok('and marks the read span on it', cov.readSpans >= 1, cov.readSpans);
  ok('and ticks both differences', cov.ticks === 2, cov.ticks);
  ok('and names the feature the read covers', cov.rows.some(r => /GFP/i.test(r[0] || '')),
    JSON.stringify(cov.rows.map(r => r[0])));

  const gfpRow = cov.rows.find(r => /GFP/i.test(r[0] || '')) || [];
  ok('it reports the residues the read actually reached', /\d+–\d+ of \d+/.test(gfpRow[4] || ''), gfpRow[4]);
  ok('and both changes by residue number', /54/.test(gfpRow[5] || '') && /105/.test(gfpRow[5] || ''), gfpRow[5]);
  ok('a residue the read never reached is not invented',
    !/all \d+ residues/.test(gfpRow[4] || ''), gfpRow[4]);

  /* The consequence calls themselves, against the engine, so the table is
     not grading its own arithmetic. */
  const eff = await page.evaluate(() => {
    const rec = window.PB.active(), rd = window.PB.state.reads[0];
    const egfp = rec.features.find(f => /EGFP|GFP/i.test(f.label || '') && window.Bio.isCodingFeature(f));
    return rd.aln.differences.map(d => {
      const e = window.Bio.diffEffect(rec, egfp, d, {});
      return e ? { kind: e.kind, res: e.residue } : null;
    });
  });
  ok('every difference inside the CDS gets a consequence', eff.length === 2 && eff.every(e => e), JSON.stringify(eff));
  ok('and the residue numbers are 54 and 105',
    eff.length === 2 && eff.map(e => e && e.res).sort((a, b) => a - b).join(',') === '54,105',
    JSON.stringify(eff));

  /* ---------- the bases behind each peptide ---------- */
  console.log('\n== copying the ORF as nucleotides ==');
  /* A peptide you can read is only half of what a verification is for: the
     other half is the DNA behind it, to put in a primer or an order. The
     thing that has to hold is that the bases handed over translate to the
     peptide shown beside them — on a minus frame that means the reverse
     complement of the read, which is the strand the gene is on and not
     what is in the trace file. */
  const ntCheck = await page.evaluate(() => {
    const out = { rows: [], api: true };
    const fr = window.PB.state.verifyFrames;
    if (!fr || !fr.length) { out.api = false; return out; }
    fr.forEach(r => {
      const pep = r.clean || r.longest, nt = r.clean ? r.cleanNt : r.longestNt;
      out.rows.push({
        label: r.label, dir: r.dir, pep: pep.length, nt: nt.length,
        backTranslates: window.Bio.translate(nt, {}) === pep,
        thirdOfIt: nt.length === pep.length * 3,
        metOk: !r.fromMet || window.Bio.translate(r.fromMetNt, {}) === r.fromMet,
        metStartsATG: !r.fromMet || /^ATG/.test(r.fromMetNt)
      });
    });
    return out;
  });
  ok('every frame carries the bases of its own peptide', ntCheck.api && ntCheck.rows.length === 6,
     ntCheck.rows.length);
  ok('and each one is exactly three bases per residue',
     ntCheck.rows.every(r => r.thirdOfIt), JSON.stringify(ntCheck.rows.filter(r => !r.thirdOfIt)));
  ok('and translates back to the peptide shown beside it',
     ntCheck.rows.every(r => r.backTranslates), JSON.stringify(ntCheck.rows.filter(r => !r.backTranslates)));
  ok('including on the minus frames, where they are the reverse complement',
     ntCheck.rows.filter(r => r.dir === -1).length === 3 &&
     ntCheck.rows.filter(r => r.dir === -1).every(r => r.backTranslates),
     JSON.stringify(ntCheck.rows.filter(r => r.dir === -1).map(r => r.label + ':' + r.backTranslates)));
  ok('the from-Met bases translate to the from-Met peptide',
     ntCheck.rows.every(r => r.metOk), JSON.stringify(ntCheck.rows.filter(r => !r.metOk)));
  ok('and begin at an ATG, which is what makes them worth cloning',
     ntCheck.rows.every(r => r.metStartsATG), JSON.stringify(ntCheck.rows.filter(r => !r.metStartsATG)));

  /* and the buttons are really there, in both tables */
  const btns = await page.evaluate(() => {
    const host = document.querySelector('#view-verify') || document.body;
    const t = Array.from(host.querySelectorAll('#trTable tbody tr'))
      .map(tr => Array.from(tr.querySelectorAll('button')).map(b => b.textContent.trim()))
      .filter(a => a.length);
    return { main: t[0] || [],
      detailNt: host.querySelectorAll('[data-fnt]').length,
      detailMet: host.querySelectorAll('[data-fntmet]').length,
      mainNt: host.querySelectorAll('[data-trnt]').length };
  });
  ok('the read table offers both copies side by side',
     btns.main.indexOf('Copy aa') >= 0 && btns.main.indexOf('Copy nt') >= 0, JSON.stringify(btns.main));
  ok('one Copy nt per read', btns.mainNt >= 1, btns.mainNt);
  ok('and the six-frame table offers one per frame', btns.detailNt === 6, btns.detailNt);
  ok('with nt from Met wherever a Met is found', btns.detailMet >= 1, btns.detailMet);

  /* what actually lands on the clipboard */
  await page.evaluate(() => navigator.clipboard && 0);
  const copied = await page.evaluate(() => {
    let grabbed = null;
    const real = window.PB.copyText;
    window.PB.copyText = (t, what) => { grabbed = { text: t, what: what }; };
    document.querySelector('[data-trnt]').click();
    window.PB.copyText = real;
    const fr = window.PB.state.trCache || null;
    return grabbed;
  });
  ok('clicking Copy nt copies bases, not residues',
     copied && /^[ACGTN]+$/.test(copied.text), copied && copied.text.slice(0, 40));
  ok('and says how many and on which frame',
     copied && /frame [+\u2212]\d/.test(copied.what) && /\d+ bp/.test(copied.what), copied && copied.what);

  /* ---------- the BLAST round trip ---------- */
  console.log('\n== BLAST results pasted back ==');
  const shapes = [
    ['Hit Table', HITTABLE, 3, 'MN908947.3'],
    ['Descriptions CSV', DESCCSV, 2, 'U55761.1'],
    ['Single-file JSON', JSON2, 1, 'U55761.1'],
    ['pairwise text', PAIRWISE, 1, 'U55761.1']
  ];
  for (const [label, text, n, firstAcc] of shapes) {
    await page.evaluate(() => { window.PB.state.blast = null; window.PB.renderViews(); });
    await page.waitForTimeout(200);
    const opened = await page.evaluate(() => {
      const b = document.querySelector('#blPaste');
      if (!b) return false;
      b.click();
      return true;
    });
    ok(label + ': the paste dialog opens', opened);
    await page.waitForTimeout(250);
    await page.evaluate(t => { document.querySelector('#blText').value = t; }, text);
    await page.evaluate(() => {
      const b = Array.from(document.querySelectorAll('#modalRoot button')).find(x => x.textContent.trim() === 'Read it');
      b.click();
    });
    await page.waitForTimeout(400);
    const got = await page.evaluate(() => {
      const b = window.PB.state.blast;
      const host = document.querySelector('#view-verify') || document.body;
      const tbl = Array.from(host.querySelectorAll('table.tbl')).find(t => /Accession/.test(t.textContent) && /Bits/.test(t.textContent));
      const rows = tbl ? Array.from(tbl.querySelectorAll('tbody tr')).map(tr =>
        Array.from(tr.children).map(td => td.textContent.trim())) : [];
      return { stored: b ? b.hits.length : 0, format: b ? b.format : '', rows: rows, open: !!document.querySelector('#blText') };
    });
    ok(label + ': the dialog closes on a good paste', !got.open);
    ok(label + ': ' + n + ' hit(s) read', got.stored === n, got.stored);
    ok(label + ': all of them reach the table', got.rows.length === n, got.rows.length);
    ok(label + ': the first accession is ' + firstAcc, (got.rows[0] || [])[0] === firstAcc, (got.rows[0] || [])[0]);
    ok(label + ': the panel says which format it read', new RegExp(got.format, 'i').test(got.format) && got.format.length > 2, got.format);
  }

  /* a column the format does not carry must stay blank, not be guessed */
  await page.evaluate(() => { window.PB.state.blast = null; window.PB.renderViews(); });
  await page.evaluate(() => document.querySelector('#blPaste').click());
  await page.waitForTimeout(200);
  await page.evaluate(t => { document.querySelector('#blText').value = t; }, HITTABLE);
  await page.evaluate(() => {
    Array.from(document.querySelectorAll('#modalRoot button')).find(x => x.textContent.trim() === 'Read it').click();
  });
  await page.waitForTimeout(350);
  const blankCover = await page.evaluate(() => {
    const host = document.querySelector('#view-verify') || document.body;
    const tbl = Array.from(host.querySelectorAll('table.tbl')).find(t => /Accession/.test(t.textContent) && /Bits/.test(t.textContent));
    const r = tbl ? Array.from(tbl.querySelectorAll('tbody tr')[0].children).map(td => td.textContent.trim()) : [];
    return { identity: r[2], cover: r[3], evalue: r[4], bits: r[5] };
  });
  ok('a hit table keeps its identity', blankCover.identity === '99.8%', blankCover.identity);
  ok('and leaves query cover blank rather than guessing it', blankCover.cover === '—', blankCover.cover);
  ok('and keeps the E value and bit score', blankCover.evalue === '0.0' && blankCover.bits === '1103',
    blankCover.evalue + ' / ' + blankCover.bits);

  /* refusals */
  console.log('\n== what it refuses to read ==');
  for (const [label, text] of [['empty paste', '   '], ['unrelated text', 'Dear sir, please find attached the results of our experiment.']]) {
    await page.evaluate(() => { window.PB.state.blast = null; window.PB.renderViews(); });
    await page.evaluate(() => document.querySelector('#blPaste').click());
    await page.waitForTimeout(180);
    await page.evaluate(t => { document.querySelector('#blText').value = t; }, text);
    await page.evaluate(() => {
      Array.from(document.querySelectorAll('#modalRoot button')).find(x => x.textContent.trim() === 'Read it').click();
    });
    await page.waitForTimeout(300);
    const st = await page.evaluate(() => ({ stored: !!window.PB.state.blast, open: !!document.querySelector('#blText') }));
    ok(label + ': nothing is stored', !st.stored);
    ok(label + ': and the dialog stays open to fix it', st.open);
  }
  await page.evaluate(() => { const b = Array.from(document.querySelectorAll('#modalRoot button')).find(x => /Cancel/.test(x.textContent)); if (b) b.click(); });
  await page.waitForTimeout(200);

  /* clearing */
  await page.evaluate(t => {
    window.PB.state.blast = null; window.PB.renderViews();
    document.querySelector('#blPaste').click();
  });
  await page.waitForTimeout(200);
  await page.evaluate(t => { document.querySelector('#blText').value = t; }, HITTABLE);
  await page.evaluate(() => Array.from(document.querySelectorAll('#modalRoot button')).find(x => x.textContent.trim() === 'Read it').click());
  await page.waitForTimeout(350);
  await page.evaluate(() => { const b = document.querySelector('#blClear'); if (b) b.click(); });
  await page.waitForTimeout(300);
  ok('Clear puts the panel back', await page.evaluate(() => !window.PB.state.blast && !!document.querySelector('#blPaste')));

  console.log('\n== console ==');
  ok('no page errors anywhere in that', errors.length === 0, errors.slice(0, 4).join(' | '));

  console.log('\n   ' + pass + ' passed, ' + fail + ' failed\n');
  await browser.close();
  process.exit(fail ? 1 : 0);
})();
