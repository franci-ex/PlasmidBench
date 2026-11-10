// Seamless In-Fusion: pick an insert, a vector and a place, and everything
// follows with no further clicks. Driven on the real pRSFDuet-1 file.
const { chromium } = require('playwright');
const path = require('path'), fs = require('fs');
const B = require('../src/engine.js'); require('../src/engine2.js');
let pass = 0, fail = 0;
function ok(name, cond, extra) {
  if (cond) { pass++; console.log('  PASS  ' + name); }
  else { fail++; console.log('  FAIL  ' + name + (extra !== undefined ? '   -> ' + extra : '')); }
}
const VEC = path.resolve(__dirname, 'fixtures/pRSFDuet-1.dna');
/* the user's own CDK1 coding sequence */
const CDK1 = 'ATGGAAGATTATACCAAAATAGAGAAAATAGGAGAAGGAACTTACGGAGTTGTGTATAAAGGACGTCATAAAACAACAGGACAGGTAGTGGCAATGAAAAAGATTCGTCTAGAAAGTGAAGAGGAAGGGGTTCCTTCAACTGCAATTCGGGAAATCTCTCTGTTAAAGGAACTTCGACATCCAAATATAGTCAGTCTTCAGGATGTGCTTATGCAGGATTCCAGGTTATATTTGATTTTTGAATTTCTTTCCATGGATCTGAAAAAATACTTGGATTCTATCCCTCCTGGCCAGTACATGGATTCTAGTTTAGTAAAAAGTTATTTGTACCAAATACTACAAGGAATTGTGTTTTGCCATTCAAGAAGAGTATTACACAGAGATCTAAAACCACAGAATCTGCTCATTGATGATAAAGGAACCATTAAACTAGCAGACTTTGGACTAGCCAGAGCTTTTGGAATACCTATTCGTGTATATACACATGAGGTAGTAACACTCTGGTATAGATCTCCAGAGGTACTTCTTGGATCTGCAAGATACTCCACACCTGTTGACATCTGGAGCATTGGTACTATATTTGCAGAATTAGCAACTAAAAAGCCTCTGTTTCATGGAGATTCTGAGATTGACCAACTCTTCCGGATTTTTCGGGCCTTAGGAACTCCAAATAATGAGGTGTGGCCTGAAGTGGAATCTTTGCAAGACTATAAAAACACATTCCCCAAATGGAAACCAGGCAGCCTTGCTTCTCATGTTAAAAACCTTGATGAAAATGGATTGGATCTTCTCTCAAAAATGCTTATTTATGATCCTGCTAAACGAATTTCTGGAAAAATGGCACTTAATCATCCTTATTTCAATGATTTGGATAACCAAATAAAGAAGATGTAG';

(async () => {
  const browser = await chromium.launch();
  const page = await browser.newPage({ viewport: { width: 1500, height: 1100 } });
  const errors = [];
  page.on('pageerror', e => errors.push(String(e)));
  page.on('console', m => { if (m.type() === 'error' && !/Failed to load resource|ERR_(TUNNEL|NAME|INTERNET|CONNECTION)/.test(m.text())) errors.push(m.text()); });
  await page.goto('file://' + path.resolve('plasmidbench.html'));
  await page.waitForTimeout(600);
  if (!fs.existsSync(VEC)) { console.log('SKIP: vector file not present'); await browser.close(); process.exit(0); }
  const b64 = fs.readFileSync(VEC).toString('base64');
  await page.evaluate(([b64, cdk]) => {
    const bin = atob(b64), buf = new Uint8Array(bin.length);
    for (let i = 0; i < bin.length; i++) buf[i] = bin.charCodeAt(i);
    const r = window.Bio.parseSnapGene(buf.buffer);
    const P = window.PB;
    P.addRecord(Object.assign({ name: 'pRSFDuet-1' }, r, { topology: 'circular' }));
    /* the insert as it would arrive: a record whose CDS is the gene */
    P.addRecord(P.makeRecord({ name: 'CDK1 source', seq: 'GGCACGAGG' + cdk + 'TTAACCGGTT', topology: 'linear',
      features: [{ start: 9, end: 9 + cdk.length, strand: 1, type: 'CDS', label: 'CDK1', quals: {} }] }));
    /* The dialog takes the construct on screen as the backbone, which is
       right — it is the one you are looking at — but adding the insert second
       leaves that one active. Go back to the vector, which is what a cloner
       would have open before reaching for the tool. */
    const vrec = P.state.records.find(x => x.seq === r.seq);
    P.state.activeId = vrec.id;
    P.renderAll();
  }, [b64, CDK1]);
  await page.waitForTimeout(400);

  const live = () => page.evaluate(() => document.querySelector('#gbLive').textContent);
  const openDlg = async () => {
    await page.evaluate(() => { const m = document.querySelector('#modalRoot'); if (!m.hidden) window.PB.closeModal(); });
    await page.evaluate(() => window.PB.openTools());
    await page.waitForTimeout(250);
    await page.evaluate(() => {
      const b = Array.from(document.querySelectorAll('#modalRoot button, #modalRoot [role=button], #modalRoot .tool'))
        .find(x => /In-Fusion/.test(x.textContent));
      b.click();
    });
    await page.waitForTimeout(700);
  };

  console.log('\n== it opens ready, with sensible choices already made ==');
  await openDlg();
  /* The dialog opens on "vector only", which is a real starting state — a
     backbone can be re-closed on itself with no insert at all. Choosing the
     record that holds the gene is the one decision it cannot make for you,
     since any record in the library might be the insert. Everything after
     that choice is what must need no further clicks. */
  const insId = await page.evaluate(() => window.PB.state.records.find(r => r.name === 'CDK1 source').id);
  const opened = await page.evaluate(() => ({
    insSrc: document.querySelector('#gbIns_src').selectedOptions[0].textContent,
    parts: Array.from(document.querySelector('#gbInsPart').options).map(o => o.textContent)
  }));
  ok('it opens with no insert chosen rather than guessing one',
     /vector only/i.test(opened.insSrc) || /vector only/i.test(opened.parts.join(' ')),
     opened.insSrc + ' | ' + JSON.stringify(opened.parts));
  await page.selectOption('#gbIns_src', insId);
  await page.waitForTimeout(600);
  const first = await page.evaluate(() => ({
    insSrc: document.querySelector('#gbIns_src').selectedOptions[0].textContent,
    parts: Array.from(document.querySelector('#gbInsPart').options).map(o => o.textContent),
    partPicked: document.querySelector('#gbInsPart').selectedOptions[0].textContent,
    vec: document.querySelector('#gbVec_src').selectedOptions[0].textContent,
    where: Array.from(document.querySelector('#gbWhere').options).map(o => o.textContent),
    name: document.querySelector('#gbInsName').value,
    live: document.querySelector('#gbLive').textContent
  }));
  ok('the insert record offers its coding features by name',
     first.parts.some(p => /^CDK1 · CDS/.test(p)), JSON.stringify(first.parts));
  ok('and a coding feature is picked for you', /^CDK1 · CDS/.test(first.partPicked), first.partPicked);
  ok('the insert is named from it', first.name === 'CDK1', first.name);
  /* The placement list is phrased as what is left out rather than what is cut
     — "everything except KanR" rather than "cut with NcoI" — because the
     backbone is opened by PCR here, not by an enzyme. */
  /* With an insert chosen the list is about where the gene goes, and it is
     specific to this vector: both of pRSFDuet-1's cloning sites, each offered
     behind its own tag in frame, untagged from the ATG, or replacing the
     whole site. */
  ok('the placement list is built from the vector itself',
     first.where.some(x => /MCS-1/.test(x)) && first.where.some(x => /MCS-2/.test(x)),
     JSON.stringify(first.where.slice(0, 6)));
  ok('and it offers each tag as an in-frame fusion point',
     first.where.some(x => /after 6xHis, in frame/.test(x)) &&
     first.where.some(x => /after S-Tag, in frame/.test(x)),
     JSON.stringify(first.where.filter(x => /in frame/.test(x))));
  ok('and an untagged placement that keeps the vector\u2019s spacing',
     first.where.some(x => /untagged/.test(x) && /RBS/.test(x)),
     JSON.stringify(first.where.filter(x => /untagged/.test(x))));

  console.log('\n== with nothing else clicked, the whole design is already there ==');
  const ins = insId;
  await page.waitForTimeout(500);
  await page.evaluate(() => {
    /* Pick by the option's own value, not by where it sits in the list. The
       values are indices into the placement list and the options are grouped
       into optgroups in a reading order, so the two do not line up — and
       setting the value to a DOM position silently chooses a different
       placement, which is how this test used to end up designing a BamHI
       digest while believing it had asked for an in-frame fusion. */
    const sel = document.querySelector('#gbWhere');
    const o = Array.from(sel.options).find(x => /in frame/.test(x.textContent));
    if (!o) throw new Error('no in-frame placement offered');
    sel.value = o.value; sel.dispatchEvent(new Event('change', { bubbles: true }));
  });
  await page.waitForTimeout(700);
  const d = await page.evaluate(() => {
    const t = document.querySelector('#gbLive').textContent;
    /* Scoped to the primer table. The live panel also holds the cycling
       conditions and the junction readout, and a selector that swept up every
       tbody row was reading thermocycler steps as if they were primers. */
    const rows = Array.from(document.querySelectorAll('#gbLive table.primers tbody tr')).map(tr =>
      Array.from(tr.children).map(td => td.textContent.trim()));
    return { t: t, rows: rows, save: !!document.querySelector('#gbSave') };
  });
  ok('the product is reported without pressing anything', /bp circular/.test(d.t), d.t.slice(0, 120));
  ok('both primer pairs are listed', d.rows.length === 4, JSON.stringify(d.rows.map(r => r[0])));
  /* The table grew a "Tm full" column between the annealing Tm and the tail:
     a primer with a 15 bp homology arm has two Tms that matter, the one it
     anneals at in the PCR and the one of the whole oligo. Tail is the last
     column. */
  ok('the insert primers carry a homology tail', d.rows[0] && +d.rows[0][6] === 15, d.rows[0] && d.rows[0][6]);
  ok('and a Tm for the whole oligo as well as for the annealing part',
     d.rows[0] && /^\d+\.\d$/.test(d.rows[0][5]), d.rows[0] && d.rows[0][5]);
  ok('the backbone primers have no tail', d.rows[2] && d.rows[2][6] === '—', d.rows[2] && d.rows[2][6]);
  ok('and so report no full-length Tm of their own', d.rows[2] && d.rows[2][5] === '—', d.rows[2] && d.rows[2][5]);
  ok('every primer has a Tm', d.rows.every(r => /^\d+\.\d$/.test(r[4])), JSON.stringify(d.rows.map(r => r[4])));
  ok('the reading frame of the fusion is reported', /in frame|out of frame/.test(d.t));
  ok('it is in frame with the vector tag', /in frame/.test(d.t), d.t.match(/.{0,60}frame.{0,60}/)[0]);
  /* Both junctions are drawn and named for what they are rather than labelled
     "junction (vector -> insert)": the panel shows the two strands pairing
     along the 15 bp arm at each end, and says which primer's tail makes each
     pairing. */
  ok('the junctions are shown',
     /Where the insert meets the backbone — 5′ end/.test(d.t) &&
     /Where the insert meets the backbone — 3′ end/.test(d.t),
     (d.t.match(/Where the insert meets[^|]{0,40}/g) || []).join(' / '));
  ok('and each one says which tail makes the pairing',
     /tail on CDK1_F/.test(d.t) && /tail on CDK1_R/.test(d.t),
     (d.t.match(/tail on \w+/g) || []).join(', '));
  ok('the insert is described', /The insert:/.test(d.t));
  ok('it notices the trailing stop codon',
     /ends in a stop codon/.test(d.t), d.t.match(/.{0,40}stop codon.{0,60}/) || '');
  ok('and offers to save', d.save);

  console.log('\n== the design really is correct ==');
  const chk = await page.evaluate(cdk => {
    const P = window.PB, B = window.Bio;
    const d = P.state.records;                       // force evaluation
    const el = document.querySelectorAll('#gbLive table.primers tbody tr');
    const seqOf = i => el[i].children[2].textContent.trim();
    const insF = seqOf(0), insR = seqOf(1), bbF = seqOf(2), bbR = seqOf(3);
    /* The SnapGene file carries its own name (pRSFDuet(tm)-1), so find the
       vector by what it is rather than by the name the test hoped for. */
    const vec = P.state.records.find(r => r.topology === 'circular' && /RSFDuet/i.test(r.name));
    /* the backbone primers must really amplify the vector */
    const bb = B.pcr(vec.seq, bbF, bbR, { circular: true });
    /* the insert primers must really amplify the insert */
    const amp = B.pcr(cdk, insF, insR, {});
    return { bb: bb.length, bbLen: bb.length ? bb[0].length : 0, vecLen: vec.seq.length,
             amp: amp.length, ampLen: amp.length ? amp[0].length : 0, cdkLen: cdk.length,
             insF: insF, insR: insR, bbF: bbF, bbR: bbR };
  }, CDK1);
  ok('the backbone primers amplify the vector', chk.bb >= 1, chk.bb);
  ok('the insert primers amplify the insert', chk.amp >= 1, chk.amp);
  ok('the insert product is the gene plus two 15 bp arms',
     chk.ampLen === chk.cdkLen + 30, chk.ampLen + ' vs ' + (chk.cdkLen + 30));
  /* the arms must match the backbone ends, which is the whole point */
  const arms = await page.evaluate(([insF, insR, bbF, bbR]) => {
    const P = window.PB, B = window.Bio;
    /* The SnapGene file carries its own name (pRSFDuet(tm)-1), so find the
       vector by what it is rather than by the name the test hoped for. */
    const vec = P.state.records.find(r => r.topology === 'circular' && /RSFDuet/i.test(r.name));
    const bb = B.pcr(vec.seq, bbF, bbR, { circular: true })[0];
    return { fwdArm: insF.slice(0, 15), bbEnd: bb.seq.slice(-15),
             revArm: B.revcomp(insR.slice(0, 15)), bbStart: bb.seq.slice(0, 15) };
  }, [chk.insF, chk.insR, chk.bbF, chk.bbR]);
  ok('the forward arm is the backbone\'s 3′ end', arms.fwdArm === arms.bbEnd, arms.fwdArm + ' vs ' + arms.bbEnd);
  ok('the reverse arm is the backbone\'s 5′ end', arms.revArm === arms.bbStart, arms.revArm + ' vs ' + arms.bbStart);

  console.log('\n== changing anything re-designs at once ==');
  const before = await live();
  await page.fill('#gbOv', '25');
  await page.waitForTimeout(500);
  const after = await page.evaluate(() => ({
    t: document.querySelector('#gbLive').textContent,
    tail: document.querySelectorAll('#gbLive table.primers tbody tr')[0].children[6].textContent.trim() }));
  ok('a longer arm is applied with no button', after.tail === '25', after.tail);
  ok('and the product changes with it', after.t !== before);
  await page.fill('#gbOv', '15'); await page.waitForTimeout(400);

  console.log('\n== an out-of-frame placement is called out ==');
  await page.evaluate(() => {
    const sel = document.querySelector('#gbWhere');
    const o = Array.from(sel.options).find(x => /a single base I type/.test(x.textContent));
    if (!o) throw new Error('no "single base" placement in: ' +
      Array.from(sel.options).map(x => x.textContent).join(' | '));
    sel.value = o.value; sel.dispatchEvent(new Event('change', { bubbles: true }));
  });
  await page.waitForTimeout(400);
  /* base 90 sits inside the His-tag reading frame of MCS-1 at a shifted phase */
  await page.fill('#gbFrom', '90');
  await page.waitForTimeout(600);
  const off = await live();
  ok('a frame problem is stated with what to do about it',
     /out of frame by \d/.test(off) && /Shift the opening|add \d/.test(off),
     (off.match(/out of frame.{0,120}/) || [''])[0]);

  console.log('\n== saving ==');
  await page.evaluate(() => {
    /* Pick by the option's own value, not by where it sits in the list. The
       values are indices into the placement list and the options are grouped
       into optgroups in a reading order, so the two do not line up — and
       setting the value to a DOM position silently chooses a different
       placement, which is how this test used to end up designing a BamHI
       digest while believing it had asked for an in-frame fusion. */
    const sel = document.querySelector('#gbWhere');
    const o = Array.from(sel.options).find(x => /in frame/.test(x.textContent));
    if (!o) throw new Error('no in-frame placement offered');
    sel.value = o.value; sel.dispatchEvent(new Event('change', { bubbles: true }));
  });
  await page.waitForTimeout(700);
  const nBefore = await page.evaluate(() => window.PB.state.records.length);
  await page.evaluate(() => document.querySelector('#gbSave').click());
  await page.waitForTimeout(800);
  const saved = await page.evaluate(() => {
    const r = window.PB.active();
    return { n: window.PB.state.records.length, name: r.name, len: r.seq.length, topo: r.topology,
             feats: r.features.map(f => f.label), hist: (r.history || []).map(h => h.op + ' ' + h.detail).join(' | ') };
  });
  ok('the construct is saved and opened', saved.n === nBefore + 1, saved.n);
  /* The vector keeps the name the SnapGene file gives it, trademark symbol
     and all, which is right: renaming someone's vector silently would be
     worse than carrying the symbol through. */
  ok('named after the vector and the insert', /pRSFDuet.{0,2}-1 \+ CDK1/.test(saved.name), saved.name);
  ok('the insert is annotated on the map', saved.feats.indexOf('CDK1') >= 0, JSON.stringify(saved.feats.slice(0, 8)));
  ok('the vector features came across', saved.feats.some(f => /His/.test(f || '')), JSON.stringify(saved.feats.slice(0, 10)));
  ok('the history records the primers', /insert primers/.test(saved.hist) && /backbone primers/.test(saved.hist));
  const prot = await page.evaluate(() => {
    const r = window.PB.active(), B = window.Bio;
    const f = r.features.find(x => x.label === 'CDK1');
    if (!f) return null;
    return { type: f.type, aa: f.quals && f.quals.translation ? f.quals.translation.slice(0, 10) : null };
  });
  ok('and annotates it as coding with its translation', prot && prot.type === 'CDS' && prot.aa === 'MEDYTKIEKI', JSON.stringify(prot));

  console.log('\n== primers go to the library ==');
  await openDlg();
  await page.waitForTimeout(400);
  await page.selectOption('#gbIns_src', ins);
  await page.waitForTimeout(600);
  const pBefore = await page.evaluate(() => window.PB.state.primers.length);
  const hasSave = await page.evaluate(() => !!document.querySelector('#gbSavePrimers'));
  if (hasSave) {
    await page.evaluate(() => document.querySelector('#gbSavePrimers').click());
    await page.waitForTimeout(400);
    const pAfter = await page.evaluate(() => window.PB.state.primers.map(p => p.name));
    ok('all four primers are stored', pAfter.length === pBefore + 4, pAfter.length - pBefore);
    ok('named after the insert', pAfter.some(n => /CDK1.*_F$/.test(n)) && pAfter.some(n => /bb_F$/.test(n)),
       JSON.stringify(pAfter.slice(-4)));
  } else ok('all four primers are stored', false, 'no save button');

  console.log('\n== a primer with an overhang keeps its annealing Tm out of the designer ==');
  /* design again from scratch, so this section does not depend on whatever
     the last one left in the dialog */
  await openDlg();
  await page.waitForTimeout(300);
  await page.selectOption('#gbIns_src', ins);
  await page.waitForTimeout(800);
  /* The whole point of the two columns in the designer is that the first
     cycles anneal on the 3' part alone. A primer saved to the library and
     then shown with its full-length Tm would undo that silently: the number
     is higher, it looks like an ordinary good primer, and a block set to it
     anneals above what actually holds. */
  const inDesigner = await page.evaluate(() =>
    Array.from(document.querySelectorAll('#gbLive table.primers tbody tr'))
      .map(tr => Array.from(tr.children).map(td => td.textContent.trim())));
  const tailed = inDesigner.filter(r => /^\d+$/.test(r[6]));
  ok('the designer made primers with overhangs', tailed.length >= 2, JSON.stringify(inDesigner.map(r => r[6])));
  ok('whose two Tms genuinely differ',
     tailed.every(r => Math.abs(+r[4] - +r[5]) > 1), JSON.stringify(tailed.map(r => r[4] + '/' + r[5])));

  await page.evaluate(() => {
    const b = Array.from(document.querySelectorAll('#modalRoot button'))
      .find(x => /primer/i.test(x.textContent) && /save|library/i.test(x.textContent));
    if (b) b.click();
  });
  await page.waitForTimeout(400);
  const stored = await page.evaluate(() => window.PB.state.primers.map(p =>
    ({ name: p.name, len: p.seq.length, anneal: (p.anneal || '').length, tail: (p.tail || '').length,
       suffix: !p.anneal || p.seq.slice(p.seq.length - p.anneal.length) === p.anneal })));
  ok('saving carries the annealing portion with the primer',
     stored.some(x => x.anneal > 0), JSON.stringify(stored));
  ok('and it really is the 3′ end of what was saved',
     stored.every(x => x.suffix), JSON.stringify(stored.filter(x => !x.suffix)));
  ok('a primer with no overhang stores none',
     stored.some(x => x.anneal === 0 && /bb/.test(x.name)), JSON.stringify(stored));

  await page.evaluate(() => { window.PB.closeModal(); window.PB.state.view = 'primers'; window.PB.renderAll(); });
  await page.waitForTimeout(600);
  const lib = await page.evaluate(() => Array.from(document.querySelectorAll('#view-primers table.tbl tbody tr'))
    .map(tr => {
      /* the cell is a button now, so that the split can be set by hand;
         its first line is the number, the rest are the notes under it */
      const c = Array.from(tr.children);
      const btn = c[3].querySelector('button');
      const lines = (btn ? btn.innerText : c[3].textContent).split('\n')
        .map(x => x.trim()).filter(Boolean);
      return { name: (tr.querySelector('input') || {}).value,
               tm: lines[0] || '',
               full: lines.find(x => /^full /.test(x)) || '',
               title: c[3].getAttribute('title') || '' };
    }));
  const libF = lib.find(r => /_F$/.test(r.name || '') && !/bb/.test(r.name || ''));
  const desF = tailed[0];
  ok('the library shows the annealing Tm, the same number the designer did',
     libF && desF && libF.tm === desF[4], (libF && libF.tm) + ' vs ' + (desF && desF[4]));
  ok('and shows the whole primer’s Tm underneath rather than dropping it',
     libF && libF.full === 'full ' + desF[5], libF && libF.full);
  ok('and explains which is which on hover',
     libF && /anneal/i.test(libF.title) && /overhang/i.test(libF.title), libF && libF.title.slice(0, 80));
  const libBB = lib.find(r => /bb/.test(r.name || ''));
  ok('a primer with no overhang shows one number and no second line',
     libBB && libBB.full === '', libBB && JSON.stringify(libBB));

  const sheet = await page.evaluate(() => {
    const rows = window.PB.primerExportRows ? window.PB.primerExportRows('all') : null;
    return rows ? rows.map(r => ({ name: r.name, len: r.len, tm: r.tm, tmFull: r.tmFull,
      annealLen: r.annealLen, tailLen: r.tailLen, basis: r.tmBasis })) : null;
  });
  if (sheet) {
    const ex = sheet.find(r => /_F$/.test(r.name) && !/bb/.test(r.name));
    ok('the export carries both temperatures, not just one',
       ex && ex.tm != null && ex.tmFull != null && Math.abs(ex.tmFull - ex.tm) > 1,
       JSON.stringify(ex));
    ok('with the annealing length beside them',
       ex && ex.annealLen > 0 && ex.tailLen > 0 && ex.annealLen + ex.tailLen === ex.len,
       JSON.stringify(ex));
    ok('and says where the split came from',
       ex && /as designed/.test(ex.basis || ''), ex && ex.basis);
  } else ok('the export rows are reachable', false, 'PB.primerExportRows missing');

  console.log('\n== a primer that could not reach the target says so ==');
  /* The annealing length has a floor, so on a GC-rich end even the shortest
     allowed stretch can melt well above the target. The designer still
     hands over the best it found — that is right — but a silent 76 reads
     as a primer that met the 60 that was asked for. */
  await openDlg();
  await page.waitForTimeout(300);
  await page.selectOption('#gbIns_src', '');
  await page.fill('#gbIns_txt',
    'ATGGCGGCCGCGGCGGCCGCGGCGGCGGCCGCGGCC' +
    'ACGATTAAACTGAATGGTGAAACCTATAAAGATATTAAAGATCTGAATGATAAAGAATATAAACTG' +
    'AATTTAAATAAATTAAATTTAAATAAATTAAATTTATAA');
  await page.fill('#gbInsName', 'Lopsided');
  await page.waitForTimeout(1000);
  const lop = await page.evaluate(() => {
    const live = document.querySelector('#gbLive');
    return {
      target: +document.querySelector('#gbAn').value,
      rows: Array.from(live.querySelectorAll('table.primers tbody tr'))
        .map(tr => Array.from(tr.children).map(td => td.textContent.trim())),
      badges: Array.from(live.querySelectorAll('table.primers .badge')).map(b => b.textContent.trim()),
      pair: Array.from(live.querySelectorAll('.note')).map(n => n.textContent.replace(/\s+/g, ' '))
        .filter(n => /apart/.test(n))
    };
  });
  ok('the target was 60', lop.target === 60, lop.target);
  ok('a primer well above it is flagged on its own row',
     lop.badges.some(b => /^\+1?\d+\.\d off 60/.test(b)), JSON.stringify(lop.badges));
  ok('and one well below it too',
     lop.badges.some(b => /^−9\.\d off 60/.test(b) || /^−\d+\.\d off 60/.test(b)),
     JSON.stringify(lop.badges));
  /* The backbone pair at this opening cannot reach 60 either — its
     shortest allowed length already melts above it — so both are flagged,
     and the two are balanced against each other rather than one of them
     being left at the target and the other stranded above it. */
  const bbRows = lop.rows.filter(r => /Backbone/.test(r[0]));
  ok('the backbone pair is balanced to within a couple of degrees',
     bbRows.length === 2 && Math.abs(parseFloat(bbRows[0][4]) - parseFloat(bbRows[1][4])) <= 2,
     JSON.stringify(bbRows.map(r => r[4])));
  ok('and both say how far off the target they are, rather than one of them hiding it',
     bbRows.every(r => /off 60/.test(r[4])), JSON.stringify(bbRows.map(r => r[4])));
  ok('the pair being far apart is called out separately',
     lop.pair.length === 1 && /insert pair are 2\d\.\d °C apart/.test(lop.pair[0]),
     JSON.stringify(lop.pair).slice(0, 140));
  ok('and it says which way to compromise rather than only that it is wrong',
     lop.pair[0] && /is the usual compromise/.test(lop.pair[0]) && /Other lengths/.test(lop.pair[0]),
     (lop.pair[0] || '').slice(-150));

  console.log('\n== a pair that is fine is left in peace ==');
  await openDlg();
  await page.waitForTimeout(300);
  await page.selectOption('#gbIns_src', ins);
  await page.waitForTimeout(900);
  const calm = await page.evaluate(() => {
    const live = document.querySelector('#gbLive');
    return { badges: Array.from(live.querySelectorAll('table.primers .badge')).map(b => b.textContent.trim()),
             rows: Array.from(live.querySelectorAll('table.primers tbody tr'))
               .map(tr => Array.from(tr.children).map(td => td.textContent.trim())),
             pair: Array.from(live.querySelectorAll('.note')).map(n => n.textContent).filter(n => /apart/.test(n)).length };
  });
  ok('the insert pair it could balance carries no badge',
     !calm.rows.filter(r => /Insert/.test(r[0])).some(r => /off 60/.test(r[4])),
     JSON.stringify(calm.rows.filter(r => /Insert/.test(r[0])).map(r => r[4])));
  ok('and no pair is far enough apart to warn about', calm.pair === 0, calm.pair);

  console.log('\n== one backbone pair across a series of inserts ==');
  /* The workflow this is for: one backbone, one opening, a different gene
     each time. The backbone primers are ordered once, so every construct
     after the first has to be designed against that same pair — and the
     arms on its insert primers have to match them. */
  const INS2 = 'ATGCGTCGCGTTAGTCAGGCTAGCGTGAAACCCTCAGGTTTAACCGGTAAAGAATTCGCGCGCTTATTTAAACCCTAA';
  let chosenLabel = '';
  const designWith = async (seq, nm, tm, reuse) => {
    await openDlg();
    await page.waitForTimeout(300);
    await page.selectOption('#gbIns_src', '');
    await page.fill('#gbIns_txt', seq);
    await page.fill('#gbInsName', nm);
    await page.fill('#gbAn', String(tm));
    await page.waitForTimeout(900);
    if (reuse) {
      const idx = await page.evaluate(() => {
        const sel = document.querySelector('#gbBB');
        const i = Array.from(sel.options).findIndex(o => /^reuse /.test(o.textContent));
        if (i < 0) return -1;
        sel.value = sel.options[i].value;
        sel.dispatchEvent(new Event('change', { bubbles: true }));
        return i;
      });
      if (idx < 0) return { none: true };
      await page.waitForTimeout(900);
      chosenLabel = await page.evaluate(() => {
        const s = document.querySelector('#gbBB');
        return s.options[s.selectedIndex].textContent;
      });
    }
    return page.evaluate(() => ({
      rows: Array.from(document.querySelectorAll('#gbLive table.primers tbody tr'))
        .map(tr => Array.from(tr.children).map(td => td.textContent.trim())),
      names: Array.from(document.querySelectorAll('#gbLive table.primers tbody tr'))
        .map(tr => [tr.children[0].textContent.trim(),
                    (tr.querySelector('input') || {}).value || tr.children[1].textContent.trim()]),
      live: (document.querySelector('#gbLive') || {}).textContent || '',
      note: (document.querySelector('#gbBBnote') || {}).textContent || '',
      options: Array.from(document.querySelectorAll('#gbBB option')).map(o => o.textContent)
    }));
  };

  const serie1 = await designWith(INS2, 'SeriesOne', 60, false);
  ok('a fresh design is always the first choice',
     /design a new pair/.test(serie1.options[0]), JSON.stringify(serie1.options.slice(0, 2)));
  ok('and each saved pair is offered once, not once per copy of it in the library',
     new Set(serie1.options).size === serie1.options.length,
     JSON.stringify(serie1.options));
  ok('while a fresh design warns that the pair depends on the Tm setting',
     /changing that between/.test(serie1.note), serie1.note.slice(0, 110));
  const bbFirst = serie1.rows.filter(r => /Backbone/.test(r[0])).map(r => r[2]);
  ok('a backbone pair was designed', bbFirst.length === 2, JSON.stringify(bbFirst));

  await page.evaluate(() => {
    const b = Array.from(document.querySelectorAll('#modalRoot button'))
      .find(x => /primer/i.test(x.textContent) && /save|library/i.test(x.textContent));
    if (b) b.click();
  });
  await page.waitForTimeout(400);

  const second = await designWith('ATGGGTTTTCTGCGTTTTCAGGCTAGCGTGCCCTCAGGTTTAACCGGTAAAGAATTCCGTACGTTACCCGGGAAATAA',
                                  'SeriesTwo', 66, true);
  ok('the saved pair is offered for the next insert', !second.none,
     JSON.stringify((second.options || []).slice(0, 3)));
  if (!second.none) {
    const bbSecond = second.rows.filter(r => /Backbone/.test(r[0])).map(r => r[2]);
    ok('and choosing it gives back exactly those two primers',
       JSON.stringify(bbSecond) === JSON.stringify(bbFirst),
       JSON.stringify(bbSecond) + ' vs ' + JSON.stringify(bbFirst));
    ok('even though the anneal Tm was changed between the two designs',
       /66/.test(await page.evaluate(() => document.querySelector('#gbAn').value)));
    ok('the panel says the primers were reused rather than chosen',
       /primers you already have/.test(second.live), (second.live.match(/.{0,30}already have.{0,40}/) || [''])[0]);
    ok('and the note explains the opening now comes from them',
       /opening comes from these two primers/.test(second.note), second.note.slice(0, 90));
    /* the thing that makes it usable: the arms must match the reused pair */
    const arms = await page.evaluate(() => {
      const rows = Array.from(document.querySelectorAll('#gbLive table.primers tbody tr'))
        .map(tr => Array.from(tr.children).map(td => td.textContent.trim()));
      const ins = rows.filter(r => /Insert/.test(r[0])).map(r => r[2]);
      const bb = rows.filter(r => /Backbone/.test(r[0])).map(r => r[2]);
      const tail = parseInt(rows.filter(r => /Insert/.test(r[0]))[0][6], 10);
      return { fArm: ins[0].slice(0, tail), rArm: ins[1].slice(0, tail), bbF: bb[0], bbR: bb[1], tail: tail };
    });
    const rc = x => x.split('').reverse().map(c => ({ A: 'T', T: 'A', G: 'C', C: 'G' })[c]).join('');
    ok('the forward arm is the end of the reused reverse primer',
       arms.fArm === rc(arms.bbR).slice(-arms.tail), arms.fArm + ' vs ' + rc(arms.bbR).slice(-arms.tail));
    ok('and the reverse arm is the start of the reused forward primer',
       arms.rArm === rc(arms.bbF.slice(0, arms.tail)), arms.rArm + ' vs ' + rc(arms.bbF.slice(0, arms.tail)));

    /* The tube in the freezer is labelled after the construct it was made
       for. Relabelling the same two primers after today's insert would
       mean two names for one tube, and the wrong one pulled. */
    const nameOf = (set, role) => (set.names.find(n => n[0] === role) || [])[1];
    ok('the reused pair carries the names the library holds it under',
       chosenLabel.indexOf(nameOf(second, 'Backbone forward')) >= 0 &&
       chosenLabel.indexOf(nameOf(second, 'Backbone reverse')) >= 0,
       nameOf(second, 'Backbone forward') + ' / ' + nameOf(second, 'Backbone reverse') +
       '  offered as  ' + chosenLabel);
    ok('which is what the tube in the freezer says, not a new label',
       await page.evaluate(n => window.PB.state.primers.some(p => p.name === n),
                           nameOf(second, 'Backbone forward')),
       nameOf(second, 'Backbone forward'));
    ok('rather than being relabelled after the insert being cloned today',
       !/SeriesTwo/.test(String(nameOf(second, 'Backbone forward'))),
       nameOf(second, 'Backbone forward'));
    ok('while the insert primers are named after this insert',
       /SeriesTwo/.test(String(nameOf(second, 'Insert forward'))),
       nameOf(second, 'Insert forward'));

    /* and saving the second design must not put a second copy of the same
       two primers in the library */
    const lib = await page.evaluate(() => {
      const before = window.PB.state.primers.length;
      const b = Array.from(document.querySelectorAll('#modalRoot button'))
        .find(x => /primer/i.test(x.textContent) && /save|library/i.test(x.textContent));
      if (b) b.click();
      const seqs = window.PB.state.primers.map(p => p.name + '|' + p.seq);
      return { before: before, after: window.PB.state.primers.length,
               dup: seqs.length - new Set(seqs).size,
               toast: (document.querySelector('.toast') || {}).textContent || '' };
    });
    ok('saving the second design adds only the two insert primers',
       lib.after - lib.before === 2, lib.before + ' -> ' + lib.after);
    ok('so the library holds no duplicate of the backbone pair', lib.dup === 0, lib.dup);
    ok('and it says so instead of claiming four were saved',
       /2 primers/.test(lib.toast) && /already there/.test(lib.toast), lib.toast);

    /* The trap this closes: with a pair chosen, the placement dropdown no
       longer decides anything. Point it somewhere else and the product is
       still the pair's — which has to be said, not left to be discovered
       on a gel. */
    const clash = await page.evaluate(() => {
      const sel = document.querySelector('#gbWhere');
      const opts = Array.from(sel.options);
      /* any placement other than the one the pair happens to match */
      const cur = sel.value;
      const other = opts.find(o => o.value !== cur && !o.disabled && /\d/.test(o.textContent));
      if (!other) return { none: true };
      sel.value = other.value; sel.dispatchEvent(new Event('change', { bubbles: true }));
      return { none: false, chose: other.textContent };
    });
    if (!clash.none) {
      await page.waitForTimeout(900);
      const warned = await page.evaluate(() => {
        const t = (document.querySelector('#gbLive') || {}).textContent || '';
        return { warn: /not what these primers do/.test(t),
                 said: (t.match(/.{0,20}primers won.{0,60}/) || [''])[0],
                 still: Array.from(document.querySelectorAll('#gbLive table.primers tbody tr'))
                   .filter(tr => /Backbone/.test(tr.children[0].textContent))
                   .map(tr => tr.children[2].textContent.trim()) };
      });
      ok('pointing the placement elsewhere while a pair is chosen is called out',
         warned.warn, warned.said || 'no warning');
      ok('and the primers still win rather than silently changing',
         JSON.stringify(warned.still) === JSON.stringify(bbFirst),
         JSON.stringify(warned.still));
    }
  }

  await page.evaluate(() => { const m = document.querySelector('#modalRoot'); if (m && !m.hidden) window.PB.closeModal(); });
  ok('no errors anywhere in the run', errors.length === 0, errors.slice(0, 3).join(' | '));

  console.log('\n   ' + pass + ' passed, ' + fail + ' failed\n');
  await browser.close();
  process.exit(fail ? 1 : 0);
})();
