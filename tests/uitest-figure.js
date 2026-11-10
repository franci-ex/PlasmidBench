/* ===========================================================================
   The alignment figure, through the Align tab as a person would reach it.

   test16.js checks the drawing. This checks the way in: that a protein
   alignment survives the FASTA box at all (it used not to — the parser
   kept only nucleotide codes, so a peptide arrived as the handful of its
   residues that happen to be IUPAC bases), that an already-aligned file is
   not silently re-aligned with a worse aligner, and that the panel offers
   the right shading modes for the kind of sequence in front of it.
   ========================================================================= */
const { chromium } = require('playwright');
const path = require('path');
let pass = 0, fail = 0;
function ok(name, cond, extra) {
  if (cond) { pass++; console.log('  PASS  ' + name); }
  else { fail++; console.log('  FAIL  ' + name + (extra !== undefined ? '   -> ' + extra : '')); }
}

const PROT = '>human\nMKVLAAGIVGLNLGGKLYEQ--WATRDF\n>mouse\nMKVLAAGIVGLNLGSKLYEQ--WATRDF\n' +
             '>yeast\nMKILSAGVVGMNIGGRLFDQKKWSTRDY\n>fly\nMRVLAAGIIGLNLGGKLYEQ--WATRNF\n';
const DNA = '>pET28a\nATGGCTAGCAAAGGAGAAGAACTTTTCACTGGAGTTGTC--CCAATTCTTGTTGAATTAGATGGT\n' +
            '>pET21b\nATGGCTAGCAAAGGGGAAGAACTCTTCACTGGAGTTGTC--CCAATTCTTGTCGAATTAGATGGT\n' +
            '>pRSF\nATGGCCAGTAAAGGAGAAGAACTTTTCACCGGAGTCGTGAACCAATTCTTGTTGAACTAGATGGT\n' +
            '>pACYC\nATGGCTAGCAAAGGAGAAGAGCTTTTCACTGGAGTTGTC--CCAATACTTGTTGAATTAGATGGT\n';
const UNALIGNED_PROT = '>a\nMKVLAAGIVGLNLGGKLYEQWATRDF\n>b\nMKVLAAGIVGLNLGSKLYEQWATRDF\n>c\nMKILSAGVVGMNIGGRLFDQKKWSTRDY\n';

(async () => {
  const browser = await chromium.launch();
  const page = await browser.newPage({ viewport: { width: 1500, height: 1150 } });
  const errors = [];
  page.on('pageerror', e => errors.push('pageerror: ' + String(e)));
  page.on('console', m => { if (m.type() === 'error' && !/Failed to load resource/.test(m.text())) errors.push('console: ' + m.text()); });
  await page.goto('file://' + path.resolve('plasmidbench.html'));
  await page.waitForTimeout(800);

  const align = async (text) => {
    await page.evaluate(() => { window.PB.state.view = 'align'; window.PB.state.alignMode = 'msa'; window.PB.renderAll(); });
    await page.waitForTimeout(250);
    await page.evaluate(t => { document.querySelector('#msaText').value = t; }, text);
    await page.click('#msaGo');
    await page.waitForTimeout(500);
    await page.evaluate(() => { const d = document.querySelector('#figDet'); if (d) d.open = true; });
    await page.waitForTimeout(250);
  };
  const result = () => page.evaluate(() => {
    const r = window.PB.state.alignState.msaResult;
    return r && r.ok ? { method: r.method, n: r.rows.length, len: r.length, seqs: r.rows.map(x => x.seq) } : null;
  });

  /* ------------------------------------------------ protein gets through */
  console.log('\n== a protein alignment survives the box ==');
  await align(PROT);
  let r = await result();
  ok('it aligns at all', !!r, 'no result');
  ok('all four sequences arrive', r && r.n === 4, r && r.n);
  ok('the residues are not filtered down to the IUPAC bases among them',
    r && r.seqs[0] === 'MKVLAAGIVGLNLGGKLYEQ--WATRDF', r && r.seqs[0]);
  ok('the gaps it came with are kept', r && /-/.test(r.seqs[0]));
  ok('and it is not re-aligned behind your back', r && /already aligned/.test(r.method), r && r.method);
  ok('which it says, naming the molecule', r && /protein/.test(r.method), r && r.method);

  console.log('\n== the figure panel, for protein ==');
  const modes = () => page.evaluate(() => Array.from(document.querySelectorAll('#figMode option')).map(o => o.value));
  ok('the protein groupings are offered', (await modes()).join(',') === 'physchem,multalin,strict,blosum62', (await modes()).join(','));
  ok('the nucleotide grouping is not', (await modes()).indexOf('purine') === -1);
  ok('no purine colouring checkbox on a protein figure',
    await page.evaluate(() => !document.querySelector('#figNuc')));
  ok('the structure box asks for DSSP letters', await page.evaluate(() =>
    /DSSP/.test(document.querySelector('.fig-wrap').textContent)));
  ok('a figure is drawn straight away', await page.evaluate(() => !!document.querySelector('#figPrev img')));

  console.log('\n== the controls change the drawing ==');
  const svgOf = () => page.evaluate(() => decodeURIComponent((document.querySelector('#figPrev img') || {}).src || ''));
  const before = await svgOf();
  await page.evaluate(() => {
    const s = document.querySelector('#figSS');
    s.value = 'HHHHHHHHHEEEEEEETTTT  EEEEEE';
    s.dispatchEvent(new Event('input'));
  });
  await page.waitForTimeout(350);
  const withSS = await svgOf();
  ok('a structure string adds a track', withSS !== before && withSS.length > before.length,
    JSON.stringify({ before: before.length, after: withSS.length }));
  await page.evaluate(() => {
    const m = document.querySelector('#figMono'); m.checked = true; m.dispatchEvent(new Event('change'));
  });
  await page.waitForTimeout(350);
  const mono = await svgOf();
  ok('monochrome takes the red out', mono.indexOf('#e02020') === -1 && withSS.indexOf('#e02020') !== -1);
  ok('and leaves the alignment in', /human/.test(mono) && /yeast/.test(mono));
  await page.evaluate(() => {
    const m = document.querySelector('#figMono'); m.checked = false; m.dispatchEvent(new Event('change'));
    const t = document.querySelector('#figThr'); t.value = '0.95'; t.dispatchEvent(new Event('input'));
  });
  await page.waitForTimeout(350);
  ok('the threshold is honoured', await page.evaluate(() => window.PB.state.alignState.fig.threshold === 0.95),
    await page.evaluate(() => window.PB.state.alignState.fig.threshold));

  /* ------------------------------------------------------- DNA and RNA */
  console.log('\n== the same panel, for DNA ==');
  await align(DNA);
  r = await result();
  ok('a nucleotide alignment is recognised as one', r && /nucleotide/.test(r.method), r && r.method);
  ok('the nucleotide groupings are offered', (await modes()).join(',') === 'strict,purine', (await modes()).join(','));
  ok('the protein groupings are not', (await modes()).indexOf('physchem') === -1);
  ok('the purine colouring checkbox appears', await page.evaluate(() => !!document.querySelector('#figNuc')));
  ok('the structure box now asks for dot-bracket', await page.evaluate(() =>
    /dot-bracket/.test(document.querySelector('.fig-wrap').textContent)));
  ok('all four molecule choices are there', await page.evaluate(() =>
    Array.from(document.querySelectorAll('#figMol option')).map(o => o.value).join(',') === 'auto,protein,dna,rna'),
    await page.evaluate(() => Array.from(document.querySelectorAll('#figMol option')).map(o => o.value).join(',')));
  ok('auto says what it decided', await page.evaluate(() =>
    /read as DNA/.test(document.querySelector('#figMol').textContent)));

  await page.evaluate(() => {
    /* the threshold is a preference and carries over from the protein
       figure above; at 0.95 only strictly conserved columns shade, and
       those are a red box rather than a coloured letter */
    const t = document.querySelector('#figThr'); t.value = '0.7'; t.dispatchEvent(new Event('input'));
    const n = document.querySelector('#figNuc'); n.checked = true; n.dispatchEvent(new Event('change'));
    const s = document.querySelector('#figSS'); s.value = '((((((((....))))))))'; s.dispatchEvent(new Event('input'));
  });
  await page.waitForTimeout(400);
  const dnaSvg = await svgOf();
  ok('a dot-bracket structure draws stems', dnaSvg.indexOf('<path') !== -1);
  ok('purines and pyrimidines get different colours',
    dnaSvg.indexOf('#c02626') !== -1 && dnaSvg.indexOf('#2b4fd8') !== -1);

  await page.evaluate(() => { const m = document.querySelector('#figMol'); m.value = 'rna'; m.dispatchEvent(new Event('change')); });
  await page.waitForTimeout(600);
  const rnaSvg = await svgOf();
  ok('switching to RNA shows uracil', rnaSvg.indexOf('>U<') !== -1 && rnaSvg.indexOf('>T<') === -1);
  ok('and the panel stayed open through the rebuild',
    await page.evaluate(() => document.querySelector('#figDet').open));

  /* -------------------------------------- genuinely unaligned sequences */
  console.log('\n== sequences that are NOT already aligned ==');
  await align(UNALIGNED_PROT);
  r = await result();
  ok('are aligned here', r && /centre-star|center-star/i.test(r.method), r && r.method);
  ok('with BLOSUM62 rather than a flat match score', r && /BLOSUM62/.test(r.method), r && r.method);
  ok('and come out padded to one width', r && r.seqs.every(s => s.length === r.len), r && r.seqs.map(s => s.length).join(','));

  /* --------------------------------------------------------- the export */
  console.log('\n== saving the figure ==');
  await page.evaluate(() => {
    window.__saved = null;
    window.PB.offerFile = async (fn, text) => { window.__saved = { fn: fn, n: text.length, text: text }; return true; };
    window.PB.offerBlob = async (fn, blob) => { window.__saved = { fn: fn, n: blob.size }; return true; };
  });
  await page.evaluate(() => { const d = document.querySelector('#figDet'); if (d) d.open = true; });
  await page.waitForTimeout(200);
  await page.click('#figSvg');
  await page.waitForTimeout(400);
  let saved = await page.evaluate(() => window.__saved);
  ok('Save SVG hands over an SVG', saved && /\.svg$/.test(saved.fn) && saved.n > 1000, JSON.stringify(saved && saved.fn));
  ok('and it is the figure, not the page', saved && saved.text.indexOf('<svg') !== -1 && /MKVL|>M</.test(saved.text));
  await page.click('#figPng');
  await page.waitForTimeout(1600);
  saved = await page.evaluate(() => window.__saved);
  ok('Save PNG hands over a PNG', saved && /\.png$/.test(saved.fn) && saved.n > 2000, JSON.stringify(saved && [saved.fn, saved.n]));

  ok('nothing threw anywhere in the run', errors.length === 0, errors.slice(0, 3).join(' | '));
  console.log('\n   ' + pass + ' passed, ' + fail + ' failed\n');
  await browser.close();
  process.exit(fail ? 1 : 0);
})().catch(e => { console.error('HARNESS ERROR', e); process.exit(2); });
