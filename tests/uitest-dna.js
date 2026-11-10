const { chromium } = require('playwright');
const path = require('path'), fs = require('fs');
let pass = 0, fail = 0;
function ok(name, cond, extra) {
  if (cond) { pass++; console.log('  PASS  ' + name); }
  else { fail++; console.log('  FAIL  ' + name + (extra !== undefined ? '   -> ' + extra : '')); }
}
function prng(seed) { let x = seed; return () => (x = (x * 1103515245 + 12345) & 0x7fffffff) / 0x7fffffff; }
function randDNA(n, seed) { const r = prng(seed); let s = ''; for (let i = 0; i < n; i++) s += 'ACGT'[Math.floor(r() * 4)]; return s; }
function packet(type, payload) {
  const h = Buffer.alloc(5); h.writeUInt8(type, 0); h.writeUInt32BE(payload.length, 1);
  return Buffer.concat([h, payload]);
}
function buildDna(seq, circular, featuresXml, primersXml, notesXml) {
  const header = Buffer.alloc(14);
  header.write('SnapGene', 0, 'latin1');
  header.writeUInt16BE(1, 8); header.writeUInt16BE(14, 10); header.writeUInt16BE(14, 12);
  const dna = Buffer.concat([Buffer.from([circular ? 1 : 0]), Buffer.from(seq, 'latin1')]);
  const parts = [packet(9, header), packet(0, dna)];
  if (notesXml) parts.push(packet(6, Buffer.from(notesXml, 'utf8')));
  if (primersXml) parts.push(packet(5, Buffer.from(primersXml, 'utf8')));
  if (featuresXml) parts.push(packet(10, Buffer.from(featuresXml, 'utf8')));
  return Buffer.concat(parts);
}
function wrap(body) {
  return '<!doctype html><html><head><meta charset="utf8">' +
    '<meta name="viewport" content="width=device-width,initial-scale=1,viewport-fit=cover">' +
    '<style>:root{color-scheme:light;padding-top:env(safe-area-inset-top,0px);padding-bottom:env(safe-area-inset-bottom,0px)}' +
    'body{margin:0;font:14px system-ui;background:#fafaf9}img{max-width:100%}[hidden]{display:none!important}</style>' +
    '</head><body>' + body + '</body></html>';
}

(async () => {
  fs.writeFileSync(path.join(__dirname, 'preview.html'),
    wrap(fs.readFileSync(path.join(__dirname, '..', 'plasmidbench.html'), 'utf8')));

  const seq = randDNA(5400, 4242);
  const featuresXml = '<Features nextValidID="3">' +
    '<Feature recentID="0" name="AmpR" type="CDS" directionality="1">' +
      '<Segment range="401-1261" color="#ccffcc" type="standard"/>' +
      '<Q name="label"><V text="AmpR (bla)"/></Q></Feature>' +
    '<Feature recentID="1" name="ori" type="rep_origin" directionality="1">' +
      '<Segment range="1500-2088" color="#ffef86" type="standard"/></Feature>' +
    '<Feature recentID="2" name="EGFP" type="CDS" directionality="2">' +
      '<Segment range="2500-3217" color="#b1ff67" type="standard"/></Feature>' +
    '</Features>';
  const primersXml = '<Primers><Primer name="M13 fwd" sequence="GTAAAACGACGGCCAGT"/></Primers>';
  const notesXml = '<Notes><CustomMapLabel>pMyPlasmid</CustomMapLabel></Notes>';
  const dna = buildDna(seq, true, featuresXml, primersXml, notesXml);

  const browser = await chromium.launch();
  const page = await browser.newPage({ viewport: { width: 1440, height: 900 } });
  const errors = [];
  const ignorable = t => /ERR_TUNNEL_CONNECTION_FAILED|ERR_NAME_NOT_RESOLVED|fonts\.(googleapis|gstatic)\.com|Failed to load resource/.test(t);
  page.on('console', m => { if (m.type() === 'error' && !ignorable(m.text())) errors.push(m.text()); });
  page.on('pageerror', e => errors.push('PAGEERROR: ' + e.message));
  await page.goto('file://' + path.join(__dirname, 'preview.html'));
  await page.waitForTimeout(800);

  console.log('\n== drag-and-drop a SnapGene .dna file ==');
  // exercise the real handleFiles path with a genuine File object
  await page.evaluate(bytes => {
    const u8 = new Uint8Array(bytes);
    const f = new File([u8], 'pMyPlasmid.dna', { type: 'application/octet-stream' });
    window.PB.handleFiles([f]);
  }, Array.from(dna));
  await page.waitForTimeout(900);

  const rec = await page.evaluate(() => {
    const r = window.PB.active();
    return { name: r.name, len: r.seq.length, topo: r.topology, seq: r.seq,
             feats: r.features.map(f => ({ l: f.label, s: f.start, e: f.end, st: f.strand, c: f.color, t: f.type })) };
  });
  ok('file imported, not skipped', rec.len > 0, rec.len);
  ok('sequence length is exactly 5,400 bp (no byte-soup padding)', rec.len === 5400, rec.len);
  ok('sequence matches the source exactly', rec.seq === seq, rec.seq.slice(0, 24));
  ok('name taken from the map label, not "pasted sequence"', rec.name === 'pMyPlasmid', rec.name);
  ok('topology read as circular', rec.topo === 'circular', rec.topo);
  ok('features imported', rec.feats.length >= 3, rec.feats.length);
  const amp = rec.feats.find(f => /AmpR/.test(f.l));
  ok('AmpR label from the label qualifier', amp && amp.l === 'AmpR (bla)', amp && amp.l);
  ok('AmpR coordinates correct', amp && amp.s === 400 && amp.e === 1261, amp && (amp.s + '-' + amp.e));
  ok('AmpR typed as CDS', amp && amp.t === 'CDS', amp && amp.t);
  ok('SnapGene colour preserved', amp && amp.c === '#ccffcc', amp && amp.c);
  const egfp = rec.feats.find(f => f.l === 'EGFP');
  ok('reverse-strand feature kept its strand', egfp && egfp.st === -1, egfp && egfp.st);
  const ori = rec.feats.find(f => f.l === 'ori');
  ok('origin feature imported', ori && ori.t === 'rep_origin', ori && ori.t);
  const prim = await page.evaluate(() => window.PB.state.primers.map(p => p.name + ':' + p.seq));
  ok('primer pulled out of the .dna file', prim.some(p => /M13 fwd:GTAAAACGACGGCCAGT/.test(p)), prim.join(','));

  console.log('\n== the map reflects it ==');
  await page.click('[data-v="map"]');
  await page.waitForTimeout(700);
  const mode = await page.evaluate(() => window.PB.state.mapMode);
  ok('circular plasmid shows the circular projection', mode === 'circular', mode);
  const featG = await page.locator('#view-map svg.map g.feat').count();
  ok('features drawn on the map', featG >= 3, featG);
  const centre = await page.locator('#view-map .ctrname').textContent();
  ok('map centre shows the real name', centre === 'pMyPlasmid', centre);
  const meta = await page.locator('#view-map .ctrmeta').first().textContent();
  ok('map centre shows 5,400 bp circular', /5,400 bp - circular/.test(meta), meta);
  // the imported colour must actually be used
  const strokes = await page.evaluate(() =>
    Array.from(document.querySelectorAll('#view-map svg.map g.feat path[stroke]')).map(p => p.getAttribute('stroke')));
  ok('imported colour used for the arc', strokes.some(s => s.toLowerCase() === '#ccffcc'), strokes.join(','));

  /* Zoom state is kept per construct in PB.mapView(rec), not in one global
     state.mapZoom, so that switching constructs does not drag the last one's
     zoom along. And on the circle a plain wheel turns the plasmid: it is a
     pinch — a wheel event with ctrlKey set — that zooms. */
  console.log('\n== zoom and pan ==');
  const svgBox = await page.locator('#view-map svg.map').boundingBox();
  const pinch = async (dy) => {
    await page.keyboard.down('Control');
    await page.mouse.wheel(0, dy);
    await page.keyboard.up('Control');
  };
  const before = await page.evaluate(() => JSON.parse(JSON.stringify(window.PB.mapView(window.PB.active()))));
  ok('starts at 100%', Math.abs(before.k - 1) < 1e-9, before.k);
  await page.mouse.move(svgBox.x + svgBox.width / 2, svgBox.y + svgBox.height / 2);
  await pinch(-600);
  await page.waitForTimeout(300);
  const zoomed = await page.evaluate(() => JSON.parse(JSON.stringify(window.PB.mapView(window.PB.active()))));
  ok('scrolling up zooms in', zoomed.k > 1.2, zoomed.k);
  const pctTxt = await page.locator('#mapZoomPct').textContent();
  ok('zoom readout updates', /\d+%/.test(pctTxt) && parseInt(pctTxt) > 120, pctTxt);
  const layerT = await page.evaluate(() => document.getElementById('mapZoomLayer').getAttribute('transform'));
  ok('transform applied to the zoom layer', /scale\(/.test(layerT), layerT);
  ok('no NaN in the transform', !/NaN/.test(layerT), layerT);

  // zoom toward the pointer keeps that point roughly fixed
  await page.evaluate(() => { window.PB.resetMapZoom(); });
  await page.mouse.move(svgBox.x + svgBox.width * 0.3, svgBox.y + svgBox.height * 0.3);
  await pinch(-500);
  await page.waitForTimeout(250);
  const z2 = await page.evaluate(() => JSON.parse(JSON.stringify(window.PB.mapView(window.PB.active()))));
  ok('zooming toward the pointer offsets the layer', z2.k > 1 && (z2.x !== 0 || z2.y !== 0),
     JSON.stringify(z2));

  // drag to pan
  const beforePan = await page.evaluate(() => JSON.parse(JSON.stringify(window.PB.mapView(window.PB.active()))));
  await page.mouse.move(svgBox.x + svgBox.width / 2, svgBox.y + svgBox.height / 2);
  await page.mouse.down();
  await page.mouse.move(svgBox.x + svgBox.width / 2 + 90, svgBox.y + svgBox.height / 2 + 40, { steps: 8 });
  await page.mouse.up();
  await page.waitForTimeout(250);
  const afterPan = await page.evaluate(() => JSON.parse(JSON.stringify(window.PB.mapView(window.PB.active()))));
  ok('dragging pans the map', afterPan.x !== beforePan.x || afterPan.y !== beforePan.y,
     JSON.stringify(beforePan) + ' -> ' + JSON.stringify(afterPan));
  ok('panning does not change the zoom level', Math.abs(afterPan.k - beforePan.k) < 1e-9);

  // reset
  await page.click('[data-zm="fit"]');
  await page.waitForTimeout(250);
  const reset = await page.evaluate(() => JSON.parse(JSON.stringify(window.PB.mapView(window.PB.active()))));
  ok('reset returns to 100% and no offset', reset.k === 1 && reset.x === 0 && reset.y === 0, JSON.stringify(reset));
  await page.click('[data-zm="in"]');
  await page.waitForTimeout(200);
  const btnZoom = await page.evaluate(() => window.PB.mapView(window.PB.active()).k);
  ok('the + button zooms in', btnZoom > 1, btnZoom);
  await page.click('[data-zm="fit"]');

  console.log('\n== a binary file is no longer imported as junk ==');
  const countBefore = await page.evaluate(() => window.PB.state.records.length);
  await page.evaluate(() => {
    const junk = new Uint8Array(3000);
    for (let i = 0; i < junk.length; i++) junk[i] = i % 256;
    window.PB.handleFiles([new File([junk], 'mystery.bin', { type: 'application/octet-stream' })]);
  });
  await page.waitForTimeout(700);
  const countAfter = await page.evaluate(() => window.PB.state.records.length);
  ok('a non-sequence binary file is refused, not imported', countAfter === countBefore,
     countBefore + ' -> ' + countAfter);

  console.log('\n== a linear construct opens in the linear projection ==');
  const linDna = buildDna(randDNA(1200, 77), false, '<Features></Features>', '', '');
  await page.evaluate(bytes => {
    window.PB.handleFiles([new File([new Uint8Array(bytes)], 'linearFrag.dna', { type: 'application/octet-stream' })]);
  }, Array.from(linDna));
  await page.waitForTimeout(800);
  await page.click('[data-v="map"]');
  await page.waitForTimeout(600);
  const mode2 = await page.evaluate(() => ({ m: window.PB.state.mapMode, t: window.PB.active().topology }));
  ok('linear construct uses the linear map', mode2.m === 'linear' && mode2.t === 'linear', JSON.stringify(mode2));

  await page.click('[data-v="map"]');
  await page.waitForTimeout(400);
  await page.evaluate(() => {
    const r = window.PB.state.records.find(x => x.name === 'pMyPlasmid');
    if (r) { window.PB.state.activeId = r.id; window.PB.renderAll(); }
  });
  await page.waitForTimeout(700);
  await page.screenshot({ path: path.join(__dirname, 'shot-dna.png') });
  await page.mouse.move(svgBox.x + svgBox.width / 2, svgBox.y + svgBox.height / 2);
  await page.mouse.wheel(0, -700);
  await page.waitForTimeout(400);
  await page.screenshot({ path: path.join(__dirname, 'shot-zoom.png') });

  ok('no console errors throughout', errors.length === 0, errors.slice(0, 4).join(' | '));
  await browser.close();
  console.log('\n==========================================');
  console.log('   ' + pass + ' passed, ' + fail + ' failed');
  console.log('==========================================\n');
  process.exit(fail ? 1 : 0);
})().catch(e => { console.error('HARNESS ERROR', e); process.exit(2); });
