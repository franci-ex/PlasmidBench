// Map image export: the file must stand on its own — no page classes, no
// theme variables, fixed size, the same text as the map — and the PNG must
// actually contain the drawing.
const { chromium } = require('playwright');
const path = require('path'), fs = require('fs');
let pass = 0, fail = 0;
function ok(name, cond, extra) {
  if (cond) { pass++; console.log('  PASS  ' + name); }
  else { fail++; console.log('  FAIL  ' + name + (extra !== undefined ? '   -> ' + extra : '')); }
}
const V = path.resolve(__dirname, 'fixtures/pRSFDuet-1.dna');

(async () => {
  const browser = await chromium.launch();
  const page = await browser.newPage({ viewport: { width: 1400, height: 1000 } });
  const errors = [];
  page.on('pageerror', e => errors.push(String(e)));
  page.on('console', m => {
    if (m.type() === 'error' && !/Failed to load resource|ERR_(TUNNEL|NAME|INTERNET|CONNECTION)/.test(m.text())) errors.push(m.text());
  });
  await page.goto('file://' + path.resolve('plasmidbench.html'));
  await page.waitForTimeout(600);
  const useReal = fs.existsSync(V);
  if (useReal) {
    const b64 = fs.readFileSync(V).toString('base64');
    await page.evaluate(b64 => {
      const bin = atob(b64), buf = new Uint8Array(bin.length);
      for (let i = 0; i < bin.length; i++) buf[i] = bin.charCodeAt(i);
      const rec = window.Bio.parseSnapGene(buf.buffer);
      window.PB.addRecord(Object.assign({ name: 'pRSFDuet-1' }, rec, { topology: 'circular' }));
    }, b64);
  }
  await page.waitForTimeout(400);
  console.log('\n== ' + (useReal ? 'pRSFDuet-1' : 'the demo construct') + ' ==');

  for (const mode of ['circular', 'linear']) {
    console.log('\n== ' + mode + ' map as SVG ==');
    const r = await page.evaluate(mode => {
      const P = window.PB;
      P.state.mapMode = mode; P.state.view = 'map'; P.renderTabs(); P.renderViews();
      const onScreen = Array.from(document.querySelectorAll('#mapWrap svg text')).map(t => t.textContent).filter(Boolean);
      const svg = P.mapSVGString({ mode: mode, whole: true, theme: 'light' });
      const doc = new DOMParser().parseFromString(svg, 'image/svg+xml');
      const err = doc.querySelector('parsererror');
      const root = doc.documentElement;
      const texts = Array.from(doc.querySelectorAll('text')).map(t => t.textContent).filter(Boolean);
      return {
        len: svg.length, parses: !err, rootName: root.nodeName,
        w: root.getAttribute('width'), h: root.getAttribute('height'), ns: root.getAttribute('xmlns'),
        classes: (svg.match(/\sclass="/g) || []).length,
        vars: (svg.match(/var\(--/g) || []).length,
        dataAttrs: (svg.match(/\sdata-[a-z-]+="/g) || []).length,
        texts: texts.length, screen: onScreen.length,
        /* The export is of the WHOLE construct, so it can legitimately carry
           labels the view on screen has left out at its current zoom. What
           must never happen is a label on screen that is missing from the
           file. */
        /* "+11 more sites — zoom in" is an instruction to the reader of the
           screen, not a label on the map: the exported file draws all the
           sites, so the hint has nothing to say in it and is rightly left
           out. */
        sameText: onScreen.filter(t => !/zoom in/.test(t)).every(t => texts.indexOf(t) >= 0),
        missing: onScreen.filter(t => !/zoom in/.test(t) && texts.indexOf(t) < 0),
        title: (doc.querySelector('title') || {}).textContent || '',
        bg: (doc.querySelector('svg > rect') || { getAttribute: () => '' }).getAttribute('fill'),
        /* A feature arc is a stroked path, not a filled one — the colour is
           on the stroke and the fill is none. Look at both, so the check is
           about the colour being a real resolved rgb() and not a theme
           variable, whichever property carries it. */
        featFill: (() => {
          const re = /(?:^|;)\s*(?:fill|stroke):\s*(rgb\([^)]*\))/;
          for (const x of Array.from(doc.querySelectorAll('path'))) {
            const m = re.exec(x.getAttribute('style') || '');
            if (m) return m[1];
            const f = x.getAttribute('fill') || '', st = x.getAttribute('stroke') || '';
            if (/^rgb\(/.test(f)) return f;
            if (/^rgb\(/.test(st)) return st;
          }
          return '';
        })()
      };
    }, mode);
    ok(mode + ': the file parses as XML', r.parses && r.rootName === 'svg');
    ok(mode + ': it declares the SVG namespace', r.ns === 'http://www.w3.org/2000/svg');
    ok(mode + ': it has a fixed width and height', +r.w > 100 && +r.h > 50, r.w + ' x ' + r.h);
    ok(mode + ': no page classes are left', r.classes === 0, r.classes);
    ok(mode + ': no theme variables are left', r.vars === 0, r.vars);
    ok(mode + ': no page data attributes are left', r.dataAttrs === 0, r.dataAttrs);
    ok(mode + ': every label on the map is in the file', r.sameText && r.texts > 3,
      r.texts + ' in file vs ' + r.screen + ' on screen; missing: ' + JSON.stringify(r.missing));
    ok(mode + ': it is titled with the construct', /bp/.test(r.title), r.title);
    ok(mode + ': the background is white for the light export', /^#fff/i.test(r.bg), r.bg);
    ok(mode + ': feature colours are real colours', /^rgb/.test(r.featFill), r.featFill);
  }

  console.log('\n== the page itself is left as it was ==');
  const theme = await page.evaluate(() => {
    const root = document.documentElement;
    root.setAttribute('data-theme', 'dark');
    const P = window.PB; P.state.selection = { start: 10, end: 50 };
    P.mapSVGString({ theme: 'light' });
    return { theme: root.getAttribute('data-theme'), sel: JSON.stringify(P.state.selection) };
  });
  ok('the dark theme is back after a light export', theme.theme === 'dark');
  ok('the selection is back after an export that left it out', theme.sel === '{"start":10,"end":50}');
  const dark = await page.evaluate(() => {
    const svg = window.PB.mapSVGString({ theme: 'dark' });
    return (new DOMParser().parseFromString(svg, 'image/svg+xml').querySelector('svg > rect') || { getAttribute: () => '' }).getAttribute('fill');
  });
  ok('a dark export has a dark background', dark && !/^#fff/i.test(dark), dark);
  const transp = await page.evaluate(() => {
    const svg = window.PB.mapSVGString({ theme: 'light', background: 'transparent' });
    const d = new DOMParser().parseFromString(svg, 'image/svg+xml');
    const first = d.documentElement.querySelector('rect');
    return { firstRectIsBg: !!first && first.getAttribute('width') === d.documentElement.getAttribute('width') &&
             first.getAttribute('x') === null };
  });
  ok('a transparent export has no background rectangle', !transp.firstRectIsBg);
  const withSel = await page.evaluate(() => {
    const a = window.PB.mapSVGString({ theme: 'light', withSelection: false, mode: 'circular' });
    const b = window.PB.mapSVGString({ theme: 'light', withSelection: true, mode: 'circular' });
    return { a: a.length, b: b.length };
  });
  ok('the selection is drawn only when asked for', withSel.b > withSel.a, JSON.stringify(withSel));
  await page.evaluate(() => { document.documentElement.removeAttribute('data-theme'); window.PB.state.selection = null; });

  console.log('\n== PNG ==');
  const png = await page.evaluate(async () => {
    const svg = window.PB.mapSVGString({ mode: 'circular', theme: 'light' });
    const blob = await window.PB.mapPNGBlob(svg, 2);
    const head = new Uint8Array(await blob.slice(0, 8).arrayBuffer());
    const bmp = await createImageBitmap(blob);
    const c = document.createElement('canvas'); c.width = bmp.width; c.height = bmp.height;
    const ctx = c.getContext('2d'); ctx.drawImage(bmp, 0, 0);
    const px = ctx.getImageData(0, 0, c.width, c.height).data;
    let dark = 0, coloured = 0;
    for (let i = 0; i < px.length; i += 4 * 7) {
      const r = px[i], g = px[i + 1], b = px[i + 2];
      if (r + g + b < 380) dark++;
      if (Math.max(r, g, b) - Math.min(r, g, b) > 60) coloured++;
    }
    const m = svg.match(/<svg[^>]*\swidth="([\d.]+)"/);
    return { sig: Array.from(head).join(','), w: bmp.width, h: bmp.height, svgW: m ? +m[1] : 0, dark, coloured, size: blob.size };
  });
  ok('the file is a real PNG', png.sig === '137,80,78,71,13,10,26,10', png.sig);
  ok('it is drawn at the chosen resolution', png.w === Math.round(png.svgW * 2), png.w + ' vs ' + png.svgW * 2);
  ok('it contains dark text and lines, not a blank image', png.dark > 200, png.dark);
  ok('it contains the coloured features', png.coloured > 200, png.coloured);

  console.log('\n== the dialog ==');
  await page.evaluate(() => { const P = window.PB; P.state.view = 'map'; P.renderViews(); });
  await page.evaluate(() => document.querySelector('[data-act="mapexport"]').click());
  await page.waitForTimeout(600);
  const dlg = await page.evaluate(() => {
    const img = document.querySelector('#mePrev img');
    return { prev: !!img, loaded: img ? img.complete && img.naturalWidth > 0 : false,
             scaleHidden: document.querySelector('#meScaleF').style.display === 'none' };
  });
  ok('the dialog shows a preview of the export', dlg.prev && dlg.loaded, JSON.stringify(dlg));
  ok('the PNG resolution choice is hidden while SVG is chosen', dlg.scaleHidden);
  await page.selectOption('#meFmt', 'png'); await page.waitForTimeout(300);
  ok('and shown for PNG', await page.evaluate(() => document.querySelector('#meScaleF').style.display !== 'none'));
  /* There is no host here, which is also the case on a lab machine: the
     page then uses the browser's own download rather than putting the
     picture on screen to be saved by hand. The old fallback is still there
     for an embedded view, where a click can be dropped without an error,
     but it is no longer what a plain browser gets. */
  const [dl] = await Promise.all([
    page.waitForEvent('download', { timeout: 9000 }).catch(() => null),
    page.evaluate(() => {
      const b = Array.from(document.querySelectorAll('#modalRoot button')).find(x => x.textContent.trim() === 'Save image');
      b.click();
    })
  ]);
  await page.waitForTimeout(600);
  ok('with no host, saving the PNG downloads it', dl !== null,
     dl ? dl.suggestedFilename() : 'no download event');
  ok('and it is named as a PNG', dl && /\.png$/i.test(dl.suggestedFilename()),
     dl && dl.suggestedFilename());

  ok('no errors anywhere in the run', errors.length === 0, errors.slice(0, 3).join(' | '));
  console.log('\n   ' + pass + ' passed, ' + fail + ' failed\n');
  await browser.close();
  process.exit(fail ? 1 : 0);
})();
