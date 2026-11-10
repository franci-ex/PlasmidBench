/* ===========================================================================
   What happens when the viewer's host will not write the file.

   A published page cannot write a file itself. It asks whatever is showing
   it, and that host keeps a list of file types it will accept — a list that
   differs between a browser tab, the desktop app and a phone. SVG is on the
   edge of it. The old behaviour when a save was refused was to put the text
   on the clipboard and announce that cheerfully, which for a 72,000-character
   SVG is a dead end dressed up as a success.

   These tests stand in a fake host and make it refuse in each of the ways it
   really can, then check that what the person is told is true and that the
   ways out actually lead somewhere.
   ========================================================================= */
const { chromium } = require('playwright');
const path = require('path'), fs = require('fs');
let pass = 0, fail = 0;
function ok(name, cond, extra) {
  if (cond) { pass++; console.log('  PASS  ' + name); }
  else { fail++; console.log('  FAIL  ' + name + (extra !== undefined ? '   -> ' + extra : '')); }
}

/* a host whose allowlist and failure mode the test chooses */
const FAKE = `(() => {
  window.__saves = [];
  window.__mode = { fail: null, allow: ['txt','json','png'] };
  window.claude = {
    use: async (name) => {
      if (name !== 'downloads') return null;
      return {
        save: async ({ filename, data }) => {
          const ext = (filename.split('.').pop() || '').toLowerCase();
          window.__saves.push({ filename: filename, ext: ext,
            size: typeof data === 'string' ? data.length : (data.size || 0),
            text: typeof data === 'string' ? data : null });
          if (window.__mode.fail) { const e = new Error('nope'); e.code = window.__mode.fail; throw e; }
          if (window.__mode.allow.indexOf(ext) === -1) { const e = new Error('no'); e.code = 'rejected_extension'; throw e; }
          return { status: 'saved' };
        }
      };
    }
  };
})();`;

(async () => {
  const browser = await chromium.launch();
  const page = await browser.newPage({ viewport: { width: 1400, height: 1000 } });
  const errors = [];
  page.on('pageerror', e => errors.push('pageerror: ' + String(e)));
  page.on('console', m => { if (m.type() === 'error' && !/Failed to load resource/.test(m.text())) errors.push('console: ' + m.text()); });
  await page.addInitScript(FAKE);
  await page.goto('file://' + path.resolve('plasmidbench.html'));
  await page.waitForTimeout(900);
  await page.evaluate(t => window.PB.importText(t, 'pGro7.gb'),
    fs.readFileSync(path.resolve(__dirname, 'fixtures/pGro7.gb'), 'utf8'));
  await page.waitForTimeout(500);

  const ready = await page.evaluate(() => !!window.PB.state.downloads);
  ok('the page finds the host save capability', ready);

  const dialog = () => page.evaluate(() => {
    const m = document.querySelector('#modalRoot');
    if (!m || m.hidden) return null;
    return { text: m.textContent, buttons: Array.from(m.querySelectorAll('button, a.btn')).map(b => b.textContent.trim()) };
  });
  const toastText = () => page.evaluate(() => (document.querySelector('#toast') || {}).textContent || '');
  const reset = () => page.evaluate(() => { window.__saves = []; window.__mode.fail = null; window.PB.closeModal(); });

  /* ------------------------------------------------- the save that works */
  console.log('\n== a file the host accepts ==');
  await reset();
  await page.evaluate(() => window.PB.offerFile('notes.txt', 'hello', 'Notes'));
  await page.waitForTimeout(250);
  ok('goes straight through', await page.evaluate(() => window.__saves.length === 1 && window.__saves[0].filename === 'notes.txt'));
  ok('and says it was saved', /Saved notes\.txt/.test(await toastText()), await toastText());
  ok('with no dialog in the way', (await dialog()) === null);

  /* --------------------------------------- the type the host will not take */
  console.log('\n== an SVG the host will not write ==');
  await reset();
  await page.evaluate(async () => {
    const svg = window.PB.mapSVGString({ whole: true, theme: 'light', background: 'solid' });
    window.__svgLen = svg.length;
    await window.PB.offerFile('map.svg', svg, 'Map (SVG)', {
      alt: { label: 'Save it as a PNG instead', run: () => { window.__altRan = true; } }
    });
  });
  await page.waitForTimeout(350);
  let d = await dialog();
  ok('the clipboard is not quietly called a save', !/on the clipboard/.test(await toastText()), await toastText());
  ok('a dialog explains what was refused', !!d && /will not save/.test(d.text), d && d.text.slice(0, 120));
  ok('and names the extension, not just "the file"', !!d && /\.svg/.test(d.text), d && d.text.slice(0, 120));
  ok('it offers the format that does work', !!d && d.buttons.some(b => /PNG/.test(b)), d && d.buttons.join(' | '));
  ok('it offers the same bytes under a name the host takes',
    !!d && d.buttons.some(b => /map\.svg\.txt/.test(b)), d && d.buttons.join(' | '));
  ok('it still offers the clipboard, as a choice rather than a consolation',
    !!d && d.buttons.some(b => /Copy the text/.test(b)), d && d.buttons.join(' | '));
  ok('and the text itself is on screen to take', await page.evaluate(() =>
    !!document.querySelector('#sbText') && document.querySelector('#sbText').value.length > 400));
  ok('nothing was substituted without being asked for', await page.evaluate(() =>
    window.__saves.length === 1 && window.__saves[0].ext === 'svg' && !window.__altRan));

  console.log('\n== and the ways out lead somewhere ==');
  await page.evaluate(() => {
    const b = Array.from(document.querySelectorAll('#modalRoot button')).find(x => /map\.svg\.txt/.test(x.textContent));
    b.click();
  });
  await page.waitForTimeout(350);
  const saved = await page.evaluate(() => window.__saves[window.__saves.length - 1]);
  ok('the .txt save really happens', saved && saved.filename === 'map.svg.txt', JSON.stringify(saved && saved.filename));
  ok('and carries the whole SVG, byte for byte, not a summary of it',
    await page.evaluate(() => {
      const s = window.__saves[window.__saves.length - 1];
      return s.text && s.text.length === window.__svgLen && s.text.indexOf('<svg') !== -1;
    }), JSON.stringify({ got: saved && saved.size, want: await page.evaluate(() => window.__svgLen) }));
  ok('and the message says to rename it', /rename it to \.svg/.test(await toastText()), await toastText());

  await reset();
  await page.evaluate(async () => {
    await window.PB.offerFile('map.svg', '<svg/>', 'Map (SVG)',
      { alt: { label: 'Save it as a PNG instead', run: () => { window.__altRan = true; } } });
  });
  await page.waitForTimeout(300);
  await page.evaluate(() => {
    const b = Array.from(document.querySelectorAll('#modalRoot button')).find(x => /PNG/.test(x.textContent));
    b.click();
  });
  await page.waitForTimeout(300);
  ok('the PNG offer runs the PNG export', await page.evaluate(() => window.__altRan === true));

  /* ----------------------------------------------- the other refusals */
  console.log('\n== the other ways a host says no ==');
  const cases = [
    ['declined', /cancelled/i, false],
    ['too_large', /too large/i, false],
    ['rate_limited', /already waiting/i, false],
    ['extension_not_enabled', /switched off/i, true],
    ['unavailable', /could not save/i, true]
  ];
  for (const [code, re, wantsDialog] of cases) {
    await reset();
    await page.evaluate(code => { window.__mode.fail = code; }, code);
    await page.evaluate(() => window.PB.offerFile('map.svg', '<svg/>', 'Map (SVG)'));
    await page.waitForTimeout(300);
    const dd = await dialog(), tt = await toastText();
    ok(code + ' is reported in its own words', re.test(wantsDialog ? (dd ? dd.text : '') : tt),
      wantsDialog ? (dd ? dd.text.slice(0, 110) : 'no dialog') : tt);
    ok(code + (wantsDialog ? ' offers a way on' : ' stops there, with no dialog'),
      wantsDialog ? !!dd : dd === null, dd ? dd.buttons.join(' | ') : 'none');
  }
  await page.evaluate(() => { window.__mode.fail = null; });

  /* ------------------------------------------------- no capability at all */
  console.log('\n== a view with no save capability at all ==');
  await reset();
  await page.evaluate(() => { window.PB.state.downloads = null; });
  await page.evaluate(() => window.PB.offerFile('map.svg', '<svg/>', 'Map (SVG)'));
  await page.waitForTimeout(300);
  d = await dialog();
  /* The wording changed when the real download link was added: the dialog
     now says which file types this view allows rather than that saving is
     "not available", because saving IS available — through the browser's own
     download, which is the first button offered. The clipboard is the last
     resort, not the only one. */
  ok('says so plainly', !!d && /will not save|not available in this view/.test(d.text),
    d && d.text.slice(0, 140));
  ok('and offers the browser download first of all',
    !!d && /^Download /.test((d.buttons.filter(Boolean)[0] || '')), d && d.buttons.join(' | '));
  ok('with a way to open it in a tab as well',
    !!d && d.buttons.some(b => /new tab/.test(b)), d && d.buttons.join(' | '));
  ok('but still hands over the text', !!d && d.buttons.some(b => /Copy the text/.test(b)), d && d.buttons.join(' | '));

  /* ---------------------------------------------- through the real dialog */
  console.log('\n== through the export dialog, as a person would ==');
  await page.evaluate(() => { window.PB.closeModal(); });
  await page.reload();
  await page.waitForTimeout(900);
  await page.evaluate(t => window.PB.importText(t, 'pGro7.gb'),
    fs.readFileSync(path.resolve(__dirname, 'fixtures/pGro7.gb'), 'utf8'));
  await page.waitForTimeout(400);
  await page.evaluate(() => { window.PB.state.view = 'map'; window.PB.renderAll(); window.PB.exportMapDialog(); });
  await page.waitForTimeout(500);
  await page.selectOption('#meFmt', 'svg');
  await page.waitForTimeout(400);
  await page.evaluate(() => {
    const b = Array.from(document.querySelectorAll('#modalRoot button')).find(x => x.textContent.trim() === 'Save image');
    b.click();
  });
  await page.waitForTimeout(700);
  d = await dialog();
  ok('choosing SVG and pressing Save gets the explanation, not a clipboard notice',
    !!d && /will not save/.test(d.text), d ? d.text.slice(0, 120) : 'no dialog');
  /* The real download comes first now; the PNG re-render is the next way on
     for a view that will not take an .svg at all. */
  ok('and the PNG way out is offered',
    !!d && d.buttons.some(b => /PNG/.test(b)), d && d.buttons.join(' | '));
  ok('after the browser download itself',
    !!d && /^Download /.test((d.buttons.filter(Boolean)[0] || '')), d && d.buttons.join(' | '));
  await page.evaluate(() => {
    const b = Array.from(document.querySelectorAll('#modalRoot button')).find(x => /PNG/.test(x.textContent));
    b.click();
  });
  await page.waitForTimeout(1800);
  const png = await page.evaluate(() => window.__saves.filter(s => s.ext === 'png').pop());
  ok('which really writes the PNG of the same map', !!png && png.size > 10000,
    JSON.stringify(png && { f: png.filename, n: png.size }));
  ok('nothing threw anywhere in the run', errors.length === 0, errors.slice(0, 3).join(' | '));

  console.log('\n   ' + pass + ' passed, ' + fail + ' failed\n');
  await browser.close();
  process.exit(fail ? 1 : 0);
})().catch(e => { console.error('HARNESS ERROR', e); process.exit(2); });
