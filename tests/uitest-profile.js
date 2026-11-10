/* ===========================================================================
   Making the tool this lab's.

   The project's distinctive claim is that it knows what a particular freezer
   holds and what a particular bench runs on. That only means anything if a
   second lab can take the whole setup without editing the source, and if a
   new student can be handed a file rather than a list of settings to type
   in again. So: one profile, exported, imported, and — the part that
   matters — actually changing what the page calculates.

   The two things the profile must never carry are checked here rather than
   left to whoever adds the next field: no API key, and no records. Both are
   enforced in profile.js by pattern, and a test that adds a key-shaped
   setting and tries to export it is the only way to know the enforcement
   works.
   ========================================================================= */
const { chromium } = require('playwright');
const path = require('path');
let pass = 0, fail = 0;
function ok(name, cond, extra) {
  if (cond) { pass++; console.log('  PASS  ' + name); }
  else { fail++; console.log('  FAIL  ' + name + (extra !== undefined ? '   -> ' + extra : '')); }
}

const PAGE = 'file://' + path.resolve(__dirname, '..', 'plasmidbench.html');

(async () => {
  const browser = await chromium.launch();
  const page = await browser.newPage({ viewport: { width: 1500, height: 1100 } });
  const errors = [];
  page.on('pageerror', e => errors.push('pageerror: ' + String(e)));
  page.on('console', m => {
    if (m.type() === 'error' && !/Failed to load resource|ERR_(TUNNEL|NAME|INTERNET|CONNECTION)/.test(m.text()))
      errors.push(m.text());
  });
  await page.goto(PAGE);
  await page.waitForTimeout(1100);

  console.log('\n== a profile holds configuration, and nothing else ==');
  const exported = await page.evaluate(() => {
    const P = window.PB;
    /* a lab's worth of settings, plus two things that must not travel */
    P.lsSet('lab.name', 'Tong lab');
    P.lsSet('lab.pi', 'Y. Tong');
    P.lsSet('tmCond', { oligo: 0.5, Na: 100, Mg: 2.0, dNTP: 0.3 });
    P.lsSet('suppliers', ['N', 'T']);
    P.lsSet('sanger.good', 520);
    P.lsSet('sanger.dead', 45);
    P.lsSet('sanger.facility', 'the core facility downstairs');
    P.elab.setConfig({ url: 'https://eln.tonglab.org', key: 'SUPER-SECRET-KEY', remember: true });
    P.lsSet('stock', [{ name: 'OL-1', seq: 'ACGT' }]);     /* records: must not travel */
    P.lsSet('primers', [{ name: 'p1', seq: 'ACGT' }]);
    const text = P.profile.toText();
    return { text: text, parsed: JSON.parse(text), filename: P.profile.filename() };
  });
  ok('it says what format it is', exported.parsed.format === 'plasmidbench-lab-profile',
     exported.parsed.format);
  ok('and is named after the lab', exported.parsed.name === 'Tong lab', exported.parsed.name);
  ok('the file is named after the lab too', exported.filename === 'tong_lab-plasmidbench-profile.json',
     exported.filename);
  ok('it carries the bench conditions',
     exported.parsed.values.tmCond && exported.parsed.values.tmCond.Na === 100,
     JSON.stringify(exported.parsed.values.tmCond));
  ok('the suppliers', JSON.stringify(exported.parsed.values.suppliers) === '["N","T"]',
     JSON.stringify(exported.parsed.values.suppliers));
  ok('what the sequencing facility really gives',
     exported.parsed.values['sanger.good'] === 520 && exported.parsed.values['sanger.dead'] === 45,
     JSON.stringify(exported.parsed.values));
  ok('and the eLabFTW address', exported.parsed.values['elab.url'] === 'https://eln.tonglab.org',
     exported.parsed.values['elab.url']);

  console.log('\n== what it must not carry, checked against the file itself ==');
  ok('the API key is nowhere in it', exported.text.indexOf('SUPER-SECRET-KEY') === -1,
     exported.text.length + ' characters, searched');
  ok('nor is the word key used as a field', !/"[^"]*key[^"]*"\s*:/i.test(exported.text), '');
  ok('no oligo inventory', !/OL-1/.test(exported.text), '');
  ok('no primer library', !/"p1"/.test(exported.text), '');
  ok('and it says so on its face', /No API keys, no plasmids/.test(exported.parsed.note),
     exported.parsed.note);

  console.log('\n== importing shows what would change before changing it ==');
  const incoming = {
    format: 'plasmidbench-lab-profile', version: 1, name: 'Francis lab', exported: '2026-11-01',
    values: {
      'lab.name': 'Francis lab',
      tmCond: { oligo: 0.25, Na: 50, Mg: 1.5, dNTP: 0.2 },
      'sanger.good': 850,
      suppliers: ['N'],
      'bench.kit': 'In-Fusion HD',
      /* a key and a record list, as a badly made profile might carry */
      'elab.key': 'ANOTHER-SECRET',
      stock: [{ name: 'X', seq: 'ACGT' }],
      /* and a setting from a future version */
      'future.thing': 42
    }
  };
  const rep = await page.evaluate(j => window.PB.profile.inspect(j), incoming);
  ok('it reads', rep.ok, rep.error);
  ok('named after the lab it came from', rep.name === 'Francis lab', rep.name);
  ok('it lists the settings that would change',
     rep.changes.length >= 4, JSON.stringify(rep.changes.map(c => c.key)));
  ok('each with what it is for, in words, not a key name',
     rep.changes.every(c => c.what && c.what !== c.key),
     JSON.stringify(rep.changes.map(c => c.what).slice(0, 2)));
  ok('and with the value now and the value after',
     rep.changes.some(c => c.key === 'sanger.good' && c.from === 520 && c.to === 850),
     JSON.stringify(rep.changes.find(c => c.key === 'sanger.good')));
  ok('a key in the incoming file is refused, not applied',
     rep.refused.indexOf('elab.key') >= 0, JSON.stringify(rep.refused));
  ok('so are records', rep.refused.indexOf('stock') >= 0, JSON.stringify(rep.refused));
  ok('a setting from a later version is skipped and counted',
     rep.ignored.indexOf('future.thing') >= 0, JSON.stringify(rep.ignored));
  ok('and nothing has changed yet',
     await page.evaluate(() => window.PB.lsGet('sanger.good', null)) === 520, '');

  console.log('\n== and then it changes what the page calculates ==');
  /* The arriving profile has to differ from what is running, or there is
     nothing to see. Salt five times higher moves a melting temperature by
     several degrees — a real difference, and the reason the import step
     shows what it will change before changing it. */
  incoming.values.tmCond = { oligo: 1.0, Na: 500, Mg: 4.0, dNTP: 0.4 };
  const applied = await page.evaluate(async j => {
    const P = window.PB, B = window.Bio;
    /* One oligo's melting temperature, before and after, read through
       P.conditions() — which is what the whole page calculates with, so
       this measures the path a person's results actually travel. */
    const primer = 'TTGTATTTCCAGGGCTCACCC';
    const condBefore = P.conditions();
    const before = B.tmNN(primer, condBefore);
    const rep2 = P.profile.inspect(j);
    const n = P.profile.apply(rep2);
    const after = B.tmNN(primer, P.conditions());
    return {
      n: n, before: before, after: after,
      cond: P.conditions(), condBefore: condBefore,
      labName: P.lsGet('lab.name', ''),
      sangerGood: P.lsGet('sanger.good', null),
      /* the key the refused field tried to overwrite is untouched */
      key: P.elab.config().key,
      url: P.elab.config().url
    };
  }, incoming);
  ok('the settings are taken', applied.n >= 4, applied.n);
  ok('the lab name follows', applied.labName === 'Francis lab', applied.labName);
  /* conditions() hands the engine molar monovalent salt and molar oligo but
     millimolar magnesium and dNTP; that mixture is the documented contract
     at engine.js:174, so the expected values are converted the same way
     here rather than assumed to be the stored ones. */
  ok('the bench conditions follow',
     applied.cond.Na === 0.5 && applied.cond.oligo === 1.0e-6 && applied.cond.Mg === 4.0,
     JSON.stringify(applied.cond) + ' was ' + JSON.stringify(applied.condBefore));
  ok('and a melting temperature actually changes with them',
     Math.abs(applied.before - applied.after) > 1,
     applied.before.toFixed(2) + ' → ' + applied.after.toFixed(2));
  ok('the sequencing read length follows', applied.sangerGood === 850, applied.sangerGood);
  ok('the API key was not touched by the import', applied.key === 'SUPER-SECRET-KEY', applied.key);
  ok('and neither was the address, which the profile did not set',
     applied.url === 'https://eln.tonglab.org', applied.url);

  console.log('\n== the suppliers setting reaches the enzyme list ==');
  const enz = await page.evaluate(() => {
    const P = window.PB, B = window.Bio;
    const count = () => {
      P.state.view = 'map';
      P.renderEnzymeRail();
      return document.querySelectorAll('#enzList .enz-row').length;
    };
    P.state.enzSet = 'all'; P.lsSet('enzSet', 'all');
    P.lsSet('suppliers', []);
    const all = count();
    /* one small supplier only */
    P.lsSet('suppliers', ['Q']);
    const few = count();
    const name = (B.SUPPLIER_NAMES || {})['Q'] || 'Q';
    P.lsSet('suppliers', []);
    const back = count();
    return { all: all, few: few, back: back, name: name };
  });
  ok('with no suppliers set, every enzyme is offered', enz.all > 400, enz.all);
  ok('choosing one supplier cuts the list to what they sell',
     enz.few > 0 && enz.few < enz.all, enz.few + ' of ' + enz.all + ' (' + enz.name + ')');
  ok('and clearing it puts them all back', enz.back === enz.all, enz.back);

  console.log('\n== the sequencing panel opens on the facility’s real numbers ==');
  const sq = await page.evaluate(async () => {
    const P = window.PB;
    P.lsSet('sanger.good', 480); P.lsSet('sanger.dead', 60);
    const rec = P.makeRecord({ name: 'x', seq: 'ACGT'.repeat(500), topology: 'circular', features: [] });
    P.addRecord(rec);
    P.state.selection = { start: 100, end: 300 };
    P.toolSanger();
    await new Promise(r => setTimeout(r, 400));
    const g = document.querySelector('#sqGood'), d = document.querySelector('#sqDead');
    const out = { good: g && g.value, dead: d && d.value };
    P.closeModal();
    return out;
  });
  ok('the confident read length is the lab’s, not the textbook default',
     sq.good === '480', sq.good);
  ok('and so is the dead zone', sq.dead === '60', sq.dead);

  console.log('\n== a file that is not a profile ==');
  const junk = await page.evaluate(() => [
    window.PB.profile.inspect('not json at all'),
    window.PB.profile.inspect('{"hello":"world"}')
  ]);
  ok('plain text is refused, saying why', !junk[0].ok && /not even JSON/.test(junk[0].error),
     junk[0].error);
  ok('and so is unrelated JSON', !junk[1].ok && /not a PlasmidBench lab profile/.test(junk[1].error),
     junk[1].error);

  console.log('\n== console ==');
  ok('no page errors anywhere in that', errors.length === 0, errors.slice(0, 3).join(' | '));

  console.log('\n   ' + pass + ' passed, ' + fail + ' failed\n');
  await browser.close();
  process.exit(fail ? 1 : 0);
})();
