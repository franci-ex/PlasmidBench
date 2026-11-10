/* ===========================================================================
   Talking to eLabFTW.

   Section 7.3 of the project evaluation says one question decides the whole
   architecture: may a browser call eLabFTW's API directly, or does the
   server's CORS configuration forbid it? And it says to answer that before
   designing anything. The answer cannot be read off the code, and it cannot
   be got from a test suite against someone's real instance — so mock-elab.js
   reproduces the behaviours that matter and these tests drive each of them
   through the real panel.

   The cases, and why each one is here:

     cors=full      headers set properly: direct access works, and the whole
                    read-and-add cycle has to work end to end
     cors=simple    Allow-Origin set, Allow-Headers missing. This is a
                    DEFAULT INSTALL: a plain request succeeds and one
                    carrying the API key is blocked. The test must say that
                    in those words rather than "network error", because the
                    person reading it has to know whether to change the
                    server or use the file route
     cors=none      no headers at all: nothing works
     readonly key   reading works, adding is refused, and the refusal says
                    what to do about it and that nothing was changed
     bad key        the server took the request and rejected the key, which
                    is a different answer from being blocked

   And the rule the whole module exists to keep: nothing is ever replaced.
   The suite checks the requests the server actually received, so that is
   verified against the wire rather than against intent.
   ========================================================================= */
const { chromium } = require('playwright');
const path = require('path');
const mock = require('./mock-elab.js');

let pass = 0, fail = 0;
function ok(name, cond, extra) {
  if (cond) { pass++; console.log('  PASS  ' + name); }
  else { fail++; console.log('  FAIL  ' + name + (extra !== undefined ? '   -> ' + extra : '')); }
}

const PAGE_FILE = 'file://' + path.resolve(__dirname, '..', 'plasmidbench.html');

(async () => {
  const browser = await chromium.launch();
  const errors = [];

  /* The page is served over http rather than opened from file://, because a
     file:// page has a null origin and every cross-origin request from it is
     refused whatever the server says — which would make every case below
     look like the blocked one. A lab serving this from a web server is the
     deployment the evaluation describes, so that is what is tested. */
  /* the built page lives at the repo root, one level up from the suites */
  const serveDir = path.resolve(__dirname, '..');
  const http = require('http');
  const fs = require('fs');
  const pageServer = http.createServer((req, res) => {
    const p = req.url === '/' ? '/plasmidbench.html' : req.url.split('?')[0];
    const file = path.join(serveDir, path.basename(p));
    fs.readFile(file, (err, data) => {
      if (err) { res.statusCode = 404; res.end('no'); return; }
      res.setHeader('Content-Type', /\.html$/.test(file) ? 'text/html' : 'text/plain');
      res.end(data);
    });
  });
  await new Promise(r => pageServer.listen(0, '127.0.0.1', r));
  const PAGE = 'http://127.0.0.1:' + pageServer.address().port + '/plasmidbench.html';
  void PAGE_FILE;

  /* ------------------------------------------------------------------ */
  async function openPage() {
    const page = await browser.newPage({ viewport: { width: 1500, height: 1100 } });
    page.on('pageerror', e => errors.push('pageerror: ' + String(e)));
    page.on('console', m => {
      if (m.type() === 'error' &&
          !/Failed to load resource|ERR_(TUNNEL|NAME|INTERNET|CONNECTION)|CORS|Access-Control|net::ERR_FAILED/.test(m.text()))
        errors.push(m.text());
    });
    await page.goto(PAGE);
    await page.waitForTimeout(900);
    return page;
  }

  async function configure(page, url, key) {
    await page.evaluate(([u, k]) => {
      window.PB.elab.setConfig({ url: u, key: k, remember: true, scope: 3 });
    }, [url, key]);
  }

  async function runTest(page) {
    return page.evaluate(() => window.PB.elab.test());
  }

  /* ================================================================== */
  console.log('\n== a server configured for it: the test says so ==');
  let srv = await mock.start({ cors: 'full' });
  let page = await openPage();
  await configure(page, srv.url, srv.keys.rw);
  let r = await runTest(page);
  ok('the verdict is that direct access works', r.verdict === 'direct browser access works', r.verdict);
  ok('the preflight step passed',
     r.steps.some(s => /carrying the API key/.test(s.name) && s.ok === true),
     JSON.stringify(r.steps.map(s => [s.name, s.ok])));
  ok('it names whose key it is', r.user && /Tong/.test(r.user.name), JSON.stringify(r.user));
  ok('it reports the write permission', r.canWrite === true, r.canWrite);
  ok('and the version it found',
     r.steps.some(s => /version/i.test(s.name) && /5\.2/.test(s.note)),
     JSON.stringify(r.steps.filter(s => /version/i.test(s.name))));

  console.log('\n== finding a record and reading its file ==');
  const search = await page.evaluate(() => window.PB.elab.search('pTEST', 'items'));
  ok('the search finds it', search.ok && search.rows.length === 1, JSON.stringify(search).slice(0, 160));
  ok('with its title and category',
     search.ok && search.rows[0].title === 'pTEST plasmid' && search.rows[0].category === 'Plasmid',
     JSON.stringify(search.rows && search.rows[0]));
  const ups = await page.evaluate(() => window.PB.elab.uploads('items', 1));
  ok('its attachment is listed', ups.ok && ups.rows.length === 1 && ups.rows[0].name === 'pTEST.gb',
     JSON.stringify(ups).slice(0, 160));
  ok('with the checksum eLabFTW holds', ups.ok && /^[0-9a-f]{64}$/.test(ups.rows[0].hash),
     ups.ok && ups.rows[0].hash);

  const got = await page.evaluate(async () => {
    const u = await window.PB.elab.uploads('items', 1);
    const f = u.rows[0];
    const d = await window.PB.elab.download('items', 1, f.id, f.hash, f.hashAlgo);
    if (!d.ok) return { ok: false, error: d.error };
    const text = new TextDecoder().decode(d.bytes);
    const recs = window.Bio.parseGenBank(text);
    return { ok: true, verified: d.verified, name: recs[0].name, len: recs[0].seq.length,
             topo: recs[0].topology, feats: recs[0].features.length };
  });
  ok('the file comes down and parses', got.ok && got.name === 'pTEST', JSON.stringify(got));
  ok('its checksum is checked against what eLabFTW stored', got.verified === true, got.verified);
  ok('and it is the plasmid it should be',
     got.len === 60 && got.topo === 'circular' && got.feats === 2, JSON.stringify(got));

  console.log('\n== a corrupted download is refused rather than opened ==');
  const bad = await mock.start({ cors: 'full', corruptDownloads: true });
  const corrupt = await page.evaluate(async (base) => {
    window.PB.elab.setConfig({ url: base, key: 'test-key-rw', remember: true, scope: 3 });
    const u = await window.PB.elab.uploads('items', 1);
    const f = u.rows[0];
    const d = await window.PB.elab.download('items', 1, f.id, f.hash, f.hashAlgo);
    return { ok: d.ok, kind: d.kind, error: d.error || '' };
  }, bad.url);
  ok('a file that does not match its checksum is not opened',
     corrupt.ok === false && corrupt.kind === 'checksum', JSON.stringify(corrupt));
  ok('and the reason says so plainly', /does not match the checksum/.test(corrupt.error), corrupt.error);
  await bad.stop();
  await configure(page, srv.url, srv.keys.rw);

  console.log('\n== adding a result back, without replacing anything ==');
  const sent = await page.evaluate(async () => {
    const P = window.PB, B = window.Bio;
    const rec = P.makeRecord({ name: 'pTEST_SPDYC', seq: 'ATGAAACGCATTAGCACCACCATTACCACC',
                               topology: 'circular', features: [] });
    P.addRecord(rec);
    const sum = B.seqChecksum(rec.seq, { topology: rec.topology });
    const r = await P.elab.attach('items', 1, 'pTEST_SPDYC.gb', B.toGenBank(rec), {
      comment: 'GenBank from PlasmidBench', checksum: sum,
      note: 'Added from PlasmidBench. Not verified stock.'
    });
    return { ok: r.ok, error: r.error, checksum: r.checksum, sum: sum };
  });
  ok('the attachment is accepted', sent.ok, JSON.stringify(sent).slice(0, 180));
  ok('and carries the sequence checksum', sent.checksum === sent.sum, sent.checksum);
  const after = await page.evaluate(() => window.PB.elab.uploads('items', 1));
  ok('the record now has two files, not one replaced',
     after.ok && after.rows.length === 2, after.ok && after.rows.length);
  ok('the original is still there, untouched',
     after.ok && after.rows.some(f => f.name === 'pTEST.gb' &&
       f.hash === srv.state.items[0].uploads[0].hash),
     JSON.stringify(after.rows && after.rows.map(f => f.name)));
  ok('the new one says what made it and when',
     after.ok && after.rows.some(f => /PlasmidBench/.test(f.comment) && /sha256:/.test(f.comment)),
     JSON.stringify(after.rows && after.rows.map(f => f.comment)));
  ok('a comment was added to the record itself',
     srv.state.items[0].comments.length === 1 &&
     /PlasmidBench/.test(srv.state.items[0].comments[0].comment),
     JSON.stringify(srv.state.items[0].comments));

  /* the rule, checked against the wire rather than against intent */
  console.log('\n== nothing in the module can replace, patch or delete ==');
  const methods = srv.seen.map(s => s.method + ' ' + s.path);
  ok('no PATCH was sent', !srv.seen.some(s => s.method === 'PATCH'),
     methods.filter(m => /PATCH/.test(m)).join(', '));
  ok('no DELETE was sent', !srv.seen.some(s => s.method === 'DELETE'),
     methods.filter(m => /DELETE/.test(m)).join(', '));
  ok('nothing was posted to an existing upload (which would replace it)',
     !srv.seen.some(s => s.method === 'POST' && /\/uploads\/\d+$/.test(s.path)),
     methods.filter(m => /POST.*uploads\/\d/.test(m)).join(', '));
  /* Two kinds of request carry no key, both on purpose: the preflight,
     which a browser sends by itself and which cannot carry one, and the
     test's first probe, which is sent keyless so that "unreachable" can be
     told apart from "the preflight is refused". Every request that reads or
     writes anything has to carry it. */
  const dataRequests = srv.seen.filter(s => s.method !== 'OPTIONS' && s.path !== '/api/v2/info');
  ok('every request that reads or writes carried the key',
     dataRequests.length > 4 && dataRequests.every(s => s.auth),
     JSON.stringify(dataRequests.filter(s => !s.auth).slice(0, 2)));
  ok('and the only keyless request is the test’s own first probe',
     srv.seen.filter(s => s.method !== 'OPTIONS' && !s.auth).every(s => s.path === '/api/v2/info'),
     JSON.stringify(srv.seen.filter(s => s.method !== 'OPTIONS' && !s.auth).map(s => s.path)));

  console.log('\n== proposing a construct as a new record ==');
  const prop = await page.evaluate(async () => {
    const P = window.PB;
    const rec = P.active();
    return P.elab.propose(rec, { from: 'pTEST + SPDYC', how: 'In-Fusion, 15 bp arms' });
  });
  ok('a new resource is created', prop.ok && prop.id, JSON.stringify(prop));
  const created = srv.state.items.find(i => i.id === prop.id);
  ok('its title says it is a proposal, not stock',
     created && /^\[proposed\]/.test(created.title), created && created.title);
  ok('its body carries the checksum and where it came from',
     created && /sha256:/.test(created.body) && /pTEST \+ SPDYC/.test(created.body),
     created && created.body.slice(0, 200));
  ok('and it says to sequence it before trusting it',
     created && /[Ss]equence it/.test(created.body), created && created.body.slice(-120));
  ok('with the GenBank attached', created && created.uploads.length === 1,
     created && created.uploads.length);

  await page.close();
  await srv.stop();

  /* ================================================================== */
  console.log('\n== a default install: the preflight is refused ==');
  /* Allow-Origin present, Allow-Headers absent. A plain GET is permitted and
     a request carrying the API key is not. The test has to tell these apart,
     because the first says "it is reachable" and the second says "change
     the server or use files". */
  const simple = await mock.start({ cors: 'simple' });
  page = await openPage();
  await configure(page, simple.url, simple.keys.rw);
  r = await runTest(page);
  ok('the verdict is that CORS blocks it', r.verdict === 'blocked by CORS', r.verdict);
  ok('and it does not say the server is unreachable, because it is not',
     r.steps.some(s => /Something answers/.test(s.name) && s.ok === true),
     JSON.stringify(r.steps.map(s => [s.name, s.ok])));
  const corsStep = r.steps.find(s => /carrying the API key/.test(s.name));
  ok('it names the header that is missing',
     corsStep && /Access-Control-Allow-Headers/.test(corsStep.note), corsStep && corsStep.note);
  ok('and it names the way that does work meanwhile',
     corsStep && /file exchange|proxy/.test(corsStep.note), corsStep && corsStep.note);
  /* and the panel offers the file route, which needs none of this */
  const panel = await page.evaluate(() => {
    window.PB.toolElab();
    return new Promise(res => setTimeout(() => {
      const m = document.querySelector('#modalRoot');
      res({ text: m ? m.textContent : '', open: !!m && !m.hidden });
    }, 400));
  });
  ok('the panel is there', panel.open);
  ok('and describes the file route in full',
     /download the GenBank/i.test(panel.text) && /audit trail is identical/i.test(panel.text),
     panel.text.slice(0, 120));
  ok('it says the key never goes into the distributed file',
     /never in the file you downloaded/i.test(panel.text), '');
  ok('and that results are added rather than replacing anything',
     /never replaced/i.test(panel.text), '');
  await page.close();
  await simple.stop();

  console.log('\n== a server with no CORS headers at all ==');
  const none = await mock.start({ cors: 'none' });
  page = await openPage();
  await configure(page, none.url, none.keys.rw);
  r = await runTest(page);
  ok('nothing answers, and it says that rather than guessing',
     r.verdict === 'nothing answered', r.verdict);
  ok('the note points at the missing header',
     r.steps.some(s => /Access-Control-Allow-Origin/.test(s.note)),
     JSON.stringify(r.steps.map(s => s.note.slice(0, 60))));
  await page.close();
  await none.stop();

  console.log('\n== a read-only key ==');
  const ro = await mock.start({ cors: 'full', readonly: true });
  page = await openPage();
  await configure(page, ro.url, ro.keys.ro);
  r = await runTest(page);
  ok('reading still works', r.verdict === 'direct browser access works', r.verdict);
  ok('and it says the key cannot add anything', r.canWrite === false, r.canWrite);
  const refused = await page.evaluate(async () => {
    const r2 = await window.PB.elab.attach('items', 1, 'x.gb', 'LOCUS x 1 bp\n//\n', {});
    return { ok: r2.ok, kind: r2.kind, error: r2.error };
  });
  ok('adding is refused', refused.ok === false && refused.kind === 'perm', JSON.stringify(refused));
  ok('the refusal says what to do instead', /read-write key|attach it there yourself/.test(refused.error), refused.error);
  ok('and says nothing was changed', /nothing has been changed/i.test(refused.error), refused.error);
  ok('and nothing was', ro.state.items[0].uploads.length === 1, ro.state.items[0].uploads.length);
  await page.close();
  await ro.stop();

  console.log('\n== a key the server rejects ==');
  const good = await mock.start({ cors: 'full' });
  page = await openPage();
  await configure(page, good.url, 'not-a-real-key');
  r = await runTest(page);
  ok('that is reported as the key, not as CORS', r.verdict === 'the key was rejected', r.verdict);
  const step = r.steps.find(s => /carrying the API key/.test(s.name));
  ok('and it says CORS is fine', step && /CORS is fine/.test(step.note), step && step.note);

  console.log('\n== the key is kept where it was asked to be ==');
  const keyHandling = await page.evaluate(async () => {
    const P = window.PB;
    P.elab.setConfig({ url: 'https://eln.example.org', key: 'session-only-key', remember: false });
    /* session only: nothing in the stored settings */
    let stored = JSON.stringify(P.lsGet('elab', {}));
    const inSession = P.elab.config().key;
    P.elab.setConfig({ key: 'remembered-key', remember: true });
    const storedNow = JSON.stringify(P.lsGet('elab', {}));
    P.elab.forget();
    return { afterSessionOnly: stored, usableInSession: inSession,
             afterRemember: storedNow, afterForget: JSON.stringify(P.lsGet('elab', {})),
             keyAfterForget: P.elab.config().key };
  });
  ok('"this session only" keeps the key out of storage',
     !/session-only-key/.test(keyHandling.afterSessionOnly), keyHandling.afterSessionOnly);
  ok('but it still works for the session', keyHandling.usableInSession === 'session-only-key',
     keyHandling.usableInSession);
  ok('"remember it" stores it', /remembered-key/.test(keyHandling.afterRemember), '');
  ok('and forgetting clears both', keyHandling.afterForget === '{}' && !keyHandling.keyAfterForget,
     keyHandling.afterForget + ' / ' + keyHandling.keyAfterForget);

  console.log('\n== an address pasted from the browser bar still works ==');
  const urls = await page.evaluate(() => {
    const out = {};
    for (const u of ['https://eln.lab.org', 'https://eln.lab.org/', 'https://eln.lab.org/api/v2',
                     'https://eln.lab.org/experiments.php?mode=view&id=7',
                     'https://eln.lab.org/database.php']) {
      window.PB.elab.setConfig({ url: u, key: 'k', remember: false });
      out[u] = window.PB.elab.base();
    }
    return out;
  });
  const allSame = Object.keys(urls).every(k => urls[k] === 'https://eln.lab.org/api/v2');
  ok('every form of the address resolves to the API root', allSame, JSON.stringify(urls));

  await page.close();
  await good.stop();
  await new Promise(r2 => pageServer.close(r2));

  console.log('\n== console ==');
  ok('no page errors anywhere in that', errors.length === 0, errors.slice(0, 3).join(' | '));

  console.log('\n   ' + pass + ' passed, ' + fail + ' failed\n');
  await browser.close();
  process.exit(fail ? 1 : 0);
})();
