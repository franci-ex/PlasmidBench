/* ============================================================================
   A small stand-in for eLabFTW's REST API v2, for testing against.

   It exists because the one question section 7.3 says to answer first —
   whether a browser may call eLabFTW directly, which turns on the server's
   CORS configuration — cannot be answered by reading code, and cannot be
   answered against a real instance from a test suite. So the behaviours
   that matter are reproduced here and the suite drives each of them:

     cors=full     Access-Control-Allow-Origin and -Allow-Headers set, so
                   the preflight succeeds and direct access works
     cors=simple   Allow-Origin set but Allow-Headers missing: a plain GET
                   is permitted and a request carrying the API key is not.
                   This is a default install, and the case the connection
                   test has to name correctly.
     cors=none     no CORS headers at all: nothing works from a page
     readonly      the key may read but not add
     down          nothing answers

   The shapes it returns follow eLabFTW's published OpenAPI description:
   the key in an Authorization header with no prefix, 201 plus a location
   header on create, ?format=binary to fetch an attachment's bytes.

   It is a test fixture and nothing else: it stores in memory, it is not
   authentication, and it never runs as part of the page.

   Usage:  node mock-elab.js [port] [mode]
   or:     const m = require('./mock-elab.js'); const s = await m.start({...})
   ============================================================================ */
const http = require('http');
const crypto = require('crypto');

const GOOD_KEY = 'test-key-rw';
const RO_KEY = 'test-key-ro';

function makeState() {
  const gb = [
    'LOCUS       pTEST                     60 bp    ds-DNA     circular SYN 10-NOV-2026',
    'DEFINITION  a tiny plasmid for the tests.',
    'FEATURES             Location/Qualifiers',
    '     source          1..60',
    '     CDS             1..30',
    '                     /label="tinyORF"',
    'ORIGIN',
    '        1 atgaaacgca ttagcaccac cattaccacc accatcacca ttaccacagg taacggtgcg',
    '//',
    ''
  ].join('\n');
  const uploads = [{
    id: 11, real_name: 'pTEST.gb', filesize: gb.length, comment: 'the plasmid map',
    hash: crypto.createHash('sha256').update(gb).digest('hex'),
    hash_algorithm: 'sha256', created_at: '2026-11-01 09:00:00',
    _bytes: Buffer.from(gb, 'utf8')
  }];
  return {
    items: [
      { id: 1, title: 'pTEST plasmid', elabid: '20261101-aaa', category_title: 'Plasmid',
        date: '2026-11-01', fullname: 'Y. Tong', tags: 'plasmid|test', body: '', uploads: uploads,
        comments: [] },
      { id: 2, title: 'pGro7 chaperone vector', elabid: '20261101-bbb', category_title: 'Plasmid',
        date: '2026-11-02', fullname: 'J. Francis', tags: 'plasmid', body: '', uploads: [],
        comments: [] }
    ],
    experiments: [
      { id: 7, title: 'SPDYC cloning, attempt 3', elabid: '20261105-ccc', date: '2026-11-05',
        fullname: 'J. Francis', tags: 'cloning', body: '', uploads: [], comments: [] }
    ],
    nextItemId: 3,
    nextUploadId: 12
  };
}

function start(opts) {
  opts = opts || {};
  const mode = opts.cors || 'full';
  const readonly = !!opts.readonly;
  const state = makeState();
  const seen = [];                       /* every request, for the suite to read */

  const server = http.createServer((req, res) => {
    const u = new URL(req.url, 'http://localhost');
    seen.push({ method: req.method, path: u.pathname, query: u.search,
                auth: req.headers.authorization || null });

    /* ---- CORS, which is the whole point of this fixture ---- */
    const origin = req.headers.origin || '*';
    if (mode !== 'none') {
      res.setHeader('Access-Control-Allow-Origin', origin === 'null' ? '*' : origin);
      res.setHeader('Access-Control-Expose-Headers', 'location');
    }
    if (mode === 'full') {
      res.setHeader('Access-Control-Allow-Headers', 'Authorization, Content-Type');
      res.setHeader('Access-Control-Allow-Methods', 'GET, POST, PATCH, OPTIONS');
      res.setHeader('Access-Control-Max-Age', '600');
    }
    if (req.method === 'OPTIONS') {
      /* A default install answers the preflight without naming
         Authorization, so the browser refuses to send the real request.
         Answering 204 with no Allow-Headers is exactly that case. */
      res.statusCode = mode === 'none' ? 403 : 204;
      res.end();
      return;
    }

    const send = (code, body, type) => {
      res.statusCode = code;
      res.setHeader('Content-Type', type || 'application/json');
      res.end(typeof body === 'string' || Buffer.isBuffer(body) ? body : JSON.stringify(body));
    };

    const key = req.headers.authorization || '';
    const p = u.pathname.replace(/^\/api\/v2/, '');

    if (!key) return send(401, { code: '401', message: 'Unauthorized', description: 'No API key.' });
    if (key !== GOOD_KEY && key !== RO_KEY)
      return send(401, { code: '401', message: 'Unauthorized', description: 'Invalid API key.' });
    const mayWrite = !readonly && key === GOOD_KEY;

    if (p === '/info' && req.method === 'GET')
      return send(200, { elabftw_version: '5.2.0', ts_balance: 0, all_users_count: 3 });

    if (p === '/users/me' && req.method === 'GET')
      return send(200, { userid: 4, fullname: 'Yufeng Tong', email: 'tong@example.org',
                         teams: [{ id: 1, name: 'Tong lab' }] });

    if (p === '/apikeys' && req.method === 'GET')
      return send(200, [{ id: 1, name: 'PlasmidBench', can_write: mayWrite ? 1 : 0 }]);

    let m;
    /* list */
    if ((m = /^\/(items|experiments)$/.exec(p)) && req.method === 'GET') {
      const q = (u.searchParams.get('q') || '').toLowerCase();
      const limit = +(u.searchParams.get('limit') || 15);
      const rows = state[m[1]].filter(e => !q || (e.title || '').toLowerCase().includes(q));
      return send(200, rows.slice(0, limit).map(e => Object.assign({}, e, { uploads: undefined, _bytes: undefined })));
    }
    /* create */
    if ((m = /^\/items$/.exec(p)) && req.method === 'POST') {
      if (!mayWrite) return send(403, { code: '403', message: 'Forbidden',
                                        description: 'This key cannot write.' });
      const id = state.nextItemId++;
      let body = '';
      req.on('data', d => { body += d; });
      return req.on('end', () => {
        let j = {};
        try { j = JSON.parse(body || '{}'); } catch (e) { /* eLabFTW would 400 */ }
        state.items.push({ id: id, title: j.title || 'Untitled', body: j.body || '',
                           category_title: 'Plasmid', uploads: [], comments: [] });
        res.statusCode = 201;
        res.setHeader('location', 'http://localhost:' + server.address().port + '/api/v2/items/' + id);
        res.end(JSON.stringify({ id: id }));
      });
    }
    /* read one */
    if ((m = /^\/(items|experiments)\/(\d+)$/.exec(p)) && req.method === 'GET') {
      const e = state[m[1]].find(x => x.id === +m[2]);
      if (!e) return send(404, { code: '404', message: 'Not found', description: 'No such entry.' });
      return send(200, Object.assign({}, e, { uploads: undefined }));
    }
    /* uploads: list */
    if ((m = /^\/(items|experiments)\/(\d+)\/uploads$/.exec(p)) && req.method === 'GET') {
      const e = state[m[1]].find(x => x.id === +m[2]);
      if (!e) return send(404, { code: '404', message: 'Not found', description: 'No such entry.' });
      return send(200, (e.uploads || []).map(x => Object.assign({}, x, { _bytes: undefined })));
    }
    /* uploads: add */
    if ((m = /^\/(items|experiments)\/(\d+)\/uploads$/.exec(p)) && req.method === 'POST') {
      if (!mayWrite) return send(403, { code: '403', message: 'Forbidden',
                                        description: 'This key cannot write.' });
      const e = state[m[1]].find(x => x.id === +m[2]);
      if (!e) return send(404, { code: '404', message: 'Not found', description: 'No such entry.' });
      const chunks = [];
      req.on('data', d => chunks.push(d));
      return req.on('end', () => {
        const raw = Buffer.concat(chunks);
        /* enough multipart parsing to see the filename, the comment and the
           bytes, which is all the suite needs to check */
        const text = raw.toString('binary');
        const nameM = /name="file";\s*filename="([^"]*)"/.exec(text);
        const commM = /name="comment"\r\n\r\n([\s\S]*?)\r\n--/.exec(text);
        let bytes = Buffer.alloc(0);
        if (nameM) {
          const after = text.indexOf('\r\n\r\n', text.indexOf(nameM[0]));
          const bound = text.indexOf('\r\n--', after + 4);
          bytes = Buffer.from(text.slice(after + 4, bound < 0 ? undefined : bound), 'binary');
        }
        const up = { id: state.nextUploadId++, real_name: nameM ? nameM[1] : 'file',
                     filesize: bytes.length, comment: commM ? commM[1] : '',
                     hash: crypto.createHash('sha256').update(bytes).digest('hex'),
                     hash_algorithm: 'sha256', created_at: '2026-11-10 10:00:00', _bytes: bytes };
        e.uploads.push(up);
        res.statusCode = 201;
        res.setHeader('location', '/api/v2/' + m[1] + '/' + e.id + '/uploads/' + up.id);
        res.end(JSON.stringify({ id: up.id }));
      });
    }
    /* uploads: fetch one */
    if ((m = /^\/(items|experiments)\/(\d+)\/uploads\/(\d+)$/.exec(p)) && req.method === 'GET') {
      const e = state[m[1]].find(x => x.id === +m[2]);
      const up = e && (e.uploads || []).find(x => x.id === +m[3]);
      if (!up) return send(404, { code: '404', message: 'Not found', description: 'No such file.' });
      if (u.searchParams.get('format') === 'binary') {
        /* optionally corrupt it, so the checksum check can be tested */
        const body = opts.corruptDownloads ? Buffer.concat([up._bytes, Buffer.from('X')]) : up._bytes;
        return send(200, body, 'application/octet-stream');
      }
      return send(200, Object.assign({}, up, { _bytes: undefined }));
    }
    /* comments */
    if ((m = /^\/(items|experiments)\/(\d+)\/comments$/.exec(p)) && req.method === 'POST') {
      if (!mayWrite) return send(403, { code: '403', message: 'Forbidden',
                                        description: 'This key cannot write.' });
      const e = state[m[1]].find(x => x.id === +m[2]);
      if (!e) return send(404, { code: '404', message: 'Not found', description: 'No such entry.' });
      let body = '';
      req.on('data', d => { body += d; });
      return req.on('end', () => {
        let j = {};
        try { j = JSON.parse(body || '{}'); } catch (err) { /* ignore */ }
        e.comments.push({ id: e.comments.length + 1, comment: j.comment || '' });
        res.statusCode = 201;
        res.end(JSON.stringify({ id: e.comments.length }));
      });
    }

    return send(404, { code: '404', message: 'Not found', description: 'No such endpoint: ' + p });
  });

  return new Promise(resolve => {
    server.listen(opts.port || 0, '127.0.0.1', () => {
      const port = server.address().port;
      resolve({
        port: port,
        url: 'http://127.0.0.1:' + port,
        keys: { rw: GOOD_KEY, ro: RO_KEY },
        state: state,
        seen: seen,
        stop: () => new Promise(r => server.close(() => r()))
      });
    });
  });
}

module.exports = { start: start, GOOD_KEY: GOOD_KEY, RO_KEY: RO_KEY };

if (require.main === module) {
  const port = +(process.argv[2] || 8099);
  const mode = process.argv[3] || 'full';
  start({ port: port, cors: mode }).then(s => {
    console.log('mock eLabFTW on ' + s.url + '/api/v2  (cors=' + mode + ')');
    console.log('keys: read-write ' + s.keys.rw + ' / read-only ' + s.keys.ro);
  });
}
