/* ============================================================================
   PlasmidBench — eLabFTW.

   eLabFTW is the lab's record system and this is not one. The division of
   labour, from section 7 of the project evaluation:

     eLabFTW     the records: plasmids, primers, samples, experiments, and
                 the audit trail that makes them GLP records at all
     PlasmidBench  portable analysis, design and visualisation, aware of what
                 the freezer holds

   Which gives one rule that everything here follows:

     READ FREELY. WRITE ONLY BY ADDING, AND ONLY THROUGH THE API.

   Specifically, nothing in this file ever replaces an attachment, patches a
   body, or deletes anything. A result goes back as a NEW attachment plus a
   comment saying where it came from — so the original stays exactly as it
   was, the addition is attributable to the person whose key was used, and
   eLabFTW's own audit trail records it with a timestamp. A construct does
   not become verified lab stock because a design tool produced it, and the
   way to express that is to propose rather than to overwrite.

   THE API KEY
   It is the person's own, generated in their eLabFTW settings, and it is
   never embedded in the distributed page — there is nowhere in this file
   for one to live. It is typed in and kept in this browser, or for the
   session only, as they choose. A key with read-only permission is enough
   for everything except sending a result back, and the connection test says
   which one has been given.

   CORS
   The evaluation puts this first, and rightly: whether a browser may call
   eLabFTW's API at all is decided by the server's CORS configuration, and
   the answer settles whether this works with no new infrastructure or needs
   a proxy. The custom Authorization header makes every request preflighted,
   so a server that allows GET but not the preflight still fails. There is a
   test here that reports each step separately, so the answer is a reading
   rather than a guess. See PB.elab.test.

   API: eLabFTW REST v2. Base <host>/api/v2, key in the Authorization
   header with no prefix. Endpoints used, all from the published OpenAPI
   description: GET /info, GET /items, GET /experiments, GET /{type}/{id},
   GET /{type}/{id}/uploads, GET /{type}/{id}/uploads/{sub}?format=binary,
   POST /{type}/{id}/uploads (multipart, field "file"),
   POST /{type}/{id}/comments, POST /items.
   ============================================================================ */
(function (root) {
  'use strict';
  /** @type {any} */
  const PB = root.PB = root.PB || {};
  const B = root.Bio;

  /* ------------------------------------------------------------- settings -- */

  /* The key is held apart from the rest of the settings so that "remember
     it" and "this session only" are different places, and so that nothing
     that exports settings can carry a credential out with it. */
  let sessionKey = null;

  function config() {
    const c = PB.lsGet ? PB.lsGet('elab', {}) : {};
    return {
      url: (c && c.url) || '',
      key: sessionKey || (c && c.key) || '',
      remember: !!(c && c.key),
      scope: (c && c.scope) || 2,          /* 1 self, 2 team, 3 everything */
      user: (c && c.user) || null,         /* filled in by the test */
      canWrite: c ? c.canWrite : null
    };
  }

  function setConfig(next) {
    const c = PB.lsGet ? PB.lsGet('elab', {}) : {};
    const out = Object.assign({}, c);
    if (next.url !== undefined) out.url = String(next.url || '').replace(/\/+$/, '');
    if (next.scope !== undefined) out.scope = +next.scope || 2;
    if (next.user !== undefined) out.user = next.user;
    if (next.canWrite !== undefined) out.canWrite = next.canWrite;
    if (next.key !== undefined) {
      if (next.remember) { out.key = next.key; sessionKey = null; }
      else { delete out.key; sessionKey = next.key || null; }
    }
    if (PB.lsSet) PB.lsSet('elab', out);
    return config();
  }

  function forget() {
    sessionKey = null;
    if (PB.lsSet) PB.lsSet('elab', {});
  }

  /* --------------------------------------------------------------- calling -- */

  function base() {
    const c = config();
    if (!c.url) throw new Error('No eLabFTW address has been set.');
    /* someone will paste the address of the page they were looking at */
    let u = c.url.replace(/\/+$/, '');
    u = u.replace(/\/(experiments|database|items)\.php.*$/, '');
    if (!/\/api\/v2$/.test(u)) u += '/api/v2';
    return u;
  }
  PB.elabBase = base;

  /**
   * One call. Returns { ok, status, data, text, error, kind } and never
   * throws, because every caller has to be able to tell a refusal apart
   * from a network failure and say which happened.
   *
   * `kind` is the thing worth reporting: 'cors' when the browser would not
   * make the request at all, 'auth' for a rejected key, 'perm' for a key
   * without the permission, 'notfound', 'server', 'network', 'ok'.
   */
  async function call(path, opts) {
    opts = opts || {};
    const c = config();
    if (!c.key) return { ok: false, kind: 'auth', error: 'No API key has been given.' };
    let url;
    try { url = base() + path; }
    catch (e) { return { ok: false, kind: 'config', error: e.message }; }

    const headers = { Authorization: c.key };
    if (opts.json) headers['Content-Type'] = 'application/json';
    if (!opts.binary) headers.Accept = 'application/json';

    let res;
    try {
      res = await fetch(url, {
        method: opts.method || 'GET',
        headers: headers,
        body: opts.body || undefined,
        /* no cookies: the key is the whole of the authentication, and
           sending credentials would make the server's CORS requirements
           stricter for no gain */
        credentials: 'omit',
        mode: 'cors'
      });
    } catch (e) {
      /* A browser reports a blocked cross-origin request as a plain network
         failure and tells the page nothing more — by design. So this cannot
         claim it was CORS; it says what it knows and points at the test. */
      return { ok: false, kind: 'network',
               error: 'The browser could not reach ' + url + '. Either the server is not ' +
                      'answering, or it is refusing requests from this page. Run the ' +
                      'connection test, which separates the two.',
               detail: String(e && e.message || e) };
    }

    const kindOf = s => s === 401 || s === 403 ? (s === 401 ? 'auth' : 'perm')
                      : s === 404 ? 'notfound' : s >= 500 ? 'server' : 'http';
    if (!res.ok) {
      let body = '';
      try { body = await res.text(); } catch (e) { /* nothing to add */ }
      let msg = body;
      try { const j = JSON.parse(body); msg = j.description || j.message || body; } catch (e) { /* text */ }
      return { ok: false, status: res.status, kind: kindOf(res.status),
               error: msg || ('eLabFTW answered ' + res.status) };
    }
    if (opts.binary) {
      const buf = await res.arrayBuffer();
      return { ok: true, status: res.status, kind: 'ok', buffer: buf,
               type: res.headers.get('content-type') || '' };
    }
    const text = await res.text();
    let data = null;
    try { data = text ? JSON.parse(text) : null; } catch (e) { /* not JSON */ }
    return { ok: true, status: res.status, kind: 'ok', data: data, text: text,
             location: res.headers.get('location') || null };
  }
  PB.elabCall = call;

  /* ------------------------------------------------------- connection test -- */

  /**
   * The half-day question from section 7.3, as a button: can a browser talk
   * to this eLabFTW directly, and with what permissions?
   *
   * Each step is reported separately, because the useful answer is which
   * one failed. A server that allows a simple GET but not the preflight
   * that a custom Authorization header forces will pass step 1 and fail
   * step 2, and that distinction is the whole architecture decision.
   *
   * @returns {Promise<{steps:Array<{name:string, ok:boolean|null, note:string}>, verdict:string, canWrite:boolean|null, user:any}>}
   */
  async function test() {
    const steps = [];
    const add = (name, ok, note) => { steps.push({ name: name, ok: ok, note: note }); return ok; };
    const c = config();

    if (!c.url) {
      add('An address to try', false, 'Nothing has been entered yet.');
      return { steps: steps, verdict: 'nothing to test', canWrite: null, user: null };
    }
    let url;
    try { url = base(); } catch (e) { url = ''; }
    add('The address', !!url, url ? 'Calling ' + url : 'Could not be read as a URL.');

    /* 1 — is anything there at all? A request with no custom header is not
       preflighted, so this separates "the server is unreachable" from "the
       server will not accept our preflight". It is sent without the key,
       so a 401 is a success for this step: something answered. */
    let reachable = null;
    try {
      const r = await fetch(url + '/info', { method: 'GET', mode: 'cors', credentials: 'omit' });
      reachable = true;
      add('Something answers', true, 'HTTP ' + r.status +
          (r.status === 401 ? ' — which is the right answer to a request with no key.' : '.'));
    } catch (e) {
      reachable = false;
      add('Something answers', false,
          'The browser would not complete even a plain request. Either the host is wrong or ' +
          'unreachable, or the server sends no Access-Control-Allow-Origin header at all.');
    }

    /* 2 — the preflight. The Authorization header is not a CORS-safelisted
       one, so the browser sends OPTIONS first and will not proceed unless
       the server allows that header from this origin. This is the step a
       default eLabFTW install fails. */
    if (reachable) {
      const r = await call('/info', {});
      if (r.ok) {
        add('A request carrying the API key', true,
            'The server allows the Authorization header from this page. Direct browser ' +
            'access works here — no proxy is needed.');
        add('eLabFTW version', true, (r.data && (r.data.elabftw_version || r.data.version)) || 'answered');
      } else if (r.kind === 'network') {
        add('A request carrying the API key', false,
            'The plain request worked and this one did not, which leaves one cause: the ' +
            'server does not allow the Authorization header from this origin. It needs ' +
            'Access-Control-Allow-Headers to include Authorization, and ' +
            'Access-Control-Allow-Origin to include this page. Until then, use the file ' +
            'exchange route or put a small proxy in front.');
        return { steps: steps, verdict: 'blocked by CORS', canWrite: null, user: null };
      } else if (r.kind === 'auth') {
        add('A request carrying the API key', false,
            'The server took the request and rejected the key (401). CORS is fine; the key is not.');
        return { steps: steps, verdict: 'the key was rejected', canWrite: null, user: null };
      } else {
        add('A request carrying the API key', false, r.error || ('HTTP ' + r.status));
        return { steps: steps, verdict: 'the server refused', canWrite: null, user: null };
      }
    } else {
      return { steps: steps, verdict: 'nothing answered', canWrite: null, user: null };
    }

    /* 3 — who the key belongs to. Attribution comes from eLabFTW, and this
       is where it is read; nothing here invents an author. */
    let user = null;
    const me = await call('/users/me', {});
    if (me.ok && me.data) {
      user = { id: me.data.userid || me.data.id, name: me.data.fullname || me.data.firstname || '',
               email: me.data.email || '', teams: me.data.teams || null };
      add('Whose key it is', true, user.name || ('user ' + user.id) +
          ' — anything sent back will be recorded against this account.');
    } else {
      add('Whose key it is', null, 'This server did not say. Not a problem for reading.');
    }

    /* 4 — can it read records? */
    const items = await call('/items?limit=1&scope=' + (config().scope || 2), {});
    add('Reading resources', items.ok,
        items.ok ? 'Read ' + ((items.data || []).length) + ' of them, which is enough to list.'
                 : (items.error || 'refused'));

    /* 5 — can it write? Asked without writing anything: a key's write
       permission is a property of the key, and eLabFTW reports it. */
    let canWrite = null;
    const keys = await call('/apikeys', {});
    if (keys.ok && Array.isArray(keys.data)) {
      /* canWrite here means some key on this account can write; it does not
         single out the one in use, which eLabFTW does not reveal. The real
         answer comes when something is sent, and the attempt reports it. */
      const anyRW = keys.data.some(k => (k.can_write === true || k.can_write === 1));
      canWrite = anyRW;
      add('Permission to add things', anyRW,
          anyRW ? 'At least one key on this account may write. Sending a result back adds a ' +
                  'new attachment and a comment; nothing is ever replaced.'
                : 'Every key on this account is read-only. Reading works; sending results back ' +
                  'needs a read-write key, which you can make in your eLabFTW settings.');
    } else {
      add('Permission to add things', null,
          'Could not be read in advance. It will be clear the first time something is sent, ' +
          'and nothing is replaced either way.');
    }

    setConfig({ user: user, canWrite: canWrite });
    return { steps: steps, verdict: 'direct browser access works', canWrite: canWrite, user: user };
  }

  /* ------------------------------------------------------------- reading --- */

  /** Search resources ('items') or experiments. */
  async function search(q, type, limit) {
    const t = type === 'experiments' ? 'experiments' : 'items';
    const c = config();
    const qs = '?limit=' + (limit || 20) + '&scope=' + (c.scope || 2) +
               (q ? '&q=' + encodeURIComponent(q) : '');
    const r = await call('/' + t + qs, {});
    if (!r.ok) return r;
    return { ok: true, type: t, rows: (r.data || []).map(e => ({
      id: e.id, title: e.title || '(untitled)', elabid: e.elabid || null,
      category: e.category_title || e.category || '', date: e.date || e.created_at || '',
      owner: e.fullname || '', tags: e.tags || ''
    })) };
  }

  /** The files attached to one record. */
  async function uploads(type, id) {
    const r = await call('/' + type + '/' + id + '/uploads', {});
    if (!r.ok) return r;
    return { ok: true, rows: (r.data || []).map(u => ({
      id: u.id, name: u.real_name || u.long_name || ('upload ' + u.id),
      size: u.filesize || 0, comment: u.comment || '', hash: u.hash || '',
      hashAlgo: u.hash_algorithm || '', created: u.created_at || ''
    })) };
  }

  /**
   * Fetch one attachment's bytes. The checksum eLabFTW stored with the file
   * is checked against the bytes that arrived when the algorithm is one we
   * have, so a truncated download is noticed here rather than becoming a
   * strange-looking plasmid later.
   */
  async function download(type, id, uploadId, expectHash, algo) {
    const r = await call('/' + type + '/' + id + '/uploads/' + uploadId + '?format=binary',
                         { binary: true });
    if (!r.ok) return r;
    const bytes = new Uint8Array(r.buffer);
    let verified = null;
    if (expectHash && /sha ?-?256/i.test(algo || '') && B && B.sha256) {
      const got = B.sha256(bytes);
      verified = got.toLowerCase() === String(expectHash).toLowerCase();
      if (!verified) {
        return { ok: false, kind: 'checksum',
                 error: 'The file that arrived does not match the checksum eLabFTW holds for ' +
                        'it (' + got.slice(0, 12) + '… against ' +
                        String(expectHash).slice(0, 12) + '…). It was not opened.' };
      }
    }
    return { ok: true, buffer: r.buffer, bytes: bytes, type: r.type, verified: verified };
  }

  /* ------------------------------------------------------------- adding ---- */

  /**
   * Send a result back as a NEW attachment, with a comment saying what it
   * is and where it came from. Nothing is replaced and nothing is patched:
   * see the rule at the top of this file.
   *
   * @param {string} type 'items' or 'experiments'
   * @param {number|string} id
   * @param {string} filename
   * @param {string|Blob} content
   * @param {{comment?:string, checksum?:string, note?:string}} [meta]
   */
  async function attach(type, id, filename, content, meta) {
    meta = meta || {};
    if (typeof FormData === 'undefined') {
      return { ok: false, kind: 'config', error: 'This browser cannot send a file.' };
    }
    const blob = content instanceof Blob ? content
      : new Blob([content], { type: 'text/plain' });
    const fd = new FormData();
    fd.append('file', blob, filename);
    /* The comment is part of the record, so it says what made the file and
       what it is of — a file called "product.gb" with no provenance is a
       file nobody will trust in a year. */
    const stamp = new Date().toISOString().slice(0, 16).replace('T', ' ');
    const checksum = meta.checksum ||
      (typeof content === 'string' && B && B.sha256 ? 'sha256:' + B.sha256(content) : '');
    fd.append('comment', (meta.comment || 'From PlasmidBench') +
      ' · added ' + stamp + (checksum ? ' · ' + checksum : ''));

    const r = await call('/' + type + '/' + id + '/uploads', { method: 'POST', body: fd });
    if (!r.ok) {
      if (r.kind === 'perm') {
        return { ok: false, kind: 'perm',
                 error: 'That key may read but not add. Make a read-write key in your eLabFTW ' +
                        'settings, or download the file and attach it there yourself — ' +
                        'nothing has been changed.' };
      }
      return r;
    }
    /* and a comment on the record itself, so the addition is visible in the
       entry rather than only in its file list */
    if (meta.note) {
      await call('/' + type + '/' + id + '/comments', {
        method: 'POST', json: true,
        body: JSON.stringify({ comment: meta.note })
      });
    }
    return { ok: true, location: r.location, checksum: checksum };
  }

  /**
   * Propose a construct as a new resource. Deliberately a *proposal*: the
   * title is marked, the body says what made it and from what, and the
   * sequence goes on as a GenBank attachment with its checksum. Whether it
   * becomes real lab stock is a decision for a person in eLabFTW, which is
   * where that decision belongs.
   */
  async function propose(rec, opts) {
    opts = opts || {};
    const gb = B.toGenBank(rec);
    const sum = B.seqChecksum(rec.seq, { topology: rec.topology });
    const title = (opts.title || rec.name || 'construct');
    const r = await call('/items', {
      method: 'POST', json: true,
      body: JSON.stringify({
        title: '[proposed] ' + title,
        /* category is a team setting, so it is passed through rather than
           guessed; without one eLabFTW uses the team's default */
        category_id: opts.categoryId || undefined,
        body: [
          '<p><strong>Proposed by PlasmidBench</strong>, not yet verified lab stock.</p>',
          '<ul>',
          '<li>' + (rec.seq || '').length + ' bp, ' + (rec.topology || 'linear') + '</li>',
          '<li>Checksum of the sequence: <code>' + sum + '</code></li>',
          opts.from ? '<li>Built from: ' + esc(opts.from) + '</li>' : '',
          opts.how ? '<li>How: ' + esc(opts.how) + '</li>' : '',
          '</ul>',
          '<p>The GenBank file is attached. Sequence it before treating this as stock.</p>'
        ].join('')
      })
    });
    if (!r.ok) return r;
    /* the location header names the new resource */
    const m = /\/items\/(\d+)/.exec(r.location || '');
    const id = m ? +m[1] : null;
    if (id) {
      await attach('items', id, (rec.name || 'construct').replace(/[^\w.-]+/g, '_') + '.gb', gb,
        { comment: 'GenBank of the proposed construct', checksum: sum });
    }
    return { ok: true, id: id, checksum: sum };
  }

  function esc(s) {
    return String(s == null ? '' : s).replace(/&/g, '&amp;').replace(/</g, '&lt;')
      .replace(/>/g, '&gt;').replace(/"/g, '&quot;');
  }

  PB.elab = {
    config: config, setConfig: setConfig, forget: forget,
    test: test, search: search, uploads: uploads, download: download,
    attach: attach, propose: propose, call: call, base: base
  };

  if (typeof module !== 'undefined' && module.exports) module.exports = PB.elab;
})(typeof window !== 'undefined' ? window : globalThis);
