/* ============================================================================
   PlasmidBench — the eLabFTW panel.

   Three things, in the order a person needs them:

     1. Say where eLabFTW is and give a key, then TEST it. The test is first
        because the answer decides whether any of this works at all, and
        because the answer is a property of the server, not of this page: a
        default install does not allow a browser to send an Authorization
        header cross-origin, and until that is changed the file-exchange
        route is the one that works. The test says which case you are in,
        in those words.
     2. Find a record and open what is attached to it. Reading is the
        common case and it is also the safe one.
     3. Send a result back — as an added attachment and a comment, never as
        a replacement. The panel says so on the button, because "save to
        eLabFTW" would imply something this must not do.

   The key is never in the distributed file. It is typed here and kept in
   this browser, or for the session only.
   ============================================================================ */
(function (root) {
  'use strict';
  /** @type {any} */
  const PB = root.PB;
  const B = root.Bio;
  const $ = PB.$, esc = PB.esc, fmt = PB.fmt;

  /* what the panel is looking at, between renders */
  let found = [];       /* search results */
  let openOn = null;    /* { type, id, title } */
  let files = [];       /* its attachments */

  function head() {
    const c = PB.elab.config();
    const linked = c.url && c.key;
    return '<div class="row"><div class="field wide">' +
      '<label>Your eLabFTW address</label>' +
      '<input id="elUrl" type="text" spellcheck="false" placeholder="https://eln.your-lab.org" value="' +
        esc(c.url) + '">' +
      '<div class="hint">The address you log in at. The <code>/api/v2</code> part is added for you.</div>' +
      '</div></div>' +
      '<div class="row"><div class="field wide">' +
      '<label>Your personal API key</label>' +
      '<input id="elKey" type="password" spellcheck="false" autocomplete="off" placeholder="' +
        (c.key ? '•••• already set' : 'from eLabFTW → Settings → API keys') +
        '">' +
      '<div class="hint">Yours, not the lab’s. It stays in this browser and is never in the ' +
      'file you downloaded. A read-only key is enough for everything but sending results back.</div>' +
      '</div></div>' +
      '<div class="row">' +
      '<div class="field"><label>Remember it</label><select id="elRemember">' +
        '<option value="0"' + (c.remember ? '' : ' selected') + '>no — this session only</option>' +
        '<option value="1"' + (c.remember ? ' selected' : '') + '>yes — keep it in this browser</option>' +
      '</select></div>' +
      '<div class="field"><label>What to search</label><select id="elScope">' +
        '<option value="1"' + (c.scope === 1 ? ' selected' : '') + '>only mine</option>' +
        '<option value="2"' + (c.scope === 2 ? ' selected' : '') + '>my team</option>' +
        '<option value="3"' + (c.scope === 3 ? ' selected' : '') + '>everything I can see</option>' +
      '</select></div>' +
      '<div class="field"><label>&nbsp;</label>' +
        '<button class="btn primary" id="elTest">Test the connection</button></div>' +
      '<div class="field"><label>&nbsp;</label>' +
        '<button class="btn" id="elForget">Forget the key</button></div>' +
      '</div>' +
      (c.user ? '<div class="note ok">Connected as <strong>' + esc(c.user.name || ('user ' + c.user.id)) +
        '</strong>. Anything added is recorded against that account.</div>' : '') +
      '<div id="elTestOut"></div>' +
      (linked ? searchBox() : '<div class="note">Enter an address and a key, then test. ' +
        'Nothing is sent anywhere until you do.</div>') +
      fileExchangeNote();
  }

  /* The route that needs no permission, no configuration and no network
     access from this page at all. It is written down because it is the
     right first step and because it keeps working when the API route is
     blocked. */
  function fileExchangeNote() {
    return '<details style="margin-top:14px"><summary>If the API is blocked — the ' +
      'file route, which always works</summary><div class="note">' +
      '<ol style="margin:4px 0 0 18px;padding:0">' +
      '<li>In eLabFTW, download the GenBank or FASTA attached to the record.</li>' +
      '<li>Open it here with <strong>Import</strong>, or drop it on the window.</li>' +
      '<li>Do the work. Export the result — GenBank keeps the features.</li>' +
      '<li>In eLabFTW, attach that file to the record, or to a new entry, and say in a ' +
      'comment what made it.</li>' +
      '</ol>' +
      '<p style="margin:8px 0 0">Slower by a few clicks, and the audit trail is identical: ' +
      'the file arrives through eLabFTW’s own upload, attributed to you, timestamped by it. ' +
      'This is the proof of concept section 7.3 asks for, and it needs nothing installing.</p>' +
      '</div></details>';
  }

  function searchBox() {
    return '<hr style="border:0;border-top:1px solid var(--line);margin:14px 0">' +
      '<div class="row">' +
      '<div class="field wide"><label>Find a record</label>' +
        '<input id="elQ" type="text" placeholder="plasmid name, or part of a title">' +
      '</div>' +
      '<div class="field"><label>Where</label><select id="elType">' +
        '<option value="items">resources (plasmids, primers, strains)</option>' +
        '<option value="experiments">experiments</option>' +
      '</select></div>' +
      '<div class="field"><label>&nbsp;</label><button class="btn primary" id="elGo">Search</button></div>' +
      '</div><div id="elResults"></div>';
  }

  function renderResults() {
    const host = $('#elResults');
    if (!host) return;
    if (!found.length) { host.innerHTML = '<div class="hint">Nothing found yet.</div>'; return; }
    host.innerHTML = '<div class="scrollx"><table class="tbl sm"><thead><tr>' +
      '<th>Title</th><th>Category</th><th>Date</th><th>Who</th><th></th></tr></thead><tbody>' +
      found.map((r, i) => '<tr><td>' + esc(r.title) + '</td><td>' + esc(r.category) +
        '</td><td class="mono">' + esc(String(r.date).slice(0, 10)) + '</td><td>' + esc(r.owner) +
        '</td><td><button class="btn sm" data-elopen="' + i + '">Open its files</button></td></tr>').join('') +
      '</tbody></table></div>' + '<div id="elFiles"></div>';
  }

  function renderFiles() {
    const host = $('#elFiles');
    if (!host) return;
    if (!openOn) { host.innerHTML = ''; return; }
    const c = PB.elab.config();
    let html = '<div class="note"><strong>' + esc(openOn.title) + '</strong> — ' +
      files.length + ' file' + (files.length === 1 ? '' : 's') + ' attached.</div>';
    if (files.length) {
      html += '<div class="scrollx"><table class="tbl sm"><thead><tr><th>File</th>' +
        '<th class="num">Size</th><th>Comment</th><th>Checksum</th><th></th></tr></thead><tbody>' +
        files.map((f, i) => '<tr><td class="mono">' + esc(f.name) + '</td>' +
          '<td class="num">' + fmt(f.size) + '</td><td>' + esc(f.comment) + '</td>' +
          '<td class="mono" title="' + esc(f.hashAlgo + ' ' + f.hash) + '">' +
            (f.hash ? esc(f.hash.slice(0, 10)) + '…' : '—') + '</td>' +
          '<td>' + (openable(f.name)
            ? '<button class="btn sm" data-elget="' + i + '">Open here</button>'
            : '<span class="hint">not a sequence file</span>') + '</td></tr>').join('') +
        '</tbody></table></div>';
    }
    /* sending something back, named for what it does */
    const rec = PB.active && PB.active();
    html += '<hr style="border:0;border-top:1px solid var(--line);margin:12px 0">';
    if (!rec) {
      html += '<div class="hint">Open a construct here and you can add it to this record.</div>';
    } else {
      html += '<div class="row"><div class="field wide"><label>Add to this record</label>' +
        '<div class="hint">Adds <strong>' + esc(rec.name) + '</strong> as a new GenBank attachment, ' +
        'with its checksum, and a comment saying what made it. Nothing already there is changed ' +
        'or replaced.' + (c.canWrite === false ? ' <strong>The keys on this account are ' +
        'read-only, so this will be refused.</strong>' : '') + '</div></div>' +
        '<div class="field"><label>&nbsp;</label><button class="btn primary" id="elSend">' +
        'Attach it</button></div></div>';
    }
    host.innerHTML = html;
  }

  const openable = n => /\.(gb|gbk|genbank|fa|fasta|fna|seq|embl|dna|ab1|gff3?|json|xml|txt)$/i.test(n || '');

  function renderTest(result) {
    const host = $('#elTestOut');
    if (!host) return;
    const icon = s => s.ok === true ? '✓' : s.ok === false ? '✗' : '–';
    const cls = result.verdict === 'direct browser access works' ? 'ok'
      : /blocked|rejected|refused|nothing/.test(result.verdict) ? 'bad' : 'warn';
    host.innerHTML = '<div class="note ' + cls + '"><strong>' + esc(result.verdict) + '</strong></div>' +
      '<div class="scrollx"><table class="tbl sm"><tbody>' +
      result.steps.map(s => '<tr><td style="width:26px;text-align:center;font-weight:700;color:' +
        (s.ok === true ? 'var(--ok)' : s.ok === false ? 'var(--bad)' : 'var(--muted)') + '">' +
        icon(s) + '</td><td style="white-space:nowrap">' + esc(s.name) + '</td><td>' +
        esc(s.note) + '</td></tr>').join('') +
      '</tbody></table></div>' +
      (result.verdict === 'blocked by CORS'
        ? '<div class="note warn"><strong>What to change, if you administer the server.</strong> ' +
          'eLabFTW has to send <code>Access-Control-Allow-Origin</code> naming this page’s ' +
          'origin and <code>Access-Control-Allow-Headers: Authorization</code>, and answer the ' +
          '<code>OPTIONS</code> preflight. In an nginx front end that is a few lines. ' +
          'Until then the file route below works with nothing installed, and is the honest ' +
          'first proof of concept.</div>'
        : '');
  }

  /* ----------------------------------------------------------------- wiring */

  PB.toolElab = function () {
    const w = PB.modal({
      title: 'eLabFTW',
      sub: 'read records, open their files, add results back',
      wide: true,
      body: '<div class="note">eLabFTW holds the lab’s records; this holds none of them. ' +
        'Files are read from it and results are <em>added</em> to it — never replaced, ' +
        'never patched — so its audit trail records every change, attributed to you.</div>' +
        head()
    });

    const cfgFromForm = () => {
      const keyEl = /** @type {HTMLInputElement} */ ($('#elKey'));
      const next = {
        url: /** @type {HTMLInputElement} */ ($('#elUrl')).value,
        scope: +(/** @type {HTMLSelectElement} */ ($('#elScope')) || { value: 2 }).value,
        remember: (/** @type {HTMLSelectElement} */ ($('#elRemember')) || { value: '0' }).value === '1'
      };
      if (keyEl && keyEl.value) next.key = keyEl.value;
      return PB.elab.setConfig(next);
    };

    const rewire = () => {
      /* the panel is rebuilt in place after a change of state */
      const body = w.querySelector('.modal-body') || w;
      body.innerHTML = body.innerHTML;
    };
    void rewire;

    w.addEventListener('click', async ev => {
      const t = /** @type {HTMLElement} */ (ev.target);
      if (!t) return;

      if (t.id === 'elTest') {
        cfgFromForm();
        const out = $('#elTestOut');
        if (out) out.innerHTML = '<div class="note">Testing…</div>';
        const r = await PB.elab.test();
        renderTest(r);
        /* a successful test reveals the search box, which was not there */
        if (r.verdict === 'direct browser access works' && !$('#elQ')) {
          out.insertAdjacentHTML('afterend', searchBox());
        }
        return;
      }

      if (t.id === 'elForget') {
        PB.elab.forget();
        found = []; openOn = null; files = [];
        PB.toast('The key and address are forgotten.', 'ok');
        PB.closeModal(); PB.toolElab();
        return;
      }

      if (t.id === 'elGo') {
        cfgFromForm();
        const q = (/** @type {HTMLInputElement} */ ($('#elQ')) || { value: '' }).value.trim();
        const type = (/** @type {HTMLSelectElement} */ ($('#elType')) || { value: 'items' }).value;
        const host = $('#elResults');
        if (host) host.innerHTML = '<div class="note">Searching…</div>';
        const r = await PB.elab.search(q, type);
        if (!r.ok) { if (host) host.innerHTML = '<div class="note bad">' + esc(r.error) + '</div>'; return; }
        found = r.rows; openOn = null; files = [];
        found.forEach(x => { x._type = r.type; });
        renderResults();
        return;
      }

      const openIdx = t.getAttribute && t.getAttribute('data-elopen');
      if (openIdx != null) {
        const r = found[+openIdx];
        openOn = { type: r._type, id: r.id, title: r.title };
        const host = $('#elFiles');
        if (host) host.innerHTML = '<div class="note">Reading its files…</div>';
        const u = await PB.elab.uploads(openOn.type, openOn.id);
        files = u.ok ? u.rows : [];
        if (!u.ok && host) { host.innerHTML = '<div class="note bad">' + esc(u.error) + '</div>'; return; }
        renderFiles();
        return;
      }

      const getIdx = t.getAttribute && t.getAttribute('data-elget');
      if (getIdx != null) {
        const f = files[+getIdx];
        PB.toast('Fetching ' + f.name + '…');
        const d = await PB.elab.download(openOn.type, openOn.id, f.id, f.hash, f.hashAlgo);
        if (!d.ok) { PB.toast(d.error, 'bad'); return; }
        /* hand it to the ordinary import path, which already knows every
           format this reads, so there is one importer and not two */
        try {
          /* handed to the ordinary import path, which already knows every
             format this reads: one importer, not two */
          const file = new File([d.buffer], f.name);
          PB.handleFiles([file]);
          await new Promise(r => setTimeout(r, 350));
          const note = 'from eLabFTW ' + openOn.type.replace(/s$/, '') + ' ' + openOn.id +
                       (d.verified ? ', checksum verified' : '');
          const rec = PB.active && PB.active();
          if (rec) {
            /* where it came from, on the record, so the link back is not
               only in someone's memory */
            rec.source = note;
            rec.elab = { type: openOn.type, id: openOn.id, uploadId: f.id,
                         checksum: f.hash ? f.hashAlgo + ':' + f.hash : null };
            PB.persist(rec);
          }
          PB.closeModal();
          PB.toast('Opened ' + f.name + ' — ' + note, 'ok');
        } catch (e) {
          PB.toast('That file arrived but could not be read: ' + (e && e.message), 'bad');
        }
        return;
      }

      if (t.id === 'elSend') {
        const rec = PB.active && PB.active();
        if (!rec || !openOn) return;
        const gb = B.toGenBank(rec);
        const sum = B.seqChecksum(rec.seq, { topology: rec.topology });
        const fname = (rec.name || 'construct').replace(/[^\w.-]+/g, '_') + '.gb';
        PB.toast('Adding ' + fname + '…');
        const r = await PB.elab.attach(openOn.type, openOn.id, fname, gb, {
          comment: 'GenBank from PlasmidBench',
          checksum: sum,
          note: 'Added <strong>' + esc(fname) + '</strong> from PlasmidBench: ' +
                (rec.seq || '').length + ' bp ' + (rec.topology || 'linear') +
                ', checksum <code>' + sum + '</code>. Not verified stock — sequence it first.'
        });
        if (!r.ok) { PB.toast(r.error, 'bad'); return; }
        rec.elab = Object.assign({}, rec.elab, { type: openOn.type, id: openOn.id, checksum: sum });
        PB.persist(rec);
        const u = await PB.elab.uploads(openOn.type, openOn.id);
        files = u.ok ? u.rows : files;
        renderFiles();
        PB.toast('Added to ' + openOn.title + '. Nothing there was changed.', 'ok');
        return;
      }
    });

    /* Enter in the search box searches */
    w.addEventListener('keydown', ev => {
      if (ev.key === 'Enter' && /** @type {HTMLElement} */ (ev.target).id === 'elQ') {
        ev.preventDefault();
        const b = $('#elGo'); if (b) b.click();
      }
    });
  };
})(typeof window !== 'undefined' ? window : globalThis);
