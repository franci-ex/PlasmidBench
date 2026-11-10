/* ============================================================================
   PlasmidBench — the lab profile panel.

   One screen holding everything a lab has to set before this page fits the
   way they work, and one file that carries it to the next lab or the next
   student. See profile.js for what is in a profile and, more importantly,
   what is kept out of it: no API keys and no records.

   Importing is a two-step: read the file, see what it would change, then
   take it. Not because a confirmation step is good manners, but because
   these settings change results silently — a profile that sets the salt
   concentration changes every melting temperature on the page, and someone
   who did not know that happened would have no way to account for it.
   ============================================================================ */
(function (root) {
  'use strict';
  /** @type {any} */
  const PB = root.PB;
  const $ = PB.$, esc = PB.esc;

  const SUPPLIER_LABELS = () => {
    const names = root.Bio.SUPPLIER_NAMES || {};
    return Object.keys(names).sort((a, b) => names[a].localeCompare(names[b]))
      .map(code => [code, names[code]]);
  };

  function field(key, what) {
    const v = PB.lsGet(key, '');
    const id = 'pf_' + key.replace(/[^\w]/g, '_');

    /* the three that are not single values get their own control */
    if (key === 'tmCond') {
      const c = PB.lsGet('tmCond', { oligo: 0.25, Na: 50, Mg: 1.5, dNTP: 0.2 });
      return '<div class="field wide"><label>' + esc(what) + '</label><div class="row">' +
        ['oligo', 'Na', 'Mg', 'dNTP'].map(k =>
          '<div class="field narrow"><label>' + k +
          (k === 'oligo' ? ' µM' : k === 'dNTP' ? ' mM' : ' mM') + '</label>' +
          '<input class="inp mono" type="number" step="0.05" min="0" id="pf_tm_' + k +
          '" value="' + esc(String(c[k])) + '"></div>').join('') +
        '</div></div>';
    }
    if (key === 'suppliers') {
      const chosen = new Set(PB.lsGet('suppliers', []) || []);
      return '<div class="field wide"><label>' + esc(what) + '</label>' +
        '<div class="chips" id="pf_suppliers">' +
        SUPPLIER_LABELS().map(([code, name]) =>
          '<label class="chip' + (chosen.has(code) ? ' on' : '') + '">' +
          '<input type="checkbox" value="' + code + '"' + (chosen.has(code) ? ' checked' : '') +
          '> ' + esc(name) + '</label>').join('') +
        '</div><div class="hint">' + (chosen.size
          ? 'The enzyme list shows only what these companies sell. Enzymes you added by hand are always kept.'
          : 'Nothing selected, so every enzyme in REBASE is offered.') + '</div></div>';
    }
    if (key === 'customEnzymes' || /\.columns$/.test(key)) {
      const n = Array.isArray(v) ? v.length : (v && typeof v === 'object' ? Object.keys(v).length : 0);
      return '<div class="field wide"><label>' + esc(what) + '</label>' +
        '<div class="hint">' + (n ? n + ' set here, and carried in the profile.'
          : 'Nothing set. This fills in when you add an enzyme or import a stock sheet.') +
        '</div></div>';
    }
    if (key === 'theme') {
      return '<div class="field"><label>' + esc(what) + '</label><select class="inp" id="' + id + '">' +
        ['system', 'light', 'dark'].map(t => '<option value="' + t + '"' +
          (PB.lsGet('theme', 'system') === t ? ' selected' : '') + '>' + t + '</option>').join('') +
        '</select></div>';
    }

    const numeric = /good|dead|armLen/.test(key);
    return '<div class="field' + (/name|pi|contact|facility|url|polymerase|kit/.test(key) ? ' wide' : '') +
      '"><label>' + esc(what) + '</label>' +
      '<input class="inp' + (numeric ? ' mono' : '') + '" type="' + (numeric ? 'number' : 'text') +
      '" id="' + id + '" value="' + esc(String(v == null ? '' : v)) + '"' +
      (key === 'elab.url' ? ' placeholder="https://eln.your-lab.org"' : '') + '></div>';
  }

  function body() {
    const groupTitles = {
      who: 'The lab',
      bench: 'The bench — these change what every calculation returns',
      enzymes: 'Enzymes you can actually buy',
      sequencing: 'What your sequencing facility really gives you',
      display: 'How it looks',
      records: 'Records',
      stock: 'Your own spreadsheet layouts'
    };
    return '<div class="note">Everything here is this lab’s, not built into the page. ' +
      'Export it and another lab can take the whole setup in one file — ' +
      '<strong>no API keys and no records travel with it</strong>, only configuration.</div>' +
      PB.profile.FIELDS.map(([group, rows]) =>
        '<h4 style="margin:16px 0 6px">' + esc(groupTitles[group] || group) + '</h4>' +
        '<div class="row">' + rows.map(([k, what]) => field(k, what)).join('') + '</div>'
      ).join('') +
      '<hr style="border:0;border-top:1px solid var(--line);margin:16px 0">' +
      '<div class="row">' +
      '<div class="field"><label>&nbsp;</label><button class="btn primary" id="pfSave">Save</button></div>' +
      '<div class="field"><label>&nbsp;</label><button class="btn" id="pfExport">Export this profile</button></div>' +
      '<div class="field"><label>&nbsp;</label><button class="btn" id="pfImport">Import a profile…</button></div>' +
      '</div><div id="pfOut"></div>' +
      '<input type="file" id="pfFile" accept=".json,application/json" style="display:none">';
  }

  function read() {
    const out = {};
    const get = id => { const el = $('#' + id); return el ? /** @type {HTMLInputElement} */ (el).value : null; };
    PB.profile.FIELDS.forEach(([, rows]) => rows.forEach(([k]) => {
      if (k === 'tmCond') {
        const c = {};
        ['oligo', 'Na', 'Mg', 'dNTP'].forEach(x => {
          const v = parseFloat(get('pf_tm_' + x));
          if (isFinite(v)) c[x] = v;
        });
        if (Object.keys(c).length === 4) out.tmCond = c;
        return;
      }
      if (k === 'suppliers') {
        const host = $('#pf_suppliers');
        if (host) out.suppliers = Array.from(host.querySelectorAll('input:checked'))
          .map(i => /** @type {HTMLInputElement} */ (i).value);
        return;
      }
      if (k === 'customEnzymes' || /\.columns$/.test(k)) return;    /* shown, not edited here */
      const v = get('pf_' + k.replace(/[^\w]/g, '_'));
      if (v == null) return;
      out[k] = /good|dead|armLen/.test(k) ? (v === '' ? '' : +v) : v;
    }));
    return out;
  }

  function save() {
    const vals = read();
    Object.keys(vals).forEach(k => {
      if (k === 'elab.url') { if (PB.elab) PB.elab.setConfig({ url: vals[k] }); return; }
      if (vals[k] === '' ) { PB.lsSet(k, ''); return; }
      PB.lsSet(k, vals[k]);
    });
    PB.reloadSettings();
    PB.renderAll();
    PB.toast('Profile saved. Melting temperatures and the enzyme list follow it from now on.', 'ok');
  }

  function showReport(rep) {
    const host = $('#pfOut');
    if (!host) return;
    if (!rep.ok) { host.innerHTML = '<div class="note bad">' + esc(rep.error) + '</div>'; return; }
    if (!rep.changes.length) {
      host.innerHTML = '<div class="note ok">That profile (' + esc(rep.name) +
        ') matches what you already have. Nothing to change.</div>';
      return;
    }
    const show = v => v === undefined || v === null || v === ''
      ? '<span class="hint">not set</span>'
      : '<span class="mono">' + esc(typeof v === 'object' ? JSON.stringify(v) : String(v)) + '</span>';
    host.innerHTML = '<div class="note warn"><strong>' + esc(rep.name) + '</strong>' +
      (rep.exported ? ', exported ' + esc(rep.exported) : '') + ' would change ' +
      rep.changes.length + ' setting' + (rep.changes.length === 1 ? '' : 's') +
      '. Look through them before taking it: these change what the page calculates.</div>' +
      '<div class="scrollx"><table class="tbl sm"><thead><tr><th>Setting</th><th>Now</th>' +
      '<th>Would become</th></tr></thead><tbody>' +
      rep.changes.map(c => '<tr><td>' + esc(c.what) + '</td><td>' + show(c.from) +
        '</td><td>' + show(c.to) + '</td></tr>').join('') +
      '</tbody></table></div>' +
      (rep.refused.length ? '<div class="note">' + rep.refused.length +
        ' field' + (rep.refused.length === 1 ? '' : 's') + ' in that file ' +
        (rep.refused.length === 1 ? 'is' : 'are') + ' not configuration — a key, or ' +
        'somebody’s records — and ' + (rep.refused.length === 1 ? 'was' : 'were') +
        ' ignored: ' + esc(rep.refused.join(', ')) + '</div>' : '') +
      (rep.ignored.length ? '<div class="hint">' + rep.ignored.length +
        ' setting(s) this version does not know about were skipped.</div>' : '') +
      '<div class="row"><div class="field"><label>&nbsp;</label>' +
      '<button class="btn primary" id="pfApply">Take these ' + rep.changes.length +
      ' setting' + (rep.changes.length === 1 ? '' : 's') + '</button></div>' +
      '<div class="field"><label>&nbsp;</label><button class="btn" id="pfCancel">Leave mine alone</button></div>' +
      '</div>';
    PB._pendingProfile = rep;
  }

  PB.toolProfile = function () {
    const w = PB.modal({
      title: 'Lab profile',
      sub: PB.lsGet('lab.name', '') || 'not named yet',
      wide: true,
      body: body()
    });

    w.addEventListener('change', ev => {
      const t = /** @type {HTMLElement} */ (ev.target);
      if (t && t.closest && t.closest('#pf_suppliers')) {
        const lab = t.closest('.chip');
        if (lab) lab.classList.toggle('on', /** @type {HTMLInputElement} */ (t).checked);
      }
    });

    w.addEventListener('click', ev => {
      const t = /** @type {HTMLElement} */ (ev.target);
      if (!t || !t.id) return;
      if (t.id === 'pfSave') { save(); return; }
      if (t.id === 'pfExport') {
        save();
        PB.offerFile(PB.profile.filename(), PB.profile.toText(), 'lab profile');
        return;
      }
      if (t.id === 'pfImport') { const f = $('#pfFile'); if (f) f.click(); return; }
      if (t.id === 'pfApply') {
        const n = PB.profile.apply(PB._pendingProfile);
        PB._pendingProfile = null;
        PB.closeModal();
        PB.toast('Took ' + n + ' setting' + (n === 1 ? '' : 's') + ' from that profile.', 'ok');
        return;
      }
      if (t.id === 'pfCancel') {
        PB._pendingProfile = null;
        const host = $('#pfOut'); if (host) host.innerHTML = '';
        return;
      }
    });

    const fileEl = $('#pfFile');
    if (fileEl) fileEl.addEventListener('change', async ev => {
      const f = /** @type {HTMLInputElement} */ (ev.target).files[0];
      if (!f) return;
      const text = await f.text();
      showReport(PB.profile.inspect(text));
    });
  };
})(typeof window !== 'undefined' ? window : globalThis);
