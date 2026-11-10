/* ============================================================================
   PlasmidBench — the annotation panel.

   Two matchers, shown side by side, because they answer different questions
   and a person needs to know which one spoke:

     exact         this sequence is in the library, letter for letter. One
                   pass over the sequence however large the library is.
     approximate   something very like this is in the library: a variant, or
                   a fragment left behind by a cloning step. Reported with
                   the identity and how much of the library entry was found.

   Every row says which, and what it is resting on. "CMV promoter, 97%, all
   of it" and "CMV promoter, exact" are different findings, and a panel that
   printed both as "CMV promoter" would be hiding the only thing worth
   knowing about the second one.

   Nothing is added to the construct until the person chooses it. An
   annotation is a claim about someone's plasmid, and 27 built-in entries
   plus whatever library they imported is not grounds for writing it in
   without asking.
   ============================================================================ */
(function (root) {
  'use strict';
  /** @type {any} */
  const PB = root.PB;
  const B = root.Bio;
  const $ = PB.$, esc = PB.esc, fmt = PB.fmt;

  let lastResult = null;
  let index = null;          /* the k-mer index, kept between runs */
  let indexFor = null;       /* which library it was built for */

  function library() {
    const extra = PB.lsGet('featureLib', null);
    if (extra && extra.entries && extra.entries.length) {
      return { entries: B.FEATURE_LIB.concat(extra.entries),
               name: extra.source || 'imported', imported: extra.entries.length };
    }
    return { entries: B.FEATURE_LIB, name: 'built in', imported: 0 };
  }

  function ensureIndex(lib, k) {
    const key = lib.name + ':' + lib.entries.length + ':' + k;
    if (index && indexFor === key) return index;
    index = B.buildFeatureIndex(lib.entries, { k: k });
    indexFor = key;
    return index;
  }

  function row(h, i) {
    const pc = n => (n * 100).toFixed(n > 0.995 ? 0 : 1) + '%';
    const how = h.exact
      ? '<span class="badge ok">exact</span>'
      : '<span class="badge ' + (h.identity >= 0.95 ? 'neutral' : 'warn') + '">' +
        pc(h.identity) + ' identical</span>' +
        (h.gaps ? ' <span class="badge warn">' + h.gaps + ' bp indel</span>' : '');
    const cov = h.partial
      ? '<span class="badge warn" title="Only part of the library entry is here — ' +
        'often what a cloning step left behind">' + pc(h.coverage) + ' of it</span>'
      : '<span class="badge neutral">all of it</span>';
    return '<tr><td><input type="checkbox" data-anno="' + i + '" checked></td>' +
      '<td>' + esc(h.name) + (h.alsoMatched
        ? '<div class="hint">also matched ' + esc(h.alsoMatched.slice(0, 3).join(', ')) + '</div>' : '') +
      '</td><td>' + esc(h.type) + '</td>' +
      '<td class="mono num">' + fmt(h.start + 1) + '–' + fmt(h.end) + '</td>' +
      '<td class="mono">' + (h.strand < 0 ? '−' : '+') + '</td>' +
      '<td class="num">' + fmt(h.libLength) + '</td>' +
      '<td>' + how + ' ' + cov + (h.clippedByBand
        ? ' <span class="badge warn" title="The alignment reached the edge of the band, so a ' +
          'larger insertion or deletion than this may be the real explanation">at the band edge</span>'
        : '') + '</td>' +
      '<td class="hint">' + esc(h.source || 'built in') + '</td></tr>';
  }

  function render() {
    const host = $('#anOut');
    if (!host || !lastResult) return;
    const r = lastResult;
    if (!r.hits.length) {
      host.innerHTML = '<div class="note">Nothing in the library matched, at the settings above. ' +
        'A library of ' + fmt(r.searched) + ' entries is small — importing one is under ' +
        '<em>Import a library</em>.</div>';
      return;
    }
    const exact = r.hits.filter(h => h.exact).length;
    const partial = r.hits.filter(h => h.partial).length;
    host.innerHTML =
      '<div class="note">' + r.hits.length + ' found in ' + fmt(r.searched) + ' library entries, ' +
      'in ' + r.ms + ' ms. ' + exact + ' exact, ' + (r.hits.length - exact) + ' approximate' +
      (partial ? ', ' + partial + ' of them only part of the entry' : '') + '.</div>' +
      '<div class="scrollx"><table class="tbl sm"><thead><tr><th></th><th>Name</th><th>Type</th>' +
      '<th>Where</th><th>Strand</th><th class="num">Entry bp</th><th>Match</th><th>From</th>' +
      '</tr></thead><tbody>' + r.hits.map(row).join('') + '</tbody></table></div>' +
      '<div class="row"><div class="field"><label>&nbsp;</label>' +
      '<button class="btn primary" id="anAdd">Add the ticked ones to this construct</button></div>' +
      '<div class="field"><label>&nbsp;</label><button class="btn" id="anNone">Untick all</button></div>' +
      '</div>';
  }

  PB.toolAnnotate = function () {
    const rec = PB.active && PB.active();
    if (!rec) { PB.toast('Open a construct first.', 'warn'); return; }
    const lib = library();

    const w = PB.modal({
      title: 'Annotate',
      sub: rec.name + ' · ' + fmt(lib.entries.length) + ' library entries',
      wide: true,
      body:
        '<div class="note">Exact search finds what is in the library letter for letter. ' +
        'Approximate search also finds variants and fragments — a promoter one base off, ' +
        'a terminator a cloning step took the end from — and says by how much each ' +
        'differs. Nothing is written into the construct until you pick it.</div>' +
        (lib.imported ? '' : '<div class="note warn">Only the ' + lib.entries.length +
          ' built-in entries are loaded. That is a small library and it will miss most of ' +
          'what is in a real plasmid. <strong>Import one below</strong> — ApE ships a ' +
          'feature list in a format this reads, and FPbase has an open API for fluorescent ' +
          'proteins.</div>') +
        '<div class="row">' +
        '<div class="field narrow"><label>Seed length</label>' +
          '<input class="inp mono" id="anK" type="number" min="8" max="20" value="12">' +
          '<div class="hint">12 finds a variant 5% diverged; 14 is far quicker and slightly less ' +
          'sensitive. Measured, not guessed — see engine8.js.</div></div>' +
        '<div class="field narrow"><label>Least identity</label>' +
          '<input class="inp mono" id="anId" type="number" min="0.6" max="1" step="0.01" value="0.85"></div>' +
        '<div class="field narrow"><label>Least of the entry</label>' +
          '<input class="inp mono" id="anCov" type="number" min="0.1" max="1" step="0.05" value="0.5">' +
          '<div class="hint">Below 1 reports fragments.</div></div>' +
        '<div class="field narrow"><label>Largest indel</label>' +
          '<input class="inp mono" id="anBand" type="number" min="0" max="60" value="12"></div>' +
        '<div class="field"><label>&nbsp;</label><button class="btn primary" id="anGo">Search</button></div>' +
        '</div>' +
        '<div id="anOut"></div>' +
        '<hr style="border:0;border-top:1px solid var(--line);margin:14px 0">' +
        '<details><summary>Import a library</summary><div class="note">' +
        '<p style="margin:0 0 8px">The libraries worth importing already exist and are curated; ' +
        'this reads three shapes:</p><ul style="margin:0 0 8px 18px;padding:0">' +
        '<li><strong>ApE</strong> — its default feature file, tab-separated: name, sequence, type. ' +
        'The closest fit and the easiest import.</li>' +
        '<li><strong>FPbase</strong> — a JSON protein list from its open API, for fluorescent proteins.</li>' +
        '<li><strong>A PlasmidBench library</strong> — JSON with the source of every entry kept.</li>' +
        '</ul>' +
        '<p style="margin:0">Each entry keeps where it came from, and that shows in the ' +
        '<em>From</em> column above. Check the licence of a source before passing a library on: ' +
        'pLannotate’s feature set traces through GenoLIB to SnapGene’s curation, and ApE ' +
        'and FPbase have terms of their own.</p></div>' +
        '<div class="row"><div class="field wide"><label>Call this source</label>' +
        '<input class="inp" id="anSrc" placeholder="ApE default features"></div>' +
        '<div class="field wide"><label>Its licence, if you know it</label>' +
        '<input class="inp" id="anLic" placeholder="GPL-3.0, CC-BY, ask the author…"></div>' +
        '<div class="field"><label>&nbsp;</label><button class="btn" id="anPick">Choose a file…</button></div>' +
        '</div><div id="anLibOut">' + libraryState(lib) + '</div>' +
        '<input type="file" id="anFile" accept=".txt,.tsv,.tab,.json,.csv" style="display:none">' +
        '</details>'
    });

    w.addEventListener('click', async ev => {
      const t = /** @type {HTMLElement} */ (ev.target);
      if (!t) return;

      if (t.id === 'anGo') {
        const k = Math.max(8, Math.min(20, +(/** @type {HTMLInputElement} */ ($('#anK')).value) || 12));
        const opts = {
          k: k,
          minIdentity: +(/** @type {HTMLInputElement} */ ($('#anId')).value) || 0.85,
          minCoverage: +(/** @type {HTMLInputElement} */ ($('#anCov')).value) || 0.5,
          band: +(/** @type {HTMLInputElement} */ ($('#anBand')).value) || 12,
          circular: rec.topology === 'circular'
        };
        const host = $('#anOut');
        if (host) host.innerHTML = '<div class="note">Searching…</div>';
        /* let the note paint before a search that may take a moment */
        await new Promise(r => setTimeout(r, 20));
        const lib2 = library();
        const idx = ensureIndex(lib2, k);
        lastResult = B.annotateApprox(rec.seq, idx, opts);
        render();
        return;
      }

      if (t.id === 'anNone') {
        w.querySelectorAll('[data-anno]').forEach(c =>
          { /** @type {HTMLInputElement} */ (c).checked = false; });
        return;
      }

      if (t.id === 'anAdd') {
        const picked = Array.from(w.querySelectorAll('[data-anno]'))
          .filter(c => /** @type {HTMLInputElement} */ (c).checked)
          .map(c => lastResult.hits[+c.getAttribute('data-anno')]);
        if (!picked.length) { PB.toast('Nothing ticked.', 'warn'); return; }
        picked.forEach(h => {
          /* What the match was is kept on the feature, not thrown away once
             it is drawn: a feature called at 91% identity should still say
             so a year later. */
          const note = h.exact ? 'exact match to ' + (h.source || 'the built-in library')
            : (h.identity * 100).toFixed(1) + '% identical over ' +
              (h.coverage * 100).toFixed(0) + '% of ' + (h.source || 'the built-in library') +
              (h.gaps ? ', ' + h.gaps + ' bp indel' : '');
          rec.features.push({
            start: h.start, end: h.end, strand: h.strand, type: h.type,
            label: h.name,
            quals: { note: 'PlasmidBench annotation: ' + note },
            annotated: { identity: h.identity, coverage: h.coverage, exact: h.exact,
                         source: h.source || null, how: h.how }
          });
        });
        rec.features.sort((a, b) => a.start - b.start);
        /* PB.touch records the step so it can be undone and writes the
           history line, and persists; the same path every other editing
           tool uses. */
        PB.touch(rec);
        PB.closeModal();
        PB.renderAll();
        PB.toast('Added ' + picked.length + ' feature' + (picked.length === 1 ? '' : 's') +
          ', each noting how well it matched.', 'ok');
        return;
      }

      if (t.id === 'anPick') { const f = $('#anFile'); if (f) f.click(); return; }
    });

    const fileEl = $('#anFile');
    if (fileEl) fileEl.addEventListener('change', async ev => {
      const f = /** @type {HTMLInputElement} */ (ev.target).files[0];
      if (!f) return;
      const text = await f.text();
      const src = (/** @type {HTMLInputElement} */ ($('#anSrc')).value || '').trim() || f.name;
      const lic = (/** @type {HTMLInputElement} */ ($('#anLic')).value || '').trim() || null;
      const r = B.parseFeatureLibrary(text, { source: src, licence: lic });
      const host = $('#anLibOut');
      if (!r.ok) { if (host) host.innerHTML = '<div class="note bad">' + esc(r.error) + '</div>'; return; }
      PB.lsSet('featureLib', { source: src, licence: lic, entries: r.entries });
      index = null;                                  /* the index is for the old library */
      if (host) host.innerHTML =
        '<div class="note ok"><strong>' + fmt(r.count) + ' features</strong> from ' + esc(src) +
        (lic ? ' (' + esc(lic) + ')' : '') + '. Lengths ' + r.shortest + ' to ' + fmt(r.longest) +
        ' bp, median ' + r.median + '.' +
        (r.duplicates ? ' ' + r.duplicates + ' duplicate sequence' +
          (r.duplicates === 1 ? '' : 's') + ' folded into one entry with aliases.' : '') +
        (r.skipped ? ' ' + r.skipped + ' line' + (r.skipped === 1 ? '' : 's') +
          ' could not be read.' : '') +
        (r.coincidental ? '<br><strong>' + r.coincidental + ' of them are short enough to turn ' +
          'up in a 5 kb plasmid by chance</strong>, so treat those calls with suspicion — ' +
          'the length of each entry is in the table above.' : '') +
        '</div>' + libraryState(library());
      PB.toast(fmt(r.count) + ' features imported from ' + src, 'ok');
    });
  };

  function libraryState(lib) {
    const extra = PB.lsGet('featureLib', null);
    if (!extra || !extra.entries) {
      return '<div class="hint">No library imported. The ' + B.FEATURE_LIB.length +
        ' built-in entries are in use.</div>';
    }
    return '<div class="row"><div class="field wide"><div class="hint">' +
      fmt(extra.entries.length) + ' imported from <strong>' + esc(extra.source || 'a file') +
      '</strong>' + (extra.licence ? ', ' + esc(extra.licence) : '') +
      ', plus the ' + B.FEATURE_LIB.length + ' built in.</div></div>' +
      '<div class="field"><label>&nbsp;</label><button class="btn sm" id="anExport">Export it</button></div>' +
      '<div class="field"><label>&nbsp;</label><button class="btn sm" id="anClear">Remove it</button></div>' +
      '</div>';
  }

  /* export and clear are delegated from the panel, since the block is redrawn */
  document.addEventListener('click', ev => {
    const t = /** @type {HTMLElement} */ (ev.target);
    if (!t || !t.id) return;
    if (t.id === 'anExport') {
      const extra = PB.lsGet('featureLib', null);
      if (!extra) return;
      PB.offerFile((extra.source || 'library').replace(/[^\w.-]+/g, '_') + '.json',
        B.featureLibraryToJSON(extra.entries, extra), 'feature library');
    }
    if (t.id === 'anClear') {
      PB.lsSet('featureLib', null);
      index = null;
      const host = $('#anLibOut');
      if (host) host.innerHTML = '<div class="note">Removed. Back to the ' +
        B.FEATURE_LIB.length + ' built-in entries.</div>';
      PB.toast('Imported library removed.', 'ok');
    }
  });
})(typeof window !== 'undefined' ? window : globalThis);
