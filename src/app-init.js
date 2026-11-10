/* ============================================================================
   PlasmidBench — wiring and start-up
   ============================================================================ */
(function (PB) {
  'use strict';
  const B = window.Bio;
  const $ = PB.$, state = PB.state;

  function wireChrome() {
    // construct name
    const nameEl = $('#recName');
    nameEl.addEventListener('change', () => {
      const rec = PB.active();
      if (!rec) return;
      if (rec.builtin) { PB.toast('Duplicate the demo construct to rename it.', 'warn'); nameEl.value = rec.name; return; }
      const v = nameEl.value.trim();
      if (!v) { nameEl.value = rec.name; return; }
      rec.name = v;
      PB.addHistory(rec, 'Renamed', v);
      PB.touch(rec); PB.renderLibrary(); PB.renderViews();
    });

    $('#btnFind').addEventListener('click', () => PB.searchDialog());
    $('#btnUndo').addEventListener('click', () => PB.undo());
    $('#btnRedo').addEventListener('click', () => PB.redo());
    $('#btnTools').addEventListener('click', () => PB.openTools());
    $('#btnImport').addEventListener('click', () => PB.openImport());
    $('#btnExport').addEventListener('click', () => PB.openExport());
    $('#btnNew').addEventListener('click', () => {
      const nr = PB.makeRecord({ name: 'New construct', seq: '', topology: 'circular' });
      PB.addHistory(nr, 'Created', 'Empty construct');
      PB.addRecord(nr);
      PB.toast('Paste a sequence with Tools → Replace the whole sequence, or use Import.');
    });
    $('#btnDup').addEventListener('click', () => {
      const rec = PB.active();
      if (!rec) return;
      const nr = PB.makeRecord({
        name: rec.name.replace(/\s*\(demo\)\s*/, ' ') .trim() + ' copy',
        seq: rec.seq, topology: rec.topology,
        features: rec.features.map(f => Object.assign({}, f)),
        description: rec.description
      });
      PB.addHistory(nr, 'Duplicated', 'Copy of ' + rec.name);
      PB.addRecord(nr);
    });

    $('#btnEnzNone').addEventListener('click', () => {
      state.enzymes = [];
      PB.lsSet('enzymes', state.enzymes);
      PB.invalidateSites(); PB.renderAll();
    });
    $('#btnEnzCustom').addEventListener('click', () => PB.customEnzymeDialog());

    const search = $('#enzSearch');
    search.addEventListener('input', () => { state.enzQuery = search.value; PB.renderEnzymeRail(); });
    PB.$$('#enzChips .chip').forEach(c => c.addEventListener('click', () => {
      state.enzSet = c.getAttribute('data-eset');
      PB.lsSet('enzSet', state.enzSet);
      PB.renderEnzymeRail();
    }));

    // theme
    const applyTheme = t => {
      if (t === 'system') document.documentElement.removeAttribute('data-theme');
      else document.documentElement.setAttribute('data-theme', t);
    };
    /* so that applying a lab profile can set it on a page already running */
    PB.setTheme = applyTheme;
    const saved = PB.lsGet('theme', 'system');
    applyTheme(saved);
    $('#btnTheme').addEventListener('click', () => {
      const cur = PB.lsGet('theme', 'system');
      const order = ['system', 'light', 'dark'];
      const next = order[(order.indexOf(cur) + 1) % order.length];
      PB.lsSet('theme', next);
      applyTheme(next);
      PB.toast('Theme: ' + next);
      PB.renderViews();
    });

    // slide-over panels on narrow screens
    $('#btnRail').addEventListener('click', () => $('#rail').classList.toggle('open'));
    $('#btnInsp').addEventListener('click', () => $('#inspector').classList.toggle('open'));
    document.addEventListener('click', e => {
      if (window.innerWidth > 760) return;
      const rail = $('#rail');
      if (rail.classList.contains('open') && !rail.contains(e.target) && !$('#btnRail').contains(e.target)) {
        rail.classList.remove('open');
      }
    });

    // file input + drag and drop
    $('#fileInput').addEventListener('change', e => PB.handleFiles(e.target.files));
    const stop = e => { e.preventDefault(); e.stopPropagation(); };
    ['dragenter', 'dragover', 'dragleave', 'drop'].forEach(ev =>
      document.addEventListener(ev, stop, false));
    document.addEventListener('drop', e => {
      if (e.dataTransfer && e.dataTransfer.files && e.dataTransfer.files.length) {
        PB.handleFiles(e.dataTransfer.files);
      } else if (e.dataTransfer) {
        const t = e.dataTransfer.getData('text');
        if (t) PB.importText(t, 'dropped text');
      }
    });

    /* Pasting DNA while bases are selected opens the editor with the pasted
       bases in it, rather than changing anything straight away: a click
       selects a base, and silently replacing it would be a surprise. */
    document.addEventListener('paste', e => {
      const tag = (e.target.tagName || '').toLowerCase();
      if (tag === 'input' || tag === 'textarea' || tag === 'select' || e.target.isContentEditable) return;
      if (!(state.view === 'seq' || state.view === 'map')) return;
      const rec = PB.active();
      if (!rec || !state.selection) return;
      const txt = (e.clipboardData && e.clipboardData.getData('text')) || '';
      const bases = B.cleanSeq(txt.replace(/^>.*$/gm, ''));
      if (!bases.length) return;
      e.preventDefault();
      const L = rec.seq.length, n = PB.normSel(rec, state.selection);
      if (!n) return;
      PB.toolInsertDelete({ kind: n.length === 1 ? 'insert' : 'replace', from: n.start + 1,
        to: (n.end - 1) % L + 1, seq: bases });
    });

    // keyboard shortcuts
    document.addEventListener('keydown', e => {
      const tag = (e.target.tagName || '').toLowerCase();
      if (tag === 'input' || tag === 'textarea' || tag === 'select') return;
      if (e.metaKey || e.ctrlKey) {
        /* inside a text box, the browser's own text undo applies instead */
        const k = (e.key || '').toLowerCase();
        if (k === 'z' && !e.shiftKey) { PB.undo(); e.preventDefault(); return; }
        if ((k === 'z' && e.shiftKey) || k === 'y') { PB.redo(); e.preventDefault(); return; }
        if (e.key === 'c' && state.selection) {
          const rec = PB.active();
          if (rec) { PB.copyText(PB.selectionSeq(rec, state.selection), 'Selection'); e.preventDefault(); }
        }
        return;
      }
      /* Delete or Backspace removes the selected bases, in place; it is one
         undoable step like any other edit */
      if ((e.key === 'Delete' || e.key === 'Backspace') && state.selection &&
          (state.view === 'seq' || state.view === 'map')) {
        const rec = PB.active();
        if (rec) {
          e.preventDefault();
          if (rec.builtin) { PB.toast('The demo construct does not change — duplicate it to edit.', 'warn'); return; }
          const n = PB.normSel(rec, state.selection);
          if (n) PB.editInPlace(n.start, n.end, '');
        }
        return;
      }
      if (state.view === 'map' && ('+=-0'.indexOf(e.key) !== -1)) {
        const which = e.key === '0' ? 'fit' : (e.key === '-' ? 'out' : 'in');
        const btn = document.querySelector('[data-zm="' + which + '"]');
        if (btn) { btn.click(); e.preventDefault(); return; }
      }
      const map = { m: 'map', s: 'seq', f: 'features', e: 'enzymes', o: 'orfs',
                    p: 'primers', g: 'gel', v: 'verify', a: 'align', h: 'history' };
      if (map[e.key]) {
        state.view = map[e.key];
        PB.lsSet('view', state.view);
        PB.renderTabs(); PB.renderViews();
      } else if (e.key === 't') PB.openTools();
      else if (e.key === '/') { PB.searchDialog(); e.preventDefault(); }
      else if (e.key === 'Escape') {
        state.selection = null; state.selectedFeature = null;
        PB.renderHeader(); PB.renderViews(); PB.renderInspector();
      }
    });
  }

  function boot() {
    try {
      wireChrome();
      /* the demo construct opens the tool in a working state, unless it
         has been sent away: someone with their own plasmids loaded does
         not need an example cluttering the library */
      if (!PB.demoHidden()) {
        const demo = PB.buildDemo();
        state.records.push(demo);
        state.activeId = demo.id;
      }
      PB.renderAll();
      PB.initStorage();
    } catch (e) {
      document.body.innerHTML = '<div style="padding:24px;font-family:system-ui">' +
        '<h2 style="margin-bottom:8px">PlasmidBench could not start</h2>' +
        '<p style="color:#666">' + PB.esc(e && e.message ? e.message : String(e)) + '</p></div>';
      if (window.console) console.error(e);
    }
  }

  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', boot);
  else boot();
})(window.PB);
