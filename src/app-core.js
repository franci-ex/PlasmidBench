/* ============================================================================
   PlasmidBench — core: state, persistence, library rail, enzyme rail,
   inspector, helpers.
   ============================================================================ */
window.PB = window.PB || {};
(function (PB) {
  'use strict';
  const B = window.Bio;

  /* ---------------- tiny DOM helpers ---------------- */
  const $ = sel => document.querySelector(sel);
  const $$ = sel => Array.prototype.slice.call(document.querySelectorAll(sel));
  function esc(s) {
    return String(s == null ? '' : s)
      .replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;')
      .replace(/"/g, '&quot;').replace(/'/g, '&#39;');
  }
  function fmt(n) { return Number(n).toLocaleString('en-US'); }
  function f1(n) { return (Math.round(n * 10) / 10).toFixed(1); }
  function f2(n) { return (Math.round(n * 100) / 100).toFixed(2); }
  function uid() {
    return 'c' + Date.now().toString(36) + Math.random().toString(36).slice(2, 8);
  }
  let toastTimer = null;
  function toast(msg, kind) {
    const t = $('#toast');
    t.textContent = msg;
    t.className = 'toast' + (kind ? ' ' + kind : '');
    t.hidden = false;
    clearTimeout(toastTimer);
    toastTimer = setTimeout(() => { t.hidden = true; },
      Math.min(12000, Math.max(kind === 'bad' ? 5200 : (kind === 'warn' ? 4200 : 2700), String(msg).length * 55)));
  }
  async function copyText(text, what) {
    try {
      await navigator.clipboard.writeText(text);
      toast((what || 'Copied') + ' — ' + fmt(text.length) + ' characters on the clipboard', 'ok');
      return true;
    } catch (e) {
      PB.showTextFallback(what || 'Copy this text', text);
      return false;
    }
  }
  /* Settings read synchronously from the store's in-memory cache and are
     written through to IndexedDB behind the scenes; see store.js for why it
     is arranged that way. The two names are kept because around a hundred
     call sites use them, and because they still mean the same thing: a
     small value that belongs to this browser. If store.js has not loaded —
     a stripped copy of the page, a test harness — they fall back to
     localStorage so nothing breaks outright. */
  function lsGet(k, dflt) {
    if (PB.store) return PB.store.get(k, dflt);
    try {
      const v = localStorage.getItem('pb.' + k);
      return v == null ? dflt : JSON.parse(v);
    } catch (e) { return dflt; }
  }
  function lsSet(k, v) {
    if (PB.store) { PB.store.set(k, v); return; }
    try { localStorage.setItem('pb.' + k, JSON.stringify(v)); } catch (e) { /* ignore */ }
  }

  PB.$ = $; PB.$$ = $$; PB.esc = esc; PB.fmt = fmt; PB.f1 = f1; PB.f2 = f2;
  PB.uid = uid; PB.toast = toast; PB.copyText = copyText;
  PB.lsGet = lsGet; PB.lsSet = lsSet;

  /* ---------------- feature colouring ---------------- */
  const TYPE_COLOR = {
    promoter: 'f-promoter', enhancer: 'f-promoter', TATA_signal: 'f-promoter',
    CDS: 'f-cds', gene: 'f-cds', ORF: 'f-cds', mat_peptide: 'f-cds',
    terminator: 'f-terminator', polyA_signal: 'f-terminator', polyA_site: 'f-terminator',
    rep_origin: 'f-origin', origin: 'f-origin', oriT: 'f-origin',
    primer_bind: 'f-primer',
    RBS: 'f-reg', protein_bind: 'f-reg', regulatory: 'f-reg', operator: 'f-reg',
    misc_binding: 'f-reg', LTR: 'f-reg', intron: 'f-reg', exon: 'f-cds',
    recomb: 'f-recomb', misc_recomb: 'f-recomb', protein_tag: 'f-tag',
    sig_peptide: 'f-tag', misc_feature: 'f-misc', source: 'f-misc'
  };
  const MARKER_RE = /\b(ampr|bla|kanr|neor?|puror?|hygror?|bsd|cmr|tetr|smr|zeor|blastic|amp\b|kan\b|chlor|spec|gent|resist|ccdb|sacb)/i;
  const TAG_RE = /\b(his6|6xhis|his-tag|flag|ha tag|myc|gst|mbp|strep|sumo|snap|halo|v5|avi|tag)\b/i;
  function featColorVar(f) {
    const lbl = (f.label || '') + ' ' + (f.name || '');
    if (MARKER_RE.test(lbl)) return 'f-marker';
    if (TAG_RE.test(lbl)) return 'f-tag';
    if (/gfp|cherry|rfp|yfp|cfp|venus|tdtomato|luciferase|luc\b/i.test(lbl)) return 'f-cds';
    return TYPE_COLOR[f.type] || 'f-misc';
  }
  function featColor(f) {
    // a colour carried in from SnapGene wins, so imported maps keep their look
    if (f && f.color && /^#[0-9a-f]{3,8}$/i.test(f.color)) return f.color;
    return getComputedStyle(document.documentElement)
      .getPropertyValue('--' + featColorVar(f)).trim() || '#5c6a6e';
  }
  PB.featColorVar = featColorVar; PB.featColor = featColor;

  /* ---------------- state ---------------- */
  const state = {
    records: [],
    activeId: null,
    blast: null,                      // {query, format, hits[]} pasted back from NCBI
    view: lsGet('view', 'map'),
    selection: null,                  // {start, end} 0-based, end exclusive
    selectedFeature: null,            // index into active record features
    /* Several features can be selected at once, with ctrl or cmd held.
       `selectedFeature` stays the last one picked, so everything that
       reads a single selection — the inspector, Edit feature, the codon
       map — keeps working unchanged; `selectedFeatures` is the whole set
       and only the views that paint ranges look at it. Anything that
       clears `selectedFeature` clears the set with it, because
       featureSelection() treats a null primary as nothing selected. */
    selectedFeatures: [],
    enzymes: lsGet('enzymes', ['EcoRI', 'BamHI', 'HindIII', 'XhoI', 'NotI', 'NdeI', 'XbaI', 'SalI', 'KpnI', 'SacI', 'PstI', 'SphI', 'BsaI', 'EcoRV', 'NcoI', 'SpeI', 'AgeI', 'BglII']),
    enzSet: lsGet('enzSet', 'common'),
    enzQuery: '',
    customEnzymes: lsGet('customEnzymes', []),
    primers: lsGet('primers', []),
    mapMode: lsGet('mapMode', 'circular'),
    showPrimers: lsGet('showPrimers', true),
    /* which stored primers are drawn on the map; null means all of them */
    primersShown: lsGet('primersShown', null),
    /* the primer-binding sites annotated in the file itself, which on a
       commercial vector are a dozen labels nobody is designing against */
    showFilePrimers: lsGet('showFilePrimers', true),
    seqOpts: lsGet('seqOpts', { bpl: 60, bottom: true, aa: true, colorBases: false, revTrans: false }),
    tmCond: lsGet('tmCond', { oligo: 0.25, Na: 50, Mg: 1.5, dNTP: 0.2 }),
    codeTable: lsGet('codeTable', 1),
    gelLanes: null,
    db: null, user: null, downloads: null,
    canWrite: null, dbReady: false
  };
  PB.state = state;

  function conditions() {
    const c = state.tmCond;
    return { oligo: c.oligo * 1e-6, Na: c.Na / 1000, Mg: c.Mg, dNTP: c.dNTP };
  }
  PB.conditions = conditions;

  /* Re-read the settings that were copied into state at start-up. Applying a
     lab profile changes them underneath a page that is already running, and
     a profile that only takes effect after a reload is a profile people do
     not trust. Everything not listed here is read at the point of use and
     needs nothing done. */
  PB.reloadSettings = function () {
    state.enzSet = lsGet('enzSet', 'common');
    state.seqOpts = lsGet('seqOpts', { bpl: 60, bottom: true, aa: true, colorBases: false, revTrans: false });
    state.tmCond = lsGet('tmCond', { oligo: 0.25, Na: 50, Mg: 1.5, dNTP: 0.2 });
    state.codeTable = lsGet('codeTable', 1);
    state.mapMode = lsGet('mapMode', state.mapMode);
    const th = lsGet('theme', null);
    if (th && PB.setTheme) PB.setTheme(th);
    if (PB.renderEnzymeRail) PB.renderEnzymeRail();
  };

  /* A selection on a circle can run across base 1, and it arrives in two
     shapes: a search hit there ends past the last base (end > length),
     and a feature that spans base 1 ends before it starts (end <= start).
     Everything that reads a selection goes through these two, so both
     shapes give the same bases. */
  function normSel(rec, sel) {
    if (!rec || !sel) return null;
    const L = rec.seq.length;
    let a = sel.start, b = sel.end;
    if (!(a >= 0 && a < L)) return null;
    if (b <= a) { if (rec.topology !== 'circular') return null; b += L; }
    if (b > L && rec.topology !== 'circular') b = L;
    if (b - a > L) b = a + L;
    return { start: a, end: b, wraps: b > L, length: b - a };
  }
  function selectionSeq(rec, sel) {
    const n = normSel(rec, sel);
    if (!n) return '';
    return n.wraps ? rec.seq.slice(n.start) + rec.seq.slice(0, n.end - rec.seq.length) : rec.seq.slice(n.start, n.end);
  }
  /* "1,151–1,200 and 1–40" for a selection across base 1 */
  function selRangeText(rec, sel) {
    const n = normSel(rec, sel);
    if (!n) return '';
    const L = rec.seq.length;
    return n.wraps ? fmt(n.start + 1) + '–' + fmt(L) + ' and 1–' + fmt(n.end - L) + ' (across base 1)'
                   : fmt(n.start + 1) + '–' + fmt(n.end);
  }
  /* Nearest-neighbour Tm is a property of a short duplex. Below 8 nt the
     model gives numbers such as −68 °C for GGG; above 60 nt a duplex does
     not melt as a two-state unit and the figure is not a Tm in any useful
     sense. Outside that window the value is not shown at all. */
  const TM_MIN = 8, TM_MAX = 60;
  function tmShown(s) { return s.length >= TM_MIN && s.length <= TM_MAX; }
  /* ---------------- selecting more than one feature ----------------
     The set is kept normalised against the record on every read, so a
     feature deleted or an edit that shortened the list cannot leave a
     stale index behind, and so that the hundred places which clear
     `selectedFeature` do not each have to know about the set. */
  function featureSelection(rec) {
    const one = state.selectedFeature;
    if (one == null || !rec || !rec.features || !rec.features[one]) return [];
    const arr = Array.isArray(state.selectedFeatures) ? state.selectedFeatures : [];
    if (arr.indexOf(one) === -1) return [one];
    const seen = {}, out = [];
    arr.forEach(i => { if (rec.features[i] && !seen[i]) { seen[i] = 1; out.push(i); } });
    return out.length ? out : [one];
  }
  /* click: this one alone, or nothing if it was the only one already.
     ctrl/cmd-click: add it, or take it out if it is already in. */
  function toggleFeatureSelection(rec, i, additive) {
    const cur = featureSelection(rec);
    const has = cur.indexOf(i) !== -1;
    let next;
    if (additive) next = has ? cur.filter(x => x !== i) : cur.concat([i]);
    else next = (cur.length === 1 && has) ? [] : [i];
    state.selectedFeatures = next;
    state.selectedFeature = next.length ? next[next.length - 1] : null;
    const f = next.length ? rec.features[next[next.length - 1]] : null;
    state.selection = f ? { start: f.start, end: f.end } : null;
    return next;
  }
  /* every selected feature as plain top-strand pieces, joins and
     origin-wraps opened out, for the views that paint ranges */
  function featureSelectionPieces(rec) {
    const idx = featureSelection(rec);
    if (idx.length < 2) return null;         /* one or none: the ordinary path */
    const L = rec.seq.length, out = [];
    idx.forEach(i => {
      B.featurePieces(rec.features[i], L).forEach(p => out.push({ start: p.start, end: p.end }));
    });
    return out;
  }
  /* cut a plain string before it is escaped, so a slice can never land in
     the middle of an &amp; and leave half an entity on the page */
  function trim(s, n) { s = String(s || ''); return s.length > n ? s.slice(0, n - 1) + '…' : s; }
  PB.trimText = trim;
  PB.featureSelection = featureSelection;
  PB.toggleFeatureSelection = toggleFeatureSelection;
  PB.featureSelectionPieces = featureSelectionPieces;
  PB.isFeatureSelected = (rec, i) => featureSelection(rec).indexOf(i) !== -1;

  PB.selRangeText = selRangeText; PB.tmShown = tmShown; PB.TM_RANGE = [TM_MIN, TM_MAX];
  PB.normSel = normSel;
  PB.selectionSeq = selectionSeq;
  PB.selectionText = () => selectionSeq(state.records.find(x => x.id === state.activeId), state.selection);

  /* The state a construct was in when it was first looked at is the base
     the first undo returns to — every change is made to the active
     construct, so it is always captured before anything can alter it. */
  function active() {
    const r = state.records.find(x => x.id === state.activeId) || null;
    if (r && PB.undoBaseline) PB.undoBaseline(r);
    return r;
  }
  PB.active = active;

  function allEnzymes() { return B.ALL_ENZYMES.concat(state.customEnzymes); }
  PB.allEnzymes = allEnzymes;
  function enzByName(n) {
    return allEnzymes().find(e => e.name.toLowerCase() === String(n).toLowerCase());
  }
  PB.enzByName = enzByName;

  /* ---------------- records ---------------- */
  /* GenBank and EMBL require a `source` feature spanning the whole molecule,
     carrying /organism and /mol_type. It is file bookkeeping, not biology:
     drawn, it is a solid band right round the plasmid that buries every real
     feature under it, and it pads the feature count. SnapGene and Benchling
     both hide it, so it is dropped on the way in. The GenBank exporter writes
     its own source line, so a round trip still produces a valid file. */
  function isBoilerplate(f) {
    return !f || String(f.type).toLowerCase() === 'source';
  }
  PB.isBoilerplateFeature = isBoilerplate;

  function makeRecord(o) {
    return {
      id: o.id || uid(),
      name: o.name || 'Untitled',
      seq: B.cleanSeq(o.seq || ''),
      topology: o.topology === 'circular' ? 'circular' : 'linear',
      features: (o.features || []).filter(f => !isBoilerplate(f)).map(normFeature),
      description: o.description || '',
      builtin: !!o.builtin,
      history: o.history || [],
      methylation: o.methylation || null,
      createdAt: o.createdAt || Date.now(),
      updatedAt: o.updatedAt || Date.now()
    };
  }
  function normFeature(f) {
    const int = v => Math.max(0, Math.round(+v) || 0);
    const g = {
      type: f.type || 'misc_feature',
      start: int(f.start),
      end: int(f.end),
      strand: f.strand === -1 ? -1 : 1,
      label: f.label || f.name || f.type || 'feature',
      quals: f.quals || {},
      color: f.color || null,
      note: f.note || ''
    };
    if (Array.isArray(f.parts) && f.parts.length > 1) {
      g.parts = f.parts.map(p => {
        const q = { start: int(p.start), end: int(p.end) };
        if (f.mixed) q.strand = p.strand === -1 ? -1 : 1;
        return q;
      });
      if (f.mixed) g.mixed = true;
    }
    if (f.between) { g.between = true; g.end = g.start; }
    return g;
  }
  PB.makeRecord = makeRecord; PB.normFeature = normFeature;

  function addHistory(rec, op, detail) {
    rec.history = rec.history || [];
    rec.history.push({ op, detail: detail || '', at: Date.now() });
    if (rec.history.length > 200) rec.history = rec.history.slice(-200);
  }
  PB.addHistory = addHistory;

  function addRecord(rec, opts) {
    opts = opts || {};
    /* every record needs an id and normalised features; anything that did
       not come through makeRecord is passed through it here */
    if (!rec.id || !Array.isArray(rec.history)) {
      const keep = rec;
      rec = makeRecord(rec);
      if (keep.builtin) rec.builtin = true;
    }
    state.records.push(rec);
    if (opts.activate !== false) { state.activeId = rec.id; state.selection = null; state.selectedFeature = null; }
    persist(rec);
    PB.renderAll();
    return rec;
  }
  PB.addRecord = addRecord;

  /* ---------------- undo / redo ----------------
     Every change made to a construct in place is saved through touch(), so
     that is where the state before the change is captured. Each construct
     keeps its own history of steps, the way each document does in a
     desktop editor; undo never reaches into a construct you are not looking
     at. Steps are kept for this session only. The label of each step comes
     from the history entry the change wrote, so "Undo" can say what it is
     about to undo. Long sequences keep fewer steps, so a genome-sized
     record cannot exhaust memory. */
  const undoStacks = {};
  const UNDO_MAX = 100, UNDO_CHARS = 8e6;
  function snapOf(rec) {
    return {
      name: rec.name, seq: rec.seq, topology: rec.topology,
      description: rec.description || '',
      features: JSON.parse(JSON.stringify(rec.features || [])),
      methylation: rec.methylation ? JSON.parse(JSON.stringify(rec.methylation)) : null
    };
  }
  function sameSnap(a, b) {
    return a.seq === b.seq && a.name === b.name && a.topology === b.topology &&
      a.description === b.description &&
      JSON.stringify(a.features) === JSON.stringify(b.features) &&
      JSON.stringify(a.methylation) === JSON.stringify(b.methylation);
  }
  function lastHist(rec) { return rec.history && rec.history.length ? rec.history[rec.history.length - 1] : null; }
  function undoBaseline(rec) {
    if (!rec || rec.builtin || undoStacks[rec.id]) return;
    undoStacks[rec.id] = { base: snapOf(rec), hist: lastHist(rec), past: [], future: [] };
  }
  function stepCap(rec) {
    return Math.max(1, Math.min(UNDO_MAX, Math.floor(UNDO_CHARS / Math.max(1, rec.seq.length))));
  }
  function touch(rec, opts) {
    const u = undoStacks[rec.id];
    if (u && !(opts && opts.fromUndo)) {
      const now = snapOf(rec);
      if (!sameSnap(u.base, now)) {
        const h = lastHist(rec);
        const label = h && h !== u.hist ? h.op + (h.detail ? ' — ' + h.detail : '') : 'Edit';
        u.past.push({ snap: u.base, label: label });
        const cap = stepCap(rec);
        while (u.past.length > cap) u.past.shift();
        u.future = [];
        u.base = now; u.hist = h;
      }
    } else if (!u) undoBaseline(rec);
    rec.updatedAt = Date.now();
    persist(rec);
    if (PB.renderUndo) PB.renderUndo();
  }
  PB.touch = touch;

  function applySnap(rec, snap) {
    rec.name = snap.name; rec.seq = snap.seq; rec.topology = snap.topology;
    rec.description = snap.description;
    rec.features = JSON.parse(JSON.stringify(snap.features));
    rec.methylation = snap.methylation ? JSON.parse(JSON.stringify(snap.methylation)) : null;
  }
  function stepAcross(dir) {
    const rec = active();
    const u = rec && undoStacks[rec.id];
    const from = u && (dir < 0 ? u.past : u.future);
    if (!from || !from.length) {
      toast(dir < 0 ? 'Nothing to undo for this construct.' : 'Nothing to redo.', 'warn');
      return null;
    }
    const step = from.pop();
    (dir < 0 ? u.future : u.past).push({ snap: u.base, label: step.label });
    applySnap(rec, step.snap);
    addHistory(rec, dir < 0 ? 'Undone' : 'Redone', step.label);
    u.base = snapOf(rec); u.hist = lastHist(rec);
    /* a selection can point past the end of a sequence that just got shorter */
    if (state.selection && state.selection.start >= rec.seq.length) state.selection = null;
    if (state.selectedFeature != null && !rec.features[state.selectedFeature]) state.selectedFeature = null;
    touch(rec, { fromUndo: true });
    if (PB.invalidateSites) PB.invalidateSites();
    PB.renderAll();
    toast((dir < 0 ? 'Undid: ' : 'Redid: ') + step.label, 'ok');
    return step.label;
  }
  PB.undo = () => stepAcross(-1);
  PB.redo = () => stepAcross(1);
  PB.undoState = function () {
    const rec = active();
    const u = rec && undoStacks[rec.id];
    return {
      canUndo: !!(u && u.past.length), canRedo: !!(u && u.future.length),
      undoLabel: u && u.past.length ? u.past[u.past.length - 1].label : '',
      redoLabel: u && u.future.length ? u.future[u.future.length - 1].label : '',
      steps: u ? u.past.length : 0
    };
  };
  PB.forgetUndo = id => { delete undoStacks[id]; };
  PB.undoBaseline = undoBaseline;

  /* ---------------- persistence ----------------
     The library lives in IndexedDB, in this browser, through PB.store. That
     is the system of record: it works with no network, no account and no
     host, which is what a page served from a lab machine needs.

     When the page is opened somewhere that offers a host-side store — the
     claude.ai artifact runtime is the one that does — that store is still
     read at start-up and still written to, as a mirror. Not because a tool
     like this should hold authoritative records off the machine (it should
     not: attribution and contemporaneity belong to the lab's record system,
     which is eLabFTW here), but because constructs saved that way in the
     past would otherwise be stranded, and because the person who has been
     opening this on two devices should not silently lose the second one.
     Setting `pb.mirrorHost` to false turns the mirror off; nothing else
     depends on it, and the page is complete without it. */
  const LOCAL_KEY = 'library.v1';

  function serialize(rec) {
    return {
      name: rec.name, seq: rec.seq, topology: rec.topology,
      features: rec.features, description: rec.description,
      history: rec.history, createdAt: rec.createdAt, updatedAt: rec.updatedAt
    };
  }

  /* Every construct is a row of its own in IndexedDB, so a 300 kb plasmid
     no longer has to fit in a shared string budget with everything else.
     The old single-key array is still written while it fits, so that an
     older copy of the page opened on the same browser still finds the
     library; once it does not fit, IndexedDB alone holds it, which is the
     whole reason for the change. */
  function saveLocal() {
    const own = state.records.filter(r => !r.builtin);
    if (PB.store && PB.store.isAvailable()) {
      own.forEach(r => PB.store.putRecord(Object.assign({ id: r.id }, serialize(r))));
    }
    lsSet(LOCAL_KEY, own.map(r => Object.assign({ id: r.id }, serialize(r))));
  }

  let writeQueue = Promise.resolve();
  function persist(rec) {
    if (rec && rec.builtin) return;
    saveLocal();                          /* IndexedDB: the system of record */
    if (rec && PB.store && PB.store.isAvailable()) PB.store.putRecord(Object.assign({ id: rec.id }, serialize(rec)));
    if (!state.db || state.canWrite === false || !rec) return;
    if (lsGet('mirrorHost', true) === false) return;
    // one write at a time, coalesced
    writeQueue = writeQueue.then(async () => {
      try {
        const body = serialize(rec);
        /* The host store takes 200,000 bases at most. IndexedDB has already
           saved this one, so there is nothing to lose here and nothing to
           warn about beyond saying the mirror did not take it. */
        if (body.seq.length > 200000) {
          if (!(PB.store && PB.store.isAvailable()))
            toast('Sequence is too large for the account store; kept in this browser only.', 'warn');
          return;
        }
        await state.db.collection('plasmids').doc(rec.id).set(body);
      } catch (e) {
        if (e && e.code === 'invalid_argument') {
          state.canWrite = false;
          renderStoreState();
        }
      }
    }).catch(() => {});
    return writeQueue;
  }
  PB.persist = persist;

  function removeRecord(id) {
    const rec = state.records.find(r => r.id === id);
    if (!rec) return;
    /* The demo is rebuilt on every load, so removing it has to be
       remembered or it walks straight back in. It is never stored, so
       there is nothing to delete - only a note not to build it again. */
    if (rec.builtin) {
      lsSet('demoHidden', true);
      state.records = state.records.filter(r => r.id !== id);
      if (state.activeId === id) {
        state.activeId = state.records.length ? state.records[0].id : null;
        state.selection = null; state.selectedFeature = null;
      }
      PB.renderAll();
      return;
    }
    state.records = state.records.filter(r => r.id !== id);
    PB.forgetUndo(id);
    if (state.activeId === id) {
      state.activeId = state.records.length ? state.records[0].id : null;
      state.selection = null; state.selectedFeature = null;
    }
    saveLocal();
    if (PB.store && PB.store.isAvailable()) PB.store.deleteRecord(id);
    if (state.db && state.canWrite !== false && lsGet('mirrorHost', true) !== false) {
      writeQueue = writeQueue.then(() =>
        state.db.collection('plasmids').doc(id).delete().catch(() => {})).catch(() => {});
    }
    PB.renderAll();
  }
  PB.removeRecord = removeRecord;

  PB.demoHidden = () => !!lsGet('demoHidden', false);
  PB.restoreDemo = function () {
    if (state.records.some(r => r.builtin)) { toast('The demo construct is already here'); return; }
    lsSet('demoHidden', false);
    const demo = PB.buildDemo();
    state.records.push(demo);
    state.activeId = demo.id;
    PB.renderAll();
    toast('Demo construct restored', 'ok');
  };

  /* What this says has to be true, because someone decides whether to close
     the tab on it. "Saved to your account" was true only inside the
     claude.ai artifact and was shown everywhere. */
  function renderStoreState() {
    const el = $('#storeState');
    if (!el) return;
    const rep = state.storeReport;
    const idb = PB.store && PB.store.isAvailable();
    let txt;
    if (idb) {
      txt = 'Saved in this browser.';
      if (rep && rep.migratedRecords) txt += ' ' + rep.migratedRecords +
        ' construct' + (rep.migratedRecords === 1 ? '' : 's') + ' moved out of the old store.';
      if (state.db && state.canWrite !== false && lsGet('mirrorHost', true) !== false)
        txt += ' Also mirrored to your account, so it reaches your other devices.';
    } else if (rep && rep.available === false) {
      /* the one case worth a warning: a private window or blocked site data
         means what is "saved" goes when the tab does */
      txt = 'This browser will not keep data (' + (rep.why || 'storage unavailable') +
            '). Export anything you want to keep before closing.';
      el.className = 'warn';
    } else if (!state.dbReady) {
      txt = 'Opening local storage\u2026';
    } else {
      txt = 'Saved in this browser.';
    }
    el.textContent = txt;
  }
  PB.renderStoreState = renderStoreState;

  async function initStorage() {
    /* 1 — IndexedDB. Opening it also seeds the settings cache and brings
       across anything an earlier version of the page left in localStorage,
       so this has to happen before any setting is read for real. */
    if (PB.store) {
      try {
        state.storeReport = await PB.store.init();
      } catch (e) {
        state.storeReport = { available: false, why: String(e && e.message || e) };
      }
      const rows = await PB.store.allRecords();
      for (const r of rows) {
        if (!r || !r.seq) continue;
        const existing = state.records.find(x => x.id === r.id);
        if (existing) {
          if ((r.updatedAt || 0) > (existing.updatedAt || 0)) {
            Object.assign(existing, makeRecord(r));
            PB.forgetUndo(r.id);
          }
        } else state.records.push(makeRecord(r));
      }
    }

    /* 2 — anything still only under the old single-key array */
    const local = lsGet(LOCAL_KEY, []);
    for (const r of local) {
      if (!state.records.some(x => x.id === r.id)) state.records.push(makeRecord(r));
    }
    if (!state.activeId && state.records.length) state.activeId = state.records[0].id;
    renderStoreState();
    PB.renderAll();

    /* 3 — the host store, where there is one, as a mirror. See the note on
       persistence above: this is here so that constructs saved through it
       are not stranded, not because records belong off the machine. */
    if (lsGet('mirrorHost', true) === false ||
        !window.claude || typeof window.claude.use !== 'function') {
      state.dbReady = true; renderStoreState(); return;
    }
    try { state.downloads = await window.claude.use('downloads'); } catch (e) { state.downloads = null; }
    try { state.user = await window.claude.use('user'); } catch (e) { state.user = null; }
    if (state.user) {
      try { state.canWrite = await state.user.can('data.write'); } catch (e) { state.canWrite = null; }
    }
    try { state.db = await window.claude.use('db'); } catch (e) { state.db = null; }
    state.dbReady = true;
    renderStoreState();
    PB.renderExportAvailability && PB.renderExportAvailability();
    if (!state.db) return;
    try {
      const snap = await state.db.collection('plasmids').limit(500).get();
      let added = 0;
      snap.docs.forEach(d => {
        const data = d.data();
        if (!data || !data.seq) return;
        const existing = state.records.find(r => r.id === d.id);
        if (existing) {
          // remote wins when it is newer
          if ((data.updatedAt || 0) > (existing.updatedAt || 0)) {
            Object.assign(existing, makeRecord(Object.assign({ id: d.id }, data)));
            /* the newer saved copy is now the starting point; undo must not
               step back to the stale local one it replaced */
            PB.forgetUndo(d.id);
          }
        } else {
          state.records.push(makeRecord(Object.assign({ id: d.id }, data)));
          added++;
        }
      });
      // push up anything local the store has not seen
      const remoteIds = new Set(snap.docs.map(d => d.id));
      for (const r of state.records) {
        if (!r.builtin && !remoteIds.has(r.id)) persist(r);
      }
      if (!state.activeId && state.records.length) state.activeId = state.records[0].id;
      saveLocal();
      PB.renderAll();
      if (added) toast(added + ' construct' + (added === 1 ? '' : 's') + ' loaded from your account', 'ok');
    } catch (e) { /* store unavailable; local copy already rendered */ }
  }
  PB.initStorage = initStorage;

  /* ---------------- downloads / clipboard output ----------------
     A published page cannot write a file itself; the viewer's host does it,
     after asking, and only for file types on its allowlist. Which types
     those are is decided by the host and can differ between a browser tab,
     the desktop app and the phone, so the same button that writes a file on
     one machine is refused on another.

     What used to happen when it was refused was that the text went to the
     clipboard and a cheerful notice said so. For a 72,000-character SVG
     that is not a save — it is a dead end dressed up as a success. So a
     refusal now says plainly what was refused and offers the ways out that
     do work on that machine, as things to click rather than as a chain of
     silent retries: the same picture in a format the host will take, the
     same bytes under a .txt name to rename afterwards, and the text itself
     to copy. Nothing is tried without being asked for. */
  /* The save capability is resolved once at start-up, and `use()` answers
     null not only when it is unavailable but also when nothing answered
     within ten seconds — which a busy first paint can cause. Caching that
     null would disable saving for the rest of the page's life over a
     moment's timing, so it is asked for again at the moment of saving. */
  async function ensureDownloads() {
    if (state.downloads) return state.downloads;
    if (!window.claude || typeof window.claude.use !== 'function') return null;
    try { state.downloads = await window.claude.use('downloads'); }
    catch (e) { state.downloads = null; }
    if (state.downloads && PB.renderExportAvailability) PB.renderExportAvailability();
    return state.downloads;
  }
  PB.ensureDownloads = ensureDownloads;

  /* The ordinary way a web page hands over a file: a Blob and an anchor
     with a download attribute. It needs no host and no permission, and it
     is what happens when this page is opened from a lab machine or served
     from a static web server — which is the deployment the project is for.

     It is not used *instead* of the host capability where there is one: an
     embedded view may run without allow-downloads, where the click is
     dropped silently and the person is told a file was saved that was not.
     So the host is asked first where it exists, and this is the path
     everywhere else. */
  const SAVE_MIME = {
    svg: 'image/svg+xml', json: 'application/json', xml: 'application/xml',
    html: 'text/html', csv: 'text/csv', tsv: 'text/tab-separated-values',
    gb: 'text/plain', gbk: 'text/plain', fa: 'text/plain', fasta: 'text/plain',
    embl: 'text/plain', gff3: 'text/plain', md: 'text/markdown', txt: 'text/plain',
    xlsx: 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet',
    png: 'image/png', pdf: 'application/pdf'
  };
  const standalone = () => typeof window === 'undefined' || !window.claude;
  PB.standalone = standalone;
  function nativeSave(filename, data) {
    try {
      if (typeof document === 'undefined' || !document.body) return false;
      const a = document.createElement('a');
      if (!('download' in a)) return false;
      const ext = (String(filename).split('.').pop() || '').toLowerCase();
      const blob = (data instanceof Blob) ? data
        : new Blob([data], { type: SAVE_MIME[ext] || 'application/octet-stream' });
      const url = URL.createObjectURL(blob);
      a.href = url;
      a.download = String(filename).split('/').pop();
      a.style.display = 'none';
      document.body.appendChild(a);
      a.click();
      /* the object URL has to outlive the click, not the statement */
      setTimeout(() => { try { URL.revokeObjectURL(url); a.remove(); } catch (e) { /* gone */ } }, 30000);
      return true;
    } catch (e) { return false; }
  }
  PB.nativeSave = nativeSave;

  async function offerFile(filename, text, label, opts) {
    opts = opts || {};
    await ensureDownloads();
    if (!state.downloads) {
      /* Standalone means no host at all — no window.claude — and then the
         page is in an ordinary tab where the browser's own download runs
         and can be reported as done. Inside an embedded view the click can
         be dropped with no error, so there the dialog below is the honest
         answer even though the same code would "work". */
      if (standalone() && nativeSave(filename, text)) { toast('Saved ' + filename, 'ok'); return true; }
      return saveBlocked(filename, text, label, opts, null);
    }
    try {
      await state.downloads.save({ filename, data: text });
      toast('Saved ' + filename, 'ok');
      return true;
    } catch (e) {
      const code = (e && e.code) || 'unavailable';
      if (code === 'declined') { toast('Save cancelled'); return false; }
      if (code === 'too_large') { toast('That file is too large to save here.', 'bad'); return false; }
      if (code === 'rate_limited') { toast('A save is already waiting for an answer — finish that one first.', 'warn'); return false; }
      return saveBlocked(filename, text, label, opts, code);
    }
  }
  PB.offerFile = offerFile;

  /* What to say, and what to offer, when the host will not write this file.
     `code` is null when there is no save capability at all. */
  function saveBlocked(filename, text, label, opts, code) {
    const ext = (String(filename).split('.').pop() || '').toLowerCase();
    const txtName = /\.txt$/i.test(filename) ? filename : filename + '.txt';
    const why = code === null
      ? 'Saving files is not available in this view at all.'
      : (code === 'rejected_extension'
          ? 'This view will not save <b>.' + esc(ext) + '</b> files — its list of allowed file types does not include that one.'
          : (code === 'extension_not_enabled'
              ? 'This view has <b>.' + esc(ext) + '</b> files switched off. That is not the usual state; reloading the page sometimes brings it back.'
              : 'This view could not save <b>.' + esc(ext) + '</b> files just now.'));
    const acts = [];
    /* the same content in a form this host will take, offered, not assumed */
    if (opts.alt) acts.push({ label: opts.alt.label, primary: true, onClick: () => { closeModal(); opts.alt.run(); return false; } });
    if (state.downloads && ext !== 'txt') acts.push({ label: 'Save as ' + txtName.split('/').pop(), onClick: () => {
      closeModal();
      state.downloads.save({ filename: txtName, data: text })
        .then(() => toast('Saved ' + txtName + ' — rename it to .' + ext + ' and it opens as before', 'ok'))
        .catch(err => toast('That did not save either' + (err && err.message ? ': ' + err.message : '.'), 'bad'));
      return false;
    } });
    acts.push({ label: 'Copy the text', onClick: () => { closeModal(); copyText(text, label || filename); return false; } });
    acts.push({ label: 'Close' });
    /* The browser's own download, offered as a link the person clicks
       themselves. A synthetic click is ignored in an embedded view; a real
       one on a real anchor is a user gesture and is often allowed where
       the synthetic one was not. It costs nothing to offer and it is the
       thing most likely to actually produce the file. */
    let url = null;
    try { url = URL.createObjectURL(new Blob([text], { type: SAVE_MIME[ext] || 'text/plain' })); }
    catch (e) { url = null; }
    modal({
      title: label || filename,
      body: '<div class="note">' + why +
        (ext === 'svg' ? ' The picture itself is fine — only writing it to disk under that name is refused.' : '') +
        '</div>' +
        (url ? '<div class="row" style="gap:6px"><a class="btn primary" href="' + url + '" download="' + esc(filename) +
          '">Download ' + esc(filename) + '</a>' +
          '<a class="btn" href="' + url + '" target="_blank" rel="noopener">Open it in a new tab</a></div>' +
          '<div class="hint">The download link is the browser’s own, not this page’s, so it works in places the button above does not. ' +
          'If the new tab opens the picture, saving it from there works too.</div>' : '') +
        '<div class="hint">' + fmt(text.length) + ' characters. Anything saved as <b>.txt</b> holds exactly the same bytes: ' +
        'renaming the file to <b>.' + esc(ext) + '</b> afterwards gives the original file, byte for byte.</div>' +
        '<div class="row"><div class="field wide"><label>Or take it from here</label>' +
        '<textarea class="inp" rows="8" readonly spellcheck="false" id="sbText"></textarea></div></div>',
      actions: acts,
      onOpen: w => {
        const t = w.querySelector('#sbText');
        /* a megabyte of SVG in a textarea locks the page up; the file is
           still whole on the clipboard and in the save */
        t.value = text.length > 400000 ? text.slice(0, 400000) + '\n\n…truncated for display; Copy the text copies all of it.' : text;
        t.addEventListener('focus', () => t.select());
      },
      onClose: () => { if (url) setTimeout(() => URL.revokeObjectURL(url), 60000); }
    });
    return false;
  }
  PB.saveBlocked = saveBlocked;

  /* A binary file — a PNG. Where the save capability is not available it is
     shown on screen instead, where it can be saved or copied from the
     browser's own image menu. */
  async function offerBlob(filename, blob, label) {
    await ensureDownloads();
    if (state.downloads) {
      try {
        await state.downloads.save({ filename, data: blob });
        toast('Saved ' + filename, 'ok');
        return true;
      } catch (e) {
        const code = e && e.code;
        if (code === 'declined') { toast('Save cancelled'); return false; }
        if (code === 'too_large') { toast('That image is too large to save here — choose a lower resolution.', 'bad'); return false; }
      }
    }
    /* No host save, or it refused. In an ordinary tab, with no host at all,
       the browser's own download simply works and saying so is correct. */
    if (standalone() && nativeSave(filename, blob)) { toast('Saved ' + filename, 'ok'); return true; }
    /* Otherwise: it is attempted, nothing claims it worked, and a link the
       viewer clicks themselves is put in front of them, which is allowed in
       places a synthetic click is not. */
    const url = URL.createObjectURL(blob);
    try {
      const a = document.createElement('a');
      a.href = url; a.download = filename; a.rel = 'noopener';
      a.style.display = 'none';
      document.body.appendChild(a);
      a.click();
      setTimeout(() => a.remove(), 4000);
    } catch (e) { /* the dialog below is the real answer either way */ }
    const isImg = /^image\//.test(blob.type || '');
    modal({
      title: label || filename,
      body: '<div class="note">If <b>' + esc(filename) + '</b> did not download on its own, this view would not write it.' +
        (isImg ? ' The picture itself is fine — use the download link below, or right-click (long-press on a phone) the image and choose <strong>Save image</strong>.'
               : ' Use the download link below.') +
        '</div>' +
        '<div class="row"><a class="btn primary" href="' + url + '" download="' + esc(filename) + '">Download ' + esc(filename) + '</a></div>' +
        (isImg ? '<div style="margin-top:10px;text-align:center"><img alt="' + esc(filename) + '" src="' + url +
          '" style="max-width:100%;border:1px solid var(--border);border-radius:6px"></div>' : ''),
      actions: [{ label: 'Close', onClick: () => { setTimeout(() => URL.revokeObjectURL(url), 1000); } }]
    });
    return false;
  }
  PB.offerBlob = offerBlob;

  function showTextFallback(title, text) {
    PB.modal({
      title: title,
      sub: fmt(text.length) + ' characters',
      body: '<div class="note">Downloads are not available in this view. Select the text below and copy it.</div>' +
            '<textarea class="inp" id="fbText" rows="14" spellcheck="false"></textarea>',
      actions: [{ label: 'Select all', onClick: () => { const t = $('#fbText'); t.focus(); t.select(); } },
                { label: 'Close', primary: true, close: true }],
      onOpen: () => { $('#fbText').value = text; }
    });
  }
  PB.showTextFallback = showTextFallback;

  /* ---------------- modal ---------------- */
  let modalStack = [];
  function modal(cfg) {
    const root = $('#modalRoot');
    const wrap = document.createElement('div');
    wrap.className = 'modal' + (cfg.wide ? ' wide' : '');
    wrap.innerHTML =
      '<div class="modal-head"><h3>' + esc(cfg.title) +
        (cfg.sub ? ' <span class="sub">' + esc(cfg.sub) + '</span>' : '') +
      '</h3><button class="btn icon ghost" data-close aria-label="Close">' +
      '<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" style="width:16px;height:16px"><path d="M18 6 6 18M6 6l12 12"/></svg></button></div>' +
      '<div class="modal-body">' + (cfg.body || '') + '</div>' +
      (cfg.actions ? '<div class="modal-foot"></div>' : '');
    root.innerHTML = '';
    root.appendChild(wrap);
    root.hidden = false;
    if (cfg.actions) {
      const foot = wrap.querySelector('.modal-foot');
      cfg.actions.forEach(a => {
        const b = document.createElement('button');
        b.className = 'btn' + (a.primary ? ' primary' : '');
        b.textContent = a.label;
        if (a.id) b.id = a.id;
        b.addEventListener('click', () => {
          if (a.onClick) {
            const r = a.onClick(wrap);
            if (r === false) return;
          }
          if (a.close !== false) closeModal();
        });
        foot.appendChild(b);
      });
    }
    wrap.querySelectorAll('[data-close]').forEach(b => b.addEventListener('click', closeModal));
    root.onclick = e => { if (e.target === root) closeModal(); };
    modalStack.push(cfg);
    if (cfg.onOpen) cfg.onOpen(wrap);
    return wrap;
  }
  function closeModal() {
    const root = $('#modalRoot');
    /* a dialog that made an object URL or started something needs to know
       it is going away, however it was closed — button, ×, backdrop or Esc */
    modalStack.forEach(c => { if (c && typeof c.onClose === 'function') { try { c.onClose(); } catch (e) { /* closing anyway */ } } });
    root.hidden = true; root.innerHTML = ''; modalStack = [];
  }
  PB.modal = modal; PB.closeModal = closeModal;
  document.addEventListener('keydown', e => {
    if (e.key === 'Escape' && !$('#modalRoot').hidden) closeModal();
  });

  /* ---------------- DNA methylation ----------------
     Plasmid DNA grown in ordinary E. coli strains carries Dam (GATC) and Dcm
     (CCWGG) methylation; PCR products, synthetic DNA and DNA from dam-/dcm-
     strains carry none; mammalian genomic DNA is CpG-methylated. Each record
     remembers which applies; the default follows topology.             */
  const METHYL_MODES = {
    ecoli: { label: 'E. coli dam⁺ dcm⁺ (DH5α, TOP10…)', sys: { dam: true, dcm: true } },
    none:  { label: 'Unmethylated (PCR, synthetic, dam⁻ dcm⁻)', sys: {} },
    dam:   { label: 'Dam only (dcm⁻ strain)', sys: { dam: true } },
    cpg:   { label: 'CpG-methylated (mammalian genomic)', sys: { cpg: true } }
  };
  PB.METHYL_MODES = METHYL_MODES;
  PB.methylMode = rec => (rec && METHYL_MODES[rec.methylation]) ? rec.methylation
    : (rec && rec.topology === 'circular' ? 'ecoli' : 'none');
  PB.methylSystems = rec => METHYL_MODES[PB.methylMode(rec)].sys;
  let methCache = { key: null, sets: null };
  PB.methylSets = function (rec) {
    const mode = PB.methylMode(rec);
    const key = rec.id + '|' + rec.updatedAt + '|' + rec.seq.length + '|' + rec.topology + '|' + mode;
    if (methCache.key !== key) methCache = { key, sets: B.methylatedBases(rec.seq, rec.topology === 'circular', METHYL_MODES[mode].sys) };
    return methCache.sets;
  };

  /* ---------------- enzyme computation cache ---------------- */
  let siteCache = { key: null, map: null };
  function sitesFor(rec, names) {
    const key = rec.id + '|' + rec.updatedAt + '|' + rec.seq.length + '|' + rec.topology + '|' + rec.seq.slice(0, 40) +
                '|' + names.join(',') + '|' + state.customEnzymes.length + '|' + PB.methylMode(rec);
    if (siteCache.key === key) return siteCache.map;
    const map = {};
    const meth = PB.methylSets(rec);
    for (const n of names) {
      const e = enzByName(n);
      if (!e) continue;
      try {
        map[n] = B.findSites(rec.seq, e, rec.topology === 'circular');
        map[n].forEach(s => { const m = B.methylationEffect(e, s, rec.seq.length, meth); if (m) s.methyl = m; });
      }
      catch (err) { map[n] = []; }
    }
    siteCache = { key, map };
    return map;
  }
  PB.sitesFor = sitesFor;
  PB.invalidateSites = function () { siteCache = { key: null, map: null }; };

  function allCuts(rec) {
    const map = sitesFor(rec, state.enzymes);
    const out = [];
    for (const n in map) for (const s of map[n]) out.push(Object.assign({ enzyme: n }, s));
    return out.sort((a, b) => a.cutTop - b.cutTop);
  }
  PB.allCuts = allCuts;
  // "XbaI" / "XbaI ⊘Dam" (blocked) / "BsaI ?Dcm" (impaired or uncertain)
  PB.cutMark = function (c) {
    const m = c.methyl;
    if (!m) return '';
    const by = m.by.map(k => ({ dam: 'Dam', dcm: 'Dcm', cpg: 'CpG' })[k]).join('/');
    return (m.effect === 'blocked' ? ' ⊘' : ' ?') + by;
  };

  /* ---------------- library rail ---------------- */
  function renderLibrary() {
    const list = $('#libList');
    const recs = state.records.slice().sort((a, b) =>
      (a.builtin === b.builtin) ? (b.updatedAt || 0) - (a.updatedAt || 0) : (a.builtin ? 1 : -1));
    if (!recs.length) {
      list.innerHTML = '<div class="hint" style="padding:6px 2px">No constructs yet. Import a GenBank or FASTA file, or start a new sequence.' +
        (lsGet('demoHidden', false)
          ? ' <a href="#" id="libDemoBack">Bring the demo back.</a>' : '') + '</div>';
      const back = $('#libDemoBack');
      if (back) back.addEventListener('click', e => { e.preventDefault(); PB.restoreDemo(); });
    } else {
      list.innerHTML = recs.map(r =>
        '<div class="lib-item' + (r.id === state.activeId ? ' active' : '') + '" data-id="' + r.id + '">' +
          '<span class="ring' + (r.topology === 'linear' ? ' linear' : '') + '"></span>' +
          '<span class="nm" title="' + esc(r.name) + '">' + esc(r.name) + '</span>' +
          (r.builtin ? '<span class="tagx">demo</span>' : '') +
          '<span class="bp">' + fmt(r.seq.length) + '</span>' +
          '<button class="del" data-del="' + r.id + '" title="' +
            (r.builtin ? 'Hide the demo construct' : 'Remove') + '" aria-label="Remove ' + esc(r.name) + '">&times;</button>' +
        '</div>').join('');
    }
    $('#libCount').textContent = recs.length ? recs.length : '';
    list.querySelectorAll('.lib-item').forEach(it => {
      it.addEventListener('click', e => {
        if (e.target.hasAttribute('data-del')) return;
        state.activeId = it.getAttribute('data-id');
        state.selection = null; state.selectedFeature = null;
        PB.renderAll();
      });
    });
    list.querySelectorAll('[data-del]').forEach(b => {
      b.addEventListener('click', e => {
        e.stopPropagation();
        const id = b.getAttribute('data-del');
        const rec = state.records.find(r => r.id === id);
        modal({
          title: rec.builtin ? 'Hide the demo construct?' : 'Remove this construct?',
          body: rec.builtin
            ? '<p>“' + esc(rec.name) + '” will stop appearing in your library. ' +
              'It is a built-in example rather than anything of yours, so nothing is lost — ' +
              'you can bring it back from <strong>Tools → Demo construct</strong> whenever you want it.</p>'
            : '<p>“' + esc(rec.name) + '” (' + fmt(rec.seq.length) + ' bp) will be deleted from your library. This cannot be undone.</p>',
          actions: [
            { label: 'Cancel' },
            { label: rec.builtin ? 'Hide it' : 'Remove', primary: true, onClick: () => {
              removeRecord(id);
              toast(rec.builtin ? 'Demo construct hidden' : 'Removed ' + rec.name, 'ok');
            } }
          ]
        });
      });
    });
    renderStoreState();
  }
  PB.renderLibrary = renderLibrary;

  /* ---------------- enzyme rail ---------------- */
  function enzymeRows(rec) {
    const q = state.enzQuery.trim().toLowerCase();
    let list = allEnzymes();
    if (state.enzSet === 'iis') list = list.filter(e => e.typeIIS);
    else if (state.enzSet === 'blunt') list = list.filter(e => !e.typeIIS && e.top === e.bot);
    else if (state.enzSet === 'common') list = list.filter(e => COMMON.has(e.name) || state.customEnzymes.indexOf(e) !== -1);

    /* A lab that buys from two companies cannot use the other six hundred
       enzymes, and offering them is how a design ends up resting on a tube
       nobody can order. The suppliers come from the lab profile; with none
       set, nothing is filtered. Enzymes added by hand are always kept —
       they are there because someone has them. */
    const suppliers = lsGet('suppliers', null);
    if (suppliers && suppliers.length) {
      const want = new Set(suppliers);
      list = list.filter(e => state.customEnzymes.indexOf(e) !== -1 ||
        Array.from(String(e.suppliers || '')).some(code => want.has(code)));
    }
    if (q) list = list.filter(e => e.name.toLowerCase().indexOf(q) !== -1 || e.site.toLowerCase().indexOf(q) !== -1);

    const counts = {};
    if (rec && rec.seq.length) {
      const circ = rec.topology === 'circular';
      for (const e of list) {
        try { counts[e.name] = B.findSites(rec.seq, e, circ).length; }
        catch (err) { counts[e.name] = 0; }
      }
    }
    if (state.enzSet === 'unique') list = list.filter(e => counts[e.name] === 1);
    else if (state.enzSet === 'noncut') list = list.filter(e => counts[e.name] === 0);
    list = list.slice().sort((a, b) => a.name.localeCompare(b.name));
    return { list, counts };
  }
  const COMMON = new Set(['EcoRI','BamHI','HindIII','NotI','XhoI','SalI','PstI','SphI','KpnI',
    'SacI','SmaI','XmaI','SpeI','XbaI','NheI','NcoI','NdeI','BglII','EcoRV','PvuII','ScaI',
    'StuI','DraI','ClaI','MluI','AscI','FseI','PacI','SbfI','AgeI','MfeI','BspEI','BsrGI',
    'AatII','ApaI','AvrII','EagI','NsiI','PmeI','SnaBI','NarI','BsiWI','HpaI','AflII',
    'BsaI','BsmBI','BbsI','SapI','PaqCI','AarI','BtgZI','Sau3AI','MboI','DpnI','AluI','HaeIII','RsaI','TaqI','MspI']);
  PB.COMMON = COMMON;

  function renderEnzymeRail() {
    const rec = active();
    const { list, counts } = enzymeRows(rec);
    const host = $('#enzList');
    if (!list.length) {
      host.innerHTML = '<div class="hint">No enzymes match.</div>';
    } else {
      host.innerHTML = list.map(e => {
        const c = counts[e.name];
        const cls = c === undefined ? 'c0' : (c === 0 ? 'c0' : (c === 1 ? 'c1' : 'cm'));
        return '<label class="enz-row" title="' + esc(e.name + '  ' + e.site) + '">' +
          '<input type="checkbox" data-enz="' + esc(e.name) + '"' +
            (state.enzymes.indexOf(e.name) !== -1 ? ' checked' : '') + '>' +
          '<span class="en">' + esc(e.name) + '</span>' +
          '<span class="es">' + esc(e.site) + (e.typeIIS ? ' ▸' : '') + '</span>' +
          '<span class="ec ' + cls + '">' + (c === undefined ? '' : c) + '</span>' +
        '</label>';
      }).join('');
    }
    $('#enzCount').textContent = state.enzymes.length ? state.enzymes.length + ' on' : '';
    host.querySelectorAll('[data-enz]').forEach(cb => {
      cb.addEventListener('change', () => {
        const n = cb.getAttribute('data-enz');
        if (cb.checked) { if (state.enzymes.indexOf(n) === -1) state.enzymes.push(n); }
        else state.enzymes = state.enzymes.filter(x => x !== n);
        lsSet('enzymes', state.enzymes);
        PB.invalidateSites();
        PB.renderViews(); renderEnzymeRail(); PB.renderInspector();
      });
    });
    $$('#enzChips .chip').forEach(c => {
      c.classList.toggle('on', c.getAttribute('data-eset') === state.enzSet);
    });
  }
  PB.renderEnzymeRail = renderEnzymeRail;

  /* ---------------- header ---------------- */
  function renderHeader() {
    const rec = active();
    const nameEl = $('#recName'), metaEl = $('#recMeta');
    if (!rec) {
      nameEl.value = ''; nameEl.placeholder = 'No construct loaded';
      metaEl.innerHTML = '<span class="pill">import or create a sequence to begin</span>';
      return;
    }
    nameEl.value = rec.name;
    const gc = rec.seq.length ? B.gcContent(rec.seq) : 0;
    const cuts = rec.seq.length ? allCuts(rec).length : 0;
    metaEl.innerHTML =
      '<span>' + fmt(rec.seq.length) + ' bp</span>' +
      '<span class="pill">' + rec.topology + '</span>' +
      '<span>' + f1(gc) + '% GC</span>' +
      '<span class="pill">' + rec.features.length + ' feature' + (rec.features.length === 1 ? '' : 's') + '</span>' +
      (state.enzymes.length ? '<span class="pill">' + cuts + ' site' + (cuts === 1 ? '' : 's') + '</span>' : '') +
      (featureSelection(rec).length > 1
        ? '<span style="color:var(--accent)">' + featureSelection(rec).length + ' features selected (' +
          fmt(featureSelectionPieces(rec).reduce((t, p) => t + (p.end - p.start), 0)) + ' bp) · ' +
          esc(trim(featureSelection(rec).map(i => rec.features[i].label).join(', '), 90)) + '</span>'
        : normSel(rec, state.selection) ? '<span style="color:var(--accent)">sel ' + selRangeText(rec, state.selection) +
          ' (' + fmt(normSel(rec, state.selection).length) + ' bp)</span>' : '');
  }
  PB.renderHeader = renderHeader;

  /* ---------------- inspector ---------------- */
  function renderInspector() {
    const host = $('#inspector');
    const rec = active();
    if (!rec) { host.innerHTML = '<div class="insp-sec"><div class="hint">Nothing loaded.</div></div>'; return; }
    let html = '';

    // selection / feature detail
    const sel = state.selection;
    if (state.selectedFeature != null && rec.features[state.selectedFeature]) {
      const f = rec.features[state.selectedFeature];
      const sub = featSeq(rec, f);
      html += '<section class="insp-sec"><h2>Feature</h2>' +
        '<div style="display:flex;align-items:center;gap:7px;margin-bottom:7px">' +
        '<span class="sw" style="width:12px;height:12px;border-radius:3px;background:' + featColor(f) + '"></span>' +
        '<b style="font-size:13.5px">' + esc(f.label) + '</b></div>' +
        '<dl class="kv">' +
        '<dt>Type</dt><dd>' + esc(f.type) + '</dd>' +
        '<dt>Location</dt><dd>' + esc(PB.locationText(rec, f)) + '</dd>' +
        '<dt>Length</dt><dd>' + fmt(PB.featureLen(rec, f)) + ' bp' +
          (f.parts ? ' in ' + f.parts.length + ' parts' : '') + '</dd>' +
        '<dt>Strand</dt><dd>' + (f.mixed ? 'mixed (trans-spliced)' : (f.strand === -1 ? 'reverse (−)' : 'forward (+)')) + '</dd>' +
        '<dt>GC</dt><dd>' + f1(B.gcContent(sub)) + '%</dd>' +
        '</dl>';
      if (/CDS|gene|ORF|exon/i.test(f.type) && sub.length >= 3) {
        const prot = PB.codonMap(rec, state.selectedFeature).map(c => c.aa).join('');
        html += '<h2 style="margin-top:10px">Translation</h2><div class="seqblock">' + esc(prot) + '</div>' +
          '<div class="hint">' + prot.replace(/\*$/, '').length + ' aa' +
          (prot.indexOf('*') !== -1 && prot.indexOf('*') < prot.length - 1 ? ' · <span style="color:var(--bad)">internal stop codon</span>' : '') +
          (sub.length % 3 ? ' · <span style="color:var(--warn)">not a multiple of 3</span>' : '') +
          ((f.quals && f.quals.translation && f.quals.translation !== prot.replace(/\*$/, '') &&
            f.quals.translation.slice(1) !== prot.replace(/\*$/, '').slice(1))
            ? ' · <span style="color:var(--bad)">differs from the /translation stored in the file</span>' : '') +
          '</div>';
      }
      html += '<div class="row" style="margin-top:9px;gap:5px">' +
        '<button class="btn sm" data-act="selfeat">Select range</button>' +
        '<button class="btn sm" data-act="copyfeat">Copy</button>' +
        '<button class="btn sm" data-act="editfeat">Edit</button>' +
        '</div></section>';
    } else if (normSel(rec, sel)) {
      const sub = selectionSeq(rec, sel);
      html += '<section class="insp-sec"><h2>Selection</h2><dl class="kv">' +
        '<dt>Range</dt><dd>' + selRangeText(rec, sel) + '</dd>' +
        '<dt>Length</dt><dd>' + fmt(sub.length) + ' bp</dd>' +
        '<dt>GC</dt><dd>' + f1(B.gcContent(sub)) + '%</dd>' +
        '<dt>Tm</dt><dd' + (tmShown(sub) ? '' : ' title="Nearest-neighbour Tm describes an oligo of ' + TM_MIN + '–' + TM_MAX +
          ' nt; outside that range the number would not mean anything"') + '>' +
          (tmShown(sub) ? f1(B.tmNN(sub, conditions())) + ' °C' : '— <span class="hint" style="margin:0">(' + TM_MIN + '–' + TM_MAX + ' nt only)</span>') + '</dd>' +
        '<dt>MW (ds)</dt><dd>' + fmt(Math.round(B.dsDNAWeight(sub))) + ' g/mol</dd>' +
        '</dl>' +
        '<h2 style="margin-top:10px">Sequence</h2><div class="seqblock">' + esc(sub.length > 2000 ? sub.slice(0, 2000) + '…' : sub) + '</div>' +
        '<div class="row" style="margin-top:8px;gap:5px">' +
        '<button class="btn sm" data-act="copysel">Copy</button>' +
        '<button class="btn sm" data-act="copyrc">Copy rev-comp</button>' +
        '<button class="btn sm" data-act="translate">Translate</button>' +
        '<button class="btn sm" data-act="mkfeat">Make feature</button>' +
        '<button class="btn sm" data-act="mkprimer">Make primer</button>' +
        '</div></section>';
    }

    // composition
    const counts = B.baseCounts(rec.seq);
    html += '<section class="insp-sec"><h2>Composition</h2><dl class="kv">' +
      '<dt>Length</dt><dd>' + fmt(rec.seq.length) + ' bp</dd>' +
      '<dt>GC</dt><dd>' + f1(B.gcContent(rec.seq)) + '%</dd>' +
      '<dt>A / T</dt><dd>' + fmt(counts.A) + ' / ' + fmt(counts.T) + '</dd>' +
      '<dt>G / C</dt><dd>' + fmt(counts.G) + ' / ' + fmt(counts.C) + '</dd>' +
      (counts.N ? '<dt>Ambiguous</dt><dd>' + fmt(counts.N) + '</dd>' : '') +
      '<dt>MW (ds)</dt><dd>' + (rec.seq.length <= 60000 ? fmt(Math.round(B.dsDNAWeight(rec.seq) / 1000)) + ' kDa' : '—') + '</dd>' +
      '</dl></section>';

    // legend
    const used = {};
    rec.features.forEach(f => { used[featColorVar(f)] = (used[featColorVar(f)] || 0) + 1; });
    const LEG = { 'f-promoter': 'Promoter / enhancer', 'f-cds': 'CDS / gene', 'f-terminator': 'Terminator / polyA',
      'f-origin': 'Origin of replication', 'f-tag': 'Tag', 'f-marker': 'Selection marker',
      'f-reg': 'Regulatory / binding', 'f-recomb': 'Recombination site', 'f-primer': 'Primer binding', 'f-misc': 'Other' };
    const legKeys = Object.keys(used);
    if (legKeys.length) {
      html += '<section class="insp-sec"><h2>Feature key</h2><div class="legend">' +
        legKeys.map(k => '<div class="legend-item"><span class="sw" style="background:var(--' + k + ')"></span>' +
          esc(LEG[k] || k) + ' <span style="color:var(--fg-faint)">· ' + used[k] + '</span></div>').join('') +
        '</div></section>';
    }

    // what this DNA is, for methylation-sensitive digests
    const mm = PB.methylMode(rec);
    html += '<section class="insp-sec"><h2>DNA methylation</h2>' +
      '<select class="inp" id="condMeth" title="Which methylation this DNA carries">' +
      Object.keys(METHYL_MODES).map(k => '<option value="' + k + '"' + (k === mm ? ' selected' : '') + '>' + esc(METHYL_MODES[k].label) + '</option>').join('') +
      '</select><div class="hint">Sites blocked by this methylation are not cut in digests, the gel or cloning, and are flagged on the map and in the Enzymes tab. Sensitivities from the NEB methylation chart; for enzymes NEB does not sell it is shown as unknown.' +
      (rec.methylation ? '' : ' Default for a ' + rec.topology + ' record.') + '</div></section>';

    // conditions
    const c = state.tmCond;
    html += '<section class="insp-sec"><h2>Tm conditions</h2>' +
      '<div class="row" style="gap:6px">' +
      '<div class="field narrow"><label>Oligo µM</label><input class="inp mono" id="condOligo" type="number" step="0.05" min="0.01" value="' + c.oligo + '"></div>' +
      '<div class="field narrow"><label>Na⁺ mM</label><input class="inp mono" id="condNa" type="number" step="5" min="0" value="' + c.Na + '"></div>' +
      '</div><div class="row" style="gap:6px">' +
      '<div class="field narrow"><label>Mg²⁺ mM</label><input class="inp mono" id="condMg" type="number" step="0.5" min="0" value="' + c.Mg + '"></div>' +
      '<div class="field narrow"><label>dNTP mM</label><input class="inp mono" id="conddNTP" type="number" step="0.1" min="0" value="' + c.dNTP + '"></div>' +
      '</div>' +
      '<div class="hint">Nearest-neighbour thermodynamics (SantaLucia 1998) with the Owczarzy monovalent and Mg²⁺ salt corrections.</div>' +
      '<div class="field" style="margin-top:8px"><label>Genetic code</label><select class="inp" id="condCode">' +
      B.geneticCodes().map(g => '<option value="' + g.id + '"' + (g.id === state.codeTable ? ' selected' : '') + '>' + esc(g.id + ' · ' + g.name) + '</option>').join('') +
      '</select></div></section>';

    host.innerHTML = html;

    const bind = (act, fn) => host.querySelectorAll('[data-act="' + act + '"]').forEach(b => b.addEventListener('click', fn));
    bind('copysel', () => copyText(selectionSeq(rec, sel), 'Selection'));
    bind('copyrc', () => copyText(B.revcomp(selectionSeq(rec, sel)), 'Reverse complement'));
    bind('translate', () => PB.toolTranslate && PB.toolTranslate());
    bind('mkfeat', () => PB.featureDialog && PB.featureDialog(null));
    bind('mkprimer', () => PB.primerFromSelection && PB.primerFromSelection());
    bind('selfeat', () => {
      const f = rec.features[state.selectedFeature];
      state.selection = { start: f.start, end: f.end };
      PB.renderHeader(); PB.renderViews(); renderInspector();
    });
    bind('copyfeat', () => copyText(featSeq(rec, rec.features[state.selectedFeature]), 'Feature sequence'));
    bind('editfeat', () => PB.featureDialog && PB.featureDialog(state.selectedFeature));

    const onCond = () => {
      state.tmCond = {
        oligo: parseFloat($('#condOligo').value) || 0.25,
        Na: parseFloat($('#condNa').value) || 50,
        Mg: parseFloat($('#condMg').value) || 0,
        dNTP: parseFloat($('#conddNTP').value) || 0
      };
      lsSet('tmCond', state.tmCond);
      PB.renderViews(); renderInspector();
    };
    ['#condOligo', '#condNa', '#condMg', '#conddNTP'].forEach(id => {
      const el = $(id); if (el) el.addEventListener('change', onCond);
    });
    const meth = $('#condMeth');
    if (meth) meth.addEventListener('change', () => {
      rec.methylation = meth.value;
      addHistory(rec, 'Methylation set', METHYL_MODES[meth.value].label);
      PB.touch(rec); PB.invalidateSites(); PB.renderAll();
    });
    const code = $('#condCode');
    if (code) code.addEventListener('change', () => {
      state.codeTable = parseInt(code.value, 10) || 1;
      lsSet('codeTable', state.codeTable);
      PB.renderViews(); renderInspector();
    });
  }
  PB.renderInspector = renderInspector;

  function featSeq(rec, f) {
    if (!f || !rec) return '';
    return B.featureSequence(rec.seq, f);
  }
  PB.featureLen = (rec, f) => B.featureLength(f, rec.seq.length);
  PB.locationText = (rec, f) => B.formatLocation(f, rec.seq.length) || '';

  /* Codon map: for a coding feature, every codon with the genomic (top-strand)
     positions of its three bases and its amino acid. Built through joined
     parts, origin wraps and either strand, honouring /codon_start and
     /transl_table, so translations drawn anywhere in the app agree with the
     feature's real spliced sequence. Cached per record version.          */
  const codonCache = new Map();
  function featurePositions(rec, f) {
    const L = rec.seq.length;
    const out = [], strandOf = [];
    if (f.mixed && f.parts) {
      f.parts.forEach(p => {
        if (p.strand === -1) for (let x = p.end - 1; x >= p.start; x--) { out.push(x); strandOf.push(-1); }
        else for (let x = p.start; x < p.end; x++) { out.push(x); strandOf.push(1); }
      });
      return { pos: out, strand: strandOf };
    }
    B.featurePieces(f, L).forEach(p => { for (let x = p.start; x < p.end; x++) out.push(x); });
    if (f.strand === -1) out.reverse();
    return { pos: out, strand: out.map(() => f.strand === -1 ? -1 : 1) };
  }
  PB.featurePositions = featurePositions;
  PB.codonMap = function (rec, idx) {
    const f = rec.features[idx];
    if (!f) return [];
    const key = rec.id + ':' + rec.updatedAt + ':' + rec.seq.length + ':' + idx + ':' + state.codeTable +
                ':' + f.start + ':' + f.end + ':' + f.strand;
    if (codonCache.has(key)) return codonCache.get(key);
    const fp = featurePositions(rec, f);
    const fs = B.featureSequence(rec.seq, f);
    const cs = Math.max(1, Math.min(3, parseInt((f.quals && f.quals.codon_start) || '1', 10) || 1)) - 1;
    const table = parseInt((f.quals && f.quals.transl_table) || state.codeTable, 10) || state.codeTable;
    const out = [];
    for (let i = cs, k = 0; i + 3 <= fs.length; i += 3, k++) {
      out.push({ k: k, pos: [fp.pos[i], fp.pos[i + 1], fp.pos[i + 2]],
                 aa: B.translate(fs.substr(i, 3), { table: table }) });
    }
    if (codonCache.size > 400) codonCache.clear();
    codonCache.set(key, out);
    return out;
  };
  // Put an edited feature sequence (same length, e.g. after silent mutation or
  // codon optimisation) back into the record base by base.
  PB.writeFeatureSequence = function (rec, f, newFeatSeq) {
    const fp = featurePositions(rec, f);
    if (newFeatSeq.length !== fp.pos.length) throw new Error('Edited feature changed length.');
    const arr = rec.seq.split('');
    for (let i = 0; i < fp.pos.length; i++) {
      const b = newFeatSeq[i];
      arr[fp.pos[i]] = fp.strand[i] === -1 ? B.complement(b) : b;
    }
    return arr.join('');
  };
  PB.featSeq = featSeq;

  /* ---------------- view tabs ---------------- */
  const VIEWS = [
    ['map', 'Map'], ['seq', 'Sequence'], ['features', 'Features'], ['enzymes', 'Enzymes'],
    ['orfs', 'ORFs'], ['primers', 'Primers'], ['gel', 'Gel'], ['verify', 'Verify'],
    ['align', 'Align'], ['history', 'History']
  ];
  function renderTabs() {
    const host = $('#viewTabs');
    host.innerHTML = VIEWS.map(([k, label]) =>
      '<button role="tab" data-v="' + k + '"' + (state.view === k ? ' class="on" aria-selected="true"' : '') + '>' + label + '</button>').join('');
    host.querySelectorAll('[data-v]').forEach(b => {
      b.addEventListener('click', () => {
        state.view = b.getAttribute('data-v');
        lsSet('view', state.view);
        renderTabs(); PB.renderViews();
      });
    });
    $$('.view').forEach(v => v.classList.toggle('on', v.getAttribute('data-view') === state.view));
  }
  PB.renderTabs = renderTabs;

  /* ---------------- top-level render ---------------- */
  PB.renderAll = function () {
    renderHeader(); renderLibrary(); renderEnzymeRail(); renderTabs();
    PB.renderViews(); renderInspector(); renderUndo();
  };

  /* the buttons say what they would undo or redo */
  function renderUndo() {
    const u = $('#btnUndo'), r = $('#btnRedo');
    if (!u || !r) return;
    const st = PB.undoState();
    u.disabled = !st.canUndo; r.disabled = !st.canRedo;
    u.title = st.canUndo ? 'Undo: ' + st.undoLabel + '  (Ctrl/\u2318 Z)' : 'Nothing to undo for this construct';
    r.title = st.canRedo ? 'Redo: ' + st.redoLabel + '  (Ctrl/\u2318 Shift Z)' : 'Nothing to redo';
  }
  PB.renderUndo = renderUndo;

  /* ---------------- demo construct ----------------
     Built from short elements whose sequences are individually verified
     (promoters, tags, linkers, the pUC19 polylinker). It is a synthetic
     teaching construct, clearly labelled, not a real catalogue plasmid. */
  function buildDemo() {
    const T7 = 'TAATACGACTCACTATAG';
    const lacO = 'AATTGTGAGCGGATAACAATT';
    const RBS = 'AAGGAGGAGATATACC';
    const his = 'CATCATCATCATCATCAT';
    const linker = 'GGTGGTGGTGGTTCTGGTGGTGGTGGTTCT';
    const egfpN = 'ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTAC';
    const mcs = 'GAATTCGAGCTCGGTACCCGGGGATCCTCTAGAGTCGACCTGCAGGCATGCAAGCTT';
    const stop = 'TAATAA';
    const T7term = 'CTAGCATAACCCCTTGGGGCCTCTAAACGGGTCTTGAGGGGTTTTTTG';
    const flag = 'GATTACAAGGATGACGACGATAAG';
    const spacerA = 'CGTAGCTAGCTTGGATCATAGCGTACGATCAGCTACGGATCCATAGC'.slice(0, 24);
    const pad = s => s;

    const parts = [];
    let seq = '';
    const feats = [];
    function put(s, type, label, strand) {
      if (label) feats.push({ type, start: seq.length, end: seq.length + s.length, strand: strand || 1, label });
      seq += s;
    }
    put('GGTCTCAGGAGA'.slice(0, 6) + 'CTAGCCATCGGCATTACGGCA', 'misc_feature', '', 1);
    put(T7, 'promoter', 'T7 promoter', 1);
    put('CC', 'misc_feature', '', 1);
    put(lacO, 'protein_bind', 'lac operator', 1);
    put('GAATTGTGAGCGGA'.slice(0, 8), 'misc_feature', '', 1);
    put(RBS, 'RBS', 'ribosome binding site', 1);
    put('ATG', 'CDS', 'start codon', 1);
    put(his, 'protein_tag', 'His6 tag', 1);
    put(linker, 'misc_feature', 'flexible linker', 1);
    put(egfpN, 'CDS', 'EGFP', 1);
    put(flag, 'protein_tag', 'FLAG tag', 1);
    put(stop, 'misc_feature', 'stop', 1);
    put(mcs, 'misc_feature', 'pUC19 polylinker (MCS)', 1);
    put(T7term, 'terminator', 'T7 terminator', 1);
    put('GCTAGTTATTGCTCAGCGGTGGCAGCAGCCAACTCAGCTTCCTTTCGGGCTTTGTTAGCAGCCGGATC', 'misc_feature', '', 1);
    put('ATGAGTATTCAACATTTCCGTGTCGCCCTTATTCCCTTTTTTGCGGCATTTTGCCTTCCTGTTTTTGCTCACCCAGAAACGCTGGTGAAAGTAAAAGATGCTGAAGATCAGTTGGGTGCACGAGTGGGTTACATCGAACTGGATCTCAACAGCGGTAAGATCCTTGAGAGTTTTCGCCCCGAAGAACGTTTTCCAATGATGAGCACTTTTAAAGTTCTGCTATGTGGCGCGGTATTATCCCGTATTGACGCCGGGCAAGAGCAACTCGGTCGCCGCATACACTATTCTCAGAATGACTTGGTTGAGTACTCACCAGTCACAGAAAAGCATCTTACGGATGGCATGACAGTAAGAGAATTATGCAGTGCTGCCATAACCATGAGTGATAACACTGCGGCCAACTTACTTCTGACAACGATCGGAGGACCGAAGGAGCTAACCGCTTTTTTGCACAACATGGGGGATCATGTAACTCGCCTTGATCGTTGGGAACCGGAGCTGAATGAAGCCATACCAAACGACGAGCGTGACACCACGATGCCTGTAGCAATGGCAACAACGTTGCGCAAACTATTAACTGGCGAACTACTTACTCTAGCTTCCCGGCAACAATTAATAGACTGGATGGAGGCGGATAAAGTTGCAGGACCACTTCTGCGCTCGGCCCTTCCGGCTGGCTGGTTTATTGCTGATAAATCTGGAGCCGGTGAGCGTGGGTCTCGCGGTATCATTGCAGCACTGGGGCCAGATGGTAAGCCCTCCCGTATCGTAGTTATCTACACGACGGGGAGTCAGGCAACTATGGATGAACGAAATAGACAGATCGCTGAGATAGGTGCCTCACTGATTAAGCATTGGTAA', 'CDS', 'AmpR (bla)', 1);
    put('CGCCAGCAACCGCACCTGTGGCGCCGGTGATGCCGGCCACGATGCGTCCGGCGTAGAGGATC', 'misc_feature', '', 1);
    put('TTGAGATCCTTTTTTTCTGCGCGTAATCTGCTGCTTGCAAACAAAAAAACCACCGCTACCAGCGGTGGTTTGTTTGCCGGATCAAGAGCTACCAACTCTTTTTCCGAAGGTAACTGGCTTCAGCAGAGCGCAGATACCAAATACTGTTCTTCTAGTGTAGCCGTAGTTAGGCCACCACTTCAAGAACTCTGTAGCACCGCCTACATACCTCGCTCTGCTAATCCTGTTACCAGTGGCTGCTGCCAGTGGCGATAAGTCGTGTCTTACCGGGTTGGACTCAAGACGATAGTTACCGGATAAGGCGCAGCGGTCGGGCTGAACGGGGGGTTCGTGCACACAGCCCAGCTTGGAGCGAACGACCTACACCGAACTGAGATACCTACAGCGTGAGCTATGAGAAAGCGCCACGCTTCCCGAAGGGAGAAAGGCGGACAGGTATCCGGTAAGCGGCAGGGTCGGAACAGGAGAGCGCACGAGGGAGCTTCCAGGGGGAAACGCCTGGTATCTTTATAGTCCTGTCGGGTTTCGCCACCTCTGACTTGAGCGTCGATTTTTGTGATGCTCGTCAGGGGGGCGGAGCCTATGGAAAAACGCCAGCAACGCGG', 'rep_origin', 'pUC/ori (rep origin)', 1);
    const rec = makeRecord({
      id: 'demo-construct',
      name: 'pPB-T7-His6-EGFP (demo)',
      seq: seq,
      topology: 'circular',
      features: feats.filter(f => f.label),
      description: 'Synthetic demonstration construct assembled from individually verified elements (T7 promoter, lac operator, RBS, His6, EGFP, FLAG, pUC19 polylinker, T7 terminator, bla, pUC origin). Not a real catalogue plasmid — import your own GenBank files for real work.',
      builtin: true
    });
    addHistory(rec, 'Created', 'Built in demonstration construct');
    return rec;
  }
  PB.buildDemo = buildDemo;
})(window.PB);
