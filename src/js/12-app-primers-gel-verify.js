
/* ============================================================================
   PlasmidBench — views: primers, gel, verify (sequencing), align, history
   ============================================================================ */
(function (PB) {
  'use strict';
  const B = window.Bio;
  const $ = PB.$, $$ = PB.$$, esc = PB.esc, fmt = PB.fmt, f1 = PB.f1, f2 = PB.f2;
  const state = PB.state;
  const emptyState = PB.emptyState;

  /* ============================== PRIMERS ============================== */
  const TAIL_PRESETS = [
    ['', 'none'], ['GGGGAATTC', 'EcoRI (GAATTC)'], ['GGGGGATCC', 'BamHI (GGATCC)'],
    ['GGGAAGCTT', 'HindIII (AAGCTT)'], ['GGGCTCGAG', 'XhoI (CTCGAG)'],
    ['GGGGTCGAC', 'SalI (GTCGAC)'], ['GGGCATATG', 'NdeI (CATATG)'],
    ['GGGCCATGG', 'NcoI (CCATGG)'], ['GGGTCTAGA', 'XbaI (TCTAGA)'],
    ['GGGGGTCTCA', 'BsaI (GGTCTC)'], ['GGGCGTCTCA', 'BsmBI (CGTCTC)'],
    ['CACC', 'Directional TOPO (CACC)']
  ];

  function renderPrimers() {
    const host = $('#view-primers');
    const rec = PB.active();
    if (!rec) { host.innerHTML = emptyState('No construct loaded.'); return; }
    const sel = state.selection;
    const d = state.primerDesign || {};
    const start = d.start != null ? d.start : (sel ? sel.start : 0);
    const end = d.end != null ? d.end : (sel ? sel.end : Math.min(rec.seq.length, 600));

    let html = '<div class="panel"><h3>Design a primer pair ' +
      '<span class="sub">optimises length for the target Tm, then checks structure</span></h3>' +
      '<div class="row">' +
      '<div class="field narrow"><label>From (bp)</label><input class="inp mono" id="pdStart" type="number" min="1" max="' + rec.seq.length + '" value="' + (start + 1) + '"></div>' +
      '<div class="field narrow"><label>To (bp)</label><input class="inp mono" id="pdEnd" type="number" min="1" max="' + rec.seq.length + '" value="' + end + '"></div>' +
      '<div class="field narrow"><label>Target Tm °C</label><input class="inp mono" id="pdTm" type="number" step="0.5" value="' + (d.targetTm || 62) + '"></div>' +
      '<div class="field narrow"><label>Min len</label><input class="inp mono" id="pdMin" type="number" min="12" value="' + (d.minLen || 18) + '"></div>' +
      '<div class="field narrow"><label>Max len</label><input class="inp mono" id="pdMax" type="number" max="60" value="' + (d.maxLen || 30) + '"></div>' +
      (sel && sel.end > sel.start ? '<div class="field" style="flex:0 0 auto;align-self:flex-end"><button class="btn sm" id="pdUseSel">Use selection</button></div>' : '') +
      '</div><div class="row">' +
      '<div class="field"><label>Forward 5′ tail</label><select class="inp" id="pdFTail">' +
        TAIL_PRESETS.map(([v, l]) => '<option value="' + v + '"' + ((d.fTail || '') === v ? ' selected' : '') + '>' + esc(l) + '</option>').join('') +
      '</select></div>' +
      '<div class="field"><label>Reverse 5′ tail</label><select class="inp" id="pdRTail">' +
        TAIL_PRESETS.map(([v, l]) => '<option value="' + v + '"' + ((d.rTail || '') === v ? ' selected' : '') + '>' + esc(l) + '</option>').join('') +
      '</select></div>' +
      '<div class="field" style="flex:0 0 auto;align-self:flex-end"><button class="btn primary" id="pdGo">Design</button></div>' +
      '</div>';

    if (d.result) {
      const r = d.result;
      if (!r.ok) {
        html += '<div class="note bad">' + esc(r.error) + '</div>';
      } else {
        html += '<div class="row" style="margin-top:12px;align-items:stretch">' +
          primerCard(r.forward, 'Forward') + primerCard(r.reverse, 'Reverse') + '</div>' +
          '<div class="row"><div class="field wide"><dl class="kv">' +
          '<dt>Tm difference</dt><dd>' + f1(r.tmDifference) + ' °C ' +
            (r.tmDifference <= 3 ? '<span class="badge ok">good</span>' : r.tmDifference <= 5 ? '<span class="badge warn">acceptable</span>' : '<span class="badge bad">rebalance</span>') + '</dd>' +
          '<dt>Suggested Ta</dt><dd>' + r.annealingSuggestion + ' °C</dd>' +
          '<dt>Product</dt><dd>' + (r.amplicon ? fmt(r.amplicon.length) + ' bp' : '<span style="color:var(--bad)">no product predicted</span>') + '</dd>' +
          '<dt>Cross-dimer</dt><dd>' + r.heterodimer.matches + ' bp complementarity' +
            (r.heterodimer.threePrime ? ' <span class="badge bad">at a 3′ end</span>' :
             r.heterodimer.matches >= 8 ? ' <span class="badge warn">check</span>' : ' <span class="badge ok">clear</span>') + '</dd>' +
          '</dl></div></div>' +
          '<div class="row"><button class="btn sm" id="pdSave">Save both to primer library</button>' +
          '<button class="btn sm" id="pdCopy">Copy both</button>' +
          (r.amplicon ? '<button class="btn sm" id="pdAmp">Save amplicon as new sequence</button>' : '') +
          '</div>';
      }
    }
    html += '</div>';

    // library
    html += '<div class="panel"><h3>Primer library <span class="sub">' + state.primers.length + ' stored</span>' +
      '<button class="btn sm" id="pAdd">Add primer</button>' +
      (state.primers.length ? '<button class="btn sm" id="pExport" title="Write them out with their melting temperatures, for ordering or for the record">Export</button>' : '') + '</h3>';
    if (!state.primers.length) {
      html += emptyState('No primers stored yet. Design a pair above, or add one you already have.');
    } else {
      html += '<div class="scrollx"><table class="tbl"><thead><tr>' +
        '<th>Name</th><th>Sequence</th><th class="num">nt</th>' +
        '<th class="num" title="For a primer with a 5\u2032 overhang this is the Tm of the part that anneals, ' +
          'which is the temperature the first cycles need; the whole primer\u2019s Tm is shown under it.">Tm</th>' +
        '<th class="num">GC</th>' +
        '<th>Hairpin</th><th>Self-dimer</th><th>Binds here</th><th></th></tr></thead><tbody>' +
        state.primers.map((p, i) => {
          const allHits = rec ? PB.primerHits(rec) : [];
          const T = PB.primerTmOf(p, rec, allHits, i);
          const tm = T.tm;
          const hp = B.hairpin(p.seq);
          const sd = B.selfDimer(p.seq);
          const hits = allHits.filter(h => h.primer === i);
          return '<tr>' +
            '<td style="padding:2px 6px"><input class="inp" data-prename="' + i + '" value="' + esc(p.name) +
              '" spellcheck="false" aria-label="Primer name" style="font-size:12px;padding:3px 6px;min-width:110px;width:100%"></td>' +
            '<td class="seqcell" style="max-width:260px">' + esc(p.seq) + '</td>' +
            '<td class="num">' + p.seq.length + '</td>' +
            '<td class="num" title="' + esc(tmTitle(T)) + '">' +
              '<button class="btn sm ghost" data-ptm="' + i + '" style="padding:0;border:none;background:none;' +
                'font:inherit;color:inherit;cursor:pointer;display:block;width:100%;text-align:right;' +
                'line-height:1.25">' +
                (isFinite(tm) ? f1(tm) : '\u2014') +
                (T.tailed ? '<div class="hint" style="margin:0;font-size:10px;white-space:nowrap">full ' +
                   (isFinite(T.tmFull) ? f1(T.tmFull) : '\u2014') + '</div>' +
                   (T.basis === 'arm' || T.basis === 'construct'
                     ? '<div class="hint" style="margin:0;font-size:9.5px;white-space:nowrap;opacity:.75">' +
                       esc(T.tailLen) + ' nt arm, found</div>' : '') +
                   (T.basis === 'by hand'
                     ? '<div class="hint" style="margin:0;font-size:9.5px;white-space:nowrap;opacity:.75">set by you</div>' : '')
                 : '') +
              '</button>' +
            '</td>' +
            '<td class="num">' + f1(B.gcContent(p.seq)) + '</td>' +
            '<td>' + (hp.stemLength >= 4 ? '<span class="badge warn">' + hp.stemLength + ' bp stem</span>' : '<span class="badge ok">none</span>') + '</td>' +
            '<td>' + (sd.matches >= 8 ? '<span class="badge ' + (sd.threePrime ? 'bad' : 'warn') + '">' + sd.matches + ' bp</span>' : '<span class="badge ok">' + sd.matches + ' bp</span>') + '</td>' +
            '<td style="font-size:11.5px">' + (hits.length ? hits.slice(0, 4).map(h => {
              const k = allHits.indexOf(h);
              const L0 = rec.seq.length;
              const span = h.wraps ? fmt(h.start + 1) + '–' + fmt(L0) + ', 1–' + fmt(h.end) : fmt(h.start + 1) + '–' + fmt(h.end);
              return '<button class="btn sm ghost mono" style="padding:1px 5px" data-psite="' + k + '" title="' +
                esc(PB.primerSiteText(rec, h)) + '">' + (h.strand === 1 ? '\u2192 ' : '\u2190 ') + span + '</button>' +
                (h.mismatches.length ? '<span class="badge warn">' + h.mismatches.length + ' mm</span>' : '') +
                (h.tail ? '<span class="badge neutral">' + h.tail + ' nt tail</span>' : '');
            }).join(' ') + (hits.length > 4 ? ' <span class="badge warn">+' + (hits.length - 4) + ' more</span>' : '') +
              (hits.length > 1 ? ' <span class="badge warn">anneals ' + hits.length + '\u00d7</span>' : '')
              : '<span class="hint" style="margin:0">does not anneal here</span>') + '</td>' +
            '<td style="white-space:nowrap">' +
              '<button class="btn sm ghost" data-pcopy="' + i + '" title="Copy">⧉</button>' +
              '<button class="btn sm ghost" data-pdel="' + i + '" title="Delete">×</button>' +
            '</td></tr>';
        }).join('') + '</tbody></table></div>';
    }
    html += '</div>';
    host.innerHTML = html;
    wirePrimers(host, rec);
  }
  PB.renderPrimers = renderPrimers;

  function primerCard(p, which) {
    const warn = [];
    if (p.hairpin && p.hairpin.stemLength >= 5) warn.push('hairpin: ' + p.hairpin.stemLength + ' bp stem, ΔG ' + f1(p.hairpin.dG) + ' kcal/mol');
    if (p.selfDimer && p.selfDimer.matches >= 8) warn.push('self-dimer: ' + p.selfDimer.matches + ' bp' + (p.selfDimer.threePrime ? ' involving the 3′ end' : ''));
    if (p.runs >= 5) warn.push('run of ' + p.runs + ' identical bases');
    if (p.gc < 35 || p.gc > 70) warn.push('GC ' + f1(p.gc) + '% is outside 35–70%');
    return '<div class="field" style="flex:1 1 280px">' +
      '<label>' + which + ' primer' + (p.tail ? ' (with 5′ tail)' : '') + '</label>' +
      '<div class="seqblock" style="max-height:none">' +
        (p.tail ? '<span style="color:var(--accent);font-weight:600">' + esc(p.tail) + '</span>' : '') +
        esc(p.anneal) + '</div>' +
      '<dl class="kv" style="margin-top:7px">' +
      '<dt>Length</dt><dd>' + p.length + ' nt' + (p.tail ? ' (' + p.anneal.length + ' annealing)' : '') + '</dd>' +
      '<dt>Tm (annealing)</dt><dd>' + f1(p.tm) + ' °C</dd>' +
      '<dt>GC</dt><dd>' + f1(p.gc) + '%</dd>' +
      '<dt>3′ end</dt><dd>' + esc(p.gc3) + '</dd>' +
      (p.thermo ? '<dt>ΔG°37</dt><dd>' + f2(p.thermo.dG37) + ' kcal/mol</dd>' : '') +
      '</dl>' +
      (warn.length ? '<div class="note warn" style="margin:7px 0 0">' + warn.map(esc).join('<br>') + '</div>'
                   : '<div class="note ok" style="margin:7px 0 0">No structural warnings.</div>') +
      '</div>';
  }

  /* ---------------- library primers on the construct ----------------
     Every primer in the library that anneals to the active construct, found
     with Bio.primerSites (3' end seeded, internal mismatches tolerated, 5'
     tails reported as tails). They are drawn on the map and the sequence as
     display-only primer_bind features: they are never written into the
     construct, and a primer the file already annotates at the same spot is
     not drawn twice. */
  /* The two melting temperatures of a stored primer. The annealing portion
     comes from the design that made it, when it was made here; failing
     that, from how the primer binds the construct on screen, which is the
     only other thing that can say which bases pair. With neither, the
     whole primer is all there is to measure and `tailed` is false, so no
     caller mistakes a full-length figure for an annealing one. */
  /* Worked out once for the whole library: every stored primer against
     every loaded construct. Primers saved before their design was recorded
     carry no annealing portion, so this is the only thing standing between
     them and a full-length Tm that would be read as an annealing one. */
  PB.overhangScan = function () {
    const key = state.primers.map(p => p.name + '=' + p.seq + '/' + (p.anneal || '') + '/' + (p.annealLen || '')).join('|') +
      '#' + state.records.map(r => r.id + ':' + r.updatedAt + ':' + r.seq.length).join('|');
    if (state.overhangCache && state.overhangCache.key === key) return state.overhangCache.map;
    const map = {};
    state.primers.forEach((p, i) => {
      if ((typeof p.anneal === 'string' && p.anneal) || p.annealLen) return;   // already known
      if (p.noOverhang) return;                                                // said so by hand
      let r = null;
      try { r = B.findOverhang(p.seq, state.records, {}); } catch (e) { r = null; }
      if (r) map[i] = r;
    });
    state.overhangCache = { key: key, map: map };
    return map;
  };
  /* why this number and not the other one */
  function tmTitle(T) {
    if (!T.tailed) return 'Nearest-neighbour Tm of the whole primer. No overhang is known for it — ' +
      'click to set one if it has one.';
    const where = {
      design: 'as recorded when this primer was designed',
      'by hand': 'as you set it',
      construct: 'from where it anneals on a loaded construct',
      arm: 'because its first ' + T.tailLen + ' bases are sequence of ' + (T.foundOn || 'a loaded construct') +
           ', which is what a homology arm is'
    }[T.basis] || 'from the whole primer';
    return 'The 3\u2032 ' + T.annealLen + ' bases are what binds; the ' + T.tailLen +
      ' base overhang does not, in the first cycles, so set the annealing temperature from this number. ' +
      'The whole primer melts at ' + (isFinite(T.tmFull) ? T.tmFull.toFixed(1) : '?') +
      ' \u00b0C, which is what later cycles see once the overhang has been copied into the product. ' +
      'The split is ' + where +
      (T.basis === 'arm'
        ? '. Where the arm ends can be a base or two out when the insert happens to carry on matching the vector; ' +
          'that errs towards a lower temperature, never a higher one. Click to set it exactly.'
        : '. Click to change it.');
  }
  PB.primerTmOf = function (p, rec, hits, idx) {
    const cond = PB.conditions();
    /* 1 — recorded when the primer was designed: exact, and true whatever
       is loaded */
    if (p && typeof p.anneal === 'string' && p.anneal)
      return B.primerTm(p.seq, { anneal: p.anneal, conditions: cond, basis: 'design' });
    /* 2 — set by hand on this row */
    if (p && p.annealLen)
      return B.primerTm(p.seq, { annealLen: p.annealLen, conditions: cond, basis: 'by hand' });
    if (p && p.noOverhang) return B.primerTm(p.seq, { conditions: cond });
    /* 3 — where it anneals on the construct in front of you */
    const mine = (hits || (rec ? PB.primerHits(rec) : [])).filter(h => h.name === p.name && h.seq === p.seq);
    const best = mine.slice().sort((a, b) => b.annealed - a.annealed)[0];
    if (best && best.tail > 0)
      return B.primerTm(p.seq, { annealLen: best.annealed, conditions: cond, basis: 'construct' });
    /* 4 — its 5' end is sequence of one of the loaded constructs, which is
       what an assembly arm looks like against the vector it was built for */
    const scan = PB.overhangScan();
    const found = idx != null ? scan[idx] : null;
    if (found) {
      const t = B.primerTm(p.seq, { annealLen: found.annealLen, conditions: cond,
        basis: found.how === 'anneals' ? 'construct' : 'arm' });
      t.foundOn = found.on; t.how = found.how; t.disagrees = found.disagrees;
      return t;
    }
    return B.primerTm(p.seq, { conditions: cond });
  };
  PB.primerHits = function (rec) {
    if (!rec || !rec.seq.length || !state.primers.length) return [];
    const key = rec.id + ':' + rec.updatedAt + ':' + rec.seq.length + ':' + rec.topology + ':' +
      state.primers.map(p => p.name + '=' + p.seq).join('|');
    if (state.primerHitCache && state.primerHitCache.key === key) return state.primerHitCache.hits;
    const hits = [];
    state.primers.forEach((p, i) => {
      let sites = [];
      try { sites = B.primerSites(rec.seq, p.seq, { circular: rec.topology === 'circular' }); } catch (e) { sites = []; }
      sites.forEach(h => hits.push(Object.assign({ primer: i, name: p.name, seq: p.seq }, h)));
    });
    state.primerHitCache = { key: key, hits: hits };
    return hits;
  };
  /* which of the stored primers are drawn: all of them unless a choice has
     been made, and then only those */
  PB.primerShown = function (name) {
    const sel = state.primersShown;
    return !sel || sel.indexOf(name) >= 0;
  };
  PB.primerFeatures = function (rec) {
    const real = (rec.features || []).filter(f => /primer/i.test(f.type || ''));
    return PB.primerHits(rec).filter(h => PB.primerShown(h.name) &&
      !real.some(f => f.start === h.start && f.end === h.end &&
      (f.strand || 1) === h.strand)).map(h => ({
        start: h.start, end: h.end, strand: h.strand, type: 'primer_bind',
        label: h.name + (h.mismatches.length ? ' (' + h.mismatches.length + ' mismatch' + (h.mismatches.length === 1 ? '' : 'es') + ')' : ''),
        color: h.mismatches.length ? '#b45309' : '#c2410c', quals: {}, note: 'library primer',
        virtual: true, hit: h
      }));
  };
  /* the construct as drawn: its own features, plus library primers when shown */
  PB.displayRecord = function (rec) {
    if (!rec) return rec;
    let feats = rec.features;
    if (state.showFilePrimers === false) {
      const hid = feats.filter(f => !/^primer_bind$/i.test(f.type || ''));
      if (hid.length !== feats.length) feats = hid;
    }
    const extra = state.showPrimers === false ? [] : PB.primerFeatures(rec);
    if (feats === rec.features && !extra.length) return rec;
    return Object.assign({}, rec, { features: feats.concat(extra) });
  };
  PB.primerSiteText = function (rec, h) {
    const L = rec.seq.length;
    const span = h.wraps ? (h.start + 1) + '–' + L + ' and 1–' + h.end : (h.start + 1) + '–' + h.end;
    return h.name + ' anneals ' + (h.strand === 1 ? 'forward' : 'reverse') + ' at ' + span +
      ' (' + h.annealed + ' nt' + (h.tail ? ', ' + h.tail + ' nt 5′ tail' : '') + ')' +
      (h.mismatches.length ? '; mismatch at primer base ' + h.mismatches.join(', ') : '; perfect match');
  };
  /* a click on a drawn library primer selects where it anneals */
  PB.selectPrimerFeature = function (rec, f) {
    state.selectedFeature = null;
    state.selection = { start: f.start, end: f.end };
    PB.toast(PB.primerSiteText(rec, f.hit), f.hit.mismatches.length ? 'warn' : 'ok');
  };

  function wirePrimers(host, rec) {
    /* Setting the overhang by hand. Needed for a primer that was ordered
       elsewhere, or when nothing that carries its arm is loaded, or to
       correct a boundary the arm search put a base out. */
    host.querySelectorAll('[data-ptm]').forEach(b => b.addEventListener('click', () => {
      const i = +b.getAttribute('data-ptm');
      const p = state.primers[i];
      if (!p) return;
      const T = PB.primerTmOf(p, rec, null, i);
      const cond = PB.conditions();
      const line = n => {
        const t = B.primerTm(p.seq, { annealLen: n, conditions: cond });
        return isFinite(t.tm) ? t.tm.toFixed(1) : '—';
      };
      PB.modal({
        title: 'Overhang on ' + p.name,
        body: '<div class="field wide"><label>The primer</label>' +
            '<div class="seqblock" id="poSeq">' + esc(p.seq) + '</div></div>' +
          '<div class="row"><div class="field narrow"><label>Overhang (5\u2032 bases that do not bind)</label>' +
            '<input class="inp mono" id="poTail" type="number" min="0" max="' + (p.seq.length - 8) +
            '" value="' + (T.tailed ? T.tailLen : 0) + '"></div>' +
          '<div class="field"><label>What that leaves</label><div class="inp mono" id="poOut" ' +
            'style="background:var(--surface-2)">' + (T.tailed ? T.annealLen : p.seq.length) +
            ' nt anneals \u00b7 Tm ' + line(T.tailed ? T.annealLen : p.seq.length) + ' \u00b0C</div></div></div>' +
          '<div class="hint" id="poWhy">' + esc(T.tailed
            ? 'Currently ' + (T.basis === 'design' ? 'as recorded when it was designed'
              : T.basis === 'by hand' ? 'as you set it'
              : T.basis === 'arm' ? 'worked out from ' + (T.foundOn || 'a loaded construct') +
                 ', whose sequence the first ' + T.tailLen + ' bases are'
              : 'worked out from where it anneals on the construct') + '.'
            : 'No overhang is known for this primer, so the whole of it is being measured.') +
          ' Set it to 0 to say it has none; the whole primer is then what is reported.</div>',
        actions: [{ label: 'Cancel' }, { label: 'Save', primary: true, onClick: () => {
          const n = Math.max(0, Math.min(p.seq.length - 8, parseInt(document.querySelector('#poTail').value, 10) || 0));
          delete p.anneal; delete p.tail;
          if (n > 0) { p.annealLen = p.seq.length - n; delete p.noOverhang; }
          else { delete p.annealLen; p.noOverhang = true; }
          PB.lsSet('primers', state.primers);
          state.overhangCache = null;
          PB.closeModal(); PB.renderPrimers();
          PB.toast(n ? (p.name + ': ' + n + ' nt overhang, anneals over ' + (p.seq.length - n))
                     : (p.name + ': no overhang'), 'ok');
          return false;
        } }],
        onOpen: () => {
          const t = document.querySelector('#poTail'), out = document.querySelector('#poOut'),
                sq = document.querySelector('#poSeq');
          const upd = () => {
            const n = Math.max(0, Math.min(p.seq.length - 8, parseInt(t.value, 10) || 0));
            out.textContent = (p.seq.length - n) + ' nt anneals \u00b7 Tm ' + line(p.seq.length - n) + ' \u00b0C';
            sq.innerHTML = n ? '<span style="color:var(--accent)">' + esc(p.seq.slice(0, n)) + '</span>' +
              esc(p.seq.slice(n)) : esc(p.seq);
          };
          if (t) { t.addEventListener('input', upd); upd(); }
        }
      });
    }));
    host.querySelectorAll('[data-psite]').forEach(b => b.addEventListener('click', () => {
      const h = PB.primerHits(rec)[+b.getAttribute('data-psite')];
      if (!h) return;
      state.selection = { start: h.start, end: h.end }; state.selectedFeature = null;
      PB.toast(PB.primerSiteText(rec, h), h.mismatches.length ? 'warn' : 'ok');
      PB.renderHeader(); PB.renderInspector();
    }));
    const read = () => {
      const s = Math.max(1, parseInt($('#pdStart').value, 10) || 1) - 1;
      const e = Math.min(rec.seq.length, parseInt($('#pdEnd').value, 10) || rec.seq.length);
      return {
        start: s, end: e,
        targetTm: parseFloat($('#pdTm').value) || 62,
        minLen: parseInt($('#pdMin').value, 10) || 18,
        maxLen: parseInt($('#pdMax').value, 10) || 30,
        fTail: $('#pdFTail').value, rTail: $('#pdRTail').value
      };
    };
    const useSel = $('#pdUseSel');
    if (useSel) useSel.addEventListener('click', () => {
      state.primerDesign = Object.assign({}, state.primerDesign, {
        start: state.selection.start, end: state.selection.end, result: null });
      renderPrimers();
    });
    $('#pdGo').addEventListener('click', () => {
      const o = read();
      if (o.end - o.start < 40) { PB.toast('Pick a target region of at least 40 bp.', 'warn'); return; }
      const result = B.designPrimers(rec.seq, {
        start: o.start, end: o.end, targetTm: o.targetTm,
        minLen: o.minLen, maxLen: o.maxLen,
        forwardTail: o.fTail, reverseTail: o.rTail,
        conditions: PB.conditions()
      });
      state.primerDesign = Object.assign({}, o, { result });
      renderPrimers();
    });
    const save = $('#pdSave');
    if (save) save.addEventListener('click', () => {
      const r = state.primerDesign.result;
      const base = rec.name.replace(/[^\w.-]+/g, '_').slice(0, 18);
      state.primers.push({ name: base + '_F', seq: r.forward.seq });
      state.primers.push({ name: base + '_R', seq: r.reverse.seq });
      PB.lsSet('primers', state.primers);
      PB.addHistory(rec, 'Primers designed', r.forward.seq + ' / ' + r.reverse.seq);
      PB.touch(rec);
      renderPrimers();
      PB.toast('Saved 2 primers', 'ok');
    });
    const cp = $('#pdCopy');
    if (cp) cp.addEventListener('click', () => {
      const r = state.primerDesign.result;
      PB.copyText('Forward\t' + r.forward.seq + '\t' + f1(r.forward.tm) + ' C\nReverse\t' +
        r.reverse.seq + '\t' + f1(r.reverse.tm) + ' C', 'Primer pair');
    });
    const amp = $('#pdAmp');
    if (amp) amp.addEventListener('click', () => {
      const r = state.primerDesign.result;
      const nr = PB.makeRecord({
        name: rec.name + ' amplicon (' + fmt(r.amplicon.length) + ' bp)',
        seq: r.amplicon.seq, topology: 'linear'
      });
      PB.addHistory(nr, 'PCR product', 'Amplified from ' + rec.name + ' with the designed primer pair');
      PB.addRecord(nr);
      PB.toast('Amplicon saved', 'ok');
    });
    const padd = $('#pAdd');
    if (padd) padd.addEventListener('click', () => {
      PB.modal({
        title: 'Add a primer',
        body: '<div class="row"><div class="field"><label>Name</label><input class="inp" id="npName" placeholder="e.g. seq_fwd_1"></div></div>' +
          '<div class="row"><div class="field wide"><label>Sequence (5′→3′)</label>' +
          '<textarea class="inp" id="npSeq" rows="3" spellcheck="false"></textarea></div></div>',
        actions: [{ label: 'Cancel' }, { label: 'Add', primary: true, onClick: () => {
          const s = B.cleanSeq($('#npSeq').value);
          if (s.length < 8) { PB.toast('That is too short to be a primer.', 'bad'); return false; }
          state.primers.push({ name: ($('#npName').value || 'primer ' + (state.primers.length + 1)).trim(), seq: s });
          PB.lsSet('primers', state.primers);
          renderPrimers();
        } }]
      });
    });
    host.querySelectorAll('[data-pcopy]').forEach(b => b.addEventListener('click', () =>
      PB.copyText(state.primers[parseInt(b.getAttribute('data-pcopy'), 10)].seq, 'Primer')));
    const pex = host.querySelector('#pExport');
    if (pex) pex.addEventListener('click', () => PB.primerExportDialog());
    host.querySelectorAll('[data-pdel]').forEach(b => b.addEventListener('click', () => {
      state.primers.splice(parseInt(b.getAttribute('data-pdel'), 10), 1);
      PB.lsSet('primers', state.primers); renderPrimers();
    }));
    /* renaming in place: a primer ordered under one name should be called
       that here too. Typed characters are stored as they are typed, and the
       panel is only redrawn on leaving the box, so the cursor stays put. */
    host.querySelectorAll('[data-prename]').forEach(el => {
      const i = parseInt(el.getAttribute('data-prename'), 10);
      el.addEventListener('input', () => {
        if (!state.primers[i]) return;
        state.primers[i].name = el.value;
        PB.lsSet('primers', state.primers);
      });
      el.addEventListener('change', () => {
        if (!state.primers[i]) return;
        const v = el.value.trim();
        state.primers[i].name = v || 'primer ' + (i + 1);
        PB.lsSet('primers', state.primers);
        PB.renderAll ? PB.renderAll() : renderPrimers();
      });
      el.addEventListener('keydown', e => { if (e.key === 'Enter') el.blur(); });
    });
  }

  PB.primerFromSelection = function () {
    const rec = PB.active(), sel = state.selection;
    if (!rec || !sel) return;
    const s = PB.selectionSeq(rec, sel);
    if (!s) return;
    const n = PB.normSel(rec, sel);
    state.primers.push({ name: 'sel_' + (n.start + 1) + '_' + ((n.end - 1) % rec.seq.length + 1), seq: s });
    PB.lsSet('primers', state.primers);
    PB.toast('Added ' + s.length + ' nt primer to the library', 'ok');
    if (state.view === 'primers') renderPrimers();
  };

  /* ============================== GEL ============================== */
  function defaultLanes() {
    const rec = PB.active();
    return [
      { kind: 'ladder', ladder: 'NEB 1 kb Plus' },
      { kind: 'digest', recordId: rec ? rec.id : null, enzymes: state.enzymes.slice(0, 2) }
    ];
  }

  function gelBands(sizes, gelPct) {
    const total = sizes.reduce((a, b) => a + b, 0) || 1;
    const maxLen = Math.max.apply(null, sizes.concat([1]));
    const bands = sizes.map(bp => ({
      bp, migration: B.gelMigration(bp, gelPct), mass: bp / total, count: 1
    })).sort((a, b) => a.migration - b.migration);
    const merged = [];
    for (const b of bands) {
      const last = merged[merged.length - 1];
      if (last && Math.abs(last.migration - b.migration) < 0.006) {
        last.mass += b.mass; last.count++;
        last.bp = Math.round((last.bp * (last.count - 1) + b.bp) / last.count);
      } else merged.push(Object.assign({}, b));
    }
    merged.forEach(b => { b.intensity = Math.max(0.3, Math.min(1, (b.bp / maxLen) * 1.25)); });
    return merged;
  }

  /* A primer named from the library, or a sequence typed straight in. */
  function resolvePrimer(txt) {
    const t = String(txt || '').trim();
    if (!t) return null;
    const lib = (state.primers || []).filter(p => (p.name || '').toLowerCase() === t.toLowerCase())[0];
    if (lib) {
      const site = (function () {
        try { return B.primerSites(PB.active().seq, lib.seq, { circular: true })[0]; } catch (e) { return null; }
      })();
      return { name: lib.name, seq: lib.seq, tail: site ? site.tail : 0 };
    }
    const bare = t.replace(/[^ACGTacgtNn]/g, '').toUpperCase();
    if (bare.length >= 10) return { name: bare.length + ' nt typed', seq: bare, tail: 0 };
    return null;
  }

  function renderGel() {
    const host = $('#view-gel');
    const rec = PB.active();
    if (!state.gelLanes) state.gelLanes = defaultLanes();
    const gelPct = state.gelPct || 1.0;
    const lanes = state.gelLanes;

    const laneData = lanes.map(l => {
      if (l.kind === 'ladder') {
        const sizes = B.LADDERS[l.ladder] || B.LADDERS['NEB 1 kb Plus'];
        return { label: l.ladder, ladder: true,
                 bands: sizes.map(bp => ({ bp, migration: B.gelMigration(bp, gelPct), intensity: 0.8 })) };
      }
      const r = state.records.find(x => x.id === l.recordId) || rec;
      if (!r || !r.seq.length) return { label: 'empty', bands: [] };

      /* A PCR lane runs the amplification rather than a digest: what you
         want from a gel before you pour one is not only the band you are
         after but the ones you are not, and those come from the same
         engine that designs the primers. */
      if (l.kind === 'pcr') {
        const f = resolvePrimer(l.fwd), rv = resolvePrimer(l.rev);
        const lab = shortName(r.name) + ' / PCR ' +
          (f ? f.name : '?') + ' + ' + (rv ? rv.name : '?');
        if (!f || !rv)
          return { label: lab, bands: [], sizes: [], record: r, pcr: true,
                   notes: ['Name a forward and a reverse primer, from the library or typed as a sequence.'] };
        let prods = [];
        try {
          prods = B.pcr(r.seq, f.seq, rv.seq, { circular: r.topology === 'circular' });
        } catch (e) { prods = []; }
        const sizes = prods.map(p => p.length).sort((a, b) => b - a);
        const notes = [];
        if (!prods.length)
          notes.push('No product. These two primers do not both prime on this template in the right orientation \u2014 check the template, or that the reverse primer is written 5\u2032\u21923\u2032 on the bottom strand.');
        else if (prods.length > 1)
          notes.push(prods.length + ' products. The band you want is usually the one you designed for; the others come from ' +
            (prods.sites ? prods.sites.forward + ' forward and ' + prods.sites.reverse + ' reverse priming sites' : 'additional priming sites') +
            ' on this template. Gel-purify, or move a primer.');
        if (prods.truncated)
          notes.push('This template is repetitive enough that the primers prime in very many places; only the ' +
            prods.length + ' shortest products are drawn, out of ' + fmt(prods.found) + '.');
        const tail = Math.max(f.tail || 0, rv.tail || 0);
        if (prods.length === 1 && tail)
          notes.push('The product includes the 5\u2032 tails that do not anneal, so it runs ' + tail +
            ' bp longer than the stretch of template between the primers.');
        const rg0 = B.GEL_RANGES[gelPct] || B.GEL_RANGES[1.0];
        const tooSmall = sizes.filter(x => x < rg0[0]), tooBig = sizes.filter(x => x > rg0[1]);
        if (tooSmall.length) notes.push(tooSmall.map(fmt).join(', ') + ' bp is below this gel\u2019s resolving range \u2014 use a higher-percentage gel.');
        if (tooBig.length) notes.push(tooBig.map(fmt).join(', ') + ' bp is above it \u2014 bands compress near the well.');
        return { label: lab, bands: gelBands(sizes, gelPct), sizes: sizes, record: r,
                 notes: notes, pcr: true, products: prods };
      }

      const enzymes = (l.enzymes || []).filter(n => PB.enzByName(n));
      let sizes, notes = [], uncutCircle = false;
      const circ = r.topology === 'circular';
      if (!enzymes.length) {
        sizes = [r.seq.length]; uncutCircle = circ;
      } else {
        const d = B.digest(r.seq, enzymes, { circular: circ, methylation: PB.methylSystems(r) });
        sizes = d.fragments.map(f => f.length).sort((a, b) => b - a);
        uncutCircle = circ && !d.cuts.length;
        const sysName = k => ({ dam: 'Dam', dcm: 'Dcm', cpg: 'CpG' })[k];
        if (d.blocked.length) notes.push(d.blocked.length + ' site' + (d.blocked.length === 1 ? '' : 's') + ' blocked by methylation (' +
          d.blocked.map(b => b.enzyme + ' at ' + fmt(b.pos + 1) + ', ' + b.by.map(sysName).join('/')).slice(0, 4).join('; ') + ') — not cut');
        const maybe = d.caution.filter(c => c.effect !== 'unknown');
        if (maybe.length) notes.push(maybe.length + ' site' + (maybe.length === 1 ? '' : 's') + ' may cut poorly (' +
          maybe.map(b => b.enzyme + ' at ' + fmt(b.pos + 1)).slice(0, 4).join('; ') + ') — partial-digest bands possible');
      }
      const rg = B.GEL_RANGES[gelPct] || B.GEL_RANGES[1.0];
      const small = sizes.filter(x => x < rg[0]), big = sizes.filter(x => x > rg[1]);
      if (small.length) notes.push(small.map(fmt).join(', ') + ' bp below this gel’s resolving range (' + fmt(rg[0]) + ' bp): runs near the dye front and may be faint or invisible — use a higher-percentage gel to see it');
      if (big.length) notes.push(big.map(fmt).join(', ') + ' bp above the resolving range: bands compress near the well');
      if (uncutCircle) notes.push('uncut plasmid: supercoiled DNA runs faster than linear DNA of the same size, and nicked DNA slower, so its bands do not read against a linear ladder; the band is drawn at its true length, outlined');
      return {
        label: shortName(r.name) + (enzymes.length ? ' / ' + enzymes.join('+') : ' uncut'),
        bands: gelBands(sizes, gelPct), sizes: sizes, record: r, enzymes, notes, uncutCircle
      };
    });

    const laneW = 64, padX = 14, gelH = 340, topPad = 26;
    const W = padX * 2 + laneData.length * laneW + 46;
    const H = topPad + gelH + 16;
    let svg = '<svg class="gel" viewBox="0 0 ' + W + ' ' + H + '" width="' + W + '" role="img" aria-label="Simulated agarose gel">';
    svg += '<rect x="0" y="0" width="' + W + '" height="' + H + '" fill="#0d0f10" rx="6"/>';
    laneData.forEach((ld, i) => {
      const x = padX + 46 + i * laneW;
      svg += '<rect x="' + (x - laneW / 2 + 6) + '" y="' + topPad + '" width="' + (laneW - 12) + '" height="' + gelH + '" fill="#14181a"/>';
      svg += '<rect x="' + (x - laneW / 2 + 6) + '" y="' + (topPad - 4) + '" width="' + (laneW - 12) + '" height="4" fill="#2a3135"/>';
      const lbl = ld.label.length > 13 ? ld.label.slice(0, 12) + '…' : ld.label;
      svg += '<text class="lanelbl" x="' + x + '" y="' + (topPad - 9) + '" text-anchor="middle">' + esc(lbl) + '</text>';
      let lastLblY = -99;
      const rgl = B.GEL_RANGES[gelPct] || B.GEL_RANGES[1.0];
      (ld.bands || []).forEach(b => {
        const y = topPad + b.migration * gelH;
        const outside = b.bp < rgl[0] || b.bp > rgl[1];
        const showLbl = Math.abs(y - lastLblY) >= 9;
        if (showLbl) lastLblY = y;
        const h = ld.ladder ? 3 : Math.max(2.4, 3.6 * b.intensity + 0.8);
        const op = (ld.ladder ? 0.85 : (0.35 + 0.6 * b.intensity)) * (outside ? 0.35 : 1);
        svg += ld.uncutCircle
          ? '<rect x="' + (x - laneW / 2 + 9) + '" y="' + (y - h / 2).toFixed(1) + '" width="' + (laneW - 18) +
            '" height="' + h.toFixed(1) + '" fill="none" stroke="#e8f6f4" stroke-dasharray="3 2" opacity=".75" rx="1"/>'
          : '<rect x="' + (x - laneW / 2 + 9) + '" y="' + (y - h / 2).toFixed(1) + '" width="' + (laneW - 18) +
            '" height="' + h.toFixed(1) + '" fill="#e8f6f4" opacity="' + op.toFixed(2) + '" rx="1"/>';
        if (!showLbl) { /* label would collide with the one above */ }
        else if (ld.ladder) {
          svg += '<text x="' + (x - laneW / 2 + 4) + '" y="' + (y + 3).toFixed(1) + '" text-anchor="end">' + fmt(b.bp) + '</text>';
        } else if (b.count === 1) {
          svg += '<text x="' + (x + laneW / 2 - 8) + '" y="' + (y + 3).toFixed(1) + '" text-anchor="start" opacity=".8">' + fmt(b.bp) + '</text>';
        } else {
          svg += '<text x="' + (x + laneW / 2 - 8) + '" y="' + (y + 3).toFixed(1) + '" text-anchor="start" opacity=".8">' +
            fmt(b.bp) + '×' + b.count + '</text>';
        }
      });
    });
    svg += '</svg>';

    let html = '<div class="panel"><h3>Simulated agarose gel <span class="sub">' +
      gelPct + '% agarose · resolving ' + fmt(B.GEL_RANGES[gelPct][0]) + '–' + fmt(B.GEL_RANGES[gelPct][1]) + ' bp</span></h3>' +
      '<div class="gel-stage">' +
      '<div style="flex:0 0 auto;overflow-x:auto;max-width:100%">' + svg + '</div>' +
      '<div class="gel-lanes">' +
      '<div class="row"><div class="field narrow"><label>Agarose %</label><select class="inp" id="glPct">' +
        Object.keys(B.GEL_RANGES).map(p => '<option value="' + p + '"' + (+p === gelPct ? ' selected' : '') + '>' + p + '%</option>').join('') +
      '</select></div></div>';

    lanes.forEach((l, i) => {
      html += '<div class="lane-cfg"><span class="lname">Lane ' + (i + 1) + '</span>' +
        '<select class="inp" data-lk="' + i + '" style="flex:0 0 96px">' +
          '<option value="ladder"' + (l.kind === 'ladder' ? ' selected' : '') + '>Ladder</option>' +
          '<option value="digest"' + (l.kind === 'digest' ? ' selected' : '') + '>Digest</option>' +
          '<option value="pcr"' + (l.kind === 'pcr' ? ' selected' : '') + '>PCR</option>' +
        '</select>';
      if (l.kind === 'ladder') {
        html += '<select class="inp" data-ll="' + i + '">' +
          Object.keys(B.LADDERS).map(k => '<option' + (k === l.ladder ? ' selected' : '') + '>' + esc(k) + '</option>').join('') +
        '</select>';
      } else if (l.kind === 'pcr') {
        html += '<select class="inp" data-lr="' + i + '" style="flex:1 1 100px">' +
          state.records.map(r => '<option value="' + r.id + '"' + (r.id === l.recordId ? ' selected' : '') + '>' + esc(shortName(r.name)) + '</option>').join('') +
        '</select>' +
        '<input class="inp mono" list="glPrimers" data-lf="' + i + '" value="' + esc(l.fwd || '') +
          '" placeholder="forward" style="flex:1 1 105px">' +
        '<input class="inp mono" list="glPrimers" data-lv="' + i + '" value="' + esc(l.rev || '') +
          '" placeholder="reverse" style="flex:1 1 105px">';
      } else {
        html += '<select class="inp" data-lr="' + i + '" style="flex:1 1 110px">' +
          state.records.map(r => '<option value="' + r.id + '"' + (r.id === l.recordId ? ' selected' : '') + '>' + esc(shortName(r.name)) + '</option>').join('') +
        '</select>' +
        '<input class="inp mono" data-le="' + i + '" value="' + esc((l.enzymes || []).join(' ')) + '" placeholder="EcoRI HindIII" style="flex:1 1 130px">';
      }
      html += '<button class="btn sm ghost" data-lx="' + i + '" title="Remove lane">×</button></div>';
    });
    html += '<datalist id="glPrimers">' +
      (state.primers || []).map(p => '<option value="' + esc(p.name) + '">').join('') + '</datalist>';
    html += '<div class="row" style="gap:5px"><button class="btn sm" id="glAdd">Add lane</button>' +
      '<button class="btn sm" id="glDigestAll">One lane per selected enzyme</button>' +
      '<button class="btn sm" id="glPcrAdd">Add a PCR lane</button>' +
      (state.gbDesign && state.gbDesign.p
        ? '<button class="btn sm" id="glFromDesign">Lane from the primers I just designed</button>' : '') +
      '</div>' +
      '<div class="hint">Band positions follow a logarithmic mobility model over the gel’s resolving window, and band brightness scales with the mass of DNA in each fragment. Fragments that would not resolve are merged and marked ×n.</div>' +
      '</div></div>';

    // fragment table
    const digestLanes = laneData.filter(l => l.sizes);
    if (digestLanes.length) {
      html += '<div class="panel"><h3>' +
        (digestLanes.some(l => l.pcr) ? 'Fragments and products' : 'Fragments') + '</h3><div class="scrollx"><table class="tbl"><thead><tr>' +
        '<th>Lane</th><th class="num">Bands</th><th>Sizes (bp)</th><th>Notes</th></tr></thead><tbody>' +
        digestLanes.map(l => '<tr><td>' + esc(l.label) + '</td><td class="num">' + l.sizes.length + '</td>' +
          '<td class="mono" style="font-size:11.5px">' + l.sizes.map(fmt).join(' · ') + '</td>' +
          '<td style="font-size:11.5px;max-width:340px">' + (l.notes || []).map(n => '<div>' + esc(n) + '</div>').join('') + '</td></tr>').join('') +
        '</tbody></table></div></div>';
    }
    html += '</div>';
    host.innerHTML = html;

    $('#glPct').addEventListener('change', e => { state.gelPct = parseFloat(e.target.value); renderGel(); });
    const pcrAdd = $('#glPcrAdd');
    if (pcrAdd) pcrAdd.addEventListener('click', () => {
      const lib = state.primers || [];
      lanes.push({ kind: 'pcr', recordId: rec ? rec.id : null,
                   fwd: lib[0] ? lib[0].name : '', rev: lib[1] ? lib[1].name : '' });
      renderGel();
    });
    /* straight from the assembly tool: the insert amplification you just
       designed, on the template it will actually be run on */
    const fromDesign = $('#glFromDesign');
    if (fromDesign) fromDesign.addEventListener('click', () => {
      const d = state.gbDesign;
      if (!d || !d.p) return;
      const donor = (d.part && d.part.rec) || d.r || null;
      lanes.push({ kind: 'pcr', recordId: donor ? donor.id : (rec ? rec.id : null),
                   fwd: d.p.forward.seq, rev: d.p.reverse.seq });
      if (d.open && d.open.primers)
        lanes.push({ kind: 'pcr', recordId: d.r ? d.r.id : (rec ? rec.id : null),
                     fwd: d.open.primers.forward.seq, rev: d.open.primers.reverse.seq });
      PB.toast('Added the insert' + (d.open && d.open.primers ? ' and backbone' : '') + ' amplification', 'ok');
      renderGel();
    });
    host.querySelectorAll('[data-lk]').forEach(s => s.addEventListener('change', () => {
      const i = +s.getAttribute('data-lk');
      lanes[i].kind = s.value;
      if (s.value === 'ladder' && !lanes[i].ladder) lanes[i].ladder = 'NEB 1 kb Plus';
      if (s.value === 'digest' && !lanes[i].recordId) {
        lanes[i].recordId = rec ? rec.id : null; lanes[i].enzymes = state.enzymes.slice(0, 2);
      }
      if (s.value === 'pcr') {
        if (!lanes[i].recordId) lanes[i].recordId = rec ? rec.id : null;
        const lib = state.primers || [];
        if (!lanes[i].fwd && lib[0]) lanes[i].fwd = lib[0].name;
        if (!lanes[i].rev && lib[1]) lanes[i].rev = lib[1].name;
      }
      renderGel();
    }));
    host.querySelectorAll('[data-lf]').forEach(inp => inp.addEventListener('change', () => {
      lanes[+inp.getAttribute('data-lf')].fwd = inp.value; renderGel();
    }));
    host.querySelectorAll('[data-lv]').forEach(inp => inp.addEventListener('change', () => {
      lanes[+inp.getAttribute('data-lv')].rev = inp.value; renderGel();
    }));
    host.querySelectorAll('[data-ll]').forEach(s => s.addEventListener('change', () => {
      lanes[+s.getAttribute('data-ll')].ladder = s.value; renderGel();
    }));
    host.querySelectorAll('[data-lr]').forEach(s => s.addEventListener('change', () => {
      lanes[+s.getAttribute('data-lr')].recordId = s.value; renderGel();
    }));
    host.querySelectorAll('[data-le]').forEach(inp => inp.addEventListener('change', () => {
      const names = inp.value.split(/[\s,+]+/).filter(Boolean);
      const bad = names.filter(n => !PB.enzByName(n));
      if (bad.length) PB.toast('Unknown enzyme: ' + bad.join(', '), 'warn');
      lanes[+inp.getAttribute('data-le')].enzymes = names.filter(n => PB.enzByName(n)).map(n => PB.enzByName(n).name);
      renderGel();
    }));
    host.querySelectorAll('[data-lx]').forEach(b => b.addEventListener('click', () => {
      lanes.splice(+b.getAttribute('data-lx'), 1);
      if (!lanes.length) state.gelLanes = defaultLanes();
      renderGel();
    }));
    $('#glAdd').addEventListener('click', () => {
      lanes.push({ kind: 'digest', recordId: rec ? rec.id : null, enzymes: state.enzymes.slice(0, 1) });
      renderGel();
    });
    $('#glDigestAll').addEventListener('click', () => {
      state.gelLanes = [{ kind: 'ladder', ladder: 'NEB 1 kb Plus' }].concat(
        state.enzymes.slice(0, 10).map(n => ({ kind: 'digest', recordId: rec ? rec.id : null, enzymes: [n] })));
      renderGel();
    });
  }
  PB.renderGel = renderGel;

  function shortName(n) { return n.length > 22 ? n.slice(0, 21) + '…' : n; }

  /* ============================== VERIFY (sequencing) ============================== */
  function renderVerify() {
    const host = $('#view-verify');
    const rec = PB.active();
    state.reads = state.reads || [];
    /* A read aligned against a different construct carries numbers that mean
       nothing here, so it is realigned once when the reference changes. */
    if (rec && rec.seq.length) {
      const fp = rec.id + ':' + rec.updatedAt;
      state.reads.forEach(r => {
        if (r.alnFor === fp) return;
        const forAlign = r.trim && r.trim.seq.length > 40 ? r.trim.seq : r.seq;
        if (forAlign.length < 20) { r.aln = null; r.alnFor = fp; return; }
        try { r.aln = B.alignToReference(rec.seq, forAlign); } catch (e) { r.aln = null; }
        r.alnFor = fp;
      });
    }
    let html = '<div class="panel"><h3>Verify by sequencing ' +
      '<span class="sub">align Sanger reads to “' + esc(rec ? rec.name : 'nothing') + '” and list every discrepancy</span></h3>' +
      '<div class="row">' +
      '<button class="btn" id="vAb1">Load reads (.ab1, .seq, FASTA)</button>' +
      '<button class="btn" id="vPaste">Paste reads</button>' +
      (state.reads.length ? '<button class="btn sm" id="vClear">Clear reads</button>' : '') +
      '</div>' +
      '<div class="hint">Pick a whole plate at once — any mix of .ab1 traces and .seq or FASTA files; a multi-record FASTA gives one read per record. ' +
      'Files are read in your browser; nothing is uploaded. From an .ab1, base calls, quality values and the four dye traces are all parsed from the ABIF container. ' +
      'Reads without quality values (.seq, FASTA, pasted) have their N-rich ends set aside before they are judged. ' +
      'The reads can be translated without any construct loaded.</div>';

    if (!rec) { html += '<div class="note warn">No construct is loaded, so the reads cannot be aligned or checked for differences. They can still be translated and triaged below; load the construct to compare against it and to pin frames by its coding sequences.</div>'; }
    if (!state.reads.length) {
      html += emptyState('No reads loaded. Add an .ab1 file from your sequencing provider, or paste the base calls as text.');
    } else {
      html += '<div class="scrollx" style="margin-top:10px"><table class="tbl"><thead><tr>' +
        '<th>Read</th><th class="num">nt</th><th class="num">Mean Q</th><th class="num">Trimmed</th>' +
        '<th class="num">Identity</th><th class="num">Covers</th><th>Result</th><th></th></tr></thead><tbody>' +
        state.reads.map((rd, i) => {
          const a = rd.aln;
          let verdict = '<span class="badge neutral">not aligned</span>';
          if (a) {
            verdict = a.clean
              ? '<span class="badge ok">perfect match</span>'
              : '<span class="badge ' + (a.differences.length <= 2 ? 'warn' : 'bad') + '">' +
                a.differences.length + ' difference' + (a.differences.length === 1 ? '' : 's') + '</span>';
          }
          return '<tr>' +
            '<td>' + esc(rd.name) + (rd.hasTraces ? ' <span class="badge acc">trace</span>' : '') +
              (a && a.flipped ? ' <span class="badge neutral">rev</span>' : '') + '</td>' +
            '<td class="num">' + fmt(rd.seq.length) + '</td>' +
            '<td class="num">' + (rd.meanQuality != null ? Math.round(rd.meanQuality) : '—') + '</td>' +
            '<td class="num">' + (rd.trim ? fmt(rd.trim.seq.length) : '—') + '</td>' +
            '<td class="num">' + (a ? f1(a.identity) + '%' : '—') + '</td>' +
            '<td class="num">' + (a ? f1(a.coverage) + '%' : '—') + '</td>' +
            '<td>' + verdict + '</td>' +
            '<td style="white-space:nowrap">' +
              '<button class="btn sm ghost" data-vdet="' + i + '">Details</button>' +
              '<button class="btn sm ghost" data-vdel="' + i + '" title="Remove">×</button>' +
            '</td></tr>';
        }).join('') + '</tbody></table></div>';
    }
    html += '</div>';

    if (state.reads.length) html += translatePanel();

    // detail panel
    const di = state.verifyDetail;
    if (di != null && state.reads[di]) {
      const rd = state.reads[di];
      html += '<div class="panel"><h3>' + esc(rd.name) + ' <span class="sub">' +
        (rd.aln ? f1(rd.aln.identity) + '% identity over ' + fmt(rd.aln.alignedLength) + ' bp, reference ' +
          fmt(rd.aln.refStart + 1) + '–' + fmt(rd.aln.refEnd) : 'not aligned') + '</span></h3>';
      if (rd.aln) {
        html += coveragePanel(rec, rd);
        if (rd.aln.clean) {
          html += '<div class="note ok"><strong>No discrepancies.</strong> The read matches the reference exactly over the aligned region.</div>';
        } else {
          html += '<div class="scrollx"><table class="tbl"><thead><tr>' +
            '<th class="num">Position</th><th>Type</th><th>Reference</th><th>Read</th>' +
            '<th>In feature</th><th>Effect</th><th></th></tr></thead><tbody>' +
            rd.aln.differences.slice(0, 400).map((d, k) => {
              const ann = annotateDiff(rec, d);
              return '<tr>' +
                '<td class="num">' + fmt(d.refPos + 1) + '</td>' +
                '<td>' + d.type + '</td>' +
                '<td class="mono">' + esc(d.ref || d.base || '') + '</td>' +
                '<td class="mono">' + esc(d.read || (d.type === 'deletion' ? '–' : d.base) || '') + '</td>' +
                '<td>' + (ann.feature ? esc(ann.feature) : '—') + '</td>' +
                '<td>' + ann.effect + '</td>' +
                '<td>' + (rd.hasTraces ? '<button class="btn sm ghost" data-vtr="' + k + '">trace</button>' : '') + '</td>' +
                '</tr>';
            }).join('') + '</tbody></table></div>';
        }
        html += '<div class="row" style="margin-top:10px"><div class="field wide"><label>Alignment</label>' +
          alignBlocks(rd.aln.alignment.a, rd.aln.alignment.b, 'ref', 'read', 60, rd.aln.refStart) + '</div></div>';
      }
      if (rd.quality && rd.quality.length) {
        html += '<div class="row"><div class="field wide"><label>Quality (Phred)</label>' + qualitySVG(rd) + '</div></div>';
      }
      if (state.traceAt != null && rd.hasTraces) {
        html += '<div class="row"><div class="field wide"><label>Chromatogram around the difference</label>' +
          chromaSVG(rd, state.traceAt) + '</div></div>';
      }
      html += translateBlock(rd);
      html += '<div class="row"><button class="btn sm" data-vcopy="' + di + '">Copy read</button>' +
        '<button class="btn sm" data-vsave="' + di + '">Save read as sequence</button>' +
        '<button class="btn sm" data-vnblast="' + di + '">BLAST the read (blastn)</button></div>';
      html += '</div>';
    }

    if (rec && state.reads.length) html += proteinPanel(rec);
    html += blastPanel();

    host.innerHTML = html;
    wireVerify(host, rec);
  }
  PB.renderVerify = renderVerify;

  /* ---------------- translate the loaded reads ----------------
     Straight nucleotide-to-protein on the reads themselves. No reference is
     needed: paste or load reads, get the amino acids, copy the FASTA or send
     a peptide to BLAST. Quality-trimmed bases are used where they exist
     because the ragged ends of a read translate to nonsense.            */
  function framesOf(seq, table) {
    return B.frameTranslate(seq, { table: table }).map(f =>
      Object.assign({}, f, { ntStart: f.longestAt }));
  }

  /* Triage every loaded read: which frame is the protein, and is the read
     worth believing. Cached, because this runs on every render of the tab. */
  function readTriages() {
    const rec = PB.active();
    const table = state.codeTable;
    const key = (rec ? rec.id + ':' + rec.updatedAt : 'none') + ':' + table + ':' +
      state.reads.map(r => r.name + '/' + r.seq.length).join('|');
    if (state.triageCache && state.triageCache.key === key) return state.triageCache.rows;
    let lms = [];
    if (rec) { try { lms = B.landmarksFromRecord(rec, { table: table }); } catch (e) { lms = []; } }
    const rows = [];
    state.reads.forEach((rd, i) => {
      const seq = (rd.trim && rd.trim.seq.length > 40 ? rd.trim.seq : rd.seq) || '';
      if (seq.length < 30) return;
      let t;
      try { t = B.readTriage(seq, { table: table, landmarks: lms }); } catch (e) { return; }
      rows.push({ read: rd.name, idx: i, seq: seq, t: t });
    });
    let dups = [];
    try { dups = B.duplicateReads(rows.map(r => ({ seq: r.seq })), {}); } catch (e) { dups = []; }
    dups.forEach((g, gi) => g.forEach(j => { if (rows[j]) rows[j].dupGroup = gi + 1; }));
    state.triageCache = { key: key, rows: rows, dups: dups, usedConstruct: lms.length > 0 };
    return rows;
  }

  /* The whole triage as a spreadsheet: one row per read, with every signal
     the verdict was made from, so a plate can be reviewed or filed.
     RFC 4180 quoting, so a flag containing a comma stays in its cell. */
  function triageCsv() {
    const cache = state.triageCache;
    if (!cache || !cache.rows.length) return '';
    const q = v => {
      const t = String(v == null ? '' : v).replace(/\r?\n/g, ' ');
      return /[",]/.test(t) ? '"' + t.replace(/"/g, '""') + '"' : t;
    };
    const head = ['read', 'verdict', 'verdict_label', 'frame', 'clean_aa', 'no_stop_aa', 'stops',
      'ambiguous_pct', 'composition_fit', 'dominance', 'pinned_by', 'judged_from', 'judged_to',
      'read_length', 'same_dna_group', 'indel_between', 'flags', 'peptide'];
    const lines = [head.join(',')].concat(cache.rows.map(r => {
      const t = r.t, c = t.chosen;
      const pins = c ? c.landmarks.filter(l => l.pins).map(l => l.name).join('; ') : '';
      return [r.read, t.verdict, B.VERDICTS[t.verdict], c ? c.label.replace('\u2212', '-') : '',
        c ? c.cleanLen : '', c ? c.longest.length : '', c ? c.stops : '',
        (t.ambiguousFrac * 100).toFixed(1), c && c.fit !== null ? c.fit.toFixed(3) : '',
        !c ? '' : isFinite(t.dominance) ? t.dominance.toFixed(2) : 'pinned', pins,
        t.span ? t.span.start + 1 : 1, t.span ? t.span.end : t.length, t.length,
        r.dupGroup || '', t.shift ? t.shift.readFrom + '-' + t.shift.readTo : '',
        t.flags.map(f => f.level + ': ' + f.text).join(' | '),
        c ? (c.clean || c.longest) : ''].map(q).join(',');
    }));
    return lines.join('\r\n') + '\r\n';
  }
  PB.triageCsv = triageCsv;

  const VERDICT_STYLE = { 'clean': 'ok', 'partial': 'acc', 'open-only': 'warn',
                          'frameshift': 'warn', 'frame-ambiguous': 'warn',
                          'low-quality': 'bad', 'no-protein': 'bad', 'too-short': 'neutral' };
  const WORTH = { 'clean': 1, 'partial': 1, 'frameshift': 1 };

  function translatePanel() {
    const mode = state.trMode || 'best';
    const rows = readTriages();
    const cache = state.triageCache || {};
    if (!rows.length) return '';

    /* one row per read, or all six frames when asked */
    const items = [];
    rows.forEach(r => {
      if (mode === 'all') r.t.frames.forEach(f => items.push({ read: r.read, idx: r.idx, f: f, t: r.t }));
      else {
        items.push({ read: r.read, idx: r.idx, f: r.t.chosen, t: r.t, lead: true });
        /* when no frame stands out, the runner-up is shown too rather than
           quietly dropped — the whole point is not to pick one and hope */
        if (r.t.verdict === 'frame-ambiguous' && r.t.runnerUp)
          items.push({ read: r.read, idx: r.idx, f: r.t.runnerUp, t: r.t, alt: true });
      }
    });
    state.trItems = items;

    const tally = {};
    rows.forEach(r => { tally[r.t.verdict] = (tally[r.t.verdict] || 0) + 1; });
    const order = ['clean', 'partial', 'frameshift', 'open-only', 'frame-ambiguous', 'low-quality', 'no-protein', 'too-short'];
    const summary = order.filter(v => tally[v]).map(v =>
      '<span class="badge ' + VERDICT_STYLE[v] + '" title="' + esc(B.VERDICTS[v]) + '">' +
      tally[v] + ' \u00d7 ' + esc(B.VERDICTS_SHORT[v]) + '</span>').join(' ');
    const worth = rows.filter(r => WORTH[r.t.verdict]).length;

    let dupNote = '';
    if (cache.dups && cache.dups.length) {
      /* A group whose members did not all get the same verdict is the one
         worth stopping at. The same DNA cannot be both clean and frameshifted,
         so either the population is mixed, the polymerase made the difference,
         or one of the traces is not to be trusted — and a reader who takes the
         group as "duplicates, ignore the rest" would throw away the one read
         that found the problem. */
      const split = cache.dups.map((g, i) => {
        const vs = Array.from(new Set(g.map(j => rows[j] && rows[j].t ? rows[j].t.verdict : '?')));
        return vs.length > 1 ? { i: i + 1, names: g.map(j => rows[j] ? rows[j].read : '?'), verdicts: vs } : null;
      }).filter(Boolean);
      dupNote = '<div class="note warn"><strong>Some reads are the same DNA.</strong> ' +
        cache.dups.map((g, i) => 'Group ' + (i + 1) + ': ' +
          g.map(j => esc(rows[j] ? rows[j].read : '?')).join(', ')).join(' · ') +
        '. That is expected when one clone was read with two primers or read twice. ' +
        'It is not expected when the names say different clones — then one template was sequenced ' +
        'several times, the wells were swapped, or one trace came back under several names, and ' +
        'those reads are not independent evidence. Check the traces before concluding anything ' +
        'about those clones.</div>' +
        split.map(d => '<div class="note bad" data-dupsplit="' + d.i + '"><strong>Group ' + d.i +
          ' reads the same DNA but did not get the same verdict</strong> — ' +
          esc(d.names.join(', ')) + ' came back as ' + esc(d.verdicts.join(' and ')) + '. ' +
          'One template cannot be both, so either the colony was mixed, the difference was ' +
          'introduced during amplification, or one of these traces is wrong. Do not treat either ' +
          'read as confirming the other: look at both traces at the position where they differ.' +
          '</div>').join('');
    }

    const flagRow = (t, span) => t.flags.length
      ? '<tr class="flagrow"><td colspan="' + span + '" style="padding-top:0">' +
        t.flags.map(f => '<div class="note ' + (f.level === 'ok' ? 'ok' : f.level === 'bad' ? 'bad' : f.level === 'warn' ? 'warn' : '') +
          '" style="margin:2px 0;font-size:11.5px">' + esc(f.text) + '</div>').join('') + '</td></tr>'
      : '';

    const filt = state.trOnlyOk;
    return '<div class="panel"><h3>Translate the reads ' +
      '<span class="sub">nucleotides to amino acids, with the frame called for you</span></h3>' +
      dupNote +
      '<div class="row" style="margin-bottom:6px">' + summary + '</div>' +
      '<div class="chipbar">' +
      [['best', 'Called frame per read'], ['all', 'All six frames']].map(([v, l]) =>
        '<button class="chip' + (mode === v ? ' on' : '') + '" data-trm="' + v + '">' + l + '</button>').join('') +
      '<button class="chip' + (filt ? ' on' : '') + '" id="trOnlyOk">Only reads worth BLASTing (' + worth + ')</button>' +
      '</div>' +
      '<div class="scrollx scrolly"><table class="tbl" id="trTable"><thead><tr><th>Read</th><th>Verdict</th><th>Frame</th>' +
      '<th class="num">Peptide</th><th class="num">Stops</th><th class="num">N</th><th class="num">Fit</th>' +
      '<th>Amino acids</th><th></th></tr></thead><tbody>' +
      items.map((x, k) => {
        if (filt && !WORTH[x.t.verdict]) return '';
        const pep = x.f.clean || x.f.longest;
        const v = x.t.verdict;
        const main = '<tr' + (x.alt ? ' style="opacity:.75"' : '') + '>' +
          '<td>' + (x.alt ? '' : esc(x.read) +
            (x.t.dupGroup ? ' <span class="badge bad" title="same DNA as another read">dup ' + x.t.dupGroup + '</span>' : '')) + '</td>' +
          '<td>' + (x.alt ? '<span class="hint" style="margin:0">runner-up</span>'
            : '<span class="badge ' + VERDICT_STYLE[v] + '">' + esc(B.VERDICTS[v]) + '</span>') + '</td>' +
          '<td class="mono">' + x.f.label + (x.f.pinned ? ' <span class="badge ok" title="' +
            esc(x.f.landmarks.filter(l => l.pins).map(l => l.name).join(', ')) + '">pinned</span>' : '') + '</td>' +
          '<td class="num"><b>' + pep.length + '</b> aa</td>' +
          '<td class="num">' + x.f.stops + '</td>' +
          '<td class="num">' + (x.t.ambiguousFrac * 100 < 0.05 ? '—'
            : (x.t.ambiguousFrac * 100).toFixed(1) + '%') + '</td>' +
          '<td class="num" title="chi-square per residue against the Swiss-Prot average; below 0.3 is protein-like">' +
            (x.f.fit === null ? '—' : x.f.fit.toFixed(2)) + '</td>' +
          '<td class="seqcell mono" style="max-width:260px;font-size:11.5px">' +
            esc(pep.slice(0, 44)) + (pep.length > 44 ? '…' : '') + '</td>' +
          '<td style="white-space:nowrap">' +
            '<button class="btn sm ghost" data-trc="' + k + '" title="the peptide, as amino acids">Copy aa</button>' +
            '<button class="btn sm ghost" data-trnt="' + k + '" title="the bases this peptide was translated from, ' +
              'on the strand of frame ' + x.f.label + '">Copy nt</button>' +
            '<button class="btn sm ghost" data-trb="' + k + '">BLAST</button>' +
            (x.f.fromMet && x.f.fromMet.length !== pep.length
              ? '<button class="btn sm ghost" data-trm2="' + k + '">from Met</button>' : '') +
          '</td></tr>';
        return main + (mode === 'best' && x.lead ? flagRow(x.t, 9) : '');
      }).join('') +
      '</tbody></table></div>' +
      '<div class="row" style="margin-top:8px">' +
        '<button class="btn sm" id="trFasta">Copy FASTA</button>' +
        '<button class="btn sm" id="trDl">Download FASTA</button>' +
        '<button class="btn sm" id="trCsv">Download triage table (.csv)</button>' +
        '<span class="hint" style="margin:0;align-self:center">' + rows.length + ' read' +
          (rows.length === 1 ? '' : 's') + ', ' + worth + ' worth sending to BLAST' +
          (cache.usedConstruct ? '; frames also checked against the loaded construct' : '') + '</span>' +
      '</div>' +
      '<div class="hint">The frame is called from three things at once: the longest stretch with no stop ' +
      '<em>and no ambiguous residue</em>, how far the amino-acid composition sits from an ordinary protein, ' +
      'and whether a tag, protease site or a coding sequence of the loaded construct matches in that frame. ' +
      'Stops alone are not enough — an N base translates to X rather than to a stop, so a dead trace reports ' +
      'a deceptively low stop count, which is why the N column sits next to it. ' +
      'Load the construct to pin frames that carry no tag. ' +
      'BLAST copies the peptide and opens NCBI blastp; paste it if the box comes up empty. ' +
      '<em>Copy nt</em> gives the bases that peptide was translated from, on the strand of its own frame \u2014 ' +
      'so a minus frame hands back the reverse complement of the read, which is what reads 5\u2032\u21923\u2032 ' +
      'through the gene. They stop where the peptide does, so the stop codon that ended it is not included.</div></div>';
  }

  /* ---------------- translate a read, six frames ----------------
     What you BLAST is usually the longest peptide a read can produce, so
     each frame is translated and cut at its stop codons, and the longest
     piece of each is offered. Quality-trimmed bases are used where they
     exist, because the ragged ends of a read translate to nonsense.     */
  function translateBlock(rd) {
    const seq = (rd.trim && rd.trim.seq.length > 40 ? rd.trim.seq : rd.seq) || '';
    if (seq.length < 60) return '';
    const table = state.codeTable;
    const rec = PB.active();
    let lms = [];
    if (rec) { try { lms = B.landmarksFromRecord(rec, { table: table }); } catch (e) { lms = []; } }
    let tri;
    try { tri = B.readTriage(seq, { table: table, landmarks: lms }); } catch (e) { return ''; }
    const rows = tri.frames;
    state.verifyFrames = rows;
    return '<div class="row" style="margin-top:12px"><div class="field wide">' +
      '<label>Translate this read — all six frames, best call first</label>' +
      '<div class="scrollx"><table class="tbl"><thead><tr><th>Frame</th><th class="num">Clean run</th>' +
      '<th class="num">No stops</th><th class="num">From the first Met</th><th class="num">Stops</th>' +
      '<th class="num">Fit</th><th>Peptide</th><th></th></tr></thead><tbody>' +
      rows.map((r, i) => '<tr' + (i === 0 ? ' style="font-weight:500"' : '') + '>' +
        '<td class="mono">' + r.label +
          (r.pinned ? ' <span class="badge ok" title="' +
            esc(r.landmarks.filter(l => l.pins).map(l => l.name).join(', ')) + '">pinned</span>' : '') + '</td>' +
        '<td class="num"><b>' + r.cleanLen + '</b> aa</td>' +
        '<td class="num">' + r.longest.length + ' aa</td>' +
        '<td class="num">' + (r.fromMet ? r.fromMet.length + ' aa' : '—') + '</td>' +
        '<td class="num">' + r.stops + '</td>' +
        '<td class="num">' + (r.fit === null ? '—' : r.fit.toFixed(2)) + '</td>' +
        '<td class="seqcell mono" style="max-width:300px;font-size:11.5px">' +
          esc((r.clean || r.longest).slice(0, 50)) + ((r.clean || r.longest).length > 50 ? '…' : '') + '</td>' +
        '<td style="white-space:nowrap">' +
          '<button class="btn sm ghost" data-fcopy="' + i + '">Copy aa</button>' +
          '<button class="btn sm ghost" data-fnt="' + i + '" title="the bases this peptide was translated from, ' +
            'on the strand of frame ' + r.label + '">Copy nt</button>' +
          (r.fromMet && r.fromMetNt
            ? '<button class="btn sm ghost" data-fntmet="' + i + '" title="the bases from the first Met on">nt from Met</button>' : '') +
          '<button class="btn sm ghost" data-fblast="' + i + '">BLAST</button>' +
          (r.fromMet && r.fromMet.length !== r.longest.length
            ? '<button class="btn sm ghost" data-fmet="' + i + '">from Met</button>' : '') +
        '</td></tr>').join('') +
      '</tbody></table></div>' +
      '<div class="row"><span class="badge ' + VERDICT_STYLE[tri.verdict] + '">' +
        esc(B.VERDICTS[tri.verdict]) + '</span></div>' +
      tri.flags.map(f => '<div class="note ' + (f.level === 'ok' ? 'ok' : f.level === 'bad' ? 'bad' : f.level === 'warn' ? 'warn' : '') +
        '" style="font-size:11.5px">' + esc(f.text) + '</div>').join('') +
      '<div class="hint">"Clean run" is the longest stretch with no stop codon and no ambiguous residue, which is what ' +
      'the frames are ranked on — the plain no-stop length is shown next to it because an N translates to X, not to a stop. ' +
      '"Fit" is the chi-square per residue of the amino-acid composition against the Swiss-Prot average; a real frame sits below 0.3. ' +
      'BLAST copies the peptide and opens NCBI blastp — paste it if the box comes up empty. ' +
      '<em>Copy nt</em> gives that peptide’s own bases on the strand of its frame, and <em>nt from Met</em> ' +
      'the same from the first ATG on, which is usually what you want to clone.</div></div></div>';
  }

  /* ---------------- protein read off the reads ----------------
     Vote the reads onto the reference, then translate the coding regions
     from what the reads actually support. A residue whose codon no read
     covered is not a result, so coverage is reported beside every protein
     and nothing is called verified without it.                          */
  function proteinPanel(rec) {
    const key = rec.id + ':' + rec.updatedAt + ':' + state.reads.map(r => r.name + '/' + r.seq.length).join('|');
    let pile;
    if (state.pileCache && state.pileCache.key === key) pile = state.pileCache.pile;
    else {
      try {
        pile = B.readPileup(rec.seq, state.reads.map(r => ({
          name: r.name, seq: (r.trim && r.trim.seq) || r.seq, aln: r.aln || null })), {});
      } catch (e) { return ''; }
      state.pileCache = { key: key, pile: pile };
    }
    state.verifyPile = pile;
    const prots = [];
    const code = f => parseInt((f.quals && f.quals.transl_table), 10) || state.codeTable;
    const frameOf = f => { const c = f.quals && parseInt(f.quals.codon_start, 10); return c >= 1 && c <= 3 ? c - 1 : 0; };
    rec.features.forEach((f, idx) => {
      if (!/^CDS$/i.test(f.type) || f.between) return;
      const fp = PB.featurePositions(rec, f);
      if (fp.pos.length < 75) return;                  // a tag or a start codon is not a protein
      const cov = B.coverageOver(pile, fp.pos);
      if (cov.covered === 0) return;                   // the reads never reach it
      const fr = frameOf(f);
      const got = B.translate(B.featureSequence(pile.consensus, f).slice(fr), { table: code(f) }).replace(/\*$/, '');
      const want = B.translate(B.featureSequence(rec.seq, f).slice(fr), { table: code(f) }).replace(/\*$/, '');
      const diffs = [];
      for (let i = 0; i < Math.max(got.length, want.length); i++) {
        if (got[i] !== want[i]) diffs.push((want[i] || '-') + (i + 1) + (got[i] || '-'));
      }
      prots.push({ name: f.label || ('CDS ' + (idx + 1)), seq: got, expected: want, cov, diffs, source: 'annotated' });
    });
    // anything coding that the reads cover but the record does not annotate
    const covered = [];
    let runStart = -1;
    for (let i = 0; i <= pile.coverage.length; i++) {
      const c = i < pile.coverage.length ? pile.coverage[i] : 0;
      if (c > 0 && runStart === -1) runStart = i;
      else if (c === 0 && runStart !== -1) { covered.push([runStart, i]); runStart = -1; }
    }
    covered.forEach(([a, b]) => {
      if (b - a < 150) return;
      let orfs = [];
      try { orfs = B.findORFs(pile.consensus.slice(a, b), { table: state.codeTable, minAA: 50, requireStart: true }); }
      catch (e) { return; }
      orfs.forEach(o => {
        const st = a + o.start, en = a + o.end;
        if (prots.some(p => p.source === 'annotated' && p.seq && o.protein && (p.seq.indexOf(o.protein) !== -1 || o.protein.indexOf(p.seq) !== -1))) return;
        const pos = []; for (let i = st; i < en; i++) pos.push(i);
        prots.push({ name: 'unannotated ORF ' + fmt(st + 1) + '–' + fmt(en), seq: o.protein,
                     expected: null, cov: B.coverageOver(pile, pos), diffs: [], source: 'orf' });
      });
    });
    const isVerified = p => p.cov.covered >= 100 && !p.cov.conflicts;
    const nOk = prots.filter(isVerified).length;
    const verdict = p => {
      if (p.cov.conflicts) return '<span class="badge bad">reads disagree at ' + p.cov.conflicts + ' base' + (p.cov.conflicts === 1 ? '' : 's') + '</span>';
      if (p.cov.covered < 100) return '<span class="badge warn">' + f1(p.cov.covered) + '% covered</span>';
      if (p.expected === null) return '<span class="badge neutral">not in the reference</span>';
      if (!p.diffs.length) return '<span class="badge ok">matches the reference</span>';
      return '<span class="badge warn">' + p.diffs.length + ' residue change' + (p.diffs.length === 1 ? '' : 's') + '</span>';
    };
    let h = '<div class="panel"><h3>Protein from these reads ' +
      '<span class="sub">' + pile.reads.length + ' read' + (pile.reads.length === 1 ? '' : 's') + ' voted onto “' + esc(rec.name) + '”</span></h3>';
    if (pile.rejected.length) h += '<div class="note warn">Left out: ' +
      pile.rejected.map(r => esc(r.name) + ' (' + esc(r.why) + ')').join('; ') + '.</div>';
    if (!prots.length) {
      h += '<div class="hint">The reads do not cover any coding region of this construct, and hold no open reading frame of 50 aa or more.</div></div>';
      return h;
    }
    h += '<div class="scrollx"><table class="tbl"><thead><tr><th>Protein</th><th class="num">aa</th>' +
      '<th class="num">Covered</th><th class="num">Min depth</th><th>Status</th><th>Changes</th><th></th></tr></thead><tbody>' +
      prots.map((p, i) => '<tr>' +
        '<td>' + esc(p.name) + (p.source === 'orf' ? ' <span class="badge neutral">found in the reads</span>' : '') + '</td>' +
        '<td class="num">' + fmt(p.seq.length) + '</td>' +
        '<td class="num">' + f1(p.cov.covered) + '%</td>' +
        '<td class="num">' + p.cov.minDepth + '×</td>' +
        '<td>' + verdict(p) + '</td>' +
        '<td class="mono" style="font-size:11.5px">' + esc(p.diffs.slice(0, 6).join(', ') + (p.diffs.length > 6 ? ' +' + (p.diffs.length - 6) : '')) + '</td>' +
        '<td style="white-space:nowrap">' +
          '<button class="btn sm ghost" data-pcopy="' + i + '">Copy</button>' +
          '<button class="btn sm ghost" data-pblast="' + i + '">BLAST</button></td></tr>').join('') +
      '</tbody></table></div>' +
      '<div class="row" style="margin-top:8px"><button class="btn sm" id="pFasta">Copy FASTA</button>' +
      '<button class="btn sm" id="pFastaDl">Download FASTA</button>' +
      '<label class="chip' + (state.protVerifiedOnly ? ' on' : '') + '" id="pOnlyOk" style="align-self:center">verified only</label>' +
      '<span class="hint" style="margin:0;align-self:center">' + nOk + ' of ' + prots.length + ' fully covered with no conflicting reads</span></div>' +
      '<div class="hint">Each protein is translated from the bases the reads support, not from the stored construct. ' +
      '<strong>Covered</strong> is the share of its codons at least one read reached and <strong>min depth</strong> the thinnest point — anything under 100% is partly the reference’s own sequence, so it is not verified. ' +
      'Where reads disagree, no call is made; open the read’s details to see the position and its chromatogram. ' +
      'Coding features under 25 residues — tags, start codons — are left out. ' +
      'BLAST copies the sequence and opens NCBI blastp — if the box comes up empty, paste it.</div></div>';
    state.verifyProteins = prots;
    return h;
  }

  /* What one difference does, for the differences table.

     This used to work the codon out here, by slicing the feature's
     start..end out of the sequence. That is right only for a plain
     forward feature: on a spliced gene the slice includes the introns, so
     every codon index past the first intron is wrong, and on a gene that
     wraps the origin the slice returns nothing at all and no consequence
     was reported. The engine now keeps an explicit genomic-position to
     codon map that handles both, and this is a thin wrapper over it. */
  function annotateDiff(rec, d) {
    const out = { feature: null, effect: '—' };
    if (!rec) return out;
    for (const f of rec.features) {
      const pieces = B.featurePieces(f, rec.seq.length);
      if (!pieces.some(p => d.refPos >= p.start && d.refPos < p.end)) continue;
      out.feature = f.label;
      if (!B.isCodingFeature(f)) { out.effect = 'non-coding'; continue; }
      const e = B.diffEffect(rec, f, d, { table: state.codeTable });
      if (!e) { out.effect = 'non-coding'; continue; }
      const badge = { silent: 'ok', missense: 'warn', nonsense: 'bad',
                      readthrough: 'bad', frameshift: 'bad', partial: 'neutral' }[e.kind] || 'neutral';
      out.effect = '<span class="badge ' + badge + '">' + esc(e.kind) + '</span> ' +
        esc(e.text.replace(/^[a-z ]+\u00b7\s*/, ''));
      out.detail = e;
      return out;
    }
    return out;
  }

  /* ---------------- BLAST results, pasted back ----------------
     A published page cannot reach NCBI: the sandbox blocks outside
     requests, and the URL API wants a polling loop against a Request ID
     besides. So the round trip is one paste. The parser takes the four
     shapes NCBI's Download menu produces plus the pairwise text people
     copy off the screen, and what comes back is a real table in the tool
     rather than another browser tab to keep open. */
  function blastPanel() {
    const b = state.blast;
    let h = '<div class="panel"><h3>BLAST results ' +
      '<span class="sub">' + (b ? esc(b.query || 'pasted') + ' · ' + esc(b.format) + ' · ' + b.hits.length + ' hit' + (b.hits.length === 1 ? '' : 's') : 'paste what NCBI gave you') + '</span></h3>' +
      '<div class="row"><button class="btn" id="blPaste">Paste BLAST results</button>' +
      (b ? '<button class="btn sm ghost" id="blClear">Clear</button>' : '') + '</div>';
    if (!b) {
      h += '<div class="hint">Run BLAST from any of the BLAST buttons above, then on the results page use ' +
        '<strong>Download → Hit Table (text)</strong>, <strong>Descriptions (CSV)</strong> or ' +
        '<strong>Single-file JSON</strong> and paste the file in here. The pairwise text copied straight off the ' +
        'screen works too. Nothing is fetched: this page has no network access, which is also why no sequence of ' +
        'yours ever leaves it.</div></div>';
      return h;
    }
    const n = v => (v == null ? '—' : v);
    const ev = v => v == null ? '—' : (v === 0 ? '0.0' : (v < 1e-4 ? v.toExponential(0) : String(v)));
    h += '<div class="scrollx"><table class="tbl"><thead><tr><th>Accession</th><th>Description</th>' +
      '<th class="num">Identity</th><th class="num">Cover</th><th class="num">E value</th><th class="num">Bits</th>' +
      '<th class="num">HSPs</th><th></th></tr></thead><tbody>' +
      b.hits.slice(0, 200).map(hit => '<tr>' +
        '<td class="mono">' + esc(hit.accession || hit.id || '—') + '</td>' +
        '<td>' + esc(PB.trimText(hit.title || '', 90)) + '</td>' +
        '<td class="num">' + (hit.identity == null ? '—' : f1(hit.identity) + '%') + '</td>' +
        '<td class="num">' + (hit.cover == null ? '—' : f1(hit.cover) + '%') + '</td>' +
        '<td class="num mono">' + ev(hit.evalue) + '</td>' +
        '<td class="num">' + n(hit.bits) + '</td>' +
        '<td class="num">' + (hit.hsps && hit.hsps.length ? hit.hsps.length : '—') + '</td>' +
        '<td>' + (hit.accession ? '<a class="btn sm ghost" target="_blank" rel="noopener" href="https://www.ncbi.nlm.nih.gov/search/all/?term=' +
          encodeURIComponent(hit.accession) + '">NCBI</a>' : '') + '</td></tr>').join('') +
      '</tbody></table></div>' +
      '<div class="hint">Read as <strong>' + esc(b.format) + '</strong>. A column the format does not carry is left ' +
      'blank rather than guessed — a Hit Table has no query cover in it, and a Descriptions CSV has no HSP ' +
      'coordinates. The NCBI link opens the record; it is not fetched here.</div></div>';
    return h;
  }
  PB.blastPanel = blastPanel;

  function blastPasteDialog() {
    PB.modal({
      title: 'Paste BLAST results',
      body: '<div class="row"><div class="field wide"><label>Paste the file, or the text off the results page</label>' +
        '<textarea class="inp mono" id="blText" rows="12" spellcheck="false" placeholder="Hit table, Descriptions CSV, single-file JSON, or pairwise text"></textarea></div></div>' +
        '<div class="row"><div class="field"><label>What was this a BLAST of? (optional)</label>' +
        '<input class="inp" id="blQuery" placeholder="e.g. EGFP_F.ab1" value="' +
        esc(state.blastFrom || '') + '"></div></div>' +
        '<div class="hint">From the BLAST results page: <strong>Download</strong> → Hit Table (text), ' +
        'Descriptions (CSV) or Single-file JSON. All three paste in here, and so does the pairwise text.</div>',
      actions: [{ label: 'Cancel' }, { label: 'Read it', primary: true, onClick: () => {
        const txt = document.querySelector('#blText').value;
        const q = document.querySelector('#blQuery').value.trim();
        let r;
        try { r = B.parseBlast(txt); }
        catch (e) { PB.toast(e.message, 'bad'); return false; }
        if (!r.hits.length) { PB.toast('Nothing to read in that.', 'warn'); return false; }
        state.blast = { query: q, format: r.format, hits: r.hits };
        PB.closeModal(); PB.renderVerify();
        PB.toast('Read ' + r.hits.length + ' hit' + (r.hits.length === 1 ? '' : 's') + ' from the ' + r.format.toLowerCase(), 'ok');
        return false;
      } }]
    });
  }

  /* ---------------- what this read actually landed on ----------------
     The question at the bench is never "how many bases differ" but "which
     part of my protein did this read see, and is it right there". This
     draws the read's span against the construct, lists every feature it
     touched, and for the coding ones gives the residue range it reached
     and the protein that stretch of the read actually encodes. */
  function coveragePanel(rec, rd) {
    if (!rec || !rd.aln) return '';
    const aln = rd.aln, L = rec.seq.length;
    const rows = B.readOverlap(rec, aln, { table: state.codeTable });
    const W = 900, H = 58, pad = 6, trackY = 26, trackH = 9;
    const x = p => pad + (p / L) * (W - pad * 2);
    let g = '<svg class="covmap" viewBox="0 0 ' + W + ' ' + H + '" width="100%" height="' + H +
      '" role="img" aria-label="where this read lands on the construct">';
    /* the construct */
    g += '<rect x="' + pad + '" y="' + trackY + '" width="' + (W - pad * 2) + '" height="' + trackH +
      '" rx="2" fill="var(--surface-3)"/>';
    /* every feature, in its own colour */
    rec.features.forEach((f, i) => {
      if (f.virtual) return;
      B.featurePieces(f, L).forEach(pc => {
        g += '<rect x="' + x(pc.start).toFixed(1) + '" y="' + trackY + '" width="' +
          Math.max(1, x(pc.end) - x(pc.start)).toFixed(1) + '" height="' + trackH +
          '" fill="' + PB.featColor(f) + '" opacity="0.55"><title>' + esc(f.label) + '</title></rect>';
      });
    });
    /* the read on top */
    const rs = x(aln.refStart), re = x(Math.min(L, aln.refEnd));
    g += '<rect class="covread" x="' + rs.toFixed(1) + '" y="' + (trackY - 7) + '" width="' + Math.max(2, re - rs).toFixed(1) +
      '" height="' + (trackH + 14) + '" rx="2" fill="none" stroke="var(--accent)" stroke-width="1.6"/>';
    g += '<rect x="' + rs.toFixed(1) + '" y="' + (trackY - 7) + '" width="' + Math.max(2, re - rs).toFixed(1) +
      '" height="5" fill="var(--accent)"/>';
    /* each difference as a tick */
    (aln.differences || []).forEach(d => {
      g += '<line class="covdiff" x1="' + x(d.refPos).toFixed(1) + '" y1="' + (trackY - 7) + '" x2="' + x(d.refPos).toFixed(1) +
        '" y2="' + (trackY + trackH + 7) + '" stroke="var(--bad)" stroke-width="1"/>';
    });
    g += '<text x="' + pad + '" y="14" font-size="10" fill="var(--fg-faint)">1</text>' +
      '<text x="' + (W - pad) + '" y="14" font-size="10" text-anchor="end" fill="var(--fg-faint)">' + fmt(L) + '</text>' +
      '<text x="' + ((rs + re) / 2).toFixed(1) + '" y="' + (H - 4) + '" font-size="10" text-anchor="middle" fill="var(--accent)">' +
      fmt(aln.refStart + 1) + '–' + fmt(aln.refEnd) + '</text></svg>';

    let h = '<div class="row"><div class="field wide"><label>Where this read lands</label>' + g + '</div></div>';
    if (!rows.length) {
      return h + '<div class="hint">The read lands on no annotated feature of this construct.</div>';
    }
    const pct = v => f1(v * 100) + '%';
    h += '<div class="scrollx"><table class="tbl sm"><thead><tr><th>Feature</th><th>Type</th>' +
      '<th class="num">Bases read</th><th class="num">Of feature</th><th>Residues reached</th><th>Protein over that stretch</th>' +
      '</tr></thead><tbody>';
    rows.forEach(r => {
      let prot = '<span class="dim">—</span>';
      if (r.coding && r.residues) {
        const p = B.readProtein(rec, aln, rec.features[r.index], { table: state.codeTable });
        if (p) {
          if (!p.covered) prot = '<span class="badge neutral">no whole codon covered</span>';
          else if (p.identical) prot = '<span class="badge ok">' + p.covered + ' of ' + p.length +
            ' residues read, no change in them</span>';
          else prot = '<span class="badge warn">' + p.changes.length + ' residue change' +
            (p.changes.length === 1 ? '' : 's') + '</span> <span class="mono">' +
            esc(p.changes.slice(0, 5).map(c => c.from + c.residue + c.to).join(', ')) +
            (p.changes.length > 5 ? ' +' + (p.changes.length - 5) : '') + '</span>';
        }
      }
      h += '<tr><td>' + esc(r.label) + '</td><td class="mono" style="font-size:11.5px">' + esc(r.type) + '</td>' +
        '<td class="num">' + fmt(r.bases) + '</td>' +
        '<td class="num">' + pct(r.fraction) + '</td>' +
        '<td class="mono">' + (r.residues ? r.residues.from + '–' + r.residues.to + ' of ' + (r.aaTotal || '?') : '—') + '</td>' +
        '<td>' + prot + '</td></tr>';
    });
    h += '</tbody></table></div>' +
      '<div class="hint">Residue numbers are counted along the coding sequence itself, so they are right for a gene ' +
      'written as a join across an intron or across the origin of the plasmid — where slicing start to end would ' +
      'give the wrong codon, or nothing at all. <strong>Protein over that stretch</strong> is translated from the ' +
      'read’s own bases; residues the read never reached are not filled in from the reference and are not counted.</div>';
    return h;
  }

  function alignBlocks(a, b, nameA, nameB, width, refStart) {
    let out = '';
    const pad = s => (s + '      ').slice(0, 6);
    for (let i = 0; i < a.length; i += width) {
      const as = a.substr(i, width), bs = b.substr(i, width);
      let mid = '', aOut = '', bOut = '';
      for (let k = 0; k < as.length; k++) {
        const x = as[k], y = bs[k] || ' ';
        if (x === '-' || y === '-') { mid += ' '; aOut += '<span class="gp">' + esc(x) + '</span>'; bOut += '<span class="gp">' + esc(y) + '</span>'; }
        else if (x === y) { mid += '|'; aOut += esc(x); bOut += esc(y); }
        else { mid += '·'; aOut += '<span class="mm">' + esc(x) + '</span>'; bOut += '<span class="mm">' + esc(y) + '</span>'; }
      }
      const refCoord = (refStart || 0) + a.substr(0, i).replace(/-/g, '').length + 1;
      out += '<span class="lbl">' + pad(nameA) + String(refCoord).padStart(7) + ' </span>' + aOut + '\n' +
             '<span class="lbl">' + pad('') + '        </span><span class="mid">' + mid + '</span>\n' +
             '<span class="lbl">' + pad(nameB) + '        </span>' + bOut + '\n\n';
    }
    return '<div class="aln-box">' + out + '</div>';
  }
  PB.alignBlocks = alignBlocks;

  function qualitySVG(rd) {
    const q = rd.quality;
    const W = 860, H = 72, pad = 4;
    const n = q.length;
    const step = Math.max(1, Math.floor(n / 600));
    let pts = '';
    for (let i = 0; i < n; i += step) {
      const x = pad + (i / n) * (W - 2 * pad);
      const y = H - pad - Math.min(q[i], 60) / 60 * (H - 2 * pad - 8);
      pts += x.toFixed(1) + ',' + y.toFixed(1) + ' ';
    }
    let s = '<svg class="chroma" viewBox="0 0 ' + W + ' ' + H + '" preserveAspectRatio="none" role="img" aria-label="Read quality">';
    s += '<rect width="' + W + '" height="' + H + '" fill="var(--bg)"/>';
    [20, 40].forEach(t => {
      const y = H - pad - t / 60 * (H - 2 * pad - 8);
      s += '<line x1="' + pad + '" y1="' + y.toFixed(1) + '" x2="' + (W - pad) + '" y2="' + y.toFixed(1) +
        '" stroke="var(--border)" stroke-dasharray="3 3"/>';
      s += '<text x="' + (pad + 2) + '" y="' + (y - 2).toFixed(1) + '" fill="var(--fg-faint)">Q' + t + '</text>';
    });
    if (rd.trim) {
      const x0 = pad + (rd.trim.start / n) * (W - 2 * pad);
      const x1 = pad + (rd.trim.end / n) * (W - 2 * pad);
      s += '<rect x="' + x0.toFixed(1) + '" y="' + pad + '" width="' + Math.max(1, x1 - x0).toFixed(1) +
        '" height="' + (H - 2 * pad) + '" fill="var(--accent)" opacity=".10"/>';
    }
    s += '<polyline points="' + pts + '" fill="none" stroke="var(--accent)" stroke-width="1.2"/>';
    s += '</svg>';
    return s;
  }

  function chromaSVG(rd, baseIdx) {
    const locs = rd.peakLocations;
    const chans = ['A', 'C', 'G', 'T'];
    const colorOf = { A: 'var(--base-a)', C: 'var(--base-c)', G: 'var(--base-g)', T: 'var(--base-t)' };
    if (!locs || !locs.length) return '<div class="hint">This file has base calls but no peak positions, so a chromatogram cannot be aligned to the bases.</div>';
    const half = 12;
    const i0 = Math.max(0, baseIdx - half), i1 = Math.min(rd.seq.length - 1, baseIdx + half);
    const t0 = locs[i0], t1 = locs[i1];
    if (t0 == null || t1 == null || t1 <= t0) return '<div class="hint">Peak positions do not cover this base.</div>';
    const W = 860, H = 150, padL = 6, padB = 20;
    let maxV = 1;
    for (const c of chans) {
      const tr = rd.traces[c]; if (!tr) continue;
      for (let t = t0; t <= t1; t++) if (tr[t] > maxV) maxV = tr[t];
    }
    const xOf = t => padL + ((t - t0) / (t1 - t0)) * (W - 2 * padL);
    const yOf = v => (H - padB) - (v / maxV) * (H - padB - 8);
    let s = '<svg class="chroma" viewBox="0 0 ' + W + ' ' + H + '" role="img" aria-label="Chromatogram">';
    s += '<rect width="' + W + '" height="' + H + '" fill="var(--bg)"/>';
    // highlight the base in question
    if (locs[baseIdx] != null) {
      const bx = xOf(locs[baseIdx]);
      s += '<rect x="' + (bx - 7).toFixed(1) + '" y="0" width="14" height="' + (H - padB + 4) + '" fill="var(--sel-strong)" opacity=".22"/>';
    }
    for (const c of chans) {
      const tr = rd.traces[c];
      if (!tr) continue;
      let pts = '';
      for (let t = t0; t <= t1; t++) pts += xOf(t).toFixed(1) + ',' + yOf(tr[t] || 0).toFixed(1) + ' ';
      s += '<polyline points="' + pts + '" fill="none" stroke="' + colorOf[c] + '" stroke-width="1.3" opacity=".9"/>';
    }
    for (let i = i0; i <= i1; i++) {
      if (locs[i] == null) continue;
      const base = rd.seq[i] || 'N';
      s += '<text x="' + xOf(locs[i]).toFixed(1) + '" y="' + (H - 6) + '" text-anchor="middle" fill="' +
        (colorOf[base] || 'var(--fg-faint)') + '" style="font-weight:' + (i === baseIdx ? '700' : '400') + '">' + esc(base) + '</text>';
    }
    s += '</svg>';
    return s;
  }

  function blastProtein(seq, what) {
    PB.state.blastFrom = what;
    PB.copyText(seq, what + ' (also on the clipboard)');
    window.open('https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=Proteins&PROGRAM=blastp' +
      '&BLAST_PROGRAMS=blastp&DATABASE=nr&CMD=Web&QUERY=' + encodeURIComponent(seq), '_blank', 'noopener');
  }

  function wireVerify(host, rec) {
    const prots = state.verifyProteins || [];
    /* The header carries where the sequence came from and how well the reads
       covered it, so the file still says what it is once it leaves here. */
    const statusOf = p => p.cov.conflicts ? 'reads_disagree'
      : p.cov.covered < 100 ? 'partial_coverage'
      : p.expected === null ? 'unannotated_ORF'
      : p.diffs.length ? 'changed' : 'verified';
    const fastaOf = p => '>' + p.name.replace(/\s+/g, '_') + ' ' + p.seq.length + 'aa' +
      ' covered=' + p.cov.covered.toFixed(1) + '%' +
      ' min_depth=' + p.cov.minDepth + 'x ' + statusOf(p) +
      (p.diffs.length ? ' changes=' + p.diffs.join(',') : '') + '\n' +
      (p.seq.match(/.{1,60}/g) || []).join('\n') + '\n';
    const chosen = () => prots.filter(p => !state.protVerifiedOnly || (p.cov.covered >= 100 && !p.cov.conflicts));
    host.querySelectorAll('[data-pcopy]').forEach(b => b.addEventListener('click', () => {
      const p = prots[+b.getAttribute('data-pcopy')];
      if (p) PB.copyText(fastaOf(p), p.name);
    }));
    host.querySelectorAll('[data-pblast]').forEach(b => b.addEventListener('click', () => {
      const p = prots[+b.getAttribute('data-pblast')];
      if (!p) return;
      state.blastFrom = p.name;
      PB.copyText(p.seq, p.name + ' (also on the clipboard)');
      const u = 'https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=Proteins&PROGRAM=blastp' +
                '&BLAST_PROGRAMS=blastp&DATABASE=nr&CMD=Web&QUERY=' + encodeURIComponent(p.seq);
      window.open(u, '_blank', 'noopener');
    }));
    const pf = $('#pFasta');
    if (pf) pf.addEventListener('click', () => {
      const list = chosen();
      if (!list.length) { PB.toast('Nothing to copy with that filter.', 'warn'); return; }
      PB.copyText(list.map(fastaOf).join(''), list.length + ' protein' + (list.length === 1 ? '' : 's'));
    });
    const pfd = $('#pFastaDl');
    if (pfd) pfd.addEventListener('click', () => {
      const list = chosen();
      if (!list.length) { PB.toast('Nothing to export with that filter.', 'warn'); return; }
      const rec = PB.active();
      const name = (rec ? rec.name : 'reads').replace(/[^\w.-]+/g, '_').slice(0, 40);
      PB.offerFile(name + '_proteins.txt', list.map(fastaOf).join(''),
        list.length + ' protein' + (list.length === 1 ? '' : 's') + ' as FASTA');
    });
    const po = $('#pOnlyOk');
    if (po) po.addEventListener('click', () => { state.protVerifiedOnly = !state.protVerifiedOnly; renderVerify(); });
    /* translate panel */
    const tr = state.trItems || [];
    const pepOf = x => x.f.clean || x.f.longest;
    /* whichever peptide pepOf picked, these are its own bases */
    const ntOf = x => x.f.clean ? x.f.cleanNt : x.f.longestNt;
    /* The header carries the verdict, so a peptide that leaves here cannot
       be mistaken later for a confident translation when it was not one. */
    const trFasta = x => {
      const pep = pepOf(x);
      const pins = (x.f.landmarks || []).filter(l => l.pins).map(l => l.name.replace(/\s+/g, '_'));
      return '>' + x.read.replace(/\s+/g, '_') + '_frame' + x.f.label.replace('−', 'minus') +
        ' ' + pep.length + 'aa stops=' + x.f.stops +
        ' ambiguous=' + (x.t.ambiguousFrac * 100).toFixed(1) + '%' +
        (x.f.fit === null ? '' : ' composition=' + x.f.fit.toFixed(2)) +
        ' verdict=' + x.t.verdict + (pins.length ? ' pinned_by=' + pins.join(',') : '') +
        (x.alt ? ' runner_up_frame' : '') + '\n' +
        (pep.match(/.{1,60}/g) || []).join('\n') + '\n';
    };
    const trChosen = () => tr.filter(x => !state.trOnlyOk || WORTH[x.t.verdict]);
    host.querySelectorAll('[data-trm]').forEach(b => b.addEventListener('click', () => {
      state.trMode = b.getAttribute('data-trm'); renderVerify();
    }));
    const too = $('#trOnlyOk');
    if (too) too.addEventListener('click', () => { state.trOnlyOk = !state.trOnlyOk; renderVerify(); });
    host.querySelectorAll('[data-trc]').forEach(b => b.addEventListener('click', () => {
      const x = tr[+b.getAttribute('data-trc')];
      if (x) PB.copyText(pepOf(x), x.read + ' frame ' + x.f.label);
    }));
    /* The bases a peptide was translated from. On a minus frame these are
       the reverse complement of the read, which is the useful thing to hand
       over — they translate to the peptide on screen and can go straight
       into a primer or an order — but it is not what is in the trace file,
       so the message says so rather than leaving it to be discovered. */
    const ntNote = (f, nt) => 'frame ' + f.label + ', ' + nt.length + ' bp' +
      (f.dir === -1 ? ' (reverse complement of the read)' : '');
    host.querySelectorAll('[data-trnt]').forEach(b => b.addEventListener('click', () => {
      const x = tr[+b.getAttribute('data-trnt')];
      if (!x) return;
      const nt = ntOf(x);
      if (!nt) { PB.toast('That frame has no bases to copy.', 'warn'); return; }
      PB.copyText(nt, x.read + ' \u00b7 ' + ntNote(x.f, nt));
    }));
    host.querySelectorAll('[data-trb]').forEach(b => b.addEventListener('click', () => {
      const x = tr[+b.getAttribute('data-trb')];
      if (!x) return;
      if (!WORTH[x.t.verdict])
        PB.toast('This read is ' + B.VERDICTS[x.t.verdict].toLowerCase() + ' — BLAST it if you want, ' +
          'but the peptide may not mean anything.', 'warn');
      blastProtein(pepOf(x), x.read + ' frame ' + x.f.label);
    }));
    host.querySelectorAll('[data-trm2]').forEach(b => b.addEventListener('click', () => {
      const x = tr[+b.getAttribute('data-trm2')];
      if (x) blastProtein(x.f.fromMet, x.read + ' frame ' + x.f.label + ' from Met');
    }));
    const trf = $('#trFasta');
    if (trf) trf.addEventListener('click', () => {
      const list = trChosen();
      if (!list.length) { PB.toast('Nothing to copy with that filter.', 'warn'); return; }
      PB.copyText(list.map(trFasta).join(''), list.length + ' translations');
    });
    /* The whole triage as a spreadsheet-ready table: one row per read, every
       signal the verdict was made from, so a plate can be reviewed or filed. */
    const trt = $('#trCsv');
    if (trt) trt.addEventListener('click', () => {
      const csv = triageCsv();
      if (!csv) { PB.toast('No reads to export.', 'warn'); return; }
      const rec0 = PB.active();
      const nm = (rec0 ? rec0.name : 'reads').replace(/[^\w.-]+/g, '_').slice(0, 40);
      PB.offerFile(nm + '_read_triage.csv', csv,
        state.triageCache.rows.length + ' reads as a table');
    });
    const trd = $('#trDl');
    if (trd) trd.addEventListener('click', () => {
      const list = trChosen();
      if (!list.length) { PB.toast('Nothing to export with that filter.', 'warn'); return; }
      const rec0 = PB.active();
      const nm = (rec0 ? rec0.name : 'reads').replace(/[^\w.-]+/g, '_').slice(0, 40);
      PB.offerFile(nm + '_read_translations.txt', list.map(trFasta).join(''),
        list.length + ' translations as FASTA');
    });

    const frames = state.verifyFrames || [];
    const blastAA = blastProtein;
    host.querySelectorAll('[data-fcopy]').forEach(b => b.addEventListener('click', () => {
      const r = frames[+b.getAttribute('data-fcopy')];
      if (r) PB.copyText(r.clean || r.longest, 'Frame ' + r.label);
    }));
    host.querySelectorAll('[data-fnt]').forEach(b => b.addEventListener('click', () => {
      const r = frames[+b.getAttribute('data-fnt')];
      if (!r) return;
      const nt = r.clean ? r.cleanNt : r.longestNt;
      if (!nt) { PB.toast('That frame has no bases to copy.', 'warn'); return; }
      PB.copyText(nt, 'Frame ' + r.label + ', ' + nt.length + ' bp' +
        (r.dir === -1 ? ' (reverse complement of the read)' : ''));
    }));
    host.querySelectorAll('[data-fntmet]').forEach(b => b.addEventListener('click', () => {
      const r = frames[+b.getAttribute('data-fntmet')];
      if (!r || !r.fromMetNt) return;
      PB.copyText(r.fromMetNt, 'Frame ' + r.label + ' from the first Met, ' + r.fromMetNt.length + ' bp' +
        (r.dir === -1 ? ' (reverse complement of the read)' : ''));
    }));
    host.querySelectorAll('[data-fblast]').forEach(b => b.addEventListener('click', () => {
      const r = frames[+b.getAttribute('data-fblast')];
      if (r) blastAA(r.clean || r.longest, 'Frame ' + r.label);
    }));
    host.querySelectorAll('[data-fmet]').forEach(b => b.addEventListener('click', () => {
      const r = frames[+b.getAttribute('data-fmet')];
      if (r) blastAA(r.fromMet, 'Frame ' + r.label + ' from Met');
    }));
    const blp = host.querySelector('#blPaste');
    if (blp) blp.addEventListener('click', blastPasteDialog);
    const blc = host.querySelector('#blClear');
    if (blc) blc.addEventListener('click', () => { state.blast = null; renderVerify(); });
    host.querySelectorAll('[data-vnblast]').forEach(b => b.addEventListener('click', () => {
      const rd = state.reads[+b.getAttribute('data-vnblast')];
      if (!rd) return;
      const seq = (rd.trim && rd.trim.seq.length > 40 ? rd.trim.seq : rd.seq);
      state.blastFrom = rd.name;
      PB.copyText(seq, rd.name + ' (also on the clipboard)');
      window.open('https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=MegaBlast&PROGRAM=blastn' +
        '&BLAST_PROGRAMS=megaBlast&DATABASE=nt&CMD=Web&QUERY=' + encodeURIComponent(seq), '_blank', 'noopener');
    }));
    const ab1 = $('#vAb1');
    if (ab1) ab1.addEventListener('click', () => PB.pickReadFiles());
    const paste = $('#vPaste');
    if (paste) paste.addEventListener('click', () => {
      PB.modal({
        title: 'Paste sequencing reads',
        body: '<div class="row"><div class="field"><label>Name (used when there is no FASTA header)</label>' +
          '<input class="inp" id="rdName" placeholder="clone 7 fwd"></div></div>' +
          '<div class="row"><div class="field wide"><label>Base calls</label>' +
          '<textarea class="inp" id="rdSeq" rows="9" spellcheck="false" placeholder="One read, or many in FASTA:&#10;>JM02-T7&#10;ACGT…&#10;>JM03-T7&#10;ACGT…"></textarea></div></div>' +
          '<div class="hint" id="rdCount">Each FASTA record becomes its own read, named after its header.</div>',
        actions: [{ label: 'Cancel' }, { label: 'Add', primary: true, onClick: () => {
          const raw = $('#rdSeq').value;
          const got = B.parseReadsText(raw, ($('#rdName').value || 'pasted read').trim());
          const usable = got.filter(r => r.seq.length >= 20);
          if (!usable.length) { PB.toast('No read of 20 bases or more in that text.', 'bad'); return false; }
          if (usable.length < got.length)
            PB.toast((got.length - usable.length) + ' record' + (got.length - usable.length === 1 ? ' was' : 's were') +
              ' under 20 bases and left out.', 'warn');
          PB.addReads(usable.map(r => Object.assign(r, { source: 'pasted' })));
        } }],
        onOpen: w => {
          const ta = w.querySelector('#rdSeq'), cnt = w.querySelector('#rdCount');
          ta.addEventListener('input', () => {
            const n = B.parseReadsText(ta.value, 'x').length;
            cnt.textContent = n > 1 ? n + ' reads will be added, one per FASTA record.'
              : 'Each FASTA record becomes its own read, named after its header.';
          });
        }
      });
    });
    const clr = $('#vClear');
    if (clr) clr.addEventListener('click', () => {
      state.reads = []; state.verifyDetail = null; state.traceAt = null; renderVerify();
    });
    host.querySelectorAll('[data-vdet]').forEach(b => b.addEventListener('click', () => {
      state.verifyDetail = +b.getAttribute('data-vdet'); state.traceAt = null; renderVerify();
    }));
    host.querySelectorAll('[data-vdel]').forEach(b => b.addEventListener('click', () => {
      state.reads.splice(+b.getAttribute('data-vdel'), 1);
      state.verifyDetail = null; renderVerify();
    }));
    host.querySelectorAll('[data-vtr]').forEach(b => b.addEventListener('click', () => {
      const rd = state.reads[state.verifyDetail];
      const k = +b.getAttribute('data-vtr');
      const d = rd.aln.differences[k];
      // map a reference position back to the read index
      let readIdx = 0;
      const a = rd.aln.alignment.a, bb = rd.aln.alignment.b;
      let refPos = rd.aln.refStart, ri = 0;
      for (let i = 0; i < a.length; i++) {
        if (refPos === d.refPos) { readIdx = ri; break; }
        if (a[i] !== '-') refPos++;
        if (bb[i] !== '-') ri++;
      }
      state.traceAt = Math.max(0, Math.min(rd.seq.length - 1, readIdx + (rd.trimOffset || 0)));
      renderVerify();
    }));
    host.querySelectorAll('[data-vcopy]').forEach(b => b.addEventListener('click', () =>
      PB.copyText(state.reads[+b.getAttribute('data-vcopy')].seq, 'Read')));
    host.querySelectorAll('[data-vsave]').forEach(b => b.addEventListener('click', () => {
      const rd = state.reads[+b.getAttribute('data-vsave')];
      const nr = PB.makeRecord({ name: rd.name, seq: rd.seq, topology: 'linear' });
      PB.addHistory(nr, 'Imported read', rd.source || 'sequencing read');
      PB.addRecord(nr);
    }));
  }

  function prepareRead(rd, rec) {
    const read = Object.assign({ hasTraces: false, quality: null, traces: {}, peakLocations: null }, rd);
    if (read.quality && read.quality.length) {
      read.trim = B.trimByQuality(read.seq, read.quality, { window: 10, minMean: 20 });
    }
    const forAlign = read.trim && read.trim.seq.length > 40 ? read.trim.seq : read.seq;
    read.trimOffset = read.trim ? read.trim.start : 0;
    if (rec && rec.seq.length && forAlign.length >= 20) {
      try { read.aln = B.alignToReference(rec.seq, forAlign); } catch (e) { read.aln = null; }
      read.alnFor = rec.id + ':' + rec.updatedAt;
    }
    return read;
  }

  /* Many reads at once — a whole plate. One render and one summary at the
     end, instead of one per read. */
  /* A handful of reads is aligned straight away. A plate is aligned a read
     at a time, yielding between them, so the page stays responsive and says
     how far it has got — each alignment against a full vector takes a
     noticeable fraction of a second. Returns the count, or a promise of it
     for a large batch. */
  PB.addReads = function (list) {
    list = (list || []).filter(r => r && r.seq && r.seq.length >= 20);
    if (!list.length) return 0;
    const rec = PB.active();
    state.reads = state.reads || [];
    /* names must stay distinct, or two reads become indistinguishable in
       every table and in the exported FASTA */
    const taken = new Set(state.reads.map(r => r.name));
    const named = list.map(rd => {
      let nm = (rd.name || 'read').trim() || 'read', k = 2;
      while (taken.has(nm)) nm = (rd.name || 'read') + ' (' + (k++) + ')';
      taken.add(nm);
      return Object.assign({}, rd, { name: nm });
    });
    if (named.length <= 4) return finishAdd(named.map(rd => prepareRead(rd, rec)));
    state.readsLoading = { done: 0, total: named.length };
    const added = [];
    return new Promise(resolve => {
      let i = 0;
      const step = () => {
        const t0 = Date.now();
        /* as many as fit in ~60 ms, then let the page breathe */
        while (i < named.length && Date.now() - t0 < 60) added.push(prepareRead(named[i++], rec));
        state.readsLoading.done = i;
        PB.toast('Aligning reads… ' + i + ' of ' + named.length);
        if (i < named.length) setTimeout(step, 0);
        else { state.readsLoading = null; resolve(finishAdd(added)); }
      };
      setTimeout(step, 0);
    });
  };

  function finishAdd(added) {
    added.forEach(r => state.reads.push(r));
    state.verifyDetail = added.length === 1 ? state.reads.length - 1 : null;
    state.view = 'verify'; PB.lsSet('view', 'verify');
    PB.renderTabs(); PB.renderViews();
    if (added.length === 1) {
      const read = added[0];
      if (read.aln) PB.toast(read.name + ': ' + f1(read.aln.identity) + '% identity, ' +
        (read.aln.clean ? 'no differences' : read.aln.differences.length + ' difference(s)'),
        read.aln.clean ? 'ok' : 'warn');
      else PB.toast('Added ' + read.name, 'ok');
    } else {
      PB.toast('Added ' + added.length + ' reads', 'ok');
    }
    return added.length;
  };

  PB.addRead = function (rd) { return PB.addReads([rd]); };

  /* Reads from files: .ab1/.abi traces, and .seq / FASTA / plain text,
     which a sequencing facility usually sends alongside. Everything picked
     is read first and added together. */
  PB.loadReadFiles = function (files) {
    const list = Array.prototype.slice.call(files || []);
    if (!list.length) return Promise.resolve(0);
    const reads = [], problems = [];
    const one = file => new Promise(resolve => {
      const nm = file.name, low = nm.toLowerCase();
      const fr = new FileReader();
      fr.onerror = () => { problems.push(nm + ': could not be read'); resolve(); };
      if (/\.(ab1|abi)$/.test(low)) {
        fr.onload = () => {
          try {
            const p = B.parseAB1(fr.result);
            if (!p.sequence.length) throw new Error('no base calls');
            reads.push({ name: nm.replace(/\.(ab1|abi)$/i, ''), seq: p.sequence, quality: p.quality,
                         traces: p.traces, peakLocations: p.peakLocations, hasTraces: p.hasTraces,
                         meanQuality: p.meanQuality, source: nm, order: list.indexOf(file) });
          } catch (e) { problems.push(nm + ': ' + (e.message || 'not a readable trace')); }
          resolve();
        };
        fr.readAsArrayBuffer(file);
      } else if (/\.(scf)$/.test(low)) {
        problems.push(nm + ': SCF traces are not supported — export .ab1 or FASTA instead'); resolve();
      } else {
        fr.onload = () => {
          const txt = String(fr.result || '');
          if (PB.looksBinary && PB.looksBinary(txt)) { problems.push(nm + ': not a text sequence file'); resolve(); return; }
          const base = nm.replace(/\.(seq|fasta|fa|fas|fna|txt)$/i, '');
          const got = B.parseReadsText(txt, base);
          if (!got.length) problems.push(nm + ': no sequence found');
          got.forEach((r, i) => reads.push(Object.assign(r, { source: nm, order: list.indexOf(file) + i / 1000 })));
          resolve();
        };
        fr.readAsText(file);
      }
    });
    return Promise.all(list.map(one)).then(() => {
      reads.sort((a, b) => a.order - b.order);
      /* problems are reported after the reads are in, so the progress
         messages cannot hide them */
      return Promise.resolve(PB.addReads(reads)).then(n => {
        if (problems.length) PB.toast((n ? n + ' read' + (n === 1 ? '' : 's') + ' added. ' : '') +
          'Not loaded: ' + problems.slice(0, 3).join(' · ') +
          (problems.length > 3 ? ' · and ' + (problems.length - 3) + ' more' : ''), 'bad');
        return n;
      });
    });
  };

  PB.pickReadFiles = function () {
    const inp = document.createElement('input');
    inp.type = 'file'; inp.multiple = true;
    inp.accept = '.ab1,.abi,.seq,.fasta,.fa,.fas,.fna,.txt';
    inp.style.display = 'none';
    inp.addEventListener('change', () => { PB.loadReadFiles(inp.files); inp.remove(); });
    document.body.appendChild(inp);
    inp.click();
  };

  /* ============================== ALIGN ============================== */
  function renderAlign() {
    const host = $('#view-align');
    const mode = state.alignMode || 'pairwise';
    const a = state.alignState || {};
    let html = '<div class="panel"><h3>Alignment</h3>' +
      '<div class="chipbar">' +
      [['pairwise', 'Pairwise'], ['msa', 'Multiple'], ['contigs', 'Assemble reads']].map(([v, l]) =>
        '<button class="chip' + (mode === v ? ' on' : '') + '" data-am="' + v + '">' + l + '</button>').join('') +
      '</div>';

    if (mode === 'pairwise') {
      html += '<div class="row">' + seqPicker('A', a.aSrc, a.aText) + seqPicker('B', a.bSrc, a.bText) + '</div>' +
        '<div class="row">' +
        '<div class="field narrow"><label>Mode</label><select class="inp" id="alMode">' +
          '<option value="global"' + (a.mode !== 'local' ? ' selected' : '') + '>Global</option>' +
          '<option value="local"' + (a.mode === 'local' ? ' selected' : '') + '>Local</option></select></div>' +
        '<div class="field narrow"><label>Match</label><input class="inp mono" id="alMatch" type="number" value="' + (a.match != null ? a.match : 2) + '"></div>' +
        '<div class="field narrow"><label>Mismatch</label><input class="inp mono" id="alMis" type="number" value="' + (a.mismatch != null ? a.mismatch : -1) + '"></div>' +
        '<div class="field narrow"><label>Gap</label><input class="inp mono" id="alGap" type="number" value="' + (a.gap != null ? a.gap : -2) + '"></div>' +
        '<div class="field" style="flex:0 0 auto;align-self:flex-end"><button class="btn primary" id="alGo">Align</button></div>' +
        '</div>';
      if (a.result) {
        const r = a.result;
        html += '<div class="row"><div class="field wide"><dl class="kv">' +
          '<dt>Identity</dt><dd>' + f1(r.identity) + '% (' + fmt(r.matches) + ' / ' + fmt(r.length) + ')</dd>' +
          '<dt>Score</dt><dd>' + fmt(r.score) + '</dd>' +
          '<dt>Gaps</dt><dd>' + fmt((r.a.match(/-/g) || []).length + (r.b.match(/-/g) || []).length) + '</dd>' +
          '</dl></div></div>' +
          '<div class="row"><div class="field wide">' + alignBlocks(r.a, r.b, 'A', 'B', 60, 0) + '</div></div>';
      }
    } else if (mode === 'msa') {
      html += '<div class="row"><div class="field wide"><label>Sequences — paste in FASTA, or pick from the library below</label>' +
        '<textarea class="inp" id="msaText" rows="7" spellcheck="false" placeholder=">seq1&#10;ACGT…&#10;>seq2&#10;ACGT…">' + esc(a.msaText || '') + '</textarea></div></div>' +
        '<div class="row"><div class="field"><label>Add from library</label><select class="inp" id="msaPick">' +
          '<option value="">choose…</option>' +
          state.records.map(r => '<option value="' + r.id + '">' + esc(shortName(r.name)) + '</option>').join('') +
        '</select></div>' +
        '<div class="field" style="flex:0 0 auto;align-self:flex-end"><button class="btn primary" id="msaGo">Align</button></div></div>';
      if (a.msaResult) {
        const r = a.msaResult;
        if (!r.ok) html += '<div class="note bad">' + esc(r.error) + '</div>';
        else {
          html += '<div class="note">' + esc(r.method) + ' · ' + f1(r.identity) + '% of columns fully conserved · ' +
            fmt(r.length) + ' columns</div>' + msaBlocks(r, 60) + figurePanel(a, r.rows);
        }
      }
    } else {
      html += '<div class="row"><div class="field wide"><label>Reads — one per line, or FASTA</label>' +
        '<textarea class="inp" id="asmText" rows="8" spellcheck="false" placeholder=">read1&#10;ACGT…">' + esc(a.asmText || '') + '</textarea></div></div>' +
        '<div class="row"><div class="field narrow"><label>Min overlap</label><input class="inp mono" id="asmOv" type="number" value="' + (a.minOverlap || 25) + '"></div>' +
        '<div class="field narrow"><label>Min identity %</label><input class="inp mono" id="asmId" type="number" value="' + (a.minIdentity || 95) + '"></div>' +
        '<div class="field" style="flex:0 0 auto;align-self:flex-end"><button class="btn primary" id="asmGo">Assemble</button></div>' +
        '<div class="field" style="flex:0 0 auto;align-self:flex-end"><button class="btn sm" id="asmFromReads">Use loaded reads</button></div>' +
        '</div>';
      if (a.asmResult) {
        const r = a.asmResult;
        html += '<div class="note' + (r.contigCount === 1 ? ' ok' : ' warn') + '">' +
          r.contigCount + ' contig' + (r.contigCount === 1 ? '' : 's') + ' from the reads provided.</div>' +
          '<div class="scrollx"><table class="tbl"><thead><tr><th>Contig</th><th class="num">Length</th><th>Reads</th><th></th></tr></thead><tbody>' +
          r.contigs.map((c, i) => '<tr><td>' + esc(c.name.slice(0, 40)) + '</td><td class="num">' + fmt(c.length) + '</td>' +
            '<td style="font-size:11.5px">' + esc(c.reads.join(', ')) + '</td>' +
            '<td><button class="btn sm ghost" data-asmsave="' + i + '">Save</button>' +
            '<button class="btn sm ghost" data-asmcopy="' + i + '">⧉</button></td></tr>').join('') +
          '</tbody></table></div>';
      }
    }
    html += '<div class="note">Pairwise alignment is an exact Needleman–Wunsch (global) or Smith–Waterman (local) dynamic program. Multiple alignment uses centre-star progressive alignment built on the same engine — it is deterministic and correct, but it is not Clustal Omega, MAFFT, MUSCLE or T-Coffee, which are separate compiled programs.</div>';
    html += '</div>';
    host.innerHTML = html;
    wireAlign(host);
  }
  PB.renderAlign = renderAlign;

  function seqPicker(which, src, text) {
    const key = which.toLowerCase();
    return '<div class="field" style="flex:1 1 300px"><label>Sequence ' + which + '</label>' +
      '<select class="inp" data-src="' + key + '" style="margin-bottom:5px">' +
      '<option value="">paste below…</option>' +
      state.records.map(r => '<option value="' + r.id + '"' + (src === r.id ? ' selected' : '') + '>' +
        esc(shortName(r.name)) + ' (' + fmt(r.seq.length) + ' bp)</option>').join('') +
      (state.selection ? '<option value="__sel"' + (src === '__sel' ? ' selected' : '') + '>current selection</option>' : '') +
      '</select>' +
      '<textarea class="inp" data-txt="' + key + '" rows="4" spellcheck="false" placeholder="ACGT…">' + esc(text || '') + '</textarea></div>';
  }

  function resolveSeq(src, text) {
    if (src === '__sel') {
      const rec = PB.active();
      if (rec && state.selection) return PB.selectionSeq(rec, state.selection);
      return '';
    }
    if (src) {
      const r = state.records.find(x => x.id === src);
      if (r) return r.seq;
    }
    return B.cleanSeq((text || '').replace(/^>.*$/gm, ''));
  }

  /* ---------------- publication figure of the alignment ----------------
     An alignment you read on screen and one you put in a paper are not the
     same drawing. The conventions of the second have been settled since
     ESPript: white letter on red for a strictly conserved column, red
     letters inside a blue frame for a conserved run, a structure track
     above, a consensus line below. The engine draws it; this is the panel
     that sets it up and hands the file over. */
  const FIG_MODE_LABEL = {
    physchem: 'Physico-chemical equivalence',
    multalin: 'MultAlin groups',
    strict: 'Strict identity',
    blosum62: 'BLOSUM62 similarity',
    purine: 'Purine / pyrimidine (transitions)'
  };
  function figurePanel(a, rows) {
    const f = a.fig = a.fig || {};
    /* protein or nucleotide decides which shading modes are even
       meaningful: the amino-acid groupings say nothing about a base
       column, and the purine grouping says nothing about a residue */
    /* Settle the molecule on the way in. Leaving it undefined meant the
       first change to ANY control saw it go from undefined to "auto",
       read that as a molecule change and rebuilt the whole panel — which
       detached every element a caller was holding and silently dropped
       that first change. */
    if (f.molecule === undefined) f.molecule = 'auto';
    let mol = f.molecule;
    const auto = rows.some(r => B.looksProteinSeq(String(r.seq).replace(/[-.]/g, ''))) ? 'protein' : 'dna';
    const kind = (mol === 'auto' ? auto : mol) === 'protein' ? 'protein' : 'nucleotide';
    const modes = B.ALN_MODES[kind].map(m => [m, FIG_MODE_LABEL[m]]);
    let mode = f.mode || (kind === 'protein' ? 'physchem' : 'strict');
    if (B.ALN_MODES[kind].indexOf(mode) === -1) mode = B.ALN_MODES[kind][0];
    const nuc = kind === 'nucleotide';
    const sel = (id, opts, val) => '<select class="inp" id="' + id + '">' + opts.map(([v, l]) =>
      '<option value="' + v + '"' + (String(v) === String(val) ? ' selected' : '') + '>' + l + '</option>').join('') + '</select>';
    return '<details class="fig-wrap"' + (f.open ? ' open' : '') + ' id="figDet">' +
      '<summary><strong>Publication figure</strong> — shaded by conservation, with a structure track</summary>' +
      '<div class="row" style="margin-top:8px">' +
      '<div class="field"><label>Sequence</label>' + sel('figMol', [
        ['auto', 'Auto — read as ' + (auto === 'protein' ? 'protein' : 'DNA')],
        ['protein', 'Protein'], ['dna', 'DNA'], ['rna', 'RNA (T shown as U)']], mol) + '</div>' +
      '<div class="field"><label>Shade by</label>' + sel('figMode', modes, mode) + '</div>' +
      '<div class="field narrow"><label>Threshold</label><input class="inp mono" id="figThr" type="number" step="0.05" min="0" max="1" value="' +
        (f.threshold != null ? f.threshold : 0.7) + '"></div>' +
      '<div class="field narrow"><label>Columns a line</label><input class="inp mono" id="figCols" type="number" min="10" max="200" value="' +
        (f.columns || 60) + '"></div>' +
      '<div class="field"><label>Title</label><input class="inp" id="figTitle" value="' + esc(f.title || '') + '" placeholder="optional"></div>' +
      '</div><div class="row">' +
      '<label class="chk"><input type="checkbox" id="figCons"' + (f.consensus === false ? '' : ' checked') + '> Consensus line</label>' +
      '<label class="chk"><input type="checkbox" id="figRuler"' + (f.ruler === false ? '' : ' checked') + '> Numbering</label>' +
      '<label class="chk"><input type="checkbox" id="figLegend"' + (f.legend === false ? '' : ' checked') + '> Key</label>' +
      '<label class="chk"><input type="checkbox" id="figMono"' + (f.mono ? ' checked' : '') + '> Monochrome</label>' +
      (nuc ? '<label class="chk"><input type="checkbox" id="figNuc"' + (f.nucleotideColour ? ' checked' : '') +
        '> Purines red, pyrimidines blue</label>' : '') +
      '</div>' +
      '<div class="row"><div class="field wide"><label>' + (nuc
        ? 'Structure over the first sequence — dot-bracket, e.g. <span class="mono">((((....))))</span>'
        : 'Secondary structure over the first sequence — a DSSP string, or paste PSIPRED or JPred output as it comes') +
      '</label><textarea class="inp mono" id="figSS" rows="2" spellcheck="false" placeholder="' +
      (nuc ? '((((((((....))))))))' : 'HHHHHHHH  EEEEEE TT EEEE     — or the whole Conf:/Pred:/AA: block from PSIPRED') +
      '">' + esc(f.ss || '') + '</textarea>' +
      '<div class="row" style="gap:6px;margin-top:6px">' +
      '<button class="btn sm ghost" id="figClearSS">Clear</button>' +
      '<span class="hint" style="margin:0;align-self:center" id="figSSNote"></span></div>' +
      '</div></div>' +
      '<div class="hint">' + (nuc
        ? 'For nucleotides the figure follows N·ESPript: strict identity, or purine–pyrimidine equivalence so that a transition (A↔G, C↔T) counts as conserved and a transversion does not. ' +
          'Structure is given in dot-bracket notation rather than DSSP letters, and the chemical-property colouring replaces the flat red of a conserved column rather than sitting on top of it.'
        : 'The three grouping modes are ESPript’s own, implemented from its user guide: strict identity, the MultAlin groups ' +
          '(IV, LM, FY, NDQEBZ), and equivalence by physico-chemical property (HKR, DE, STNQ, AVLIM, FYW, PGC). ' +
          'ESPript’s default Risler global score is not reproduced, because its formula is not published; BLOSUM62 similarity ' +
          '— the share of residue pairs in a column scoring above zero, what a BLAST alignment marks with a + — stands in its place.') +
      '</div>' +
      '<div class="row"><button class="btn" id="figSvg">Save SVG</button>' +
      '<button class="btn" id="figPng">Save PNG</button>' +
      '<button class="btn ghost" id="figCopy">Copy SVG</button></div>' +
      '<div class="fig-prev" id="figPrev"></div>' +
      '</details>';
  }
  PB.figurePanel = figurePanel;

  function figureOpts(a) {
    const f = a.fig || {};
    return {
      molecule: f.molecule || 'auto',
      mode: f.mode || 'physchem',
      threshold: f.threshold != null ? f.threshold : 0.7,
      columns: f.columns || 60,
      consensus: f.consensus !== false,
      ruler: f.ruler !== false,
      mono: !!f.mono,
      legend: f.legend !== false,
      nucleotideColour: !!f.nucleotideColour,
      ss: f.ss || '',
      title: f.title || ''
    };
  }

  function msaBlocks(r, width) {
    const nameW = Math.min(16, Math.max.apply(null, r.rows.map(x => x.name.length).concat([4])));
    let out = '';
    for (let i = 0; i < r.length; i += width) {
      r.rows.forEach(row => {
        const seg = row.seq.substr(i, width);
        let painted = '';
        for (let k = 0; k < seg.length; k++) {
          const col = i + k;
          const allSame = r.rows.every(x => x.seq[col] === r.rows[0].seq[col]) && seg[k] !== '-';
          painted += allSame ? esc(seg[k]) :
            (seg[k] === '-' ? '<span class="gp">-</span>' : '<span class="mm">' + esc(seg[k]) + '</span>');
        }
        out += '<span class="lbl">' + esc((row.name + '                ').slice(0, nameW)) + ' </span>' + painted + '\n';
      });
      out += '<span class="lbl">' + ' '.repeat(nameW + 1) + '</span><span class="mid">' +
        esc(r.consensus.substr(i, width)) + '</span>\n\n';
    }
    return '<div class="aln-box">' + out + '</div>';
  }

  /* The figure panel: every control writes to state and redraws the
     preview in place, so nothing re-renders the whole view and loses the
     open/closed state of the panel. */
  function wireFigure(host, st) {
    const det = host.querySelector('#figDet');
    if (!det || !st.msaResult || !st.msaResult.ok) return;
    const a = st, f = a.fig = a.fig || {};
    det.addEventListener('toggle', () => { f.open = det.open; });
    const prev = host.querySelector('#figPrev');
    const build = () => B.alnFigureSVG(st.msaResult.rows, figureOpts(a));
    const draw = () => {
      try {
        const svg = build();
        prev.innerHTML = '<img alt="Preview of the alignment figure" src="data:image/svg+xml;charset=utf-8,' +
          encodeURIComponent(svg) + '">';
      } catch (e) { prev.textContent = e.message; }
    };
    const read = () => {
      const molEl = host.querySelector('#figMol');
      const wasMol = f.molecule;
      if (molEl) f.molecule = molEl.value;
      f.mode = host.querySelector('#figMode').value;
      f.threshold = Math.max(0, Math.min(1, parseFloat(host.querySelector('#figThr').value)));
      if (!isFinite(f.threshold)) f.threshold = 0.7;
      f.columns = Math.max(10, Math.min(200, parseInt(host.querySelector('#figCols').value, 10) || 60));
      f.title = host.querySelector('#figTitle').value;
      f.consensus = host.querySelector('#figCons').checked;
      f.ruler = host.querySelector('#figRuler').checked;
      f.mono = host.querySelector('#figMono').checked;
      const legEl = host.querySelector('#figLegend');
      if (legEl) f.legend = legEl.checked;
      const nucEl = host.querySelector('#figNuc');
      if (nucEl) f.nucleotideColour = nucEl.checked;
      f.ss = host.querySelector('#figSS').value;
      /* say what the pasted text was understood as, or that it was not */
      const note = host.querySelector('#figSSNote');
      if (note) {
        const got = B.parseSSInput(f.ss);
        note.textContent = !f.ss.trim() ? 'No structure track. A DSSP string, PSIPRED or JPred output, or dot-bracket all work.'
          : got ? ('read as ' + got.source + ' \u2014 ' + got.text.replace(/\s/g, '').length + ' assigned positions')
                : 'not recognised \u2014 nothing is drawn. Expected DSSP letters, H/E/C, or dot-bracket.';
      }
      /* changing protein to DNA changes which shading modes exist and
         which structure notation is asked for, so the panel is rebuilt */
      if (wasMol !== f.molecule) { f.open = true; renderAlign(); return; }
      draw();
    };
    const clearBtn = host.querySelector('#figClearSS');
    if (clearBtn) clearBtn.addEventListener('click', () => {
      f.ss = '';
      host.querySelector('#figSS').value = '';
      draw();
    });
    ['#figMol', '#figMode', '#figThr', '#figCols', '#figTitle', '#figCons', '#figRuler', '#figLegend', '#figMono', '#figNuc', '#figSS']
      .forEach(id => {
        const el = host.querySelector(id);
        if (!el) return;
        el.addEventListener('change', read);
        el.addEventListener('input', read);
      });
    const name = () => (f.title || 'alignment').replace(/[^\w.-]+/g, '_').slice(0, 50) + '_figure';
    const svgBtn = host.querySelector('#figSvg');
    if (svgBtn) svgBtn.addEventListener('click', async () => {
      let svg; try { svg = build(); } catch (e) { PB.toast(e.message, 'bad'); return; }
      await PB.offerFile(name() + '.svg', svg, 'Alignment figure (SVG)', {
        alt: { label: 'Save it as a PNG instead', run: () => savePng() }
      });
    });
    const savePng = async () => {
      let svg; try { svg = build(); } catch (e) { PB.toast(e.message, 'bad'); return; }
      let blob;
      try { blob = await PB.mapPNGBlob(svg, 3); } catch (e) { PB.toast(e.message, 'bad'); return; }
      await PB.offerBlob(name() + '.png', blob, 'Alignment figure (PNG)');
    };
    const pngBtn = host.querySelector('#figPng');
    if (pngBtn) pngBtn.addEventListener('click', savePng);
    const copyBtn = host.querySelector('#figCopy');
    if (copyBtn) copyBtn.addEventListener('click', () => {
      let svg; try { svg = build(); } catch (e) { PB.toast(e.message, 'bad'); return; }
      PB.copyText(svg, 'Alignment figure (SVG)');
    });
    draw();
  }

  function wireAlign(host) {
    host.querySelectorAll('[data-am]').forEach(b => b.addEventListener('click', () => {
      state.alignMode = b.getAttribute('data-am'); renderAlign();
    }));
    const st = state.alignState = state.alignState || {};
    host.querySelectorAll('[data-src]').forEach(s => s.addEventListener('change', () => {
      st[s.getAttribute('data-src') + 'Src'] = s.value; renderAlign();
    }));
    host.querySelectorAll('[data-txt]').forEach(t => t.addEventListener('input', () => {
      st[t.getAttribute('data-txt') + 'Text'] = t.value;
    }));
    const go = $('#alGo');
    if (go) go.addEventListener('click', () => {
      const A = resolveSeq(st.aSrc, st.aText), Bq = resolveSeq(st.bSrc, st.bText);
      if (A.length < 2 || Bq.length < 2) { PB.toast('Provide two sequences.', 'warn'); return; }
      if (A.length * Bq.length > 36e6) { PB.toast('Those sequences are too long to align in the browser (limit about 6 kb × 6 kb).', 'bad'); return; }
      st.mode = $('#alMode').value;
      st.match = parseFloat($('#alMatch').value); st.mismatch = parseFloat($('#alMis').value);
      st.gap = parseFloat($('#alGap').value);
      st.result = B.align(A, Bq, { mode: st.mode, match: st.match, mismatch: st.mismatch, gap: st.gap });
      renderAlign();
    });
    const mp = $('#msaPick');
    if (mp) mp.addEventListener('change', () => {
      const r = state.records.find(x => x.id === mp.value);
      if (!r) return;
      const ta = $('#msaText');
      ta.value = (ta.value ? ta.value.replace(/\s*$/, '') + '\n' : '') + '>' + r.name.replace(/\s+/g, '_') + '\n' + r.seq + '\n';
      st.msaText = ta.value;
      mp.value = '';
    });
    const mgo = $('#msaGo');
    if (mgo) mgo.addEventListener('click', () => {
      st.msaText = $('#msaText').value;
      const p = parseAlnFasta(st.msaText);
      if (p.rows.length < 2) { PB.toast('Provide at least two sequences in FASTA form.', 'warn'); return; }
      /* Sequences that arrive already aligned — same length, with gaps —
         are drawn as they came. Re-aligning somebody's MAFFT output with
         a centre-star progressive aligner would throw away the better
         alignment and silently replace it with a worse one. */
      if (p.aligned) {
        const seqs = p.rows.map(r => r.seq);
        st.msaResult = {
          ok: true, rows: p.rows, length: seqs[0].length,
          consensus: B.alnConsensus(Array.from({ length: seqs[0].length }, (_, i) => seqs.map(s => s[i]))),
          identity: B.alnAnalyse(p.rows, {}).identity,
          method: 'taken as already aligned (' + (p.protein ? 'protein' : 'nucleotide') + ', ' +
                  p.rows.length + ' sequences, gaps kept as given)'
        };
        st.msaProtein = p.protein;
        renderAlign();
        return;
      }
      const tot = p.rows.reduce((a, r) => a + r.seq.length, 0);
      if (tot > 30000) { PB.toast('Too much sequence for in-browser multiple alignment (about 30 kb total).', 'bad'); return; }
      /* protein scored +2/-1 is a poor alignment; BLOSUM62 is used instead */
      const clean = p.rows.map(r => ({ name: r.name, seq: r.seq.replace(/[-.]/g, '') }));
      st.msaResult = p.protein ? B.msaProtein(clean) : B.msa(clean);
      if (st.msaResult.ok && p.protein) st.msaResult.method += ', BLOSUM62';
      st.msaProtein = p.protein;
      renderAlign();
    });
    wireFigure(host, st);
    const ago = $('#asmGo');
    if (ago) ago.addEventListener('click', () => {
      st.asmText = $('#asmText').value;
      st.minOverlap = parseInt($('#asmOv').value, 10) || 25;
      st.minIdentity = parseFloat($('#asmId').value) || 95;
      const reads = parseFastaish(st.asmText);
      if (reads.length < 2) { PB.toast('Provide at least two reads.', 'warn'); return; }
      st.asmResult = B.assembleContigs(reads, { minOverlap: st.minOverlap, minIdentity: st.minIdentity });
      renderAlign();
    });
    const afr = $('#asmFromReads');
    if (afr) afr.addEventListener('click', () => {
      if (!state.reads || !state.reads.length) { PB.toast('No reads loaded — add them in the Verify tab.', 'warn'); return; }
      const txt = state.reads.map(r => '>' + r.name.replace(/\s+/g, '_') + '\n' + (r.trim ? r.trim.seq : r.seq)).join('\n');
      $('#asmText').value = txt; st.asmText = txt;
    });
    host.querySelectorAll('[data-asmsave]').forEach(b => b.addEventListener('click', () => {
      const c = st.asmResult.contigs[+b.getAttribute('data-asmsave')];
      const nr = PB.makeRecord({ name: 'contig (' + fmt(c.length) + ' bp)', seq: c.seq, topology: 'linear' });
      PB.addHistory(nr, 'Contig assembled', 'From reads: ' + c.reads.join(', '));
      PB.addRecord(nr);
    }));
    host.querySelectorAll('[data-asmcopy]').forEach(b => b.addEventListener('click', () =>
      PB.copyText(st.asmResult.contigs[+b.getAttribute('data-asmcopy')].seq, 'Contig')));
  }

  function parseFastaish(text) {
    text = text || '';
    if (/^\s*>/.test(text)) {
      return B.parseFASTA(text).filter(r => r.seq.length).map((r, i) => ({ name: r.name || ('seq' + (i + 1)), seq: r.seq }));
    }
    return text.split(/\n+/).map(l => B.cleanSeq(l)).filter(s => s.length)
      .map((s, i) => ({ name: 'seq' + (i + 1), seq: s }));
  }
  PB.parseFastaish = parseFastaish;

  /* The alignment view needs something parseFastaish cannot give it. That
     one runs every sequence through cleanSeq, which keeps only nucleotide
     codes — so a protein arrives as the handful of its residues that
     happen to be IUPAC bases, and an already-aligned file loses the gaps
     that are the whole point of it. This keeps both, works out whether
     the set is protein or nucleotide, and says whether it arrived
     aligned. */
  function parseAlnFasta(text) {
    text = String(text || '');
    const rows = [];
    if (/^\s*>/.test(text)) {
      text.replace(/\r\n?/g, '\n').split(/^>/m).forEach(block => {
        if (!block.trim()) return;
        const nl = block.indexOf('\n');
        const header = (nl === -1 ? block : block.slice(0, nl)).trim();
        const body = nl === -1 ? '' : block.slice(nl + 1);
        const seq = body.toUpperCase().replace(/[^A-Z*.\-]/g, '');
        if (seq) rows.push({ name: header.split(/\s+/)[0] || ('seq' + (rows.length + 1)), seq: seq });
      });
    } else {
      text.split(/\n+/).forEach(line => {
        const seq = line.toUpperCase().replace(/[^A-Z*.\-]/g, '');
        if (seq) rows.push({ name: 'seq' + (rows.length + 1), seq: seq });
      });
    }
    const protein = rows.some(r => B.looksProteinSeq(r.seq));
    const widths = rows.map(r => r.seq.length);
    const aligned = rows.length > 1 && widths.every(w => w === widths[0]) &&
      rows.some(r => /[-.]/.test(r.seq));
    return { rows: rows, protein: protein, aligned: aligned };
  }
  PB.parseAlnFasta = parseAlnFasta;

  /* ============================== HISTORY ============================== */
  function renderHistory() {
    const host = $('#view-history');
    const rec = PB.active();
    if (!rec) { host.innerHTML = emptyState('No construct loaded.'); return; }
    const h = (rec.history || []).slice().reverse();
    let html = '<div class="panel"><h3>History <span class="sub">' + h.length + ' recorded step' + (h.length === 1 ? '' : 's') + ' for “' + esc(rec.name) + '”</span>' +
      (h.length ? '<button class="btn sm" id="hCopy">Copy as text</button>' : '') + '</h3>';
    if (!h.length) {
      html += emptyState('Nothing recorded yet. Every clone, digest, mutagenesis and annotation step you run from the Tools menu is logged here.');
    } else {
      html += '<div class="hist">' + h.map((e, i) =>
        '<div class="hist-item"><div class="hist-dot' + (i === h.length - 1 ? ' src' : '') + '"></div><div class="hist-body">' +
        '<div class="hist-op">' + esc(e.op) + '</div>' +
        (e.detail ? '<div class="hist-detail">' + esc(e.detail) + '</div>' : '') +
        '<div class="hist-when">' + new Date(e.at).toLocaleString() + '</div>' +
        '</div></div>').join('') + '</div>';
    }
    html += '<dl class="kv" style="margin-top:14px;padding-top:12px;border-top:1px solid var(--border)">' +
      '<dt>Created</dt><dd>' + new Date(rec.createdAt).toLocaleString() + '</dd>' +
      '<dt>Last change</dt><dd>' + new Date(rec.updatedAt).toLocaleString() + '</dd>' +
      '<dt>Identifier</dt><dd>' + esc(rec.id) + '</dd>' +
      '</dl>';
    if (rec.description) {
      html += '<div class="note" style="margin-top:12px">' + esc(rec.description) + '</div>';
    }
    html += '</div>';
    host.innerHTML = html;
    const c = $('#hCopy');
    if (c) c.addEventListener('click', () => PB.copyText(
      rec.name + '\n' + (rec.history || []).map(e =>
        new Date(e.at).toISOString() + '\t' + e.op + '\t' + e.detail).join('\n'), 'History'));
  }
  PB.renderHistory = renderHistory;

  /* ---------------- view dispatcher ---------------- */
  PB.renderViews = function () {
    const v = state.view;
    try {
      if (v === 'map') PB.renderMap();
      else if (v === 'seq') PB.renderSeq();
      else if (v === 'features') PB.renderFeatures();
      else if (v === 'enzymes') PB.renderEnzymes();
      else if (v === 'orfs') PB.renderORFs();
      else if (v === 'primers') PB.renderPrimers();
      else if (v === 'gel') PB.renderGel();
      else if (v === 'verify') PB.renderVerify();
      else if (v === 'align') PB.renderAlign();
      else if (v === 'history') PB.renderHistory();
    } catch (e) {
      const host = document.getElementById('view-' + v);
      if (host) host.innerHTML = '<div class="note bad"><strong>That view hit an error.</strong><br>' +
        esc(e && e.message ? e.message : String(e)) + '</div>';
      if (window.console) console.error(e);
    }
  };

  /* ======================= PRIMERS, READY TO ORDER =======================
     Everything measured about each primer in one file. The melting
     temperature depends on the salt it is calculated for, so the
     conditions travel with the numbers rather than being left implicit in
     a column heading — a Tm without them is not a figure anyone can
     check or reproduce.
  ------------------------------------------------------------------ */
  function primerRows(which) {
    const rec = PB.active();
    const cond = PB.conditions();
    const c = state.tmCond;
    const condText = c.oligo + ' µM oligo, ' + c.Na + ' mM Na+, ' + c.Mg + ' mM Mg2+, ' + c.dNTP + ' mM dNTP';
    const list = which === 'shown'
      ? state.primers.filter(p => PB.primerShown(p.name))
      : state.primers;
    const hits = rec ? PB.primerHits(rec) : [];
    return list.map(p => {
      const seq = p.seq || '';
      /* The exported Tm is the one you would set the block to: for a primer
         with an overhang that is the annealing part's, not the whole
         primer's. Both go out, with the annealing length beside them, so a
         row read away from this program still says which is which. */
      const T = PB.primerTmOf(p, rec, hits);
      const tm = T.tm;
      const hp = B.hairpin(seq);
      const sd = B.selfDimer(seq);
      const mine = hits.filter(h => h.name === p.name);
      return {
        name: p.name, seq: seq, len: seq.length,
        tm: isFinite(tm) ? tm : null,
        tmFull: isFinite(T.tmFull) ? T.tmFull : null,
        annealLen: T.annealLen,
        tailLen: T.tailLen,
        tmBasis: T.tailed ? ('anneals over ' + T.annealLen + ' of ' + seq.length +
                             ' nt (' + (T.basis === 'design' ? 'as designed' : 'against ' + (rec ? rec.name : 'the construct')) + ')')
                          : 'whole primer',
        gc: B.gcContent(seq),
        hairpin: hp.stemLength >= 4 ? hp.stemLength : 0,
        hairpinDG: hp.stemLength >= 4 ? hp.dG : 0,
        dimer: sd.matches, dimer3: !!sd.threePrime,
        sites: mine.length,
        where: mine.slice(0, 4).map(h => (h.strand === 1 ? '+' : '-') + (h.start + 1) + '..' + h.end).join(' '),
        tail: mine.length ? (mine[0].tail || 0) : 0,
        cond: condText,
        on: rec ? rec.name : ''
      };
    });
  }
  function csvEsc(v) {
    const t = String(v == null ? '' : v);
    return /[",\n]/.test(t) ? '"' + t.replace(/"/g, '""') + '"' : t;
  }
  function primerCsv(rows, kind) {
    if (kind === 'order') {
      /* the four columns a synthesis order asks for; scale follows length,
         and desalted is right for anything used in a PCR or an assembly */
      const head = ['Name', 'Sequence', 'Scale', 'Purification'];
      return [head.join(',')].concat(rows.map(r => [
        csvEsc(r.name), csvEsc(r.seq),
        r.len > 60 ? '100nm' : '25nm',
        'STD'
      ].join(','))).join('\n') + '\n';
    }
    if (kind === 'fasta') {
      return rows.map(r => '>' + r.name + '\n' + r.seq).join('\n') + '\n';
    }
    const head = ['name', 'sequence', 'length_nt', 'tm_c', 'gc_percent', 'hairpin_stem_bp', 'hairpin_dg',
                  'self_dimer_bp', 'self_dimer_3prime', 'sites_on_construct', 'where', 'tail_nt',
                  'construct', 'tm_conditions'];
    return [head.join(',')].concat(rows.map(r => [
      csvEsc(r.name), csvEsc(r.seq), r.len,
      r.tm == null ? '' : r.tm.toFixed(1),
      r.gc.toFixed(1),
      r.hairpin || '', r.hairpin ? r.hairpinDG.toFixed(2) : '',
      r.dimer, r.dimer3 ? 'yes' : 'no',
      r.sites, csvEsc(r.where), r.tail || '',
      csvEsc(r.on), csvEsc(r.cond)
    ].join(','))).join('\n') + '\n';
  }
  /* the workbook: the table on one sheet, and what the numbers mean and
     what they were computed under on another, so the file explains itself
     to whoever opens it next */
  PB.primerExportRows = primerRows;

  function primerSheets(rows) {
    const c = state.tmCond;
    const r1 = n => n == null ? '' : Math.round(n * 10) / 10;
    const head = ['Name', 'Sequence', 'Length (nt)', 'Tm anneal (°C)', 'Tm whole primer (°C)',
                  'Anneals over (nt)', 'Overhang (nt)', 'Tm measured on', 'GC (%)', 'Hairpin stem (bp)',
                  'Hairpin dG', 'Self-dimer (bp)', '3-prime end involved', 'Sites here', 'Where', 'Tail (nt)'];
    const body = rows.map(r => [
      r.name, r.seq, r.len,
      r1(r.tm), r1(r.tmFull),
      r.annealLen == null ? '' : r.annealLen,
      r.tailLen == null ? '' : r.tailLen,
      r.tmBasis || '',
      Math.round(r.gc * 10) / 10,
      r.hairpin || 0,
      r.hairpin ? Math.round(r.hairpinDG * 100) / 100 : 0,
      r.dimer, r.dimer3 ? 'yes' : 'no',
      r.sites, r.where || '', r.tail || 0
    ]);
    const about = [
      ['Item', 'Value'],
      ['Construct open when exported', rows.length && rows[0].on ? rows[0].on : ''],
      ['Primers written', rows.length],
      ['', ''],
      ['Melting temperature', 'Nearest-neighbour (SantaLucia 1998) with the Owczarzy monovalent and Mg corrections'],
      ['Primers with an overhang', 'Only the 3-prime part of such a primer anneals in the first cycles, so ' +
        '"Tm anneal" is measured on that part alone and is the temperature to set the block to. ' +
        '"Tm whole primer" is what later cycles see, once the overhang has been copied into the product. ' +
        'For a primer with no overhang the two are the same number.'],
      ['Oligo concentration', c.oligo + ' µM'],
      ['Na+', c.Na + ' mM'],
      ['Mg2+', c.Mg + ' mM'],
      ['dNTP', c.dNTP + ' mM'],
      ['', ''],
      ['Hairpin stem', 'Longest self-complementary stem, in base pairs; 0 below 4 bp'],
      ['Hairpin dG', 'Free energy of that stem at 37 C, kcal/mol — more negative is more stable'],
      ['Self-dimer', 'Longest run where the primer pairs with another copy of itself'],
      ['3-prime end involved', 'Whether that run reaches the 3-prime end, which is the end that extends'],
      ['Sites here', 'How many places it anneals on the construct named above'],
      ['Where', 'Those positions; + is the top strand, - the bottom'],
      ['Tail (nt)', 'Bases at the 5-prime end that do not match the template — a homology arm or an added site']
    ];
    return [
      { name: 'Primers', widths: [18, 46, 11, 9, 9, 15, 11, 14, 18, 10, 22, 10],
        rows: [head].concat(body) },
      { name: 'About', widths: [32, 74], rows: about }
    ];
  }
  PB.primerExportDialog = function () {
    if (!state.primers.length) { PB.toast('No primers stored yet.', 'warn'); return; }
    const shownCount = state.primers.filter(p => PB.primerShown(p.name)).length;
    const w = PB.modal({
      title: 'Export primers',
      sub: state.primers.length + ' in the library',
      wide: true,
      body:
        '<div class="row">' +
          '<div class="field" style="flex:1 1 320px"><label>What to write</label><select class="inp" id="peKind">' +
            '<option value="xlsx">A spreadsheet — everything measured, in columns, with the conditions on a second sheet</option>' +
            '<option value="full">The same as a CSV file</option>' +
            '<option value="order">An ordering list — name, sequence, scale, purification</option>' +
            '<option value="fasta">FASTA — names and sequences only</option>' +
          '</select></div>' +
          '<div class="field" style="flex:1 1 200px"><label>Which primers</label><select class="inp" id="peWhich">' +
            '<option value="all">All ' + state.primers.length + '</option>' +
            '<option value="shown">The ' + shownCount + ' shown on the map</option>' +
          '</select></div>' +
        '</div>' +
        '<div class="row"><div class="field wide"><label>Preview</label>' +
          '<textarea class="inp mono" id="pePrev" rows="11" spellcheck="false" readonly></textarea></div></div>' +
        '<div class="hint" id="peHint"></div>',
      actions: [
        { label: 'Copy', close: false, onClick: () => {
            if (w.querySelector('#peKind').value === 'xlsx') {
              PB.toast('A spreadsheet cannot be copied as text — use Download.', 'warn'); return;
            }
            PB.copyText(w.querySelector('#pePrev').value, 'Primer table');
          } },
        { label: 'Download', primary: true, close: false, onClick: () => {
            const kind = w.querySelector('#peKind').value;
            const rows = primerRows(w.querySelector('#peWhich').value);
            if (kind === 'xlsx') {
              PB.offerBlob('primers.xlsx', PB.xlsxBlob(primerSheets(rows)), 'Primers');
              return;
            }
            const text = w.querySelector('#pePrev').value;
            const ext = kind === 'fasta' ? '.fasta' : '.csv';
            const base = kind === 'order' ? 'primer-order' : (kind === 'fasta' ? 'primers' : 'primers-full');
            PB.offerBlob(base + ext, new Blob([text], { type: kind === 'fasta' ? 'text/plain' : 'text/csv' }),
                         'Primers');
          } }
      ]
    });
    const build = () => {
      const kind = w.querySelector('#peKind').value;
      const rows = primerRows(w.querySelector('#peWhich').value);
      w.querySelector('#pePrev').value = kind === 'xlsx'
        ? primerSheets(rows)[0].rows.map(r => r.join('\t')).join('\n') + '\n'
        : primerCsv(rows, kind);
      const c = state.tmCond;
      w.querySelector('#peHint').innerHTML = kind === 'order'
        ? 'Scale is 25 nmol up to 60 bases and 100 nmol above that; desalted is the usual purification for a primer ' +
          'used in a PCR or an assembly. Change either in the file if your supplier wants something else — ' +
          'the column names match what the bulk-entry boxes expect.'
        : (kind === 'fasta' ? 'Names and sequences, nothing else.'
          : 'Melting temperatures are nearest-neighbour values at ' + c.oligo + ' µM oligo, ' + c.Na +
            ' mM Na⁺, ' + c.Mg + ' mM Mg²⁺ and ' + c.dNTP + ' mM dNTP — the conditions set in the ' +
            'side panel, and written into the file so the numbers can be checked later. ' +
            '“Where” is where each primer anneals on the construct that is open.');
    };
    w.querySelector('#peKind').addEventListener('change', build);
    w.querySelector('#peWhich').addEventListener('change', build);
    build();
  };

})(window.PB);

