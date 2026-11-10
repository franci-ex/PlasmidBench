/* ============================================================================
   PlasmidBench — views: map, sequence, features, enzymes, ORFs
   ============================================================================ */
(function (PB) {
  'use strict';
  const B = window.Bio;
  const $ = PB.$, $$ = PB.$$, esc = PB.esc, fmt = PB.fmt, f1 = PB.f1;
  const state = PB.state;

  function emptyState(msg, action) {
    return '<div class="empty">' +
      '<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.4"><circle cx="12" cy="12" r="9"/><path d="M12 7v5l3 2"/></svg>' +
      '<div>' + msg + '</div>' + (action || '') + '</div>';
  }
  PB.emptyState = emptyState;

  /* ---------------- shared geometry helpers ---------------- */
  // Every feature as drawable geometry: its pieces (non-wrapping, top-strand
  // order), its spliced length, its full span (start..end along the circle)
  // and a display track chosen so no two features' pieces overlap.
  function featIntervals(rec) {
    const L = rec.seq.length;
    const out = rec.features.map((f, i) => {
      const pieces = f.between ? [{ start: f.start, end: Math.min(L, f.start + 1) }] : B.featurePieces(f, L);
      const bpLen = f.between ? 0 : B.featureLength(f, L);
      const span = f.between ? 1 : (f.end > f.start ? f.end - f.start : (L - f.start) + f.end);
      return { idx: i, f: f, pieces: pieces, bpLen: bpLen,
               start: f.start, len: Math.max(span, 1), end: f.start + Math.max(span, 1) };
    });
    out.sort((a, b) => a.start - b.start || b.len - a.len);
    // occupied ranges include the whole span (introns too) so joined features
    // read as one unit; a wrapped span also occupies [0, end) at the origin
    const occ = [];
    for (const it of out) {
      const ranges = it.end <= L ? [[it.start, it.end]] : [[it.start, L], [0, it.end - L]];
      let t = 0;
      while ((occ[t] || []).some(r => ranges.some(q => q[0] < r[1] && r[0] < q[1]))) t++;
      it.track = t;
      (occ[t] = occ[t] || []).push.apply(occ[t], ranges);
    }
    return out;
  }
  PB.featIntervals = featIntervals;

  /* niceStep and layoutLabels used to live here, and were still here after
     the map viewer moved into app-map.js, which kept the copies that are
     called. Nothing in this file has called either since. layoutLabels in
     particular had gone stale: app-map.js grew the version that fans a
     crowded side into several columns — which is what makes the cloning
     site of a plasmid readable — while this copy stayed single-column. Two
     copies of a drawing routine, one of them improved, is precisely how a
     fix comes to apply in one view and not the other, so the dead pair is
     gone rather than kept in step.
     Removed 2026-11; both are in app-map.js, which is the only caller. */

  /* ============================== MAP ============================== */
  /* The map viewer lives in its own module (app-map.js) — it owns
     renderMap, the semantic-zoom linear view and the circular view. */

  /* ============================== SEQUENCE ============================== */
  const SEQ_PAGE = 30000;

  function renderSeq() {
    const host = $('#view-seq');
    const rec = PB.active();
    if (!rec) { host.innerHTML = emptyState('No construct loaded.'); return; }
    const o = state.seqOpts;
    const L = rec.seq.length;
    const pages = Math.max(1, Math.ceil(L / SEQ_PAGE));
    if (state.seqPage == null || state.seqPage >= pages) state.seqPage = 0;

    host.innerHTML =
      '<div class="seq-toolbar">' +
        '<div class="field narrow" style="flex:0 0 86px"><label>bp / line</label>' +
        '<select class="inp" id="sqBpl">' + [30, 60, 90, 120].map(n =>
          '<option value="' + n + '"' + (o.bpl === n ? ' selected' : '') + '>' + n + '</option>').join('') + '</select></div>' +
        '<label class="chip' + (o.bottom ? ' on' : '') + '" data-tg="bottom">Complement</label>' +
        '<label class="chip' + (o.aa ? ' on' : '') + '" data-tg="aa">Translation</label>' +
        '<label class="chip' + (o.colorBases ? ' on' : '') + '" data-tg="colorBases">Colour bases</label>' +
        '<span style="width:6px"></span>' +
        '<div class="field" style="flex:0 1 190px"><label>Find</label>' +
        '<input class="inp mono" id="sqFind" placeholder="ACGT… or enzyme" spellcheck="false"></div>' +
        '<button class="btn sm" id="sqFindBtn">Find next</button>' +
        '<button class="btn sm ghost" id="sqFindAll" title="Search every construct, peptides and labels too">Search library</button>' +
        '<div class="field narrow" style="flex:0 0 94px"><label>Go to bp</label>' +
        '<input class="inp mono" id="sqGoto" type="number" min="1" max="' + L + '"></div>' +
        '<span class="spacer" style="flex:1 1 auto"></span>' +
        (pages > 1 ? '<button class="btn sm" id="sqPrev"' + (state.seqPage === 0 ? ' disabled' : '') + '>&larr;</button>' +
          '<span class="mono" style="font-size:11px">' + fmt(state.seqPage * SEQ_PAGE + 1) + '–' +
          fmt(Math.min(L, (state.seqPage + 1) * SEQ_PAGE)) + '</span>' +
          '<button class="btn sm" id="sqNext"' + (state.seqPage >= pages - 1 ? ' disabled' : '') + '>&rarr;</button>' : '') +
        '<button class="btn sm" id="sqCopy">Copy</button>' +
      '</div>' +
      '<div class="seq-scroll" id="seqScroll"><div class="seq-inner" id="seqInner"></div></div>' +
      '<div class="seq-status" id="seqStatus"></div>';

    renderSeqLines(rec);
    wireSeq(host, rec);
    updateSeqStatus(rec);
  }
  PB.renderSeq = renderSeq;

  function renderSeqLines(rec) {
    const inner = $('#seqInner');
    if (!inner) return;
    const o = state.seqOpts;
    const L = rec.seq.length;
    if (!L) { inner.innerHTML = '<div class="hint">No sequence. Use Tools → Replace the whole sequence to paste one.</div>'; return; }
    const from = state.seqPage * SEQ_PAGE;
    const to = Math.min(L, from + SEQ_PAGE);
    const bpl = o.bpl;
    const colorOn = o.colorBases && (to - from) <= 12000;
    if (o.colorBases && !colorOn) {
      PB.toast('Base colouring is off above 12,000 bp on screen for speed.', 'warn');
    }
    const ivs = featIntervals(PB.displayRecord(rec));
    const cuts = state.enzymes.length ? PB.allCuts(rec) : [];
    const cutsByLine = {};
    cuts.forEach(c => {
      const ln = Math.floor(c.cutTop / bpl);
      (cutsByLine[ln] = cutsByLine[ln] || []).push(c);
    });
    const sel = selPieces(rec, state.selection);
    const out = [];
    for (let s = from; s < to; s += bpl) {
      const n = Math.min(bpl, to - s);
      const top = rec.seq.substr(s, n);
      const ln = Math.floor(s / bpl);
      let body = '';

      // enzyme row
      const lineCuts = cutsByLine[ln] || [];
      if (lineCuts.length) {
        const rows = assignRowsCh(lineCuts.map(c => ({ col: c.cutTop - s, text: c.enzyme + PB.cutMark(c) })), n);
        const rowCount = Math.max.apply(null, rows.map(r => r.row)) + 1;
        for (let r = 0; r < rowCount; r++) {
          body += '<div class="sl-enz">' + rows.filter(x => x.row === r).map(x =>
            '<span class="em" style="left:' + x.col + 'ch">' + esc(x.text) + '</span>').join('') + '</div>';
        }
      }

      // strands
      body += '<div class="sl-strand" data-ls="' + s + '" data-ln="' + n + '">' + strandHtml(top, s, sel, colorOn) + '</div>';
      if (o.bottom) {
        body += '<div class="sl-strand bottom">' + strandHtml(B.complement(top), s, sel, colorOn) + '</div>';
      }

      // translations, from each coding feature's codon map
      if (o.aa) {
        for (const it of ivs) {
          if (!/CDS|ORF/i.test(it.f.type) && !(it.f.type === 'gene' && !rec.features.some(g => g.type === 'CDS'))) continue;
          if (!it.pieces.some(p => p.end > s && p.start < s + n)) continue;
          const codons = PB.codonMap(rec, it.idx);
          const row = [];
          for (const c of codons) {
            const mid = c.pos[1];
            if (mid < s || mid >= s + n) continue;
            row.push('<span class="aa' + (c.aa === '*' ? ' stop' : (c.k === 0 ? ' start' : '')) +
              '" style="left:' + (mid - 1 - s) + 'ch">' + esc(c.aa) + '</span>');
          }
          if (row.length) body += '<div class="sl-aa" title="' + esc(it.f.label) + '">' + row.join('') + '</div>';
        }
      }

      // feature bars: one bar per piece, a thin line across introns
      const lineFeats = ivs.filter(it =>
        it.pieces.some(p => p.end > s && p.start < s + n) ||
        (it.pieces.length > 1 && gapsOf(it, L).some(g => g.end > s && g.start < s + n)));
      if (lineFeats.length) {
        const maxTrack = Math.max.apply(null, lineFeats.map(i => i.track));
        for (let t = 0; t <= maxTrack; t++) {
          const row = lineFeats.filter(i => i.track === t);
          if (!row.length) continue;
          let html = '';
          row.forEach(it => {
            const color = PB.featColor(it.f);
            gapsOf(it, L).forEach(g => {
              const a = Math.max(g.start, s), b = Math.min(g.end, s + n);
              if (b > a) html += '<span class="sl-featgap" style="left:' + (a - s) + 'ch;width:' + (b - a) +
                'ch;border-color:' + color + '"></span>';
            });
            let labelled = false;
            it.pieces.forEach(p => {
              const a = Math.max(p.start, s), b = Math.min(p.end, s + n);
              if (b <= a) return;
              const w = b - a;
              const showText = !labelled && w >= it.f.label.length + 1;
              if (showText) labelled = true;
              html += '<span class="sl-featbar" data-fi="' + it.idx + '" title="' + esc(it.f.label + ' — ' + PB.locationText(rec, it.f)) +
                '" style="left:' + (a - s) + 'ch;width:' + w + 'ch;background:' + color + '">' +
                (showText ? esc(it.f.label) : '') + '</span>';
            });
          });
          body += '<div class="sl-feats" style="height:15px">' + html + '</div>';
        }
      }

      out.push('<div class="seqline"><span class="sl-num">' + fmt(s + 1) + '</span><div class="sl-body">' + body + '</div></div>');
    }
    inner.innerHTML = out.join('');
  }

  /* The selection as plain top-strand ranges: two when it runs across
     base 1, and one per piece of every feature when several features are
     selected at once. strandHtml already paints a list of ranges, so the
     multiple case needs nothing further downstream. */
  function selPieces(rec, sel) {
    const multi = PB.featureSelectionPieces(rec);
    if (multi) return multi;
    const n = PB.normSel(rec, sel);
    if (!n) return [];
    const L = rec.seq.length;
    return n.wraps ? [{ start: n.start, end: L }, { start: 0, end: n.end - L }] : [{ start: n.start, end: n.end }];
  }

  function strandHtml(str, lineStart, sel, color) {
    /* the parts of this line that are selected, in line coordinates */
    const spans = (sel || []).map(p => [Math.max(p.start - lineStart, 0), Math.min(p.end - lineStart, str.length)])
      .filter(([a, b]) => b > a).sort((x, y) => x[0] - y[0]);
    if (!color) {
      if (!spans.length) return esc(str);
      let out = '', at = 0;
      spans.forEach(([a, b]) => {
        out += esc(str.slice(at, a)) + '<span style="background:var(--sel)">' + esc(str.slice(a, b)) + '</span>';
        at = b;
      });
      return out + esc(str.slice(at));
    }
    let out = '';
    for (let i = 0; i < str.length; i++) {
      const inSel = spans.some(([a, b]) => i >= a && i < b);
      out += '<span class="b-' + str[i] + '"' + (inSel ? ' style="background:var(--sel)"' : '') + '>' + str[i] + '</span>';
    }
    return out;
  }

  // the stretches between consecutive pieces of a joined feature (introns),
  // in top-strand order; a hop across the origin is split at it
  function gapsOf(it, L) {
    const f = it.f;
    if (!f.parts || f.parts.length < 2) return [];
    L = L || Infinity;
    const ps = f.mixed ? f.parts.slice().sort((a, b) => a.start - b.start) : f.parts;
    const out = [];
    for (let i = 1; i < ps.length; i++) {
      const prev = ps[i - 1], cur = ps[i];
      const prevEnd = prev.end <= prev.start ? prev.end : prev.end;
      if (cur.start >= prevEnd) { if (cur.start > prevEnd) out.push({ start: prevEnd, end: cur.start }); }
      else {
        if (prevEnd < L) out.push({ start: prevEnd, end: L });
        if (cur.start > 0) out.push({ start: 0, end: cur.start });
      }
    }
    return out;
  }
  PB.gapsOf = gapsOf;

  function assignRowsCh(items, lineLen) {
    items = items.slice().sort((a, b) => a.col - b.col);
    const rowEnd = [];
    items.forEach(it => {
      let r = 0;
      while (r < rowEnd.length && rowEnd[r] > it.col) r++;
      it.row = r;
      rowEnd[r] = it.col + it.text.length + 1;
    });
    return items;
  }

  let charW = 0;
  function measureCharW() {
    const inner = $('#seqInner');
    if (!inner) return 7.8;
    const probe = document.createElement('div');
    probe.className = 'seqline';
    probe.style.cssText = 'position:absolute;visibility:hidden;white-space:pre';
    probe.textContent = '0'.repeat(100);
    inner.appendChild(probe);
    const w = probe.getBoundingClientRect().width / 100;
    probe.remove();
    return w > 0.5 ? w : 7.8;
  }

  function wireSeq(host, rec) {
    const o = state.seqOpts;
    const bplSel = $('#sqBpl');
    if (bplSel) bplSel.addEventListener('change', () => {
      o.bpl = parseInt(bplSel.value, 10) || 60;
      PB.lsSet('seqOpts', o); renderSeq();
    });
    host.querySelectorAll('[data-tg]').forEach(el => el.addEventListener('click', () => {
      const k = el.getAttribute('data-tg');
      o[k] = !o[k];
      PB.lsSet('seqOpts', o); renderSeq();
    }));
    const copyBtn = $('#sqCopy');
    if (copyBtn) copyBtn.addEventListener('click', () => {
      const sel = state.selection;
      const text = PB.normSel(rec, sel) ? PB.selectionSeq(rec, sel) : rec.seq;
      PB.copyText(text, sel ? 'Selection' : 'Full sequence');
    });
    const prev = $('#sqPrev'), next = $('#sqNext');
    if (prev) prev.addEventListener('click', () => { state.seqPage--; renderSeq(); });
    if (next) next.addEventListener('click', () => { state.seqPage++; renderSeq(); });

    const goto = $('#sqGoto');
    if (goto) goto.addEventListener('change', () => {
      const p = Math.max(1, Math.min(rec.seq.length, parseInt(goto.value, 10) || 1)) - 1;
      scrollToPos(p);
    });
    const findAll = $('#sqFindAll');
    if (findAll) findAll.addEventListener('click', () => PB.searchDialog(($('#sqFind').value || '').trim()));
    const findBox = $('#sqFind'), findBtn = $('#sqFindBtn');
    const doFind = () => {
      const q = (findBox.value || '').trim();
      if (!q) return;
      const enz = PB.enzByName(q);
      if (enz) {
        const sites = B.findSites(rec.seq, enz, rec.topology === 'circular');
        if (!sites.length) { PB.toast('No ' + enz.name + ' sites in this sequence.', 'warn'); return; }
        if (state.enzymes.indexOf(enz.name) === -1) {
          state.enzymes.push(enz.name); PB.lsSet('enzymes', state.enzymes); PB.invalidateSites();
        }
        const fromPos = state.selection ? state.selection.end : 0;
        const nextSite = sites.find(s => s.pos >= fromPos) || sites[0];
        state.selection = { start: nextSite.pos, end: nextSite.pos + enz.site.length };
        PB.renderHeader(); renderSeq(); PB.renderInspector(); PB.renderEnzymeRail();
        scrollToPos(nextSite.pos);
        PB.toast(enz.name + ': ' + sites.length + ' site' + (sites.length === 1 ? '' : 's') + ', cuts at ' +
          sites.map(s => fmt(s.cutTop)).join(', '));
        return;
      }
      const pat = B.cleanSeq(q);
      if (!pat) { PB.toast('Enter bases (ACGT…) or an enzyme name.', 'warn'); return; }
      const hay = rec.topology === 'circular' ? rec.seq + rec.seq.slice(0, pat.length - 1) : rec.seq;
      const fromPos = state.selection ? state.selection.start + 1 : 0;
      let i = hay.indexOf(pat, fromPos);
      let rc = false;
      if (i === -1) i = hay.indexOf(pat, 0);
      if (i === -1) {
        const rcPat = B.revcomp(pat);
        i = hay.indexOf(rcPat, fromPos);
        if (i === -1) i = hay.indexOf(rcPat, 0);
        rc = i !== -1;
      }
      if (i === -1) { PB.toast('“' + pat + '” not found on either strand.', 'warn'); return; }
      state.selection = { start: i % rec.seq.length, end: (i % rec.seq.length) + pat.length };
      PB.renderHeader(); renderSeq(); PB.renderInspector();
      scrollToPos(i);
      PB.toast('Found at ' + fmt(i + 1) + (rc ? ' (reverse strand)' : ''));
    };
    if (findBtn) findBtn.addEventListener('click', doFind);
    if (findBox) findBox.addEventListener('keydown', e => { if (e.key === 'Enter') { e.preventDefault(); doFind(); } });

    // selection by dragging across the top strand
    const inner = $('#seqInner');
    charW = measureCharW();
    let dragging = false, anchor = 0;
    function idxFromEvent(e) {
      const el = e.target.closest ? e.target.closest('.sl-strand') : null;
      if (!el || el.classList.contains('bottom')) return null;
      const ls = parseInt(el.getAttribute('data-ls'), 10);
      const n = parseInt(el.getAttribute('data-ln'), 10);
      const rect = el.getBoundingClientRect();
      let col = Math.floor((e.clientX - rect.left) / charW);
      col = Math.max(0, Math.min(n - 1, col));
      return ls + col;
    }
    inner.addEventListener('mousedown', e => {
      if (e.target.classList && e.target.classList.contains('sl-featbar')) return;
      const i = idxFromEvent(e);
      if (i == null) return;
      e.preventDefault();
      dragging = true; anchor = i;
      state.selection = { start: i, end: i + 1 };
      state.selectedFeature = null;
      paintSelection(rec);
    });
    inner.addEventListener('mousemove', e => {
      if (!dragging) return;
      const i = idxFromEvent(e);
      if (i == null) return;
      state.selection = { start: Math.min(anchor, i), end: Math.max(anchor, i) + 1 };
      paintSelection(rec);
    });
    window.addEventListener('mouseup', () => {
      if (!dragging) return;
      dragging = false;
      PB.renderHeader(); PB.renderInspector();
    });
    inner.querySelectorAll('.sl-featbar').forEach(bar => {
      bar.addEventListener('click', e => {
        e.stopPropagation();
        const i = parseInt(bar.getAttribute('data-fi'), 10);
        const f = PB.displayRecord(rec).features[i];
        if (!f) return;
        if (f.virtual) { PB.selectPrimerFeature(rec, f); PB.renderHeader(); renderSeq(); PB.renderInspector(); return; }
        state.selectedFeature = i;
        state.selection = { start: f.start, end: f.end };
        PB.renderHeader(); renderSeq(); PB.renderInspector();
      });
    });
  }

  // repaint only the strand rows + status while dragging (cheap)
  function paintSelection(rec) {
    const o = state.seqOpts;
    const colorOn = o.colorBases && rec.seq.length <= 12000;
    const sel = selPieces(rec, state.selection);
    $$('#seqInner .sl-strand').forEach(el => {
      const isBottom = el.classList.contains('bottom');
      const holder = isBottom ? el.previousElementSibling : el;
      const ls = parseInt((holder && holder.getAttribute('data-ls')) || el.getAttribute('data-ls'), 10);
      const n = parseInt((holder && holder.getAttribute('data-ln')) || el.getAttribute('data-ln'), 10);
      if (isNaN(ls)) return;
      const str = rec.seq.substr(ls, n);
      el.innerHTML = strandHtml(isBottom ? B.complement(str) : str, ls, sel, colorOn);
    });
    updateSeqStatus(rec);
  }

  function scrollToPos(p) {
    const bpl = state.seqOpts.bpl;
    const page = Math.floor(p / SEQ_PAGE);
    if (page !== state.seqPage) { state.seqPage = page; renderSeq(); }
    const inner = $('#seqInner'), scroll = $('#seqScroll');
    if (!inner || !scroll) return;
    const lineIdx = Math.floor((p - state.seqPage * SEQ_PAGE) / bpl);
    const line = inner.children[lineIdx];
    if (line) scroll.scrollTop = Math.max(0, line.offsetTop - 60);
  }

  PB.scrollSeqTo = scrollToPos;

  function updateSeqStatus(rec) {
    const el = $('#seqStatus');
    if (!el) return;
    const sel = state.selection;
    const bits = ['' + fmt(rec.seq.length) + ' bp', f1(B.gcContent(rec.seq)) + '% GC'];
    if (PB.normSel(rec, sel)) {
      const sub = PB.selectionSeq(rec, sel);
      bits.push('selected ' + PB.selRangeText(rec, sel) + ' (' + fmt(sub.length) + ' bp)');
      bits.push(f1(B.gcContent(sub)) + '% GC');
      if (PB.tmShown(sub)) {
        const tm = B.tmNN(sub, PB.conditions());
        if (isFinite(tm)) bits.push('Tm ' + f1(tm) + ' °C');
      }
      if (sub.length % 3 === 0 && sub.length <= 3000) {
        bits.push(B.translate(sub, { table: state.codeTable }).slice(0, 40));
      }
    } else {
      bits.push('drag across the sequence to select');
    }
    el.innerHTML = bits.map(b => '<span>' + esc(b) + '</span>').join('');
  }

  /* ============================== FEATURES ============================== */
  function renderFeatures() {
    const host = $('#view-features');
    const rec = PB.active();
    if (!rec) { host.innerHTML = emptyState('No construct loaded.'); return; }
    const ivs = featIntervals(rec).slice().sort((a, b) => a.start - b.start);
    let html = '<div class="panel"><h3>Features <span class="sub">' + rec.features.length + ' annotated</span>' +
      '<span style="display:flex;gap:6px">' +
      '<button class="btn sm" data-act="detect">Detect common features</button>' +
      '<button class="btn sm" data-act="add">Add feature</button>' +
      (rec.features.length ? '<button class="btn sm" data-act="clear">Remove all</button>' : '') +
      '</span></h3>';
    if (!ivs.length) {
      html += emptyState('No features yet. “Detect common features” scans for promoters, tags, markers and origins from the built-in library, or add one from a selection.');
    } else {
      html += '<div class="scrollx"><table class="tbl"><thead><tr>' +
        '<th>Label</th><th>Type</th><th class="num">Start</th><th class="num">End</th>' +
        '<th class="num">Length</th><th>Strand</th><th class="num">GC</th><th>Product</th><th></th>' +
        '</tr></thead><tbody>' +
        ivs.map(it => {
          const f = it.f;
          const sub = PB.featSeq(rec, f);
          const isCds = /CDS|gene|ORF|exon/i.test(f.type);
          let prod = '';
          if (isCds && sub.length >= 3) {
            const p = B.translate(sub, { table: state.codeTable });
            const aaLen = p.replace(/\*+$/, '').length;
            const internalStop = p.slice(0, -1).indexOf('*') !== -1;
            prod = aaLen + ' aa' +
              (sub.length % 3 ? ' <span class="badge warn">frame</span>' : '') +
              (internalStop ? ' <span class="badge bad">stop</span>' : '');
          }
          return '<tr data-fi="' + it.idx + '">' +
            '<td><span style="display:inline-block;width:9px;height:9px;border-radius:2px;background:' +
              PB.featColor(f) + ';margin-right:6px"></span>' + esc(f.label) + '</td>' +
            '<td class="mono" style="font-size:11.5px">' + esc(f.type) + '</td>' +
            '<td class="num">' + fmt(f.start + 1) + '</td>' +
            '<td class="num">' + fmt(f.end) + '</td>' +
            '<td class="num">' + fmt(it.bpLen) + (f.parts ? ' <span class="badge neutral" title="' +
              esc(PB.locationText(rec, f)) + '">' + f.parts.length + ' parts</span>' : '') +
              (f.end <= f.start && !f.between ? ' <span class="badge neutral" title="crosses the origin">wraps</span>' : '') + '</td>' +
            '<td>' + (f.strand === -1 ? '−' : '+') + '</td>' +
            '<td class="num">' + f1(B.gcContent(sub)) + '</td>' +
            '<td>' + prod + '</td>' +
            '<td style="white-space:nowrap">' +
              '<button class="btn sm ghost" data-sel="' + it.idx + '" title="Select">◧</button>' +
              '<button class="btn sm ghost" data-edit="' + it.idx + '" title="Edit">✎</button>' +
              '<button class="btn sm ghost" data-del="' + it.idx + '" title="Delete">×</button>' +
            '</td></tr>';
        }).join('') + '</tbody></table></div>';
    }
    html += '</div>';
    host.innerHTML = html;

    host.querySelectorAll('[data-act="detect"]').forEach(b => b.addEventListener('click', () => PB.detectFeatures()));
    host.querySelectorAll('[data-act="add"]').forEach(b => b.addEventListener('click', () => PB.featureDialog(null)));
    host.querySelectorAll('[data-act="clear"]').forEach(b => b.addEventListener('click', () => {
      PB.modal({
        title: 'Remove all features?',
        body: '<p>All ' + rec.features.length + ' annotations on “' + esc(rec.name) + '” will be deleted. The sequence itself is unchanged.</p>',
        actions: [{ label: 'Cancel' }, { label: 'Remove all', primary: true, onClick: () => {
          rec.features = []; PB.addHistory(rec, 'Features cleared', 'All annotations removed');
          PB.touch(rec); PB.renderAll();
        } }]
      });
    }));
    host.querySelectorAll('[data-sel]').forEach(b => b.addEventListener('click', ev => {
      const i = parseInt(b.getAttribute('data-sel'), 10);
      const additive = ev.ctrlKey || ev.metaKey;
      PB.toggleFeatureSelection(rec, i, additive);
      PB.renderHeader(); PB.renderInspector();
      /* ctrl-clicking down a list to build a selection should not throw
         you into another tab on the first click */
      if (!additive) { state.view = 'seq'; PB.lsSet('view', 'seq'); PB.renderTabs(); }
      PB.renderViews();
    }));
    host.querySelectorAll('[data-edit]').forEach(b => b.addEventListener('click', () =>
      PB.featureDialog(parseInt(b.getAttribute('data-edit'), 10))));
    host.querySelectorAll('[data-del]').forEach(b => b.addEventListener('click', () => {
      const i = parseInt(b.getAttribute('data-del'), 10);
      const f = rec.features[i];
      rec.features.splice(i, 1);
      state.selectedFeature = null;
      PB.addHistory(rec, 'Feature removed', f.label);
      PB.touch(rec); PB.renderAll();
    }));
  }
  PB.renderFeatures = renderFeatures;

  /* ============================== ENZYMES ============================== */
  function renderEnzymes() {
    const host = $('#view-enzymes');
    const rec = PB.active();
    if (!rec || !rec.seq.length) { host.innerHTML = emptyState('No sequence loaded.'); return; }
    const circ = rec.topology === 'circular';
    const list = PB.allEnzymes();
    const rows = [];
    const meth = PB.methylSets(rec), methSys = PB.methylSystems(rec);
    const methOn = Object.keys(methSys).length > 0;
    for (const e of list) {
      let sites;
      try { sites = B.findSites(rec.seq, e, circ); } catch (err) { continue; }
      const eff = methOn ? sites.map(st => B.methylationEffect(e, st, rec.seq.length, meth)) : [];
      const blocked = eff.filter(m => m && m.effect === 'blocked');
      const doubt = eff.filter(m => m && m.effect !== 'blocked');
      rows.push({ e, sites, blocked, doubt });
    }
    const q = (state.enzViewQuery || '').toLowerCase();
    const mode = state.enzViewMode || 'cutters';
    let shown = rows.filter(r => {
      if (q && r.e.name.toLowerCase().indexOf(q) === -1 && r.e.site.toLowerCase().indexOf(q) === -1) return false;
      if (mode === 'cutters') return r.sites.length > 0;
      if (mode === 'unique') return r.sites.length === 1;
      if (mode === 'uniqueCut') return r.sites.length === 1 && !r.blocked.length && !r.doubt.length;
      if (mode === 'meth') return r.blocked.length + r.doubt.length > 0;
      if (mode === 'double') return r.sites.length === 2;
      if (mode === 'none') return r.sites.length === 0;
      return true;
    });
    shown.sort((a, b) => a.sites.length - b.sites.length || a.e.name.localeCompare(b.e.name));

    const uniqueCount = rows.filter(r => r.sites.length === 1).length;
    const nonCount = rows.filter(r => r.sites.length === 0).length;

    const sysName = k => ({ dam: 'Dam', dcm: 'Dcm', cpg: 'CpG' })[k] || k;
    function methBadge(r) {
      if (!r.sites.length || (!r.blocked.length && !r.doubt.length)) return '';
      const by = arr => Array.from(new Set([].concat.apply([], arr.map(m => m.by)))).map(sysName).join('/');
      let h = '';
      if (r.blocked.length) h += '<span class="badge bad" title="The NEB chart lists this enzyme as blocked when its site overlaps this methylation">' +
        (r.blocked.length === r.sites.length ? 'all ' : r.blocked.length + ' of ' + r.sites.length + ' ') + 'blocked · ' + by(r.blocked) + '</span> ';
      const maybe = r.doubt.filter(m => m.effect !== 'unknown'), unk = r.doubt.filter(m => m.effect === 'unknown');
      if (maybe.length) h += '<span class="badge warn" title="Impaired, or blocked only by some overlapping combinations — expect partial cutting">' +
        maybe.length + ' impaired/possible · ' + by(maybe) + '</span> ';
      if (unk.length) h += '<span class="badge neutral" title="Not in the NEB chart; check the supplier">' + unk.length + ' overlap ' + by(unk) + ', sensitivity unknown</span>';
      return h;
    }
    host.innerHTML = '<div class="panel"><h3>Restriction sites ' +
      '<span class="sub">' + rows.length + ' enzymes screened · ' + uniqueCount + ' unique cutters · ' + nonCount + ' non-cutters · methylation: ' +
        esc(PB.METHYL_MODES[PB.methylMode(rec)].label.split(' (')[0]) + '</span></h3>' +
      '<div class="row">' +
      '<div class="field" style="flex:0 1 180px"><label>Filter</label><input class="inp mono" id="evQ" value="' + esc(state.enzViewQuery || '') + '" placeholder="name or site"></div>' +
      '<div class="field narrow" style="flex:0 0 230px"><label>Show</label><select class="inp" id="evMode">' +
        [['cutters', 'Cutters'], ['unique', 'Unique (1 site)'], ['uniqueCut', 'Unique, cuts despite methylation'], ['double', 'Two sites'], ['none', 'Non-cutters'], ['meth', 'Affected by methylation'], ['all', 'All enzymes']]
          .map(([v, l]) => '<option value="' + v + '"' + (mode === v ? ' selected' : '') + '>' + l + '</option>').join('') +
      '</select></div>' +
      '<div class="field" style="flex:1 1 auto;align-self:flex-end">' +
      '<button class="btn sm" id="evSelAll">Show all listed on map</button></div>' +
      '</div>' +
      '<div class="scrollx scrolly"><table class="tbl"><thead><tr>' +
      '<th style="width:34px"></th><th>Enzyme</th><th>Recognition</th><th class="num">Sites</th>' +
      '<th>Cut positions</th><th>Ends</th><th>Fragments</th><th>Methylation</th></tr></thead><tbody>' +
      shown.map(r => {
        const on = state.enzymes.indexOf(r.e.name) !== -1;
        const positions = r.sites.slice(0, 9).map(s => fmt(s.cutTop)).join(', ') +
          (r.sites.length > 9 ? ' +' + (r.sites.length - 9) + ' more' : '');
        let ends = '—';
        if (r.sites.length) {
          const s0 = r.sites[0];
          const d = s0.cutBottom - s0.cutTop;
          if (d === 0) ends = 'blunt';
          else {
            const L = rec.seq.length;
            const oh = d > 0
              ? rec.seq.substr(((s0.cutTop % L) + L) % L, d)
              : rec.seq.substr(((s0.cutBottom % L) + L) % L, -d);
            ends = (d > 0 ? "5′ " : "3′ ") + oh;
          }
        }
        let frags = '—';
        if (r.sites.length) {
          const sizes = B.digestSizes(rec.seq, [r.e.name], { circular: circ, methylation: methSys });
          frags = sizes.slice(0, 6).map(fmt).join(' / ') + (sizes.length > 6 ? ' …' : '');
        }
        return '<tr>' +
          '<td><input type="checkbox" data-ev="' + esc(r.e.name) + '"' + (on ? ' checked' : '') + ' style="accent-color:var(--accent)"></td>' +
          '<td><b>' + esc(r.e.name) + '</b>' + (r.e.typeIIS ? ' <span class="badge neutral">IIS</span>' : '') +
            (r.e.custom ? ' <span class="badge acc">custom</span>' : '') + '</td>' +
          '<td class="mono">' + esc(r.e.site) + '</td>' +
          '<td class="num"><b style="color:' + (r.sites.length === 1 ? 'var(--good)' : 'inherit') + '">' + r.sites.length + '</b></td>' +
          '<td class="mono" style="font-size:11.5px">' + esc(positions) + '</td>' +
          '<td class="mono" style="font-size:11.5px">' + esc(ends) + '</td>' +
          '<td class="mono" style="font-size:11.5px">' + esc(frags) + '</td>' +
          '<td>' + methBadge(r) + '</td>' +
        '</tr>';
      }).join('') + '</tbody></table></div>' +
      (shown.length ? '' : '<div class="hint">Nothing matches this filter.</div>') +
      '</div>';

    const qEl = $('#evQ');
    qEl.addEventListener('input', () => { state.enzViewQuery = qEl.value; renderEnzymes(); });
    $('#evMode').addEventListener('change', e => { state.enzViewMode = e.target.value; renderEnzymes(); });
    $('#evSelAll').addEventListener('click', () => {
      shown.forEach(r => { if (state.enzymes.indexOf(r.e.name) === -1) state.enzymes.push(r.e.name); });
      PB.lsSet('enzymes', state.enzymes); PB.invalidateSites(); PB.renderAll();
    });
    host.querySelectorAll('[data-ev]').forEach(cb => cb.addEventListener('change', () => {
      const n = cb.getAttribute('data-ev');
      if (cb.checked) { if (state.enzymes.indexOf(n) === -1) state.enzymes.push(n); }
      else state.enzymes = state.enzymes.filter(x => x !== n);
      PB.lsSet('enzymes', state.enzymes); PB.invalidateSites();
      PB.renderHeader(); PB.renderEnzymeRail();
    }));
  }
  PB.renderEnzymes = renderEnzymes;

  /* ============================== ORFs ============================== */
  function renderORFs() {
    const host = $('#view-orfs');
    const rec = PB.active();
    if (!rec || !rec.seq.length) { host.innerHTML = emptyState('No sequence loaded.'); return; }
    const minAA = state.orfMin == null ? 50 : state.orfMin;
    const requireStart = state.orfStart !== false;
    let orfs = [];
    try {
      orfs = B.findORFs(rec.seq, { table: state.codeTable, minAA, requireStart, circular: rec.topology === 'circular' });
    } catch (e) { orfs = []; }
    orfs.sort((a, b) => b.length - a.length);
    const shown = orfs.slice(0, 300);

    host.innerHTML = '<div class="panel"><h3>Open reading frames ' +
      '<span class="sub">' + orfs.length + ' found in 6 frames, genetic code ' + state.codeTable + '</span></h3>' +
      '<div class="row">' +
      '<div class="field narrow"><label>Min length (aa)</label><input class="inp mono" id="orfMin" type="number" min="1" value="' + minAA + '"></div>' +
      '<div class="field" style="flex:0 1 190px;align-self:flex-end">' +
      '<label class="chip' + (requireStart ? ' on' : '') + '" id="orfStart">Require ATG start</label></div>' +
      '</div>' +
      (shown.length ? '<div class="scrollx scrolly"><table class="tbl"><thead><tr>' +
        '<th class="num">Start</th><th class="num">End</th><th>Frame</th><th class="num">aa</th>' +
        '<th class="num">bp</th><th>Protein</th><th></th></tr></thead><tbody>' +
        shown.map((o, i) => '<tr>' +
          '<td class="num">' + fmt(o.start + 1) + '</td>' +
          '<td class="num">' + fmt(o.end) + '</td>' +
          '<td class="mono">' + (o.strand === 1 ? '+' : '−') + Math.abs(o.frame) +
            (o.wrapped ? ' <span class="badge neutral" title="runs across the origin">wraps</span>' : '') + '</td>' +
          '<td class="num"><b>' + fmt(o.length) + '</b></td>' +
          '<td class="num">' + fmt(o.ntLength) + '</td>' +
          '<td class="seqcell" style="max-width:300px">' + esc(o.protein.slice(0, 70)) + (o.protein.length > 70 ? '…' : '') + '</td>' +
          '<td style="white-space:nowrap">' +
            '<button class="btn sm ghost" data-osel="' + i + '" title="Select">◧</button>' +
            '<button class="btn sm ghost" data-oadd="' + i + '" title="Annotate as CDS">+</button>' +
            '<button class="btn sm ghost" data-ocopy="' + i + '" title="Copy protein">⧉</button>' +
          '</td></tr>').join('') +
        '</tbody></table></div>' +
        (orfs.length > 300 ? '<div class="hint">Showing the 300 longest of ' + fmt(orfs.length) + '. Raise the minimum length to narrow the list.</div>' : '')
        : emptyState('No ORFs of at least ' + minAA + ' aa. Lower the minimum length, or untick “Require ATG start”.')) +
      '</div>';

    $('#orfMin').addEventListener('change', e => {
      state.orfMin = Math.max(1, parseInt(e.target.value, 10) || 50); renderORFs();
    });
    $('#orfStart').addEventListener('click', () => { state.orfStart = !requireStart; renderORFs(); });
    host.querySelectorAll('[data-osel]').forEach(b => b.addEventListener('click', () => {
      const o = shown[parseInt(b.getAttribute('data-osel'), 10)];
      state.selection = { start: o.start, end: o.end };
      state.selectedFeature = null;
      PB.renderHeader(); PB.renderInspector();
      state.view = 'seq'; PB.lsSet('view', 'seq'); PB.renderTabs(); PB.renderViews();
    }));
    host.querySelectorAll('[data-oadd]').forEach(b => b.addEventListener('click', () => {
      const o = shown[parseInt(b.getAttribute('data-oadd'), 10)];
      rec.features.push(PB.normFeature({
        type: 'CDS', start: o.start, end: o.end, strand: o.strand,
        label: 'ORF ' + fmt(o.start + 1) + '–' + fmt(o.end) + ' (' + o.length + ' aa)'
      }));
      PB.addHistory(rec, 'ORF annotated', o.length + ' aa at ' + fmt(o.start + 1) + '–' + fmt(o.end));
      PB.touch(rec); PB.renderAll();
      PB.toast('Annotated as CDS', 'ok');
    }));
    host.querySelectorAll('[data-ocopy]').forEach(b => b.addEventListener('click', () => {
      const o = shown[parseInt(b.getAttribute('data-ocopy'), 10)];
      PB.copyText(o.protein, 'Protein sequence');
    }));
  }
  PB.renderORFs = renderORFs;
})(window.PB);
