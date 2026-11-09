
/* ============================================================================
   PlasmidBench — map viewer with semantic zoom.

   Linear view: zooming changes how many bases span the viewport, not the size
   of the glyphs. Features, labels and ruler text keep a constant screen size,
   only the bases visible change, and past ~7 px per base the sequence itself
   is drawn. Circular view: the ring scales geometrically while every label is
   emitted at font-size/k so text stays the same size on screen.
   ============================================================================ */
(function (PB) {
  'use strict';
  const B = window.Bio;
  const $ = PB.$, esc = PB.esc, fmt = PB.fmt, f1 = PB.f1;
  const state = PB.state;

  const CIRC_MIN = 0.35, CIRC_MAX = 120;
  const MIN_SPAN = 24;              // closest the linear view may zoom, in bp
  const SEQ_PX = 4.5;               // px per base at which letters appear
  const COMP_PX = 6;                // px per base at which the complement appears
  const AA_PX = 7;                  // px per base at which translation appears

  function niceStep(span) {
    const raw = span / 10;
    const mag = Math.pow(10, Math.floor(Math.log10(Math.max(raw, 1))));
    const n = raw / mag;
    const m = n <= 1 ? 1 : n <= 2 ? 2 : n <= 5 ? 5 : 10;
    return Math.max(1, m * mag);
  }

  /* Labels round a circular map pile up wherever the features do, and on a
     plasmid everything interesting sits in the cloning site: twenty labels
     out of forty degrees of arc, stacked into one dense run that no leader
     line can be followed through. So a crowded side fans into more than one
     column. Alternate labels take the far column, which doubles the room
     each one has without moving it away from its own angle. */
  /* Labels round a circular map pile up wherever the features do, and on a
     plasmid everything interesting sits in the cloning site: twenty-five
     labels out of forty degrees of arc, stacked into one tall run with
     leaders crossing each other all the way down.

     Two things fix that. A side holding more labels than reads comfortably
     fans into columns, each taking a contiguous block of the arc, so its
     leaders stay short and roughly parallel instead of raking across the
     whole map. And where there is room to spare the lines are spread to
     use it, because a label with air around it is read at a glance and one
     in a dense run is not. */
  const LBL_PER_COL = 13;
  function layoutLabels(items, cx, cy, R, lineH, hLimit, room, charW) {
    room = room || 0; charW = charW || 5.9;
    items.forEach(it => {
      const a = it.ang;
      it.lx = cx + Math.cos(a) * R;
      it.ly = cy + Math.sin(a) * R;
      it.side = Math.cos(a) >= 0 ? 1 : -1;
      it.col = 0;
    });
    const top = lineH, bottom = hLimit - lineH;
    const avail = Math.max(lineH * 4, bottom - top);
    [1, -1].forEach(side => {
      const all = items.filter(i => i.side === side).sort((a, b) => a.ly - b.ly);
      if (!all.length) return;
      /* a column has to be as wide as its widest label, or the text of one
         column runs into the next; grouped enzyme names are long */
      const colW = all.reduce((w, it) => Math.max(w, (it.text || '').length * charW + 14), 40);
      const fit = Math.max(1, Math.floor(room / colW));
      const cols = Math.max(1, Math.min(fit, 3,
        Math.max(Math.ceil(all.length / LBL_PER_COL), Math.ceil(all.length * lineH / avail))));
      /* contiguous blocks, not alternating: each column then serves one
         stretch of the arc and its leaders do not rake across the others */
      const per = Math.ceil(all.length / cols);
      all.forEach((it, i) => { it.col = Math.min(cols - 1, Math.floor(i / per)); });
      for (let c = 0; c < cols; c++) {
        const grp = all.filter(i => i.col === c);
        if (!grp.length) continue;
        /* Only as far apart as they need to be. Spreading them to fill the
           canvas looks orderly on its own but drags every label away from
           the feature it names, and the leader then has to rake back across
           the map — which is what makes a turned plasmid look a mess. */
        const step = lineH * 1.08;
        for (let i = 1; i < grp.length; i++) {
          if (grp[i].ly - grp[i - 1].ly < step) grp[i].ly = grp[i - 1].ly + step;
        }
        const over = grp[grp.length - 1].ly - bottom;
        if (over > 0) {
          for (let i = grp.length - 1; i >= 0; i--) {
            grp[i].ly -= over;
            if (i > 0 && grp[i].ly - grp[i - 1].ly >= step) break;
          }
        }
        for (let i = 0; i < grp.length; i++) if (grp[i].ly < top) grp[i].ly = top + i * step;
      }
      all.forEach(it => {
        const dy = it.ly - cy;
        const inner = R * R - dy * dy;
        const dx = inner > 0 ? Math.sqrt(inner) : 0;
        it.lx = cx + side * (Math.max(dx, R * 0.46) + it.col * colW);
      });
    });
    return items;
  }

  /* ---------------- per-construct view state ---------------- */
  function views() {
    if (!state.mapViews) state.mapViews = {};
    return state.mapViews;
  }
  PB.mapView = function (rec, mode) {
    rec = rec || PB.active();
    mode = mode || state.mapMode;
    if (!rec) return null;
    const all = views();
    const key = rec.id + '|' + mode;
    if (!all[key]) {
      all[key] = mode === 'linear'
        ? { mode: 'linear', start: 0, span: Math.max(rec.seq.length, 1) }
        : { mode: 'circular', k: 1, x: 0, y: 0 };
    }
    return all[key];
  };
  PB.resetMapZoom = function () {
    const rec = PB.active();
    if (!rec) return;
    const v = PB.mapView(rec);
    if (!v) return;
    if (v.mode === 'linear') { v.start = 0; v.span = Math.max(rec.seq.length, 1); }
    else { v.k = 1; v.x = 0; v.y = 0; }
  };

  function clampLinear(rec, v) {
    const L = Math.max(rec.seq.length, 1);
    v.span = Math.max(Math.min(MIN_SPAN, L), Math.min(v.span, L));
    v.start = Math.max(0, Math.min(v.start, L - v.span));
  }

  function zoomPercent(rec, v) {
    if (!v) return 100;
    if (v.mode === 'circular') return Math.round(v.k * 100);
    const L = Math.max(rec.seq.length, 1);
    return Math.round(L / Math.max(v.span, 1) * 100);
  }

  /* ---------------- toolbar ---------------- */
  function mapToolbar(rec) {
    const v = PB.mapView(rec);
    const pct = zoomPercent(rec, v);
    const range = v && v.mode === 'linear'
      ? fmt(Math.round(v.start) + 1) + '–' + fmt(Math.min(rec.seq.length, Math.round(v.start + v.span))) +
        ' of ' + fmt(rec.seq.length) + ' bp'
      : fmt(rec.seq.length) + ' bp';
    return '<div class="map-toolbar">' +
      '<button class="btn sm' + (state.mapMode === 'circular' ? ' primary' : '') + '" data-mm="circular">Circular</button>' +
      '<button class="btn sm' + (state.mapMode === 'linear' ? ' primary' : '') + '" data-mm="linear">Linear</button>' +
      '<span style="width:8px"></span>' +
      '<div class="zoomctl">' +
        '<button class="btn sm icon" data-zm="out" title="Zoom out  (−)" aria-label="Zoom out">−</button>' +
        '<button class="btn sm" id="mapZoomPct" data-zm="fit" title="Zoom to fit  (0)">' + pct + '%</button>' +
        '<button class="btn sm icon" data-zm="in" title="Zoom in  (+)" aria-label="Zoom in">+</button>' +
      '</div>' +
      (state.mapMode === 'linear' && state.selection
        ? '<button class="btn sm" data-zm="sel" title="Zoom to the selection">Fit selection</button>' : '') +
      '<span style="width:8px"></span>' +
      '<button class="btn sm" data-act="detect">Detect features</button>' +
      '<button class="btn sm" data-act="setorigin"' + (rec.topology === 'circular' ? '' : ' disabled title="Circular sequences only"') + '>Set origin</button>' +
      '<button class="btn sm" data-act="flip">Reverse complement</button>' +
      '<button class="btn sm" data-act="mapexport" title="Save the map as an SVG or PNG image">Export image</button>' +
      (rec.topology === 'circular'
        ? '<button class="btn sm" data-act="mapreset" title="Back to the whole plasmid, upright">Reset view</button>'
        : '') +
      (state.primers.length ? (function () {
        const n = state.showPrimers === false ? 0 : PB.primerFeatures(rec).length;
        const sel = state.primersShown;
        const some = !!sel && sel.length !== state.primers.length;
        return '<button class="btn sm' + (n ? ' primary' : '') + '" data-act="mapprimers" title="' +
          (n ? n + ' site' + (n === 1 ? '' : 's') + ' drawn \u2014 click to choose which primers to show'
             : 'Click to choose which primers to show') + '">' +
          'Primers' + (n ? ' (' + n + ')' : '') + (some ? ' \u2022' : '') + '</button>';
      })() : '') +
      '<span class="spacer"></span>' +
      (rec.topology === 'circular'
        ? '<span class="maphint" style="opacity:.7">scroll to turn it \u00b7 pinch to zoom in to the bases</span>'
        : '') +
      '<span class="mono maphint">' + esc(range) + '</span>' +
      '</div>';
  }

  /* Which primers belong on the map. A library grows as a project does, and
     by the tenth primer the map is a thicket of sites from designs that are
     finished. The ones worth seeing are chosen here, and the choice is
     remembered. The same panel hides the primer-binding sites annotated in
     the file itself, which on a commercial backbone are a dozen labels for
     oligos nobody is using. */
  PB.primerShowDialog = function () {
    const rec = PB.active();
    if (!rec) return;
    const hits = PB.primerHits(rec);
    const count = nm => hits.filter(h => h.name === nm).length;
    const sel = state.primersShown;
    const on = nm => !sel || sel.indexOf(nm) >= 0;
    const fileSites = (rec.features || []).filter(f => /^primer_bind$/i.test(f.type || '')).length;
    const w = PB.modal({
      title: 'Primers on the map',
      sub: state.primers.length + ' in the library',
      body:
        '<div class="row" style="gap:6px;margin-bottom:8px">' +
          '<button class="btn sm" id="pmAll">All</button>' +
          '<button class="btn sm" id="pmNone">None</button>' +
          '<button class="btn sm" id="pmHere">Only those that anneal here</button>' +
        '</div>' +
        '<div class="scrollx" style="max-height:320px;overflow-y:auto"><table class="tbl"><thead><tr>' +
        '<th></th><th>Name</th><th>Sequence</th><th class="num">Sites here</th></tr></thead><tbody>' +
        state.primers.map((p, i) =>
          '<tr><td><label class="chk"><input type="checkbox" data-pm="' + i + '"' + (on(p.name) ? ' checked' : '') +
            '></label></td>' +
          '<td>' + esc(p.name) + '</td>' +
          '<td class="mono" style="font-size:11px;max-width:260px;overflow:hidden;text-overflow:ellipsis">' +
            esc(p.seq.length > 34 ? p.seq.slice(0, 34) + '\u2026' : p.seq) + '</td>' +
          '<td class="num">' + (count(p.name) || '\u2014') + '</td></tr>').join('') +
        '</tbody></table></div>' +
        (fileSites ? '<div class="row" style="margin-top:10px"><div class="field wide">' +
          '<label class="chk"><input type="checkbox" id="pmFile"' +
          (state.showFilePrimers === false ? '' : ' checked') + '> Also show the ' + fileSites +
          ' primer-binding site' + (fileSites === 1 ? '' : 's') + ' annotated in the file</label></div></div>' : '') +
        '<div class="hint">Unticked primers stay in the library and keep working everywhere else \u2014 ' +
        'they are simply not drawn on the map or the sequence.</div>',
      actions: [{ label: 'Done', primary: true }]
    });
    const apply = () => {
      const boxes = Array.from(w.querySelectorAll('[data-pm]'));
      const names = boxes.filter(b => b.checked).map(b => state.primers[+b.getAttribute('data-pm')].name);
      state.primersShown = names.length === state.primers.length ? null : names;
      state.showPrimers = names.length > 0;
      PB.lsSet('primersShown', state.primersShown);
      PB.lsSet('showPrimers', state.showPrimers);
      const f = w.querySelector('#pmFile');
      if (f) { state.showFilePrimers = f.checked; PB.lsSet('showFilePrimers', state.showFilePrimers); }
      PB.renderViews();
    };
    w.querySelectorAll('[data-pm]').forEach(b => b.addEventListener('change', apply));
    const fb = w.querySelector('#pmFile');
    if (fb) fb.addEventListener('change', apply);
    const setAll = v => { w.querySelectorAll('[data-pm]').forEach(b => { b.checked = v; }); apply(); };
    w.querySelector('#pmAll').addEventListener('click', () => setAll(true));
    w.querySelector('#pmNone').addEventListener('click', () => setAll(false));
    w.querySelector('#pmHere').addEventListener('click', () => {
      w.querySelectorAll('[data-pm]').forEach(b => {
        b.checked = count(state.primers[+b.getAttribute('data-pm')].name) > 0;
      });
      apply();
    });
  };

  /* ---------------- render ---------------- */
  let rafPending = false;
  function scheduleMapRender() {
    if (rafPending) return;
    rafPending = true;
    requestAnimationFrame(() => { rafPending = false; drawMap(); });
  }
  PB.scheduleMapRender = scheduleMapRender;

  function canvasWidth() {
    const vp = $('#viewport');
    const w = vp ? vp.clientWidth : 900;
    return Math.max(320, Math.min(1600, w - 24));
  }

  function renderMap() {
    const host = $('#view-map');
    const rec = PB.active();
    if (!rec) { host.innerHTML = PB.emptyState('No construct loaded.'); return; }
    // a construct opens in the projection that suits its topology
    if (state.lastMapRecord !== rec.id) {
      state.lastMapRecord = rec.id;
      state.mapMode = rec.topology === 'circular' ? 'circular' : 'linear';
    }
    if (!rec.seq.length) {
      host.innerHTML = mapToolbar(rec) +
        PB.emptyState('This construct has no sequence yet. Paste one in the Sequence tab.');
      wireToolbar(host, rec);
      return;
    }
    host.innerHTML = mapToolbar(rec) + '<div class="map-wrap" id="mapWrap"></div>';
    wireToolbar(host, rec);
    drawMap();
  }
  PB.renderMap = renderMap;

  // redraw only the drawing, keeping the toolbar; cheap enough for live zoom
  function drawMap() {
    /* drawn with library primers added, when they are shown */
    const rec = PB.displayRecord(PB.active());
    const wrap = $('#mapWrap');
    if (!rec || !wrap || !rec.seq.length) return;
    const v = PB.mapView(rec);
    const W = canvasWidth();
    wrap.innerHTML = v.mode === 'linear' ? linearMap(rec, v, W) : circularMap(rec, v, W);
    const pct = $('#mapZoomPct');
    if (pct) pct.textContent = zoomPercent(rec, v) + '%';
    const hint = document.querySelector('#view-map .maphint');
    if (hint) {
      hint.textContent = v.mode === 'linear'
        ? fmt(Math.round(v.start) + 1) + '–' + fmt(Math.min(rec.seq.length, Math.round(v.start + v.span))) +
          ' of ' + fmt(rec.seq.length) + ' bp'
        : fmt(rec.seq.length) + ' bp';
    }
    wireCanvas(wrap, rec, v, W);
  }

  function wireToolbar(host, rec) {
    host.querySelectorAll('[data-mm]').forEach(b => b.addEventListener('click', () => {
      state.mapMode = b.getAttribute('data-mm');
      PB.lsSet('mapMode', state.mapMode);
      renderMap();
    }));
    host.querySelectorAll('[data-act="detect"]').forEach(b => b.addEventListener('click', () => PB.detectFeatures()));
    host.querySelectorAll('[data-act="setorigin"]').forEach(b => b.addEventListener('click', () => PB.setOriginDialog()));
    host.querySelectorAll('[data-act="flip"]').forEach(b => b.addEventListener('click', () => PB.flipDialog()));
    host.querySelectorAll('[data-act="mapexport"]').forEach(b => b.addEventListener('click', () => PB.exportMapDialog()));
    host.querySelectorAll('[data-act="mapreset"]').forEach(b => b.addEventListener('click', () => {
      const v = PB.mapView(rec);
      v.k = 1; v.x = 0; v.y = 0; v.rot = 0;
      drawMap();
    }));
    host.querySelectorAll('[data-act="mapprimers"]').forEach(b => b.addEventListener('click', () => PB.primerShowDialog()));
    host.querySelectorAll('[data-zm]').forEach(b => b.addEventListener('click', () => {
      const v = PB.mapView(rec);
      const m = b.getAttribute('data-zm');
      if (m === 'fit') PB.resetMapZoom();
      else if (m === 'sel' && state.selection) {
        const pad = Math.max(12, (state.selection.end - state.selection.start) * 0.15);
        v.start = state.selection.start - pad;
        v.span = (state.selection.end - state.selection.start) + pad * 2;
        clampLinear(rec, v);
      } else if (v.mode === 'linear') {
        const mid = v.start + v.span / 2;
        v.span = v.span * (m === 'in' ? 1 / 1.6 : 1.6);
        clampLinear(rec, v);
        v.start = mid - v.span / 2;
        clampLinear(rec, v);
      } else {
        const W = canvasWidth(), H = 640;
        zoomCircular(v, v.k * (m === 'in' ? 1.5 : 1 / 1.5), { x: W / 2, y: H / 2 }, circGeom(W));
      }
      drawMap();
    }));
  }

  /* Zooming a circle at its centre magnifies empty space: the sequence is
     all out at the rim. So an anchor that falls inside the ring is moved
     out to the ring first, in the direction it already points, and zooming
     in then does what it is asked to do — bring the sequence closer. */
  /* Zoom holds the point under the cursor still, which is the one thing a
     zoom must do: whatever is being looked at stays where it is and grows.
     Nothing is snapped to the ring on the way — that made every gesture
     drift to the same place however carefully it was aimed. The only
     correction is a rescue when the ring has gone off screen entirely,
     where there would otherwise be nothing to see at all. */
  function zoomCircular(v, k2, c, geom) {
    k2 = Math.max(CIRC_MIN, Math.min(CIRC_MAX, k2));
    const m = { x: (c.x - v.x) / v.k, y: (c.y - v.y) / v.k };
    if (geom) {
      /* What is held still is the ring, at the angle the cursor points to,
         pinned exactly where it already sits on screen - not the cursor
         position itself. The sequence lives on a thin curve, and holding a
         point of empty space still sends that curve off the canvas within
         a step or two, which is what the rescue then had to haul back:
         the jump you feel when zooming. Pointing at the ring behaves as
         before, since there the two anchors are the same point; pointing
         anywhere else now simply magnifies the arc nearest the cursor
         instead of throwing it away. */
      const dx = m.x - geom.cx, dy = m.y - geom.cy;
      const d = Math.sqrt(dx * dx + dy * dy);
      const a = d > 1e-6 ? Math.atan2(dy, dx) : (-Math.PI / 2 + (v.rot || 0));
      const ax = geom.cx + Math.cos(a) * geom.R, ay = geom.cy + Math.sin(a) * geom.R;
      const sx = ax * v.k + v.x, sy = ay * v.k + v.y;    // where it is now
      v.k = k2;
      v.x = sx - ax * k2;
      v.y = sy - ay * k2;
      rescueCircular(v, geom);
      return;
    }
    v.k = k2;
    v.x = c.x - m.x * k2;
    v.y = c.y - m.y * k2;
  }
  /* Whatever the zoom, some of the ring has to be on screen: the sequence
     lives on it, and a view of the empty middle is a view of nothing. The
     part nearest where the eye already is gets brought back. */
  /* Only when none of the ring is on screen: bring the nearest piece of it
     back to the middle. While any of it is in view the view is left exactly
     where the person put it. */
  /* Does any part of the ring cross the canvas? Asked exactly, rather than
     by testing the single point of the ring nearest the middle of the
     screen. That point is a poor stand-in once the circle is large: the
     visible arc can run right across the canvas while the nearest point
     sits just past an edge, and the rescue would then haul the view back
     on a zoom that needed no rescuing - which is the jump you feel when
     zooming in. The circle crosses the rectangle when its radius is no
     smaller than the distance to the nearest point of the rectangle and
     no larger than the distance to its furthest corner. */
  function ringCrossesCanvas(v, geom) {
    const k = v.k, W = geom.W, H = geom.H;
    const C = { x: geom.cx * k + v.x, y: geom.cy * k + v.y };
    const Rk = geom.R * k;
    const dx = Math.max(0 - C.x, 0, C.x - W);
    const dy = Math.max(0 - C.y, 0, C.y - H);
    const near = Math.sqrt(dx * dx + dy * dy);
    const far = Math.max(
      Math.sqrt(C.x * C.x + C.y * C.y),
      Math.sqrt((C.x - W) * (C.x - W) + C.y * C.y),
      Math.sqrt(C.x * C.x + (C.y - H) * (C.y - H)),
      Math.sqrt((C.x - W) * (C.x - W) + (C.y - H) * (C.y - H)));
    return Rk >= near && Rk <= far;
  }
  function rescueCircular(v, geom) {
    if (ringCrossesCanvas(v, geom)) return;      /* something to look at */
    const k = v.k, W = geom.W, H = geom.H;
    const Ck = { x: geom.cx * k + v.x, y: geom.cy * k + v.y };
    const P = { x: W / 2, y: H / 2 };
    const d = Math.sqrt((P.x - Ck.x) * (P.x - Ck.x) + (P.y - Ck.y) * (P.y - Ck.y));
    const a = d > 0.5 ? Math.atan2(P.y - Ck.y, P.x - Ck.x) : (-Math.PI / 2 + (v.rot || 0));
    const rx = geom.cx + Math.cos(a) * geom.R, ry = geom.cy + Math.sin(a) * geom.R;
    const sx = rx * k + v.x, sy = ry * k + v.y;
    v.x += (W / 2 - sx);
    v.y += (H * 0.46 - sy);
  }
  function circGeom(W) {
    const H = 640;
    return { cx: W / 2, cy: H / 2, W: W, H: H,
             R: Math.max(90, Math.min(200, Math.min(W, H) / 2 - 70)) };
  }

  /* ---------------- interaction ---------------- */
  let drag = null;
  let installed = false;
  function installGlobal() {
    if (installed) return;
    installed = true;
    window.addEventListener('mousemove', e => {
      if (!drag) return;
      const rec = PB.active();
      if (!rec) return;
      const v = PB.mapView(rec);
      const dx = e.clientX - drag.lx, dy = e.clientY - drag.ly;
      drag.moved += Math.abs(dx) + Math.abs(dy);
      drag.lx = e.clientX; drag.ly = e.clientY;
      if (drag.kind === 'linear') {
        v.start -= dx * (v.span / drag.plotW);
        clampLinear(rec, v);
      } else if (drag.kind === 'mini') {
        const L = rec.seq.length;
        v.start = (e.clientX - drag.rectLeft - drag.padL) / drag.plotW * L - v.span / 2;
        clampLinear(rec, v);
      } else {
        v.x += dx * drag.unit; v.y += dy * drag.unit;
      }
      scheduleMapRender();
    });
    window.addEventListener('mouseup', () => {
      if (!drag) return;
      const moved = drag.moved;
      drag = null;
      const wrap = $('#mapWrap');
      if (wrap) {
        const svg = wrap.querySelector('svg.map');
        if (svg) svg.style.cursor = 'grab';
      }
      if (moved > 4) {
        // a drag must not also read as a click on a feature
        const swallow = e => { e.stopPropagation(); e.preventDefault(); };
        document.addEventListener('click', swallow, true);
        setTimeout(() => document.removeEventListener('click', swallow, true), 0);
      }
    });
  }

  function wireCanvas(wrap, rec, v, W) {
    installGlobal();
    const svg = wrap.querySelector('svg.map');
    if (!svg) return;
    const vbW = W;
    const padL = 10, padR = 10;
    const plotW = W - padL - padR;

    function pxPerUnit() {
      const r = svg.getBoundingClientRect();
      return r.width ? r.width / vbW : 1;
    }
    function bpAt(clientX) {
      const r = svg.getBoundingClientRect();
      const u = (clientX - r.left) / (r.width || 1) * vbW;       // viewBox units
      return v.start + (u - padL) / plotW * v.span;
    }

    svg.addEventListener('wheel', e => {
      e.preventDefault();
      /* A trackpad pinch sends a handful of units per event; a mouse wheel
         sends a hundred or more. One exponent cannot serve both, so the
         rate follows the size of the step and either feels immediate. */
      const dRaw = e.deltaY * (e.deltaMode === 1 ? 16 : 1);
      const rate = Math.abs(dRaw) < 25 ? 0.035 : 0.005;
      const factor = Math.pow(1.0022, -dRaw);
      /* and one flick of the wheel cannot be worth ten-fold */
      const zoomFactor = Math.min(3, Math.max(1 / 3, Math.exp(-dRaw * rate)));
      if (v.mode === 'linear') {
        const anchor = bpAt(e.clientX);
        const frac = (anchor - v.start) / v.span;
        v.span = v.span / factor;
        clampLinear(rec, v);
        v.start = anchor - frac * v.span;
        clampLinear(rec, v);
      } else if (e.ctrlKey || e.metaKey) {
        /* a trackpad pinch arrives as a wheel event with ctrlKey set */
        const r = svg.getBoundingClientRect();
        const c = {
          x: (e.clientX - r.left) / (r.width || 1) * vbW,
          y: (e.clientY - r.top) / (r.height || 1) * (svg.viewBox && svg.viewBox.baseVal ? svg.viewBox.baseVal.height : 640)
        };
        zoomCircular(v, v.k * zoomFactor, c, circGeom(vbW));
      } else {
        /* Two fingers up and down turn the plasmid. The step shrinks as the
           circle grows: an angle that moves the rim a few pixels when the
           whole plasmid is in view moves it most of a screen at a hundred
           times that, and what should be a smooth turn becomes a series of
           jumps. Scaling by the square root keeps the motion continuous
           without making a zoomed-in turn crawl. */
        const d = e.deltaY * (e.deltaMode === 1 ? 16 : 1);
        const turn = 0.004 / Math.sqrt(Math.max(1, v.k));
        v.rot = ((v.rot || 0) - d * turn) % (Math.PI * 2);
      }
      scheduleMapRender();
    }, { passive: false });

    svg.style.cursor = 'grab';
    svg.addEventListener('mousedown', e => {
      if (e.button !== 0) return;
      const mini = e.target.closest && e.target.closest('[data-mini]');
      const r = svg.getBoundingClientRect();
      if (mini && v.mode === 'linear') {
        const L = rec.seq.length;
        v.start = (e.clientX - r.left) / (r.width || 1) * vbW;
        v.start = (v.start - padL) / plotW * L - v.span / 2;
        clampLinear(rec, v);
        drag = { kind: 'mini', lx: e.clientX, ly: e.clientY, moved: 0,
                 rectLeft: r.left, padL: padL * (r.width / vbW), plotW: plotW * (r.width / vbW) };
        scheduleMapRender();
        e.preventDefault();
        return;
      }
      svg.style.cursor = 'grabbing';
      drag = v.mode === 'linear'
        ? { kind: 'linear', lx: e.clientX, ly: e.clientY, moved: 0, plotW: plotW * pxPerUnit() }
        : { kind: 'circular', lx: e.clientX, ly: e.clientY, moved: 0, unit: 1 / pxPerUnit() };
      e.preventDefault();
    });

    // double-click zooms in on that point
    svg.addEventListener('dblclick', e => {
      e.preventDefault();
      if (v.mode === 'linear') {
        const anchor = bpAt(e.clientX);
        v.span = v.span / 2.5;
        clampLinear(rec, v);
        v.start = anchor - v.span / 2;
        clampLinear(rec, v);
      } else {
        const r = svg.getBoundingClientRect();
        zoomCircular(v, v.k * 1.8, {
          x: (e.clientX - r.left) / (r.width || 1) * vbW,
          y: (e.clientY - r.top) / (r.height || 1) * 640
        }, circGeom(vbW));
      }
      scheduleMapRender();
    });

    // pinch
    let pinch = null;
    svg.addEventListener('touchstart', e => {
      if (e.touches.length === 2) {
        const [a, b] = [e.touches[0], e.touches[1]];
        pinch = { d: Math.hypot(a.clientX - b.clientX, a.clientY - b.clientY),
                  span: v.span, k: v.k,
                  mid: (a.clientX + b.clientX) / 2,
                  anchor: v.mode === 'linear' ? bpAt((a.clientX + b.clientX) / 2) : null };
      }
    }, { passive: true });
    svg.addEventListener('touchmove', e => {
      if (!pinch || e.touches.length !== 2) return;
      e.preventDefault();
      const [a, b] = [e.touches[0], e.touches[1]];
      const d = Math.hypot(a.clientX - b.clientX, a.clientY - b.clientY);
      if (!pinch.d) return;
      const ratio = d / pinch.d;
      if (v.mode === 'linear') {
        const frac = (pinch.anchor - v.start) / v.span;
        v.span = pinch.span / ratio;
        clampLinear(rec, v);
        v.start = pinch.anchor - frac * v.span;
        clampLinear(rec, v);
      } else {
        const r = svg.getBoundingClientRect();
        zoomCircular(v, pinch.k * ratio, {
          x: (pinch.mid - r.left) / (r.width || 1) * vbW,
          y: 320
        }, circGeom(vbW));
      }
      scheduleMapRender();
    }, { passive: false });
    svg.addEventListener('touchend', () => { pinch = null; }, { passive: true });

    // feature selection
    wrap.querySelectorAll('[data-fi]').forEach(g => g.addEventListener('click', ev => {
      const i = parseInt(g.getAttribute('data-fi'), 10);
      const f = rec.features[i];
      if (f && f.virtual) { PB.selectPrimerFeature(rec, f); PB.renderHeader(); PB.renderInspector(); drawMap(); return; }
      /* ctrl, or cmd on a Mac, adds to the selection instead of replacing
         it — and takes a feature back out if it is already in */
      PB.toggleFeatureSelection(rec, i, ev.ctrlKey || ev.metaKey);
      PB.renderHeader(); PB.renderInspector(); drawMap();
      if (state.view === 'seq') PB.renderViews();
    }));
  }

  /* ============================== CIRCULAR ============================== */
  function circularMap(rec, v, W) {
    const L = rec.seq.length;
    const H = 640, cx = W / 2, cy = H / 2;
    const R = Math.max(90, Math.min(200, Math.min(W, H) / 2 - 70));
    const k = v.k, inv = 1 / k;
    const TAU = Math.PI * 2;
    /* the whole circle can be turned, so that whatever is being worked on
       can be brought round to the top instead of read sideways */
    const rot = v.rot || 0;
    const ang = p => (p / L) * TAU - Math.PI / 2 + rot;
    const pt = (r, a) => ({ x: cx + r * Math.cos(a), y: cy + r * Math.sin(a) });
    function arcD(r, a0, a1) {
      let d = a1 - a0;
      if (d >= TAU) { a1 = a0 + TAU - 0.0001; d = a1 - a0; }
      const p0 = pt(r, a0), p1 = pt(r, a1);
      return 'M' + p0.x.toFixed(2) + ' ' + p0.y.toFixed(2) +
        ' A' + r.toFixed(2) + ' ' + r.toFixed(2) + ' 0 ' + (Math.abs(d) > Math.PI ? 1 : 0) +
        ' ' + (d > 0 ? 1 : 0) + ' ' + p1.x.toFixed(2) + ' ' + p1.y.toFixed(2);
    }
    // every text size and hairline is divided by k, so it holds its size on screen
    const fs = n => (n * inv).toFixed(2);
    const sw = n => (n * inv).toFixed(3);
    let defs = '', s = '';

    /* How much of the screen one base gets. Past a couple of pixels the
       backbone stops being a line and becomes the sequence itself, which is
       what zooming in is for: keep going and the bases are there. */
    const pxPerBase = (R * k * TAU) / L;
    const basesDrawn = pxPerBase >= 1.6;
    /* Far enough in, coloured blocks stop being the clearest thing to draw:
       the sequence itself is, in both strands, with a ruler under it. The
       blocks are for the zooms in between, where a base is a few pixels
       wide and a letter would not fit. */
    const lettersOnly = pxPerBase >= 13;
    const baseBandPx = lettersOnly ? 0 : Math.min(26, 7 + (pxPerBase - 1.6) * 1.1);
    /* once the sequence itself is drawn on the backbone it takes up room:
       everything outside moves out by half the band, everything inside
       moves in by the same, so nothing is drawn over the bases */
    const halfBand = basesDrawn ? baseBandPx / 2 : 0;
    const outR = R + halfBand * inv;
    const inR = R - halfBand * inv;
    if (basesDrawn) {
      /* only the arc on screen is drawn, or a 10 kb plasmid would be ten
         thousand coloured slivers most of which are outside the canvas */
      const vis = { l: -v.x / k, r: (W - v.x) / k, t: -v.y / k, b: (H - v.y) / k };
      const bw = TAU / L;
      const seq = rec.seq;
      /* Only the bases on screen are worth building, and only they should be
         walked over to find out: testing all ten thousand of them each time
         the wheel moves is what makes a deep zoom stutter. The visible arc
         is worked out from the geometry and the loop runs over that. */
      const viewC = { x: (W / 2 - v.x) / k, y: (H / 2 - v.y) / k };
      const halfDiag = Math.sqrt(W * W + H * H) / 2 / k;
      let i0 = 0, i1 = L - 1;
      if (halfDiag < R * 0.98) {
        const aC = Math.atan2(viewC.y - cy, viewC.x - cx);
        const halfAng = Math.asin(Math.min(1, halfDiag / R)) * 1.25 + 0.01;
        const baseAt = a => ((a + Math.PI / 2 - rot) / TAU) * L;
        i0 = Math.floor(baseAt(aC - halfAng)) - 2;
        i1 = Math.ceil(baseAt(aC + halfAng)) + 2;
        if (i1 - i0 > L) { i0 = 0; i1 = L - 1; }
      }
      const marg = 30 * inv;
      /* One path per base is a thousand stroked elements to rasterise, and
         on a wide screen at two device pixels to one that is what makes a
         zoom crawl. The bases are gathered into one path per colour
         instead: the same picture out of four elements. */
      const byBase = { A: '', C: '', G: '', T: '' };
      /* the band grows with the zoom, so that bases are a readable height
         rather than a coloured hairline, and the letters go in once each
         base has the width to hold one */
      const bandPx = baseBandPx;
      const letters = pxPerBase >= 9;
      const ltrPx = Math.min(15, Math.max(9, pxPerBase * 0.5));
      const comp = { A: 'T', T: 'A', G: 'C', C: 'G' };
      let bases = '', text = '';
      for (let j = i0; j <= i1; j++) {
        const i = ((j % L) + L) % L;
        const a0 = ang(j), a1 = ang(j + 1);
        const mid = (a0 + a1) / 2;
        const p = pt(R, mid);
        if (p.x < vis.l - marg || p.x > vis.r + marg || p.y < vis.t - marg || p.y > vis.b + marg) continue;
        const deg = (mid * 180 / Math.PI + 90).toFixed(2);
        if (lettersOnly) {
          /* both strands, each base in its own colour, nothing behind them */
          const top = pt(R + 13 * inv, mid), bot = pt(R - 12 * inv, mid);
          text += '<text class="b-' + seq[i] + ' mapseqltr" style="font-size:' + fs(ltrPx) + 'px" x="' +
            top.x.toFixed(1) + '" y="' + top.y.toFixed(1) + '" text-anchor="middle" dominant-baseline="central"' +
            ' transform="rotate(' + deg + ' ' + top.x.toFixed(1) + ' ' + top.y.toFixed(1) + ')">' + seq[i] + '</text>' +
            '<text class="b-' + (comp[seq[i]] || 'N') + ' mapseqltr" style="font-size:' + fs(ltrPx) + 'px" x="' +
            bot.x.toFixed(1) + '" y="' + bot.y.toFixed(1) + '" text-anchor="middle" dominant-baseline="central"' +
            ' transform="rotate(' + deg + ' ' + bot.x.toFixed(1) + ' ' + bot.y.toFixed(1) + ')">' +
            (comp[seq[i]] || 'N') + '</text>';
          continue;
        }
        const d = arcD(R, a0, a1 + bw * 0.08);
        if (byBase[seq[i]] != null) byBase[seq[i]] += d; else bases += '<path class="b-N" stroke="currentColor" stroke-width="' +
          sw(bandPx) + '" fill="none" d="' + d + '"/>';
        if (letters) {
          text += '<text class="mapbaseltr" style="font-size:' + fs(Math.min(15, bandPx * 0.62)) + 'px" x="' +
            p.x.toFixed(1) + '" y="' + p.y.toFixed(1) + '" text-anchor="middle" dominant-baseline="central"' +
            ' transform="rotate(' + deg + ' ' + p.x.toFixed(1) + ' ' + p.y.toFixed(1) + ')">' + seq[i] + '</text>';
        }
      }
      Object.keys(byBase).forEach(b => {
        if (byBase[b]) bases += '<path class="b-' + b + '" stroke="currentColor" stroke-width="' +
          sw(bandPx) + '" fill="none" d="' + byBase[b] + '"/>';
      });
      s += '<g class="mapbases" style="stroke-linecap:butt">' + bases + '</g>' + text;
    }
    /* with both strands written out there is nothing for a backbone line to
       say, and it would only run between them */
    if (!lettersOnly) {
      s += '<circle class="backbone" cx="' + cx + '" cy="' + cy + '" r="' + R +
           '" stroke-width="' + sw(basesDrawn ? 0.8 : 2.5) + '"' + (basesDrawn ? ' opacity="0.35"' : '') + '/>';
    }
    /* ticks spaced by what the screen can hold, not by the plasmid's size:
       every 500 bases on the whole circle, every 5 among the letters */
    const step = basesDrawn ? Math.max(1, niceStep(900 / pxPerBase)) : niceStep(L / Math.max(k, 1));
    /* Zoomed in among the bases the ruler is marked every few bases, and
       stepping that way round the whole plasmid is a thousand ticks for the
       forty that can be seen. The loop is held to the arc on screen. */
    const tickVis = { l: -v.x / k, r: (W - v.x) / k, t: -v.y / k, b: (H - v.y) / k };
    const tickMarg = 40 * inv;
    for (let p = 0; p < L; p += step) {
      const a = ang(p);
      if (k > 2) {
        const tp = pt(R, a);
        if (tp.x < tickVis.l - tickMarg || tp.x > tickVis.r + tickMarg ||
            tp.y < tickVis.t - tickMarg || tp.y > tickVis.b + tickMarg) continue;
      }
      const major = basesDrawn ? true : (Math.round(p / step) % 5 === 0);
      const rulerR = lettersOnly ? R - 30 * inv : outR + 1 * inv;
      const i0 = lettersOnly ? pt(rulerR - 7 * inv, a)
                             : (basesDrawn ? pt(outR + (major ? 8 : 4) * inv, a) : pt(R - (major ? 6 : 3) * inv, a));
      const i1 = lettersOnly ? pt(rulerR, a) : (basesDrawn ? pt(outR + 1 * inv, a) : pt(R, a));
      s += '<line class="' + (lettersOnly ? 'rulertick' : (major ? 'tickmaj' : 'tick')) + '" stroke-width="' +
           sw(lettersOnly ? 1.4 : 1) + '" x1="' +
           i0.x.toFixed(1) + '" y1="' + i0.y.toFixed(1) + '" x2="' + i1.x.toFixed(1) + '" y2="' + i1.y.toFixed(1) + '"/>';
      if (major) {
        /* with the bases drawn the coordinate sits outside, clear of them */
        const t = lettersOnly ? pt(R - 46 * inv, a)
                              : (basesDrawn ? pt(outR + 15 * inv, a) : pt(R - 15 * inv, a));
        s += '<text class="ticklbl' + (lettersOnly ? ' rulerlbl' : '') + '" style="font-size:' + fs(lettersOnly ? 10 : 9) +
             'px" x="' + t.x.toFixed(1) + '" y="' + t.y.toFixed(1) +
             '" text-anchor="middle" dominant-baseline="middle">' + (p === 0 ? '1' : fmt(p)) + '</text>';
      }
    }
    /* every selected feature gets its own arc when several are picked,
       so ctrl-clicking three of them marks three stretches of the ring
       rather than only the last one */
    const multi = PB.featureSelectionPieces(rec);
    if (multi) {
      multi.forEach(p => {
        if (p.end <= p.start) return;
        s += '<path class="selarc" stroke-width="' + sw(basesDrawn ? baseBandPx + 4 : 7) + '" d="' +
          arcD(R, ang(p.start), ang(p.end)) + '"/>';
      });
    } else {
      const nsel = PB.normSel(rec, state.selection);
      if (nsel) {
        /* a selection across base 1 simply carries on past the top */
        s += '<path class="selarc" stroke-width="' + sw(basesDrawn ? baseBandPx + 4 : 7) + '" d="' +
          arcD(R, ang(nsel.start), ang(nsel.end)) + '"/>';
      }
    }

    let outer = [];
    const ivs = PB.featIntervals(rec);
    /* feature bands hold their thickness on screen, like the text does:
       scaled with the circle they would be a thousand pixels deep at the
       zoom where the bases become readable */
    const BAND = 12 * inv, GAP = 2.5 * inv;
    ivs.forEach(it => {
      const r = (lettersOnly ? R - 64 * inv : inR - 11 * inv) - it.track * (BAND + GAP);
      if (r < R * 0.2) return;
      const color = PB.featColor(it.f);
      const dir = it.f.strand === -1 ? -1 : 1;
      const isSel = PB.isFeatureSelected(rec, it.idx);
      s += '<g class="feat' + (isSel ? ' sel' : '') + '" data-fi="' + it.idx + '">' +
        '<title>' + esc(it.f.label + ' - ' + it.f.type + ' - ' + PB.locationText(rec, it.f) +
          ' - ' + fmt(it.bpLen) + ' bp - ' + (it.f.mixed ? 'mixed strands' : (dir === -1 ? 'reverse' : 'forward'))) + '</title>';
      // introns: a hairline arc joining the pieces
      PB.gapsOf(it, L).forEach(gp => {
        if (gp.end > gp.start) s += '<path d="' + arcD(r, ang(gp.start), ang(gp.end)) + '" stroke="' + color +
          '" stroke-width="' + sw(1.6) + '" fill="none" opacity=".75"/>';
      });
      // the arrowhead sits on the piece that holds the feature's 3' end
      const arrowAt = it.f.mixed ? -1 : (dir === 1 ? it.pieces.length - 1 : 0);
      let best = null;
      it.pieces.forEach((p, pi) => {
        const a0 = ang(p.start), a1 = ang(p.end);
        const pdir = it.f.mixed ? ((it.f.parts[pi] || {}).strand === -1 ? -1 : 1) : dir;
        const arrow = it.f.mixed || pi === arrowAt;
        const sweep = Math.abs(a1 - a0);
        const headLen = arrow ? Math.min(11, r * sweep * 0.4) : 0;
        const shrink = headLen / r;
        const b0 = pdir === 1 ? a0 : a0 + shrink;
        const b1 = pdir === 1 ? a1 - shrink : a1;
        if (b1 > b0) s += '<path d="' + arcD(r, b0, b1) + '" stroke="' + color + '" stroke-width="' + BAND + '" fill="none"/>';
        if (arrow && r * sweep > 3) {
          const tipA = pdir === 1 ? a1 : a0;
          const baseA = pdir === 1 ? a1 - shrink : a0 + shrink;
          const tip = pt(r, tipA), h1 = pt(r + BAND / 2, baseA), h2 = pt(r - BAND / 2, baseA);
          s += '<polygon points="' + tip.x.toFixed(1) + ',' + tip.y.toFixed(1) + ' ' +
            h1.x.toFixed(1) + ',' + h1.y.toFixed(1) + ' ' + h2.x.toFixed(1) + ',' + h2.y.toFixed(1) +
            '" fill="' + color + '"/>';
        }
        if (!best || sweep > best.sweep) best = { a0: a0, a1: a1, sweep: sweep };
      });
      if (best) {
        const txtW = it.f.label.length * 5.9 + 10;          // screen px
        if (r * best.sweep * k > txtW && BAND * k > 9) {
          const mid = (best.a0 + best.a1) / 2;
          const flip = Math.sin(mid) > 0;
          const id = 'fp' + String(rec.id).replace(/[^a-z0-9]/gi, '') + it.idx;
          defs += '<path id="' + id + '" d="' + (flip ? arcD(r, best.a1, best.a0) : arcD(r, best.a0, best.a1)) + '" fill="none"/>';
          s += '<text class="featlbl" fill="#fff" style="font-size:' + fs(10.5) + 'px" dy="' + (3.7 * inv).toFixed(2) +
            '"><textPath href="#' + id + '" xlink:href="#' + id + '" startOffset="50%" text-anchor="middle">' +
            esc(it.f.label) + '</textPath></text>';
        } else {
          outer.push({ ang: (best.a0 + best.a1) / 2, text: it.f.label, color: color, srcR: r + BAND / 2, idx: it.idx, kind: 'feat',
                       pri: (PB.isFeatureSelected(rec, it.idx) ? 1e12 : 0) + it.bpLen });
        }
      }
      s += '</g>';
    });

    const cuts = state.enzymes.length ? PB.allCuts(rec) : [];
    const tooMany = cuts.length > 70;
    cuts.forEach(c => {
      const a = ang(c.cutTop);
      if (k > 2) {
        const cp = pt(R, a);
        if (cp.x < tickVis.l - tickMarg || cp.x > tickVis.r + tickMarg ||
            cp.y < tickVis.t - tickMarg || cp.y > tickVis.b + tickMarg) return;
      }
      const p0 = pt(outR + (basesDrawn ? 20 : 1) * inv, a), p1 = pt(outR + (basesDrawn ? 28 : 9) * inv, a);
      s += '<line class="enztick" stroke-width="' + sw(1) + '" x1="' + p0.x.toFixed(1) + '" y1="' + p0.y.toFixed(1) +
           '" x2="' + p1.x.toFixed(1) + '" y2="' + p1.y.toFixed(1) + '"/>';
    });
    /* Sites a few bases apart cannot each have a label: on a 5 kb plasmid a
       polylinker puts six of them inside two degrees of arc, and six leader
       lines to six labels is unreadable. They are named together instead,
       the way a printed map does it, and zooming in separates them because
       the grouping is measured on screen. */
    if (!tooMany && cuts.length) {
      const sorted = cuts.slice().sort((x, y) => x.cutTop - y.cutTop);
      const minSep = 13 / Math.max(1, R * k / 1);            // in radians at this zoom
      const groups = [];
      sorted.forEach(c => {
        const a = ang(c.cutTop);
        const g = groups[groups.length - 1];
        if (g && Math.abs(a - g.lastAng) < minSep) { g.cuts.push(c); g.lastAng = a; }
        else groups.push({ cuts: [c], lastAng: a });
      });
      /* the first and last may meet across the origin */
      if (groups.length > 1) {
        const first = groups[0], last = groups[groups.length - 1];
        if (Math.abs((ang(first.cuts[0].cutTop) + TAU) - last.lastAng) < minSep) {
          last.cuts = last.cuts.concat(first.cuts);
          groups.shift();
        }
      }
      groups.forEach(g => {
        const mid = g.cuts.reduce((t, c) => t + ang(c.cutTop), 0) / g.cuts.length;
        const names = g.cuts.map(c => c.enzyme + PB.cutMark(c));
        const blocked = g.cuts.every(c => c.methyl && c.methyl.effect === 'blocked');
        let text;
        if (g.cuts.length === 1) text = names[0] + ' (' + fmt(g.cuts[0].cutTop) + ')';
        else if (g.cuts.length <= 4) text = names.join(', ') + ' (' + fmt(g.cuts[0].cutTop) + '\u2013' + fmt(g.cuts[g.cuts.length - 1].cutTop) + ')';
        else text = names.slice(0, 3).join(', ') + ' +' + (g.cuts.length - 3) + ' (' +
          fmt(g.cuts[0].cutTop) + '\u2013' + fmt(g.cuts[g.cuts.length - 1].cutTop) + ')';
        outer.push({ ang: mid, text: text, color: blocked ? 'var(--fg-faint)' : null,
                     srcR: outR + (basesDrawn ? 28 : 9) * inv, kind: 'enz', pri: 1e9 + g.cuts.length });
      });
    }

    // More labels than a side can hold would be pushed off the canvas, so keep
    // the most important (selected feature, enzyme sites, then the longest
    // features) and leave the rest for zooming in or for hovering the arc.
    /* ---- zoomed in, the arc is nearly a straight line ----
       Columns down the side of the canvas are right for a whole circle and
       wrong for a shallow arc: the labels end up a screen away from the
       features they name. Once only a slice of the plasmid is in view they
       go above the arc instead, stacked into as few rows as they fit, each
       on a short upright leader over its own site — which is how a
       sequence viewer shows them, and how they stay readable. */
    const visArc = W / Math.max(1, R * k);
    const flatMode = visArc < 1.2 && outer.length;
    let hiddenLabels = 0;
    if (flatMode) {
      const vis = { l: -v.x / k, r: (W - v.x) / k, t: -v.y / k, b: (H - v.y) / k };
      const charW = 5.9 * inv, lineH = 13.5 * inv;
      const keep = [];
      outer.forEach(e => {
        const p = pt(e.srcR, e.ang);
        /* Both axes have to be tested, not just the horizontal one. Zoomed
           in, the ring is enormous and the canvas cuts a narrow slice out
           of it: for any x on the visible arc there is a second point with
           the same x on the far side of the circle, thousands of units
           above or below. Testing x alone let those through, so a label
           from the opposite end of the plasmid would appear over the arc
           with a leader running off the canvas to a feature nowhere near. */
        if (p.x < vis.l - 10 * inv || p.x > vis.r + 10 * inv ||
            p.y < vis.t - 10 * inv || p.y > vis.b + 10 * inv) { hiddenLabels++; return; }
        e.ax = p.x; e.ay = p.y;
        e.w = ((e.text || '').length) * charW + 12 * inv;
        keep.push(e);
      });
      keep.sort((a, b) => a.ax - b.ax);
      const rowEnd = [];
      keep.forEach(e => {
        let r = 0;
        while (rowEnd[r] != null && rowEnd[r] > e.ax - e.w / 2) r++;
        rowEnd[r] = e.ax + e.w / 2 + 8 * inv;
        e.row = r;
      });
      const topY = keep.reduce((m, e) => Math.min(m, e.ay), Infinity) - 26 * inv;
      keep.forEach(e => {
        const ly = topY - e.row * lineH;
        s += '<path class="leader" fill="none" stroke-width="' + sw(0.8) + '" d="M' +
             e.ax.toFixed(1) + ' ' + e.ay.toFixed(1) + ' L' + e.ax.toFixed(1) + ' ' + (ly + 5 * inv).toFixed(1) + '"/>';
        /* centred over its own site, except at the edges of the canvas,
           where a centred name would hang half off: there it reads inward
           from the site instead, which keeps the whole of it legible and
           still plainly attached to the right place */
        const half = ((e.text || '').length) * 6.4 * inv / 2 + 6 * inv;
        const anchor = (e.ax + half > vis.r - 4 * inv) ? 'end'
                     : (e.ax - half < vis.l + 4 * inv) ? 'start' : 'middle';
        s += '<text class="' + (e.kind === 'enz' ? 'enzlbl' : 'featlbl') + '" style="font-size:' +
             fs(e.kind === 'enz' ? 9.5 : 10.5) + 'px' + (e.idx != null ? ';cursor:pointer' : '') + '"' +
             (e.color ? ' fill="' + e.color + '"' : '') +
             (e.idx != null ? ' data-fi="' + e.idx + '"' : '') +
             ' x="' + e.ax.toFixed(1) + '" y="' + ly.toFixed(1) +
             '" text-anchor="' + anchor + '" dominant-baseline="middle">' + esc(e.text) + '</text>';
      });
      outer = [];
    }
    /* Zoomed in but not far enough for the flat rows above, the labels are
       still laid out around the ring - and that layout drew every one of
       them, wherever its feature sat. On a view holding a sixth of the
       plasmid that meant names from the opposite side stacked into the
       gutter with leaders running clear off the canvas to features
       nowhere near. A label is only worth drawing when the thing it names
       is somewhere on screen. */
    if (outer.length && k > 1.05) {
      const vis = { l: -v.x / k, r: (W - v.x) / k, t: -v.y / k, b: (H - v.y) / k };
      const marg = 10 * inv;
      const kept = [];
      outer.forEach(e => {
        const p = pt(e.srcR, e.ang);
        if (p.x < vis.l - marg || p.x > vis.r + marg ||
            p.y < vis.t - marg || p.y > vis.b + marg) { hiddenLabels++; return; }
        kept.push(e);
      });
      outer = kept;
    }
    if (outer.length) {
      const cap = Math.max(4, (Math.floor(H / (12.4 * inv)) - 2) * 3);
      [1, -1].forEach(side => {
        const grp = outer.filter(e => (Math.cos(e.ang) >= 0 ? 1 : -1) === side);
        if (grp.length <= cap) return;
        grp.sort((a, b) => b.pri - a.pri).slice(cap).forEach(e => { e.drop = true; hiddenLabels++; });
      });
      if (hiddenLabels) outer = outer.filter(e => !e.drop);
    }
    if (outer.length) {
      // ring offset and line spacing are in map units, so divide by k to hold
      // a constant gap on screen as the circle grows
      const room = Math.max(60 * inv, (W / 2) - R * 0.46 - 10 * inv);
      layoutLabels(outer, cx, cy, R + 26 * inv, 12.4 * inv, H, room, 5.9 * inv);
      if (k > 1.05) {
        /* the layout can push a row past the top or bottom of the canvas,
           where it is drawn and never seen; counting those as hidden keeps
           the "labels outside this view" note honest */
        const lv = { l: -v.x / k, r: (W - v.x) / k, t: -v.y / k, b: (H - v.y) / k };
        const pad = 8 * inv;
        const seen = outer.filter(e => {
          /* a row past the top or bottom, or a column past either side */
          const vis2 = e.ly > lv.t - pad && e.ly < lv.b + pad &&
                       e.lx > lv.l - pad && e.lx < lv.r + pad;
          if (!vis2) hiddenLabels++;
          return vis2;
        });
        /* The column sits at a fixed distance from the centre of the
           circle, which once zoomed can put it hard against the edge of
           the canvas with the text running off it. Half a name is worse
           than none, so each column is pulled back in far enough for its
           own text to fit. */
        seen.forEach(e => {
          /* a little wider than the nominal character cell: the leader
             meets the text 4 units out, and some glyphs run over */
          const w = ((e.text || '').length) * 6.4 * inv + 14 * inv;
          if (e.side === 1) e.lx = Math.min(e.lx, lv.r - w - pad);
          else e.lx = Math.max(e.lx, lv.l + w + pad);
        });
        outer = seen;
      }
      outer.forEach(e => {
        const anchor = e.side === 1 ? 'start' : 'end';
        const ex = e.side === 1 ? e.lx + 4 * inv : e.lx - 4 * inv;
        const src = pt(e.srcR, e.ang);
        /* a short radial stub before the run to the label: two straight
           pieces are far easier to follow through a crowd than one long
           diagonal crossing a dozen others */
        const knee = pt(e.srcR + 10 * inv, e.ang);
        s += '<path class="leader" fill="none" stroke-width="' + sw(0.8) + '" d="M' + src.x.toFixed(1) + ' ' + src.y.toFixed(1) +
             ' L' + knee.x.toFixed(1) + ' ' + knee.y.toFixed(1) +
             ' L' + ex.toFixed(1) + ' ' + e.ly.toFixed(1) + '"/>';
        s += '<text class="' + (e.kind === 'enz' ? 'enzlbl' : 'featlbl') + '" style="font-size:' +
             fs(e.kind === 'enz' ? 9.5 : 10.5) + 'px' + (e.idx != null ? ';cursor:pointer' : '') + '"' +
             (e.color ? ' fill="' + e.color + '"' : '') +
             (e.idx != null ? ' data-fi="' + e.idx + '"' : '') +
             ' x="' + ex.toFixed(1) + '" y="' + e.ly.toFixed(1) +
             '" text-anchor="' + anchor + '" dominant-baseline="middle">' + esc(e.text) + '</text>';
      });
    }

    s += '<text class="ctrname" style="font-size:' + fs(15) + 'px" x="' + cx + '" y="' + (cy - 5 * inv) + '">' +
         esc(rec.name) + '</text>';
    s += '<text class="ctrmeta" style="font-size:' + fs(10.5) + 'px" x="' + cx + '" y="' + (cy + 13 * inv) + '">' +
         fmt(L) + ' bp - ' + rec.topology + '</text>';
    s += '<text class="ctrmeta" style="font-size:' + fs(10.5) + 'px" x="' + cx + '" y="' + (cy + 27 * inv) + '">' +
         f1(B.gcContent(rec.seq)) + '% GC</text>';

    let head = '';
    if (tooMany) head = '<text class="enzlbl" x="' + (W - 8) + '" y="16" text-anchor="end">' +
      cuts.length + ' sites - open the Enzymes tab for the list</text>';
    if (hiddenLabels) head += '<text class="maplegend" x="8" y="' + (H - 10) + '" text-anchor="start">' +
      hiddenLabels + (flatMode ? ' label' + (hiddenLabels === 1 ? '' : 's') + ' outside this view - pan or zoom out'
                               : ' smaller feature labels hidden - zoom in or hover the map') + '</text>';

    return '<svg class="map" viewBox="0 0 ' + W + ' ' + H + '" xmlns:xlink="http://www.w3.org/1999/xlink" ' +
      'role="img" aria-label="Plasmid map of ' + esc(rec.name) + '"><defs>' + defs + '</defs>' + head +
      '<g id="mapZoomLayer" transform="translate(' + v.x.toFixed(2) + ' ' + v.y.toFixed(2) + ') scale(' + v.k.toFixed(4) + ')">' +
      s + '</g></svg>';
  }

  /* ============================== LINEAR ============================== */
  function linearMap(rec, v, W) {
    const L = rec.seq.length;
    clampLinear(rec, v);
    const padL = 10, padR = 10;
    const plotW = W - padL - padR;
    const pxPerBp = plotW / v.span;
    const x = bp => padL + (bp - v.start) / v.span * plotW;
    const from = Math.max(0, Math.floor(v.start) - 1);
    const to = Math.min(L, Math.ceil(v.start + v.span) + 1);

    const showSeq = pxPerBp >= SEQ_PX;
    const showComp = pxPerBp >= COMP_PX;
    const showAA = pxPerBp >= AA_PX;
    // a monospace glyph advances about 0.6 em, so size the type to the space
    // each base actually has rather than hiding the sequence altogether
    const seqFont = Math.max(7, Math.min(13, pxPerBp * 1.55));

    // enzyme labels for sites in view only
    const cuts = (state.enzymes.length ? PB.allCuts(rec) : []).filter(c => c.cutTop >= from && c.cutTop <= to);
    const allEnzRows = assignRows(cuts.map(c => ({ pos: c.cutTop, text: c.enzyme + PB.cutMark(c) + ' (' + fmt(c.cutTop) + ')' })), x, 6.1);
    const ENZ_ROWS = 2, SEQ_BAND = 46, MINI = 26;
    const enzRows = allEnzRows.filter(r => r.row < ENZ_ROWS);
    const hiddenEnz = allEnzRows.length - enzRows.length;

    // The panel keeps a constant height whatever the zoom, so the drawing never
    // shifts out from under the pointer mid-gesture.
    const enzTop = MINI + 8;
    const axisY = enzTop + ENZ_ROWS * 13 + 12;
    const seqY = axisY + 22;
    const featTop = seqY + SEQ_BAND;
    const allIvs = PB.featIntervals(rec);
    const inView = it => it.pieces.some(p => p.end > from && p.start < to) ||
      PB.gapsOf(it, L).some(g => g.end > from && g.start < to);
    const ivs = allIvs.filter(inView);
    // Rows are packed by the room a feature needs INCLUDING its name, so a
    // label never lands on a neighbour; a feature that wraps the origin
    // occupies both ends of the strip. Height is reserved from the whole-
    // construct view, the most crowded case, so the panel never changes size.
    function rangesOf(it, xf) {
      const txtW = it.f.label.length * 5.9 + 10;
      if (it.end <= L) {
        const a = xf(it.start), c = xf(it.end);
        return [[a, (c - a) > txtW ? c : c + txtW + 8]];
      }
      const a = xf(it.start), c = xf(it.end - L);
      return [[a, xf(L)], [xf(0), c + ((c - xf(0)) > txtW ? 0 : txtW + 8)]];
    }
    function packRows(list, xf) {
      const sorted = list.slice().sort((a, c) => a.start - c.start);
      const rows = [];
      const out = new Map();
      sorted.forEach(it => {
        const rs = rangesOf(it, xf);
        let r = 0;
        while ((rows[r] || []).some(q => rs.some(z => z[0] < q[1] + 4 && q[0] < z[1] + 4))) r++;
        out.set(it.idx, r);
        (rows[r] = rows[r] || []).push.apply(rows[r], rs);
      });
      return { rows: out, count: Math.max(rows.length, 1) };
    }
    const fullX = bp => padL + (bp / Math.max(L, 1)) * plotW;
    const maxRows = packRows(allIvs, fullX).count;
    const packed = packRows(ivs, x);
    const H = featTop + maxRows * 20 + 16;

    let s = '<svg class="map" viewBox="0 0 ' + W + ' ' + H + '" role="img" aria-label="Linear map of ' + esc(rec.name) + '">';

    /* --- overview strip --- */
    s += '<g data-mini="1" style="cursor:pointer">';
    s += '<rect x="' + padL + '" y="6" width="' + plotW + '" height="' + (MINI - 12) + '" rx="3" fill="var(--surface-2)"/>';
    allIvs.forEach(it => it.pieces.forEach(p => {
      const mx = padL + (p.start / L) * plotW;
      const mw = Math.max(1, ((p.end - p.start) / L) * plotW);
      s += '<rect x="' + mx.toFixed(1) + '" y="8" width="' + mw.toFixed(1) + '" height="' + (MINI - 16) +
           '" fill="' + PB.featColor(it.f) + '" opacity=".75"/>';
    }));
    const wx = padL + (v.start / L) * plotW;
    const ww = Math.max(2, (v.span / L) * plotW);
    s += '<rect x="' + wx.toFixed(1) + '" y="4" width="' + ww.toFixed(1) + '" height="' + (MINI - 8) +
         '" fill="var(--accent)" opacity=".20" stroke="var(--accent)" stroke-width="1" rx="2"/>';
    s += '</g>';

    /* --- enzyme ticks --- */
    enzRows.forEach(e => {
      const px = x(e.pos);
      const ly = enzTop + e.row * 13;
      s += '<line class="enztick" x1="' + px.toFixed(1) + '" y1="' + ly + '" x2="' + px.toFixed(1) +
           '" y2="' + axisY + '" stroke-dasharray="2 2"/>';
      s += '<text class="enzlbl" x="' + px.toFixed(1) + '" y="' + (ly - 3) + '" text-anchor="middle">' + esc(e.text) + '</text>';
    });
    if (hiddenEnz > 0) {
      s += '<text class="enzlbl" x="' + (W - padR) + '" y="' + (enzTop + 8) + '" text-anchor="end">+' +
        hiddenEnz + ' more site' + (hiddenEnz === 1 ? '' : 's') + ' — zoom in</text>';
    }

    /* --- selection --- */
    const nsl = PB.normSel(rec, state.selection);
    const multiL = PB.featureSelectionPieces(rec);
    if (multiL || nsl) {
      /* on the linear view a selection across base 1 is two pieces, and
         several selected features are several more */
      const L0 = rec.seq.length;
      const pieces = multiL ? multiL.map(p => [p.start, p.end])
        : (nsl.wraps ? [[nsl.start, L0], [0, nsl.end - L0]] : [[nsl.start, nsl.end]]);
      pieces.forEach(([p0, p1]) => {
        const sx = x(p0), sw2 = Math.max(1.5, x(p1) - sx);
        s += '<rect x="' + sx.toFixed(1) + '" y="' + (axisY - 8) + '" width="' + sw2.toFixed(1) +
             '" height="' + (H - axisY + 2) + '" fill="var(--sel-strong)" opacity=".18"/>';
      });
    }

    /* --- axis --- */
    s += '<line class="backbone" x1="' + padL + '" y1="' + axisY + '" x2="' + (W - padR) + '" y2="' + axisY + '"/>';
    const step = niceStep(v.span);
    const first = Math.ceil(from / step) * step;
    for (let p = first; p <= to; p += step) {
      const px = x(p);
      if (px < padL - 2 || px > W - padR + 2) continue;
      s += '<line class="tick" x1="' + px.toFixed(1) + '" y1="' + axisY + '" x2="' + px.toFixed(1) + '" y2="' + (axisY + 5) + '"/>';
      s += '<text class="ticklbl" x="' + px.toFixed(1) + '" y="' + (axisY + 15) + '" text-anchor="middle">' +
           fmt(p === 0 ? 1 : p) + '</text>';
    }

    /* --- bases, once there is room for them --- */
    if (showSeq) {
      const top = rec.seq.slice(from, to);
      let tRow = '', cRow = '';
      for (let i = 0; i < top.length; i++) {
        const bx = x(from + i + 0.5);
        if (bx < -10 || bx > W + 10) continue;
        const base = top[i];
        tRow += '<tspan class="b-' + base + '" x="' + bx.toFixed(1) + '">' + base + '</tspan>';
        if (showComp) cRow += '<tspan class="b-' + (B.complement(base)) + '" x="' + bx.toFixed(1) + '">' + B.complement(base) + '</tspan>';
      }
      s += '<text class="seqrow" style="font-size:' + seqFont.toFixed(1) + 'px" y="' + (seqY + 12) +
           '" text-anchor="middle">' + tRow + '</text>';
      if (showComp) s += '<text class="seqrow dim" style="font-size:' + seqFont.toFixed(1) + 'px" y="' +
           (seqY + 26) + '" text-anchor="middle">' + cRow + '</text>';
      if (showAA) {
        const aaY = seqY + 26 + (showComp ? 15 : 0);
        ivs.forEach(it => {
          if (!/CDS|ORF/i.test(it.f.type) && !(it.f.type === 'gene' && !rec.features.some(g => g.type === 'CDS'))) return;
          let row = '';
          PB.codonMap(rec, it.idx).forEach(c => {
            const mid = c.pos[1];
            if (mid < from || mid >= to) return;
            row += '<tspan x="' + x(mid + 0.5).toFixed(1) + '"' +
              (c.aa === '*' ? ' fill="var(--bad)"' : '') + '>' + esc(c.aa) + '</tspan>';
          });
          if (row) s += '<text class="aarow" style="font-size:' + Math.max(7, Math.min(11, seqFont)).toFixed(1) +
            'px" y="' + aaY + '" text-anchor="middle">' + row + '</text>';
        });
      }
    }

    if (!showSeq) {
      s += '<text class="ticklbl" x="' + (W / 2) + '" y="' + (seqY + 26) +
        '" text-anchor="middle" opacity=".7">zoom in to read the sequence (' +
        fmt(Math.round(v.span)) + ' bp in view)</text>';
    }

    /* --- features ---
       Shapes first, then every label in a second pass, so a later feature's
       fill can never paint over an earlier feature's name. */
    let shapes = '', labels = '';
    const lbls = [];
    // occupancy per track, so no two names ever sit on top of one another
    const occ = {};
    const free = (t, a, b) => (occ[t] || []).every(r => a >= r[1] + 3 || b + 3 <= r[0]);
    const mark = (t, a, b) => { (occ[t] = occ[t] || []).push([a, b]); };
    ivs.forEach(it => {
      const row = packed.rows.get(it.idx) || 0;
      const y = featTop + row * 20;
      const color = PB.featColor(it.f);
      const dir = it.f.strand === -1 ? -1 : 1;
      const arrowAt = it.f.mixed ? -1 : (dir === 1 ? it.pieces.length - 1 : 0);
      let g = '<g class="feat' + (PB.isFeatureSelected(rec, it.idx) ? ' sel' : '') + '" data-fi="' + it.idx + '">' +
        '<title>' + esc(it.f.label + ' - ' + it.f.type + ' - ' + PB.locationText(rec, it.f) + ' - ' + fmt(it.bpLen) + ' bp') + '</title>';
      PB.gapsOf(it, L).forEach(gp => {
        const a = Math.max(x(gp.start), -20), b = Math.min(x(gp.end), W + 20);
        if (b > a) g += '<line x1="' + a.toFixed(1) + '" y1="' + (y + 6.5) + '" x2="' + b.toFixed(1) + '" y2="' + (y + 6.5) +
          '" stroke="' + color + '" stroke-width="1.5" opacity=".8"/>';
      });
      let best = null;
      it.pieces.forEach((p, pi) => {
        if (!(p.end > from && p.start < to)) return;
        const x0 = Math.max(x(p.start), -20), x1 = Math.min(x(p.end), W + 20);
        const w = Math.max(2, x1 - x0);
        const pdir = it.f.mixed ? ((it.f.parts[pi] || {}).strand === -1 ? -1 : 1) : dir;
        const head = (it.f.mixed || pi === arrowAt) ? Math.min(8, w * 0.4) : 0;
        const d = pdir === 1
          ? 'M' + x0 + ' ' + y + ' H' + (x0 + w - head) + ' L' + (x0 + w) + ' ' + (y + 6.5) +
            ' L' + (x0 + w - head) + ' ' + (y + 13) + ' H' + x0 + ' Z'
          : 'M' + (x0 + w) + ' ' + y + ' H' + (x0 + head) + ' L' + x0 + ' ' + (y + 6.5) +
            ' L' + (x0 + head) + ' ' + (y + 13) + ' H' + (x0 + w) + ' Z';
        g += '<path d="' + d + '" fill="' + color + '"/>';
        if (!best || w > best.w) best = { x0: x0, x1: x1, w: w };
      });
      shapes += g + '</g>';
      if (best) lbls.push({ it: it, y: y, row: row, x0: best.x0, x1: best.x1, w: best.w, color: color,
                            txtW: it.f.label.length * 5.9 + 10 });
    });
    // every visible piece reserves its own span before any name is placed outside one
    ivs.forEach(it => {
      const row = packed.rows.get(it.idx) || 0;
      it.pieces.forEach(p => { if (p.end > from && p.start < to) mark(row, Math.max(x(p.start), -20), Math.min(x(p.end), W + 20)); });
    });
    // names that fit inside their own box are always drawn
    lbls.forEach(L => {
      if (L.w <= L.txtW) return;
      labels += '<text class="featlbl" data-fi="' + L.it.idx + '" style="cursor:pointer" fill="#fff" x="' +
        (L.x0 + L.w / 2).toFixed(1) + '" y="' + (L.y + 9.8) + '" text-anchor="middle">' +
        esc(L.it.f.label) + '</text>';
    });
    // then names that must sit outside, wherever the track is still clear
    lbls.forEach(L => {
      if (L.w > L.txtW) return;
      const rA = L.x1 + 4, rB = rA + L.txtW;
      if (rB < W && free(L.row, rA, rB)) {
        labels += '<text class="featlbl" data-fi="' + L.it.idx + '" style="cursor:pointer" fill="' + L.color +
          '" x="' + rA.toFixed(1) + '" y="' + (L.y + 9.8) + '">' + esc(L.it.f.label) + '</text>';
        mark(L.row, rA, rB);
        return;
      }
      const lB = L.x0 - 4, lA = lB - L.txtW;
      if (lA > 0 && free(L.row, lA, lB)) {
        labels += '<text class="featlbl" data-fi="' + L.it.idx + '" style="cursor:pointer" fill="' + L.color +
          '" x="' + lB.toFixed(1) + '" y="' + (L.y + 9.8) + '" text-anchor="end">' + esc(L.it.f.label) + '</text>';
        mark(L.row, lA, lB);
      }
    });
    s += shapes + labels;
    if (!ivs.length) {
      s += '<text class="ticklbl" x="' + (W / 2) + '" y="' + (featTop + 14) + '" text-anchor="middle">no features in view</text>';
    }
    return s + '</svg>';
  }

  function assignRows(items, x, charPx) {
    items = items.slice().sort((a, b) => a.pos - b.pos);
    const rowEnd = [];
    items.forEach(it => {
      const px = x(it.pos);
      const halfW = (it.text.length * charPx) / 2 + 3;
      let r = 0;
      while (r < rowEnd.length && rowEnd[r] > px - halfW) r++;
      it.row = r;
      rowEnd[r] = px + halfW;
    });
    return items;
  }

  // keep the drawing in step with the window
  let resizeTimer = null;
  window.addEventListener('resize', () => {
    if (state.view !== 'map') return;
    clearTimeout(resizeTimer);
    resizeTimer = setTimeout(() => drawMap(), 140);
  });
  /* ---------------- map image export ----------------
     The map on screen is drawn with the page's stylesheet and theme
     colours, which a file on its own would not have. So the export draws
     the map off screen, in the light theme unless asked otherwise, writes
     every computed style onto each element, and strips everything that only
     makes sense in the page (classes, handlers, hover states, the
     selection). The SVG then looks the same in Illustrator, Inkscape, a
     browser or PowerPoint, and the PNG is drawn from that SVG. */
  const STYLE_PROPS = ['fill', 'fill-opacity', 'stroke', 'stroke-width', 'stroke-opacity',
    'stroke-dasharray', 'stroke-linecap', 'stroke-linejoin', 'opacity', 'font-family', 'font-size',
    'font-weight', 'font-style', 'text-anchor', 'dominant-baseline', 'letter-spacing',
    'paint-order'];

  /* --- monochrome ---
     The mono theme takes the colour out of everything the stylesheet draws,
     but a feature that came in from SnapGene or ApE carries its own hex and
     never touches the stylesheet. Those are converted here, by luminance,
     and then snapped to the same spread of tones the theme uses, so two
     imported features do not land on the same grey. Anything already grey
     is left exactly as it is, which keeps the theme's own tones, the white
     labels and the black text untouched. */
  const MONO_STEPS = [0x1c, 0x2f, 0x42, 0x55, 0x68, 0x7b, 0x8e];
  function toMono(css) {
    if (!css) return css;
    const m = /^rgba?\(\s*([\d.]+)[,\s]+([\d.]+)[,\s]+([\d.]+)\s*(?:[,/]\s*([\d.]+)\s*)?\)$/i.exec(String(css).trim());
    if (!m) return css;
    const r = +m[1], g = +m[2], b = +m[3], a = m[4] === undefined ? 1 : +m[4];
    if (Math.max(r, g, b) - Math.min(r, g, b) <= 8) return css;      // already grey
    const y = 0.2126 * r + 0.7152 * g + 0.0722 * b;
    let best = MONO_STEPS[0];
    for (const s of MONO_STEPS) if (Math.abs(s - y) < Math.abs(best - y)) best = s;
    return a < 1 ? 'rgba(' + best + ',' + best + ',' + best + ',' + a + ')'
                 : 'rgb(' + best + ',' + best + ',' + best + ')';
  }
  PB.toMono = toMono;
  /* the same conversion for a hex written as a presentation attribute */
  function monoHex(v) {
    const s = String(v).trim();
    let r, g, b;
    let m = /^#([0-9a-f]{3})$/i.exec(s);
    if (m) [r, g, b] = m[1].split('').map(c => parseInt(c + c, 16));
    else if ((m = /^#([0-9a-f]{6})(?:[0-9a-f]{2})?$/i.exec(s)))
      [r, g, b] = [0, 2, 4].map(i => parseInt(m[1].substr(i, 2), 16));
    else {
      const out = toMono(s);
      return out === s ? s : out;
    }
    const out = toMono('rgb(' + r + ', ' + g + ', ' + b + ')');
    const q = /^rgb\((\d+)/.exec(out);
    if (!q) return s;
    const h = (+q[1]).toString(16).padStart(2, '0');
    return '#' + h + h + h;
  }
  PB.monoHex = monoHex;
  PB.mapSVGString = function (opts) {
    opts = opts || {};
    const rec = PB.displayRecord(opts.rec || PB.active());
    if (!rec || !rec.seq.length) throw new Error('No construct with a sequence to draw.');
    const mode = opts.mode || state.mapMode || (rec.topology === 'circular' ? 'circular' : 'linear');
    const cur = PB.mapView(rec, mode);
    const v = opts.whole === false ? cur
      : (mode === 'linear' ? { mode: 'linear', start: 0, span: Math.max(rec.seq.length, 1) }
                           : { mode: 'circular', k: 1, x: 0, y: 0 });
    /* "as zoomed now" has to be drawn in the same frame it is being viewed
       in: the pan offsets are in canvas coordinates, so rendering them at a
       different width frames somewhere else entirely */
    const asShown = opts.whole === false;
    const W = opts.width || (asShown ? canvasWidth() : (mode === 'linear' ? 1600 : 1100));
    const root = document.documentElement;
    const prevTheme = root.getAttribute('data-theme');
    const keepSel = state.selection, keepFeat = state.selectedFeature, keepFeats = state.selectedFeatures;
    let svgText;
    const host = document.createElement('div');
    /* drawn off screen but NOT hidden: a hidden ancestor would be inherited
       by every element and copied into the file, making it all invisible */
    host.style.cssText = 'position:absolute;left:-99999px;top:0;width:' + W + 'px;pointer-events:none';
    host.setAttribute('aria-hidden', 'true');
    document.body.appendChild(host);
    try {
      if (opts.theme !== 'screen') root.setAttribute('data-theme',
        opts.theme === 'dark' ? 'dark' : (opts.theme === 'mono' ? 'mono' : 'light'));
      if (!opts.withSelection) { state.selection = null; state.selectedFeature = null; state.selectedFeatures = []; }
      host.innerHTML = mode === 'linear' ? linearMap(rec, Object.assign({}, v), W) : circularMap(rec, Object.assign({}, v), W);
      const svg = host.querySelector('svg');
      if (!svg) throw new Error('The map did not draw.');
      /* a note telling the reader to pan or zoom belongs on screen, not in
         a picture they are going to put in a slide */
      Array.from(svg.querySelectorAll('.maplegend')).forEach(e => e.remove());
      const mono = opts.theme === 'mono';
      const live = [svg].concat(Array.from(svg.querySelectorAll('*')));
      const styles = live.map(el => {
        const cs = getComputedStyle(el), o = {};
        STYLE_PROPS.forEach(p => { const val = cs.getPropertyValue(p); if (val) o[p] = val; });
        if (mono) { if (o.fill) o.fill = toMono(o.fill); if (o.stroke) o.stroke = toMono(o.stroke); }
        o.display = cs.display;
        return o;
      });
      const bgCol = getComputedStyle(document.body).backgroundColor;
      /* the drawing's own extent, labels included, so the file is not mostly
         empty margin; taken before the copy, while it is laid out */
      let crop = null;
      if (asShown) {
        /* what is on screen, and nothing else. Cropping to the drawing's
           own extent would take in the whole circle, which at a hundred
           times magnification is thirty thousand pixels across and almost
           entirely empty. */
        const vb = (svg.getAttribute('viewBox') || '').split(/\s+/).map(Number);
        if (vb.length === 4 && vb[2] > 0) crop = vb;
      } else if (opts.crop !== false) {
        try {
          const bb = svg.getBBox(), pad = 18;
          if (bb.width > 0 && bb.height > 0)
            crop = [Math.floor(bb.x - pad), Math.floor(bb.y - pad), Math.ceil(bb.width + 2 * pad), Math.ceil(bb.height + 2 * pad)];
        } catch (e) { crop = null; }
      }
      const out = svg.cloneNode(true);
      const copies = [out].concat(Array.from(out.querySelectorAll('*')));
      copies.forEach((el, i) => {
        const st = styles[i];
        if (st.display === 'none') { el.remove(); return; }
        const css = STYLE_PROPS.filter(p => st[p] !== undefined && st[p] !== '' &&
          !(p === 'paint-order' && st[p] === 'normal') && !(p === 'letter-spacing' && st[p] === 'normal'))
          .map(p => p + ':' + st[p]).join(';');
        Array.from(el.attributes).forEach(a => {
          if (a.name === 'class' || a.name === 'style' || a.name === 'role' || a.name === 'tabindex' ||
              a.name.indexOf('data-') === 0 || a.name.indexOf('on') === 0 || a.name === 'cursor') el.removeAttribute(a.name);
          else if (/var\(/.test(a.value)) el.removeAttribute(a.name);   /* resolved into the style below */
          /* A colour carried in from SnapGene or ApE is written as a
             presentation attribute. The style set below already overrides it
             when the file is drawn, but leaving the hex in the file would
             put the colour back the moment anyone stripped the styles —
             an editor, an Illustrator import, a journal's converter. So it
             is rewritten, not merely overridden. */
          else if (mono && (a.name === 'fill' || a.name === 'stroke') && /^#|^rgb/i.test(a.value))
            el.setAttribute(a.name, monoHex(a.value));
        });
        if (css) el.setAttribute('style', css);
      });
      /* a fixed size, so the image is not 0 x 0 outside the page */
      if (crop) out.setAttribute('viewBox', crop.join(' '));
      const vb = (out.getAttribute('viewBox') || '0 0 ' + W + ' ' + W).split(/\s+/).map(Number);
      out.setAttribute('width', String(vb[2])); out.setAttribute('height', String(vb[3]));
      out.setAttribute('xmlns', 'http://www.w3.org/2000/svg');
      out.setAttribute('version', '1.1');
      out.style.overflow = 'visible';
      if (opts.background !== 'transparent') {
        const bg = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
        bg.setAttribute('x', String(vb[0])); bg.setAttribute('y', String(vb[1]));
        bg.setAttribute('width', String(vb[2])); bg.setAttribute('height', String(vb[3]));
        bg.setAttribute('fill', mono ? '#ffffff'
          : (opts.background === 'page' ? bgCol : (opts.theme === 'dark' ? '#111618' : '#ffffff')));
        out.insertBefore(bg, out.firstChild);
      }
      const title = document.createElementNS('http://www.w3.org/2000/svg', 'title');
      title.textContent = rec.name + ' — ' + rec.seq.length + ' bp ' + rec.topology;
      out.insertBefore(title, out.firstChild);
      svgText = '<?xml version="1.0" encoding="UTF-8"?>\n' + new XMLSerializer().serializeToString(out);
    } finally {
      host.remove();
      if (opts.theme !== 'screen') {
        if (prevTheme === null) root.removeAttribute('data-theme'); else root.setAttribute('data-theme', prevTheme);
      }
      state.selection = keepSel; state.selectedFeature = keepFeat; state.selectedFeatures = keepFeats;
    }
    return svgText;
  };

  /* PNG drawn from the standalone SVG at the chosen scale */
  PB.mapPNGBlob = function (svgText, scale) {
    return new Promise((resolve, reject) => {
      const m = svgText.match(/<svg[^>]*\swidth="([\d.]+)"[^>]*\sheight="([\d.]+)"/) ||
                svgText.match(/<svg[^>]*\sheight="([\d.]+)"[^>]*\swidth="([\d.]+)"/);
      const w = m ? parseFloat(m[1]) : 1000, h = m ? parseFloat(m[2]) : 1000;
      const img = new Image();
      img.onload = () => {
        try {
          const c = document.createElement('canvas');
          c.width = Math.round(w * scale); c.height = Math.round(h * scale);
          const ctx = c.getContext('2d');
          ctx.scale(scale, scale);
          ctx.drawImage(img, 0, 0, w, h);
          c.toBlob(b => b ? resolve(b) : reject(new Error('The browser could not encode the PNG.')), 'image/png');
        } catch (e) { reject(e); }
      };
      img.onerror = () => reject(new Error('The map image could not be rendered.'));
      img.src = 'data:image/svg+xml;charset=utf-8,' + encodeURIComponent(svgText);
    });
  };

  PB.exportMapDialog = function () {
    const rec = PB.active();
    if (!rec || !rec.seq.length) { PB.toast('Load a construct first.', 'warn'); return; }
    const safe = rec.name.replace(/[^\w.-]+/g, '_').slice(0, 60) || 'map';
    const mode0 = state.mapMode || (rec.topology === 'circular' ? 'circular' : 'linear');
    const sel = (id, opts, val) => '<select class="inp" id="' + id + '">' + opts.map(([v, l]) =>
      '<option value="' + v + '"' + (v === val ? ' selected' : '') + '>' + l + '</option>').join('') + '</select>';
    PB.modal({
      title: 'Export map image',
      body: '<div class="row">' +
        '<div class="field narrow"><label>Format</label>' + sel('meFmt', [['svg', 'SVG (vector)'], ['png', 'PNG']], 'svg') + '</div>' +
        '<div class="field narrow"><label>Map</label>' + sel('meMode', [['circular', 'Circular'], ['linear', 'Linear']], mode0) + '</div>' +
        '<div class="field narrow"><label>Extent</label>' + sel('meWhole', [['whole', 'Whole construct'], ['view', 'As zoomed now']], 'whole') + '</div>' +
        '</div><div class="row">' +
        '<div class="field narrow"><label>Colours</label>' + sel('meTheme', [['light', 'Light (papers, slides)'], ['mono', 'Monochrome (print)'], ['dark', 'Dark'], ['screen', 'As on screen']], 'light') + '</div>' +
        '<div class="field narrow"><label>Background</label>' + sel('meBg', [['solid', 'Solid'], ['transparent', 'Transparent']], 'solid') + '</div>' +
        '<div class="field narrow" id="meScaleF"><label>PNG resolution</label>' + sel('meScale', [['2', '2× (screen)'], ['3', '3×'], ['4', '4× (print)']], '3') + '</div>' +
        '</div>' +
        '<div class="row"><label class="chk"><input type="checkbox" id="meSel"> Show the current selection</label></div>' +
        '<div class="hint">The image shows exactly what the map shows — features, labels and the enzymes ticked in the rail. ' +
        'SVG stays sharp at any size and every label stays editable text. Fonts are named rather than embedded, so on a ' +
        'computer without them a close substitute is used.<br>' +
        '<b>Monochrome</b> replaces every colour with a tone, including colours carried in from SnapGene or ApE, ' +
        'so the figure survives a greyscale printer or a journal that charges for colour. Each tone stays dark enough ' +
        'for the white label drawn along a feature to read on it.<br>' +
        'Which file types can be written is decided by whatever is showing this page, and SVG is not always one of them. ' +
        'If it is refused you are offered the same drawing as a PNG, or the SVG under a .txt name to rename — nothing is ' +
        'substituted without being asked for.</div>' +
        '<div class="map-preview" id="mePrev" style="margin-top:10px;border:1px solid var(--border);border-radius:8px;overflow:hidden;max-height:340px;display:flex;justify-content:center;background:var(--surface-2)"></div>',
      actions: [{ label: 'Cancel' }, { label: 'Save image', primary: true, onClick: () => { go(); return false; } }],
      onOpen: w => {
        const upd = () => {
          w.querySelector('#meScaleF').style.display = $('#meFmt').value === 'png' ? '' : 'none';
          try {
            const svg = build();
            w.querySelector('#mePrev').innerHTML = '<img alt="Preview of the exported map" style="max-width:100%;max-height:340px;display:block" src="data:image/svg+xml;charset=utf-8,' +
              encodeURIComponent(svg) + '">';
          } catch (e) { w.querySelector('#mePrev').textContent = e.message; }
        };
        ['#meFmt', '#meMode', '#meWhole', '#meTheme', '#meBg', '#meScale', '#meSel'].forEach(id =>
          w.querySelector(id).addEventListener('change', upd));
        upd();
      }
    });
    function build() {
      return PB.mapSVGString({ rec: rec, mode: $('#meMode').value, whole: $('#meWhole').value === 'whole',
        theme: $('#meTheme').value, background: $('#meBg').value === 'transparent' ? 'transparent' : 'solid',
        withSelection: $('#meSel').checked });
    }
    async function go() {
      let svg;
      try { svg = build(); } catch (e) { PB.toast(e.message, 'bad'); return; }
      const name = safe + '_' + $('#meMode').value + '_map';
      const scale = parseFloat($('#meScale').value) || 3;
      const asPng = async () => {
        let blob;
        try { blob = await PB.mapPNGBlob(svg, scale); }
        catch (e) { PB.toast(e.message, 'bad'); return; }
        await PB.offerBlob(name + '.png', blob, 'Map (PNG)');
      };
      if ($('#meFmt').value === 'svg') {
        PB.closeModal();
        /* If this view will not write a .svg, the picture is not lost: the
           same drawing goes out as a PNG, which every host takes, and the
           vector text is still there to save under another name. Both are
           offered; neither happens on its own. */
        await PB.offerFile(name + '.svg', svg, 'Map (SVG)',
          { alt: { label: 'Save it as a PNG instead', run: asPng } });
        return;
      }
      PB.closeModal();
      await asPng();
    }
  };

})(window.PB);

