/* ============================================================================
   PlasmidBench — sequencing from the freezer

   Two things live here: the primer stock, which is the inventory a lab has
   accumulated and which is far too long to show on a map, and the tool that
   searches it for the tubes that can read a chosen region.

   The stock is kept apart from the working primer library on purpose. A
   library of six primers is something you look at; an inventory of nine
   hundred is something you search.
   ============================================================================ */
(function (PB) {
  'use strict';
  const B = window.Bio;
  const $ = PB.$, esc = PB.esc, fmt = PB.fmt, f1 = PB.f1;
  const state = PB.state;

  /* ---------------------------------------------------------------- stock */
  state.stock = PB.lsGet('stock', []) || [];
  state.stockName = PB.lsGet('stockName', '') || '';
  function saveStock() { PB.lsSet('stock', state.stock); PB.lsSet('stockName', state.stockName); }
  PB.stockCount = () => (state.stock || []).length;

  /* every oligo the tool can search: the inventory, plus the working
     primers, which are often the ones just designed and not yet ordered */
  function allOligos(which) {
    const stock = (state.stock || []).map(p => Object.assign({ source: 'stock' }, p));
    const lib = (state.primers || []).map(p => Object.assign({ source: 'library' }, p));
    if (which === 'stock') return stock;
    if (which === 'library') return lib;
    /* an oligo in both places is one oligo */
    const seen = {};
    return stock.concat(lib).filter(p => {
      const k = (p.seq || '').toUpperCase();
      if (seen[k]) return false; seen[k] = 1; return true;
    });
  }
  PB.allOligos = allOligos;

  /* =======================================================================
     IMPORTING A STOCK LIST
  ======================================================================= */
  PB.stockImportDialog = function () {
    const n = PB.stockCount();
    PB.modal({
      title: 'Primer stock',
      sub: n ? n + ' oligos' + (state.stockName ? ' from ' + state.stockName : '') : 'nothing imported yet',
      wide: true,
      body:
        '<div class="note">Paste the inventory as it comes: a spreadsheet selection, a CSV, or text ' +
        'copied out of a printed list. Columns are recognised by what they hold rather than by their ' +
        'order, sequences typed in groups (<span class="mono">TTG TAT TTC</span>) are closed up, and a ' +
        'record that wrapped onto a second line is still read as one oligo. Nothing is sent anywhere — ' +
        'it stays in this browser.</div>' +
        '<div class="row"><div class="field"><label>Where this came from (optional)</label>' +
        '<input class="inp" id="skName" placeholder="Primer_Stocks_-20" value="' + esc(state.stockName) + '"></div></div>' +
        '<div class="field"><label>Paste the list</label>' +
        '<textarea class="inp mono" id="skText" rows="12" spellcheck="false" ' +
        'placeholder="OL-001  L3  WTOL001  A01  2019-09-25 KPRB:S425F  TTG TAT TTC CAG GGC ATC GAG  62"></textarea></div>' +
        '<div id="skPrev"></div>',
      actions: [
        { label: 'Close' },
        n ? { label: 'Clear the stock', onClick: () => {
          state.stock = []; state.stockName = ''; saveStock();
          PB.toast('Stock cleared', 'ok'); PB.stockImportDialog(); return false;
        } } : null,
        { label: 'Check it', onClick: () => { previewStock(); return false; }, close: false },
        { label: 'Import', primary: true, onClick: () => {
          const r = parsePasted();
          if (!r) return false;
          if (!r.primers.length) { PB.toast('Nothing in that paste looked like an oligo.', 'bad'); return false; }
          /* replacing wholesale is what people mean by importing their
             stock again, but the ones already there are kept if the
             sequence is new */
          const seen = {};
          const merged = r.primers.slice();
          merged.forEach(p => { seen[p.seq] = 1; });
          (state.stock || []).forEach(p => { if (!seen[p.seq]) merged.push(p); });
          state.stock = merged;
          state.stockName = ($('#skName') || {}).value || '';
          saveStock();
          PB.toast(r.primers.length + ' oligos imported; the stock now holds ' + merged.length, 'ok');
          return true;
        } }
      ].filter(Boolean)
    });
    const ta = $('#skText');
    if (ta) ta.addEventListener('input', () => { const p = $('#skPrev'); if (p) p.innerHTML = ''; });
  };

  function parsePasted() {
    const ta = $('#skText');
    if (!ta || !ta.value.trim()) { PB.toast('Paste the list first.', 'warn'); return null; }
    try { return B.parseStockTable(ta.value); }
    catch (e) { PB.toast('That list could not be read: ' + e.message, 'bad'); return null; }
  }

  function previewStock() {
    const r = parsePasted();
    const host = $('#skPrev');
    if (!r || !host) return;
    const rows = r.primers.slice(0, 12).map(p =>
      '<tr><td class="mono">' + esc(p.id || '—') + '</td><td>' + esc(p.name || '') + '</td>' +
      '<td class="mono">' + esc(p.loc || '—') + '</td>' +
      '<td class="num">' + (p.tm == null ? '—' : p.tm) + '</td>' +
      '<td class="num">' + p.seq.length + '</td>' +
      '<td class="mono">' + esc(p.seq.length > 30 ? p.seq.slice(0, 29) + '…' : p.seq) + '</td></tr>').join('');
    host.innerHTML =
      '<div class="note ' + (r.primers.length ? 'ok' : 'bad') + '">Read <strong>' + r.primers.length +
      '</strong> oligos' + (r.mode === 'header' ? ' using the column headings' : ' by the shape of each row') +
      (r.skipped.length ? ', and skipped ' + r.skipped.length + ' row' + (r.skipped.length === 1 ? '' : 's') +
        ' with no sequence in them' : '') + '.</div>' +
      (r.skipped.length ? '<div class="hint">Skipped: ' +
        esc(r.skipped.slice(0, 4).map(s => (s.line || '').trim().slice(0, 48)).join(' · ')) +
        (r.skipped.length > 4 ? ' …' : '') + '</div>' : '') +
      '<div class="scrollx"><table class="tbl sm"><thead><tr><th>ID</th><th>Name</th><th>Where</th>' +
      '<th>Tm</th><th>nt</th><th>Sequence</th></tr></thead><tbody>' + rows + '</tbody></table></div>' +
      (r.primers.length > 12 ? '<div class="hint">… and ' + (r.primers.length - 12) + ' more.</div>' : '');
  }

  /* =======================================================================
     THE SEQUENCING TOOL
  ======================================================================= */
  let last = null;        /* the plan on screen, for export and gap filling */

  PB.toolSanger = function () {
    const rec = PB.active();
    if (!rec) { PB.toast('Open a construct first.', 'warn'); return; }
    const sel = state.selection;
    const feats = (rec.features || []).filter(f => f.start != null && f.end != null);
    /* The published defaults, overridden by whatever the lab profile says
       this facility actually gives. A plan drawn for 700 confident bases
       when the facility returns 500 sends someone to the bench with a gap
       they will only find after paying for the reads. */
    const S = Object.assign({}, B.SANGER, {
      dead: PB.lsGet('sanger.dead', B.SANGER.dead),
      good: PB.lsGet('sanger.good', B.SANGER.good),
      max: Math.max(PB.lsGet('sanger.good', B.SANGER.good) + 100, B.SANGER.max)
    });
    const n = PB.stockCount();

    PB.modal({
      title: 'Sanger sequencing from stock',
      sub: rec.name,
      wide: true,
      body:
        (n ? '' : '<div class="note warn">No primer stock has been imported yet, so only the ' +
          (state.primers || []).length + ' primers in the working library will be searched. ' +
          'Import your inventory from <strong>Primer stock</strong> to search the freezer.</div>') +
        '<div class="row">' +
          '<div class="field"><label>Region to read</label><select class="inp" id="sqWhat">' +
            (sel ? '<option value="sel">the selection (' + fmt(Math.abs(sel.end - sel.start)) + ' bp)</option>' : '') +
            '<option value="feat">a feature…</option>' +
            '<option value="manual">bases I type</option>' +
            '<option value="all">the whole construct</option>' +
          '</select></div>' +
          '<div class="field narrow" id="sqFromF"><label>From</label><input class="inp mono" id="sqFrom" type="number" min="1" value="' +
            ((sel ? Math.min(sel.start, sel.end) : 0) + 1) + '"></div>' +
          '<div class="field narrow" id="sqToF"><label>To</label><input class="inp mono" id="sqTo" type="number" min="1" value="' +
            (sel ? Math.max(sel.start, sel.end) : Math.min(1000, rec.seq.length)) + '"></div>' +
        '</div>' +
        '<div class="row" id="sqFeatRow" style="display:none"><div class="field"><label>Feature</label>' +
          '<select class="inp" id="sqFeat">' + feats.map((f, i) =>
            '<option value="' + i + '">' + esc((f.label || f.type || 'feature')) + ' · ' +
            fmt(f.start + 1) + '–' + fmt(f.end) + '</option>').join('') + '</select></div></div>' +
        '<div class="row">' +
          '<div class="field narrow"><label>Unreadable at the start</label>' +
            '<input class="inp mono" id="sqDead" type="number" min="0" max="200" value="' + S.dead + '"></div>' +
          '<div class="field narrow"><label>Confident read</label>' +
            '<input class="inp mono" id="sqGood" type="number" min="100" max="1200" value="' + S.good + '"></div>' +
          '<div class="field narrow"><label>Furthest reach</label>' +
            '<input class="inp mono" id="sqMax" type="number" min="100" max="1500" value="' + S.max + '"></div>' +
          '<div class="field"><label>Search</label><select class="inp" id="sqSrc">' +
            '<option value="all">stock and working primers (' + allOligos('all').length + ')</option>' +
            '<option value="stock">only the stock (' + n + ')</option>' +
            '<option value="library">only the working primers (' + (state.primers || []).length + ')</option>' +
          '</select></div>' +
        '</div>' +
        '<div class="hint">Defaults follow what sequencing services publish: the first 25–35 bases after a ' +
        'primer are expected to be poor, 700 bases is a safe confident read, and 800–900 is reachable on clean ' +
        'template with anything beyond read cautiously. Set them to whatever your facility actually delivers.</div>' +
        '<div id="sqOut"></div>',
      actions: [
        { label: 'Close' },
        { label: 'Import stock…', onClick: () => { PB.stockImportDialog(); return false; } },
        { label: 'Find primers', primary: true, onClick: () => { runSanger(rec); return false; }, close: false }
      ]
    });

    const what = $('#sqWhat');
    const sync = () => {
      const v = what.value;
      $('#sqFeatRow').style.display = v === 'feat' ? '' : 'none';
      const manual = v === 'manual';
      $('#sqFromF').style.display = manual ? '' : 'none';
      $('#sqToF').style.display = manual ? '' : 'none';
    };
    what.addEventListener('change', sync); sync();
    if (!sel) what.value = feats.length ? 'feat' : 'manual', sync();
  };

  function regionFrom(rec) {
    const v = $('#sqWhat').value;
    const L = rec.seq.length;
    if (v === 'sel' && state.selection)
      return { start: Math.min(state.selection.start, state.selection.end),
               end: Math.max(state.selection.start, state.selection.end) };
    if (v === 'feat') {
      const f = (rec.features || []).filter(x => x.start != null && x.end != null)[+$('#sqFeat').value];
      if (f) return { start: f.start, end: f.end };
    }
    if (v === 'all') return { start: 0, end: L };
    const a = Math.max(1, parseInt($('#sqFrom').value, 10) || 1) - 1;
    const b = Math.max(1, parseInt($('#sqTo').value, 10) || L);
    return { start: a, end: b };
  }

  function runSanger(rec) {
    const host = $('#sqOut');
    const region = regionFrom(rec);
    const circ = rec.topology === 'circular';
    const dead = Math.max(0, parseInt($('#sqDead').value, 10) || 0);
    const good = Math.max(50, parseInt($('#sqGood').value, 10) || 700);
    const max = Math.max(good, parseInt($('#sqMax').value, 10) || 900);
    const oligos = allOligos($('#sqSrc').value);
    if (!oligos.length) { host.innerHTML = '<div class="note bad">There are no primers to search.</div>'; return; }
    host.innerHTML = '<div class="note">Searching ' + oligos.length + ' oligos…</div>';
    /* let the note paint before the work starts */
    setTimeout(() => {
      let plan;
      try {
        plan = B.sangerPlan(rec.seq, region, oligos,
          { circular: circ, dead: dead, good: good, max: max });
      } catch (e) {
        host.innerHTML = '<div class="note bad">That search failed: ' + esc(e.message) + '</div>';
        return;
      }
      /* pairs that would amplify the region, for sequencing a PCR product
         rather than the plasmid itself */
      let amps = [];
      try {
        amps = B.sangerAmplicons(rec.seq, region, oligos,
          { circular: circ, dead: dead }).slice(0, 8);
      } catch (e) { amps = []; }
      last = { plan: plan, rec: rec, region: region, circ: circ,
               dead: dead, good: good, max: max, oligos: oligos, amps: amps };
      host.innerHTML = renderPlan(plan, rec, circ) + renderAmps(amps, rec);
      wirePlan(rec);
    }, 10);
  }

  /* -------------------------------------------------------------- drawing */
  function coverageBar(plan, rec) {
    const W = 760, H = 104, pad = 8;
    const len = plan.region.length;
    const L = rec.seq.length;
    const x = off => pad + (off / len) * (W - 2 * pad);
    const rel = p => ((p - plan.region.start) % L + L) % L;
    const lane = (cands, y, colour) => cands.map(c => {
      const w = c.window.good;
      let from = rel(w.start), to = from + w.length;
      /* a read that starts before the region, or runs past it, is clipped
         to the picture rather than drawn off the end of it */
      if (from > len) { to -= L; from -= L; }
      const a = Math.max(0, from), b = Math.min(len, to);
      if (b <= a) return '';
      return '<rect x="' + x(a).toFixed(1) + '" y="' + y + '" width="' + Math.max(1.5, x(b) - x(a)).toFixed(1) +
        '" height="12" rx="3" fill="' + colour + '" opacity="0.85"></rect>' +
        /* the direction it reads */
        '<path d="M' + (c.strand === 1 ? x(b).toFixed(1) : x(a).toFixed(1)) + ' ' + (y + 6) +
        ' l' + (c.strand === 1 ? '-5 -4 0 8' : '5 -4 0 8') + 'z" fill="' + colour + '"></path>';
    }).join('');
    const gapMarks = (gaps, y, colour) => gaps.map(g =>
      '<rect x="' + x(g.from).toFixed(1) + '" y="' + y + '" width="' + Math.max(1.5, x(g.to) - x(g.from)).toFixed(1) +
      '" height="12" rx="3" fill="none" stroke="' + colour + '" stroke-dasharray="3 2"></rect>').join('');
    return '<svg viewBox="0 0 ' + W + ' ' + H + '" class="covbar" style="width:100%;height:auto;display:block">' +
      '<rect x="' + pad + '" y="50" width="' + (W - 2 * pad) + '" height="14" rx="3" fill="var(--panel2, #e5e7eb)"></rect>' +
      '<text x="' + pad + '" y="13" class="cbl">forward reads →</text>' +
      lane(plan.plan.forward, 20, '#2563eb') + gapMarks(plan.gaps.forward, 20, '#dc2626') +
      '<text x="' + pad + '" y="46" class="cbl">the region, ' + fmt(plan.region.length) + ' bp</text>' +
      '<text x="' + pad + '" y="99" class="cbl">← reverse reads</text>' +
      lane(plan.plan.reverse, 72, '#059669') + gapMarks(plan.gaps.reverse, 72, '#dc2626') +
      '</svg>';
  }

  function pill(c, i) {
    if (c.warnings.length === 0) return '<span class="badge ok">clean</span>';
    const bad = c.warnings.some(w => /mixed|unreadable|slip/.test(w));
    /* a count of notes is no use on its own: the notes are what decides
       whether you trust the read, so the badge opens them */
    return '<button class="badge ' + (bad ? 'bad' : 'warn') + '" data-notes="' + i + '" ' +
      'title="' + esc(c.warnings.join(' — ')) + '">' + c.warnings.length +
      ' note' + (c.warnings.length === 1 ? '' : 's') + '</button>';
  }

  function primerRow(c, rec, i, picked) {
    const L = rec.seq.length;
    const site = c.site.wraps ? (c.site.start + 1) + '–' + L + ' + 1–' + c.site.end
                              : (c.site.start + 1) + '–' + c.site.end;
    return '<tr data-cand="' + i + '">' +
      '<td><input type="checkbox" class="sqpick" data-pick="' + i + '"' + (picked ? ' checked' : '') +
        ' title="Include this primer when the button below writes them onto the construct"></td>' +
      '<td class="mono">' + esc(c.id || '—') + '</td>' +
      '<td>' + esc(c.name || '') + (c.primer.source === 'library' ? ' <span class="badge neutral">library</span>' : '') + '</td>' +
      '<td class="mono">' + esc(c.loc || '—') + '</td>' +
      '<td>' + (c.strand === 1 ? 'forward →' : '← reverse') + '</td>' +
      '<td class="mono">' + site + '</td>' +
      '<td class="num">' + f1(c.tm) + '</td>' +
      '<td class="num">' + fmt(c.covers) + '</td>' +
      '<td class="num">' + (c.inside ? 'inside' : fmt(c.gap)) + '</td>' +
      '<td>' + pill(c, i) + '</td>' +
      '<td><button class="btn sm ghost" data-show="' + i + '">show</button>' +
      '<button class="btn sm ghost" data-annot="' + i + '" title="Add this primer to the construct as a labelled feature">annotate</button></td></tr>' +
      (c.warnings.length
        ? '<tr class="noterow" id="nr' + i + '" hidden><td colspan="11"><ul class="notelist">' +
          c.warnings.map(w => '<li>' + esc(w) + '</li>').join('') + '</ul></td></tr>'
        : '');
  }

  function renderPlan(plan, rec, circ) {
    const cov = plan.coverage;
    const pc = v => Math.round((v / cov.length) * 100) + '%';
    const complete = cov.both === cov.length;
    const halfway = cov.either === cov.length;
    const head =
      '<div class="note ' + (complete ? 'ok' : halfway ? 'warn' : 'bad') + '">' +
      (complete
        ? 'Your stock covers the whole region on <strong>both strands</strong> — ' + fmt(cov.length) + ' bp, confirmed each way.'
        : halfway
          ? 'Your stock reads every base of the region, but only ' + pc(cov.both) + ' of it on both strands. ' +
            'Single-strand sequence is enough to see what is there; both strands is what makes it confirmed.'
          : 'Your stock reads ' + pc(cov.either) + ' of the region. ' + fmt(cov.length - cov.either) +
            ' bp cannot be reached by any primer you have.') +
      '</div>' +
      '<div class="chipbar">' +
        '<span class="chip">' + fmt(cov.length) + ' bp region</span>' +
        '<span class="chip">forward ' + pc(cov.forward) + '</span>' +
        '<span class="chip">reverse ' + pc(cov.reverse) + '</span>' +
        '<span class="chip">both strands ' + pc(cov.both) + '</span>' +
        '<span class="chip">' + plan.candidates.length + ' usable primer' + (plan.candidates.length === 1 ? '' : 's') + '</span>' +
        (plan.unusable.length ? '<span class="chip bad">' + plan.unusable.length + ' refused</span>' : '') +
      '</div>' +
      coverageBar(plan, rec);

    const setTable = (which, label) => {
      const rows = plan.plan[which];
      if (!rows.length) return '<div class="note warn">Nothing in your stock reads this region ' +
        (which === 'forward' ? 'forward' : 'in reverse') + '.</div>';
      return '<h4>' + label + '</h4><div class="scrollx"><table class="tbl sm"><thead><tr>' +
        '<th title="Include when annotating">\u2713</th>' +
        '<th>ID</th><th>Name</th><th>Where it is</th><th>Reads</th><th>Anneals at</th><th>Tm</th>' +
        '<th>Covers</th><th>Run-up</th><th></th><th></th></tr></thead><tbody>' +
        rows.map(c => primerRow(c, rec, plan.candidates.indexOf(c), true)).join('') +
        '</tbody></table></div>';
    };

    /* the reactions to actually set up */
    const reactions = plan.plan.forward.concat(plan.plan.reverse);
    const order = reactions.length
      ? '<h4>The reactions this comes to</h4><div class="note">' + reactions.length +
        ' sequencing reaction' + (reactions.length === 1 ? '' : 's') + ' on <strong>' + esc(rec.name) +
        '</strong>, pulling ' + reactions.length + ' tube' + (reactions.length === 1 ? '' : 's') + ': ' +
        reactions.map(c => '<span class="mono">' + esc(c.id || c.name) + '</span>' +
          (c.loc ? ' (' + esc(c.loc) + ')' : '')).join(', ') + '.</div>'
      : '';

    const gapList = ['forward', 'reverse'].map(k => {
      const g = plan.gaps[k];
      if (!g.length) return '';
      return '<div class="note warn"><strong>' + (k === 'forward' ? 'Forward' : 'Reverse') +
        ' strand not covered:</strong> ' + g.map(x =>
          'bases ' + fmt(x.from + 1) + '–' + fmt(x.to) + ' of the region (' +
          fmt(x.start + 1) + '–' + fmt(x.end) + ' on the construct)').join('; ') +
        ' <button class="btn sm" data-fill="' + k + '">Design a primer for this</button></div>';
    }).join('');

    const refused = plan.unusable.length
      ? '<h4>Refused, and why</h4><div class="scrollx"><table class="tbl sm"><thead><tr><th>ID</th><th>Name</th>' +
        '<th>Where it is</th><th>Why not</th></tr></thead><tbody>' +
        plan.unusable.slice(0, 20).map(c => '<tr><td class="mono">' + esc(c.id || '—') + '</td><td>' +
          esc(c.name || '') + '</td><td class="mono">' + esc(c.loc || '—') + '</td><td>' +
          esc(c.warnings[0] || '') + '</td></tr>').join('') + '</tbody></table></div>'
      : '';

    const others = plan.candidates.filter(c => plan.plan.forward.indexOf(c) < 0 && plan.plan.reverse.indexOf(c) < 0);
    const rest = others.length
      ? '<h4>Other primers that would also work</h4><div class="scrollx"><table class="tbl sm"><thead><tr>' +
        '<th title="Include when annotating">\u2713</th>' +
        '<th>ID</th><th>Name</th><th>Where it is</th><th>Reads</th><th>Anneals at</th><th>Tm</th>' +
        '<th>Covers</th><th>Run-up</th><th></th><th></th></tr></thead><tbody>' +
        others.slice(0, 40).map(c => primerRow(c, rec, plan.candidates.indexOf(c), false)).join('') +
        '</tbody></table></div>' +
        (others.length > 40 ? '<div class="hint">… and ' + (others.length - 40) + ' more.</div>' : '')
      : '';

    return head + order + setTable('forward', 'Forward reads to set up') +
      setTable('reverse', 'Reverse reads to set up') + gapList +
      '<div class="row" style="margin-top:10px">' +
      '<button class="btn" id="sqXls">Export the plan as a spreadsheet</button>' +
      '<button class="btn ghost" id="sqCopy">Copy as text</button>' +
      '<button class="btn" id="sqAnnot">Add the ticked primers as features</button>' +
      '<button class="btn sm ghost" id="sqPickAll">Tick all</button>' +
      '<button class="btn sm ghost" id="sqPickNone">Tick none</button>' +
      '<button class="btn ghost" id="sqMark">Draw these on the map</button></div>' +
      '<div class="hint">Every primer in the tables above has a tick box, and only the ticked ones are written ' +
      'onto the construct. The ones the plan chose start ticked; the rest do not. Annotating writes them onto ' +
      'the construct itself, so they show on the map, sit in the Features tab and travel with the file when ' +
      'you export it. Drawing them on the map instead adds them to your primer library, which is kept in this ' +
      'browser and not in the construct.</div>' +
      rest + refused +
      '<div class="hint">A read is counted only where the confident window reaches, so a primer whose ' +
      'window ends before the region is listed as within reach but not as covering it.</div>';
  }

  /* pairs in the stock that would amplify the region instead */
  function renderAmps(amps, rec) {
    if (!amps.length) return '';
    return '<h4>Or amplify it first</h4>' +
      '<div class="note">These pairs are already in your stock and would give a product with the region ' +
      'well inside it, for sequencing an amplicon rather than the plasmid. Both primers of a pair anneal ' +
      'in exactly one place.</div>' +
      '<div class="scrollx"><table class="tbl sm"><thead><tr><th>Forward</th><th>Where it is</th>' +
      '<th>Reverse</th><th>Where it is</th><th>Product</th><th>Tm f / r</th><th>Clear of the ends</th>' +
      '<th>Notes</th></tr></thead><tbody>' +
      amps.map(a => '<tr>' +
        '<td class="mono">' + esc(a.forward.id || a.forward.name || '') + '</td>' +
        '<td class="mono">' + esc(a.forward.loc || '\u2014') + '</td>' +
        '<td class="mono">' + esc(a.reverse.id || a.reverse.name || '') + '</td>' +
        '<td class="mono">' + esc(a.reverse.loc || '\u2014') + '</td>' +
        '<td class="num">' + fmt(a.length) + ' bp</td>' +
        '<td class="num">' + f1(a.forward.tm) + ' / ' + f1(a.reverse.tm) + '</td>' +
        '<td class="num">' + fmt(a.inset.fromStart) + ' / ' + fmt(a.inset.fromEnd) + '</td>' +
        '<td>' + (a.warnings.length
          ? '<span class="badge warn" title="' + esc(a.warnings.join(' \u2014 ')) + '">' +
            a.warnings.length + '</span>'
          : '<span class="badge ok">clean</span>') + '</td></tr>').join('') +
      '</tbody></table></div>';
  }

  /* -------------------------------------------------------------- wiring */
  function wirePlan(rec) {
    const plan = last.plan;
    (PB.$$ || (s => Array.from(document.querySelectorAll(s))))('#sqOut [data-show]').forEach(b =>
      b.addEventListener('click', () => {
        const c = plan.candidates[+b.getAttribute('data-show')];
        if (!c) return;
        state.selection = { start: c.site.start, end: c.site.end };
        PB.renderAll();
        PB.toast((c.id ? c.id + ' · ' : '') + (c.name || '') + ' anneals ' +
          (c.strand === 1 ? 'forward' : 'reverse') + ' at ' + (c.site.start + 1) + '–' + c.site.end +
          (c.warnings.length ? ' — ' + c.warnings[0] : ''), c.warnings.length ? 'warn' : 'ok');
      }));
    Array.from(document.querySelectorAll('#sqOut [data-fill]')).forEach(b =>
      b.addEventListener('click', () => fillGap(b.getAttribute('data-fill'))));
    Array.from(document.querySelectorAll('#sqOut [data-notes]')).forEach(b =>
      b.addEventListener('click', () => {
        const r = document.getElementById('nr' + b.getAttribute('data-notes'));
        if (r) r.hidden = !r.hidden;
      }));
    Array.from(document.querySelectorAll('#sqOut [data-annot]')).forEach(b =>
      b.addEventListener('click', () => {
        const c = plan.candidates[+b.getAttribute('data-annot')];
        if (c) annotatePrimers([c]);
      }));
    const picks = () => Array.from(document.querySelectorAll('#sqOut .sqpick'));
    const chosen = () => picks().filter(b => b.checked)
      .map(b => plan.candidates[+b.getAttribute('data-pick')]).filter(Boolean);
    const an = $('#sqAnnot');
    const relabel = () => {
      if (!an) return;
      const n = picks().filter(b => b.checked).length;
      an.textContent = n ? 'Add the ' + n + ' ticked primer' + (n === 1 ? '' : 's') + ' as features'
                         : 'Add the ticked primers as features';
      an.disabled = !n;
    };
    picks().forEach(b => b.addEventListener('change', relabel));
    const setAll = v => { picks().forEach(b => { b.checked = v; }); relabel(); };
    const pa = $('#sqPickAll'); if (pa) pa.addEventListener('click', () => setAll(true));
    const pn = $('#sqPickNone'); if (pn) pn.addEventListener('click', () => setAll(false));
    relabel();
    if (an) an.addEventListener('click', () => {
      const list = chosen();
      if (!list.length) { PB.toast('Tick the primers you want written onto the construct first.', 'warn'); return; }
      annotatePrimers(list);
    });
    const x = $('#sqXls'); if (x) x.addEventListener('click', sequencingExportPlan);
    const cp = $('#sqCopy'); if (cp) cp.addEventListener('click', () => PB.copyText(sequencingPlanText()));
    const mk = $('#sqMark'); if (mk) mk.addEventListener('click', markOnMap);
  }

  /* Write a primer onto the construct as a real, labelled feature, so it
     is part of the map and travels with the file on export - as opposed
     to the library primers, which are drawn over the map but belong to
     the browser rather than to the construct. */
  function annotatePrimers(cands) {
    const rec = PB.active();
    if (!rec) return;
    if (rec.builtin) {
      PB.toast('The demo construct is read-only — duplicate it first.', 'warn');
      return;
    }
    let added = 0, already = 0;
    cands.forEach(c => {
      const label = ((c.id ? c.id + ' ' : '') + (c.name || 'primer')).trim();
      const dup = (rec.features || []).some(f =>
        f.start === c.site.start && f.end === c.site.end &&
        (f.strand || 1) === c.strand && (f.label || '') === label);
      if (dup) { already++; return; }
      rec.features.push(PB.normFeature({
        type: 'primer_bind', start: c.site.start, end: c.site.end, strand: c.strand,
        label: label,
        quals: {
          note: 'sequencing primer' + (c.loc ? ', stock location ' + c.loc : '') +
                '; anneals ' + c.annealed + ' nt, Tm ' + c.tm.toFixed(1) + ' \u00b0C' +
                (c.mismatches.length ? ', ' + c.mismatches.length + ' mismatch' +
                  (c.mismatches.length === 1 ? '' : 'es') : ''),
          PCR_primer_seq: (c.seq || '').toUpperCase()
        }
      }));
      added++;
    });
    if (added) {
      PB.addHistory(rec, 'Sequencing primers annotated',
        added + ' primer' + (added === 1 ? '' : 's') + ' added as features');
      PB.touch(rec);
      PB.renderAll();
    }
    PB.toast(added
      ? added + ' primer' + (added === 1 ? '' : 's') + ' added to ' + rec.name + ' as features'
      : (already ? 'Already annotated' : 'Nothing to add'), added ? 'ok' : 'warn');
  }
  PB.annotateSeqPrimers = annotatePrimers;

  /* put the planned primers into the working library so the map draws them */
  function markOnMap() {
    const plan = last.plan;
    const add = plan.plan.forward.concat(plan.plan.reverse);
    let n = 0;
    add.forEach(c => {
      const nm = (c.id ? c.id + ' ' : '') + (c.name || 'primer');
      if ((state.primers || []).some(p => (p.seq || '').toUpperCase() === c.seq.toUpperCase())) return;
      state.primers.push({ name: nm, seq: c.seq });
      n++;
    });
    PB.lsSet('primers', state.primers);
    state.primersShown = null;
    PB.renderAll();
    PB.toast(n ? n + ' primer' + (n === 1 ? '' : 's') + ' added to the library and drawn on the map'
               : 'Those primers are already in the library', 'ok');
  }

  function fillGap(which) {
    const plan = last.plan, rec = last.rec;
    const gaps = plan.gaps[which];
    if (!gaps.length) return;
    const strand = which === 'forward' ? 1 : -1;
    const made = [];
    gaps.forEach(g => {
      /* read into the hole from before it on the forward strand, and from
         after it on the reverse */
      const target = strand === 1 ? g.start : g.end;
      const d = B.designSeqPrimer(rec.seq, target, strand,
        { circular: last.circ, dead: last.dead });
      if (d) made.push({ gap: g, primer: d });
    });
    if (!made.length) {
      PB.toast('No primer with a single site could be found for that stretch.', 'warn');
      return;
    }
    PB.modal({
      title: 'Primers to order for the ' + which + ' gap' + (made.length === 1 ? '' : 's'),
      wide: true,
      body: '<div class="note">These are not in your stock — each is a new oligo to order. Every one ' +
        'anneals in exactly one place on this construct and sits clear of the unreadable window before ' +
        'the stretch it is meant to read.</div>' +
        '<div class="scrollx"><table class="tbl sm"><thead><tr><th>For bases</th><th>Sequence 5′→3′</th>' +
        '<th>nt</th><th>Tm</th><th>GC</th><th>Run-up</th></tr></thead><tbody>' +
        made.map(m => '<tr><td class="mono">' + fmt(m.gap.start + 1) + '–' + fmt(m.gap.end) + '</td>' +
          '<td class="mono">' + m.primer.seq + '</td><td class="num">' + m.primer.length + '</td>' +
          '<td class="num">' + f1(m.primer.tm) + '</td><td class="num">' + Math.round(m.primer.gc * 100) + '%</td>' +
          '<td class="num">' + m.primer.lead + '</td></tr>').join('') +
        '</tbody></table></div>',
      actions: [
        { label: 'Close' },
        { label: 'Copy', onClick: () => { PB.copyText(made.map(m =>
            'seq-' + which + '-' + (m.gap.start + 1) + '\t' + m.primer.seq).join('\n')); return false; } },
        { label: 'Add to my primers', primary: true, onClick: () => {
          made.forEach((m, i) => state.primers.push({
            name: 'seq ' + which + ' ' + (m.gap.start + 1) + '–' + m.gap.end, seq: m.primer.seq }));
          PB.lsSet('primers', state.primers);
          state.primersShown = null;
          PB.renderAll();
          PB.toast(made.length + ' primer' + (made.length === 1 ? '' : 's') + ' added to the library', 'ok');
        } }
      ]
    });
  }

  /* ---------------------------------------------------------------- export */
  /* See the note in app-express.js: this and the expression panel's plan
     writer had the same two names in two separate scopes and are not the
     same function. */
  function sequencingPlanText() {
    const plan = last.plan, rec = last.rec, cov = plan.coverage;
    const lines = [];
    lines.push('Sanger sequencing plan — ' + rec.name);
    lines.push('Region: ' + (plan.region.start + 1) + '–' + plan.region.end + '  (' + plan.region.length + ' bp)');
    lines.push('Read model: first ' + last.dead + ' bases unreadable, confident to ' + last.good +
               ', reaching ' + last.max);
    lines.push('Coverage: forward ' + cov.forward + '/' + cov.length + ', reverse ' + cov.reverse + '/' +
               cov.length + ', both strands ' + cov.both + '/' + cov.length);
    lines.push('');
    ['forward', 'reverse'].forEach(k => {
      lines.push(k === 'forward' ? 'Forward reads:' : 'Reverse reads:');
      if (!plan.plan[k].length) lines.push('  (nothing in stock)');
      plan.plan[k].forEach(c => {
        lines.push('  ' + (c.id ? c.id + '  ' : '') + (c.name || '') +
          (c.loc ? '  [' + c.loc + ']' : '') +
          '  anneals ' + (c.site.start + 1) + '-' + c.site.end +
          '  Tm ' + c.tm.toFixed(1) + '  covers ' + c.covers + ' bp');
        c.warnings.forEach(w => lines.push('      - ' + w));
      });
      const g = plan.gaps[k];
      if (g.length) lines.push('  not covered: ' + g.map(x =>
        (x.start + 1) + '-' + x.end).join(', '));
      lines.push('');
    });
    return lines.join('\n');
  }

  function sequencingExportPlan() {
    const plan = last.plan, rec = last.rec, cov = plan.coverage;
    if (!PB.xlsxBlob) { PB.copyText(sequencingPlanText()); PB.toast('Copied as text instead.', 'warn'); return; }
    const head = ['Use', 'Stock ID', 'Name', 'Where it is', 'Reads', 'Anneals from', 'Anneals to',
                  'Tm (C)', 'GC (%)', 'Length (nt)', 'Anneals (nt)', '5′ tail (nt)',
                  'Covers (bp)', 'Run-up (bp)', 'Mismatches', 'Notes', 'Sequence 5′→3′'];
    const row = (c, use) => [use, c.id || '', c.name || '', c.loc || '',
      c.strand === 1 ? 'forward' : 'reverse', c.site.start + 1, c.site.end,
      +c.tm.toFixed(1), Math.round(c.gc * 100), (c.seq || '').length, c.annealed, c.tail,
      c.covers, c.inside ? 0 : c.gap, c.mismatches.length,
      c.warnings.join('; '), (c.seq || '').toUpperCase()];
    const rows = [head];
    plan.plan.forward.forEach(c => rows.push(row(c, 'plan — forward')));
    plan.plan.reverse.forEach(c => rows.push(row(c, 'plan — reverse')));
    plan.candidates.forEach(c => {
      if (plan.plan.forward.indexOf(c) >= 0 || plan.plan.reverse.indexOf(c) >= 0) return;
      rows.push(row(c, 'also possible'));
    });
    plan.unusable.forEach(c => rows.push(row(c, 'refused')));

    const about = [['Sanger sequencing plan', ''],
      ['Construct', rec.name],
      ['Length', rec.seq.length + ' bp' + (last.circ ? ', circular' : ', linear')],
      ['Region read', (plan.region.start + 1) + '–' + plan.region.end],
      ['Region length', plan.region.length + ' bp'],
      ['', ''],
      ['Unreadable bases after the primer', last.dead],
      ['Confident read length', last.good],
      ['Furthest reach', last.max],
      ['', ''],
      ['Forward coverage', cov.forward + ' / ' + cov.length + ' bp'],
      ['Reverse coverage', cov.reverse + ' / ' + cov.length + ' bp'],
      ['Covered on both strands', cov.both + ' / ' + cov.length + ' bp'],
      ['Reactions in the plan', plan.plan.forward.length + plan.plan.reverse.length],
      ['Oligos searched', last.oligos.length],
      ['', ''],
      ['Gaps, forward', plan.gaps.forward.map(g => (g.start + 1) + '–' + g.end).join(', ') || 'none'],
      ['Gaps, reverse', plan.gaps.reverse.map(g => (g.start + 1) + '–' + g.end).join(', ') || 'none'],
      ['', ''],
      ['Tm model', 'nearest-neighbour, SantaLucia & Hicks 2004'],
      ['Salt', '50 mM Na+, 1.5 mM Mg2+, 0.2 mM dNTP, 0.25 uM oligo'],
      ['Coverage counted', 'only where the confident read reaches'],
      ['Made by', 'PlasmidBench']];

    const name = (rec.name || 'construct').replace(/[^\w.-]+/g, '_') + '_sequencing_plan.xlsx';
    PB.offerBlob(name, PB.xlsxBlob([
      { name: 'Plan', rows: rows, widths: [16, 11, 26, 20, 9, 13, 12, 8, 8, 11, 12, 11, 12, 11, 12, 54, 46] },
      { name: 'About', rows: about, widths: [36, 52] }
    ]), 'Sequencing plan');
  }

  /* the built-in example is useful on day one and clutter afterwards */
  PB.demoDialog = function () {
    const here = (PB.state.records || []).some(r => r.builtin);
    PB.modal({
      title: 'Demo construct',
      body: '<div class="note">' + (here
        ? 'The library holds a built-in example, <strong>pPB-T7-His6-EGFP (demo)</strong>. ' +
          'It is there so the tool opens on something rather than nothing, and it is read-only. ' +
          'Hide it once you have your own constructs loaded.'
        : 'The demo construct is hidden. Bring it back if you want something to try a tool on ' +
          'without touching your own plasmids.') + '</div>',
      actions: [
        { label: 'Close' },
        here
          ? { label: 'Hide it', primary: true, onClick: () => {
              const d = (PB.state.records || []).find(r => r.builtin);
              if (d) { PB.removeRecord(d.id); PB.toast('Demo construct hidden', 'ok'); }
            } }
          : { label: 'Bring it back', primary: true, onClick: () => PB.restoreDemo() }
      ]
    });
  };

  PB.sangerLast = () => last;
})(window.PB);
