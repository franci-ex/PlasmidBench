
/* ============================================================================
   PlasmidBench — expression strategy

   Paste a sequence, and this says how hard it is likely to be to get out in
   soluble form, then picks the construct to try from what the lab actually
   holds. Co-expression is the second half: several proteins in one cell is
   governed by hard constraints — origin families and selection markers —
   and those are checked rather than suggested.
   ============================================================================ */
(function (PB) {
  'use strict';
  const B = window.Bio;
  const $ = PB.$, esc = PB.esc, fmt = PB.fmt, f1 = PB.f1;
  const state = PB.state;

  /* --------------------------------------------------- the plasmid stock */
  state.plasmidStock = PB.lsGet('plasmidStock', []) || [];
  function savePlasmids() { PB.lsSet('plasmidStock', state.plasmidStock); }
  PB.plasmidCount = () => (state.plasmidStock || []).length;

  state.cellStock = PB.lsGet('cellStock', []) || [];
  function saveCells() { PB.lsSet('cellStock', state.cellStock); }
  PB.cellCount = () => (state.cellStock || []).length;

  PB.cellImportDialog = function () {
    const n = PB.cellCount();
    PB.modal({
      title: 'Competent cells',
      sub: n ? n + ' strains' : 'nothing imported yet',
      wide: true,
      body:
        '<div class="note">Paste the strain sheet. It needs a <strong>Strain</strong> column, and it will ' +
        'use an <strong>Existing plasmid</strong> column and a <strong>Resistance</strong> column if they ' +
        'are there \u2014 those decide what a strain has already spent before your vector goes in. A ' +
        '<strong>stock?</strong> column marks the ones you can actually use today.</div>' +
        '<div class="field"><label>Paste the sheet</label>' +
        '<textarea class="inp mono" id="csText" rows="10" spellcheck="false" ' +
        'placeholder="Location&#9;Strain&#9;Existing plasmid&#9;Use Case&#9;Resistance&#9;stock?"></textarea></div>' +
        '<div id="csPrev"></div>',
      actions: [
        { label: 'Close' },
        n ? { label: 'Clear', onClick: () => { state.cellStock = []; saveCells();
              PB.toast('Strain list cleared', 'ok'); PB.cellImportDialog(); return false; } } : null,
        { label: 'Check it', close: false, onClick: () => { previewCells(); return false; } },
        { label: 'Import', primary: true, onClick: () => {
            const ta = $('#csText');
            if (!ta || !ta.value.trim()) { PB.toast('Paste the sheet first.', 'warn'); return false; }
            let r; try { r = B.parseCellStock(ta.value); }
            catch (e) { PB.toast('That could not be read: ' + e.message, 'bad'); return false; }
            if (!r.strains.length) { PB.toast('No strains were found in that.', 'bad'); return false; }
            state.cellStock = r.strains; saveCells();
            PB.toast(r.strains.length + ' strains imported', 'ok');
            return true;
          } }
      ].filter(Boolean)
    });
  };

  function previewCells() {
    const ta = $('#csText'), host = $('#csPrev');
    if (!ta || !host) return;
    let r; try { r = B.parseCellStock(ta.value); }
    catch (e) { host.innerHTML = '<div class="note bad">' + esc(e.message) + '</div>'; return; }
    if (!r.strains.length) {
      host.innerHTML = '<div class="note bad">No header row with a <strong>Strain</strong> column was found.</div>';
      return;
    }
    host.innerHTML = '<div class="note ok">Read <strong>' + r.strains.length + '</strong> strains, ' +
      r.strains.filter(s => s.inStock).length + ' of them made up as competent cells.</div>' +
      '<div class="scrollx"><table class="tbl sm"><thead><tr><th>Strain</th><th>Carries</th>' +
      '<th>Spoken for</th><th>Use</th><th>Stock</th></tr></thead><tbody>' +
      r.strains.slice(0, 16).map(st => {
        const occ = B.hostCarries(st);
        const oris = occ.filter(c => c.ori).map(c => c.ori);
        return '<tr><td><strong>' + esc(st.name) + '</strong></td>' +
          '<td class="mono">' + esc(st.carries.map(c => c.plasmid).join(', ') || '\u2014') + '</td>' +
          '<td class="mono">' + esc((oris.join(' + ') || '\u2014') + '  /  ' + (st.markers.join(' + ') || '\u2014')) + '</td>' +
          '<td>' + esc(st.use || '') + '</td>' +
          '<td>' + (st.inStock ? '<span class="badge ok">yes</span>' : '<span class="hint">no</span>') + '</td></tr>';
      }).join('') + '</tbody></table></div>' +
      (r.strains.length > 16 ? '<div class="hint">\u2026 and ' + (r.strains.length - 16) + ' more.</div>' : '');
  }

  PB.plasmidImportDialog = function () {
    const n = PB.plasmidCount();
    PB.modal({
      title: 'Plasmid stock',
      sub: n ? n + ' plasmids' : 'nothing imported yet',
      wide: true,
      body:
        '<div class="note">Paste the sheet straight out of Excel. It needs a header row naming at ' +
        'least a <strong>Vector</strong> column, and it will use whatever else it recognises: construct, ' +
        'target, tag, resistance, shelf/box/position, and any note about how the thing expressed. ' +
        'Nothing leaves this browser.</div>' +
        '<div class="field"><label>Paste the stock</label>' +
        '<textarea class="inp mono" id="psText" rows="10" spellcheck="false" ' +
        'placeholder="Plasmid ID&#9;Construct ID&#9;Vector&#9;Resistance&#9;Tag&#9;Expression?"></textarea></div>' +
        '<div id="psPrev"></div>',
      actions: [
        { label: 'Close' },
        n ? { label: 'Clear', onClick: () => {
            state.plasmidStock = []; savePlasmids();
            PB.toast('Plasmid stock cleared', 'ok'); PB.plasmidImportDialog(); return false; } } : null,
        { label: 'Check it', close: false, onClick: () => { previewPlasmids(); return false; } },
        { label: 'Import', primary: true, onClick: () => {
            const r = parsePlasmidPaste();
            if (!r) return false;
            if (!r.plasmids.length) { PB.toast('No rows with a vector in them were found.', 'bad'); return false; }
            state.plasmidStock = r.plasmids;
            savePlasmids();
            PB.toast(r.plasmids.length + ' plasmids imported', 'ok');
            return true;
          } }
      ].filter(Boolean)
    });
  };

  function parsePlasmidPaste() {
    const ta = $('#psText');
    if (!ta || !ta.value.trim()) { PB.toast('Paste the sheet first.', 'warn'); return null; }
    try { return B.parsePlasmidStock(ta.value); }
    catch (e) { PB.toast('That could not be read: ' + e.message, 'bad'); return null; }
  }

  function previewPlasmids() {
    const r = parsePlasmidPaste(); const host = $('#psPrev');
    if (!r || !host) return;
    if (!r.plasmids.length) {
      host.innerHTML = '<div class="note bad">No header row with a <strong>Vector</strong> column was found, ' +
        'so nothing could be read. Include the header line when you copy.</div>';
      return;
    }
    const known = r.plasmids.filter(p => p.kb).length;
    const sum = B.vectorSummary(r.plasmids);
    host.innerHTML =
      '<div class="note ok">Read <strong>' + r.plasmids.length + '</strong> plasmids across <strong>' +
      sum.length + '</strong> vector families. ' + known + ' of them (' +
      Math.round(known / r.plasmids.length * 100) + '%) are backbones the tool knows the origin and marker for' +
      (r.skipped.length ? '; ' + r.skipped.length + ' rows had no vector named and were left out' : '') + '.</div>' +
      '<div class="scrollx"><table class="tbl sm"><thead><tr><th>Vector</th><th>What it is</th>' +
      '<th>Origin / marker</th><th>Copies</th><th class="num">Have</th><th>Your own records</th></tr></thead><tbody>' +
      sum.slice(0, 14).map(v => '<tr><td class="mono">' + esc(Object.keys(v.spellings)[0]) + '</td>' +
        '<td>' + esc(v.kb ? v.kb.name : 'not recognised') + '</td>' +
        '<td class="mono">' + (v.kb ? esc(v.kb.ori + ' / ' + v.kb.marker) : '—') + '</td>' +
        '<td class="num">' + (v.kb && v.kb.copy ? '~' + v.kb.copy : '—') + '</td>' +
        '<td class="num">' + v.count + '</td>' +
        '<td>' + recordCell(v) + '</td></tr>').join('') +
      '</tbody></table></div>' +
      (sum.length > 14 ? '<div class="hint">… and ' + (sum.length - 14) + ' more families.</div>' : '');
  }

  function recordCell(v) {
    const n = v.soluble + v.insoluble + v.noExpression;
    if (!n) return '<span class="hint">nothing recorded</span>';
    const bits = [];
    if (v.soluble) bits.push('<span class="badge ok">' + v.soluble + ' soluble</span>');
    if (v.insoluble) bits.push('<span class="badge bad">' + v.insoluble + ' insoluble</span>');
    if (v.noExpression) bits.push('<span class="badge warn">' + v.noExpression + ' no expression</span>');
    return bits.join(' ');
  }

  /* ===================================================================== */
  let lastPlan = null;

  PB.toolExpression = function () {
    const n = PB.plasmidCount();
    const sel = state.selection, rec = PB.active();
    PB.modal({
      title: 'Expression strategy',
      sub: n ? n + ' plasmids on the shelf' : 'no plasmid stock imported',
      wide: true,
      body:
        (n ? '' : '<div class="note warn">No plasmid stock has been imported, so there is nothing to ' +
          'recommend <em>from</em>. The reading of the sequence still works. Import your collection from ' +
          '<strong>Plasmid stock</strong> to have it choose among what you actually have.</div>') +
        '<div class="row">' +
          '<div class="field"><label>What are you expressing</label>' +
            '<select class="inp" id="exSrc">' +
              (rec ? '<option value="rec">' + esc(rec.name) + (sel ? ' — the selection' : ' — the whole thing') + '</option>' : '') +
              '<option value="paste">a sequence I paste</option>' +
            '</select></div>' +
          '<div class="field"><label>Host</label><select class="inp" id="exHost">' +
            '<option value="auto">let the tool rank my strains</option>' + hostOptions() +
            '</select></div>' +
          '<div class="field narrow"><label>Buffer pH</label>' +
            '<input class="inp mono" id="exPH" type="number" step="0.1" min="4" max="10" value="7.5"></div>' +
        '</div>' +
        '<div class="field" id="exPasteF" style="display:none"><label>Sequence — DNA or protein</label>' +
          '<textarea class="inp mono" id="exPaste" rows="4" spellcheck="false" ' +
          'placeholder="ATG… or MTYKLILNGKTLKGETTTEAV…"></textarea></div>' +
        '<div class="field"><label>Anything co-expressed with it (optional, one per line: name then sequence, or just names)</label>' +
          '<textarea class="inp mono" id="exPartners" rows="3" spellcheck="false" ' +
          'placeholder="CDK2&#10;cyclin A"></textarea></div>' +
        '<div class="row"><label class="chk"><input type="checkbox" id="exChap"> ' +
          'I want to add a chaperone plasmid (GroEL/ES)</label></div>' +
        '<div id="exOut"></div>',
      actions: [
        { label: 'Close' },
        { label: 'Import plasmid stock…', onClick: () => { PB.plasmidImportDialog(); return false; } },
        { label: 'Work it out', primary: true, close: false, onClick: () => { runExpression(); return false; } }
      ]
    });
    const src = $('#exSrc');
    const sync = () => { $('#exPasteF').style.display = src.value === 'paste' ? '' : 'none'; };
    src.addEventListener('change', sync); sync();
    if (!rec) { src.value = 'paste'; sync(); }
  };

  /* the strains they actually have, in-stock expression strains first */
  function hostOptions() {
    const cells = (state.cellStock || []).filter(s => s.expression !== false);
    if (!cells.length)
      return ['BL21(DE3)', 'LOBSTR', 'LOBSTR-RIL', 'Rosetta2(DE3)', 'SHuffle T7', 'Origami',
              'ArcticExpress', 'C41(DE3)'].map(h => '<option value="n:' + esc(h) + '">' + esc(h) + '</option>').join('');
    const rank = cells.map((s, i) => ({ s: s, i: i }))
      .sort((a, b) => (b.s.inStock ? 1 : 0) - (a.s.inStock ? 1 : 0));
    return rank.map(x => '<option value="i:' + x.i + '">' + esc(x.s.name) +
      (x.s.carries.length ? ' \u2014 carries ' + esc(x.s.carries.map(c => c.plasmid).join(', ')) : '') +
      (x.s.inStock ? '' : '  (not made up)') + '</option>').join('');
  }
  function hostChosen() {
    const v = $('#exHost').value || '';
    if (v === 'auto') return null;                 /* the tool picks */
    if (v.indexOf('i:') === 0) return (state.cellStock || [])[+v.slice(2)] || null;
    return v.slice(2) || v;
  }

  function sequenceIn() {
    const how = $('#exSrc').value;
    if (how === 'paste') return ($('#exPaste').value || '').trim();
    const rec = PB.active();
    if (!rec) return '';
    const sel = state.selection;
    if (sel) {
      const a = Math.min(sel.start, sel.end), b = Math.max(sel.start, sel.end);
      return rec.seq.slice(a, b);
    }
    return rec.seq;
  }

  function runExpression() {
    const host = $('#exOut');
    const seq = sequenceIn();
    if (!seq) { host.innerHTML = '<div class="note bad">There is no sequence to read.</div>'; return; }
    const hostPick = hostChosen();
    const hostName = !hostPick ? null : (typeof hostPick === 'string' ? hostPick : hostPick.name);
    const pH = parseFloat($('#exPH').value) || 7.5;
    const chap = $('#exChap').checked;
    const partnersRaw = ($('#exPartners').value || '').split('\n').map(s => s.trim()).filter(Boolean);
    host.innerHTML = '<div class="note">Working…</div>';
    setTimeout(() => {
      let outlook;
      try { outlook = B.solubilityOutlook(seq, { bufferPH: pH }); }
      catch (e) { host.innerHTML = '<div class="note bad">That sequence could not be read: ' + esc(e.message) + '</div>'; return; }
      if (!outlook) { host.innerHTML = '<div class="note bad">That did not look like a sequence.</div>'; return; }
      const inv = state.plasmidStock || [];
      const partnerNames = partnersRaw.map(l => l.split(/[\s,;:\t]/)[0]);
      const strat = B.strategyFor(outlook, inv, { host: hostName, partners: partnerNames });
      const strains = (state.cellStock || []).length ? B.pickStrains(outlook, state.cellStock) : null;
      /* the cross-product: every backbone against every strain, scored */
      let ranked = null;
      try {
        ranked = B.rankCombinations(outlook, inv, state.cellStock || [], { partners: partnerNames });
        if (hostPick && typeof hostPick !== 'string') {
          /* a strain was insisted on, so only pairings in it are offered */
          const only = ranked.everyPairing.filter(c => c.host.name === hostPick.name);
          if (only.length) ranked.combinations = only;
        }
      } catch (e) { ranked = null; }
      /* the partners: a name, and a sequence if one was given */
      const partners = partnersRaw.map(line => {
        const m = line.match(/^(.*?)[\s,;:\t]+([ACGTUacgtu]{30,}|[ACDEFGHIKLMNPQRSTVWY]{20,})$/);
        if (m) return { name: m[1].trim(), seq: m[2] };
        return { name: line };
      });
      let coex = null;
      if (partners.length || chap) {
        const all = [{ name: outlook.props.length + ' aa target', want: 1 }]
          .concat(partners.map(p => ({ name: p.name, want: 1 })));
        coex = B.coexpressionPlan(all, inv, { host: hostPick, chaperones: chap });
        coex.partnerOutlooks = partners.filter(p => p.seq).map(p => ({
          name: p.name, outlook: B.solubilityOutlook(p.seq, { bufferPH: pH })
        }));
      }
      lastPlan = { outlook: outlook, strategy: strat, coex: coex, host: hostName,
                   hostPick: hostPick, strains: strains, ranked: ranked, pH: pH, inv: inv };
      host.innerHTML = render(lastPlan);
      wire();
    }, 10);
  }

  /* --------------------------------------------------------------- render */
  const VERDICT = { good: 'ok', watch: 'warn', bad: 'bad' };

  function bandNote(o) {
    const cls = o.band === 'straightforward' ? 'ok' : o.band === 'tricky' ? 'warn' : 'bad';
    const head = o.band === 'wrong target' ? 'This is not a soluble-expression problem'
      : o.band === 'hard' ? 'Expect a fight'
      : o.band === 'tricky' ? 'Worth taking seriously'
      : 'Should be straightforward';
    return '<div class="note ' + cls + '"><strong>' + head + '.</strong> ' + esc(o.summary) + '</div>';
  }

  function factorTable(o) {
    return '<div class="scrollx"><table class="tbl sm"><thead><tr><th></th><th>Reading</th><th>Value</th>' +
      '<th>What it means</th></tr></thead><tbody>' +
      o.factors.map(f => '<tr><td><span class="badge ' + VERDICT[f.verdict] + '">' +
        (f.verdict === 'good' ? 'ok' : f.verdict === 'watch' ? 'watch' : 'problem') + '</span></td>' +
        '<td>' + esc(f.name) + '</td><td class="mono">' + esc(f.value) + '</td>' +
        '<td>' + esc(f.why) + '</td></tr>').join('') +
      '</tbody></table></div>';
  }

  /* group the matching plasmids by the vector they sit on, since a list of
     six hundred constructs is no use to anyone */
  function shelfFor(have) {
    const by = {};
    have.forEach(p => {
      const k = p.vector;
      (by[k] = by[k] || { vector: k, n: 0, where: null, kb: p.kb }).n++;
      if (!by[k].where && p.loc) by[k].where = p.loc;
    });
    return Object.keys(by).map(k => by[k]).sort((a, b) => b.n - a.n);
  }

  function optionCard(opt, i) {
    if (!opt.tag) return '<div class="note bad">' + esc(opt.why) + '</div>';
    const shelf = shelfFor(opt.have);
    const info = opt.info || {};
    return '<tr><td><span class="badge ' + (i === 0 ? 'ok' : 'neutral') + '">' + (i === 0 ? 'first choice' : '#' + (i + 1)) + '</span></td>' +
      '<td><strong>' + esc(info.label || opt.tag) + '</strong>' +
      (info.size > 1 ? ' <span class="hint">' + info.size + ' aa</span>' : '') + '</td>' +
      '<td>' + esc(opt.why) + (info.note ? ' <span class="hint">' + esc(info.note) + '</span>' : '') + '</td>' +
      '<td>' + (shelf.length
        ? shelf.slice(0, 3).map(s => '<div class="mono">' + esc(s.vector) + ' <span class="hint">×' + s.n +
            (s.where ? ', e.g. ' + esc(s.where) : '') + '</span></div>').join('') +
          (shelf.length > 3 ? '<div class="hint">+ ' + (shelf.length - 3) + ' more</div>' : '')
        : '<span class="badge warn">none on the shelf</span>') + '</td></tr>';
  }

  function render(p) {
    const o = p.outlook, s = p.strategy;
    let h = bandNote(o);
    h += '<div class="chipbar">' +
      '<span class="chip">' + fmt(o.props.length) + ' residues</span>' +
      '<span class="chip">' + f1(o.props.mwKDa) + ' kDa</span>' +
      '<span class="chip">pI ' + o.props.pI.toFixed(2) + '</span>' +
      '<span class="chip">GRAVY ' + o.props.gravy.toFixed(2) + '</span>' +
      '<span class="chip">' + o.props.cys + ' Cys</span>' +
      (o.tm.length ? '<span class="chip bad">' + o.tm.length + ' hydrophobic stretch' + (o.tm.length === 1 ? '' : 'es') + '</span>' : '') +
      '</div>';
    h += '<h4>What the sequence says</h4>' + factorTable(o);

    if (p.ranked && p.ranked.combinations.length) h += combinationBlock(p.ranked);
    else {
      h += '<h4>What to build, from what you have</h4>';
      const real = s.options.filter(x => x.tag);
      h += real.length
        ? '<div class="scrollx"><table class="tbl sm"><thead><tr><th></th><th>Fusion</th><th>Why this one</th>' +
          '<th>On the shelf</th></tr></thead><tbody>' +
          real.map((opt, i) => optionCard(opt, i)).join('') + '</tbody></table></div>'
        : s.options.map(x => optionCard(x, 0)).join('');
    }

    h += '<h4>How to run it</h4><div class="scrollx"><table class="tbl sm"><tbody>' +
      s.conditions.map(c => '<tr><td style="white-space:nowrap"><strong>' + esc(c.what) + '</strong></td>' +
        '<td>' + esc(c.detail) + '</td></tr>').join('') + '</tbody></table></div>';

    if (p.strains && p.strains.length) h += strainBlock(p);
    if (p.coex) h += coexBlock(p.coex);

    h += '<div class="row" style="margin-top:10px">' +
      '<button class="btn" id="exXls">Export the plan as a spreadsheet</button>' +
      '<button class="btn ghost" id="exCopy">Copy as text</button></div>';
    h += '<div class="hint">These readings come from composition alone, and composition cannot tell you ' +
      'whether a protein folds. They are here to point at the specific thing most likely to go wrong, so ' +
      'the first construct is a considered one rather than a guess. The origin and marker conflicts above ' +
      'are not predictions — those are hard constraints.</div>';
    return h;
  }

  const short = (t, n) => { t = String(t || ''); return t.length > n ? t.slice(0, n - 1) + '\u2026' : t; };

  function elementBar(inv) {
    const row = (label, m) => {
      const keys = Object.keys(m).sort((a, b) => m[b] - m[a]);
      if (!keys.length) return '';
      return '<tr><td style="white-space:nowrap"><strong>' + label + '</strong></td><td>' +
        keys.map(k => '<span class="chip">' + esc(k) + ' <span class="hint">\u00d7' + m[k] + '</span></span>').join(' ') +
        '</td></tr>';
    };
    return '<h4>What your collection contains</h4>' +
      '<div class="scrollx"><table class="tbl sm"><tbody>' +
      row('Origins', inv.origins) + row('Promoters', inv.promoters) +
      row('Markers', inv.markers) + row('Fusion tags', inv.tags) +
      '</tbody></table></div>';
  }

  function combinationBlock(r) {
    const rows = r.combinations.slice(0, 10);
    return elementBar(r.inventory) +
      '<h4>Every combination, ranked</h4>' +
      '<div class="note">' + r.counts.considered + ' backbone-and-strain pairings were considered, ' +
      r.counts.valid + ' of them workable across ' + r.counts.backbones + ' distinct backbones. ' +
      (r.rejected.length ? r.rejected.length + ' were struck out for a reason you can read below. ' : '') +
      'One line per backbone, each with the strain that suits it best.</div>' +
      '<div class="scrollx"><table class="tbl sm"><thead><tr><th></th><th>Backbone</th><th>Strain</th>' +
      '<th>Element</th><th>Why this one</th><th>Where</th></tr></thead><tbody>' +
      rows.map((c, i) => '<tr><td>' + (i === 0 ? '<span class="badge ok">best fit</span>'
          : '<span class="badge neutral">#' + (i + 1) + '</span>') + '</td>' +
        '<td><strong class="mono">' + esc(short(c.spelling, 26)) + '</strong> <span class="hint">\u00d7' + c.have + '</span>' +
        '<div class="hint">' + esc(c.kb.name) + '</div></td>' +
        '<td>' + esc(c.host.name) + (c.host.inStock ? '' : ' <span class="badge warn">not made up</span>') + '</td>' +
        '<td class="mono">' + esc(c.kb.ori + ' / ' + c.kb.marker) + '<div class="hint">' + esc(c.kb.promoter || '') + '</div></td>' +
        '<td>' + esc(c.why.join('; ')) + '</td>' +
        '<td class="mono">' + esc(short(c.where || '\u2014', 24)) + '</td></tr>').join('') +
      '</tbody></table></div>' +
      (r.rejected.length
        ? '<h4>Struck out, and why</h4><div class="scrollx"><table class="tbl sm"><thead><tr>' +
          '<th>Backbone</th><th>Strain</th><th>Why not</th></tr></thead><tbody>' +
          r.rejected.slice(0, 8).map(c => '<tr><td class="mono">' + esc(short(c.spelling, 24)) + '</td>' +
            '<td>' + esc(c.host.name) + '</td><td>' + esc(c.against[0]) + '</td></tr>').join('') +
          '</tbody></table></div>'
        : '');
  }

  function strainBlock(p) {
    const ranked = p.strains.filter(x => x.strain.expression !== false).slice(0, 6);
    if (!ranked.length) return '';
    return '<h4>Which cells to use</h4>' +
      '<div class="scrollx"><table class="tbl sm"><thead><tr><th></th><th>Strain</th>' +
      '<th>Already spent</th><th>Why</th></tr></thead><tbody>' +
      ranked.map((x, i) => {
        const st = x.strain;
        const occ = B.hostCarries(st);
        const oris = occ.filter(c => c.ori).map(c => c.ori);
        return '<tr><td>' + (i === 0 ? '<span class="badge ok">first choice</span>'
                : st.inStock ? '<span class="badge neutral">#' + (i + 1) + '</span>'
                : '<span class="badge warn">not made up</span>') + '</td>' +
          '<td><strong>' + esc(st.name) + '</strong>' +
          (st.loc ? '<div class="hint">' + esc(st.loc) + '</div>' : '') + '</td>' +
          '<td class="mono">' + esc(oris.length || st.markers.length
            ? (oris.join(' + ') || '\u2014') + ' / ' + (st.markers.join(' + ') || '\u2014')
            : 'nothing \u2014 every origin and drug still free') + '</td>' +
          '<td>' + esc(x.why.join('; ') || 'a plain expression strain') + '</td></tr>';
      }).join('') + '</tbody></table></div>' +
      '<div class="hint">A strain that already carries a plasmid has spent that origin and that drug ' +
      'before your vector goes in, which is why it is listed here rather than left to be discovered later.</div>';
  }

  function coexBlock(c) {
    let h = '<h4>Co-expression</h4>';
    if (c.compatibility.problems.length)
      h += c.compatibility.problems.map(x =>
        '<div class="note bad"><strong>' + esc(x.what) + '</strong><br>' + esc(x.why) + '</div>').join('');
    if (c.assigned.length)
      h += '<div class="scrollx"><table class="tbl sm"><thead><tr><th>Partner</th><th>Backbone to use</th>' +
        '<th>Origin / marker</th><th class="num">Copies</th><th>On the shelf</th></tr></thead><tbody>' +
        c.assigned.map(a => '<tr><td><strong>' + esc(a.partner.name || '?') + '</strong></td>' +
          '<td>' + esc(a.slot.vector) + '</td>' +
          '<td class="mono">' + esc(a.slot.ori + ' / ' + a.slot.marker) + '</td>' +
          '<td class="num">~' + a.slot.copy + '</td>' +
          '<td>' + (a.haveCount ? a.haveCount + ' plasmid' + (a.haveCount === 1 ? '' : 's') + ' already on this backbone'
                                : '<span class="badge warn">nothing on this backbone</span>') + '</td></tr>').join('') +
        '</tbody></table></div>';
    (c.advice || []).forEach(a => {
      h += '<div class="note ' + (a.kind === 'blocking' ? 'bad' : '') + '">' + esc(a.text) + '</div>';
    });
    (c.compatibility.notes || []).forEach(t => { h += '<div class="hint">' + esc(t) + '</div>'; });
    if (c.partnerOutlooks && c.partnerOutlooks.length) {
      h += '<h4>The partners, read the same way</h4><div class="scrollx"><table class="tbl sm"><thead><tr>' +
        '<th>Partner</th><th>Size</th><th>pI</th><th>GRAVY</th><th>Reading</th></tr></thead><tbody>' +
        c.partnerOutlooks.filter(x => x.outlook).map(x => '<tr><td>' + esc(x.name) + '</td>' +
          '<td class="num">' + f1(x.outlook.props.mwKDa) + ' kDa</td>' +
          '<td class="num">' + x.outlook.props.pI.toFixed(2) + '</td>' +
          '<td class="num">' + x.outlook.props.gravy.toFixed(2) + '</td>' +
          '<td><span class="badge ' + (x.outlook.band === 'straightforward' ? 'ok' :
            x.outlook.band === 'tricky' ? 'warn' : 'bad') + '">' + esc(x.outlook.band) + '</span></td></tr>').join('') +
        '</tbody></table></div>';
    }
    return h;
  }

  function wire() {
    const x = $('#exXls'); if (x) x.addEventListener('click', exportPlan);
    const c = $('#exCopy'); if (c) c.addEventListener('click', () => PB.copyText(planText()));
  }

  function planText() {
    const p = lastPlan, o = p.outlook;
    const L = [];
    L.push('Expression strategy');
    L.push('Protein: ' + o.props.length + ' aa, ' + o.props.mwKDa.toFixed(1) + ' kDa, pI ' +
           o.props.pI.toFixed(2) + ', GRAVY ' + o.props.gravy.toFixed(3) + ', ' + o.props.cys + ' Cys');
    L.push('Host: ' + p.host + '   buffer pH ' + p.pH);
    L.push('Reading: ' + o.band.toUpperCase() + ' — ' + o.summary);
    L.push('');
    o.factors.forEach(f => L.push('  [' + f.verdict + '] ' + f.name + ': ' + f.value));
    L.push('');
    L.push('Suggested constructs:');
    p.strategy.options.filter(x => x.tag).forEach((x, i) => {
      L.push('  ' + (i + 1) + '. ' + ((x.info && x.info.label) || x.tag) + ' — ' + x.why);
      shelfFor(x.have).slice(0, 3).forEach(s => L.push('       ' + s.vector + ' x' + s.n + (s.where ? '  [' + s.where + ']' : '')));
    });
    L.push('');
    L.push('Conditions:');
    p.strategy.conditions.forEach(c => L.push('  - ' + c.what + ': ' + c.detail));
    if (p.coex) {
      L.push('');
      L.push('Co-expression:');
      p.coex.compatibility.problems.forEach(x => L.push('  ! ' + x.what + ' ' + x.why));
      p.coex.assigned.forEach(a => L.push('  ' + (a.partner.name || '?') + ' -> ' + a.slot.vector +
        ' (' + a.slot.ori + '/' + a.slot.marker + ', ~' + a.slot.copy + ' copies)'));
      p.coex.advice.forEach(a => L.push('  - ' + a.text));
    }
    return L.join('\n');
  }

  function exportPlan() {
    const p = lastPlan, o = p.outlook;
    if (!PB.xlsxBlob) { PB.copyText(planText()); PB.toast('Copied as text instead.', 'warn'); return; }
    const reading = [['Reading', 'Value', 'Verdict', 'What it means']]
      .concat(o.factors.map(f => [f.name, f.value, f.verdict, f.why]));
    const opts = [['Rank', 'Fusion', 'Size (aa)', 'Why', 'Vector on the shelf', 'How many', 'Where']];
    p.strategy.options.filter(x => x.tag).forEach((x, i) => {
      const shelf = shelfFor(x.have);
      if (!shelf.length) opts.push([i + 1, (x.info && x.info.label) || x.tag, (x.info && x.info.size) || '', x.why, 'none', 0, '']);
      shelf.slice(0, 5).forEach(s => opts.push([i + 1, (x.info && x.info.label) || x.tag,
        (x.info && x.info.size) || '', x.why, s.vector, s.n, s.where || '']));
    });
    const cond = [['Step', 'Detail']].concat(p.strategy.conditions.map(c => [c.what, c.detail]));
    const about = [['Expression strategy', ''],
      ['Protein length', o.props.length + ' aa'],
      ['Molecular weight', +o.props.mwKDa.toFixed(2) + ' kDa'],
      ['Isoelectric point', +o.props.pI.toFixed(2)],
      ['GRAVY', +o.props.gravy.toFixed(3)],
      ['Instability index', +o.props.instability.toFixed(1)],
      ['Aliphatic index', +o.props.aliphatic.toFixed(1)],
      ['Cysteines', o.props.cys],
      ['Predicted membrane helices', o.tm.length],
      ['Host', p.host], ['Buffer pH', p.pH],
      ['Overall reading', o.band],
      ['', ''],
      ['Hydropathy scale', 'Kyte & Doolittle 1982'],
      ['Instability index', 'Guruprasad, Reddy & Pandit 1990'],
      ['Composition solubility index', 'Wilkinson & Harrison 1991'],
      ['Weight, pI, GRAVY and instability', 'checked against Biopython ProtParam to the last decimal'],
      ['Caution', 'Composition cannot decide whether a protein folds. Origin and marker conflicts are hard constraints; the rest is guidance.'],
      ['Made by', 'PlasmidBench']];
    const sheets = [
      { name: 'Reading', rows: reading, widths: [30, 34, 10, 86] },
      { name: 'Constructs', rows: opts, widths: [7, 20, 10, 70, 24, 11, 22] },
      { name: 'Conditions', rows: cond, widths: [28, 96] },
      { name: 'About', rows: about, widths: [36, 60] }
    ];
    if (p.coex) {
      const cx = [['Partner', 'Backbone', 'Origin', 'Marker', 'Copies', 'On the shelf']];
      p.coex.assigned.forEach(a => cx.push([a.partner.name || '', a.slot.vector, a.slot.ori,
        a.slot.marker, a.slot.copy, a.haveCount]));
      cx.push(['', '', '', '', '', '']);
      p.coex.compatibility.problems.forEach(x => cx.push(['CONFLICT', x.what, x.why, '', '', '']));
      p.coex.advice.forEach(a => cx.push(['note', a.text, '', '', '', '']));
      sheets.splice(3, 0, { name: 'Co-expression', rows: cx, widths: [22, 30, 12, 10, 9, 14] });
    }
    PB.offerBlob('expression_strategy.xlsx', PB.xlsxBlob(sheets), 'Expression strategy');
  }

  PB.expressionLast = () => lastPlan;
})(window.PB);

