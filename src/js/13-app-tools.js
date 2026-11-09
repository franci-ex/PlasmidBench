
/* ============================================================================
   PlasmidBench — tools: cloning simulators, mutagenesis, design, import/export
   ============================================================================ */
(function (PB) {
  'use strict';
  const B = window.Bio;
  const $ = PB.$, esc = PB.esc, fmt = PB.fmt, f1 = PB.f1, f2 = PB.f2;
  const state = PB.state;

  /* ---------------- shared dialog helpers ---------------- */
  function recOptions(selId, allowSel) {
    let h = state.records.map(r => '<option value="' + r.id + '"' + (r.id === selId ? ' selected' : '') + '>' +
      esc(r.name.length > 30 ? r.name.slice(0, 29) + '…' : r.name) + ' · ' + fmt(r.seq.length) + ' bp' +
      (r.topology === 'circular' ? ' ○' : '') + '</option>').join('');
    if (allowSel && state.selection) h += '<option value="__sel">current selection (' +
      fmt(state.selection.end - state.selection.start) + ' bp)</option>';
    return h;
  }

  function seqField(id, label, opts) {
    opts = opts || {};
    const sel = opts.defaultId || (PB.active() ? PB.active().id : '');
    return '<div class="field" style="flex:1 1 ' + (opts.flex || '290px') + '"><label>' + esc(label) + '</label>' +
      '<select class="inp" id="' + id + '_src" style="margin-bottom:5px">' +
        recOptions(sel, opts.allowSel !== false) +
        '<option value="">paste below…</option>' +
      '</select>' +
      '<textarea class="inp" id="' + id + '_txt" rows="' + (opts.rows || 3) + '" spellcheck="false" placeholder="' +
        esc(opts.placeholder || 'or paste sequence / GenBank here') + '"></textarea></div>';
  }

  function seqValue(id) {
    const srcEl = $('#' + id + '_src'), txtEl = $('#' + id + '_txt');
    const src = srcEl ? srcEl.value : '';
    if (src === '__sel') {
      const rec = PB.active();
      if (rec && state.selection) {
        // a selection is a linear piece of the record: it is not the record,
        // but it carries the record's features and methylation
        return { seq: PB.selectionSeq(rec, state.selection),
                 rec: null, from: rec, circular: false, name: rec.name + ' selection' };
      }
    } else if (src) {
      const r = state.records.find(x => x.id === src);
      if (r) return { seq: r.seq, rec: r, circular: r.topology === 'circular', name: r.name };
    }
    const raw = txtEl ? txtEl.value : '';
    if (/LOCUS|ORIGIN/.test(raw)) {
      const recs = B.parseGenBank(raw);
      if (recs.length && recs[0].seq) {
        return { seq: recs[0].seq, rec: PB.makeRecord(recs[0]), circular: recs[0].topology === 'circular', name: recs[0].name };
      }
    }
    const s = B.cleanSeq(raw.replace(/^>.*$/gm, ''));
    return { seq: s, rec: null, circular: false, name: 'pasted sequence' };
  }

  function enzField(id, label, dflt, multi) {
    return '<div class="field" style="flex:1 1 170px"><label>' + esc(label) + '</label>' +
      '<input class="inp mono" id="' + id + '" value="' + esc(dflt || '') + '" placeholder="' +
      (multi ? 'EcoRI HindIII' : 'BsaI') + '" spellcheck="false"></div>';
  }
  function enzNames(id) {
    return ($('#' + id).value || '').split(/[\s,+]+/).filter(Boolean)
      .map(n => { const e = PB.enzByName(n); return e ? e.name : null; }).filter(Boolean);
  }

  function result(html) {
    const host = $('#toolResult');
    if (host) host.innerHTML = html;
  }
  function resultSlot() { return '<div id="toolResult"></div>'; }

  /* ---------------- feature transfer ---------------- */
  // The contiguous stretch of the source that holds a feature (exons and the
  // introns between them), with each piece's offset inside that stretch.
  function featureBlock(rec, f) {
    const L = rec.seq.length, circ = rec.topology === 'circular';
    let pieces, s0;
    if (f.mixed && f.parts) {
      // keep biological order; the block starts at the leftmost part
      pieces = f.parts.map(p => ({ start: p.start, end: p.end, strand: p.strand }));
      s0 = Math.min.apply(null, pieces.map(p => p.start));
    } else {
      pieces = B.featurePieces(f, L).map(p => ({ start: p.start, end: p.end, strand: f.strand }));
      if (pieces.length) s0 = pieces[0].start;
    }
    if (!pieces.length) return null;
    const off = x => ((x - s0) % L + L) % L;
    let len = 0;
    const rel = pieces.map(p => {
      const o = off(p.start), n = p.end - p.start;
      len = Math.max(len, o + n);
      return { off: o, len: n, strand: p.strand };
    });
    if (len > L) return null;
    const e0 = s0 + len;
    const text = e0 <= L ? rec.seq.slice(s0, e0) : (circ ? rec.seq.slice(s0) + rec.seq.slice(0, e0 - L) : null);
    return text ? { text, rel, start: s0, L: L, circ: circ } : null;
  }
  function transferFeatures(productSeq, sources, circular) {
    const out = [];
    const placed = new Set();          // so two identical copies get two slots
    const P = productSeq.length;
    const hay = circular ? productSeq + productSeq : productSeq;
    for (const src of sources) {
      if (!src || !src.features) continue;
      for (const f of src.features) {
        if (f.type === 'source' || f.between) continue;
        const blk = featureBlock(src, f);
        if (!blk) continue;
        const sub = blk.text;
        if (!sub.length || sub.length > P) continue;
        /* A short feature — a ribosome binding site, a start codon, a 4 nt
           overhang — cannot be placed by its own sequence: six bases turn up
           all over a plasmid. It is placed by the sequence AROUND it instead,
           which is unambiguous, and the flank is shortened until it fits
           inside the product.                                             */
        /* Every feature is placed by the sequence AROUND it first. That is
           what tells two identical copies apart — a vector with two T7
           promoters has two of them in the product, not one — and it is the
           only way to place something as short as a six-base ribosome
           binding site. Only when the context straddles a junction, and the
           feature is long enough to be unambiguous on its own, does the bare
           sequence decide.                                                */
        let i = -1, rev = false;
        for (let fl = 12; fl >= 4 && i === -1; fl -= 4) {
          if (sub.length + 2 * fl > P) continue;
          if (!blk.circ && (blk.start - fl < 0 || blk.start + sub.length + fl > blk.L)) continue;
          const a = ((blk.start - fl) % blk.L + blk.L) % blk.L;
          const probe = (src.seq + src.seq).substr(a, sub.length + 2 * fl);
          if (probe.length < sub.length + 2 * fl) continue;
          let j = hay.indexOf(probe), r2 = false;
          if (j === -1) { j = hay.indexOf(B.revcomp(probe)); r2 = j !== -1; }
          while (j !== -1) {
            const at = (r2 ? j + fl : j + fl) % P;
            const key = at + ':' + f.type + ':' + (f.label || '');
            if (!placed.has(key)) { i = j + fl; rev = r2; break; }
            j = hay.indexOf(r2 ? B.revcomp(probe) : probe, j + 1);
          }
        }
        if (i === -1 && sub.length >= 12) {
          let j = hay.indexOf(sub), r2 = false;
          if (j === -1) { j = hay.indexOf(B.revcomp(sub)); r2 = j !== -1; }
          while (j !== -1) {
            const key = (j % P) + ':' + f.type + ':' + (f.label || '');
            if (!placed.has(key)) { i = j; rev = r2; break; }
            j = hay.indexOf(r2 ? B.revcomp(sub) : sub, j + 1);
          }
        }
        if (i >= P) i -= P;
        if (i === -1 || i >= P) continue;
        // map each piece into product coordinates (split at the origin if needed)
        let mapped = blk.rel.map(r => {
          const st = rev ? i + (sub.length - r.off - r.len) : i + r.off;
          return { start: st, end: st + r.len, strand: rev ? -(r.strand || 1) : (r.strand || 1) };
        });
        // ordinary parts are stored in top-strand order, which reverses on the
        // other strand; mixed parts are in biological order, which does not
        if (rev && !f.mixed) mapped.reverse();
        const pieces = [];
        for (const m of mapped) {
          const a = m.start % P, b = a + (m.end - m.start);
          if (b <= P) pieces.push({ start: a, end: b, strand: m.strand });
          else { pieces.push({ start: a, end: P, strand: m.strand }); pieces.push({ start: 0, end: b - P, strand: m.strand }); }
        }
        const base = Object.assign({}, f, { strand: rev ? -(f.strand || 1) : f.strand });
        let g;
        if (f.mixed) {
          // biological order: reverse-complementing flips the order too
          g = B.assembleFeature(base, pieces, P);
        } else {
          g = B.assembleFeature(base, pieces.map(p => ({ start: p.start, end: p.end })), P);
        }
        if (!g) continue;
        // never keep a feature whose bases do not match the source exactly
        if (B.featureSequence(productSeq, g) !== B.featureSequence(src.seq, f)) continue;
        placed.add((i % P) + ':' + f.type + ':' + (f.label || ''));
        out.push(PB.normFeature(g));
      }
    }
    const seen = new Set(), dedup = [];
    for (const f of out) {
      const k = (B.formatLocation(f, P) || f.start + ':' + f.end) + ':' + f.type + ':' + f.label;
      if (seen.has(k)) continue;
      seen.add(k); dedup.push(f);
    }
    return dedup;
  }
  PB.transferFeatures = transferFeatures;

  /* A coding region to recode: its sequence in reading direction, the frame
     offset (/codon_start), the genetic code, and how to write a same-length
     edit back into the construct. Joined and origin-wrapping features are
     written base by base through their own geometry.                      */
  function codingSource(rec, which) {
    if (which === '__sel' && state.selection) {
      const a = state.selection.start, b = state.selection.end;
      return { cds: rec.seq.slice(a, b), frame: 0, table: state.codeTable, isFeature: false,
               apply: s => rec.seq.slice(0, a) + s + rec.seq.slice(a + s.length) };
    }
    if (which === '__all') return { cds: rec.seq, frame: 0, table: state.codeTable, isFeature: false, apply: s => s };
    const f = rec.features[+which];
    if (!f) return null;
    const q = f.quals || {};
    const cs = parseInt(q.codon_start, 10);
    return {
      cds: B.featureSequence(rec.seq, f), feature: f, isFeature: true,
      frame: cs >= 1 && cs <= 3 ? cs - 1 : 0,
      table: parseInt(q.transl_table, 10) || state.codeTable,
      apply: s => PB.writeFeatureSequence(rec, f, s)
    };
  }
  PB.codingSource = codingSource;
  // After a recoding edit: every coding feature must still give the same
  // protein, and site counts for the named enzymes before and after.
  function verifyRecode(rec, newSeq, enzymeNames) {
    const circ = rec.topology === 'circular';
    const proteinChanged = [];
    rec.features.forEach(f => {
      if (!/^CDS$/i.test(f.type) || f.between) return;
      const q = f.quals || {}, cs = parseInt(q.codon_start, 10), fr = cs >= 1 && cs <= 3 ? cs - 1 : 0;
      const t = parseInt(q.transl_table, 10) || state.codeTable;
      const a = B.translate(B.featureSequence(rec.seq, f).slice(fr), { table: t });
      const b = B.translate(B.featureSequence(newSeq, f).slice(fr), { table: t });
      if (a !== b) proteinChanged.push(f.label || 'a CDS');
    });
    const sites = {};
    (enzymeNames || []).forEach(n => {
      const e = PB.enzByName(n);
      if (e) sites[e.name] = sites[n] = { before: B.findSites(rec.seq, e, circ).length, after: B.findSites(newSeq, e, circ).length };
    });
    return { proteinChanged, sites };
  }
  PB.verifyRecode = verifyRecode;

  function saveProduct(o) {
    const nr = PB.makeRecord({
      name: o.name, seq: o.seq,
      topology: o.circular ? 'circular' : 'linear',
      features: o.features || [],
      description: o.description || ''
    });
    (o.history || []).forEach(h => PB.addHistory(nr, h[0], h[1]));
    PB.addRecord(nr);
    PB.closeModal();
    PB.toast('Created “' + o.name + '” · ' + fmt(nr.seq.length) + ' bp', 'ok');
    return nr;
  }
  PB.saveProduct = saveProduct;

  /* What an in-place swap does to the reading frame.

     Replacing a stretch in the middle of a gene - a solubility tag, a
     domain - is the one edit where the plasmid map can look perfect and
     the protein still be ruined: one base out and everything behind the
     join is read in the wrong frame, with the protease site gone and the
     protein nonsense. So the result is translated both ways and compared,
     rather than being taken on trust. */
  function frameBox(check) {
    if (!check) return '';
    if (!check.applies) {
      return check.why && /Nothing coding/.test(check.why) ? '' :
        '<div class="note">Reading frame not checked: ' + esc(check.why || '') + '</div>';
    }
    const j = check.join;
    const arrow = ' <span style="color:var(--accent);font-weight:700">|</span> ';
    const faint = t => '<span style="color:var(--fg-faint)">' + esc(t) + '</span>';
    const line = (label, cut) => '<div class="mono" style="font-size:12px;margin:2px 0">' +
      '<span style="display:inline-block;min-width:86px;color:var(--muted)">' + label + '</span>' +
      faint(cut.left) + arrow + esc(cut.right) + '</div>';
    const head = check.ok
      ? '<strong>In frame.</strong> The protein goes from ' + fmt(check.residuesBefore) + ' to ' +
        fmt(check.residuesAfter) + ' residues' +
        (check.protease && check.protease.kept
          ? ', and the ' + esc(check.protease.name) + ' site behind it survives' : '') + '.'
      : '<strong>This breaks the protein.</strong>';
    return '<div class="note ' + (check.ok ? 'ok' : 'bad') + '">' + head +
      (check.notes.length
        ? '<ul style="margin:6px 0 0 16px;padding:0">' +
          check.notes.map(n => '<li>' + esc(n) + '</li>').join('') + '</ul>'
        : '') + '</div>' +
      '<div class="note" style="padding-top:8px">' +
      '<div style="color:var(--muted);font-size:11.5px;margin-bottom:4px">Where the new piece joins what was already there, read as protein</div>' +
      line('was', j.before.end) +
      line('now', j.after.end) +
      (j.insertResidues
        ? '<div class="hint" style="margin-top:4px">The piece going in is ' + fmt(j.insertResidues) +
          ' residues; only its far end is shown here.</div>'
        : '') +
      '</div>';
  }

  function junctionBox(label, text) {
    return '<div class="field wide" style="margin-top:8px"><label>' + esc(label) + '</label>' +
      '<div class="seqblock">' + text + '</div></div>';
  }

  /* ============================== TOOLS MENU ============================== */
  const TOOLS = [
    ['Cloning', 'Restriction cloning', 'Digest a vector and insert, check the junctions, ligate.', () => toolRestriction()],

    ['Cloning', 'Seamless assembly', 'In-Fusion, Gibson or NEBuilder HiFi \u2014 the same design, with the arm length set by the kit.', () => toolGibson('infusion')],

    ['Cloning', 'Gateway BP / LR', 'att-site recombination into entry and destination vectors.', () => toolGateway()],
    ['Cloning', 'TOPO / TA / GC', 'Topoisomerase and T-tailed vector cloning.', () => toolTopo()],
    ['PCR', 'PCR simulation', 'Amplify with any primer pair, tails included.', () => toolPCR()],
    ['PCR', 'Site-directed mutagenesis', 'Substitute, insert or delete, with primers.', () => toolSDM()],
    ['PCR', 'Overlap extension', 'Fuse two fragments through a shared overlap.', () => toolOverlap()],
    ['Design', 'Silent mutations', 'Add or remove a restriction site without changing the protein.', () => toolSilent()],
    ['Design', 'Codon optimisation', 'Re-code a CDS for a host’s codon bias.', () => toolCodon()],
    ['Design', 'Back-translate', 'Turn a protein into a host-optimised coding sequence.', () => toolBackTranslate()],


    ['Design', 'Intron splicing', 'Strip introns to make the processed mRNA.', () => toolIntrons()],
    ['Design', 'Translate', 'Translate the selection in any frame or genetic code.', () => PB.toolTranslate()],
    ['Edit', 'Edit sequence', 'Replace the whole sequence, or paste a new one.', () => toolEditSeq()],
    ['Edit', 'Edit sequence', 'Insert, delete or replace bases in place, features moving with the edit. Undoable.', () => toolInsertDelete()],
    ['Edit', 'Set origin', 'Rotate a circular construct to a new base 1.', () => PB.setOriginDialog()],
    ['Edit', 'Reverse complement', 'Flip the whole construct and its features.', () => PB.flipDialog()],
    ['Edit', 'Topology', 'Switch between circular and linear.', () => toolTopology()],
    ['Design', 'Read a fragment', 'What in this sequence is protein, and what was added around it.', () => toolAnatomy()],
    ['Edit', 'Demo construct', 'Hide the built-in example, or bring it back.', () => PB.demoDialog()],
    ['Expression', 'Expression strategy', 'Read a sequence for how hard it will be to get soluble, and pick the construct from what you have.', () => PB.toolExpression()],
    ['Expression', 'Competent cells', 'Import the strain list, so the plan knows what each host has already spent.', () => PB.cellImportDialog()],
    ['Expression', 'Plasmid stock', 'Import the lab\u2019s plasmid collection once, and plan against it from then on.', () => PB.plasmidImportDialog()],
    ['Sequencing', 'Sanger primers from stock', 'Search your whole primer inventory for the tubes that can read the selected region, on both strands.', () => PB.toolSanger()],
    ['Sequencing', 'Primer stock', 'Import the lab\u2019s oligo inventory once, and search it from then on.', () => PB.stockImportDialog()],
    ['Bench', 'Reaction setup', 'Master mixes, molar ratios, cycling and conversions \u2014 filled in from the last design.', () => toolBench()]
  ];

  PB.openTools = function () {
    const groups = {};
    TOOLS.forEach((t, i) => { (groups[t[0]] = groups[t[0]] || []).push([t, i]); });
    let body = '';
    Object.keys(groups).forEach(g => {
      body += '<h4 style="font-size:10.5px;text-transform:uppercase;letter-spacing:.08em;color:var(--fg-faint);margin:14px 0 7px">' +
        esc(g) + '</h4><div class="tool-grid">' +
        groups[g].map(([t, i]) => '<button class="tool-card" data-tool="' + i + '">' +
          '<span class="tn">' + esc(t[1]) + '</span><span class="td">' + esc(t[2]) + '</span></button>').join('') +
        '</div>';
    });
    const w = PB.modal({ title: 'Tools', sub: TOOLS.length + ' operations', wide: true, body: body });
    w.querySelectorAll('[data-tool]').forEach(b => b.addEventListener('click', () => {
      const i = parseInt(b.getAttribute('data-tool'), 10);
      PB.closeModal();
      TOOLS[i][3]();
    }));
  };

  /* ============================== RESTRICTION CLONING ============================== */
  function toolRestriction() {
    const rec = PB.active();
    PB.modal({
      title: 'Restriction cloning', wide: true,
      sub: 'digest, check ends, ligate',
      body: '<div class="row">' +
          seqField('rcVec', 'Vector (backbone)', { defaultId: rec ? rec.id : '', allowSel: false }) +
          seqField('rcIns', 'Insert (source of the fragment)', { allowSel: true }) +
        '</div><div class="row">' +
          enzField('rcVecEnz', 'Vector enzymes', (state.enzymes.slice(0, 2).join(' ') || 'EcoRI HindIII'), true) +
          enzField('rcInsEnz', 'Insert enzymes', (state.enzymes.slice(0, 2).join(' ') || 'EcoRI HindIII'), true) +
          '<div class="field" style="flex:0 0 auto;align-self:flex-end"><button class="btn primary" id="rcGo">Simulate</button></div>' +
        '</div>' +
        '<div class="hint">The backbone defaults to the largest vector fragment and the insert to the largest insert fragment. Compatible ends are matched by overhang sequence, so BamHI/BglII or SalI/XhoI pairs work as they do at the bench.</div>' +
        resultSlot(),
      onOpen: w => {
        w.querySelector('#rcGo').addEventListener('click', () => {
          const v = seqValue('rcVec'), ins = seqValue('rcIns');
          if (!v.seq.length || !ins.seq.length) { PB.toast('Provide both a vector and an insert.', 'warn'); return; }
          const ve = enzNames('rcVecEnz'), ie = enzNames('rcInsEnz');
          if (!ve.length || !ie.length) { PB.toast('Name at least one enzyme for each.', 'warn'); return; }
          const vSrc = v.rec || v.from, iSrc = ins.rec || ins.from;
          const vm = vSrc ? PB.methylSystems(vSrc) : {}, im = iSrc ? PB.methylSystems(iSrc) : {};
          const vd = B.digestWithEnds(v.seq, ve, { circular: v.circular, methylation: vm });
          const idg = B.digestWithEnds(ins.seq, ie, { circular: ins.circular, methylation: im });
          let html = '<div class="row" style="margin-top:12px">' +
            fragTable('Vector fragments', vd.fragments, 'rcv') +
            fragTable('Insert fragments', idg.fragments, 'rci') + '</div>' +
            methylNote('vector', v, vd) + methylNote('insert', ins, idg);
          const r = B.restrictionClone({
            vector: v.seq, insert: ins.seq, vectorEnzymes: ve, insertEnzymes: ie,
            vectorCircular: v.circular, insertCircular: ins.circular,
            vectorMethylation: vm, insertMethylation: im
          });
          if (!r.ok) {
            html += '<div class="note bad"><strong>' + esc(r.error) + '</strong>' +
              (r.detail ? '<br>Vector backbone ends: ' + endStr(r.detail.vectorLeft) + ' … ' + endStr(r.detail.vectorRight) +
                '<br>Insert ends: ' + endStr(r.detail.insertLeft) + ' … ' + endStr(r.detail.insertRight) +
                '<br>Overhangs must be reverse complements of one another, and both ends must be the same type.' : '') +
              '</div>';
            result(html);
            return;
          }
          html += '<div class="note ok"><strong>Ligation works.</strong> ' +
            r.products.length + ' product' + (r.products.length === 1 ? '' : 's') +
            ' (' + r.products.map(p => p.orientation).join(' and ') + ').</div>';
          r.products.forEach((p, i) => {
            html += '<div class="panel" style="margin:10px 0 0"><h3>' +
              (p.orientation === 'forward' ? 'Forward orientation' : 'Reverse orientation') +
              ' <span class="sub">' + fmt(p.length) + ' bp circular</span>' +
              '<button class="btn sm primary" data-rcsave="' + i + '">Save construct</button></h3>' +
              '<dl class="kv"><dt>5′ junction</dt><dd>' + esc(p.junction5) + '</dd>' +
              '<dt>3′ junction</dt><dd>' + esc(p.junction3) + '</dd>' +
              '<dt>Backbone</dt><dd>' + fmt(r.vectorFragment.length) + ' bp</dd>' +
              '<dt>Insert</dt><dd>' + fmt(r.insertFragment.length) + ' bp</dd></dl></div>';
          });
          result(html);
          $('#toolResult').querySelectorAll('[data-rcsave]').forEach(b => b.addEventListener('click', () => {
            const p = r.products[+b.getAttribute('data-rcsave')];
            const sources = [v.rec || v.from, ins.rec || ins.from].filter(Boolean);
            saveProduct({
              name: (v.rec ? v.rec.name : 'vector') + ' + insert (' + fmt(p.length) + ' bp)',
              seq: p.seq, circular: true,
              features: transferFeatures(p.seq, sources, true),
              history: [['Restriction cloning',
                ve.join('+') + ' vector × ' + ie.join('+') + ' insert, ' + p.orientation +
                ' orientation; junctions ' + p.junction5 + ' and ' + p.junction3]]
            });
          }));
        });
      }
    });
  }

  // Methylation consequences of a digest, in words, for the cloning tools
  function methylNote(what, src, d) {
    const sysName = k => ({ dam: 'Dam', dcm: 'Dcm', cpg: 'CpG' })[k];
    const owner = src.rec || src.from;
    const mode = owner ? PB.METHYL_MODES[PB.methylMode(owner)].label.split(' (')[0] : 'unmethylated (pasted sequence)';
    const parts = [];
    if (d.blocked && d.blocked.length) parts.push('<strong>' + d.blocked.map(b => b.enzyme + ' at ' + fmt(b.pos + 1) +
      ' (' + b.by.map(sysName).join('/') + ')').join(', ') + ' will not cut</strong> — treated as uncut above. ' +
      (d.blocked.some(b => b.by.indexOf('cpg') < 0) ? 'Grow the plasmid in a dam⁻/dcm⁻ strain (e.g. NEB C2925) or choose another enzyme.' : ''));
    const maybe = (d.caution || []).filter(c => c.effect !== 'unknown');
    if (maybe.length) parts.push(maybe.map(b => b.enzyme + ' at ' + fmt(b.pos + 1)).join(', ') +
      ' may cut incompletely (' + Array.from(new Set([].concat.apply([], maybe.map(m => m.by)))).map(sysName).join('/') + ' overlap).');
    const unk = (d.caution || []).filter(c => c.effect === 'unknown');
    if (unk.length) parts.push(unk.map(b => b.enzyme).filter((x, i, a) => a.indexOf(x) === i).join(', ') +
      ': site overlaps methylation and this enzyme is not in the NEB chart — check the supplier’s data.');
    if (!parts.length) return '';
    return '<div class="note ' + (d.blocked && d.blocked.length ? 'bad' : 'warn') + '"><strong>Methylation, ' + what + '</strong> (' + esc(mode) + '): ' + parts.join(' ') + '</div>';
  }
  PB.methylNote = methylNote;

  function endStr(e) {
    if (!e) return '?';
    return e.type === 'blunt' ? 'blunt' : e.type + ' ' + e.seq;
  }
  function fragTable(title, frags, key) {
    return '<div class="field" style="flex:1 1 320px"><label>' + esc(title) + ' (' + frags.length + ')</label>' +
      '<div class="scrollx" style="margin:0"><table class="tbl"><thead><tr><th class="num">bp</th><th>Left end</th><th>Right end</th><th>Cut by</th></tr></thead><tbody>' +
      frags.slice().sort((a, b) => b.length - a.length).slice(0, 8).map(f =>
        '<tr><td class="num">' + fmt(f.length) + '</td>' +
        '<td class="mono" style="font-size:11px">' + esc(endStr(f.left)) + '</td>' +
        '<td class="mono" style="font-size:11px">' + esc(endStr(f.right)) + '</td>' +
        '<td style="font-size:11px">' + esc([f.leftEnzyme, f.rightEnzyme].filter(x => x && x !== 'end').join(' / ') || '—') + '</td></tr>').join('') +
      '</tbody></table></div></div>';
  }

  /* ============================== GIBSON / HiFi / IN-FUSION ============================== */
  /* ===================== SEAMLESS ASSEMBLY CLONING =====================
     In-Fusion, Gibson and HiFi are the same design problem: put an insert
     into a backbone and let terminal homology do the joining. What a cloner
     actually decides is three things — what goes in, what it goes into, and
     where. Everything after that is arithmetic, so nothing here waits for a
     button: pick the three and the opened vector, both primer pairs, the
     junctions, the reading frame and the finished plasmid appear at once and
     follow every change.
  ==================================================================== */
  function toolGibson(method) {
    const isInf = method === 'infusion';
    const rec = PB.active();
    const defOv = isInf ? 15 : 20;
    PB.modal({
      title: 'Seamless assembly', wide: true,
      sub: isInf ? 'Takara seamless cloning, 15 bp homology arms' : 'exonuclease-based seamless assembly',
      body: '<div class="chipbar"><button class="chip on" data-gm="design">Design primers</button>' +
          '<button class="chip" data-gm="assemble">Assemble existing fragments</button>' +
          '<span style="flex:1 1 auto"></span>' +
          '<span class="hint" style="align-self:center">Kit</span>' +
          '<select class="inp" id="gbKit" style="flex:0 0 auto;width:auto;min-width:230px">' +
            '<option value="infusion"' + (isInf ? ' selected' : '') + '>In-Fusion \u2014 15 bp arms</option>' +
            '<option value="gibson"' + (isInf ? '' : ' selected') + '>Gibson / NEBuilder HiFi \u2014 20 bp arms</option>' +
          '</select></div>' +
        '<div id="gbDesign">' +
          '<div class="row" style="align-items:flex-start">' +
            seqField('gbIns', 'Insert — what goes in', { allowSel: true, rows: 2, flex: '260px' }) +
            '<div class="field" style="flex:0 1 250px"><label>Which part of it</label>' +
              '<select class="inp" id="gbInsPart"></select>' +
              '<input class="inp" id="gbInsName" placeholder="name it, e.g. CDK1" spellcheck="false" style="margin-top:5px">' +
              '<div class="row" style="gap:4px;margin:5px 0 0;flex-wrap:nowrap;align-items:stretch">' +
                '<input class="inp mono" id="gbInsFrom" type="number" min="1" placeholder="from" style="flex:1 1 0;min-width:0">' +
                '<input class="inp mono" id="gbInsTo" type="number" min="1" placeholder="to" style="flex:1 1 0;min-width:0">' +
                '<button class="chip" id="gbInsUnit" title="count in bases or in amino acids" style="flex:0 0 auto">bp</button>' +
              '</div>' +
              '<div class="hint" id="gbInsSpan" style="margin:3px 0 0"></div></div>' +
            seqField('gbVec', 'Vector — what it goes into', { defaultId: rec ? rec.id : '', allowSel: true, rows: 2, flex: '260px' }) +
          '</div>' +
          '<div class="row" style="align-items:flex-end">' +
            '<div class="field" style="flex:1 1 420px"><label id="gbWhereLabel">Where it goes</label>' +
              '<select class="inp" id="gbWhere"></select></div>' +
            '<div class="field narrow" id="gbFromWrap" hidden><label>From (bp)</label><input class="inp mono" id="gbFrom" type="number" min="1"></div>' +
            '<div class="field narrow" id="gbToWrap" hidden><label>To (bp)</label><input class="inp mono" id="gbTo" type="number" min="1" placeholder="blank"></div>' +
            '<div class="field narrow"><label>Arms (bp)</label><input class="inp mono" id="gbOv" type="number" value="' + defOv + '" min="15" max="60"></div>' +
            '<div class="field narrow"><label>Anneal Tm (°C)</label><input class="inp mono" id="gbAn" type="number" value="60" min="50" max="72" step="1"></div>' +
          '</div>' +
          '<div id="gbLive"></div>' +
        '</div>' +
        '<div id="gbAssemble" hidden>' +
          '<div class="hint">Every fragment here is treated as linear, exactly as written. A circular backbone must be opened first — use the <em>Design primers</em> tab, or paste it already linearised.</div>' +
          '<div class="row"><div class="field wide"><label>Fragments in FASTA (each must share terminal homology with the next)</label>' +
          '<textarea class="inp" id="gbFrags" rows="8" spellcheck="false" placeholder="&gt;backbone&#10;ACGT…&#10;&gt;insert1&#10;ACGT…"></textarea></div></div>' +
          '<div class="row"><div class="field narrow"><label>Min overlap</label><input class="inp mono" id="gbMinOv" type="number" value="15"></div>' +
          '<div class="field narrow" style="flex:0 0 130px;align-self:flex-end"><label class="chip on" id="gbCirc">Circularise</label></div>' +
          '<div class="field" style="flex:0 0 auto;align-self:flex-end"><button class="btn primary" id="gbAsmGo">Assemble</button></div>' +
          '<div class="field" style="flex:0 0 auto;align-self:flex-end"><button class="btn sm" id="gbFromLib">Add from library</button></div>' +
          '</div>' + resultSlot() +
        '</div>',
      onOpen: w => {
        /* Amplifying the vector on its own is a real job, so it is offered
           outright rather than left to be discovered by clearing the box.
           It becomes the default only when the insert would otherwise be the
           same record as the vector, which is never a design anyone means. */
        const insSel = w.querySelector('#gbIns_src');
        insSel.insertAdjacentHTML('afterbegin',
          '<option value="__none">nothing \u2014 just amplify the vector</option>');
        let insChosen = false;
        if (insSel.value === (w.querySelector('#gbVec_src') || {}).value)
          insSel.value = '__none';
        /* handed over from "Read a fragment": the coding part, already trimmed */
        if (state.gbPrefillInsert) {
          insSel.value = ''; insChosen = true;
          w.querySelector('#gbIns_txt').value = state.gbPrefillInsert.seq;
          const nm = w.querySelector('#gbInsName');
          if (nm) { nm.value = state.gbPrefillInsert.name; nm.dataset.touched = '1'; }
          state.gbPrefillInsert = null;
        }

        let circ = true;
        w.querySelectorAll('[data-gm]').forEach(c => c.addEventListener('click', () => {
          w.querySelectorAll('[data-gm]').forEach(x => x.classList.toggle('on', x === c));
          const d = c.getAttribute('data-gm') === 'design';
          $('#gbDesign').hidden = !d; $('#gbAssemble').hidden = d;
          if (d) run();
        }));
        const cEl = w.querySelector('#gbCirc');
        cEl.addEventListener('click', () => { circ = !circ; cEl.classList.toggle('on', circ); });
        /* the kit decides the arm length and the wording, and nothing else:
           the design is the same homology-arm design either way */
        const kitSel = w.querySelector('#gbKit');
        if (kitSel) kitSel.addEventListener('change', () => {
          method = kitSel.value;
          const ovEl = $('#gbOv');
          if (ovEl) ovEl.value = method === 'infusion' ? 15 : 20;
          if (typeof run === 'function') run();
        });
        w.querySelector('#gbFromLib').addEventListener('click', () => {
          const ta = $('#gbFrags');
          ta.value = state.records.slice(0, 4).map(r => '>' + r.name.replace(/\s+/g, '_') + '\n' + r.seq).join('\n');
        });

        /* ---------- what goes in ----------
           A record usually holds more than the thing being cloned, so its
           coding features are offered by name: picking "CDK1" clones the
           gene, not the plasmid it currently sits in. */
        function insertParts() {
          /* "nothing" means nothing only while nothing has been pasted: text
             in the box below is an insert, whatever the dropdown still says */
          const pasted = ($('#gbIns_txt') ? $('#gbIns_txt').value : '').trim();
          if ($('#gbIns_src') && $('#gbIns_src').value === '__none' && !pasted)
            return [{ label: 'nothing \u2014 vector only', seq: '', name: 'vector only' }];
          const v = seqValue('gbIns');
          const r = v.rec || v.from;
          const out = [{ label: 'the whole sequence (' + fmt(v.seq.length) + ' bp)', seq: v.seq, name: v.name }];
          if (r && r.features && v.rec && $('#gbIns_src').value !== '__none') {
            r.features.forEach(f => {
              if (f.between || !/^(CDS|gene|ORF|misc_feature)$/i.test(f.type || '')) return;
              const len = B.featureLength(f, r.seq.length);
              if (len < 30) return;
              const s = B.featureSequence(r.seq, f);
              out.push({ label: (f.label || f.type) + ' · ' + f.type + ', ' + fmt(len) + ' bp' +
                           (f.strand === -1 ? ', reverse strand' : ''),
                         seq: s, name: f.label || f.type, feature: f });
            });
          }
          return out;
        }
        let parts = [];
        function refreshParts(keepName) {
          const sel = $('#gbInsPart');
          if (!sel) return;
          parts = insertParts();
          const prev = sel.value;
          sel.innerHTML = parts.map((p, i) => '<option value="' + i + '">' + esc(p.label) + '</option>').join('');
          if (prev && parts[+prev]) sel.value = prev;
          /* a coding feature is almost always what someone means */
          else { const cds = parts.findIndex(p => p.feature && /^CDS$/i.test(p.feature.type)); if (cds > 0) sel.value = String(cds); }
          if (!keepName) autoName();
        }
        /* Only part of a gene is often what is wanted: a kinase domain, a
           construct truncated before its C-terminal tail. The range trims
           whatever part was picked, counted in bases or in residues — a
           domain is remembered in residues, and converting by hand is where
           an off-by-one creeps in. */
        function insertPart() {
          const base = parts[+$('#gbInsPart').value] || parts[0] || { seq: '', name: 'insert' };
          const L = base.seq.length;
          const fEl = $('#gbInsFrom'), tEl = $('#gbInsTo');
          const fv = fEl ? (fEl.value || '').trim() : '', tv = tEl ? (tEl.value || '').trim() : '';
          if (!fv && !tv) return Object.assign({}, base, { from: 0, to: L, whole: L, trimmed: false });
          const aa = (state.gbUnit || 'bp') === 'aa';
          let a, b;
          if (aa) {
            const f = Math.max(1, parseInt(fv, 10) || 1);
            const t = tv ? Math.max(f, parseInt(tv, 10) || f) : Math.floor(L / 3);
            a = (f - 1) * 3; b = Math.min(L, t * 3);
          } else {
            a = Math.max(0, (parseInt(fv, 10) || 1) - 1);
            b = tv ? Math.min(L, Math.max(a + 1, parseInt(tv, 10) || L)) : L;
          }
          if (!(b > a) || a >= L) return Object.assign({}, base, { from: 0, to: L, whole: L, trimmed: false, bad: true });
          return Object.assign({}, base, { seq: base.seq.slice(a, b), from: a, to: b, whole: L, trimmed: a > 0 || b < L });
        }
        function autoName() {
          const n = $('#gbInsName');
          if (n && !n.dataset.touched) n.value = (insertPart().name || 'insert').replace(/\s+/g, ' ').trim();
        }
        $('#gbInsName').addEventListener('input', e => { e.target.dataset.touched = '1'; });

        /* ---------- where it goes ----------
           One list of real choices, built from the vector itself: the
           placements its own annotation implies, every enzyme that cuts it
           once, whatever is selected on the map, and a free range. */
        let places = [], placeCache = '', whereChosen = false;
        /* A primer keeps the name it is ordered under, so the auto-generated
           one is a starting point, not a decision. Names are held by role,
           so they survive every recompute; clearing a box restores the
           automatic name. */
        const primerNames = {};
        /* lengths chosen by hand in the tuning panel; cleared when the
           placement changes, since they refer to that opening */
        let tuned = {};
        function primerName(role, auto) {
          const v = primerNames[role];
          return (v && v.trim()) ? v.trim() : auto;
        }
        function nameCell(role, name) {
          return '<td style="padding:2px 6px"><input class="inp mono pnm" data-role="' + esc(role) + '" ' +
            'value="' + esc(name) + '" spellcheck="false" aria-label="Name for ' + esc(role) + '" ' +
            'style="font-size:11.5px;padding:3px 6px;min-width:120px;width:100%"></td>';
        }
        function wireNames(host, rows) {
          host.querySelectorAll('.pnm').forEach(el => {
            const role = el.getAttribute('data-role');
            el.addEventListener('input', () => {
              primerNames[role] = el.value;
              const row = rows.find(x => x[0] === role);
              if (row) row[1] = el.value.trim() || row[1];
            });
            el.addEventListener('change', () => {
              if (!el.value.trim()) { delete primerNames[role]; run(); }
            });
          });
        }
        function placements(rec, seq, circular) {
          const out = [];
          if (!circular) return out;
          cloningSites(rec).forEach(x => out.push({ kind: 'ipcr', group: 'Into a cloning site', label: x.label, ex: x.ex, from: x.from, to: x.to }));
          /* Swapping one gene for another is the commonest job of all: take
             the construct you have, amplify everything except the insert
             currently in it — tags and backbone kept — and put the next gene
             in its place. Each coding feature big enough to be an insert is
             offered by name. */
          if (rec && rec.features) {
            const L0 = seq.length;
            const swaps = [];
            rec.features.forEach(f => {
              if (f.between || !/^(CDS|gene|ORF)$/i.test(f.type || '')) return;
              const len = B.featureLength(f, L0);
              if (len < 60 || len > L0 - 200) return;
              const pieces = B.featurePieces(f, L0);
              if (pieces.length > 1) return;                  // a spliced feature is not a cassette
              const ess = essentialPart(f);
              swaps.push({ kind: 'ipcr', a0: f.start, b0: f.end, group: 'Swap out a gene it already carries',
                essential: ess, label: 'replace ' + (f.label || f.type) + ' (' + f.type + ', ' + fmt(len) + ' bp)' +
                  (ess ? ' \u2014 ' + ess + ', the plasmid needs this' : ' \u2014 keeps everything else'),
                ex: 'everything except ' + (f.label || f.type) + ' (' + f.type + ', ' + fmt(len) + ' bp)' +
                  (ess ? ' \u2014 ' + ess + ', the plasmid needs this' : '') });
            });
            /* what the plasmid needs to survive goes last: nobody opens this
               list meaning to delete the resistance marker */
            swaps.sort((a, b) => (a.essential ? 1 : 0) - (b.essential ? 1 : 0));
            swaps.forEach(x => out.push(x));
          }
          const meth = rec ? B.methylatedBases(seq, true, PB.methylSystems(rec)) : null;
          const cutters = [];
          PB.allEnzymes().forEach(e => {
            if (!PB.COMMON.has(e.name) && !e.custom) return;
            let sites = [];
            try { sites = B.findSites(seq, e, true); } catch (err) { return; }
            if (sites.length !== 1) return;
            const m = meth ? B.methylationEffect(e, sites[0], seq.length, meth) : null;
            if (m && m.effect === 'blocked') return;
            cutters.push({ kind: 'enz', name: e.name, pos: sites[0].cutTop, group: 'Cut once with',
                           label: 'cut with ' + e.name + ' (single site, cuts at ' + fmt(sites[0].cutTop) + ')' +
                                  (m && m.effect ? ' ⚠ ' + m.effect : '') });
          });
          cutters.sort((a, b) => a.pos - b.pos).forEach(c => out.push(c));
          const sel = PB.normSel(rec, state.selection);
          if (rec && sel && PB.active() && PB.active().id === rec.id) {
            const selWhere = PB.selRangeText(rec, state.selection) + ', ' + fmt(sel.length) + ' bp';
            /* a0/b0 are the half-open bounds read round the circle, which is
               how inverse PCR is asked to remove a region. They matter for a
               selection that crosses base 1: given only from and to, a
               wrapping range arrives with to below from, which reads as an
               empty range — the opening was then made with nothing removed
               and the gene the selection was pointing at stayed in the
               product. The gene a person selects to replace is very often
               the one that wraps, because that is where the last assembly
               put it. */
            out.push({ kind: 'ipcr', group: 'Somewhere you choose',
              label: 'replace what is selected on the map (' + selWhere + ')',
              ex: 'everything except what is selected on the map (' + selWhere + ')',
              selection: true,
              a0: sel.start % seq.length, b0: sel.end % seq.length,
              wholeThing: sel.length >= seq.length,
              from: sel.start + 1, to: (sel.end - 1) % seq.length + 1 });
          }
          out.push({ kind: 'custom', group: 'Somewhere you choose', label: 'a range I type — inverse PCR, removing it',
                     ex: 'everything except a range I type' });
          out.push({ kind: 'pos', group: 'Somewhere you choose', label: 'a single base I type — blunt opening, nothing removed',
                     ex: 'the whole plasmid, opened at a single base I type' });
          return out;
        }
        /* With nothing going in, the same list means the opposite thing: it is
           no longer where an insert lands but what the amplified vector leaves
           out. The choices are identical, so only the wording changes. */
        function vectorOnly() {
          try { return !insertPart().seq.length; } catch (e) { return false; }
        }
        const EX_GROUP = {
          'Into a cloning site': 'Leave out a cloning site',
          'Swap out a gene it already carries': 'Leave out a gene it carries',
          'Somewhere you choose': 'Leave out a region you choose'
        };
        function refreshPlaces() {
          const v = seqValue('gbVec');
          const r = v.rec || v.from;
          const vo = vectorOnly();
          const key = (v.rec ? v.rec.id + ':' + v.rec.updatedAt : 'x') + ':' + v.seq.length + ':' +
                      (v.circular ? 'c' : 'l') + ':' + JSON.stringify(state.selection) + ':' + (vo ? 'v' : 'i');
          const sel = $('#gbWhere');
          if (!sel || key === placeCache) return;
          placeCache = key;
          const lbl = w.querySelector('#gbWhereLabel');
          if (lbl) lbl.textContent = vo ? 'What to leave out' : 'Where it goes';
          const prev = sel.value;
          places = placements(r, v.seq, v.circular);
          if (!v.circular) {
            sel.innerHTML = '<option value="lin">this vector is already linear — use it as it is</option>';
            places = [{ kind: 'linear', label: 'already linear' }];
            sel.value = '0';
            sel.innerHTML = '<option value="0">this vector is already linear — use it as it is</option>';
          } else {
            /* grouped, because a 5 kb plasmid offers thirty single cutters
               and a flat list of them buries the two choices that matter */
            /* an enzyme cut is a digest, not an amplification: with no insert
               there is nothing for it to design, so it is not offered */
            const order = vo
              ? ['Somewhere you choose', 'Swap out a gene it already carries', 'Into a cloning site']
              : ['Into a cloning site', 'Swap out a gene it already carries', 'Somewhere you choose', 'Cut once with'];
            const seen = {};
            places.forEach((pl, i) => { const g = pl.group || 'Other'; (seen[g] = seen[g] || []).push(i); });
            sel.innerHTML = order.filter(g => seen[g]).map(g =>
              '<optgroup label="' + esc(vo ? (EX_GROUP[g] || g) : g) + '">' +
              seen[g].map(i => '<option value="' + i + '">' +
                esc((vo && places[i].ex) || places[i].label) + '</option>').join('') +
              '</optgroup>').join('');
            /* a placement the person picked is kept; one that was only the
               default is picked again, since the list now means something
               else and the old index is not the same choice */
            const keep = whereChosen && prev && places[+prev] &&
                         sel.querySelector('option[value="' + prev + '"]');
            if (keep) sel.value = prev;
            else sel.value = String(defaultPlace(places, vo));
          }
          applyPlace();
        }
        /* Opened cold, the dialog has to pick something, and what it picks is
           what a hurried person will use. So it is the cloning site the vector
           declares — in frame with its tag where there is one — never a gene
           replacement, which is destructive and is only ever chosen on purpose. */
        function defaultPlace(list, vo) {
          const score = p => {
            /* amplifying the vector alone: what is selected on the map is an
               explicit choice already made, so it leads */
            if (vo) {
              if (p.selection) return 0;
              if (p.kind === 'custom') return 1;
              if (p.kind === 'pos') return 2;
              if (p.group === 'Into a cloning site') return 3;
              if (p.kind === 'enz') return 7;
              return p.essential ? 9 : 4;
            }
            if (p.group === 'Into a cloning site') return /in frame/.test(p.label) ? 0 : (/from the ATG/.test(p.label) ? 1 : 2);
            if (p.kind === 'enz') return 5;
            if (p.kind === 'custom' || p.kind === 'pos') return 6;
            return p.essential ? 9 : 8;          // a gene swap is never the default
          };
          let best = 0;
          list.forEach((p, i) => { if (score(p) < score(list[best])) best = i; });
          return best;
        }
        function place() { const s = $('#gbWhere'); return (s && places[+s.value]) || places[defaultPlace(places, vectorOnly())] || null; }
        function applyPlace() {
          if (!$('#gbWhere') || !$('#gbFromWrap')) return;
          const p = place();
          const needFrom = p && (p.kind === 'custom' || p.kind === 'pos');
          $('#gbFromWrap').hidden = !needFrom;
          $('#gbToWrap').hidden = !(p && p.kind === 'custom');
          if (p && p.kind === 'ipcr' && $('#gbFrom').value === '') {
            $('#gbFrom').value = p.from; $('#gbTo').value = p.to === '' ? '' : p.to;
          }
        }

        /* ---------- opening the vector ---------- */
        function openVector(r, seq, circular) {
          if (!circular) return { ok: true, seq: seq, how: 'already linear \u2014 used as supplied, both ends kept',
                                  points: [], primers: null, preLinear: true };
          const p = place();
          if (!p) return { ok: false, error: 'Choose where the insert goes.' };
          const L = seq.length;
          if (p.kind === 'enz') {
            const meth = r ? PB.methylSystems(r) : {};
            const d = B.digestWithEnds(seq, [p.name], { circular: true, methylation: meth });
            if (!d.cuts.length) return { ok: false, error: p.name + ' does not cut this plasmid.' };
            if (d.fragments.length !== 1) return { ok: false, error: p.name + ' cuts ' + d.cuts.length + ' times here — pick a single cutter.' };
            return { ok: true, seq: d.fragments[0].seq,
                     how: 'cut with ' + p.name + ' at ' + fmt(d.cuts[0].top),
                     points: [d.cuts[0].top], primers: null, enzyme: p.name, openAt: d.cuts[0].top };
          }
          if (p.a0 != null) {
            /* A selection of the whole molecule leaves no backbone to
               amplify, and because the bounds meet it would otherwise read
               as removing nothing at all. */
            if (p.wholeThing) return { ok: false,
              error: 'That selection is the whole plasmid — there would be no backbone left to amplify. ' +
                     'Select the part to replace, not everything.' };
            /* bounds given outright; inversePCR reads them round the circle,
               so a cassette that crosses base 1 works too */
            const r0 = B.inversePCR(seq, p.a0, p.b0, { targetTm: 60, forceLen: tuned.bb, conditions: PB.conditions ? PB.conditions() : undefined });
            if (!r0.ok) return { ok: false, error: r0.error };
            const span0 = p.b0 >= p.a0 ? fmt(p.a0 + 1) + '\u2013' + fmt(p.b0)
              : fmt(p.a0 + 1) + '\u2013' + fmt(L) + ' and 1\u2013' + fmt(p.b0);
            return { ok: true, seq: r0.backbone.seq,
                     how: 'inverse PCR, removing ' + fmt(r0.removed) + ' bp (' + span0 + ')',
                     points: [p.a0], removed: { a: p.a0, b: p.b0 }, primers: r0.primers,
                     openAt: p.a0, ipcr: true };
          }
          let from, to;
          if (p.kind === 'pos') { from = clampPos($('#gbFrom').value, L); to = from; }
          else if (p.kind === 'custom') {
            from = clampPos($('#gbFrom').value, L);
            const t = ($('#gbTo').value || '').trim();
            to = t ? Math.max(0, Math.min(L, parseInt(t, 10) || 0)) : from;
            /* Everywhere else in this tool a From above a To runs the range
               across base 1 of a circle. It has to mean the same here, or a
               range typed that way quietly becomes an opening that removes
               nothing. On a linear construct it is simply backwards. */
            if (to < from) {
              if (circular) {
                const r1 = B.inversePCR(seq, from - 1, to, { targetTm: 60, forceLen: tuned.bb,
                  conditions: PB.conditions ? PB.conditions() : undefined });
                if (!r1.ok) return { ok: false, error: r1.error };
                return { ok: true, seq: r1.backbone.seq,
                         how: 'inverse PCR, removing ' + fmt(r1.removed) + ' bp (' +
                              fmt(from) + '–' + fmt(L) + ' and 1–' + fmt(to) + ', across base 1)',
                         points: [from - 1], removed: { a: from - 1, b: to },
                         primers: r1.primers, openAt: from - 1, ipcr: true };
              }
              to = from;
            }
          } else { from = p.from; to = p.to === '' ? from - 1 : p.to; }
          const a0 = from - 1, b0 = to >= from ? to : from - 1;
          if (p.kind === 'pos' || b0 <= a0) {
            const r2 = B.inversePCR(seq, a0, a0, { targetTm: 60, forceLen: tuned.bb, conditions: PB.conditions ? PB.conditions() : undefined });
            if (!r2.ok) return { ok: false, error: r2.error };
            return { ok: true, seq: r2.backbone.seq, how: 'opened before base ' + fmt(from) + ', nothing removed',
                     points: [a0], removed: null, primers: r2.primers, openAt: a0, ipcr: true };
          }
          const r3 = B.inversePCR(seq, a0, b0, { targetTm: 60, forceLen: tuned.bb, conditions: PB.conditions ? PB.conditions() : undefined });
          if (!r3.ok) return { ok: false, error: r3.error };
          return { ok: true, seq: r3.backbone.seq,
                   how: 'inverse PCR, removing ' + fmt(r3.removed) + ' bp (' + fmt(from) + '–' + fmt(to) + ')',
                   points: [a0], removed: { a: a0, b: b0 }, primers: r3.primers, openAt: a0, ipcr: true };
        }
        function clampPos(v, L) { return Math.max(1, Math.min(L, parseInt(v, 10) || 1)); }

        /* ---------- the fusion, read as protein ----------
           What a tagged vector is for: the tag has to run into the insert in
           one frame. The junction is translated and shown, so an off-by-one
           is visible as residues rather than as a number. */
        function fusionReport(r, open, insSeq, ov) {
          if (!r || !r.features || open.openAt == null) return '';
          const L = r.seq.length, out = [];
          r.features.forEach(f => {
            if (!/^CDS$/i.test(f.type) || f.strand === -1 || f.between) return;
            const cs = Math.max(1, Math.min(3, parseInt((f.quals && f.quals.codon_start) || '1', 10) || 1)) - 1;
            const org = f.start + cs;
            const d = ((open.openAt - org) % L + L) % L;
            if (d <= 0 || d > 1200) return;
            const inFrame = d % 3 === 0;
            let peptide = '';
            if (inFrame) {
              const up = open.seq.slice(-Math.min(18, open.seq.length));
              const keep = up.length - (up.length % 3);
              const upAA = B.translate(up.slice(up.length - keep), { table: state.codeTable });
              const dnAA = B.translate(insSeq.slice(0, Math.min(24, insSeq.length - insSeq.length % 3)), { table: state.codeTable });
              peptide = '<div class="seqblock" style="max-height:none;margin-top:4px">' +
                '<span style="color:var(--fg-faint)">' + esc(upAA.split('').join(' ')) + '</span>' +
                ' <span style="color:var(--accent);font-weight:700">|</span> ' +
                '<span style="color:var(--accent);font-weight:600">' + esc(dnAA.split('').join(' ')) + '</span></div>';
            }
            out.push({ name: f.label || 'CDS', d: d, inFrame: inFrame, peptide: peptide });
          });
          if (!out.length) return '';
          /* features upstream in the same frame all say the same thing:
             report the answer once, naming them, and show the junction once */
          const groups = {};
          out.forEach(o => { const k = o.d % 3; (groups[k] = groups[k] || []).push(o); });
          const poiName = ($('#gbInsName').value || 'insert');
          return Object.keys(groups).sort().map(k => {
            const g = groups[k], inFrame = +k === 0;
            const names = g.map(o => esc(o.name)).join(', ');
            return '<div class="note ' + (inFrame ? 'ok' : 'bad') + '"><strong>' +
              names + ' \u2192 ' + esc(poiName) + ': ' +
              (inFrame ? 'in frame' : 'out of frame by ' + k + ' base' + (+k === 1 ? '' : 's')) +
              '.</strong> ' + (inFrame
                ? 'The sequence upstream reads straight into the insert:'
                : 'Shift the opening by ' + (3 - +k) + ' base' + (3 - +k === 1 ? '' : 's') +
                  ', or add ' + (3 - +k) + ' to the forward tail, or the fusion is nonsense.') +
              (inFrame ? g[0].peptide : '') + '</div>';
          }).join('');
        }

        /* What the insert itself carries into the fusion. */
        function insertNotes(insSeq, inFrameFusion, part) {
          const n = [];
          if (part && part.trimmed) {
            const aaFrom = Math.floor(part.from / 3) + 1;
            n.push('Only part of it is being cloned: bases ' + fmt(part.from + 1) + '\u2013' + fmt(part.to) +
              ' of ' + fmt(part.whole) + (part.from % 3 === 0
                ? ', which is residues ' + fmt(aaFrom) + ' onwards of its own reading frame'
                : ', starting mid-codon of its own reading frame'));
          }
          const prot = B.translate(insSeq, { table: state.codeTable });
          const trail = /\*$/.test(prot);
          const internal = prot.replace(/\*$/, '').indexOf('*') >= 0;
          if (insSeq.length % 3 !== 0)
            n.push('The insert is ' + fmt(insSeq.length) + ' bp, not a whole number of codons (' +
              (insSeq.length % 3) + ' over) — anything fused after it is out of frame.');
          if (internal) {
            /* where translation really stops decides what follows: saying
               "no stop at the end, so it reads into the vector" while an
               in-frame stop sits before that end is simply wrong */
            const at = prot.replace(/\*$/, '').indexOf('*');
            const after = insSeq.length - (at + 1) * 3;
            const codonsAfter = Math.floor(after / 3);
            if (codonsAfter <= 6) {
              /* a stop a few codons from the end is how a synthesised
                 fragment is normally ordered: stop, then a cloning site */
              n.push('It ends properly: the protein runs to residue ' + fmt(at) + ', then a stop codon, then ' +
                fmt(after) + ' bp of trailing sequence that is not translated.');
            } else {
              n.push('Translation stops inside it, at residue ' + fmt(at + 1) + ' of ' + fmt(prot.replace(/\*$/, '').length) +
                ' (bases ' + fmt(at * 3 + 1) + '\u2013' + fmt(at * 3 + 3) + '), so the last ' + fmt(after) +
                ' bp are not translated and nothing on the vector downstream is either');
            }
          } else if (trail) n.push('It ends in a stop codon, so a C-terminal tag on the vector will not be translated.');
          else n.push('It has no stop codon at the end, so translation runs on into the vector — fine for a C-terminal tag, not for an untagged protein.');
          if (/^ATG/.test(insSeq) && inFrameFusion)
            n.push('It starts with ATG; after an N-terminal tag that is an extra methionine, which is usual and harmless.');
          return n.length ? '<div class="note"><strong>The insert:</strong> ' + n.map(esc).join(' ') + '</div>' : '';
        }

        /* ---------- what each primer is made of ----------
           The table gives four sequences; this says what each half of each
           one is for. A tailed primer does two jobs at two different times:
           its 3' part sits on the template from the first cycle, and its 5'
           part does nothing until the assembly, where it is the piece that
           finds the neighbouring fragment. Both are spelled out against the
           sequence they pair with, so neither has to be taken on trust. */
        function primerRoles(rows, open, ov, vecSeq, insSeq, v) {
          if (!rows || !rows.length) return '';
          const tr = [];
          rows.forEach(x => {
            const role = x[0], name = x[1], seq = x[2], tail = x[5] || 0;
            const isIns = /Insert/.test(role);
            const fwd = /forward/i.test(role);
            const anneal = seq.slice(tail);
            if (tail) {
              tr.push({ name: name, part: '5′ tail, ' + tail + ' nt', seq: seq.slice(0, tail),
                        pairs: fwd ? 'the backbone’s last ' + tail + ' bases'
                                   : 'the backbone’s first ' + tail + ' bases, reverse complemented',
                        pairSeq: fwd ? vecSeq.slice(-tail) : B.revcomp(vecSeq.slice(0, tail)),
                        when: 'in the assembly, not in the PCR', hot: true, first: true });
            }
            tr.push({ name: tail ? '' : name, part: (tail ? '3′ ' : '') + 'anneals, ' + anneal.length + ' nt',
                      seq: anneal,
                      pairs: isIns ? 'the insert template' : 'the vector, ' +
                        (open && open.primers ? 'reading outward from the opening' : 'as supplied'),
                      pairSeq: isIns
                        ? (fwd ? insSeq.slice(0, anneal.length) : B.revcomp(insSeq.slice(insSeq.length - anneal.length)))
                        : '',
                      when: 'from the first cycle', hot: false, first: !tail });
          });
          return '<div class="field wide" style="margin-top:10px"><label>What each primer is made of</label>' +
            '<div class="scrollx"><table class="tbl"><thead><tr><th>Primer</th><th>Part</th>' +
            '<th>Sequence (5′→3′)</th><th>Pairs with</th><th>Working</th></tr></thead><tbody>' +
            tr.map(r => '<tr' + (r.first ? ' style="border-top:1px solid var(--line)"' : '') + '>' +
              '<td>' + esc(r.name) + '</td><td>' + esc(r.part) + '</td>' +
              '<td class="mono" style="font-size:11.5px' + (r.hot ? ';color:var(--accent);font-weight:600' : '') + '">' +
                esc(r.seq) + '</td>' +
              '<td>' + esc(r.pairs) +
                (r.pairSeq ? '<div class="mono" style="font-size:11px;color:var(--fg-faint)">' + esc(r.pairSeq) + '</div>' : '') +
              '</td>' +
              '<td style="font-size:11.5px' + (r.hot ? ';color:var(--accent)' : '') + '">' + esc(r.when) + '</td></tr>').join('') +
            '</tbody></table></div>' +
            '<div class="hint">Coloured rows are the homology tails. They hang off the template in the first cycles — ' +
            'which is why the annealing temperature is set from the 3′ part alone — and they are what the assembly ' +
            'uses to find the backbone afterwards. The backbone primers have no tail: the backbone’s own ends already ' +
            'are the homology.</div></div>';
        }

        /* ---------- what came on the ends of the fragment ----------
           A synthesised fragment carries whatever its designer put on it.
           When the vector supplies the same thing, including it again is a
           silent duplication in the protein, so it is named and the trim
           is offered as one click. */
        function addonNote(insSeq, vecSeq, part) {
          if (!insSeq || insSeq.length < 30) return '';
          const common = PB.allEnzymes().filter(e => PB.COMMON.has(e.name));
          let a;
          try {
            a = B.insertAddons(insSeq, { vectorUpstream: vecSeq.slice(-33), enzymes: common,
                                         table: state.codeTable });
          } catch (e) { return ''; }
          if (!a || !a.ok) return '';
          const base = part ? part.from : 0;
          const msgs = [];
          a.leading.forEach(x => {
            msgs.push('Bases ' + fmt(base + x.from + 1) + '\u2013' + fmt(base + x.to) + ' encode <strong>' +
              esc(x.aa) + '</strong> \u2014 ' + esc(x.name) +
              (x.inVector ? ', and the vector already supplies one immediately upstream, so cloning these gives the protein two'
                          : ', which was put there when the fragment was designed'));
          });
          if (a.trailing && a.trailing.bp > 0) {
            msgs.push('Bases ' + fmt(base + a.trailing.from + 1) + '\u2013' + fmt(base + a.trailing.to) + ' come after the stop codon' +
              (a.trailing.sites.length ? ' and carry a ' + esc(a.trailing.sites.map(z => z.name).join(' and ')) + ' site' : '') +
              ' \u2014 not translated either way');
          }
          if (!msgs.length) return '';
          const f = base + a.suggestFrom + 1, t = base + a.suggestTo;
          const worth = a.suggestFrom > 0 || (a.trailing && a.trailing.bp >= 6);
          return '<div class="note warn"><strong>What is on the ends of this fragment:</strong> ' +
            msgs.join('. ') + '.' +
            (worth ? ' <button class="btn sm" data-trimins="' + f + ',' + t + '" style="margin-left:6px">' +
              'clone bases ' + fmt(f) + '\u2013' + fmt(t) + ' only</button>' : '') +
            '</div>';
        }

        /* ---------- the alternatives ----------
           One primer is a choice among many, and the choice is a trade: a
           few degrees against a hairpin, a longer arm against a dimer. The
           candidates are laid out with what each costs so the trade is made
           knowingly rather than silently. */
        function flagBadges(fl) {
          if (!fl || !fl.length) return '<span class="badge ok">clean</span>';
          return fl.map(f => '<span class="badge ' + (f.level === 'bad' ? 'bad' : 'warn') + '">' +
            esc(f.text) + '</span>').join(' ');
        }
        function tunePanel(v, open, p, ov, anTm, insSeq, vecSeq) {
          const cond = PB.conditions ? PB.conditions() : undefined;
          let h = '<div class="field wide" style="margin-top:10px"><label>Other primers that would also work</label>';
          /* --- the backbone pair --- */
          if (open.primers && open.openAt != null && open.ipcr) {
            const L = v.seq.length;
            const a0 = open.removed ? open.removed.a : open.openAt;
            const b0 = open.removed ? open.removed.b : open.openAt;
            const tp = B.tunePair(v.seq, { fAt: b0, rAt: a0, targetTm: anTm, circular: true,
                                           conditions: cond, keep: 8, tmBand: 5 });
            if (tp.ok) {
              h += '<div class="hint" style="margin:0 0 4px">Backbone pair — the two ends are fixed by what is being removed, ' +
                'so only the length can move. The product is the same for every row.</div>' +
                '<div class="scrollx"><table class="tbl"><thead><tr><th></th><th class="num">fwd nt</th>' +
                '<th class="num">Tm</th><th class="num">rev nt</th><th class="num">Tm</th><th class="num">ΔTm</th>' +
                '<th>What it costs</th><th></th></tr></thead><tbody>' +
                tp.options.map(o => {
                  const cur = o.forward.length === open.primers.forward.length &&
                              o.reverse.length === open.primers.reverse.length;
                  return '<tr' + (cur ? ' style="background:var(--bg-raise)"' : '') + '><td>' +
                    (cur ? '<span class="badge ok">in use</span>' : '') + '</td>' +
                    '<td class="num mono">' + o.forward.length + '</td><td class="num mono">' + f1(o.forward.tm) + '</td>' +
                    '<td class="num mono">' + o.reverse.length + '</td><td class="num mono">' + f1(o.reverse.tm) + '</td>' +
                    '<td class="num mono">' + f1(o.spread) + '</td>' +
                    '<td>' + flagBadges(o.forward.flags.concat(o.reverse.flags, o.flags)) +
                      (o.outOfBand ? ' <span class="badge bad">far from target Tm</span>' : '') + '</td>' +
                    '<td>' + (cur ? '' : '<button class="btn sm" data-tunebb="' + o.forward.length + ',' + o.reverse.length + '">use</button>') +
                    '</td></tr>';
                }).join('') + '</tbody></table></div>';
            }
          }
          /* --- the insert pair: the arm is free as well as the annealing part --- */
          if (p && insSeq && insSeq.length > 40) {
            const f = B.tuneAssemblyPrimer({ insert: insSeq, neighbour: vecSeq, end: '5', targetTm: anTm,
                                             conditions: cond, armMin: 15, armMax: 30, keep: 6, template: insSeq });
            const r = B.tuneAssemblyPrimer({ insert: insSeq, neighbour: vecSeq, end: '3', targetTm: anTm,
                                             conditions: cond, armMin: 15, armMax: 30, keep: 6, template: insSeq });
            if (f.ok && r.ok) {
              const rowsFor = (t, which, cur) => t.options.map(o => {
                const isCur = o.arm === ov && o.annealLen === cur;
                return '<tr' + (isCur ? ' style="background:var(--bg-raise)"' : '') + '><td>' +
                  (isCur ? '<span class="badge ok">in use</span>' : '') + '</td>' +
                  '<td class="num mono">' + o.arm + '</td><td class="num mono">' + o.annealLen + '</td>' +
                  '<td class="num mono">' + o.length + '</td><td class="num mono">' + f1(o.tmAnneal) + '</td>' +
                  '<td>' + flagBadges(o.flags) + (o.outOfBand ? ' <span class="badge bad">far from target Tm</span>' : '') + '</td>' +
                  '<td>' + (isCur ? '' : '<button class="btn sm" data-tuneins="' + which + ',' + o.arm + ',' + o.annealLen + '">use</button>') +
                  '</td></tr>';
              }).join('');
              h += '<div class="hint" style="margin:10px 0 4px">Insert primers — the arm may be any length that still matches the ' +
                'backbone exactly, so both halves can move. The finished plasmid is the same for every row.</div>' +
                '<div class="scrollx"><table class="tbl"><thead><tr><th></th><th class="num">arm</th>' +
                '<th class="num">anneal nt</th><th class="num">total nt</th><th class="num">Tm anneal</th>' +
                '<th>What it costs</th><th></th></tr></thead><tbody>' +
                '<tr><td colspan="7" style="color:var(--fg-faint);font-size:11px">forward</td></tr>' +
                rowsFor(f, 'f', p.forward.anneal.length) +
                '<tr><td colspan="7" style="color:var(--fg-faint);font-size:11px">reverse</td></tr>' +
                rowsFor(r, 'r', p.reverse.anneal.length) +
                '</tbody></table></div>';
            }
          }
          h += '<div class="hint">Rows are ordered by how much trouble each is likely to give: the melting temperature against ' +
            'the target, hairpins and self-dimers (worse when they involve the 3′ end, which is the end that extends), ' +
            'runs of one base, GC at the 3′ end, and whether it anneals in more than one place. ' +
            'A row marked <em>far from target Tm</em> fixes a structural problem at a cost that is usually not worth paying.</div>';
          return h + '</div>';
        }

        /* ---------- the whole design, recomputed on any change ---------- */
        let lastDesign = null;
        function run() {
          /* the debounce can fire after the dialog is closed */
          const panel = $('#gbDesign'), host = $('#gbLive');
          if (!panel || !host || panel.hidden) return;
          refreshParts(true); refreshPlaces();
          const v = seqValue('gbVec');
          const part = insertPart();
          const insSeq = B.cleanSeq(part.seq || '');
          lastDesign = null;
          if (!v.seq.length) {
            host.innerHTML = '<div class="hint">Pick a vector. With an insert as well you get both primer pairs and the finished plasmid; with the vector alone you get the two primers that amplify it, leaving out whatever <em>Where it goes</em> removes.</div>';
            return;
          }
          /* No insert is a job in its own right: amplify the backbone
             without the gene that is in it now, to get the empty tagged
             vector. Everything the insert half would say is simply left out. */
          const noInsert = !insSeq.length;
          const r = v.rec || v.from;
          const ov = Math.max(15, Math.min(60, parseInt($('#gbOv').value, 10) || defOv));
          const anTm = Math.max(50, Math.min(72, parseInt($('#gbAn').value, 10) || 60));
          const poi = ($('#gbInsName').value || '').trim() || part.name || 'insert';
          const span = $('#gbInsSpan');
          if (span) {
            if (part.bad) span.innerHTML = '<span style="color:var(--bad)">that range is outside the sequence</span>';
            else if (part.trimmed) {
              const aaFrom = Math.floor(part.from / 3) + 1, aaTo = Math.floor(part.to / 3);
              span.textContent = 'using bases ' + fmt(part.from + 1) + '\u2013' + fmt(part.to) + ' of ' + fmt(part.whole) +
                (part.from % 3 === 0 ? ' (residues ' + fmt(aaFrom) + '\u2013' + fmt(aaTo) + ')' : '') +
                ' \u2014 ' + fmt(insSeq.length) + ' bp';
            } else span.textContent = 'all ' + fmt(part.whole) + ' bp; type a range to clone only part of it';
          }

          const open = openVector(r, v.seq, v.circular);
          if (open.ok) open.intended = place() && place().kind === 'ipcr' ? /MCS|multiple cloning|polylinker/i : null;
          if (!open.ok) { host.innerHTML = '<div class="note bad">' + esc(open.error) + '</div>'; return; }
          const vecSeq = open.seq;

          if (noInsert) { backboneOnly(host, r, v, open, vecSeq); return; }

          const des = B.designAssemblyPrimers(vecSeq, [{ name: poi, seq: insSeq }],
            { overlap: (tuned.ins && tuned.ins.arm) || ov, annealTm: anTm, method: method,
              forceAnneal: tuned.ins, conditions: PB.conditions ? PB.conditions() : undefined });
          const p = des.primers[0];
          const anc = Math.min(18, p.forward.anneal.length, p.reverse.anneal.length);
          const amps = B.pcr(insSeq, p.forward.seq, p.reverse.seq, { anchor: anc });
          /* The insert is almost never a bare fragment in a tube: it is a
             gene sitting in another plasmid, and that whole plasmid is the
             template. Checking the primers against the gene alone would miss
             a second site in the donor's backbone, which is where a spurious
             band actually comes from. */
          const iv = seqValue('gbIns');
          const donor = iv.rec || iv.from || null;
          const donorAmps = donor && donor.seq && donor.seq.length > insSeq.length
            ? B.pcr(donor.seq, p.forward.seq, p.reverse.seq,
                { anchor: anc, circular: donor.topology === 'circular' })
            : null;
          const want = insSeq.length + p.forward.tail.length + p.reverse.tail.length;
          const main = amps.find(x => x.length === want) || amps.find(x => x.seq.indexOf(insSeq) !== -1) || amps[0];
          const asm = main ? B.gibsonAssemble([{ name: 'vector', seq: vecSeq }, { name: 'insert', seq: main.seq }],
            { minOverlap: Math.min(ov, 15), circular: true }) : { ok: false, error: 'The primers give no product on the insert.' };

          let html = '';
          /* 1 — the plasmid you end up with */
          if (asm.ok) {
            html += '<div class="note ok" style="display:flex;gap:14px;align-items:center;flex-wrap:wrap">' +
              '<strong>' + esc((r ? r.name : 'vector')) + ' + ' + esc(poi) + ' — ' + fmt(asm.product.length) + ' bp circular.</strong>' +
              '<span>' + esc(open.how) + '</span><span>' + ov + ' bp arms</span></div>';
          } else {
            html += '<div class="note bad"><strong>This does not assemble.</strong> ' + esc(asm.error) + '</div>';
          }
          /* 2 — frame and what the insert carries */
          const p0 = place();
          if (p0 && p0.essential) html += '<div class="note bad"><strong>This deletes the ' + esc(p0.essential) +
            '.</strong> ' + esc(p0.label.replace(/ \u2014.*/, '')) + ' removes it from the plasmid, so after transformation ' +
            (/resistance/.test(p0.essential) ? 'nothing will grow on that antibiotic'
              : 'the plasmid will not be copied in the cell') +
            '. Unless you mean to swap it, choose a cloning site instead.</div>';
          const frame = fusionReport(r, open, insSeq, ov);
          html += frame;
          html += insertNotes(insSeq, /in frame<\/strong>/.test(frame), part);
          html += addonNote(insSeq, vecSeq, part);
          /* with a truncation, the ends of the protein are what you check */
          if (part && part.trimmed && insSeq.length >= 24) {
            const prot = B.translate(insSeq.slice(0, insSeq.length - insSeq.length % 3), { table: state.codeTable });
            if (prot.length > 8) html += '<div class="row"><div class="field wide">' +
              '<label>What this piece codes for, read from its first base</label>' +
              '<div class="seqblock" style="max-height:none">' +
              '<span style="color:var(--accent);font-weight:600">' + esc(prot.slice(0, 10)) + '</span>' +
              ' \u2026 ' + fmt(prot.length) + ' residues \u2026 ' +
              '<span style="color:var(--accent);font-weight:600">' + esc(prot.slice(-10)) + '</span></div></div></div>';
          }
          /* 3 — the primers to order */
          /* [role, name, sequence, Tm of the annealing part, Tm of the whole
             primer, tail length] — a tailed primer has two useful melting
             temperatures: the first cycles anneal on the 3' part alone, later
             cycles on the whole thing. */
          const rows = [];
          rows.push(['Insert forward', primerName('Insert forward', poi + '_F'), p.forward.seq, p.forward.tmAnneal, p.forward.tmFull, p.forward.tail.length]);
          rows.push(['Insert reverse', primerName('Insert reverse', poi + '_R'), p.reverse.seq, p.reverse.tmAnneal, p.reverse.tmFull, p.reverse.tail.length]);
          if (open.primers) {
            rows.push(['Backbone forward', primerName('Backbone forward', poi + '_bb_F'), open.primers.forward.seq, open.primers.forward.tm, open.primers.forward.tm, 0]);
            rows.push(['Backbone reverse', primerName('Backbone reverse', poi + '_bb_R'), open.primers.reverse.seq, open.primers.reverse.tm, open.primers.reverse.tm, 0]);
          }
          const tcell = t => (isFinite(t) ? f1(t) : '\u2014');
          /* A primer the designer could not bring to the target is still the
             best one available, and it is handed over either way — but it
             has to say so. The annealing length has a floor, and on a
             GC-rich end even the shortest allowed stretch melts well above
             60, so a silent 66 reads as a primer that met the target when
             nothing in the range could. The number is only usable if it is
             the number the block will be set to. */
          const offBy = t => (isFinite(t) ? t - anTm : 0);
          const tmBadge = t => {
            const d = offBy(t);
            if (!isFinite(t) || Math.abs(d) <= 3) return '';
            return ' <span class="badge ' + (Math.abs(d) > 5 ? 'bad' : 'warn') + '" title="' +
              'You asked for ' + f1(anTm) + ' \u00b0C. Nothing in the allowed range of lengths came closer ' +
              'than this; open \u201cOther lengths\u201d to see what was available.">' +
              (d > 0 ? '+' : '\u2212') + f1(Math.abs(d)) + ' off ' + f1(anTm) + '</span>';
          };
          html += '<div class="scrollx" style="margin-top:8px"><table class="tbl primers"><thead><tr><th>Primer</th><th>Name</th>' +
            '<th>Sequence (5\u2032\u21923\u2032)</th><th class="num">nt</th><th class="num">Tm anneal</th>' +
            '<th class="num">Tm full</th><th class="num">Tail</th></tr></thead><tbody>' +
            rows.map(x => '<tr><td>' + esc(x[0]) + '</td>' + nameCell(x[0], x[1]) +
              '<td class="seqcell mono" style="font-size:11.5px;max-width:320px">' +
                (x[5] ? '<span style="color:var(--accent)">' + esc(x[2].slice(0, x[5])) + '</span>' + esc(x[2].slice(x[5])) : esc(x[2])) +
              '</td><td class="num">' + x[2].length + '</td><td class="num">' + tcell(x[3]) + tmBadge(x[3]) + '</td>' +
              '<td class="num">' + (x[5] ? tcell(x[4]) : '\u2014') + '</td>' +
              '<td class="num">' + (x[5] || '\u2014') + '</td></tr>').join('') +
            '</tbody></table></div>' +
            /* Two primers that go in one tube have to anneal at one
               temperature. A pair several degrees apart means whichever
               number the block is set to is wrong for one of them: too high
               and the cooler primer does not bind, too low and the warmer
               one binds where it should not. Worth saying out loud, because
               each primer on its own can look perfectly reasonable. */
            /* A backbone primer that matches the plasmid twice gives two
               products from one inverse PCR, and the gel then shows a band
               that looks right next to one that is not. The engine already
               lengthens a primer to escape a repeat where it can; where it
               cannot, that has to be said rather than left in the data. */
            (function () {
              if (!open.primers) return '';
              const bad = [['forward', open.primers.forward], ['reverse', open.primers.reverse]]
                .filter(([, pr]) => pr && pr.unique === false);
              if (!bad.length) return '';
              return '<div class="note bad"><strong>The backbone ' +
                bad.map(([w]) => w).join(' and ') + ' primer' + (bad.length > 1 ? 's match' : ' matches') +
                ' the plasmid in more than one place.</strong> No length of ' +
                (bad.length > 1 ? 'either' : 'it') + ' binds in one place only, so the opening sits in a ' +
                'repeat \u2014 two linkers, two tags, a duplicated site. The reaction will give more than one ' +
                'product. Move the opening out of the repeat, or cut the backbone with an enzyme instead of ' +
                'amplifying it.</div>';
            })() +
            (function () {
              const pairs = [];
              if (rows.length >= 2) pairs.push(['The insert pair', rows[0], rows[1]]);
              if (rows.length >= 4) pairs.push(['The backbone pair', rows[2], rows[3]]);
              return pairs.map(([what, a, b]) => {
                if (!isFinite(a[3]) || !isFinite(b[3])) return '';
                const d = Math.abs(a[3] - b[3]);
                if (d <= 4) return '';
                const lo = a[3] < b[3] ? a : b, hi = a[3] < b[3] ? b : a;
                return '<div class="note ' + (d > 6 ? 'bad' : 'warn') + '"><strong>' + what +
                  ' are ' + f1(d) + ' \u00b0C apart</strong> \u2014 ' + esc(lo[1]) + ' anneals at ' +
                  f1(lo[3]) + ' and ' + esc(hi[1]) + ' at ' + f1(hi[3]) + '. They go in the same tube, ' +
                  'so one annealing temperature has to serve both: set it for the cooler one and the ' +
                  'warmer primer will also bind where it half-matches; set it for the warmer and the ' +
                  'cooler one will not bind at all. Around ' + f1(Math.min(a[3], b[3]) - 2) +
                  ' \u00b0C is the usual compromise, with the length of the warmer primer trimmed in ' +
                  '\u201cOther lengths\u201d if you would rather close the gap.</div>';
              }).join('');
            })() +
            '<div class="hint">Coloured bases are the homology tail. It does not anneal to the template in the first cycles, so set the annealing temperature from <em>Tm anneal</em>; later cycles run on the whole primer. ' +
            (open.ipcr ? 'Amplify the backbone with its pair and DpnI the template before assembly. ' : '') +
            (open.preLinear
              ? 'The vector is already linear, so it needs no primers of its own \u2014 these two are for the insert, and their tails are its two ends. Mix the insert PCR product with the linear vector and the assembly closes the circle. '
              : 'Order all four and the reaction is one tube. ') +
            'Names are yours to change \u2014 type over one and it is used everywhere below; empty it to get the automatic name back.</div>';
          /* 4 — the junctions */
          if (asm.ok) {
            const fName = rows[0] ? rows[0][1] : 'insert F', rName = rows[1] ? rows[1][1] : 'insert R';
            /* the two ends, drawn as the strands that actually pair */
            html += annealBlock('Where the insert meets the backbone \u2014 5\u2032 end',
              'backbone', vecSeq.slice(-14 - ov, -ov), vecSeq.slice(-ov),
              'insert  ', insSeq.slice(0, 14),
              'The backbone\u2019s own strand pairs with the complement of the tail on <strong>' + esc(fName) +
              '</strong>. That tail is the backbone\u2019s last ' + ov + ' bases copied exactly, which is why they pair along their whole length.');
            html += annealBlock('Where the insert meets the backbone \u2014 3\u2032 end',
              'insert  ', insSeq.slice(-14 - ov, -ov) || insSeq.slice(0, Math.max(0, insSeq.length - ov)).slice(-14),
              insSeq.slice(-ov),
              'backbone', vecSeq.slice(ov, ov + 14),
              'Here the insert\u2019s strand pairs with the complement of the tail on <strong>' + esc(rName) +
              '</strong>, which is the backbone\u2019s first ' + ov + ' bases, reverse complemented.');
            html += primerRoles(rows, open, ov, vecSeq, insSeq, v);
            /* a fragment that already carries the arm */
            [['start', '5\u2032'], ['end', '3\u2032']].forEach(([which, lbl]) => {
              const hit = armAlreadyThere(insSeq, vecSeq, ov, which);
              if (!hit) return;
              const exact = hit.diff === 0;
              html += '<div class="note ' + (exact ? 'warn' : 'warn') + '"><strong>The insert already carries ' +
                (exact ? 'this arm' : 'something very like this arm') + ' at its ' + lbl + ' end.</strong> Its ' +
                (which === 'start' ? 'first ' : 'last ') + hit.k + ' bases are ' +
                (exact ? 'identical to' : f1(hit.id * 100) + '% identical to (' + hit.diff + ' difference' +
                  (hit.diff === 1 ? '' : 's') + ' from)') + ' the backbone\u2019s ' +
                (which === 'start' ? 'last ' : 'first ') + hit.k + ' bases. ' +
                (exact
                  ? 'The primer adds the same arm again, so the finished plasmid would carry it twice. Amplify from base ' +
                    (hit.k + 1) + ' instead \u2014 put that in the <em>from</em> box \u2014 or set Arms to 0 and use the fragment as it is.'
                  : 'Too close to be a coincidence and too far to anneal: homology arms have to match base for base. ' +
                    'The primer supplies a correct arm, so the assembly still works, but the finished plasmid carries both copies. ' +
                    'If those bases were meant as the overlap, amplify from base ' + (hit.k + 1) + ' instead.') +
                '</div>';
            });
          }
          /* 5 — what could go wrong */
          if (amps.length > 1) html += '<div class="note warn"><strong>The insert primers prime at ' + amps.length +
            ' places on the insert</strong> (' + amps.map(x => fmt(x.length)).join(', ') + ' bp). Gel-purify the ' +
            fmt(main.length) + ' bp band, or move the primers.</div>';
          html += donorNote(donor, donorAmps, main, r, insSeq);
          if (main && main.length !== want) html += '<div class="note bad">The full-length product (' + fmt(want) +
            ' bp) was not predicted; what is shown uses a ' + fmt(main.length) + ' bp product and is probably not what you want.</div>';
          if (open.primers) html += primerPairNote(v.seq, open.primers);
          /* an in-place replacement has to keep the frame of whatever is
             behind it: this is the check that a vector swap does not need
             and a tag swap lives or dies by */
          if (open.removed) {
            let fc = null;
            try {
              fc = B.swapFrameCheck(v.seq, open.removed.a, open.removed.b, insSeq,
                { circular: v.circular, features: (r && r.features) || [] });
            } catch (e) { fc = null; }
            html += frameBox(fc);
          }
          html += tunePanel(v, open, p, ov, anTm, insSeq, vecSeq);
          html += openingReport(r, open.points, open.removed, open.intended);
          /* 6 — take it away */
          if (asm.ok) {
            html += '<div class="row" style="margin-top:10px">' +
              '<button class="btn primary" id="gbSave">Save the construct</button>' +
              '<button class="btn" id="gbSavePrimers">Save the primers</button>' +
              '<button class="btn" id="gbCopyPrimers">Copy primers</button></div>';
          }
          host.innerHTML = html;

          lastDesign = { v: v, r: r, open: open, vecSeq: vecSeq, insSeq: insSeq, p: p, asm: asm, ov: ov, poi: poi, rows: rows, part: part };
          state.gbDesign = lastDesign;     // what the dialog currently shows, for inspection
          wireNames(host, rows);
          host.querySelectorAll('[data-trimins]').forEach(b => b.addEventListener('click', () => {
            const [f, t] = b.getAttribute('data-trimins').split(',').map(Number);
            if (state.gbUnit === 'aa') { state.gbUnit = 'bp'; $('#gbInsUnit').textContent = 'bp'; }
            $('#gbInsFrom').value = f; $('#gbInsTo').value = t;
            run();
          }));
          host.querySelectorAll('[data-tunebb]').forEach(b => b.addEventListener('click', () => {
            const [f, r2] = b.getAttribute('data-tunebb').split(',').map(Number);
            tuned.bb = { f: f, r: r2 }; run();
          }));
          host.querySelectorAll('[data-tuneins]').forEach(b => b.addEventListener('click', () => {
            const [which, arm, an] = b.getAttribute('data-tuneins').split(',');
            tuned.ins = Object.assign({}, tuned.ins, { arm: +arm });
            tuned.ins[which] = +an;
            run();
          }));
          const sv = $('#gbSave');
          if (sv) sv.addEventListener('click', saveIt);
          const sp = $('#gbSavePrimers');
          if (sp) sp.addEventListener('click', savePrimers);
          const cp = $('#gbCopyPrimers');
          if (cp) cp.addEventListener('click', () => PB.copyText(rows.map(x => x[1] + '\t' + x[2]).join('\n'), rows.length + ' primers'));
          /* the annealing part is sized to the target, so a primer that still
             falls short has simply run out of template to grow into */
          const warnTm = rows.filter(x => isFinite(x[3]) && x[3] < anTm - 8);
          if (warnTm.length) host.insertAdjacentHTML('beforeend', '<div class="note warn"><strong>' +
            warnTm.map(x => esc(x[0]) + ' reaches only ' + f1(x[3]) + ' \u00b0C').join('; ') +
            '.</strong> That end of the template is AT-rich and the primer is already at its longest \u2014 ' +
            'run the annealing step cooler, or move the end a few bases.</div>');
        }

        /* The vector on its own: the two primers that amplify it without
           whatever was removed, and the plasmid you get by closing it. */
        function backboneOnly(host, r, v, open, vecSeq) {
          const gone = open.removed ? ((open.removed.b - open.removed.a + v.seq.length) % v.seq.length) || 0 : 0;
          let html = '<div class="note ok" style="display:flex;gap:14px;align-items:center;flex-wrap:wrap">' +
            '<strong>' + esc(r ? r.name : 'vector') + ' backbone \u2014 ' + fmt(vecSeq.length) + ' bp linear.</strong>' +
            '<span>' + esc(open.how) + '</span></div>';
          if (!open.primers) {
            html += '<div class="note">That placement cuts with an enzyme rather than amplifying, so there are no primers to design. ' +
              'Choose an inverse-PCR placement \u2014 <em>replace\u2026</em>, a range, or a single base \u2014 to get a primer pair.</div>';
            host.innerHTML = html; lastDesign = null; state.gbDesign = null; return;
          }
          const rows = [
            ['Backbone forward', primerName('Backbone forward', 'bb_F'), open.primers.forward.seq, open.primers.forward.tm, open.primers.forward.tm, 0],
            ['Backbone reverse', primerName('Backbone reverse', 'bb_R'), open.primers.reverse.seq, open.primers.reverse.tm, open.primers.reverse.tm, 0]
          ];
          const tcell = t => (isFinite(t) ? f1(t) : '\u2014');
          html += '<div class="scrollx" style="margin-top:8px"><table class="tbl primers"><thead><tr><th>Primer</th><th>Name</th>' +
            '<th>Sequence (5\u2032\u21923\u2032)</th><th class="num">nt</th><th class="num">Tm</th></tr></thead><tbody>' +
            rows.map(x => '<tr><td>' + esc(x[0]) + '</td>' + nameCell(x[0], x[1]) +
              '<td class="seqcell mono" style="font-size:11.5px;max-width:340px">' + esc(x[2]) + '</td>' +
              '<td class="num">' + x[2].length + '</td><td class="num">' + tcell(x[3]) + '</td></tr>').join('') +
            '</tbody></table></div>';
          html += '<div class="hint">' + (gone ? 'These two amplify everything except the ' + fmt(gone) + ' bp removed. ' : 'These two amplify the whole plasmid. ') +
            'DpnI the template afterwards, or the original plasmid will transform straight through. ' +
            'This product is what goes into the assembly with your insert \u2014 add one above and the insert primers appear with it.</div>';
          /* what the two ends run into when the product is closed on itself */
          html += junctionBox('Where the ends meet if you close it on itself',
            '<span style="color:var(--fg-faint)">' + esc(vecSeq.slice(-18)) + '</span>' +
            '<span style="color:var(--accent);font-weight:700"> | </span>' +
            '<span style="color:var(--fg-faint)">' + esc(vecSeq.slice(0, 18)) + '</span>');
          if (open.removed) {
            let fc = null;
            try {
              fc = B.swapFrameCheck(v.seq, open.removed.a, open.removed.b, '',
                { circular: v.circular, features: (r && r.features) || [] });
            } catch (e) { fc = null; }
            html += frameBox(fc);
          }
          html += '<div class="note">Inverse PCR leaves blunt ends, so closing this on itself needs a blunt ligation (KLD) rather than ' +
            'Gibson, which has no homology to work with. Most of the time you do not close it at all \u2014 you assemble it with the insert.</div>';
          html += primerPairNote(v.seq, open.primers);
          html += openingReport(r, open.points, open.removed, open.intended);
          const p0 = place();
          if (p0 && p0.essential) html += '<div class="note bad"><strong>This deletes the ' + esc(p0.essential) +
            '.</strong> Unless you mean to swap it, choose a cloning site instead.</div>';
          /* The PCR product itself, kept linear, is the thing that goes in the
             tube. Saving it means the insert can be designed against it later
             — in another sitting, or against a gene that is not loaded yet —
             and its two ends are then what the insert's arms are built from. */
          html += '<div class="row" style="margin-top:10px">' +
            '<button class="btn primary" id="gbSaveLin">Save the linear backbone</button>' +
            '<button class="btn" id="gbSaveBB">Save it closed on itself</button>' +
            '<button class="btn" id="gbSavePrimers">Save the primers</button>' +
            '<button class="btn" id="gbCopyPrimers">Copy primers</button></div>' +
            '<div class="hint">Save it linear to use it as the vector for an insert later: pick it in the ' +
            '<em>Vector</em> box and the insert primers are built to match these two ends. ' +
            'Closing it on itself is only for when you want the empty plasmid as a reagent.</div>';
          host.innerHTML = html;
          lastDesign = { v: v, r: r, open: open, vecSeq: vecSeq, insSeq: '', rows: rows,
                         poi: (r ? r.name : 'vector') + ' backbone', backboneOnly: true };
          state.gbDesign = lastDesign;
          wireNames(host, rows);
          const sl = $('#gbSaveLin');
          if (sl) sl.addEventListener('click', () => saveProduct({
            name: (r ? r.name : 'vector') + ' backbone, linear (' + fmt(vecSeq.length) + ' bp)',
            seq: vecSeq, circular: false,
            features: transferFeatures(vecSeq, [r].filter(Boolean), false),
            description: 'Inverse-PCR product' + (gone ? ', ' + fmt(gone) + ' bp removed' : '') +
              '. Use as the vector for an assembly: its ends are the homology arms.',
            history: [[isInf ? 'In-Fusion cloning' : 'Gibson / HiFi assembly',
              'Linear backbone: ' + open.how + '; primers ' + open.primers.forward.seq + ' / ' + open.primers.reverse.seq]]
          }));
          const sb = $('#gbSaveBB');
          if (sb) sb.addEventListener('click', () => saveProduct({
            name: (r ? r.name : 'vector') + ' minus ' + (gone ? fmt(gone) + ' bp' : 'nothing') + ' (' + fmt(vecSeq.length) + ' bp)',
            seq: vecSeq, circular: true,
            features: transferFeatures(vecSeq, [r].filter(Boolean), true),
            history: [[isInf ? 'In-Fusion cloning' : 'Gibson / HiFi assembly',
              'Backbone only: ' + open.how + '; primers ' + open.primers.forward.seq + ' / ' + open.primers.reverse.seq]]
          }));
          const sp2 = $('#gbSavePrimers');
          if (sp2) sp2.addEventListener('click', savePrimers);
          const cp2 = $('#gbCopyPrimers');
          if (cp2) cp2.addEventListener('click', () => PB.copyText(rows.map(x => x[1] + '\t' + x[2]).join('\n'), rows.length + ' primers'));
        }

        function insertFeatureOf(d, productSeq) {
          const P = productSeq.length, hay = productSeq + productSeq;
          const i = hay.indexOf(d.insSeq);
          if (i === -1 || i >= P) return null;
          const end = i + d.insSeq.length;
          const prot = B.translate(d.insSeq, { table: state.codeTable });
          const clean = /^M/.test(prot) && /\*$/.test(prot) && prot.slice(0, -1).indexOf('*') === -1 &&
                        d.insSeq.length % 3 === 0;
          return PB.normFeature({
            type: clean ? 'CDS' : 'misc_feature',
            start: i, end: end > P ? end - P : end, strand: 1, label: d.poi, color: '#b1ff67',
            quals: clean ? { codon_start: '1', transl_table: String(state.codeTable),
                             translation: prot.replace(/\*$/, ''), product: d.poi } : { note: d.poi }
          });
        }
        function saveIt() {
          const d = lastDesign;
          if (!d || !d.asm.ok) return;
          saveProduct({
            name: (d.r ? d.r.name : 'vector') + ' + ' + d.poi + ' (' + fmt(d.asm.product.length) + ' bp)',
            seq: d.asm.product.seq, circular: true,
            features: (f => { const g = insertFeatureOf(d, d.asm.product.seq); return g ? [g].concat(f) : f; })(
              transferFeatures(d.asm.product.seq, [d.r, d.v.rec || d.v.from].filter(Boolean), true)),
            history: [[isInf ? 'In-Fusion cloning' : 'Gibson / HiFi assembly',
              'Vector ' + d.open.how + '. ' + d.ov + ' bp homology arms; insert primers ' +
              d.p.forward.seq + ' / ' + d.p.reverse.seq +
              (d.open.primers ? '; backbone primers ' + d.open.primers.forward.seq + ' / ' + d.open.primers.reverse.seq : '')]]
          });
        }
        function savePrimers() {
          const d = lastDesign;
          if (!d) return;
          const tag = d.poi.replace(/[^\w.-]+/g, '_').slice(0, 14);
          /* a name typed by hand is kept exactly; only an automatic one gets
             the construct tag put in front of it */
          /* The annealing portion goes with the primer. Which part of a
             tailed primer binds is a fact about the primer and its
             template, and the designer is the only thing that knows it for
             certain — recovered later from whatever construct happens to be
             open, it would be a guess, and against the assembled product
             there is no tail at all because the tail is in the product by
             then. Saved here, the library can show the temperature the
             first cycles actually need instead of the whole primer's. */
          d.rows.forEach(x => state.primers.push({
            name: ((d.backboneOnly && !primerNames[x[0]]) ? tag + '_' + x[1] : x[1]).replace(/[^\w.-]+/g, '_'),
            seq: x[2],
            anneal: x[5] ? x[2].slice(x[5]) : undefined,
            tail: x[5] ? x[2].slice(0, x[5]) : undefined }));
          PB.lsSet('primers', state.primers);
          PB.toast('Saved ' + d.rows.length + ' primers to the library', 'ok');
        }

        /* every input feeds the same recompute */
        let t = null;
        const soon = () => { clearTimeout(t); t = setTimeout(run, 110); };
        ['#gbIns_src', '#gbIns_txt', '#gbVec_src', '#gbVec_txt', '#gbInsPart', '#gbWhere',
         '#gbFrom', '#gbTo', '#gbOv', '#gbAn', '#gbInsName', '#gbInsFrom', '#gbInsTo'].forEach(id => {
          const el = w.querySelector(id);
          if (!el) return;
          el.addEventListener('input', soon);
          el.addEventListener('change', soon);
        });
        w.querySelector('#gbIns_src').addEventListener('change', () => {
          insChosen = true;
          $('#gbInsPart').value = ''; $('#gbInsFrom').value = ''; $('#gbInsTo').value = '';
          refreshParts(false); soon();
        });
        /* Naming the vector resolves the collision that forced "vector only":
           whatever else is loaded is now a plausible insert, so point at it. */
        /* Pasting is itself the choice. The box below the dropdown exists to
           be pasted into, so text in it wins over whatever the dropdown
           still says — otherwise a pasted fragment is silently ignored in
           favour of a record nobody picked. Emptying it hands control back. */
        w.querySelector('#gbIns_txt').addEventListener('input', () => {
          const i = $('#gbIns_src');
          if (!i) return;
          const has = !!$('#gbIns_txt').value.trim();
          if (has && i.value !== '') {
            i.value = ''; insChosen = true;
            $('#gbInsPart').value = ''; $('#gbInsFrom').value = ''; $('#gbInsTo').value = '';
            if (state.gbUnit === 'aa') { state.gbUnit = 'bp'; $('#gbInsUnit').textContent = 'bp'; }
            refreshParts(false);
          } else if (!has && i.value === '') {
            i.value = '__none'; refreshParts(false);
          }
        });
        w.querySelector('#gbVec_txt').addEventListener('input', () => {
          const v = $('#gbVec_src');
          if (v && $('#gbVec_txt').value.trim() && v.value !== '') { v.value = ''; refreshPlaces(); }
        });
        w.querySelector('#gbVec_src').addEventListener('change', () => {
          const i = $('#gbIns_src'), v = $('#gbVec_src');
          if (!i || !v || insChosen || i.value !== '__none') return;
          const other = Array.from(i.options)
            .find(o => o.value === String(state.activeId) && o.value !== v.value);
          if (other) {
            i.value = other.value;
            $('#gbInsPart').value = ''; $('#gbInsFrom').value = ''; $('#gbInsTo').value = '';
            refreshParts(false);
          }
        });
        w.querySelector('#gbInsPart').addEventListener('change', () => {
          $('#gbInsFrom').value = ''; $('#gbInsTo').value = ''; soon();
        });
        w.querySelector('#gbInsUnit').addEventListener('click', () => {
          /* keep the same stretch of DNA when the unit changes, so the
             numbers move but the construct does not */
          const was = insertPart();
          state.gbUnit = (state.gbUnit || 'bp') === 'bp' ? 'aa' : 'bp';
          $('#gbInsUnit').textContent = state.gbUnit;
          if (was.trimmed) {
            if (state.gbUnit === 'aa') {
              $('#gbInsFrom').value = Math.floor(was.from / 3) + 1;
              $('#gbInsTo').value = Math.floor(was.to / 3);
            } else {
              $('#gbInsFrom').value = was.from + 1;
              $('#gbInsTo').value = was.to;
            }
          }
          soon();
        });
        w.querySelector('#gbVec_src').addEventListener('change', () => { placeCache = ''; $('#gbFrom').value = ''; soon(); });
        w.querySelector('#gbWhere').addEventListener('change', () => {
          whereChosen = true; tuned = {}; $('#gbFrom').value = ''; applyPlace(); soon();
        });
        refreshParts(false); refreshPlaces(); run();

        w.querySelector('#gbAsmGo').addEventListener('click', () => {
          const frags = PB.parseFastaish($('#gbFrags').value);
          if (frags.length < 1) { PB.toast('Paste at least one fragment.', 'warn'); return; }
          const minOv = Math.max(10, parseInt($('#gbMinOv').value, 10) || 15);
          const asm = B.gibsonAssemble(frags, { minOverlap: minOv, circular: circ });
          let html;
          if (!asm.ok) {
            html = '<div class="note bad"><strong>Assembly failed.</strong> ' + esc(asm.error) + '</div>' +
              (asm.order ? '<div class="hint">Chained: ' + esc(asm.order.join(' → ')) + '</div>' : '');
          } else {
            html = '<div class="note ok"><strong>' + fmt(asm.product.length) + ' bp ' +
              (asm.product.circular ? 'circular' : 'linear') + ' product.</strong><br>Order: ' +
              esc(asm.order.join(' → ')) + '<br>Overlaps: ' + asm.overlaps.join(', ') + ' bp' +
              (asm.closingOverlap ? ' (closing ' + asm.closingOverlap + ' bp)' : '') + '</div>' +
              '<button class="btn primary" id="gbAsmSave">Save construct</button>';
          }
          result(html);
          const s = $('#gbAsmSave');
          if (s) s.addEventListener('click', () => saveProduct({
            name: 'assembly (' + fmt(asm.product.length) + ' bp)',
            seq: asm.product.seq, circular: asm.product.circular,
            features: transferFeatures(asm.product.seq, state.records, asm.product.circular),
            history: [['Gibson / HiFi assembly', frags.length + ' fragments: ' + asm.order.join(' → ')]]
          }));
        });
      }
    });
  }

  /* Cloning sites read off the vector's own annotation. For each multiple
     cloning site the vector declares, work out the placements a cloner
     actually wants, from the coding features inside it:
       - fuse after the N-terminal tag, in frame, keeping the tag
       - replace from the start codon, keeping the RBS and its spacing
       - replace the whole site
     Positions come from the annotation, never from a guess, so a vector with
     no MCS feature simply offers nothing here. */
  /* The parts a plasmid cannot lose and still work: its selection marker and
     its origin of replication. Recognised from the label, so a vector that
     does not annotate them simply gets no warning. */
  /* What the insert PCR actually runs on. Amplifying a gene out of one
     plasmid to put it into another is the usual shape of this job, and it
     brings two hazards the gene's own sequence cannot show: a primer that
     also lands somewhere in the donor's backbone, and donor plasmid riding
     through the assembly into the transformation. The second only matters
     when the donor would survive the same plate as the construct. */
  /* Which antibiotics the plasmid actually grows on. The gene has to be
     there for that, so a promoter, terminator or primer-binding site that
     merely carries the name does not count: pRSFDuet-1 keeps an "AmpR
     promoter" annotation from its ancestry and confers kanamycin only. */
  function markersOf(rec) {
    const out = {};
    if (!rec || !rec.features) return out;
    rec.features.forEach(f => {
      const t = (f.type || '').toLowerCase();
      if (t && !/^(cds|gene|orf|misc_feature|misc_marker)$/.test(t)) return;
      if (/promoter|terminator|primer|operator|\brbs\b|ribosome|\bori\b|origin/i.test(f.label || '')) return;
      const m = PB.essentialPart(f);
      if (m && /resistance|selection marker/.test(m)) out[m] = (f.label || f.type || m);
    });
    return out;
  }
  function donorNote(donor, donorAmps, main, destRec, insSeq) {
    if (!donor || !donorAmps) return '';
    let html = '';
    const wanted = main ? main.length : 0;
    const extras = donorAmps.filter(x => x.length !== wanted);
    const hasWanted = donorAmps.some(x => x.length === wanted);
    const name = donor.name || 'the donor plasmid';
    if (!hasWanted) {
      html += '<div class="note bad"><strong>On ' + esc(name) + ' itself these primers do not give the ' +
        fmt(wanted) + ' bp product</strong>' +
        (donorAmps.length ? ' — they give ' + donorAmps.map(x => fmt(x.length)).join(', ') + ' bp instead' : ' — no product at all') +
        '. The insert sequence and the plasmid you would amplify it from do not agree; check you pasted the right one.</div>';
    } else if (extras.length) {
      html += '<div class="note warn"><strong>Amplified from ' + esc(name) + ', these primers also give ' +
        (extras.length === 1 ? 'a second product' : extras.length + ' further products') + '</strong> (' +
        extras.map(x => fmt(x.length)).join(', ') + ' bp alongside the ' + fmt(wanted) +
        ' bp you want). Gel-purify the right band before the assembly.</div>';
    } else {
      html += '<div class="note ok"><strong>Checked against ' + esc(name) + ':</strong> one product, ' +
        fmt(wanted) + ' bp. That is the band to gel-purify, or column-clean after DpnI.</div>';
    }
    /* carryover: uncut donor transforming alongside the assembly */
    const dm = markersOf(donor), vm = markersOf(destRec);
    const shared = Object.keys(dm).filter(k => vm[k]);
    if (shared.length) {
      html += '<div class="note warn"><strong>' + esc(name) + ' carries the same ' + esc(shared.join(' and ')) +
        ' as the construct.</strong> Any donor plasmid surviving the PCR grows on the same plate, so colonies can be donor rather than assembly. ' +
        'DpnI the insert PCR as well as the backbone, and gel-purify rather than column-cleaning it.</div>';
    } else if (Object.keys(dm).length && Object.keys(vm).length) {
      html += '<div class="note"><strong>Selection:</strong> ' + esc(name) + ' is ' + esc(Object.keys(dm).join(' and ')) +
        ' and the construct is ' + esc(Object.keys(vm).join(' and ')) +
        ', so leftover donor will not grow on your plate.</div>';
    }
    return html;
  }

  function essentialPart(f) {
    const s = ((f.label || '') + ' ' + (f.type || '')).toLowerCase();
    if (/rep_origin/.test((f.type || '').toLowerCase()) || /\bori\b|origin of replication/.test(s)) return 'origin of replication';
    if (/kanr|neo\b|nptii/.test(s)) return 'kanamycin resistance';
    if (/ampr|\bbla\b|beta-lactamase/.test(s)) return 'ampicillin resistance';
    if (/cmr|cat\b|chloramphenicol/.test(s)) return 'chloramphenicol resistance';
    if (/tetr|tetracycline/.test(s)) return 'tetracycline resistance';
    if (/specr|smr|spectinomycin|streptomycin|aada/.test(s)) return 'spectinomycin resistance';
    if (/geneticin|hygr|hygromycin|puror|puromycin|bsd|blasticidin|zeo/.test(s)) return 'selection marker';
    if (/\brop\b|\brepa\b|\brep\b(?!ress)/.test(s)) return 'replication protein';
    return null;
  }
  PB.essentialPart = essentialPart;

  function cloningSites(rec) {
    if (!rec || !rec.features) return [];
    const out = [];
    const isCDS = f => /^CDS$/i.test(f.type) && f.strand !== -1 && !f.between;
    rec.features.forEach(f => {
      if (!/MCS|multiple cloning|polylinker/i.test((f.label || '') + ' ' + ((f.quals && f.quals.note) || ''))) return;
      if (f.end <= f.start) return;                       // not handling a wrapped MCS here
      const name = f.label || 'cloning site';
      const inside = rec.features.filter(g => isCDS(g) && g.start >= f.start - 3 && g.end <= f.end + 3 && g.end > g.start)
                                 .sort((a, b) => a.start - b.start);
      const atg = inside.find(g => rec.seq.substr(g.start, 3) === 'ATG');
      const org = atg ? atg.start : (inside.length ? inside[0].start : null);
      if (inside.length && org != null) {
        const last = inside[inside.length - 1];
        let p = last.end;
        while ((p - org) % 3 !== 0 && p < f.end) p++;      // step to the next in-frame base
        if ((p - org) % 3 === 0 && p < f.end) {
          out.push({ label: name + ' · after ' + (last.label || 'the tag') + ', in frame (keeps it as an N-terminal fusion)',
                     ex: name + ' · opened after ' + (last.label || 'the tag') + ' \u2014 nothing left out',
                     from: p + 1, to: '', keeps: inside.map(g => g.label || g.type).join(', ') });
        }
      }
      if (atg) {
        out.push({ label: name + ' · replace from the ATG (untagged, keeps the vector’s RBS spacing)',
                   ex: name + ' · everything except the ATG onwards (' + fmt(atg.start + 1) + '\u2013' + fmt(f.end) + ')',
                   from: atg.start + 1, to: f.end, keeps: 'RBS and promoter' });
      }
      out.push({ label: name + ' · replace the whole site (' + fmt(f.start + 1) + '–' + fmt(f.end) + ')',
                 ex: name + ' · everything except the site itself (' + fmt(f.start + 1) + '\u2013' + fmt(f.end) + ')',
                 from: f.start + 1, to: f.end, keeps: 'RBS and promoter' });
    });
    return out;
  }
  PB.cloningSites = cloningSites;

  /* What an opening lands in. A vector is opened to receive an insert, so the
     one thing worth saying out loud is which annotated elements the opening
     interrupts. This is what catches opening a plasmid at base 1. */
  function openingReport(rec, points, removed, intended) {
    if (!rec || !rec.features || !rec.features.length) return '';
    const L = rec.seq.length;
    const hit = [];
    rec.features.forEach(f => {
      if (f.type === 'source' || f.between) return;
      /* opening the cloning site is the point of the exercise, not damage */
      if (intended && intended.test((f.label || '') + ' ' + ((f.quals && f.quals.note) || ''))) return;
      const pieces = B.featurePieces(f, L);
      let why = null;
      pieces.forEach(pc => {
        (points || []).forEach(pt => { if (pt > pc.start && pt < pc.end) why = 'split'; });
        if (removed && removed.b > removed.a) {
          const ov = (x0, x1) => x0 < pc.end && x1 > pc.start;
          const gone = removed.b >= removed.a ? ov(removed.a, removed.b) : (ov(removed.a, L) || ov(0, removed.b));
          if (gone && why !== 'split') why = 'cut into';
        }
      });
      if (why) hit.push(why + ' <b>' + esc(f.label || f.type) + '</b> (' + esc(B.formatLocation(f, L) || '') + ')');
    });
    if (!hit.length) return '';
    return '<div class="note bad"><strong>This opening interrupts annotated elements:</strong> ' +
      hit.join(', ') + '. Check this is where you mean to put the insert — a promoter, RBS or tag broken here will not work in the cell.</div>';
  }

  /* A backbone primer pair has to carry several kb, so a weak or badly
     balanced pair is worth saying out loud, as is one whose 3' end matches
     more than one place on the plasmid. */
  function primerPairNote(seqCirc, pr) {
    const msgs = [];
    const lo = Math.min(pr.forward.tm, pr.reverse.tm), d = Math.abs(pr.forward.tm - pr.reverse.tm);
    if (lo < 55) msgs.push('the weaker primer melts at ' + f1(lo) + ' °C — low for amplifying a whole backbone; raise the target or move the opening a few bases');
    if (d > 5) msgs.push('the two differ by ' + f1(d) + ' °C, so one will mis-prime before the other anneals');
    const hits = sq => {
      let n = 0;
      try { n = B.primerSites(seqCirc, sq, { circular: true }).length; } catch (e) { n = 0; }
      return n;
    };
    const hf = hits(pr.forward.seq), hr = hits(pr.reverse.seq);
    if (hf > 1 || hr > 1) msgs.push('a primer anneals to the plasmid in more than one place (' + hf + ' / ' + hr +
      ') — expect a second product. Move the opening a few bases, or shorten what is removed, so the primer reaches sequence the copies do not share');
    const warn = msgs.length
      ? '<div class="note warn"><strong>Backbone primers:</strong> ' + msgs.join('; ') + '.</div>' : '';
    /* a primer that had to grow to stay unique is worth saying out loud: the
       length looks odd next to the target Tm until you know why */
    const grew = [];
    if (pr.forward.lengthened) grew.push('forward to ' + pr.forward.length + ' nt');
    if (pr.reverse.lengthened) grew.push('reverse to ' + pr.reverse.length + ' nt');
    const note = grew.length
      ? '<div class="note"><strong>One site only:</strong> the sequence beside this opening is repeated elsewhere on the plasmid, so the ' +
        grew.join(' and the ') + ' \u2014 long enough to reach past the repeat. ' +
        'That is why the melting temperature runs above the target.</div>' : '';
    return warn + note;
  }


  /* ---------- what actually anneals ----------
     An assembly works because two single strands find each other. The
     tail is written as identity to the backbone, which is what confuses
     people reading the primer table: the bases that pair are the
     backbone's own strand and the complement of that tail. Drawn out,
     there is nothing left to take on trust.                            */
  function compInOrder(s) { return B.revcomp(s).split('').reverse().join(''); }
  function annealBlock(title, leftName, leftCtx, shared, rightName, rightCtx, caption) {
    /* The three lines have to line up base for base. Padding them with
       spaces does not do it: the prime in 5' is not one character wide in
       most monospace faces, so every line drifts by a fraction and the
       pairing bars end up over the wrong bases. Each prefix is therefore a
       box of an exact number of character widths, which the font cannot
       argue with. */
    const faint = t => '<span style="color:var(--fg-faint)">' + esc(t) + '</span>';
    const hot = t => '<span style="color:var(--accent);font-weight:700">' + esc(t) + '</span>';
    const lab = Math.max(leftName.trim().length, rightName.trim().length) + 4;   // name + " 5' "
    const box = (w, inner) => '<span style="display:inline-block;width:' + w +
      'ch;flex:0 0 auto">' + (inner || '') + '</span>';
    const prefix = (name, end) => box(lab, esc(name.trim()) + ' ' + end + '\u2032');
    const line = inner => '<div style="display:flex;white-space:pre">' + inner + '</div>';
    return '<div class="field wide" style="margin-top:8px"><label>' + esc(title) + '</label>' +
      '<div class="seqblock" style="max-height:none;line-height:1.5;overflow-x:auto">' +
      line(prefix(leftName, '5') + faint(leftCtx) + hot(shared)) +
      line(box(lab) + box(leftCtx.length) +
        '<span style="color:var(--accent)">' + '|'.repeat(shared.length) + '</span>' +
           '<span style="color:var(--fg-faint)">  ' + shared.length + ' bp</span>') +
      line(prefix(rightName, '3') + box(leftCtx.length) +
           hot(compInOrder(shared)) + faint(compInOrder(rightCtx))) +
      '</div>' +
      (caption ? '<div class="hint" style="margin:4px 0 0">' + caption + '</div>' : '') + '</div>';
  }

  /* A fragment ordered with its own homology arms already on it — a gBlock
     usually is — would carry them twice once the primers add them again.
     Near-identity is the dangerous case: close enough to look designed,
     not close enough to anneal. */
  function armAlreadyThere(insSeq, vecSeq, ov, which) {
    const ins = which === 'start' ? insSeq : B.revcomp(insSeq);
    const vec = which === 'start' ? vecSeq : B.revcomp(vecSeq);
    const vEnd = which === 'start' ? vec.slice(-1) : vec.slice(-1);
    let best = null;
    const top = Math.min(60, insSeq.length, vecSeq.length);
    for (let k = top; k >= 10; k--) {
      const a = which === 'start' ? insSeq.slice(0, k) : insSeq.slice(-k);
      const b = which === 'start' ? vecSeq.slice(-k) : vecSeq.slice(0, k);
      if (a.length < k || b.length < k) continue;
      let same = 0;
      for (let i = 0; i < k; i++) if (a[i] === b[i]) same++;
      const id = same / k;
      if (id >= 0.8 && (!best || k > best.k)) best = { k: k, id: id, a: a, b: b, diff: k - same };
    }
    return best;
  }

  function primerBox(label, p) {
    return '<div class="field" style="flex:1 1 300px"><label>' + esc(label) + ' · ' + p.length + ' nt</label>' +
      '<div class="seqblock" style="max-height:none"><span style="color:var(--accent);font-weight:600">' +
      esc(p.tail) + '</span>' + esc(p.anneal) + '</div>' +
      '<dl class="kv" style="margin-top:6px"><dt>Annealing Tm</dt><dd>' + f1(p.tmAnneal) + ' °C</dd>' +
      '<dt>Tail</dt><dd>' + p.tail.length + ' bp homology</dd></dl></div>';
  }

  function simplePrimerBox(label, seq, tm) {
    return '<div class="field" style="flex:1 1 300px"><label>' + esc(label) + ' · ' + seq.length + ' nt</label>' +
      '<div class="seqblock" style="max-height:none">' + esc(seq) + '</div>' +
      '<dl class="kv" style="margin-top:6px"><dt>Tm</dt><dd>' + f1(tm) + ' °C</dd></dl></div>';
  }

  /* ============================== GOLDEN GATE ============================== */
  function toolGoldenGate() {
    PB.modal({
      title: 'Golden Gate assembly', wide: true, sub: 'Type IIS, 4 nt fusion overhangs',
      body: '<div class="row">' + enzField('ggEnz', 'Type IIS enzyme', 'BsaI') +
          '<div class="field" style="flex:0 0 auto;align-self:flex-end"><button class="btn primary" id="ggGo">Assemble</button></div>' +
          '<div class="field" style="flex:0 0 auto;align-self:flex-end"><button class="btn sm" id="ggLib">Add from library</button></div>' +
        '</div>' +
        '<div class="row"><div class="field wide"><label>Parts in FASTA — each flanked by inward-facing enzyme sites</label>' +
        '<textarea class="inp" id="ggParts" rows="9" spellcheck="false" placeholder=">partA&#10;AAGGTCTCAAATG…GCTTAGAGACCAA"></textarea></div></div>' +
        '<div class="hint">Each part is digested, the released fragment is read off with its 4 nt overhangs, and the parts are chained by matching overhangs into a circle. Duplicate overhangs are reported — they are the usual cause of a scrambled assembly.</div>' +
        resultSlot(),
      onOpen: w => {
        w.querySelector('#ggLib').addEventListener('click', () => {
          $('#ggParts').value = state.records.slice(0, 6)
            .map(r => '>' + r.name.replace(/\s+/g, '_') + '\n' + r.seq).join('\n');
        });
        w.querySelector('#ggGo').addEventListener('click', () => {
          const names = enzNames('ggEnz');
          if (!names.length) { PB.toast('Name a Type IIS enzyme (BsaI, BsmBI, BbsI, SapI, PaqCI…).', 'warn'); return; }
          const parts = PB.parseFastaish($('#ggParts').value);
          if (parts.length < 2) { PB.toast('Paste at least two parts.', 'warn'); return; }
          const gg = B.goldenGate(parts.map(p => ({ name: p.name, seq: p.seq })), names[0]);
          let html = '<div class="scrollx" style="margin-top:12px"><table class="tbl"><thead><tr>' +
            '<th>Part</th><th class="num">Released bp</th><th>Left overhang</th><th>Right overhang</th><th class="num">Cuts</th>' +
            '</tr></thead><tbody>' + gg.pieces.map(p => '<tr><td>' + esc(p.name) + '</td>' +
              '<td class="num">' + fmt(p.length) + '</td>' +
              '<td class="mono">' + esc(endStr(p.left)) + '</td>' +
              '<td class="mono">' + esc(endStr(p.right)) + '</td>' +
              '<td class="num">' + p.nCuts + '</td></tr>').join('') + '</tbody></table></div>';
          if (gg.overhangClashes.length) {
            html += '<div class="note bad"><strong>Duplicate overhangs:</strong> ' + esc(gg.overhangClashes.join('; ')) +
              '<br>Give each junction a unique 4 nt overhang.</div>';
          }
          if (!gg.ok) {
            html += '<div class="note bad"><strong>' + esc(gg.error || 'Assembly did not close.') + '</strong>' +
              (gg.order ? '<br>Chained: ' + esc(gg.order.join(' → ')) : '') + '</div>';
          } else {
            html += '<div class="note ok"><strong>Assembles in one order:</strong> ' + esc(gg.order.join(' → ')) +
              '<br>' + fmt(gg.product.length) + ' bp circular product, enzyme sites removed.</div>' +
              '<button class="btn primary" id="ggSave">Save construct</button>';
          }
          result(html);
          const s = $('#ggSave');
          if (s) s.addEventListener('click', () => saveProduct({
            name: 'Golden Gate assembly (' + fmt(gg.product.length) + ' bp)',
            seq: gg.product.seq, circular: true,
            features: transferFeatures(gg.product.seq, state.records, true),
            history: [['Golden Gate assembly', names[0] + '; ' + gg.pieces.length + ' parts in order ' + gg.order.join(' → ')]]
          }));
        });
      }
    });
  }

  /* ============================== GATEWAY ============================== */
  function toolGateway() {
    PB.modal({
      title: 'Gateway recombination', wide: true, sub: 'BP and LR reactions',
      body: '<div class="chipbar"><button class="chip on" data-gw="lr">LR (attL × attR)</button>' +
          '<button class="chip" data-gw="bp">BP (attB × attP)</button></div>' +
        '<div class="row">' +
          seqField('gwA', 'Entry clone (attL1/attL2)', { allowSel: false, rows: 2 }) +
          seqField('gwB', 'Destination vector (attR1/attR2)', { allowSel: false, rows: 2 }) +
        '</div>' +
        '<div class="row"><div class="field" style="flex:0 0 auto"><button class="btn primary" id="gwGo">React</button></div>' +
        '<div class="field" style="flex:0 0 auto"><button class="btn sm" id="gwRef">Show att site sequences</button></div></div>' +
        '<div class="hint" id="gwHint">The cassette between attL1 and attL2 moves into the destination backbone, regenerating attB1 and attB2 in the expression clone.</div>' +
        '<div class="note">No att-site sequences are stored in this page. Both molecules are searched for their att sites, the region the two partners share base for base is found, and the product is cut and joined there — so every base of the product comes from the two sequences you supply. Annotate your att sites (a feature named attL1, attR2 and so on) and they are used directly.</div>' +
        resultSlot(),
      onOpen: w => {
        let mode = 'lr';
        w.querySelectorAll('[data-gw]').forEach(c => c.addEventListener('click', () => {
          w.querySelectorAll('[data-gw]').forEach(x => x.classList.toggle('on', x === c));
          mode = c.getAttribute('data-gw');
          const labels = w.querySelectorAll('label');
          $('#gwHint').textContent = mode === 'lr'
            ? 'The cassette between attL1 and attL2 moves into the destination backbone, regenerating attB1 and attB2 in the expression clone.'
            : 'An attB-flanked PCR product recombines with an attP donor to give the entry clone, flanked by attL1 and attL2.';
          result('');
        }));
        w.querySelector('#gwRef').addEventListener('click', () => {
          const a = seqValue('gwA'), b = seqValue('gwB');
          const show = (what, v) => {
            const rec = v.rec || v.from;
            const sites = v.seq ? B.findAttSites(v.seq, { features: rec ? rec.features : null, circular: v.circular }) : [];
            if (!sites.length) return what + ': no att sites found.\n';
            return what + ':\n' + sites.map(st => '  att' + (st.kind || '?') + st.family + '  ' + (st.start + 1) + '..' + st.end +
              (st.strand === -1 ? ' (reverse)' : '') + '  ' + (st.source === 'annotation' ? 'from your annotation "' + st.label + '"' : 'found by sequence motif') +
              '\n    ' + v.seq.slice(st.start, st.end) + '\n').join('');
          };
          PB.showTextFallback('att sites in your sequences',
            show('First molecule', a) + '\n' + show('Second molecule', b) +
            '\nPlasmidBench does not keep a table of att-site sequences. It finds the region the two partners share exactly ' +
            '(the recombination core) and cuts there, so every base of the product comes from the sequences you supplied. ' +
            'The motifs used to locate candidate sites are ' + B.ATT_PROBES[1] + ' (att1) and ' + B.ATT_PROBES[2] + ' (att2); ' +
            'annotating features named attL1, attR2 and so on makes the search use your annotation instead.');
        });
        w.querySelector('#gwGo').addEventListener('click', () => {
          const a = seqValue('gwA'), b = seqValue('gwB');
          if (!a.seq.length || !b.seq.length) { PB.toast('Provide both inputs.', 'warn'); return; }
          const aRec = a.rec || a.from, bRec = b.rec || b.from;
          const opts = { aFeatures: aRec ? aRec.features : null, bFeatures: bRec ? bRec.features : null,
                         aCircular: mode === 'lr' ? a.circular : false, bCircular: b.circular };
          const r = mode === 'lr' ? B.gatewayLR(a.seq, b.seq, opts) : B.gatewayBP(a.seq, b.seq, opts);
          if (!r.ok) { result('<div class="note bad"><strong>' + esc(r.error) + '</strong></div>'); return; }
          const sources = [aRec, bRec].filter(Boolean);
          result('<div class="note ok"><strong>' + mode.toUpperCase() + ' reaction gives a ' + fmt(r.product.length) +
            ' bp circular product.</strong><br>Recombination cores: ' +
            r.cores.map((c, i) => 'site ' + (i + 1) + ' ' + c.length + ' bp identical in both') .join(', ') + '.' +
            (r.byproduct ? ' By-product ' + fmt(r.byproduct.length) + ' bp.' : '') + '</div>' +
            '<div class="hint">' + esc(r.note) + '</div>' +
            (r.notes.length ? r.notes.map(n => '<div class="note warn">' + esc(n) + '</div>').join('') : '') +
            '<div class="row"><button class="btn primary" id="gwSave">Save ' + (mode === 'lr' ? 'expression' : 'entry') + ' clone</button>' +
            (r.byproduct ? '<button class="btn" id="gwSaveBy">Save by-product</button>' : '') + '</div>');
          $('#gwSave').addEventListener('click', () => saveProduct({
            name: (mode === 'lr' ? 'expression clone' : 'entry clone') + ' (' + fmt(r.product.length) + ' bp)',
            seq: r.product.seq, circular: true,
            features: transferFeatures(r.product.seq, sources, true),
            history: [['Gateway ' + mode.toUpperCase(), r.note + ' Cores ' + r.cores.map(c => c.length + ' bp').join(' and ') +
              (r.flipped ? '; the second molecule was used on its other strand' : '')]]
          }));
          const by = $('#gwSaveBy');
          if (by) by.addEventListener('click', () => saveProduct({
            name: (mode === 'lr' ? 'LR' : 'BP') + ' by-product (' + fmt(r.byproduct.length) + ' bp)',
            seq: r.byproduct.seq, circular: true,
            features: transferFeatures(r.byproduct.seq, sources, true),
            history: [['Gateway ' + mode.toUpperCase() + ' by-product', 'The reciprocal product of the same reaction.']]
          }));
        });
      }
    });
  }

  /* ============================== TOPO / TA ============================== */
  function toolTopo() {
    PB.modal({
      title: 'TOPO / TA / GC cloning', wide: true,
      body: '<div class="row">' +
          '<div class="field"><label>Method</label><select class="inp" id="tpMode">' +
          '<option value="blunt">Blunt TOPO</option><option value="ta">TA cloning (3′-A overhang)</option>' +
          '<option value="directional">Directional TOPO (CACC)</option></select></div>' +
        '</div><div class="row">' +
          seqField('tpIns', 'PCR product (insert)', { allowSel: true, rows: 2 }) +
          seqField('tpVec', 'Vector (optional)', { allowSel: false, rows: 2 }) +
        '</div>' +
        '<div class="row"><button class="btn primary" id="tpGo">Simulate</button></div>' + resultSlot(),
      onOpen: w => {
        w.querySelector('#tpGo').addEventListener('click', () => {
          const ins = seqValue('tpIns'), vec = seqValue('tpVec');
          if (!ins.seq.length) { PB.toast('Provide the PCR product.', 'warn'); return; }
          const r = B.topoClone(ins.seq, { mode: $('#tpMode').value, vectorSeq: vec.seq });
          if (!r.ok) { result('<div class="note bad"><strong>' + esc(r.error) + '</strong></div>'); return; }
          result('<div class="note ok"><strong>' + fmt(r.product.length) + ' bp ' +
            (r.product.circular ? 'circular' : 'linear') + ' product.</strong><br>Orientation: ' + esc(r.orientation) + '</div>' +
            '<div class="hint">' + r.notes.map(esc).join('<br>') + '</div>' +
            '<button class="btn primary" id="tpSave">Save construct</button>');
          $('#tpSave').addEventListener('click', () => saveProduct({
            name: 'TOPO product (' + fmt(r.product.length) + ' bp)',
            seq: r.product.seq, circular: r.product.circular,
            features: transferFeatures(r.product.seq, [vec.rec || vec.from, ins.rec || ins.from].filter(Boolean), r.product.circular),
            history: [['TOPO / TA cloning', r.mode + '; ' + r.orientation]]
          }));
        });
      }
    });
  }

  /* ============================== PCR ============================== */
  function toolPCR() {
    const rec = PB.active();
    PB.modal({
      title: 'PCR simulation', wide: true,
      body: '<div class="row">' + seqField('pcTpl', 'Template', { defaultId: rec ? rec.id : '', allowSel: true, rows: 2 }) + '</div>' +
        '<div class="row">' +
          '<div class="field"><label>Forward primer (5′→3′)</label><textarea class="inp" id="pcF" rows="2" spellcheck="false"></textarea></div>' +
          '<div class="field"><label>Reverse primer (5′→3′)</label><textarea class="inp" id="pcR" rows="2" spellcheck="false"></textarea></div>' +
        '</div><div class="row">' +
          '<div class="field narrow"><label>3′ anchor (bp)</label><input class="inp mono" id="pcAnchor" type="number" value="15" min="8" max="30"></div>' +
          (state.primers.length ? '<div class="field"><label>From library</label><select class="inp" id="pcLib"><option value="">choose…</option>' +
            state.primers.map((p, i) => '<option value="' + i + '">' + esc(p.name) + '</option>').join('') + '</select></div>' : '') +
          '<div class="field" style="flex:0 0 auto;align-self:flex-end"><button class="btn primary" id="pcGo">Amplify</button></div>' +
        '</div>' +
        '<div class="hint">The 3′ anchor is how many bases at each primer’s 3′ end must match the template exactly. 5′ tails that do not match (restriction sites, homology arms, CACC) are carried into the product.</div>' +
        resultSlot(),
      onOpen: w => {
        const lib = w.querySelector('#pcLib');
        if (lib) lib.addEventListener('change', () => {
          const p = state.primers[+lib.value];
          if (!p) return;
          const f = $('#pcF'), r = $('#pcR');
          if (!f.value.trim()) f.value = p.seq; else r.value = p.seq;
          lib.value = '';
        });
        w.querySelector('#pcGo').addEventListener('click', () => {
          const t = seqValue('pcTpl');
          const fwd = B.cleanSeq($('#pcF').value), rev = B.cleanSeq($('#pcR').value);
          if (!t.seq.length) { PB.toast('Provide a template.', 'warn'); return; }
          if (fwd.length < 10 || rev.length < 10) { PB.toast('Both primers need at least 10 bases.', 'warn'); return; }
          const anchor = Math.max(8, parseInt($('#pcAnchor').value, 10) || 15);
          const prods = B.pcr(t.seq, fwd, rev, { anchor, circular: t.circular });
          const cond = PB.conditions();
          let html = '<div class="row" style="margin-top:10px"><div class="field wide"><dl class="kv">' +
            '<dt>Forward Tm</dt><dd>' + f1(B.tmNN(fwd, cond)) + ' °C · ' + fwd.length + ' nt · ' + f1(B.gcContent(fwd)) + '% GC</dd>' +
            '<dt>Reverse Tm</dt><dd>' + f1(B.tmNN(rev, cond)) + ' °C · ' + rev.length + ' nt · ' + f1(B.gcContent(rev)) + '% GC</dd>' +
            '<dt>Cross-dimer</dt><dd>' + B.heteroDimer(fwd, rev).matches + ' bp complementarity</dd>' +
            '</dl></div></div>';
          if (!prods.length) {
            html += '<div class="note bad"><strong>No product.</strong> Neither primer pair orientation anchors on this template. ' +
              'Check the primer directions, or lower the 3′ anchor length.</div>';
            result(html); return;
          }
          html += '<div class="note ' + (prods.length === 1 ? 'ok' : 'warn') + '"><strong>' + prods.length +
            ' product' + (prods.length === 1 ? '' : 's') + ' predicted.</strong>' +
            (prods.length > 1 ? ' More than one product means the primers bind in several places.' : '') + '</div>' +
            '<div class="scrollx"><table class="tbl"><thead><tr><th class="num">Size</th><th class="num">From</th><th class="num">To</th><th>Product</th><th></th></tr></thead><tbody>' +
            prods.slice(0, 12).map((p, i) => '<tr><td class="num"><b>' + fmt(p.length) + '</b></td>' +
              '<td class="num">' + fmt(p.start + 1) + '</td><td class="num">' + fmt(p.end) + '</td>' +
              '<td class="seqcell" style="max-width:260px">' + esc(p.seq.slice(0, 48)) + '…</td>' +
              '<td><button class="btn sm ghost" data-pcsave="' + i + '">Save</button></td></tr>').join('') +
            '</tbody></table></div>';
          result(html);
          $('#toolResult').querySelectorAll('[data-pcsave]').forEach(b => b.addEventListener('click', () => {
            const p = prods[+b.getAttribute('data-pcsave')];
            saveProduct({
              name: (t.rec ? t.rec.name : 'template') + ' PCR ' + fmt(p.length) + ' bp',
              seq: p.seq, circular: false,
              features: transferFeatures(p.seq, [t.rec || t.from].filter(Boolean), false),
              history: [['PCR', 'Amplified ' + fmt(p.start + 1) + '–' + fmt(p.end) + ' from ' +
                (t.rec ? t.rec.name : 'pasted template') + '; primers ' + fwd + ' / ' + rev]]
            });
          }));
        });
      }
    });
  }

  /* ============================== SDM ============================== */
  function toolSDM() {
    const rec = PB.active();
    if (!rec) { PB.toast('Load a construct first.', 'warn'); return; }
    const pos = state.selection ? state.selection.start : 0;
    PB.modal({
      title: 'Site-directed mutagenesis', wide: true,
      body: '<div class="row">' +
          '<div class="field narrow"><label>Position (bp)</label><input class="inp mono" id="sdPos" type="number" min="1" max="' + rec.seq.length + '" value="' + (pos + 1) + '"></div>' +
          '<div class="field narrow"><label>Change</label><select class="inp" id="sdKind">' +
            '<option value="substitute">Substitute</option><option value="insert">Insert</option><option value="delete">Delete</option></select></div>' +
          '<div class="field"><label>New bases / length to delete</label><input class="inp mono" id="sdNew" placeholder="G, or GGATCC, or 3"></div>' +
          '<div class="field narrow"><label>Flank (bp)</label><input class="inp mono" id="sdFlank" type="number" value="15" min="8" max="30"></div>' +
          '<div class="field" style="flex:0 0 auto;align-self:flex-end"><button class="btn primary" id="sdGo">Design</button></div>' +
        '</div>' +
        '<div class="hint" id="sdCtx"></div>' + resultSlot(),
      onOpen: w => {
        const showCtx = () => {
          const p = Math.max(1, Math.min(rec.seq.length, parseInt($('#sdPos').value, 10) || 1)) - 1;
          const a = Math.max(0, p - 12), b = Math.min(rec.seq.length, p + 13);
          $('#sdCtx').innerHTML = 'Context: <span class="mono">' + esc(rec.seq.slice(a, p)) +
            '<b style="color:var(--bad)">' + esc(rec.seq[p] || '') + '</b>' + esc(rec.seq.slice(p + 1, b)) + '</span>';
        };
        w.querySelector('#sdPos').addEventListener('input', showCtx);
        showCtx();
        w.querySelector('#sdGo').addEventListener('click', () => {
          const p = Math.max(1, parseInt($('#sdPos').value, 10) || 1) - 1;
          const kind = $('#sdKind').value;
          const raw = $('#sdNew').value.trim();
          const opts = { position: p, kind, flank: Math.max(8, parseInt($('#sdFlank').value, 10) || 15) };
          if (kind === 'delete') opts.deleteLength = Math.max(1, parseInt(raw, 10) || 1);
          else {
            opts.newBases = B.cleanSeq(raw);
            if (!opts.newBases.length) { PB.toast('Enter the new bases.', 'warn'); return; }
          }
          opts.circular = rec.topology === 'circular';
          const r = B.sdm(rec.seq, opts);
          if (!r.ok) { result('<div class="note bad">' + esc(r.error) + '</div>'); return; }
          // the edit as "replace [ea, eb) with en bases", used to move every feature
          const L = rec.seq.length;
          const ea = p, eb = kind === 'insert' ? p : p + (kind === 'delete' ? opts.deleteLength : opts.newBases.length);
          const en = kind === 'delete' ? 0 : opts.newBases.length;
          const insB = kind === 'delete' ? '' : opts.newBases;
          /* features move through the same edit routine as every other edit,
             which handles base 1 of a circle correctly */
          const moveBy = list => B.editSequence(rec.seq, list, ea, eb, insB, opts.circular).features;
          const movedFeats = moveBy(rec.features);
          // protein-level consequence in every coding feature the edit touches
          const touches = f => B.featurePieces(f, L).some(pc => kind === 'insert'
            ? (ea > pc.start && ea < pc.end) : (ea < pc.end && eb > pc.start));
          const codes = (f) => { const t = f.quals && f.quals.transl_table; return t ? parseInt(t, 10) : state.codeTable; };
          const frame = f => { const c = f.quals && parseInt(f.quals.codon_start, 10); return c >= 1 && c <= 3 ? c - 1 : 0; };
          const effects = [];
          rec.features.forEach((f, i) => {
            if (!/^(CDS|ORF)$/i.test(f.type) || f.between || !touches(f)) return;
            const nf = moveBy([f])[0];
            const before = B.translate(B.featureSequence(rec.seq, f).slice(frame(f)), { table: codes(f) });
            const after = nf ? B.translate(B.featureSequence(r.mutated.seq, nf).slice(frame(f)), { table: codes(f) }) : '';
            const same = before === after;
            let firstDiff = -1;
            for (let k = 0; k < Math.max(before.length, after.length); k++) {
              if (before[k] !== after[k]) { firstDiff = k; break; }
            }
            const stopAt = after.indexOf('*');
            const premature = !same && stopAt >= 0 && stopAt < after.length - 1;
            const shift = kind !== 'substitute' && (eb - ea - en) % 3 !== 0;
            effects.push('<div class="note ' + (same ? 'ok' : 'warn') + '"><strong>In ' + esc(f.label) + ':</strong> ' +
              (same ? 'protein unchanged (silent).' :
                (firstDiff >= 0 ? 'first protein change at residue ' + (firstDiff + 1) + ': ' +
                  esc(before[firstDiff] || '–') + ' → ' + esc(after[firstDiff] || '–') : 'protein length changes') +
                (shift ? ' <strong>This shifts the reading frame.</strong>' : '') +
                (premature ? ' <strong>Creates a premature stop at residue ' + (stopAt + 1) + '.</strong>' : '')) +
              '</div>');
          });
          const effect = effects.join('');
          result('<div class="note">' + esc(r.change) + ' · product ' + fmt(r.mutated.length) + ' bp</div>' + effect +
            '<div class="row" style="margin-top:10px">' +
            '<div class="field" style="flex:1 1 300px"><label>Mutagenic forward primer · ' + r.primers.forward.length + ' nt</label>' +
            '<div class="seqblock" style="max-height:none">' + esc(r.primers.forward.seq) + '</div>' +
            '<div class="hint">Tm ' + f1(r.primers.forward.tm) + ' °C, mismatch centred</div></div>' +
            '<div class="field" style="flex:1 1 300px"><label>Mutagenic reverse primer · ' + r.primers.reverse.length + ' nt</label>' +
            '<div class="seqblock" style="max-height:none">' + esc(r.primers.reverse.seq) + '</div>' +
            '<div class="hint">Tm ' + f1(r.primers.reverse.tm) + ' °C</div></div></div>' +
            '<div class="row"><button class="btn primary" id="sdSave">Save mutant</button>' +
            '<button class="btn" id="sdSavePrim">Save primers</button></div>');
          $('#sdSave').addEventListener('click', () => {
            saveProduct({
              name: rec.name + ' ' + (kind === 'substitute' ? '' : kind + ' ') + '@' + (p + 1),
              seq: r.mutated.seq, circular: rec.topology === 'circular',
              features: movedFeats,
              history: [['Site-directed mutagenesis', r.change + '; primers ' + r.primers.forward.seq]]
            });
          });
          $('#sdSavePrim').addEventListener('click', () => {
            state.primers.push({ name: 'SDM_' + (p + 1) + '_F', seq: r.primers.forward.seq });
            state.primers.push({ name: 'SDM_' + (p + 1) + '_R', seq: r.primers.reverse.seq });
            PB.lsSet('primers', state.primers);
            PB.toast('Saved 2 primers', 'ok');
          });
        });
      }
    });
  }

  /* ============================== OVERLAP EXTENSION ============================== */
  function toolOverlap() {
    PB.modal({
      title: 'Overlap-extension fusion', wide: true,
      body: '<div class="row">' + seqField('oeA', 'Fragment A (upstream)', { allowSel: true, rows: 2 }) +
          seqField('oeB', 'Fragment B (downstream)', { allowSel: true, rows: 2 }) + '</div>' +
        '<div class="row"><div class="field narrow"><label>Overlap (bp)</label><input class="inp mono" id="oeOv" type="number" value="20" min="12" max="60"></div>' +
        '<div class="field" style="flex:0 0 auto;align-self:flex-end"><button class="btn primary" id="oeGo">Fuse</button></div></div>' +
        '<div class="hint">The chimeric inner primers carry the junction: A’s reverse primer gains B’s 5′ end, and B’s forward primer gains A’s 3′ end. The two products then anneal through the shared overlap and are extended into the fusion.</div>' +
        resultSlot(),
      onOpen: w => {
        w.querySelector('#oeGo').addEventListener('click', () => {
          const a = seqValue('oeA'), b = seqValue('oeB');
          if (!a.seq.length || !b.seq.length) { PB.toast('Provide both fragments.', 'warn'); return; }
          const ov = Math.max(12, parseInt($('#oeOv').value, 10) || 20);
          const fused = a.seq + b.seq;
          const halfA = a.seq.slice(-ov), halfB = b.seq.slice(0, ov);
          const innerR = B.revcomp(halfB + halfA.slice(-0 || undefined));
          const aRev = B.revcomp(halfB + a.seq.slice(-20));
          const bFwd = a.seq.slice(-ov) + b.seq.slice(0, 20);
          const cond = PB.conditions();
          result('<div class="note ok"><strong>Fusion product ' + fmt(fused.length) + ' bp.</strong> Junction overlap ' + ov + ' bp.</div>' +
            junctionBox('Junction',
              '<span style="color:var(--fg-faint)">' + esc(a.seq.slice(-16)) + '</span>' +
              '<span style="color:var(--accent);font-weight:600">' + esc(b.seq.slice(0, 16)) + '</span>') +
            '<div class="row" style="margin-top:10px">' +
            '<div class="field" style="flex:1 1 300px"><label>A reverse (chimeric) · ' + aRev.length + ' nt</label>' +
            '<div class="seqblock" style="max-height:none">' + esc(aRev) + '</div>' +
            '<div class="hint">Tm ' + f1(B.tmNN(aRev, cond)) + ' °C</div></div>' +
            '<div class="field" style="flex:1 1 300px"><label>B forward (chimeric) · ' + bFwd.length + ' nt</label>' +
            '<div class="seqblock" style="max-height:none">' + esc(bFwd) + '</div>' +
            '<div class="hint">Tm ' + f1(B.tmNN(bFwd, cond)) + ' °C</div></div></div>' +
            '<button class="btn primary" id="oeSave">Save fusion</button>');
          $('#oeSave').addEventListener('click', () => saveProduct({
            name: (a.name || 'A') + '–' + (b.name || 'B') + ' fusion (' + fmt(fused.length) + ' bp)',
            seq: fused, circular: false,
            features: transferFeatures(fused, [a.rec || a.from, b.rec || b.from].filter(Boolean), false),
            history: [['Overlap-extension fusion', ov + ' bp junction overlap']]
          }));
        });
      }
    });
  }

  /* ============================== SILENT MUTATIONS ============================== */
  function toolSilent() {
    const rec = PB.active();
    if (!rec) { PB.toast('Load a construct first.', 'warn'); return; }
    const cdsFeats = rec.features.filter(f => /CDS|gene|ORF/i.test(f.type));
    PB.modal({
      title: 'Silent mutations', wide: true, sub: 'change a site, keep the protein',
      body: '<div class="row">' +
          '<div class="field"><label>Coding region</label><select class="inp" id="smCds">' +
            (state.selection ? '<option value="__sel">current selection (' + fmt(state.selection.end - state.selection.start) + ' bp)</option>' : '') +
            cdsFeats.map((f, i) => '<option value="' + rec.features.indexOf(f) + '">' + esc(f.label) + ' · ' +
              fmt(PB.featureLen(rec, f)) + ' bp' + (f.parts ? ' (joined)' : '') + '</option>').join('') +
            '<option value="__all">whole sequence (frame 1)</option>' +
          '</select></div>' +
          enzField('smEnz', 'Enzyme', 'EcoRI') +
          '<div class="field narrow"><label>Action</label><select class="inp" id="smMode">' +
            '<option value="remove">Remove site</option><option value="add">Add site</option></select></div>' +
          '<div class="field" style="flex:0 0 auto;align-self:flex-end"><button class="btn primary" id="smGo">Search</button></div>' +
        '</div>' +
        '<div class="hint">Searches windows of one, two and three consecutive codons, so sites that straddle a codon boundary are found too. Every candidate is verified to leave the translated protein identical.</div>' +
        resultSlot(),
      onOpen: w => {
        w.querySelector('#smGo').addEventListener('click', () => {
          const names = enzNames('smEnz');
          if (!names.length) { PB.toast('Name an enzyme.', 'warn'); return; }
          const which = $('#smCds').value;
          const src = codingSource(rec, which);
          if (!src) { PB.toast('Pick a coding region.', 'warn'); return; }
          const cds = src.cds;
          if (cds.length < 6) { PB.toast('That region is too short.', 'warn'); return; }
          const mode = $('#smMode').value;
          const r = B.silentMutations(cds.slice(src.frame), names[0], { mode, table: src.table });
          if (!r.ok) { result('<div class="note bad">' + esc(r.error) + '</div>'); return; }
          let html;
          if (mode === 'remove') {
            if (!r.siteCount) {
              html = '<div class="note">There are no ' + esc(names[0]) + ' sites in that region, so nothing to remove.</div>';
            } else {
              html = '<div class="scrollx"><table class="tbl"><thead><tr><th class="num">Site at</th><th>Codons</th>' +
                '<th>Change</th><th>Residue</th><th class="num">Bases</th><th></th></tr></thead><tbody>' +
                r.results.map((x, i) => x.error
                  ? '<tr><td class="num">' + fmt(x.sitePos + 1) + '</td><td colspan="5"><span class="badge bad">' + esc(x.error) + '</span></td></tr>'
                  : '<tr><td class="num">' + fmt(x.sitePos + 1) + '</td>' +
                    '<td class="num">' + x.aaPosition + (x.codonSpan > 1 ? '–' + (x.aaPosition + x.codonSpan - 1) : '') + '</td>' +
                    '<td class="mono">' + esc(x.from) + ' → ' + esc(x.to) + '</td>' +
                    '<td class="mono">' + esc(x.aa) + '</td>' +
                    '<td class="num">' + x.changed.length + '</td>' +
                    '<td><button class="btn sm ghost" data-smap="' + i + '">Apply</button></td></tr>').join('') +
                '</tbody></table></div>';
            }
          } else {
            html = r.results.length
              ? '<div class="note ok">' + r.results.length + ' way' + (r.results.length === 1 ? '' : 's') +
                ' to introduce ' + esc(names[0]) + ' without changing the protein.</div>' +
                '<div class="scrollx scrolly"><table class="tbl"><thead><tr><th class="num">Residue</th><th>Change</th>' +
                '<th>Residues</th><th class="num">Bases changed</th><th></th></tr></thead><tbody>' +
                r.results.slice(0, 40).map((x, i) => '<tr><td class="num">' + x.aaPosition + '</td>' +
                  '<td class="mono">' + esc(x.from) + ' → ' + esc(x.to) + '</td>' +
                  '<td class="mono">' + esc(x.aa) + '</td>' +
                  '<td class="num">' + x.changed.length + '</td>' +
                  '<td><button class="btn sm ghost" data-smap="' + i + '">Apply</button></td></tr>').join('') +
                '</tbody></table></div>'
              : '<div class="note warn">No single-, double- or triple-codon synonymous change creates an ' + esc(names[0]) +
                ' site in that region.</div>';
          }
          result(html);
          $('#toolResult').querySelectorAll('[data-smap]').forEach(b => b.addEventListener('click', () => {
            const x = r.results[+b.getAttribute('data-smap')];
            if (!x || !x.seq) return;
            // write the edited bases back through the feature's own geometry
            const newSeq = src.apply(cds.slice(0, src.frame) + x.seq);
            const chk = verifyRecode(rec, newSeq, [names[0]]);
            if (chk.proteinChanged.length) {
              PB.toast('Refused: this edit would change ' + chk.proteinChanged.join(', ') + '.', 'bad'); return;
            }
            const site = chk.sites[names[0]];
            saveProduct({
              name: rec.name + ' ' + (mode === 'remove' ? '−' : '+') + names[0],
              seq: newSeq, circular: rec.topology === 'circular',
              features: rec.features.map(f => Object.assign({}, f)),
              history: [['Silent mutation',
                (mode === 'remove' ? 'Removed ' : 'Added ') + names[0] + ' site: ' + x.from + ' → ' + x.to +
                ' (' + x.aa + ' at residue ' + x.aaPosition + '), ' + x.changed.length + ' base(s) changed, protein unchanged; ' +
                names[0] + ' sites in construct ' + site.before + ' → ' + site.after]]
            });
            if (mode === 'remove' && site.after > 0) PB.toast(names[0] + ' still cuts the new construct ' + site.after + '× elsewhere.', 'warn');
          }));
        });
      }
    });
  }

  /* ============================== CODON OPTIMISATION ============================== */
  function toolCodon() {
    const rec = PB.active();
    if (!rec) { PB.toast('Load a construct first.', 'warn'); return; }
    const cdsFeats = rec.features.filter(f => /CDS|gene|ORF/i.test(f.type));
    PB.modal({
      title: 'Codon optimisation', wide: true,
      body: '<div class="row">' +
          '<div class="field"><label>Coding region</label><select class="inp" id="coCds">' +
            (state.selection ? '<option value="__sel">current selection</option>' : '') +
            cdsFeats.map(f => '<option value="' + rec.features.indexOf(f) + '">' + esc(f.label) + '</option>').join('') +
          '</select></div>' +
          '<div class="field narrow"><label>Host</label><select class="inp" id="coHost">' +
            '<option value="ecoli">E. coli K-12</option><option value="human">Human</option></select></div>' +
          enzField('coAvoid', 'Do not create sites for', 'BsaI, BsmBI, EcoRI, BamHI, XhoI, NotI') +
          '<div class="field" style="flex:0 0 auto;align-self:flex-end"><button class="btn primary" id="coGo">Optimise</button></div>' +
        '</div>' + resultSlot(),
      onOpen: w => {
        w.querySelector('#coGo').addEventListener('click', () => {
          const which = $('#coCds').value;
          const src = codingSource(rec, which);
          if (!src) { PB.toast('Pick a coding region.', 'warn'); return; }
          const cds = src.cds;
          if (cds.length < 6) { PB.toast('That region is too short.', 'warn'); return; }
          const usage = $('#coHost').value === 'human' ? B.HUMAN_USAGE : B.ECOLI_USAGE;
          const avoid = enzNames('coAvoid');
          // the frame offset (/codon_start) is kept verbatim in front
          const head = cds.slice(0, src.frame), body = cds.slice(src.frame);
          const opt = head + B.codonOptimize(body, usage, src.table, { keepFirst: src.isFeature && src.frame === 0, avoid });
          const protIn = B.translate(body, { table: src.table });
          const protOut = B.translate(opt.slice(src.frame), { table: src.table });
          let changed = 0;
          for (let i = 0; i < Math.min(cds.length, opt.length); i++) if (cds[i] !== opt[i]) changed++;
          const newSeq = opt.length === cds.length ? src.apply(opt) : null;
          const chk = newSeq ? verifyRecode(rec, newSeq, avoid) : null;
          const same = protIn === protOut && newSeq && !chk.proteinChanged.length;
          const created = chk ? avoid.filter(n => chk.sites[n].after > chk.sites[n].before) : [];
          result('<div class="note ' + (same ? 'ok' : 'bad') + '"><strong>' +
            (same ? 'Protein preserved exactly' + (src.isFeature ? ' (checked on every coding feature in the construct).' : '.') : 'Protein changed — do not use this result.') + '</strong><br>' +
            changed + ' of ' + cds.length + ' bases recoded (' + f1(changed / cds.length * 100) + '%). ' +
            'GC ' + f1(B.gcContent(cds)) + '% → ' + f1(B.gcContent(opt)) + '%.' +
            (src.isFeature && src.frame === 0 ? ' Start codon left unchanged.' : '') + '</div>' +
            (avoid.length ? '<div class="note ' + (created.length ? 'warn' : 'ok') + '">' +
              avoid.map(n => esc(n) + ': ' + chk.sites[n].before + ' → ' + chk.sites[n].after + ' site' + (chk.sites[n].after === 1 ? '' : 's')).join(' · ') +
              (created.length ? ' — a new site sits where the region meets its flanking sequence; adjust by hand.' : '') + '</div>' : '') +
            '<div class="hint">Greedy most-frequent codon for each residue. Usage tables are genome-wide averages (Kazusa); this does not model mRNA structure, repeats or rare-codon clusters.</div>' +
            '<div class="row"><div class="field wide"><label>Optimised sequence</label><div class="seqblock">' + esc(opt) + '</div></div></div>' +
            '<div class="row"><div class="field wide"><label>Protein</label><div class="seqblock">' + esc(protOut) + '</div></div></div>' +
            (same ? '<button class="btn primary" id="coSave">Save optimised construct</button>' : ''));
          const s = $('#coSave');
          if (s) s.addEventListener('click', () => {
            saveProduct({
              name: rec.name + ' codon-optimised',
              seq: newSeq, circular: rec.topology === 'circular',
              features: rec.features.map(f => Object.assign({}, f)),
              history: [['Codon optimisation', 'Host ' + $('#coHost').value + '; ' + changed + ' bases recoded, protein identical' +
                (avoid.length ? '; avoided ' + avoid.join(', ') : '')]]
            });
          });
        });
      }
    });
  }

  /* ============================== BACK-TRANSLATION ============================== */
  function toolBackTranslate() {
    PB.modal({
      title: 'Back-translate a protein', wide: true,
      body: '<div class="row"><div class="field wide"><label>Protein (one-letter code)</label>' +
          '<textarea class="inp mono" id="btProt" rows="5" spellcheck="false" placeholder="MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"></textarea></div></div>' +
        '<div class="row">' +
          '<div class="field"><label>Name</label><input class="inp" id="btName" value="back-translated CDS"></div>' +
          '<div class="field narrow"><label>Host</label><select class="inp" id="btHost"><option value="ecoli">E. coli K-12</option><option value="human">Human</option></select></div>' +
          enzField('btAvoid', 'Do not create sites for', 'BsaI, BsmBI, EcoRI, BamHI, XhoI, NotI') +
          '<div class="field narrow"><label>Stop codon</label><select class="inp" id="btStop"><option value="1">add TAA</option><option value="0">none</option></select></div>' +
        '</div>' +
        '<div class="hint">Uses the standard genetic code (or the code chosen in Conditions). Every residue gets the host’s most-used codon, then codons are swapped wherever that would create a site for the listed enzymes. The result is translated back and must match your protein exactly before it can be saved.</div>' +
        resultSlot(),
      actions: [{ label: 'Close' }, { label: 'Back-translate', primary: true, onClick: () => {
        const raw = ($('#btProt').value || '').replace(/^>.*$/gm, '').replace(/[\s\d]/g, '').toUpperCase();
        const prot = raw.replace(/\*+$/, '');
        const table = state.codeTable;
        const bad = Array.from(new Set(prot.replace(/[ACDEFGHIKLMNPQRSTVWY]/g, '').split(''))).join('');
        if (!prot.length) { PB.toast('Paste a protein sequence.', 'warn'); return false; }
        if (bad) { result('<div class="note bad">Not a standard amino-acid code: <span class="mono">' + esc(bad) + '</span>. Remove or replace these first.</div>'); return false; }
        let draft = '';
        for (const aa of prot) {
          const syn = B.synonymousCodons(aa, table);
          if (!syn.length) { result('<div class="note bad">Residue ' + esc(aa) + ' has no codon in genetic code ' + table + '.</div>'); return false; }
          draft += syn[0];
        }
        const usage = $('#btHost').value === 'human' ? B.HUMAN_USAGE : B.ECOLI_USAGE;
        const avoid = enzNames('btAvoid');
        let cds = B.codonOptimize(draft, usage, table, { avoid });
        if (prot[0] === 'M' && cds.slice(0, 3) !== 'ATG') cds = 'ATG' + cds.slice(3);
        if ($('#btStop').value === '1') cds += 'TAA';
        const back = B.translate(cds, { table }).replace(/\*$/, '');
        if (back !== prot) { result('<div class="note bad">Internal check failed: the DNA does not translate back to your protein. Nothing was saved.</div>'); return false; }
        const left = avoid.map(n => { const e = PB.enzByName(n); return e ? [e.name, B.findSites(cds, e, false).length] : null; })
          .filter(x => x && x[1]);
        const name = ($('#btName').value || 'back-translated CDS').trim();
        const nr = saveProduct({
          name, seq: cds, circular: false,
          features: [PB.normFeature({ type: 'CDS', start: 0, end: cds.length, strand: 1, label: name,
            quals: { translation: prot, codon_start: '1', transl_table: String(table) } })],
          history: [['Back-translated', prot.length + ' aa for ' + $('#btHost').value + ' codon usage' +
            (avoid.length ? '; avoided ' + avoid.join(', ') : '') + '; translation verified identical']]
        });
        if (left.length) PB.toast('Could not avoid every site: ' + left.map(x => x[0] + ' ×' + x[1]).join(', ') + '.', 'warn');
        return nr ? undefined : false;
      } }]
    });
  }

  /* ============================== CRISPR ============================== */
  function toolCRISPR() {
    const rec = PB.active();
    if (!rec) { PB.toast('Load a construct first.', 'warn'); return; }
    PB.modal({
      title: 'CRISPR guide design', wide: true,
      body: '<div class="row">' +
          '<div class="field narrow"><label>PAM</label><select class="inp" id="crPam">' +
            '<option value="NGG">NGG — SpCas9</option><option value="NAG">NAG — SpCas9 (weak)</option>' +
            '<option value="NGA">NGA — SpCas9-VQR</option><option value="NNGRRT">NNGRRT — SaCas9</option>' +
            '<option value="TTTV">TTTV — Cas12a (5′ PAM)</option></select></div>' +
          '<div class="field narrow"><label>Guide length</label><input class="inp mono" id="crLen" type="number" value="20" min="16" max="24"></div>' +
          '<div class="field narrow"><label>Region</label><select class="inp" id="crRegion">' +
            (state.selection ? '<option value="sel">selection</option>' : '') +
            '<option value="all">whole construct</option></select></div>' +
          '<div class="field" style="flex:0 0 auto;align-self:flex-end"><button class="btn primary" id="crGo">Find guides</button></div>' +
        '</div>' +
        '<div class="note">This finds every protospacer with the chosen PAM and reports GC content and whether the guide would be cut by the usual cloning enzymes. It does <strong>not</strong> score off-targets against a genome — that needs a genome index this page does not have. Check candidates in a genome-aware tool before ordering.</div>' +
        resultSlot(),
      onOpen: w => {
        w.querySelector('#crGo').addEventListener('click', () => {
          const pam = $('#crPam').value;
          const len = Math.max(16, parseInt($('#crLen').value, 10) || 20);
          const region = $('#crRegion').value;
          const seq = (region === 'sel' && state.selection)
            ? PB.selectionSeq(rec, state.selection) : rec.seq;
          const base = (region === 'sel' && state.selection) ? state.selection.start : 0;
          if (seq.length < len + pam.length) { PB.toast('That region is too short.', 'warn'); return; }
          let guides = B.findGuides(seq, { pam, length: len, circular: region !== 'sel' && rec.topology === 'circular' });
          guides = guides.filter(g => g.gc >= 20 && g.gc <= 85);
          guides.sort((a, b) => Math.abs(a.gc - 50) - Math.abs(b.gc - 50));
          if (!guides.length) { result('<div class="note warn">No guides with that PAM in this region.</div>'); return; }
          result('<div class="note ok">' + guides.length + ' candidate guide' + (guides.length === 1 ? '' : 's') +
            ', sorted by how close the GC content is to 50%.</div>' +
            '<div class="scrollx scrolly"><table class="tbl"><thead><tr><th class="num">Position</th><th>Strand</th>' +
            '<th>Protospacer</th><th>PAM</th><th class="num">GC</th><th>Notes</th><th></th></tr></thead><tbody>' +
            guides.slice(0, 200).map((g, i) => '<tr>' +
              '<td class="num">' + fmt((base + g.start) % rec.seq.length + 1) + '</td>' +
              '<td>' + (g.strand === 1 ? '+' : '−') + '</td>' +
              '<td class="mono">' + esc(g.protospacer) + '</td>' +
              '<td class="mono">' + esc(g.pam) + '</td>' +
              '<td class="num">' + f1(g.gc) + '</td>' +
              '<td>' + (g.wrapped ? '<span class="badge neutral">crosses the origin</span> ' : '') +
                (g.hasBsaI ? '<span class="badge warn">contains BsaI site</span>' : '') +
                (/^G/.test(g.protospacer) ? '' : ' <span class="badge neutral">no 5′ G</span>') + '</td>' +
              '<td><button class="btn sm ghost" data-crcopy="' + i + '">⧉</button>' +
                '<button class="btn sm ghost" data-crfeat="' + i + '">+</button></td>' +
              '</tr>').join('') + '</tbody></table></div>');
          $('#toolResult').querySelectorAll('[data-crcopy]').forEach(b => b.addEventListener('click', () =>
            PB.copyText(guides[+b.getAttribute('data-crcopy')].protospacer, 'Protospacer')));
          $('#toolResult').querySelectorAll('[data-crfeat]').forEach(b => b.addEventListener('click', () => {
            const g = guides[+b.getAttribute('data-crfeat')];
            rec.features.push(PB.normFeature({
              type: 'misc_feature', start: base + g.start, end: base + g.start + len,
              strand: g.strand, label: 'gRNA ' + g.protospacer.slice(0, 8) + '…'
            }));
            PB.addHistory(rec, 'gRNA annotated', g.protospacer + ' ' + g.pam);
            PB.touch(rec); PB.renderAll();
            PB.toast('Guide annotated', 'ok');
          }));
        });
      }
    });
  }

  /* ============================== RNA STRUCTURE ============================== */
  function toolRNA() {
    const rec = PB.active();
    const sel = state.selection;
    PB.modal({
      title: 'RNA secondary structure', wide: true,
      body: '<div class="row">' + seqField('rnSeq', 'RNA / DNA region', { allowSel: true, rows: 3, defaultId: sel ? undefined : (rec ? rec.id : '') }) + '</div>' +
        '<div class="row"><div class="field narrow"><label>Min loop</label><input class="inp mono" id="rnLoop" type="number" value="3" min="3" max="8"></div>' +
        '<div class="field" style="flex:0 0 auto;align-self:flex-end"><button class="btn primary" id="rnGo">Fold</button></div></div>' +
        '<div class="note">Nussinov base-pair maximisation. It finds a maximum-pairing structure, which is <strong>not</strong> the same as a minimum-free-energy prediction — use ViennaRNA or mfold when the energetics matter. Limit 1,200 nt.</div>' +
        resultSlot(),
      onOpen: w => {
        w.querySelector('#rnGo').addEventListener('click', () => {
          const v = seqValue('rnSeq');
          if (v.seq.length < 8) { PB.toast('Provide at least 8 bases.', 'warn'); return; }
          if (v.seq.length > 1200) { PB.toast('Limit is 1,200 nt — select a smaller region.', 'bad'); return; }
          const r = B.rnaFold(v.seq, { minLoop: parseInt($('#rnLoop').value, 10) || 3 });
          if (r.error) { result('<div class="note bad">' + esc(r.error) + '</div>'); return; }
          const blocks = [];
          for (let i = 0; i < r.sequence.length; i += 60) {
            blocks.push('<span class="lbl">' + String(i + 1).padStart(6) + ' </span>' + esc(r.sequence.substr(i, 60)) + '\n' +
              '<span class="lbl">       </span>' + esc(r.structure.substr(i, 60)) + '\n\n');
          }
          result('<div class="note ok"><strong>' + r.pairCount + ' base pairs</strong> over ' + r.sequence.length +
            ' nt (' + f1(r.pairCount * 2 / r.sequence.length * 100) + '% of bases paired).</div>' +
            '<div class="aln-box">' + blocks.join('') + '</div>' +
            '<div class="hint">' + esc(r.method) + '</div>');
        });
      }
    });
  }

  /* ============================== INTRONS ============================== */
  function toolIntrons() {
    const rec = PB.active();
    if (!rec) { PB.toast('Load a construct first.', 'warn'); return; }
    const introns = rec.features.filter(f => /intron/i.test(f.type) || /intron/i.test(f.label));
    PB.modal({
      title: 'Intron splicing', wide: true,
      body: '<div class="hint">Removes the regions listed below and returns the processed transcript. Annotate introns as features (type <span class="mono">intron</span>), or enter ranges by hand.</div>' +
        '<div class="row"><div class="field wide"><label>Intron ranges (start-end, one per line, 1-based inclusive)</label>' +
        '<textarea class="inp" id="inRanges" rows="5" spellcheck="false">' +
        esc(introns.map(f => (f.start + 1) + '-' + f.end).join('\n')) + '</textarea></div></div>' +
        '<div class="row"><button class="btn primary" id="inGo">Splice</button>' +
        '<button class="btn sm" id="inScan">Suggest GT…AG introns</button></div>' + resultSlot(),
      onOpen: w => {
        w.querySelector('#inScan').addEventListener('click', () => {
          // canonical donor GT ... acceptor AG, 60-10000 nt
          const s = rec.seq;
          const found = [];
          for (let i = 0; i + 60 < s.length && found.length < 40; i++) {
            if (s[i] !== 'G' || s[i + 1] !== 'T') continue;
            for (let j = i + 60; j < Math.min(s.length, i + 10000); j++) {
              if (s[j - 1] === 'A' && s[j] === 'G') { found.push([i + 1, j + 1]); i = j; break; }
            }
          }
          $('#inRanges').value = found.map(r => r[0] + '-' + r[1]).join('\n');
          PB.toast(found.length + ' candidate GT…AG intervals — these are sequence motifs only, not predictions.', 'warn');
        });
        w.querySelector('#inGo').addEventListener('click', () => {
          const ranges = ($('#inRanges').value || '').split(/\n+/).map(l => {
            const m = l.match(/(\d+)\s*[-–.]+\s*(\d+)/);
            return m ? { start: +m[1] - 1, end: +m[2] } : null;
          }).filter(Boolean).filter(r => r.end > r.start);
          if (!ranges.length) { PB.toast('Enter at least one range.', 'warn'); return; }
          const mrna = B.spliceIntrons(rec.seq, ranges);
          const removed = rec.seq.length - mrna.length;
          result('<div class="note ok"><strong>' + ranges.length + ' intron' + (ranges.length === 1 ? '' : 's') +
            ' removed</strong> (' + fmt(removed) + ' bp). Processed length ' + fmt(mrna.length) + ' bp.</div>' +
            '<button class="btn primary" id="inSave">Save processed transcript</button>');
          $('#inSave').addEventListener('click', () => saveProduct({
            name: rec.name + ' mRNA',
            seq: mrna, circular: false,
            history: [['Intron splicing', ranges.length + ' introns removed (' + fmt(removed) + ' bp): ' +
              ranges.map(r => (r.start + 1) + '-' + r.end).join(', ')]]
          }));
        });
      }
    });
  }

  /* ============================== TRANSLATE ============================== */
  PB.toolTranslate = function () {
    const rec = PB.active();
    const sel = state.selection;
    const seq = (rec && PB.normSel(rec, sel)) ? PB.selectionSeq(rec, sel) : (rec ? rec.seq.slice(0, 3000) : '');
    PB.modal({
      title: 'Translate', wide: true,
      body: '<div class="row"><div class="field wide"><label>Sequence</label>' +
        '<textarea class="inp" id="trSeq" rows="4" spellcheck="false">' + esc(seq) + '</textarea></div></div>' +
        '<div class="row"><div class="field narrow"><label>Frame</label><select class="inp" id="trFrame">' +
          '<option value="0">+1</option><option value="1">+2</option><option value="2">+3</option>' +
          '<option value="-0">−1</option><option value="-1">−2</option><option value="-2">−3</option></select></div>' +
        '<div class="field"><label>Genetic code</label><select class="inp" id="trCode">' +
          B.geneticCodes().map(g => '<option value="' + g.id + '"' + (g.id === state.codeTable ? ' selected' : '') + '>' +
            esc(g.id + ' · ' + g.name) + '</option>').join('') + '</select></div>' +
        '<div class="field narrow" style="flex:0 0 120px;align-self:flex-end"><label class="chip" id="trStop">Stop at first *</label></div>' +
        '<div class="field" style="flex:0 0 auto;align-self:flex-end"><button class="btn primary" id="trGo">Translate</button></div></div>' +
        resultSlot(),
      onOpen: w => {
        let toStop = false;
        const sEl = w.querySelector('#trStop');
        sEl.addEventListener('click', () => { toStop = !toStop; sEl.classList.toggle('on', toStop); });
        const go = () => {
          const raw = B.cleanSeq($('#trSeq').value);
          if (raw.length < 3) { PB.toast('Need at least 3 bases.', 'warn'); return; }
          const fv = $('#trFrame').value;
          const rev = fv.charAt(0) === '-';
          const frame = Math.abs(parseInt(rev ? fv.slice(1) : fv, 10)) || 0;
          const table = parseInt($('#trCode').value, 10) || 1;
          const src = rev ? B.revcomp(raw) : raw;
          const prot = B.translate(src, { table, frame, toStop });
          const mw = prot.replace(/\*/g, '').length * 110;
          result('<div class="row"><div class="field wide"><label>Protein · ' +
            prot.replace(/\*+$/, '').length + ' aa, about ' + fmt(Math.round(mw / 100) / 10) + ' kDa</label>' +
            '<div class="seqblock">' + esc(prot) + '</div></div></div>' +
            '<div class="row"><div class="field wide"><label>Three-letter</label><div class="seqblock">' +
            esc(prot.split('').map(a => B.AA3[a] || 'Xaa').join(' ')) + '</div></div></div>' +
            (prot.slice(0, -1).indexOf('*') !== -1 ? '<div class="note warn">There is a stop codon inside this reading frame.</div>' : '') +
            (raw.length % 3 ? '<div class="note warn">The sequence length is not a multiple of 3; the last partial codon is ignored.</div>' : '') +
            '<button class="btn sm" id="trCopy">Copy protein</button>');
          $('#trCopy').addEventListener('click', () => PB.copyText(prot, 'Protein'));
        };
        w.querySelector('#trGo').addEventListener('click', go);
        go();
      }
    });
  };

  /* ============================== EDIT TOOLS ============================== */
  function toolEditSeq() {
    const rec = PB.active();
    PB.modal({
      title: rec ? 'Edit sequence' : 'New sequence', wide: true,
      body: '<div class="row"><div class="field"><label>Name</label><input class="inp" id="esName" value="' +
          esc(rec ? rec.name : '') + '"></div>' +
          '<div class="field narrow"><label>Topology</label><select class="inp" id="esTopo">' +
          '<option value="circular"' + (rec && rec.topology === 'circular' ? ' selected' : '') + '>Circular</option>' +
          '<option value="linear"' + (!rec || rec.topology === 'linear' ? ' selected' : '') + '>Linear</option></select></div></div>' +
        '<div class="row"><div class="field wide"><label>Sequence, FASTA or GenBank</label>' +
        '<textarea class="inp" id="esSeq" rows="10" spellcheck="false">' + esc(rec ? rec.seq : '') + '</textarea></div></div>' +
        '<div class="hint" id="esCount"></div>',
      actions: [{ label: 'Cancel' }, { label: rec ? 'Replace sequence' : 'Create', primary: true, onClick: () => {
        const raw = $('#esSeq').value;
        let seq, feats = null, topo = $('#esTopo').value;
        if (/LOCUS|ORIGIN/.test(raw)) {
          const recs = B.parseGenBank(raw);
          if (!recs.length || !recs[0].seq) { PB.toast('Could not read that GenBank text.', 'bad'); return false; }
          seq = recs[0].seq; feats = recs[0].features; topo = recs[0].topology;
        } else {
          seq = B.cleanSeq(raw.replace(/^>.*$/gm, ''));
        }
        if (!seq.length) { PB.toast('No valid bases found.', 'bad'); return false; }
        const name = ($('#esName').value || 'Untitled').trim();
        if (rec && !rec.builtin) {
          rec.name = name; rec.seq = seq; rec.topology = topo;
          if (feats) rec.features = feats.map(PB.normFeature);
          rec.features = rec.features.filter(f => f.end <= seq.length);
          PB.addHistory(rec, 'Sequence replaced', fmt(seq.length) + ' bp');
          PB.touch(rec); PB.invalidateSites(); PB.renderAll();
        } else {
          const nr = PB.makeRecord({ name, seq, topology: topo, features: feats ? feats.map(PB.normFeature) : [] });
          PB.addHistory(nr, 'Created', fmt(seq.length) + ' bp ' + topo);
          PB.addRecord(nr);
        }
      } }],
      onOpen: w => {
        const ta = w.querySelector('#esSeq'), cnt = w.querySelector('#esCount');
        const upd = () => {
          const s = B.cleanSeq(ta.value.replace(/^>.*$/gm, '').replace(/^(LOCUS|DEFINITION|FEATURES|ORIGIN|\/\/).*$/gm, ''));
          cnt.textContent = fmt(s.length) + ' bases recognised · ' + f1(B.gcContent(s)) + '% GC';
        };
        ta.addEventListener('input', upd); upd();
      }
    });
  }

  /* Replace bases a..b of the active construct with `ins`, in place, as one
     undoable step. b may exceed the length on a circle (the range carries on
     across base 1). With opts.asCopy the result is saved as a new construct
     and the original is left alone. Returns a short description, or null
     when nothing was done. */
  PB.editInPlace = function (a, b, ins, opts) {
    opts = opts || {};
    const rec = PB.active();
    if (!rec) { PB.toast('Load a construct first.', 'warn'); return null; }
    const L = rec.seq.length, circ = rec.topology === 'circular';
    ins = B.cleanSeq(ins || '');
    if (b === a && !ins.length) { PB.toast('Nothing to insert.', 'warn'); return null; }
    let ed;
    try { ed = B.editSequence(rec.seq, rec.features, a, b, ins, circ); }
    catch (e) { PB.toast(e.message, 'bad'); return null; }
    const span = b > L ? (a + 1) + '–' + L + ' and 1–' + (b - L) : (a + 1) + '–' + b;
    const detail = b === a ? 'Inserted ' + fmt(ins.length) + ' bp before base ' + (a + 1)
      : !ins.length ? 'Deleted ' + fmt(b - a) + ' bp at ' + span
      : 'Replaced ' + fmt(b - a) + ' bp at ' + span + ' with ' + fmt(ins.length) + ' bp';
    const lost = rec.features.length - ed.features.length;
    const note = (ed.wrapped ? ' The edit crossed base 1, so the circle now starts at the first base after it.' : '') +
      (lost > 0 ? ' ' + lost + ' feature' + (lost === 1 ? ' lay' : 's lay') + ' entirely inside the removed bases and ' +
        (lost === 1 ? 'was' : 'were') + ' dropped.' : '');
    if (opts.asCopy || rec.builtin) {
      const nr = PB.makeRecord({ name: rec.name + ' edited', seq: ed.seq, topology: rec.topology,
                                 features: ed.features, methylation: rec.methylation });
      PB.addHistory(nr, 'Sequence edited', detail + ' (from ' + rec.name + ')');
      PB.addRecord(nr);
      PB.toast(detail + ' — saved as a new construct.' + note, lost > 0 ? 'warn' : 'ok');
      return detail;
    }
    rec.seq = ed.seq; rec.features = ed.features;
    PB.addHistory(rec, 'Sequence edited', detail);
    /* show what changed: the new bases, if any */
    const at = ed.wrapped ? ed.seq.length - ins.length : a;
    state.selection = ins.length ? { start: at, end: at + ins.length } : null;
    state.selectedFeature = null;
    PB.touch(rec); PB.invalidateSites(); PB.renderAll();
    PB.toast(detail + '.' + note + ' Ctrl/⌘ Z undoes it.', lost > 0 ? 'warn' : 'ok');
    return detail;
  };

  function toolInsertDelete(prefill) {
    const rec = PB.active();
    if (!rec) { PB.toast('Load a construct first.', 'warn'); return; }
    prefill = prefill || {};
    const L = rec.seq.length, circ = rec.topology === 'circular';
    const sel = PB.normSel(rec, state.selection);
    const from = prefill.from || ((sel ? sel.start : 0) + 1);
    const to = prefill.to || (sel ? (sel.end - 1) % L + 1 : from);
    PB.modal({
      title: 'Edit the sequence', wide: true,
      body: '<div class="row">' +
          '<div class="field narrow"><label>Action</label><select class="inp" id="idKind">' +
          ['insert', 'delete', 'replace'].map(k => '<option value="' + k + '"' + (prefill.kind === k ? ' selected' : '') + '>' +
            { insert: 'Insert before', delete: 'Delete range', replace: 'Replace range' }[k] + '</option>').join('') +
          '</select></div>' +
          '<div class="field narrow"><label>From (bp)</label><input class="inp mono" id="idFrom" type="number" min="1" max="' + (L + 1) + '" value="' + from + '"></div>' +
          '<div class="field narrow" id="idToWrap"><label>To (bp)</label><input class="inp mono" id="idTo" type="number" min="1" max="' + L + '" value="' + to + '"></div>' +
        '</div><div class="row" id="idSeqRow"><div class="field wide"><label>New bases</label>' +
        '<textarea class="inp" id="idSeq" rows="3" spellcheck="false">' + esc(prefill.seq || '') + '</textarea></div></div>' +
        '<div class="row"><label class="chk"><input type="checkbox" id="idCopy"' + (rec.builtin ? ' checked disabled' : '') +
          '> Keep the original — save the result as a new construct</label></div>' +
        '<div class="hint" id="idPreview"></div>' +
        '<div class="hint">Insert puts the new bases <em>before</em> base “From”; use one past the last base to add them at the end. ' +
        (circ ? 'On this circular construct, a From greater than To runs the range across base 1. ' : '') +
        'Features after the edit move with it. A feature the edit runs through takes the new bases; ' +
        'one the edit only clips gets shorter; one inside the removed bases is dropped. ' +
        'The edit is made in place and Ctrl/⌘ Z undoes it.' +
        (rec.builtin ? ' The demo construct itself never changes, so the result is saved as a new construct.' : '') + '</div>',
      actions: [{ label: 'Cancel' }, { label: 'Apply', primary: true, onClick: () => {
        const r = readRange();
        if (r.error) { PB.toast(r.error, 'bad'); return false; }
        const done = PB.editInPlace(r.a, r.b, r.ins, { asCopy: $('#idCopy').checked });
        if (!done) return false;
      } }],
      onOpen: w => {
        const upd = () => {
          const kind = $('#idKind').value;
          w.querySelector('#idToWrap').style.display = kind === 'insert' ? 'none' : '';
          w.querySelector('#idSeqRow').style.display = kind === 'delete' ? 'none' : '';
          const r = readRange();
          w.querySelector('#idPreview').innerHTML = r.error ? '<span style="color:var(--bad)">' + esc(r.error) + '</span>'
            : esc(r.text);
        };
        ['#idKind', '#idFrom', '#idTo', '#idSeq'].forEach(id => {
          w.querySelector(id).addEventListener('input', upd); w.querySelector(id).addEventListener('change', upd);
        });
        upd();
      }
    });
    function readRange() {
      const kind = $('#idKind').value;
      const f = parseInt($('#idFrom').value, 10), t = parseInt($('#idTo').value, 10);
      const ins = B.cleanSeq($('#idSeq').value || '');
      if (!(f >= 1 && f <= L + (kind === 'insert' ? 1 : 0))) return { error: 'From must be between 1 and ' + (kind === 'insert' ? L + 1 : L) + '.' };
      const a = f - 1;
      if (kind === 'insert') {
        if (!ins.length) return { error: 'Enter the bases to insert.' };
        return { a: a, b: a, ins: ins, text: fmt(ins.length) + ' bp go in before base ' + f +
          (f === L + 1 ? ' (after the last base)' : '') + '; the construct becomes ' + fmt(L + ins.length) + ' bp.' };
      }
      if (!(t >= 1 && t <= L)) return { error: 'To must be between 1 and ' + L + '.' };
      let b = t;
      if (t < f) {
        if (!circ) return { error: 'To is before From. On a linear construct the range cannot cross base 1.' };
        b = t + L;
      }
      if (b - a >= L) return { error: 'That range is the whole construct.' };
      if (kind === 'replace' && !ins.length) return { error: 'Enter the replacement bases, or choose Delete.' };
      const len = b - a, newL = L - len + (kind === 'delete' ? 0 : ins.length);
      return { a: a, b: b, ins: kind === 'delete' ? '' : ins,
        text: (kind === 'delete' ? 'Removes ' : 'Replaces ') + fmt(len) + ' bp (' + f + '–' + t +
          (b > L ? ', across base 1' : '') + ')' + (kind === 'replace' ? ' with ' + fmt(ins.length) + ' bp' : '') +
          '; the construct becomes ' + fmt(newL) + ' bp.' };
    }
  }

  /* ---------------- search the whole library ----------------
     One box for everything: bases (IUPAC codes allowed, both strands, and
     across base 1 of a circle), a peptide, an enzyme by name or by its
     recognition site, a feature label or note, a construct name, a primer.
     Every reading of the query that finds something is shown, grouped by
     construct, and a hit can be opened where it sits. */
  PB.searchDialog = function (prefill) {
    const q0 = prefill != null ? prefill : (state.lastQuery || '');
    PB.modal({
      title: 'Search the library', wide: true,
      body: '<div class="row"><div class="field wide"><label>Find</label>' +
        '<input class="inp mono" id="lsQ" spellcheck="false" autocomplete="off" ' +
        'placeholder="bases (GGWTCC), a peptide (MEDYTK), an enzyme (BamHI), a label (kanR), a primer name" value="' +
        esc(q0) + '"></div>' +
        '<div class="field" style="flex:0 0 200px"><label>Look in</label><select class="inp" id="lsScope">' +
        '<option value="all">Every construct</option><option value="one">This construct only</option>' +
        '</select></div></div>' +
        '<div id="lsKinds" class="chipbar" style="margin-top:2px"></div>' +
        '<div id="lsOut" style="margin-top:8px;max-height:min(58vh,520px);overflow:auto"></div>',
      actions: [{ label: 'Close' }],
      onOpen: w => {
        const box = w.querySelector('#lsQ');
        const run = () => {
          const q = box.value.trim();
          state.lastQuery = q;
          const rec = PB.active();
          const recs = $('#lsScope').value === 'one' ? (rec ? [rec] : []) : state.records;
          const kindsEl = w.querySelector('#lsKinds'), out = w.querySelector('#lsOut');
          if (q.length < 2) { kindsEl.innerHTML = ''; out.innerHTML = '<div class="hint">Type at least two characters.</div>'; return; }
          let res;
          try {
            res = B.searchLibrary(recs, q, { enzymes: PB.allEnzymes(), primers: state.primers,
                                             table: state.codeTable, perRecord: 300 });
          } catch (e) { out.innerHTML = '<div class="note bad">' + esc(e.message) + '</div>'; return; }
          state.lastSearch = res;
          const label = { dna: 'bases', protein: 'peptide', enzyme: 'enzyme', text: 'names and labels' };
          kindsEl.innerHTML = res.kinds.length
            ? '<span class="hint" style="margin:0;align-self:center">read as ' +
              res.kinds.map(k => label[k] || k).join(', ') + '</span>'
            : '<span class="hint" style="margin:0">that query is too short to read as anything</span>';
          if (!res.total) { out.innerHTML = '<div class="note">Nothing found' +
            ($('#lsScope').value === 'one' ? ' in this construct' : ' anywhere in the library') + '.</div>'; return; }
          /* grouped by construct, in library order */
          const groups = [];
          state.records.forEach(r => {
            const hs = res.hits.filter(h => h.recId === r.id);
            if (hs.length) groups.push({ rec: r, hits: hs });
          });
          let html = '<div class="hint" style="margin-top:0">' + fmt(res.total) + ' match' +
            (res.total === 1 ? '' : 'es') + (groups.length > 1 ? ' in ' + groups.length + ' constructs' : '') + '</div>';
          groups.forEach(g => {
            html += '<h3 style="margin:10px 0 4px;font-size:13px">' + esc(g.rec.name) +
              ' <span class="hint" style="margin:0">' + fmt(g.rec.seq.length) + ' bp · ' + g.hits.length + ' match' +
              (g.hits.length === 1 ? '' : 'es') + '</span></h3>' +
              '<div class="scrollx"><table class="tbl"><tbody>' +
              g.hits.slice(0, 60).map(h => {
                const i = res.hits.indexOf(h);
                const span = h.kind === 'name' ? '—'
                  : (h.wraps || h.end <= h.start
                      ? fmt(h.start + 1) + '–' + fmt(g.rec.seq.length) + ', 1–' + fmt(h.end)
                      : fmt(h.start + 1) + '–' + fmt(h.end));
                return '<tr><td style="width:84px"><span class="badge ' +
                  ({ dna: 'acc', protein: 'acc', enzyme: 'neutral', feature: 'ok', name: 'neutral' }[h.kind] || 'neutral') +
                  '">' + h.kind + '</span></td>' +
                  '<td class="mono" style="width:130px;font-size:11.5px">' + span + '</td>' +
                  '<td style="font-size:12px">' + esc(h.what) + '</td>' +
                  '<td style="width:74px"><button class="btn sm ghost" data-lsgo="' + i + '">Show</button></td></tr>';
              }).join('') +
              (g.hits.length > 60 ? '<tr><td colspan="4" class="hint">and ' + (g.hits.length - 60) + ' more</td></tr>' : '') +
              '</tbody></table></div>';
          });
          if (res.primers.length) {
            html += '<h3 style="margin:10px 0 4px;font-size:13px">Primer library <span class="hint" style="margin:0">' +
              res.primers.length + ' match' + (res.primers.length === 1 ? '' : 'es') + '</span></h3>' +
              '<div class="scrollx"><table class="tbl"><tbody>' +
              res.primers.map(p => '<tr><td style="width:84px"><span class="badge neutral">primer</span></td>' +
                '<td style="font-size:12px">' + esc(p.name) + ' <span class="hint" style="margin:0">' + esc(p.what) + '</span></td>' +
                '<td class="seqcell mono" style="font-size:11.5px;max-width:260px">' + esc(p.seq) + '</td></tr>').join('') +
              '</tbody></table></div>';
          }
          out.innerHTML = html;
          out.querySelectorAll('[data-lsgo]').forEach(b => b.addEventListener('click', () => {
            const h = state.lastSearch.hits[+b.getAttribute('data-lsgo')];
            if (h) PB.goToHit(h);
          }));
        };
        let t = null;
        box.addEventListener('input', () => { clearTimeout(t); t = setTimeout(run, 160); });
        box.addEventListener('keydown', e => { if (e.key === 'Enter') { clearTimeout(t); run(); } });
        w.querySelector('#lsScope').addEventListener('change', run);
        box.focus(); box.select();
        if (q0) run();
      }
    });
  };

  /* open a hit where it sits: switch to its construct, select it, show it */
  PB.goToHit = function (h) {
    const rec = state.records.find(r => r.id === h.recId);
    if (!rec) { PB.toast('That construct is no longer in the library.', 'warn'); return; }
    if (state.activeId !== rec.id) {
      state.activeId = rec.id; state.selectedFeature = null;
      PB.invalidateSites();
    }
    if (h.kind === 'name') { state.selection = null; }
    else {
      state.selection = { start: h.start, end: h.end };
      if (h.featureIndex != null) state.selectedFeature = h.featureIndex;
    }
    PB.closeModal();
    state.view = 'seq'; PB.lsSet('view', 'seq');
    PB.renderAll();
    if (state.selection && PB.scrollSeqTo) PB.scrollSeqTo(h.start);
  };

  PB.toolInsertDelete = toolInsertDelete;

  function toolTopology() {
    const rec = PB.active();
    if (!rec) return;
    if (rec.builtin) { PB.toast('Duplicate the demo construct before editing it.', 'warn'); return; }
    rec.topology = rec.topology === 'circular' ? 'linear' : 'circular';
    PB.addHistory(rec, 'Topology changed', 'Now ' + rec.topology);
    PB.touch(rec); PB.invalidateSites(); PB.renderAll();
    PB.toast('Now ' + rec.topology, 'ok');
  }

  PB.setOriginDialog = function () {
    const rec = PB.active();
    if (!rec) return;
    if (rec.topology !== 'circular') { PB.toast('Only circular sequences can be rotated.', 'warn'); return; }
    const sel = state.selection;
    PB.modal({
      title: 'Set origin',
      body: '<div class="row"><div class="field narrow"><label>New base 1 is current bp</label>' +
        '<input class="inp mono" id="soPos" type="number" min="1" max="' + rec.seq.length + '" value="' + ((sel ? sel.start : 0) + 1) + '"></div></div>' +
        '<div class="hint">Rotates the circle so this base becomes position 1. Features move with it.</div>',
      actions: [{ label: 'Cancel' }, { label: 'Rotate', primary: true, onClick: () => {
        const p = Math.max(1, Math.min(rec.seq.length, parseInt($('#soPos').value, 10) || 1)) - 1;
        const rot = B.setOrigin(rec, p);
        const target = rec.builtin
          ? PB.makeRecord({ name: rec.name + ' rotated', seq: rot.seq, topology: 'circular', features: rot.features })
          : rec;
        if (rec.builtin) {
          PB.addHistory(target, 'Origin set', 'Base ' + (p + 1) + ' of ' + rec.name + ' is now base 1');
          PB.addRecord(target);
        } else {
          rec.seq = rot.seq; rec.features = rot.features;
          PB.addHistory(rec, 'Origin set', 'Base ' + (p + 1) + ' is now base 1');
          PB.touch(rec); state.selection = null; PB.invalidateSites(); PB.renderAll();
        }
      } }]
    });
  };

  PB.flipDialog = function () {
    const rec = PB.active();
    if (!rec) return;
    PB.modal({
      title: 'Reverse complement',
      body: '<p>Replace “' + esc(rec.name) + '” with its reverse complement? Features are moved and their strands flipped.</p>',
      actions: [{ label: 'Cancel' }, { label: 'Flip', primary: true, onClick: () => {
        const fl = B.flipRecord(rec);
        if (rec.builtin) {
          const nr = PB.makeRecord({ name: rec.name + ' (rev-comp)', seq: fl.seq, topology: rec.topology, features: fl.features });
          PB.addHistory(nr, 'Reverse complement', 'From ' + rec.name);
          PB.addRecord(nr);
        } else {
          rec.seq = fl.seq; rec.features = fl.features;
          PB.addHistory(rec, 'Reverse complement', 'Whole construct flipped');
          PB.touch(rec); state.selection = null; PB.invalidateSites(); PB.renderAll();
        }
      } }]
    });
  };

  /* ============================== FEATURES ============================== */
  const FEAT_TYPES = ['CDS', 'gene', 'promoter', 'terminator', 'RBS', 'rep_origin', 'protein_bind',
    'primer_bind', 'regulatory', 'polyA_signal', 'intron', 'exon', 'protein_tag', 'misc_recomb',
    'LTR', 'sig_peptide', 'misc_feature'];

  PB.featureDialog = function (idx) {
    const rec = PB.active();
    if (!rec) return;
    const L = rec.seq.length, circ = rec.topology === 'circular';
    const f = idx != null ? rec.features[idx] : null;
    const sel = state.selection;
    const complex = !!(f && (f.parts || f.mixed || f.between));
    const start = f ? f.start : (sel ? sel.start : 0);
    const end = f ? f.end : (sel ? sel.end : Math.min(L, 100));
    const q = (f && f.quals) || {};
    PB.modal({
      title: f ? 'Edit feature' : 'Add feature',
      body: '<div class="row"><div class="field"><label>Label</label><input class="inp" id="ftLabel" value="' +
          esc(f ? f.label : '') + '" placeholder="e.g. EGFP"></div>' +
        '<div class="field"><label>Type</label><select class="inp" id="ftType">' +
          FEAT_TYPES.concat(f && FEAT_TYPES.indexOf(f.type) < 0 ? [f.type] : [])
            .map(t => '<option' + (f && f.type === t ? ' selected' : '') + '>' + esc(t) + '</option>').join('') +
        '</select></div></div>' +
        (complex ? '' :
        '<div class="row"><div class="field narrow"><label>Start (bp)</label><input class="inp mono" id="ftStart" type="number" min="1" max="' + L + '" value="' + (start + 1) + '"></div>' +
        '<div class="field narrow"><label>End (bp)</label><input class="inp mono" id="ftEnd" type="number" min="1" max="' + L + '" value="' + (end || L) + '"></div>' +
        '<div class="field narrow"><label>Strand</label><select class="inp" id="ftStrand">' +
          '<option value="1"' + (f && f.strand === -1 ? '' : ' selected') + '>Forward (+)</option>' +
          '<option value="-1"' + (f && f.strand === -1 ? ' selected' : '') + '>Reverse (−)</option></select></div></div>' +
        (circ ? '<div class="hint">On a circular construct an end before the start makes a feature that runs across the origin.</div>' : '')) +
        '<div class="row"><div class="field wide"><label>' + (complex ? 'Location' : 'Or a GenBank location (overrides start/end/strand)') + '</label>' +
          '<input class="inp mono" id="ftLoc" spellcheck="false" placeholder="join(101..250,401..600) or complement(30..90)" value="' +
          (complex ? esc(B.formatLocation(f, L) || '') : '') + '"></div></div>' +
        (complex ? '<div class="hint">This feature has ' + (f.between ? 'a between-base location' : (f.mixed ? 'parts on both strands' : f.parts.length + ' joined parts')) +
          '. Edit the location text to change it; GenBank syntax with 1-based positions.</div>' : '') +
        '<div class="row"><div class="field wide"><label>Note</label><input class="inp" id="ftNote" value="' + esc(q.note || (f && f.note) || '') + '"></div></div>' +
        (f && Object.keys(q).filter(k => k !== 'note' && k !== 'label').length
          ? '<div class="hint">Other qualifiers are kept: ' + esc(Object.keys(q).filter(k => k !== 'note' && k !== 'label').slice(0, 8).join(', ')) + '.</div>' : ''),
      actions: [{ label: 'Cancel' }].concat(f ? [{ label: 'Delete', onClick: () => {
          rec.features.splice(idx, 1);
          PB.addHistory(rec, 'Feature removed', f.label);
          PB.touch(rec); state.selectedFeature = null; PB.renderAll();
        } }] : []).concat([{ label: f ? 'Save' : 'Add', primary: true, onClick: () => {
        const locText = ($('#ftLoc').value || '').trim();
        let geom;
        if (locText) {
          geom = B.parseLocation(locText, L, circ);
          if (!geom) { PB.toast('That location is not valid for a ' + fmt(L) + ' bp ' + rec.topology + ' sequence.', 'bad'); return false; }
        } else if (complex) {
          PB.toast('Enter a location.', 'bad'); return false;
        } else {
          const s0 = parseInt($('#ftStart').value, 10), e0 = parseInt($('#ftEnd').value, 10);
          if (!(s0 >= 1 && s0 <= L && e0 >= 1 && e0 <= L)) { PB.toast('Start and end must be between 1 and ' + fmt(L) + '.', 'bad'); return false; }
          if (e0 < s0 && !circ) { PB.toast('The end must be after the start on a linear sequence.', 'bad'); return false; }
          geom = { start: s0 - 1, end: e0 < s0 ? e0 : e0, strand: parseInt($('#ftStrand').value, 10) };
          if (e0 < s0) geom.end = e0;            // wraps the origin
        }
        const quals = Object.assign({}, q);
        const note = ($('#ftNote').value || '').trim();
        if (note) quals.note = note; else delete quals.note;
        const label = ($('#ftLabel').value || 'feature').trim();
        if (quals.label != null) quals.label = label;
        const nf = PB.normFeature(Object.assign({}, f || {}, geom, {
          label, type: $('#ftType').value, quals, color: f ? f.color : null,
          parts: geom.parts, mixed: geom.mixed, between: geom.between
        }));
        if (!B.featureLength(nf, L) && !nf.between) { PB.toast('That location is empty.', 'bad'); return false; }
        const where = B.formatLocation(nf, L);
        if (f) { rec.features[idx] = nf; PB.addHistory(rec, 'Feature edited', nf.label + ' ' + where); }
        else { rec.features.push(nf); PB.addHistory(rec, 'Feature added', nf.label + ' ' + where); }
        PB.touch(rec); PB.renderAll();
      } }])
    });
  };

  PB.detectFeatures = function () {
    const rec = PB.active();
    if (!rec || !rec.seq.length) return;
    const found = B.detectFeatures(rec.seq, null, { circular: rec.topology === 'circular' });
    const existing = new Set(rec.features.map(f => f.start + ':' + f.end + ':' + f.label));
    const added = [];
    for (const f of found) {
      const nf = PB.normFeature({ type: f.type, start: f.start, end: f.end, strand: f.strand, label: f.name });
      const k = nf.start + ':' + nf.end + ':' + nf.label;
      if (existing.has(k)) continue;
      existing.add(k);
      nf.confident = f.confident; nf.expected = f.expected;
      added.push(nf);
    }
    if (!added.length) {
      PB.toast('No new features from the built-in library. Known elements may already be annotated.', 'warn');
      return;
    }
    const sure = added.filter(f => f.confident), maybe = added.filter(f => !f.confident);
    const row = (f, i) => '<tr><td><input type="checkbox" data-df="' + i + '"' + (f.confident ? ' checked' : '') +
        ' style="accent-color:var(--accent)"></td>' +
      '<td>' + esc(f.label) + '</td><td class="mono" style="font-size:11.5px">' + esc(f.type) + '</td>' +
      '<td class="num">' + fmt(f.start + 1) + '</td><td class="num">' + fmt(f.end) + '</td>' +
      '<td>' + (f.strand === -1 ? '−' : '+') + '</td>' +
      '<td>' + (f.confident
        ? '<span class="badge ok">specific</span>'
        : '<span class="badge warn">could be chance</span> <span style="color:var(--fg-faint);font-size:11px">' +
          (f.expected >= 0.01 ? f.expected.toFixed(2) : '&lt;0.01') + ' expected at random</span>') + '</td></tr>';
    PB.modal({
      title: 'Detected features', sub: sure.length + ' specific, ' + maybe.length + ' uncertain', wide: true,
      body: '<div class="scrollx scrolly"><table class="tbl"><thead><tr><th style="width:34px"></th><th>Label</th><th>Type</th>' +
        '<th class="num">Start</th><th class="num">End</th><th>Strand</th><th>How specific is the match</th></tr></thead><tbody>' +
        added.map(row).join('') + '</tbody></table></div>' +
        (maybe.length ? '<div class="note warn"><strong>' + maybe.length + ' match' + (maybe.length === 1 ? ' is' : 'es are') +
          ' too short to mean much.</strong> The number is how many times a signature that short would turn up in a sequence this size by chance alone, counting both strands. A six-base ribosome binding site inside a gene is a coincidence, not a finding — these are left unticked.</div>' : '') +
        '<div class="hint">Matches come from a curated library of common plasmid elements (promoters, tags, markers, origins, recombination sites). It is a starter set, not a complete database — anything it does not know, add by hand. Short functional motifs are also checked for context: a ribosome binding site is only reported when a start codon follows it at a workable distance.</div>',
      onOpen: w => {
        w.querySelectorAll('[data-df]').forEach(cb => cb.addEventListener('change', () => {
          added[+cb.getAttribute('data-df')].pick = cb.checked;
        }));
        added.forEach(f => { f.pick = !!f.confident; });
      },
      actions: [{ label: 'Cancel' }, { label: 'Add ticked', primary: true, onClick: () => {
        const take = added.filter(f => f.pick).map(f => {
          const g = Object.assign({}, f); delete g.confident; delete g.expected; delete g.pick; return g;
        });
        if (!take.length) { PB.toast('Nothing ticked.', 'warn'); return false; }
        rec.features = rec.features.concat(take);
        PB.addHistory(rec, 'Features detected', take.length + ' element(s) annotated from the built-in library' +
          (take.length < added.length ? '; ' + (added.length - take.length) + ' uncertain match(es) left out' : ''));
        PB.touch(rec); PB.renderAll();
      } }]
    });
  };

  PB.customEnzymeDialog = function () {
    PB.modal({
      title: 'Add a custom enzyme',
      body: '<div class="row"><div class="field"><label>Name</label><input class="inp" id="ceName" placeholder="MyEnzI"></div>' +
        '<div class="field"><label>Recognition site (IUPAC)</label><input class="inp mono" id="ceSite" placeholder="GGTCTC"></div></div>' +
        '<div class="row"><div class="field narrow"><label>Type</label><select class="inp" id="ceType">' +
        '<option value="normal">Cuts inside the site</option><option value="iis">Type IIS (cuts outside)</option></select></div>' +
        '<div class="field narrow"><label>Top cut</label><input class="inp mono" id="ceTop" type="number" value="1"></div>' +
        '<div class="field narrow"><label>Bottom cut</label><input class="inp mono" id="ceBot" type="number" value="5"></div></div>' +
        '<div class="hint">Cut positions are counted from the start of the recognition site for a normal enzyme, and from the <em>end</em> of the site for a Type IIS enzyme. For example EcoRI G^AATT_C is 1 and 5; BsaI GGTCTC(1/5) is 1 and 5 as Type IIS.</div>' +
        '<div class="hint" id="cePreview"></div>',
      actions: [{ label: 'Cancel' }, { label: 'Add', primary: true, onClick: () => {
        const name = ($('#ceName').value || '').trim();
        const site = B.cleanSeq($('#ceSite').value);
        if (!name) { PB.toast('Give it a name.', 'bad'); return false; }
        if (PB.enzByName(name)) { PB.toast('An enzyme with that name already exists.', 'bad'); return false; }
        if (site.length < 3) { PB.toast('The recognition site is too short.', 'bad'); return false; }
        const iis = $('#ceType').value === 'iis';
        const top = parseInt($('#ceTop').value, 10) || 0, bot = parseInt($('#ceBot').value, 10) || 0;
        const e = iis
          ? { name, site, o1: top, o2: bot, typeIIS: true, custom: true }
          : { name, site, top, bot, typeIIS: false, custom: true };
        state.customEnzymes.push(e);
        PB.lsSet('customEnzymes', state.customEnzymes);
        if (state.enzymes.indexOf(name) === -1) state.enzymes.push(name);
        PB.lsSet('enzymes', state.enzymes);
        PB.invalidateSites(); PB.renderAll();
        PB.toast('Added ' + name, 'ok');
      } }]
    });
  };

  /* ============================== IMPORT / EXPORT ============================== */
  PB.pickFiles = function (accept) {
    const inp = $('#fileInput');
    if (accept) inp.setAttribute('accept', accept);
    inp.value = '';
    inp.click();
  };

  PB.openImport = function () {
    PB.modal({
      title: 'Import', wide: true,
      body: '<div class="row"><button class="btn primary" id="imFile">Choose files…</button>' +
        '<span class="hint" style="margin:0;align-self:center">SnapGene <b>.dna</b> (features, colours and primers preserved), ' +
        'GenBank (.gb, .gbk, .genbank, .ape), EMBL, FASTA, GFF3, Benchling and openVectorEditor JSON, ' +
        'NCBI INSDSeq/TinySeq XML, plain sequence, .ab1 chromatograms, or a PlasmidBench .json backup.</span></div>' +
        '<div class="row"><div class="field wide"><label>…or paste here</label>' +
        '<textarea class="inp" id="imText" rows="11" spellcheck="false" placeholder="Paste GenBank, FASTA, EMBL, GFF3, JSON, XML, or just the bases.&#10;&#10;From Addgene or NCBI: open the full GenBank text, select all, and paste it here."></textarea></div></div>' +
        '<div class="note">The format is worked out from the content, not the file name, so a GenBank file saved as .txt or .ape still reads as GenBank.<br><br>' +
        'Anything imported without its own annotations is scanned against the built-in element library automatically, so promoters, tags, markers and origins get labelled on arrival. ' +
        'A GFF3 file with no bases in it is applied as annotations to the construct already open, and only if every one of them fits inside it.<br><br>' +
        'Formats that are a private binary layout — CLC .clc, Sequencher .sbd, Vector NTI .ma4, Serial Cloner .xdna — cannot be read by anything but the program that wrote them; ' +
        'exporting GenBank from that program gives a file that carries the sequence and every annotation.<br><br>' +
        'Downloading directly from NCBI or Addgene is not possible from inside a published page — the sandbox blocks outside network requests. Copy the GenBank text from their site and paste it here; everything, features included, is preserved.</div>',
      actions: [{ label: 'Close' }, { label: 'Import pasted text', primary: true, onClick: () => {
        const n = PB.importText($('#imText').value, 'pasted');
        if (!n) { PB.toast('Could not find a sequence in that text.', 'bad'); return false; }
      } }],
      onOpen: w => {
        w.querySelector('#imFile').addEventListener('click', () => PB.pickFiles());
      }
    });
  };

  // A binary file read as text yields stray ACGT letters out of its byte soup,
  // which would otherwise import as a plausible-looking but meaningless
  // sequence. Refuse it instead.
  function looksBinary(text) {
    const sample = text.slice(0, 4000);
    if (!sample.length) return false;
    let ctrl = 0;
    for (let i = 0; i < sample.length; i++) {
      const c = sample.charCodeAt(i);
      if (c === 9 || c === 10 || c === 13) continue;
      if (c < 32 || c === 0xFFFD) ctrl++;
    }
    return ctrl / sample.length > 0.02;
  }
  PB.looksBinary = looksBinary;

  PB.importText = function (text, source) {
    if (!text || !text.trim()) return 0;
    if (looksBinary(text)) {
      PB.toast('That looks like a binary file, not sequence text. PlasmidBench reads SnapGene .dna and .ab1 directly, ' +
        'and GenBank, EMBL, FASTA, GFF3, JSON and XML as text.', 'bad');
      return 0;
    }
    /* One reader picks the format from the content, not the extension,
       because the extension lies: ApE writes GenBank into .ape, people
       rename .gbk to .txt to email it, and .json is four schemas. */
    let recs = [];
    try { recs = B.parseAnySequence(text, source); }
    catch (e) { PB.toast(e.message || 'That text could not be read.', 'bad'); return 0; }
    const format = recs.format;

    /* GFF with no bases in it: annotations looking for a construct */
    if (recs.annotationsOnly) return applyGFFAnnotations(recs, source);

    // protein records would be misread as DNA (M, K, R, S, ... are IUPAC
    // nucleotide codes too), so they are refused rather than imported
    const proteins = recs.filter(r => r.molecule === 'protein');
    /* a GenBank CONTIG record points at pieces held in other files and has
       no bases of its own; noticed before the empty records are dropped */
    const contigOnly = recs.filter(r => r.contigOnly).length;
    recs = recs.filter(r => r.molecule !== 'protein' && r.seq && r.seq.length);
    const notes = [];
    if (proteins.length) {
      notes.push((proteins.length === 1 ? '“' + proteins[0].name + '” looks like a protein sequence' : proteins.length + ' records look like protein') +
        ' and ' + (proteins.length === 1 ? 'was' : 'were') + ' not imported (Tools → Back-translate makes a coding sequence from a protein)');
    }
    if (contigOnly) notes.push(contigOnly + ' CONTIG record' + (contigOnly === 1 ? '' : 's') +
      ' name pieces held in other files and carry no bases — download the full GenBank record instead');
    if (!recs.length) {
      PB.toast(notes.length ? notes.join('. ') + '.'
        : 'Nothing to import: that ' + (format || 'text') + ' holds no bases.', 'bad');
      return proteins.length || contigOnly ? -1 : 0;
    }
    const dropped = recs.reduce((a, r) => a + (r.dropped || 0), 0);
    if (dropped) notes.push('ignored ' + dropped + ' character' + (dropped === 1 ? '' : 's') + ' that are not nucleotide codes (' +
      recs.map(r => r.droppedChars || '').join('').slice(0, 12) + ') — check the sequence is complete');
    if (recs.some(r => r.hadU)) notes.push('RNA bases (U) were converted to T');
    const skipped = recs.reduce((a, r) => a.concat(r.skippedFeatures || []), []);
    if (skipped.length) notes.push(skipped.length + ' feature' + (skipped.length === 1 ? '' : 's') +
      ' with a location outside this sequence (e.g. ' + skipped[0] + ') could not be placed and ' + (skipped.length === 1 ? 'was' : 'were') + ' left out');
    let last = null, autoTotal = 0;
    recs.forEach(r => {
      const nr = PB.makeRecord(r);
      PB.addHistory(nr, 'Imported', (source || 'text') + (format ? ' · ' + format : '') +
        ' · ' + fmt(nr.seq.length) + ' bp, ' + nr.features.length + ' features' +
        (r.dropped ? '; ' + r.dropped + ' non-nucleotide characters ignored' : '') +
        (r.skippedFeatures && r.skippedFeatures.length ? '; not placed: ' + r.skippedFeatures.join(', ') : ''));
      autoTotal += PB.autoAnnotate(nr, true);
      state.records.push(nr);
      PB.persist(nr);
      last = nr;
    });
    if (last) { state.activeId = last.id; state.selection = null; state.selectedFeature = null; }
    PB.invalidateSites(); PB.closeModal(); PB.renderAll();
    PB.toast('Imported ' + recs.length + ' sequence' + (recs.length === 1 ? '' : 's') +
      (format && format !== 'raw sequence' ? ' from ' + format : '') +
      (autoTotal ? ' · ' + autoTotal + ' features auto-detected' : '') +
      (notes.length ? '. Note: ' + notes.join('; ') + '.' : ''), notes.length ? 'warn' : 'ok');
    return recs.length;
  };

  /* A GFF3 file is often annotations on their own — Prokka, bakta and
     Geneious all write one — and they are only meaningful against the
     construct they were called on. Rather than invent a sequence of Ns to
     hang them on, they are applied to the construct already open, and only
     if every one of them fits inside it. A feature silently moved to the
     nearest legal position would be worse than no feature at all: you would
     design primers against it. */
  function applyGFFAnnotations(recs, source) {
    const rec = PB.active();
    const feats = recs.reduce((a, r) => a.concat(r.features || []), []);
    if (!feats.length) { PB.toast('That GFF file has no features in it.', 'bad'); return 0; }
    if (!rec || !rec.seq.length) {
      PB.toast('That GFF file holds ' + feats.length + ' annotation' + (feats.length === 1 ? '' : 's') +
        ' but no sequence. Open the construct they belong to first, then import it again — ' +
        'or import a GFF with its bases appended after a ##FASTA line.', 'bad');
      return 0;
    }
    const L = rec.seq.length;
    const outside = feats.filter(f => f.end > L || f.start >= L);
    if (outside.length) {
      PB.toast(feats.length + ' annotation' + (feats.length === 1 ? '' : 's') + ' in that GFF run past the end of “' +
        rec.name + '” (' + fmt(L) + ' bp) — for example ' + outside[0].label + ' at ' + fmt(outside[0].end) +
        '. They belong to a different sequence, so nothing was changed.', 'bad');
      return 0;
    }
    rec.features = rec.features.concat(feats.map(PB.normFeature));
    PB.addHistory(rec, 'Annotations from GFF3', (source || 'GFF') + ' · ' + feats.length +
      ' feature' + (feats.length === 1 ? '' : 's') + ' added');
    PB.touch(rec);
    PB.invalidateSites(); PB.closeModal(); PB.renderAll();
    PB.toast('Added ' + feats.length + ' annotation' + (feats.length === 1 ? '' : 's') + ' from GFF3 to “' + rec.name + '”', 'ok');
    return feats.length;
  }

  /* Annotate an imported record from the built-in element library when it
     arrived without any features of its own (plain FASTA, bare sequence, or a
     GenBank file stripped of its FEATURES block). */
  PB.autoAnnotate = function (rec, quiet) {
    if (!rec || !rec.seq.length || (rec.features && rec.features.length)) return 0;
    let found = [];
    try { found = B.detectFeatures(rec.seq, null, { circular: rec.topology === 'circular' }); } catch (e) { return 0; }
    if (!found.length) return 0;
    const seen = new Set();
    const added = [];
    let skipped = 0;
    for (const f of found) {
      // on import, only annotate matches long enough to be specific; a short
      // motif that would turn up by chance is left for Detect features
      if (!f.confident) { skipped++; continue; }
      const nf = PB.normFeature({ type: f.type, start: f.start, end: f.end, strand: f.strand, label: f.name });
      const k = nf.start + ':' + nf.end + ':' + nf.label;
      if (seen.has(k)) continue;
      seen.add(k); added.push(nf);
    }
    if (!added.length) return 0;
    rec.features = rec.features.concat(added);
    PB.addHistory(rec, 'Features detected on import', added.length + ' elements annotated from the built-in library' +
      (skipped ? '; ' + skipped + ' match(es) too short to be specific were left out — run Detect features to see them' : ''));
    if (!quiet) PB.toast('Annotated ' + added.length + ' element' + (added.length === 1 ? '' : 's') + ' automatically', 'ok');
    return added.length;
  };

  function importSnapGene(buffer, fileName) {
    const p = B.parseSnapGene(buffer);
    const nr = PB.makeRecord({
      name: p.name || fileName.replace(/\.dna$/i, ''),
      seq: p.seq, topology: p.topology,
      features: p.features.map(PB.normFeature)
    });
    PB.addHistory(nr, 'Imported from SnapGene', fileName + ' · ' + fmt(nr.seq.length) + ' bp ' +
      nr.topology + ', ' + nr.features.length + ' features' +
      (p.primers.length ? ', ' + p.primers.length + ' primers' : ''));
    // anything the file did not annotate, fill in from the built-in library
    const extra = PB.autoAnnotate(nr, true);
    let addedPrimers = 0;
    (p.primers || []).forEach(pr => {
      if (!pr.seq || pr.seq.length < 8) return;
      if (state.primers.some(x => x.seq === pr.seq)) return;
      state.primers.push({ name: pr.name, seq: pr.seq }); addedPrimers++;
    });
    if (addedPrimers) PB.lsSet('primers', state.primers);
    state.records.push(nr);
    PB.persist(nr);
    state.activeId = nr.id; state.selection = null; state.selectedFeature = null;
    PB.invalidateSites(); PB.closeModal(); PB.renderAll();
    PB.toast('Imported ' + nr.name + ' — ' + nr.features.length + ' feature' +
      (nr.features.length === 1 ? '' : 's') +
      (extra ? ' (' + extra + ' auto-detected)' : '') +
      (addedPrimers ? ', ' + addedPrimers + ' primers' : ''), 'ok');
    return nr;
  }

  PB.handleFiles = function (files) {
    if (!files || !files.length) return;
    /* sequencing reads — traces and the .seq text files that come with
       them — go to the verify view together, in one batch */
    const all = Array.prototype.slice.call(files);
    const readFiles = all.filter(f => /\.(ab1|abi|seq)$/i.test(f.name));
    if (readFiles.length && PB.loadReadFiles) {
      PB.closeModal();
      PB.loadReadFiles(readFiles);
    }
    all.filter(f => readFiles.indexOf(f) < 0).forEach(file => {
      const nm = file.name.toLowerCase();
      const ext = (/\.([a-z0-9_]+)$/.exec(nm) || [, ''])[1];
      if (/\.prot$/.test(nm)) {
        PB.toast(file.name + ' is a SnapGene protein file — PlasmidBench works on DNA. Use Tools → Back-translate to make a coding sequence.', 'bad');
        return;
      }
      /* Formats that are somebody else's private binary layout, or a
         picture. Saying plainly that the file cannot be read, and naming
         the export that would work, beats reading half of it wrongly. */
      if (B.UNREADABLE[ext]) {
        PB.toast(file.name + ': ' + B.UNREADABLE[ext], 'bad');
        return;
      }
      if (/\.dna$/.test(nm)) {
        const fr = new FileReader();
        fr.onload = () => {
          try { importSnapGene(fr.result, file.name); }
          catch (e) { PB.toast(file.name + ': ' + (e.message || 'could not be read'), 'bad'); }
        };
        fr.onerror = () => PB.toast('Could not read ' + file.name, 'bad');
        fr.readAsArrayBuffer(file);
        return;
      }
      if (/\.(ab1|abi)$/.test(nm)) {
        const fr = new FileReader();
        fr.onload = () => {
          try {
            const p = B.parseAB1(fr.result);
            if (!p.sequence.length) throw new Error('No base calls found in that file.');
            PB.addRead({
              name: file.name.replace(/\.(ab1|abi)$/i, ''),
              seq: p.sequence, quality: p.quality, traces: p.traces,
              peakLocations: p.peakLocations, hasTraces: p.hasTraces,
              meanQuality: p.meanQuality, source: file.name
            });
            PB.closeModal();
          } catch (e) {
            PB.toast(file.name + ': ' + (e.message || 'could not be read'), 'bad');
          }
        };
        fr.onerror = () => PB.toast('Could not read ' + file.name, 'bad');
        fr.readAsArrayBuffer(file);
        return;
      }
      const fr = new FileReader();
      fr.onload = () => {
        const txt = String(fr.result || '');
        // a SnapGene file renamed or saved without its extension still starts
        // with packet type 9 and the word SnapGene
        if (txt.charCodeAt(0) === 9 && txt.slice(0, 40).indexOf('SnapGene') !== -1) {
          const br = new FileReader();
          br.onload = () => {
            try { importSnapGene(br.result, file.name); }
            catch (e) { PB.toast(file.name + ': ' + (e.message || 'could not be read'), 'bad'); }
          };
          br.readAsArrayBuffer(file);
          return;
        }
        if (/\.json$/.test(nm)) {
          /* PlasmidBench's own backup first: it carries ids, history and
             half-open coordinates, so it is restored rather than re-read.
             Any other JSON is a sequence file from somewhere else and goes
             to the JSON reader, which works out which schema it is. */
          let data = null;
          try { data = JSON.parse(txt); }
          catch (e) { PB.toast(file.name + ' is not valid JSON: ' + (e.message || 'could not be parsed'), 'bad'); return; }
          const list = Array.isArray(data) ? data : ((data && data.records) || []);
          const isBackup = list.length && list.every(r => r && typeof r.seq === 'string' &&
            r.seq.length && Array.isArray(r.history));
          if (isBackup) {
            let n = 0;
            list.forEach(r => {
              const nr = PB.makeRecord(r);
              nr.id = PB.uid();
              PB.addHistory(nr, 'Restored from backup', file.name);
              state.records.push(nr); PB.persist(nr); n++;
            });
            if (state.records.length) state.activeId = state.records[state.records.length - 1].id;
            PB.closeModal(); PB.renderAll();
            PB.toast('Restored ' + n + ' construct' + (n === 1 ? '' : 's'), 'ok');
            return;
          }
          if (!PB.importText(txt, file.name)) PB.toast(file.name + ': no sequence found in that JSON.', 'bad');
          return;
        }
        if (!PB.importText(txt, file.name)) {
          PB.toast(file.name + ': no sequence recognised.', 'bad');
        }
      };
      fr.onerror = () => PB.toast('Could not read ' + file.name, 'bad');
      fr.readAsText(file);
    });
  };

  PB.openExport = function () {
    const rec = PB.active();
    if (!rec) { PB.toast('Nothing to export.', 'warn'); return; }
    const gb = B.toGenBank(rec);
    const fa = B.toFASTA(rec.name.replace(/\s+/g, '_'), rec.seq);
    const safe = rec.name.replace(/[^\w.-]+/g, '_').slice(0, 40) || 'construct';
    PB.modal({
      title: 'Export', sub: rec.name, wide: true,
      body: '<div class="row">' +
        '<div class="field" style="flex:1 1 220px"><label>This construct</label>' +
        '<div class="row" style="gap:5px;margin:0">' +
        '<button class="btn" id="exGb">GenBank</button>' +
        '<button class="btn" id="exFa">FASTA</button>' +
        '<button class="btn" id="exSeq">Bases only</button>' +
        '<button class="btn" id="exImg">Map image (SVG / PNG)</button>' +
        '</div></div>' +
        '<div class="field" style="flex:1 1 220px"><label>Whole library</label>' +
        '<div class="row" style="gap:5px;margin:0">' +
        '<button class="btn" id="exJson">JSON backup</button>' +
        '<button class="btn" id="exAllFa">All as FASTA</button>' +
        '</div></div></div>' +
        '<div class="hint" id="exNote"></div>' +
        '<div class="row"><div class="field wide"><label>GenBank preview</label>' +
        '<textarea class="inp" rows="12" readonly spellcheck="false" id="exPrev"></textarea></div></div>',
      onOpen: w => {
        w.querySelector('#exPrev').value = gb;
        const note = w.querySelector('#exNote');
        note.textContent = state.downloads
          ? 'Files are saved through your browser after you confirm. GenBank is written as .txt because that is what this view can save — rename it to .gb if your other software wants that extension.'
          : 'Saving files is not available in this view, so these buttons copy the text to your clipboard instead.';
        w.querySelector('#exGb').addEventListener('click', () => PB.offerFile(safe + '_genbank.txt', gb, 'GenBank'));
        w.querySelector('#exImg').addEventListener('click', () => { PB.closeModal(); PB.exportMapDialog(); });
        w.querySelector('#exFa').addEventListener('click', () => PB.offerFile(safe + '.txt', fa, 'FASTA'));
        w.querySelector('#exSeq').addEventListener('click', () => PB.copyText(rec.seq, 'Sequence'));
        w.querySelector('#exJson').addEventListener('click', () => {
          const data = JSON.stringify({
            app: 'PlasmidBench', exported: new Date().toISOString(),
            records: state.records.filter(r => !r.builtin).map(r => ({
              name: r.name, seq: r.seq, topology: r.topology, features: r.features,
              description: r.description, history: r.history, createdAt: r.createdAt, updatedAt: r.updatedAt
            })),
            primers: state.primers, customEnzymes: state.customEnzymes
          }, null, 1);
          PB.offerFile('plasmidbench_backup.json', data, 'Library backup');
        });
        w.querySelector('#exAllFa').addEventListener('click', () => {
          const all = state.records.map(r => B.toFASTA(r.name.replace(/\s+/g, '_'), r.seq)).join('');
          PB.offerFile('plasmidbench_library.txt', all, 'All sequences');
        });
      }
    });
  };

  /* ======================= REACTION SETUP =======================
     The arithmetic between a finished design and a rack of tubes: how
     much of each stock goes in, what the block runs, and what the molar
     ratio really is once the fragments' own base composition is counted.
     Whatever was last designed in the assembly dialog is filled in, so
     the sizes come from the construct rather than from memory.
  ------------------------------------------------------------------ */
  function toolBench(startMode) {
    const d = state.gbDesign || null;
    const rec = PB.active();
    /* two decimals at most, and none when the number is round: a mix that
       says 5 µL is read faster than one saying 5.00 µL */
    const n2 = x => (x == null ? '—' : (Math.round(x * 100) / 100).toLocaleString('en-US',
      { minimumFractionDigits: 0, maximumFractionDigits: 2 }));
    const n1 = x => (x == null ? '—' : (Math.round(x * 10) / 10).toString());
    const num = (id, v, step, min) => '<input class="inp mono" id="' + id + '" type="number" value="' + v +
      '" step="' + (step || 1) + '"' + (min != null ? ' min="' + min + '"' : '') + '>';

    const w = PB.modal({
      title: 'Reaction setup',
      sub: 'master mixes, molar ratios and cycling',
      wide: true,
      body:
        '<div class="row" style="gap:6px;margin-bottom:10px">' +
          '<label class="chip on" data-bm="pcr">PCR</label>' +
          '<label class="chip" data-bm="asm">Assembly</label>' +
          '<label class="chip" data-bm="conv">Conversions</label>' +
        '</div>' +
        /* ---- PCR ---- */
        '<div id="bnPcr">' +
          '<div class="row">' +
            '<div class="field" style="flex:1 1 260px"><label>Polymerase</label><select class="inp" id="bnPol">' +
              Object.keys(B.POLYMERASES).map(k => '<option value="' + k + '">' + esc(B.POLYMERASES[k].name) + '</option>').join('') +
            '</select></div>' +
            '<div class="field narrow"><label>Reaction (µL)</label>' + num('bnVol', 50, 5, 5) + '</div>' +
            '<div class="field narrow"><label>Reactions</label>' + num('bnN', 1, 1, 1) + '</div>' +
            '<div class="field narrow"><label>Overage (%)</label>' + num('bnOver', 10, 5, 0) + '</div>' +
          '</div>' +
          '<div class="row">' +
            '<div class="field narrow"><label>Primer stock (µM)</label>' + num('bnPstock', 10, 1, 0.1) + '</div>' +
            '<div class="field narrow"><label>Primer final (µM)</label>' + num('bnPfinal', 0.5, 0.1, 0.05) + '</div>' +
            '<div class="field narrow"><label>Template (ng)</label>' + num('bnTng', 1, 0.5, 0) + '</div>' +
            '<div class="field narrow"><label>Template stock (ng/µL)</label>' + num('bnTstock', 10, 1, 0) + '</div>' +
            '<div class="field narrow" style="align-self:flex-end"><label class="chip" id="bnOpt">add enhancer</label></div>' +
          '</div>' +
          '<div class="row">' +
            '<div class="field narrow"><label>Forward Tm (°C)</label>' + num('bnTmF', d && d.rows && d.rows[0] ? n1(d.rows[0][3]) : 60, 0.1) + '</div>' +
            '<div class="field narrow"><label>Reverse Tm (°C)</label>' + num('bnTmR', d && d.rows && d.rows[1] ? n1(d.rows[1][3]) : 60, 0.1) + '</div>' +
            '<div class="field narrow"><label>Forward nt</label>' + num('bnLenF', d && d.rows && d.rows[0] ? d.rows[0][2].length : 20, 1, 6) + '</div>' +
            '<div class="field narrow"><label>Reverse nt</label>' + num('bnLenR', d && d.rows && d.rows[1] ? d.rows[1][2].length : 20, 1, 6) + '</div>' +
            '<div class="field narrow"><label>Product (bp)</label>' + num('bnProd', d && d.vecSeq ? d.vecSeq.length : (rec ? rec.seq.length : 1000), 100, 1) + '</div>' +
          '</div>' +
          (d && d.rows && d.rows.length ? '<div class="row" style="gap:6px"><div class="field" style="flex:0 0 auto">' +
            '<label>Take the primers from the last design</label>' +
            d.rows.map((r, i) => '<button class="btn sm" data-bnpair="' + i + '">' + esc(r[0]) + '</button>').join(' ') +
            '</div></div>' : '') +
          '<div id="bnPcrOut"></div>' +
        '</div>' +
        /* ---- assembly ---- */
        '<div id="bnAsm" hidden>' +
          '<div class="row">' +
            '<div class="field" style="flex:1 1 260px"><label>Kit</label><select class="inp" id="bnKit">' +
              Object.keys(B.ASSEMBLY_KITS).map(k => '<option value="' + k + '">' + esc(B.ASSEMBLY_KITS[k].name) + '</option>').join('') +
            '</select></div>' +
            '<div class="field narrow"><label>Reaction (µL)</label>' + num('bnAvol', 20, 5, 5) + '</div>' +
            '<div class="field narrow"><label>Vector (ng)</label>' + num('bnVecNg', 100, 10, 1) + '</div>' +
            '<div class="field narrow"><label>Insert : vector</label>' + num('bnRatio', 2, 0.5, 0.1) + '</div>' +
            '<div class="field narrow"><label>Reactions</label>' + num('bnAn', 1, 1, 1) + '</div>' +
          '</div>' +
          '<div class="row"><div class="field wide"><label>Fragments — one per line: name, length in bp, stock ng/µL</label>' +
            '<textarea class="inp mono" id="bnFrags" rows="4" spellcheck="false"></textarea>' +
            '<div class="hint" id="bnFragHint" style="margin:4px 0 0"></div></div></div>' +
          '<div id="bnAsmOut"></div>' +
        '</div>' +
        /* ---- conversions ---- */
        '<div id="bnConv" hidden>' +
          '<div class="row">' +
            '<div class="field narrow"><label>Length (bp)</label>' + num('bnCbp', d && d.vecSeq ? d.vecSeq.length : 1000, 100, 1) + '</div>' +
            '<div class="field narrow"><label>Amount (ng)</label>' + num('bnCng', 100, 10, 0) + '</div>' +
            '<div class="field narrow"><label>or (pmol)</label>' + num('bnCpmol', 0, 0.01, 0) + '</div>' +
          '</div>' +
          '<div id="bnConvOut"></div>' +
          '<div class="row" style="margin-top:12px">' +
            '<div class="field narrow"><label>Stock (ng/µL or µM)</label>' + num('bnDstock', 100, 10, 0) + '</div>' +
            '<div class="field narrow"><label>Wanted</label>' + num('bnDwant', 10, 1, 0) + '</div>' +
            '<div class="field narrow"><label>Final volume (µL)</label>' + num('bnDvol', 50, 5, 1) + '</div>' +
          '</div>' +
          '<div id="bnDilOut"></div>' +
          '<div class="row" style="margin-top:12px">' +
            '<div class="field narrow"><label>Oligo delivered (nmol)</label>' + num('bnOnmol', 25, 1, 0) + '</div>' +
            '<div class="field narrow"><label>Stock wanted (µM)</label>' + num('bnOuM', 100, 10, 1) + '</div>' +
          '</div>' +
          '<div id="bnOligoOut"></div>' +
        '</div>'
    });

    const $$ = id => w.querySelector(id);
    const val = (id, dflt) => { const e = $$(id); const v = e ? parseFloat(e.value) : NaN; return isFinite(v) ? v : dflt; };
    let useOpt = false;

    w.querySelectorAll('[data-bm]').forEach(c => c.addEventListener('click', () => {
      w.querySelectorAll('[data-bm]').forEach(x => x.classList.toggle('on', x === c));
      const m = c.getAttribute('data-bm');
      $$('#bnPcr').hidden = m !== 'pcr'; $$('#bnAsm').hidden = m !== 'asm'; $$('#bnConv').hidden = m !== 'conv';
      run();
    }));
    const oc = $$('#bnOpt');
    oc.addEventListener('click', () => { useOpt = !useOpt; oc.classList.toggle('on', useOpt); run(); });
    w.querySelectorAll('[data-bnpair]').forEach(b => b.addEventListener('click', () => {
      const i = parseInt(b.getAttribute('data-bnpair'), 10);
      const rows = (state.gbDesign || {}).rows || [];
      const pair = i % 2 === 0 ? [rows[i], rows[i + 1]] : [rows[i - 1], rows[i]];
      if (!pair[0] || !pair[1]) return;
      $$('#bnTmF').value = n1(pair[0][3]); $$('#bnTmR').value = n1(pair[1][3]);
      $$('#bnLenF').value = pair[0][2].length; $$('#bnLenR').value = pair[1][2].length;
      run();
    }));

    /* fragments start from the design: the backbone and what goes into it */
    if (d && d.vecSeq) {
      const lines = ['backbone, ' + d.vecSeq.length + ', 50'];
      if (d.insSeq && d.insSeq.length) lines.push((d.poi || 'insert') + ', ' + d.insSeq.length + ', 25');
      $$('#bnFrags').value = lines.join('\n');
      $$('#bnFragHint').textContent = 'Filled in from the last design. Lengths are exact; change the ng/µL to your own measurements.';
    } else {
      $$('#bnFrags').value = 'backbone, 4000, 50\ninsert, 900, 25';
      $$('#bnFragHint').textContent = 'Design an assembly first and these are filled in for you.';
    }

    function parseFrags() {
      return $$('#bnFrags').value.split(/\n+/).map((ln, i) => {
        const p = ln.split(/\s*,\s*/);
        if (p.length < 2 || !p[0].trim()) return null;
        const bp = parseInt(p[1], 10), st = parseFloat(p[2]);
        if (!(bp > 0)) return null;
        /* the real sequence when it is the one just designed, so the molar
           mass comes from its own bases rather than an average */
        let seq = null;
        const dd = state.gbDesign;
        if (dd) {
          if (i === 0 && dd.vecSeq && dd.vecSeq.length === bp) seq = dd.vecSeq;
          if (i === 1 && dd.insSeq && dd.insSeq.length === bp) seq = dd.insSeq;
        }
        return { name: p[0].trim(), bp: bp, seq: seq, stockNgUl: isFinite(st) ? st : 0,
                 role: i === 0 ? 'vector' : 'insert' };
      }).filter(Boolean);
    }

    function pcrOut() {
      const mix = B.pcrMix({
        polymerase: $$('#bnPol').value, volume: val('#bnVol', 50),
        reactions: val('#bnN', 1), overage: val('#bnOver', 10) / 100,
        primerStockUM: val('#bnPstock', 10), primerFinalUM: val('#bnPfinal', 0.5),
        templateNg: val('#bnTng', 1), templateStockNgUl: val('#bnTstock', 10),
        extras: useOpt ? (B.POLYMERASES[$$('#bnPol').value].optional || []).map(x => x.name) : []
      });
      const prog = B.pcrProgram({
        polymerase: $$('#bnPol').value, tmF: val('#bnTmF', 60), tmR: val('#bnTmR', 60),
        lenF: val('#bnLenF', 20), lenR: val('#bnLenR', 20), productBp: val('#bnProd', 1000)
      });
      let h = '';
      if (!mix.ok) h += '<div class="note bad">' + esc(mix.error) + '</div>';
      h += '<div class="scrollx" style="margin-top:8px"><table class="tbl"><thead><tr><th>Component</th>' +
        '<th class="num">1 reaction (µL)</th><th class="num">×' + n2(mix.factor) + ' (µL)</th><th>Final</th></tr></thead><tbody>' +
        mix.rows.map(r => '<tr' + (r.mix ? '' : ' style="opacity:.75"') + '><td>' + esc(r.name) +
          (r.mix ? '' : ' <span class="badge neutral">per tube</span>') + '</td>' +
          '<td class="num mono">' + (r.perRxn == null ? '—' : n2(r.perRxn)) + '</td>' +
          '<td class="num mono">' + (r.mix && r.total != null ? n2(r.total) : '—') + '</td>' +
          '<td>' + esc(r.final || '') + '</td></tr>').join('') +
        '</tbody></table></div>';
      h += '<div class="note ok"><strong>Master mix:</strong> ' + n2(mix.masterMixTotal) + ' µL in all — ' +
        n2(mix.masterMixPerRxn) + ' µL per reaction, so put ' + n2(mix.masterMixPerRxn) +
        ' µL in each tube and add the primers and template to make ' + n2(mix.volume) + ' µL. ' +
        'The overage covers what the tips keep.</div>';
      h += '<div class="scrollx"><table class="tbl"><thead><tr><th>Step</th><th class="num">°C</th><th class="num">Time</th><th></th></tr></thead><tbody>' +
        '<tr><td>Initial denaturation</td><td class="num mono">' + prog.denature + '</td><td class="num mono">' + prog.initialDenatureS + ' s</td><td></td></tr>' +
        '<tr><td>Denature</td><td class="num mono">' + prog.denature + '</td><td class="num mono">' + prog.denatureS + ' s</td>' +
          '<td rowspan="3" style="vertical-align:middle">' + prog.cycles[0] + '–' + prog.cycles[1] + ' cycles</td></tr>' +
        '<tr><td>Anneal</td><td class="num mono"><strong>' + n1(prog.anneal) + '</strong></td><td class="num mono">20 s</td></tr>' +
        '<tr><td>Extend</td><td class="num mono">' + prog.extendC + '</td><td class="num mono">' +
          prog.extendS[0] + (prog.extendS[1] !== prog.extendS[0] ? '–' + prog.extendS[1] : '') + ' s</td></tr>' +
        '<tr><td>Final extension</td><td class="num mono">' + prog.extendC + '</td><td class="num mono">' + prog.finalExtendMin + ' min</td><td></td></tr>' +
        '</tbody></table></div>' +
        '<div class="hint">Annealing is ' + esc(prog.annealNote) + '. Extension is ' +
        prog.polymerase.extendSPerKb[0] + '–' + prog.polymerase.extendSPerKb[1] + ' s/kb over ' +
        n2(val('#bnProd', 1000) / 1000) + ' kb; take the longer end for a GC-rich or difficult template.</div>';
      if (prog.spread > 5) h += '<div class="note warn"><strong>The two primers are ' + n1(prog.spread) +
        ' °C apart.</strong> One anneals well before the other, which is where non-specific product comes from. ' +
        'Re-design to a common Tm, or run a gradient.</div>';
      return h;
    }

    function asmOut() {
      const frags = parseFrags();
      if (!frags.length) return '<div class="note">Add at least one fragment.</div>';
      const s = B.assemblySetup({
        kit: $$('#bnKit').value, fragments: frags, vectorNg: val('#bnVecNg', 100),
        ratio: val('#bnRatio', 2), totalVol: val('#bnAvol', 20)
      });
      if (!s.ok) return '<div class="note bad">' + esc(s.error) + '</div>';
      const n = Math.max(1, val('#bnAn', 1));
      let h = '<div class="scrollx" style="margin-top:8px"><table class="tbl"><thead><tr><th>Fragment</th>' +
        '<th class="num">bp</th><th class="num">pmol</th><th class="num">ng</th><th class="num">stock ng/µL</th>' +
        '<th class="num">µL</th></tr></thead><tbody>' +
        s.rows.map(r => '<tr><td>' + esc(r.name) + (r.role === 'vector' ? ' <span class="badge neutral">vector</span>' : '') + '</td>' +
          '<td class="num mono">' + fmt(r.bp) + '</td><td class="num mono">' + n2(r.pmol) + '</td>' +
          '<td class="num mono">' + n2(r.ng) + '</td><td class="num mono">' + (r.stockNgUl || '—') + '</td>' +
          '<td class="num mono"><strong>' + (r.vol == null ? '—' : n2(r.vol)) + '</strong></td></tr>').join('') +
        '<tr><td>' + esc(s.kit.mixName) + '</td><td class="num">—</td><td class="num">—</td><td class="num">—</td>' +
          '<td class="num">—</td><td class="num mono"><strong>' + n2(s.mixVol) + '</strong></td></tr>' +
        (s.enzVol ? '<tr><td>' + esc(s.kit.enzyme.name) + '</td><td class="num">—</td><td class="num">—</td>' +
          '<td class="num">—</td><td class="num">—</td><td class="num mono"><strong>' + n2(s.enzVol) + '</strong></td></tr>' : '') +
        '<tr><td>Nuclease-free water</td><td class="num">—</td><td class="num">—</td><td class="num">—</td>' +
          '<td class="num">—</td><td class="num mono">' + n2(s.waterVol) + '</td></tr>' +
        '</tbody></table></div>';
      h += '<div class="note ok"><strong>' + n2(s.totalPmol) + ' pmol of DNA in ' + n2(s.totalVol) + ' µL</strong> — ' +
        'insert at ' + n2(s.ratio) + ':1 to the vector by moles, counted from each fragment’s own base composition. ' +
        'Incubate ' + s.kit.tempC + ' °C for ' + (s.kit.minutes[0] === s.kit.minutes[1] ? s.kit.minutes[0] + ' min' :
        s.kit.minutes[0] + '–' + s.kit.minutes[1] + ' min') + ', then transform 2 µL into competent cells.</div>';
      if (n > 1) h += '<div class="note"><strong>For ' + n + ' reactions</strong> make ' +
        n2(s.mixVol * n * 1.1) + ' µL of ' + esc(s.kit.mixName) + ' and ' + n2(s.waterVol * n * 1.1) +
        ' µL of water as one mix (10% over), then add each fragment to its own tube.</div>';
      s.warnings.forEach(x => { h += '<div class="note warn">' + esc(x.charAt(0).toUpperCase() + x.slice(1)) + '.</div>'; });
      h += '<div class="hint">A no-insert control — vector and master mix, no insert — tells you how much of the ' +
        'background is uncut or re-closed vector. It is the one control worth running every time.</div>';
      return h;
    }

    function convOut() {
      const bp = Math.max(1, val('#bnCbp', 1000));
      const mw = B.dnaMW(bp);
      const ng = val('#bnCng', 0), pmol = val('#bnCpmol', 0);
      let h = '<div class="scrollx"><table class="tbl"><thead><tr><th>From</th><th class="num">Value</th><th>Is</th><th class="num">Value</th></tr></thead><tbody>' +
        '<tr><td>' + fmt(bp) + ' bp duplex</td><td class="num mono">' + n2(mw / 1000) + ' kg/mol</td>' +
          '<td>per base pair</td><td class="num mono">' + n2(mw / bp) + ' g/mol</td></tr>' +
        '<tr><td>' + n2(ng) + ' ng</td><td class="num mono">' + n2(B.ngToPmol(ng, bp)) + ' pmol</td>' +
          '<td>or</td><td class="num mono">' + n2(B.ngToFmol(ng, bp)) + ' fmol</td></tr>' +
        '<tr><td>' + n2(pmol) + ' pmol</td><td class="num mono">' + n2(B.pmolToNg(pmol, bp)) + ' ng</td>' +
          '<td>molecules</td><td class="num mono">' + (pmol * 6.02214076e11).toExponential(2) + '</td></tr>' +
        '</tbody></table></div>' +
        '<div class="hint">Molar mass is counted base by base when the fragment is one just designed, and from the ' +
        'average pair (' + B.BP_MW_AVG + ' g/mol) otherwise. The 650 g/mol rule of thumb would read ' +
        n2(ng * 1000 / (bp * 650)) + ' pmol here instead.</div>';
      return h;
    }
    function dilOut() {
      const r = B.dilute(val('#bnDstock', 0), val('#bnDwant', 0), val('#bnDvol', 0));
      if (!r) return '';
      if (!r.ok) return '<div class="note warn">' + esc(r.error) + '</div>';
      return '<div class="note ok"><strong>' + n2(r.stockVol) + ' µL of stock + ' + n2(r.diluentVol) +
        ' µL of diluent</strong> gives ' + n2(val('#bnDvol', 0)) + ' µL at the concentration asked for — a ' +
        n2(r.fold) + '-fold dilution.</div>';
    }
    function oligoOut() {
      const r = B.resuspend(val('#bnOnmol', 0), val('#bnOuM', 0));
      if (!r) return '';
      return '<div class="note ok"><strong>Add ' + n2(r.volUl) + ' µL</strong> of TE or water to the dried pellet for ' +
        n2(val('#bnOuM', 0)) + ' µM. A 10 µM working stock is then 1 part in ' +
        n2(val('#bnOuM', 0) / 10) + '.</div>';
    }

    function run() {
      if (!$$('#bnPcrOut')) return;
      if (!$$('#bnPcr').hidden) $$('#bnPcrOut').innerHTML = pcrOut();
      if (!$$('#bnAsm').hidden) $$('#bnAsmOut').innerHTML = asmOut();
      if (!$$('#bnConv').hidden) {
        $$('#bnConvOut').innerHTML = convOut();
        $$('#bnDilOut').innerHTML = dilOut();
        $$('#bnOligoOut').innerHTML = oligoOut();
      }
    }
    let t = null;
    const soon = () => { clearTimeout(t); t = setTimeout(run, 90); };
    w.querySelectorAll('input, select, textarea').forEach(el => {
      el.addEventListener('input', soon); el.addEventListener('change', soon);
    });
    if (startMode === 'asm') { const c = w.querySelector('[data-bm="asm"]'); if (c) c.click(); }
    run();
  }
  PB.toolBench = toolBench;


  /* ======================= READ A FRAGMENT =======================
     Paste a sequence and see what is in it: the protein, and everything
     somebody put around the protein. Nothing here changes a construct —
     it answers "what am I holding?" before anything is cloned.
  ------------------------------------------------------------------ */
  const ANAT_COLOR = {
    coding: 'var(--accent)', 'fusion partner': '#4db6ac', 'protease site': '#f0a742',
    tag: '#64b5f6', linker: '#b388ff', other: 'var(--fg-faint)'
  };
  /* How much of each region the plasmid already has. A region is reported
     whole, partly or not at all, because a region that is half in the
     vector is the interesting case — it is where the boundary falls, and it
     is where the insert really begins. */
  function segCoverage(v, segments) {
    const overlap = (a, b, c, d) => Math.max(0, Math.min(b, d) - Math.max(a, c));
    return segments.map(g => {
      let vec = 0, amb = 0, diff = 0;
      v.blocks.forEach(b => {
        const o = overlap(g.from, g.to, b.from, b.to);
        if (!o) return;
        if (b.kind === 'vector') vec += o;
        else if (b.kind === 'unreadable') amb += o;
        else if (b.kind === 'differs') diff += o;
      });
      const known = g.bp - amb;
      let cell;
      if (!known) cell = '<span class="badge neutral">not read</span>';
      else if (vec === 0) cell = '<span class="badge acc">new — not in it</span>';
      else if (vec >= known) cell = '<span class="badge ok">already there</span>' +
        (diff ? ' <span class="badge warn">' + diff + ' bp differ</span>' : '');
      else cell = '<span class="badge warn">' + Math.round(vec / known * 100) + '% of it</span>';
      return { vec: vec, ambiguous: amb, differs: diff, known: known, cell: cell,
        all: known > 0 && vec >= known, none: vec === 0 };
    });
  }

  function toolAnatomy(prefill) {
    const rec = PB.active();
    const w = PB.modal({
      title: 'Read a fragment',
      sub: 'what is protein and what was added around it',
      wide: true,
      body:
        '<div class="row">' +
          seqField('anSrc', 'The fragment', { rows: 4, flex: '100%', allowSel: true,
            placeholder: 'paste the sequence — a gBlock, an old insert, anything' }) +
        '</div>' +
        '<div class="row">' +
          '<div class="field narrow" style="align-self:flex-end"><label class="chip on" id="anAuto">find the frame</label></div>' +
          '<div class="field narrow"><label id="anFrameLab">Frame found</label><select class="inp" id="anFrame" disabled>' +
            '<option value="0">1</option><option value="1">2</option><option value="2">3</option></select></div>' +
          '<div class="field narrow"><label id="anStrandLab">Strand found</label><select class="inp" id="anStrand" disabled>' +
            '<option value="1">forward</option><option value="-1">reverse</option></select></div>' +
          '<div class="field" style="flex:1 1 230px;align-self:flex-end"><div class="hint" id="anFrameWhy" style="margin:0"></div></div>' +
        '</div>' +
        '<div id="anOut"></div>'
    });
    const $$ = id => w.querySelector(id);
    let auto = true;
    const ac = $$('#anAuto');
    ac.addEventListener('click', () => {
      auto = !auto; ac.classList.toggle('on', auto);
      /* while it is finding the frame these two report an answer; they only
         take one once that has been turned off */
      $$('#anFrame').disabled = auto; $$('#anStrand').disabled = auto;
      $$('#anFrameLab').textContent = auto ? 'Frame found' : 'Frame';
      $$('#anStrandLab').textContent = auto ? 'Strand found' : 'Strand';
      run();
    });
    if (prefill) { $$('#anSrc_src').value = ''; $$('#anSrc_txt').value = prefill; }

    function run() {
      const host = $$('#anOut');
      if (!host) return;
      const v = seqValue('anSrc');
      if (!v.seq || v.seq.length < 24) {
        host.innerHTML = '<div class="hint">Paste at least 24 bases.</div>'; return;
      }
      const common = PB.allEnzymes().filter(e => PB.COMMON.has(e.name));
      let a;
      try {
        a = B.fragmentAnatomy(v.seq, Object.assign({ enzymes: common, table: state.codeTable },
          auto ? {} : { frame: +$$('#anFrame').value, strand: +$$('#anStrand').value }));
      } catch (e) { host.innerHTML = '<div class="note bad">' + esc(String(e.message || e)) + '</div>'; return; }
      if (!a.ok) { host.innerHTML = '<div class="note">' + esc(a.error) + '</div>'; return; }
      if (auto) { $$('#anFrame').value = String(a.frame); $$('#anStrand').value = String(a.strand); }
      /* why this frame and not another: the evidence, not the number */
      const why = $$('#anFrameWhy');
      if (why) {
        const bits = [];
        if (a.startsWithMet) bits.push('starts at an ATG');
        if (a.endsWithStop) bits.push('ends at a stop');
        bits.push(fmt(a.core.residues) + ' residues unbroken');
        if (a.motifs.length) bits.push('and ' + esc(a.motifs[0].name.replace(/ \(partial\)/, '')) + ' reads out of it');
        why.innerHTML = auto
          ? 'Chosen because it ' + bits.join(', ') + '. Turn the chip off to read it in another frame.'
          : 'Reading it as you asked. Turn the chip on to let it choose.';
      }

      /* Matched against the construct on screen, when there is one. A read
         off a sequencing primer begins in the vector, so a good part of
         what this table describes is already in the plasmid and is not
         something to clone into it again. */
      let vecMatch = null;
      if (rec && rec.seq && rec.seq.length > 40) {
        try { vecMatch = B.readInVector(a.work, rec, {}); } catch (e) { vecMatch = null; }
      }

      const added = a.segments.filter(g => g.kind !== 'coding').reduce((x, g) => x + g.bp, 0);
      const main = a.main || a.core;
      const fusion = a.partners && a.partners.length > 0;
      let h = '<div class="note ' + (added ? 'warn' : 'ok') + '"><strong>' + fmt(a.length) + ' bp: ' +
        (fusion
          ? 'a fusion — ' + a.partners.map(p => esc(p.name) + ' (' + fmt(p.residues) + ' residues)').join(' and ') +
            ' in front of a protein of ' + fmt(Math.floor(main.bp / 3)) + ' residues'
          : fmt(main.bp) + ' bp of protein (' + fmt(a.core.residues) + ' residues)' +
            (added ? ' and ' + fmt(added) + ' bp of other things around it' : ' and nothing else')) +
        '.</strong> ' +
        (a.strand === -1 ? 'It reads on the reverse strand — the sequence as pasted is the complement of the gene. ' : '') +
        'The reading frame ' + (a.startsWithMet ? 'starts at an ATG' : 'does not start at an ATG') + ' and ' +
        (a.endsWithStop ? 'ends at a stop codon' : 'has no stop codon at its end') + '.</div>';
      if (fusion) h += '<div class="note"><strong>These are two proteins in one reading frame.</strong> ' +
        a.partners.map(p => esc(p.name) + ' is there for ' + esc(p.why) +
          (p.mismatches ? ' (' + p.mismatches + ' residue' + (p.mismatches === 1 ? '' : 's') + ' differ from the reference)' : '')).join('; ') +
        '. ' + (a.segments.some(g => g.kind === 'protease site')
          ? 'The protease site between them is where they come apart after purification.'
          : 'There is no protease site between them, so they stay joined — check that is what you want.') +
        ' “Clone just the protein” takes the passenger and leaves the partner behind.</div>';
      if (!fusion && a.proteinParts === 1 && a.core.residues > 150 && !a.motifs.some(m => m.inside))
        h += '<div class="hint">It reads as one protein. If something is fused in front of it — a solubility ' +
          'partner, a domain from another gene — it is not one this tool recognises; the ones it knows are ' +
          esc(B.FUSION_PARTNERS.map(f => f.name.replace(/ \(.*/, '')).join(', ')) +
          '. Read the first residues of the translation below and see whether they are yours.</div>';

      /* ---- what of this the plasmid already has ---- */
      if (vecMatch) {
        const v = vecMatch;
        const pc = n => Math.round(n / v.length * 100);
        const seg = v.newSpan ? a.segments.filter(g => g.to > v.newSpan.from && g.from < v.newSpan.to) : [];
        h += '<div class="field wide" style="margin-top:10px"><label>What of this is already in ' +
          esc(v.vectorName) + '</label>' +
          '<div style="display:flex;height:22px;border-radius:5px;overflow:hidden;border:1px solid var(--border)">' +
          v.blocks.map(b => '<div title="' + esc((b.from + 1) + '–' + b.to + ' (' + b.bp + ' bp) ' +
              (b.kind === 'vector' ? 'matches ' + v.vectorName + ' ' + (b.vecFrom + 1) + '–' + b.vecTo +
                (b.strand === -1 ? ' on the bottom strand' : '')
               : b.kind === 'new' ? 'not in ' + v.vectorName
               : b.kind === 'differs' ? 'differs from ' + v.vectorName
               : b.kind === 'edge' ? 'ragged end of the read' : 'too ambiguous to tell')) +
            '" style="width:' + (b.bp / v.length * 100).toFixed(3) + '%;background:' +
            ({ vector: 'var(--good)', new: 'var(--f-promoter)', differs: 'var(--bad)',
               edge: 'var(--fg-faint)', unreadable: 'var(--surface-3)' }[b.kind]) +
            ';opacity:' + (b.kind === 'unreadable' || b.kind === 'edge' ? '.4' : '.85') + '"></div>').join('') +
          '</div><div class="hint" style="margin:4px 0 0">' +
          [['var(--good)', 'already in the plasmid'], ['var(--f-promoter)', 'not in it'],
           ['var(--bad)', 'differs'], ['var(--surface-3)', 'too ambiguous to tell']]
            .filter(x => v.blocks.some(b => ({ 'var(--good)': 'vector', 'var(--f-promoter)': 'new',
              'var(--bad)': 'differs', 'var(--surface-3)': 'unreadable' })[x[0]] === b.kind))
            .map(x => '<span style="display:inline-block;width:9px;height:9px;border-radius:2px;background:' +
              x[0] + ';margin:0 4px 0 10px"></span>' + x[1]).join('') + '</div></div>';

        if (!v.newSpan) {
          h += '<div class="note ok"><strong>All of this is already in ' + esc(v.vectorName) + '.</strong> ' +
            fmt(v.inVector) + ' of ' + fmt(v.length) + ' bases match it' +
            (v.differs ? ', apart from ' + fmt(v.differs) + ' bp that differ' : '') +
            (v.unreadable ? ' (' + fmt(v.unreadable) + ' bp too ambiguous to tell)' : '') +
            '. There is nothing here to add to it — this read is of the construct you already have.</div>';
        } else {
          const s = v.newSpan;
          h += '<div class="note ' + (v.inVector ? 'warn' : '') + '"><strong>Bases ' + fmt(s.from + 1) + '–' +
            fmt(s.to) + ' (' + fmt(s.bp) + ' bp) are not in ' + esc(v.vectorName) + '.</strong> ' +
            (v.inVector
              ? 'The ' + fmt(v.inVector) + ' bases in front of that (' + pc(v.inVector) + '% of the read) are ' +
                (v.features.length
                  ? 'the plasmid’s own ' + esc(v.features.slice(0, 5).join(', ')) +
                    (v.features.length > 5 ? ' and others' : '')
                  : 'already in the plasmid') +
                ' — the read starts inside the vector, as a read off a sequencing primer does, so that part is ' +
                'not something to put into it again. '
              : 'None of this read matches the plasmid on screen. ') +
            (seg.length
              ? 'What is new covers ' + seg.map(g => esc(g.name)).filter((x, i, arr) => arr.indexOf(x) === i).join(', ') + '. '
              : '') +
            (!s.ambiguous
              ? 'Every base of it was read.'
              : s.tailAmbiguous
                ? '<strong>Only the first ' + fmt(s.cleanBp) + ' bases of it are read cleanly.</strong> ' +
                  'After base ' + fmt(s.from + s.cleanBp + 1) + ' of the fragment there are ' + fmt(s.ambiguous) +
                  ' N in the remaining ' + fmt(s.bp - s.cleanBp) + ' bp, so the far end of the new part is not ' +
                  'established by this read at all — a missing stop codon or a changed last codon would not show ' +
                  'here. Read it from the other end before designing anything against that end.'
                : '<strong>' + fmt(s.ambiguous) + ' bases of it (' + Math.round(s.ambiguous / s.bp * 100) +
                  '%) are N</strong>, scattered through it — this read does not establish what they are, so the ' +
                  'new part is not fully determined by it. Read it again before ordering anything built on it.') +
            '</div>';
          h += '<div class="row" style="margin-top:8px">' +
            '<button class="btn primary" id="anNewClone">Clone the part it does not have (' + fmt(s.bp) + ' bp)</button>' +
            '<button class="btn" id="anNewCopy">Copy those bases</button></div>' +
            '<div class="hint">That opens the assembly designer with bases ' + fmt(s.from + 1) + '–' + fmt(s.to) +
            ' as the insert and ' + esc(v.vectorName) + ' as the backbone, where the placement list offers its ' +
            'cloning sites. It is the span to add; where in the plasmid it should go is still yours to choose.</div>';
        }
        if (v.differs) h += '<div class="note warn"><strong>' + fmt(v.differs) + ' bp differ from ' +
          esc(v.vectorName) + '</strong> inside stretches that otherwise match it. A base that is read ' +
          'differently blanks out the window around it, so a single change shows as a gap of up to ' + v.k +
          ' bp here. Use Verify by sequencing to see the changes base by base.</div>';
      } else if (rec && rec.seq) {
        h += '<div class="hint" style="margin-top:8px">Load a construct to see which of this the plasmid ' +
          'already has.</div>';
      }

      /* the track */
      h += '<div class="field wide" style="margin-top:8px"><label>What is where</label>' +
        '<div style="display:flex;height:26px;border-radius:5px;overflow:hidden;border:1px solid var(--line)">' +
        a.segments.map(g => '<div title="' + esc(g.name + ' — ' + (g.from + 1) + '–' + g.to + ' (' + g.bp + ' bp)') +
          '" style="width:' + (g.bp / a.length * 100).toFixed(3) + '%;background:' +
          (ANAT_COLOR[g.kind] || ANAT_COLOR.other) + ';opacity:' + (g.kind === 'other' ? '.35' : '.85') +
          ';border-right:1px solid rgba(0,0,0,.35)"></div>').join('') +
        '</div><div class="hint" style="margin:4px 0 0">' +
        Object.keys(ANAT_COLOR).filter(k => a.segments.some(g => g.kind === k)).map(k =>
          '<span style="display:inline-block;width:9px;height:9px;border-radius:2px;background:' + ANAT_COLOR[k] +
          ';margin:0 4px 0 10px"></span>' + esc(k)).join('') + '</div></div>';

      /* What each region is, and — on asking — the bases it is made of.
         A table that says a region is 168 bp of GB1 and then shows twelve
         residues of it is a description; the thing a person came for is the
         sequence itself, to paste into an order or a primer. Every row can
         open to its own DNA and its own translation, each copyable on its
         own, so nothing has to be counted out of the pasted fragment by
         hand. */
      const inVec = vecMatch ? segCoverage(vecMatch, a.segments) : null;
      h += '<div class="scrollx" style="margin-top:8px"><table class="tbl" id="anTable"><thead><tr><th>Bases</th>' +
        '<th class="num">bp</th><th>What it is</th><th>Reads as</th>' +
        (inVec ? '<th>In ' + esc(PB.trimText(vecMatch.vectorName, 22)) + '?</th>' : '') +
        '<th></th></tr></thead><tbody>' +
        a.segments.map((g, gi) => '<tr data-anrow="' + gi + '"><td class="mono">' + fmt(g.from + 1) + '–' + fmt(g.to) + '</td>' +
          '<td class="num mono">' + fmt(g.bp) + '</td>' +
          '<td>' + esc(g.name) + (g.kind !== 'other' && g.kind !== 'coding' ?
            ' <span class="badge neutral">' + esc(g.kind) + '</span>' : '') + '</td>' +
          '<td class="mono" style="font-size:11.5px">' + esc(!g.aa ? '—'
            : (g.aa.length > 20 ? g.aa.slice(0, 12) + '…' + g.aa.slice(-6) : g.aa)) + '</td>' +
          (inVec ? '<td style="font-size:11.5px">' + inVec[gi].cell + '</td>' : '') +
          '<td style="white-space:nowrap">' +
            '<button class="btn sm ghost" data-anshow="' + gi + '">show</button>' +
            (g.kind === 'coding' ? '' :
            '<button class="btn sm ghost" data-ancut="' + g.from + ',' + g.to + '">drop</button>') +
          '</td></tr>' +
          '<tr class="anseq" data-anseq="' + gi + '" hidden><td colspan="' + (inVec ? 6 : 5) + '">' +
            '<div class="field wide"><label>' + esc(g.name) + ' · bases ' + fmt(g.from + 1) + '–' + fmt(g.to) +
              ' (' + fmt(g.bp) + ' bp)</label>' +
              '<div class="seqblock" style="max-height:120px">' + esc(g.dna || '') + '</div></div>' +
            (g.aa ? '<div class="field wide" style="margin-top:6px"><label>as protein (' +
              fmt(g.aa.replace(/\*/g, '').length) + ' residues)</label>' +
              '<div class="seqblock" style="max-height:90px">' + esc(g.aa) + '</div></div>' : '') +
            '<div class="row" style="margin-top:6px">' +
              '<button class="btn sm" data-ancopy="' + gi + '">Copy the bases</button>' +
              (g.aa ? '<button class="btn sm" data-ancopyaa="' + gi + '">Copy the protein</button>' : '') +
              '<button class="btn sm ghost" data-anonly="' + gi + '">Clone just this part</button>' +
            '</div></td></tr>').join('') +
        '</tbody></table></div>' +
        '<div class="hint">“Show” opens a region’s own bases and its own translation. ' +
        '“Clone just this part” sends that region alone to the assembly designer.</div>';

      /* motifs sitting inside the protein are worth knowing about too */
      const inside = a.motifs.filter(m => m.inside);
      if (inside.length) h += '<div class="note"><strong>Inside the protein itself:</strong> ' +
        inside.map(m => esc(m.name) + ' (' + esc(m.aa) + ') at residue ' + fmt(m.aaFrom + 1)).join(', ') +
        '. These are part of the sequence as designed — a tag meant to stay, or a coincidence of the protein.</div>';

      const outSites = a.sites.filter(x => !x.inCore);
      if (outSites.length) h += '<div class="note"><strong>Restriction sites outside the protein:</strong> ' +
        outSites.map(x => esc(x.name) + ' at ' + fmt(x.at + 1)).join(', ') +
        '. Usually how the fragment was meant to be cloned before.</div>';

      h += '<div class="field wide" style="margin-top:8px"><label>' +
        (fusion ? 'The passenger protein (' : 'The protein it codes for (') +
        fmt(Math.floor(main.bp / 3)) + ' residues)</label><div class="seqblock" style="max-height:140px">' +
        esc(main.aa || a.core.aa) + '</div></div>' +
        (fusion ? '<div class="field wide"><label>The partner in front of it (' +
          fmt(a.partners[0].residues) + ' residues)</label><div class="seqblock" style="max-height:110px">' +
          esc(a.protein.slice(a.partners[0].aaFrom, a.partners[0].aaTo)) + '</div></div>' : '');

      h += '<div class="row" style="margin-top:10px">' +
        '<button class="btn primary" id="anClone">Clone just the protein</button>' +
        '<button class="btn" id="anSave">Save it annotated</button>' +
        '<button class="btn" id="anCopy">Copy the protein’s DNA</button></div>' +
        '<div class="hint">“Clone just the protein” opens the assembly designer with bases ' +
        fmt(main.from + 1) + '–' + fmt(main.to) + ' as the insert, leaving the rest behind.</div>';
      host.innerHTML = h;

      const core = a.work.slice(main.from, main.to);
      host.querySelectorAll('[data-ancut]').forEach(b => b.addEventListener('click', () => {
        const [f, t] = b.getAttribute('data-ancut').split(',').map(Number);
        $$('#anSrc_src').value = ''; $$('#anSrc_txt').value = a.work.slice(0, f) + a.work.slice(t);
        run();
      }));
      /* a region's own bases, on asking */
      host.querySelectorAll('[data-anshow]').forEach(b => b.addEventListener('click', () => {
        const gi = b.getAttribute('data-anshow');
        const row = host.querySelector('[data-anseq="' + gi + '"]');
        if (!row) return;
        row.hidden = !row.hidden;
        b.textContent = row.hidden ? 'show' : 'hide';
      }));
      host.querySelectorAll('[data-ancopy]').forEach(b => b.addEventListener('click', () => {
        const g = a.segments[+b.getAttribute('data-ancopy')];
        if (g) PB.copyText(g.dna, g.name + ' (' + g.bp + ' bp)');
      }));
      host.querySelectorAll('[data-ancopyaa]').forEach(b => b.addEventListener('click', () => {
        const g = a.segments[+b.getAttribute('data-ancopyaa')];
        if (g && g.aa) PB.copyText(g.aa, g.name + ' as protein');
      }));
      host.querySelectorAll('[data-anonly]').forEach(b => b.addEventListener('click', () => {
        const g = a.segments[+b.getAttribute('data-anonly')];
        if (!g) return;
        state.gbPrefillInsert = { seq: g.dna, name: (g.name === 'protein' ? (v.name || 'fragment') : g.name) };
        PB.closeModal(); toolGibson('infusion');
      }));
      if (vecMatch && vecMatch.newSpan) {
        const s = vecMatch.newSpan;
        const nc = host.querySelector('#anNewClone');
        if (nc) nc.addEventListener('click', () => {
          state.gbPrefillInsert = { seq: s.seq, name: (v.name || 'fragment') + ' insert' };
          PB.closeModal(); toolGibson('infusion');
        });
        const ncp = host.querySelector('#anNewCopy');
        if (ncp) ncp.addEventListener('click', () =>
          PB.copyText(s.seq, 'The part ' + vecMatch.vectorName + ' does not have (' + s.bp + ' bp)'));
      }
      $$('#anClone').addEventListener('click', () => {
        state.gbPrefillInsert = { seq: core, name: (v.name || 'fragment') + ' CDS' };
        PB.closeModal(); toolGibson('infusion');
      });
      $$('#anCopy').addEventListener('click', () => PB.copyText(core, 'Coding sequence'));
      $$('#anSave').addEventListener('click', () => {
        PB.saveProduct({
          name: (v.name || 'fragment') + ' read',
          seq: a.work, circular: false,
          features: a.segments.filter(g => g.kind !== 'other').map(g => PB.normFeature({
            start: g.from, end: g.to, strand: 1,
            type: g.kind === 'coding' ? 'CDS' : 'misc_feature',
            label: g.kind === 'coding' ? 'protein' : g.name,
            color: g.kind === 'coding' ? '#b1ff67' : '#f0a742',
            quals: g.kind === 'coding' ? { codon_start: '1', translation: a.core.aa } : { note: g.kind }
          })),
          history: [['Fragment read', fmt(a.core.residues) + ' residues in frame ' + (a.frame + 1) +
            (a.strand === -1 ? ' on the reverse strand' : '') + '; ' +
            (a.segments.length - 1) + ' other segment(s) identified']]
        });
      });
    }
    let t = null;
    const soon = () => { clearTimeout(t); t = setTimeout(run, 110); };
    w.querySelectorAll('input, select, textarea').forEach(el => {
      el.addEventListener('input', soon); el.addEventListener('change', soon);
    });
    run();
  }
  PB.toolAnatomy = toolAnatomy;

})(window.PB);

