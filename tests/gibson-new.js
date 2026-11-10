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
      title: isInf ? 'In-Fusion cloning' : 'Gibson / NEBuilder HiFi assembly', wide: true,
      sub: isInf ? 'Takara seamless cloning, 15 bp homology arms' : 'exonuclease-based seamless assembly',
      body: '<div class="chipbar"><button class="chip on" data-gm="design">Design primers</button>' +
          '<button class="chip" data-gm="assemble">Assemble existing fragments</button></div>' +
        '<div id="gbDesign">' +
          '<div class="row">' +
            seqField('gbIns', 'Insert — what goes in', { allowSel: true, rows: 2, flex: '260px' }) +
            '<div class="field" style="flex:0 1 210px"><label>Which part of it</label>' +
              '<select class="inp" id="gbInsPart"></select>' +
              '<input class="inp" id="gbInsName" placeholder="name it, e.g. CDK1" spellcheck="false" style="margin-top:5px"></div>' +
            seqField('gbVec', 'Vector — what it goes into', { defaultId: rec ? rec.id : '', allowSel: true, rows: 2, flex: '260px' }) +
          '</div>' +
          '<div class="row" style="align-items:flex-end">' +
            '<div class="field" style="flex:1 1 420px"><label>Where it goes</label>' +
              '<select class="inp" id="gbWhere"></select></div>' +
            '<div class="field narrow" id="gbFromWrap" hidden><label>From (bp)</label><input class="inp mono" id="gbFrom" type="number" min="1"></div>' +
            '<div class="field narrow" id="gbToWrap" hidden><label>To (bp)</label><input class="inp mono" id="gbTo" type="number" min="1" placeholder="blank"></div>' +
            '<div class="field narrow"><label>Arms (bp)</label><input class="inp mono" id="gbOv" type="number" value="' + defOv + '" min="15" max="60"></div>' +
            '<div class="field narrow"><label>Anneal (bp)</label><input class="inp mono" id="gbAn" type="number" value="20" min="15" max="40"></div>' +
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
        let circ = true;
        w.querySelectorAll('[data-gm]').forEach(c => c.addEventListener('click', () => {
          w.querySelectorAll('[data-gm]').forEach(x => x.classList.toggle('on', x === c));
          const d = c.getAttribute('data-gm') === 'design';
          $('#gbDesign').hidden = !d; $('#gbAssemble').hidden = d;
          if (d) run();
        }));
        const cEl = w.querySelector('#gbCirc');
        cEl.addEventListener('click', () => { circ = !circ; cEl.classList.toggle('on', circ); });
        w.querySelector('#gbFromLib').addEventListener('click', () => {
          const ta = $('#gbFrags');
          ta.value = state.records.slice(0, 4).map(r => '>' + r.name.replace(/\s+/g, '_') + '\n' + r.seq).join('\n');
        });

        /* ---------- what goes in ----------
           A record usually holds more than the thing being cloned, so its
           coding features are offered by name: picking "CDK1" clones the
           gene, not the plasmid it currently sits in. */
        function insertParts() {
          const v = seqValue('gbIns');
          const r = v.rec || v.from;
          const out = [{ label: 'the whole sequence (' + fmt(v.seq.length) + ' bp)', seq: v.seq, name: v.name }];
          if (r && r.features && v.rec) {
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
          parts = insertParts();
          const sel = $('#gbInsPart');
          const prev = sel.value;
          sel.innerHTML = parts.map((p, i) => '<option value="' + i + '">' + esc(p.label) + '</option>').join('');
          if (prev && parts[+prev]) sel.value = prev;
          /* a coding feature is almost always what someone means */
          else { const cds = parts.findIndex(p => p.feature && /^CDS$/i.test(p.feature.type)); if (cds > 0) sel.value = String(cds); }
          if (!keepName) autoName();
        }
        function insertPart() { return parts[+$('#gbInsPart').value] || parts[0] || { seq: '', name: 'insert' }; }
        function autoName() {
          const n = $('#gbInsName');
          if (n && !n.dataset.touched) n.value = (insertPart().name || 'insert').replace(/\s+/g, ' ').trim();
        }
        $('#gbInsName').addEventListener('input', e => { e.target.dataset.touched = '1'; });

        /* ---------- where it goes ----------
           One list of real choices, built from the vector itself: the
           placements its own annotation implies, every enzyme that cuts it
           once, whatever is selected on the map, and a free range. */
        let places = [], placeCache = '';
        function placements(rec, seq, circular) {
          const out = [];
          if (!circular) return out;
          cloningSites(rec).forEach(x => out.push({ kind: 'ipcr', label: x.label, from: x.from, to: x.to }));
          const meth = rec ? B.methylatedBases(seq, true, PB.methylSystems(rec)) : null;
          const cutters = [];
          PB.allEnzymes().forEach(e => {
            if (!PB.COMMON.has(e.name) && !e.custom) return;
            let sites = [];
            try { sites = B.findSites(seq, e, true); } catch (err) { return; }
            if (sites.length !== 1) return;
            const m = meth ? B.methylationEffect(e, sites[0], seq.length, meth) : null;
            if (m && m.effect === 'blocked') return;
            cutters.push({ kind: 'enz', name: e.name, pos: sites[0].cutTop,
                           label: 'cut with ' + e.name + ' (single site, cuts at ' + fmt(sites[0].cutTop) + ')' +
                                  (m && m.effect ? ' ⚠ ' + m.effect : '') });
          });
          cutters.sort((a, b) => a.pos - b.pos).forEach(c => out.push(c));
          const sel = PB.normSel(rec, state.selection);
          if (rec && sel && PB.active() && PB.active().id === rec.id) {
            out.push({ kind: 'ipcr', label: 'replace what is selected on the map (' +
              PB.selRangeText(rec, state.selection) + ', ' + fmt(sel.length) + ' bp)',
              from: sel.start + 1, to: (sel.end - 1) % seq.length + 1 });
          }
          out.push({ kind: 'custom', label: 'a range I type — inverse PCR, removing it' });
          out.push({ kind: 'pos', label: 'a single base I type — blunt opening, nothing removed' });
          return out;
        }
        function refreshPlaces() {
          const v = seqValue('gbVec');
          const r = v.rec || v.from;
          const key = (v.rec ? v.rec.id + ':' + v.rec.updatedAt : 'x') + ':' + v.seq.length + ':' +
                      (v.circular ? 'c' : 'l') + ':' + JSON.stringify(state.selection);
          const sel = $('#gbWhere');
          if (key === placeCache) return;
          placeCache = key;
          const prev = sel.value;
          places = placements(r, v.seq, v.circular);
          if (!v.circular) {
            sel.innerHTML = '<option value="lin">this vector is already linear — use it as it is</option>';
            places = [{ kind: 'linear', label: 'already linear' }];
            sel.value = '0';
            sel.innerHTML = '<option value="0">this vector is already linear — use it as it is</option>';
          } else {
            sel.innerHTML = places.map((p, i) => '<option value="' + i + '">' + esc(p.label) + '</option>').join('');
            if (prev && places[+prev]) sel.value = prev;
          }
          applyPlace();
        }
        function place() { return places[+$('#gbWhere').value] || places[0] || null; }
        function applyPlace() {
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
          if (!circular) return { ok: true, seq: seq, how: 'supplied linear', points: [], primers: null };
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
          let from, to;
          if (p.kind === 'pos') { from = clampPos($('#gbFrom').value, L); to = from; }
          else if (p.kind === 'custom') {
            from = clampPos($('#gbFrom').value, L);
            const t = ($('#gbTo').value || '').trim();
            to = t ? Math.max(0, Math.min(L, parseInt(t, 10) || 0)) : from;
            if (to < from) to = from;
          } else { from = p.from; to = p.to === '' ? p.from : p.to; }
          const a0 = from - 1, b0 = to >= from ? to : from - 1;
          if (p.kind === 'pos' || b0 <= a0) {
            const r2 = B.inversePCR(seq, a0, a0, { targetTm: 60, conditions: PB.conditions ? PB.conditions() : undefined });
            if (!r2.ok) return { ok: false, error: r2.error };
            return { ok: true, seq: r2.backbone.seq, how: 'opened before base ' + fmt(from) + ', nothing removed',
                     points: [a0], removed: null, primers: r2.primers, openAt: a0, ipcr: true };
          }
          const r3 = B.inversePCR(seq, a0, b0, { targetTm: 60, conditions: PB.conditions ? PB.conditions() : undefined });
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
            out.push('<div class="note ' + (inFrame ? 'ok' : 'bad') + '"><strong>' +
              esc(f.label || 'CDS') + ' → ' + esc(($('#gbInsName').value || 'insert')) + ': ' +
              (inFrame ? 'in frame' : 'out of frame by ' + (d % 3) + ' base' + (d % 3 === 1 ? '' : 's')) +
              '.</strong> ' + (inFrame
                ? 'The tag reads straight into the insert:'
                : 'Shift the opening by ' + (3 - d % 3) + ' base' + (3 - d % 3 === 1 ? '' : 's') +
                  ', or add ' + (3 - d % 3) + ' to the forward tail, or the fusion is nonsense.') +
              peptide + '</div>');
          });
          return out.slice(0, 3).join('');
        }

        /* What the insert itself carries into the fusion. */
        function insertNotes(insSeq, inFrameFusion) {
          const n = [];
          const prot = B.translate(insSeq, { table: state.codeTable });
          const trail = /\*$/.test(prot);
          const internal = prot.replace(/\*$/, '').indexOf('*') >= 0;
          if (insSeq.length % 3 !== 0)
            n.push('The insert is ' + fmt(insSeq.length) + ' bp, not a whole number of codons (' +
              (insSeq.length % 3) + ' over) — anything fused after it is out of frame.');
          if (internal) n.push('It contains a stop codon in its own frame before the end.');
          if (trail) n.push('It ends in a stop codon, so a C-terminal tag on the vector will not be translated.');
          else n.push('It has no stop codon at the end, so translation runs on into the vector — fine for a C-terminal tag, not for an untagged protein.');
          if (/^ATG/.test(insSeq) && inFrameFusion)
            n.push('It starts with ATG; after an N-terminal tag that is an extra methionine, which is usual and harmless.');
          return n.length ? '<div class="note"><strong>The insert:</strong> ' + n.map(esc).join(' ') + '</div>' : '';
        }

        /* ---------- the whole design, recomputed on any change ---------- */
        let lastDesign = null;
        function run() {
          if ($('#gbDesign').hidden) return;
          refreshParts(true); refreshPlaces();
          const host = $('#gbLive');
          const v = seqValue('gbVec');
          const part = insertPart();
          const insSeq = B.cleanSeq(part.seq || '');
          lastDesign = null;
          if (!v.seq.length || !insSeq.length) {
            host.innerHTML = '<div class="hint">Pick an insert and a vector. Everything below — the opened backbone, both primer pairs, the junctions and the finished plasmid — follows from those and from where it goes.</div>';
            return;
          }
          const r = v.rec || v.from;
          const ov = Math.max(15, Math.min(60, parseInt($('#gbOv').value, 10) || defOv));
          const an = Math.max(15, Math.min(40, parseInt($('#gbAn').value, 10) || 20));
          const poi = ($('#gbInsName').value || '').trim() || part.name || 'insert';

          const open = openVector(r, v.seq, v.circular);
          if (!open.ok) { host.innerHTML = '<div class="note bad">' + esc(open.error) + '</div>'; return; }
          const vecSeq = open.seq;

          const des = B.designAssemblyPrimers(vecSeq, [{ name: poi, seq: insSeq }], { overlap: ov, anneal: an, method: method });
          const p = des.primers[0];
          const amps = B.pcr(insSeq, p.forward.seq, p.reverse.seq, { anchor: Math.min(18, an) });
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
          const frame = fusionReport(r, open, insSeq, ov);
          html += frame;
          html += insertNotes(insSeq, /in frame<\/strong>/.test(frame));
          /* 3 — the primers to order */
          const rows = [];
          rows.push(['Insert forward', poi + '_F', p.forward.seq, p.forward.tm, p.forward.tail.length]);
          rows.push(['Insert reverse', poi + '_R', p.reverse.seq, p.reverse.tm, p.reverse.tail.length]);
          if (open.primers) {
            rows.push(['Backbone forward', poi + '_bb_F', open.primers.forward.seq, open.primers.forward.tm, 0]);
            rows.push(['Backbone reverse', poi + '_bb_R', open.primers.reverse.seq, open.primers.reverse.tm, 0]);
          }
          html += '<div class="scrollx" style="margin-top:8px"><table class="tbl"><thead><tr><th>Primer</th><th>Name</th>' +
            '<th>Sequence (5′→3′)</th><th class="num">nt</th><th class="num">Tm</th><th class="num">Tail</th></tr></thead><tbody>' +
            rows.map(x => '<tr><td>' + esc(x[0]) + '</td><td class="mono" style="font-size:11.5px">' + esc(x[1]) + '</td>' +
              '<td class="seqcell mono" style="font-size:11.5px;max-width:340px">' +
                (x[4] ? '<span style="color:var(--accent)">' + esc(x[2].slice(0, x[4])) + '</span>' + esc(x[2].slice(x[4])) : esc(x[2])) +
              '</td><td class="num">' + x[2].length + '</td><td class="num">' + (isFinite(x[3]) ? f1(x[3]) : '—') + '</td>' +
              '<td class="num">' + (x[4] || '—') + '</td></tr>').join('') +
            '</tbody></table></div>' +
            '<div class="hint">Coloured bases are the homology tail, which does not anneal to the template in the first cycles; the Tm shown is of the annealing part. ' +
            (open.ipcr ? 'Amplify the backbone with its pair and DpnI the template before assembly. ' : '') +
            'Order all four and the reaction is one tube.</div>';
          /* 4 — the junctions */
          if (asm.ok) {
            html += junctionBox('5′ junction (vector → insert)',
                '<span style="color:var(--fg-faint)">' + esc(vecSeq.slice(-15)) + '</span>' +
                '<span style="color:var(--accent);font-weight:600">' + esc(insSeq.slice(0, 15)) + '</span>') +
              junctionBox('3′ junction (insert → vector)',
                '<span style="color:var(--accent);font-weight:600">' + esc(insSeq.slice(-15)) + '</span>' +
                '<span style="color:var(--fg-faint)">' + esc(vecSeq.slice(0, 15)) + '</span>');
          }
          /* 5 — what could go wrong */
          if (amps.length > 1) html += '<div class="note warn"><strong>The insert primers prime at ' + amps.length +
            ' places on the insert</strong> (' + amps.map(x => fmt(x.length)).join(', ') + ' bp). Gel-purify the ' +
            fmt(main.length) + ' bp band, or move the primers.</div>';
          if (main && main.length !== want) html += '<div class="note bad">The full-length product (' + fmt(want) +
            ' bp) was not predicted; what is shown uses a ' + fmt(main.length) + ' bp product and is probably not what you want.</div>';
          if (open.primers) html += primerPairNote(v.seq, open.primers);
          html += openingReport(r, open.points, open.removed);
          /* 6 — take it away */
          if (asm.ok) {
            html += '<div class="row" style="margin-top:10px">' +
              '<button class="btn primary" id="gbSave">Save the construct</button>' +
              '<button class="btn" id="gbSavePrimers">Save the primers</button>' +
              '<button class="btn" id="gbCopyPrimers">Copy primers</button></div>';
          }
          host.innerHTML = html;

          lastDesign = { v: v, r: r, open: open, vecSeq: vecSeq, insSeq: insSeq, p: p, asm: asm, ov: ov, poi: poi, rows: rows, part: part };
          const sv = $('#gbSave');
          if (sv) sv.addEventListener('click', saveIt);
          const sp = $('#gbSavePrimers');
          if (sp) sp.addEventListener('click', savePrimers);
          const cp = $('#gbCopyPrimers');
          if (cp) cp.addEventListener('click', () => PB.copyText(rows.map(x => x[1] + '\t' + x[2]).join('\n'), rows.length + ' primers'));
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
          d.rows.forEach(x => state.primers.push({ name: x[1].replace(/[^\w.-]+/g, '_'), seq: x[2] }));
          PB.lsSet('primers', state.primers);
          PB.toast('Saved ' + d.rows.length + ' primers to the library', 'ok');
        }

        /* every input feeds the same recompute */
        let t = null;
        const soon = () => { clearTimeout(t); t = setTimeout(run, 110); };
        ['#gbIns_src', '#gbIns_txt', '#gbVec_src', '#gbVec_txt', '#gbInsPart', '#gbWhere',
         '#gbFrom', '#gbTo', '#gbOv', '#gbAn', '#gbInsName'].forEach(id => {
          const el = w.querySelector(id);
          if (!el) return;
          el.addEventListener('input', soon);
          el.addEventListener('change', soon);
        });
        w.querySelector('#gbIns_src').addEventListener('change', () => { $('#gbInsPart').value = ''; refreshParts(false); soon(); });
        w.querySelector('#gbVec_src').addEventListener('change', () => { placeCache = ''; $('#gbFrom').value = ''; soon(); });
        w.querySelector('#gbWhere').addEventListener('change', () => { $('#gbFrom').value = ''; applyPlace(); soon(); });
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
                     from: p + 1, to: '', keeps: inside.map(g => g.label || g.type).join(', ') });
        }
      }
      if (atg) {
        out.push({ label: name + ' · replace from the ATG (untagged, keeps the vector’s RBS spacing)',
                   from: atg.start + 1, to: f.end, keeps: 'RBS and promoter' });
      }
      out.push({ label: name + ' · replace the whole site (' + fmt(f.start + 1) + '–' + fmt(f.end) + ')',
                 from: f.start + 1, to: f.end, keeps: 'RBS and promoter' });
    });
    return out;
  }
  PB.cloningSites = cloningSites;

  /* What an opening lands in. A vector is opened to receive an insert, so the
     one thing worth saying out loud is which annotated elements the opening
     interrupts. This is what catches opening a plasmid at base 1. */
  function openingReport(rec, points, removed) {
    if (!rec || !rec.features || !rec.features.length) return '';
    const L = rec.seq.length;
    const hit = [];
    rec.features.forEach(f => {
      if (f.type === 'source' || f.between) return;
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
    if (hf > 1 || hr > 1) msgs.push('a primer anneals to the plasmid in more than one place (' + hf + ' / ' + hr + ') — expect a second product');
    if (!msgs.length) return '';
    return '<div class="note warn"><strong>Backbone primers:</strong> ' + msgs.join('; ') + '.</div>';
  }
