/* ===========================================================================
   SEQUENCING FROM AN EXISTING PRIMER STOCK

   A lab accumulates hundreds of oligos. Most of them were ordered for
   something else, but a good number of them happen to sit in the right
   place to read a region you now care about. This file answers two
   questions: which tubes in the freezer can sequence this region, and
   what is the smallest set of them that reads it on both strands.

   The geometry of a Sanger read, measured from the primer's 3' end,
   because that is where the polymerase starts:

       primer      dead zone          confident read        ragged end
     ===========>..............-------------------------············
                 ^ first bases are unreadable (dye front, poor
                   resolution): the first 25-35 bases are expected to
                   be poor, so nothing inside that window is trusted
                                  ^ good data out to ~700 bases
                                                          ^ 800-900 is
                   reachable on clean template; beyond ~900, read with
                   caution

   Sources for the defaults are quoted in SANGER below. They are only
   defaults: every facility is different, and all three numbers are
   settable, because a plan is worthless if its assumptions do not match
   the service that will run the reaction.
   ========================================================================= */
(function (root) {
  'use strict';
  const E = (typeof require !== 'undefined' && typeof module !== 'undefined')
    ? require('./engine.js') : root.Bio;
  if (typeof require !== 'undefined' && typeof module !== 'undefined') require('./engine2.js');
  /* engine2 adds to the same object, so requiring it here is enough to be
     sure primerSites exists when this file is used from Node */
  
  /* ---------------------------------------------------------------------
     What a Sanger read looks like, in bases from the primer's 3' end.

     dead  the first bases after the primer that are not trustworthy.
           "The quality of the first 25 to 35 bases of the sequence is
           expected to be poor" (Michigan State RTSF Sanger guide, 2024),
           so 35 is the cautious edge of the published range.
     good  how far the confident read reaches. Facilities quote 700-900
           high-quality bases on clean plasmid.
     max   the furthest anything should be believed at all: "sequence
           lengths of 800 to 900 bases are possible for high-quality
           templates", and beyond about 900 bases the same guide says to
           interpret with caution. Some services quote up to 1100.

     A sequencing primer itself: 18-24 nt, Tm 50-60 C, GC 45-55% is the
     range the big providers ask for, and is what gradePrimer scores
     against.
  --------------------------------------------------------------------- */
  E.SANGER = { dead: 35, good: 700, max: 900, minLen: 18, maxLen: 24,
               tmLo: 50, tmHi: 60, gcLo: 0.40, gcHi: 0.60 };

  const clean = s => E.cleanSeq(String(s || '')).replace(/U/g, 'T');

  /* =======================================================================
     READING A STOCK LIST

     Inventories arrive as whatever the lab keeps them in: a spreadsheet
     paste, a CSV, or text copied out of a printed PDF where one record
     wraps over two lines and the sequence is typed in groups of three.
     All of that has to come in without the person tidying it up first,
     so the parser works by recognising the *shape* of each field rather
     than by trusting a column order.
  ======================================================================= */

  /* a run of bases, allowing the spaces and hyphens people type for
     legibility: "TTG TAT TTC CAG GGC" is one sequence, not five */
  const SEQ_RE = /(?:^|[\s,;|"'])((?:[ACGTacgtNn]{1,}[ \-]?){6,})(?=$|[\s,;|"'])/;
  function seqFrom(text) {
    let best = '';
    const t = ' ' + String(text || '') + ' ';
    /* Every candidate run, longest wins: a name like "GB1" must not beat
       the real sequence, and a stray "TE" must never win at all. Only a
       single space or hyphen may sit between bases, because people type
       "TTG TAT TTC" but columns are separated by wider gaps - allowing a
       run of spaces would let the tail of the column before bleed into
       the sequence and silently corrupt it. */
    const re = /[ACGTNacgtn](?:[ \-]?[ACGTNacgtn])*/g;
    let m;
    while ((m = re.exec(t))) {
      let groups = m[0].split(/[ \-]/).filter(Boolean);
      /* a lone base or two before the real thing is the last word of the
         previous column ("...N-BirA  ctcgag..."), not part of the oligo */
      while (groups.length > 1 && groups[0].length < 3 && groups[1].length >= 8) groups.shift();
      const c = groups.join('');
      if (c.length > best.length) best = c;
    }
    return best.toUpperCase();
  }

  const isSeqish = s => s && s.length >= 12 && /^[ACGTN]+$/.test(s);

  /* What a primer is called, with the order code and the date that often
     sit in the same cell taken off. These sheets write the two run
     together - "2019-09-25KPRB:S425F" - and what follows the date is the
     name, so that is preferred over what comes before it. */
  function tidyName(s) {
    let t = String(s || '').replace(/\s+/g, ' ').trim();
    const d = t.match(/\d{4}[-/.]\d{1,2}[-/.]\d{1,2}/);
    if (d) {
      const after = t.slice(t.indexOf(d[0]) + d[0].length).trim();
      if (/[A-Za-z]{2}/.test(after)) t = after;
      else t = t.slice(0, t.indexOf(d[0])).trim();
    }
    t = t.replace(/^[A-Z]{2,4}[-_]\d{1,4}\s+/, '').trim();
    return t;
  }

  /* a header row, if there is one: which column holds what */
  /** @type {Array<[string, RegExp]>} column name patterns, first match wins */
  const HEAD = [
    ['seq', /^(seq|sequence|oligo|oligo ?seq|primer ?seq|5'?-?>?-?3'?|bases)$/i],
    ['name', /^(name|primer|primer ?name|description|desc|oligo ?name|label|id ?name)$/i],
    ['id', /^(id|oligo ?id|primer ?id|stock ?id|ol|catalog|cat|#|no\.?|number)$/i],
    ['tm', /^(tm|tm ?\(?c?\)?|melting|melting ?temp\w*)$/i],
    ['conc', /^(conc|conc\.?|concentration|conc ?\(?u?m?\)?|um|stock ?conc\w*)$/i],
    ['well', /^(well|position|pos|slot)$/i],
    ['plate', /^(plate|box|rack|tube ?box)$/i],
    ['shelf', /^(shelf|freezer|rm|room|location|loc|where|where ?it ?is|stored ?(at|in)?)$/i],
    ['note', /^(note|notes|comment|comments|remark\w*)$/i],
    ['date', /^(date|ordered|order ?date|arrived)$/i],
    ['method', /^(clone ?method|method|use|purpose|application)$/i]
  ];

  function splitRow(line) {
    if (line.indexOf('\t') >= 0) return line.split('\t');
    /* a comma split is only right when it does not cut a sequence apart;
       two or more spaces is how aligned text and PDF text separate fields */
    if (/,/.test(line) && !/,\s*[ACGTacgt]{1,4}\s*,/.test(line)) {
      const cells = [];
      let cur = '', q = false;
      for (let i = 0; i < line.length; i++) {
        const c = line[i];
        if (c === '"') { q = !q; continue; }
        if (c === ',' && !q) { cells.push(cur); cur = ''; continue; }
        cur += c;
      }
      cells.push(cur);
      return cells;
    }
    return line.split(/\s{2,}/);
  }

  /* Read an inventory into primer records.
     Returns { primers, skipped, columns, mode } where each primer is
     { name, seq, id, tm, conc, loc, note } and `loc` is the human string
     that says which tube to pull. */
  E.parseStockTable = function (text, opts) {
    opts = opts || {};
    const raw = String(text || '').replace(/\r\n?/g, '\n').split('\n');
    const lines = raw.filter(l => l.trim().length > 0);
    const out = [], skipped = [];
    let columns = null, mode = 'shape';

    /* --- is there a usable header? --- */
    for (let i = 0; i < Math.min(5, lines.length); i++) {
      const cells = splitRow(lines[i]).map(c => c.trim());
      if (cells.length < 2) continue;
      const map = {};
      let hits = 0;
      cells.forEach((c, j) => {
        const key = c.replace(/\s+/g, ' ').trim();
        /* Real inventories label the column "Sequence 5'->3'", "Tm (C)",
           "Conc. (uM)". The word is the header; what follows it is units, an
           arrow or a prime mark. Strip a trailing run that carries no letters
           and try again, which recognises those without letting "Length (nt)"
           or "GC (%)" pass as a sequence column. */
        const bare = key.replace(/[^A-Za-z]+$/, '');
        for (const [field, re] of HEAD) {
          if ((re.test(key) || (bare && bare !== key && re.test(bare))) && map[field] == null) {
            map[field] = j; hits++; break;
          }
        }
      });
      /* a header is only believable if it names a sequence column plus
         something to call the primer, and the row itself is not data */
      const looksLikeData = cells.some(c => isSeqish(seqFrom(c)));
      if (map.seq != null && (map.name != null || map.id != null) && hits >= 2 && !looksLikeData) {
        /* and it has to actually work: text pulled out of a PDF has a
           header whose columns no longer line up with the rows under it,
           which would quietly drop most of the list. Try it on the rows
           that follow and keep it only if the sequence column really
           holds sequences. */
        const probe = lines.slice(i + 1, i + 14);
        let hit = 0, seen = 0;
        probe.forEach(l => {
          const c = splitRow(l);
          if (c.length < 2) return;
          seen++;
          if (isSeqish(seqFrom(c[map.seq] || ''))) hit++;
        });
        /* The sequence column lining up is not enough. A header pulled out
           of a PDF can have run-together labels — "ShelfPlate",
           "PositionOrderID" — so that the count of columns still matches
           while their meanings are shifted one across. The sequence may
           then land in the right place by luck while the freezer location
           reads from the room-number column: the primer is found and the
           person is sent to the wrong box, which is worse than not being
           told where it is. A column with a type has to hold that type. */
        const typed = (field, test) => {
          if (map[field] == null) return true;
          let good = 0, n = 0;
          probe.forEach(l => {
            const c = splitRow(l);
            if (c.length < 2) return;
            const v = (c[map[field]] || '').trim();
            if (!v) return;
            n++;
            if (test(v)) good++;
          });
          return n < 2 || good / n >= 0.6;
        };
        const number = v => /^[\d.]+$/.test(v.replace(/\s*[^\d.]*$/, '')) && isFinite(parseFloat(v));
        /* And the header must have as many columns as its rows do. Two of
           this header's labels have run together — "PositionOrderID Date"
           is one cell where the rows have two — so it describes eight
           columns of a nine-column table and everything after the join is
           named by the wrong label. Counting is the cheapest way to see
           that, and it catches the shift whether or not the misplaced
           columns happen to hold a type that can be probed. */
        const widths = {};
        probe.forEach(l => {
          const n = splitRow(l).length;
          if (n >= 2) widths[n] = (widths[n] || 0) + 1;
        });
        const modal = Object.keys(widths).sort((a, b) => widths[b] - widths[a])[0];
        const sameWidth = modal == null || +modal === cells.length;
        const aligned = sameWidth &&
                        typed('tm', v => number(v) && parseFloat(v) >= 20 && parseFloat(v) <= 95) &&
                        typed('conc', number) &&
                        typed('date', v => /^\d{4}-\d{2}-\d{2}$|^\d{1,2}[/.-]\d{1,2}[/.-]\d{2,4}$|^\d{1,2} \w{3,} \d{4}$/.test(v));
        if (seen >= 2 && hit / seen >= 0.6 && aligned) {
          columns = map; mode = 'header';
          lines.splice(0, i + 1);
        }
        break;
      }
    }

    const push = (rec, srcLine) => {
      if (!isSeqish(rec.seq)) { skipped.push({ line: srcLine, why: rec.seq ? 'sequence too short' : 'no sequence found' }); return; }
      if (!rec.name) rec.name = rec.id || ('primer ' + (out.length + 1));
      out.push(rec);
    };

    if (mode === 'header') {
      lines.forEach(line => {
        const cells = splitRow(line);
        const g = f => (columns[f] != null && cells[columns[f]] != null) ? String(cells[columns[f]]).trim() : '';
        const seq = seqFrom(g('seq'));
        const num = v => { const m = String(v).match(/-?\d+(\.\d+)?/); return m ? parseFloat(m[0]) : null; };
        const loc = [g('shelf'), g('plate'), g('well')].filter(Boolean).join(' / ');
        push({ name: tidyName(g('name')) || g('id'), seq: seq, id: g('id'), tm: num(g('tm')),
               conc: num(g('conc')), loc: loc, note: g('note'), stock: true,
               date: g('date'), method: g('method') }, line);
      });
      return { primers: out, skipped: skipped, columns: columns, mode: mode };
    }

    /* --- no header: go by shape, and allow a record to wrap over the
           lines that follow it, which is what copying out of a PDF gives --- */
    const idRe = /(^|\s)([A-Za-z]{1,6}[-_ ]?\d{2,5}[A-Za-z]?)(\s|$)/;
    const groups = [];
    lines.forEach(line => {
      const hasId = idRe.test(line.slice(0, 24));
      if (hasId || !groups.length) groups.push([line]);
      else groups[groups.length - 1].push(line);
    });

    groups.forEach(gl => {
      const blob = gl.join(' ');
      const seq = seqFrom(blob);
      if (!isSeqish(seq)) { skipped.push({ line: gl[0], why: 'no sequence found' }); return; }
      const idm = gl[0].match(idRe);
      const id = idm ? idm[2].trim() : '';
      /* the name is the most word-like field left once the id, the dates,
         the sequence and the bare numbers are taken out */
      let rest = blob;
      if (id) rest = rest.replace(id, ' ');
      const seqSpaced = blob.match(new RegExp(seq.split('').join('[ \\-]?'), 'i'));
      /* What the primer is called sits before its sequence, and the
         solution and the notes come after: every inventory reads
         description-then-sequence, so the text ahead of the sequence is
         where the name is, and "TE pH 8.0" never wins. */
      let ahead = rest, behind = '';
      if (seqSpaced) {
        const at = rest.indexOf(seqSpaced[0]);
        if (at >= 0) { ahead = rest.slice(0, at); behind = rest.slice(at + seqSpaced[0].length); }
        rest = rest.replace(seqSpaced[0], ' ');
      }
      const date = (rest.match(/\d{4}[-/.]\d{1,2}[-/.]\d{1,2}/) || [''])[0];
      /* What follows the date, in the same cell, is the description: these
         sheets read "2019-09-25KPRB:S425F" with the date run together
         with the name, and that name beats anything else in the row. */
      let dateName = '';
      if (date) {
        const after = (ahead.split(date)[1] || '').split(/\s{2,}/)[0] || '';
        if (/[A-Za-z]{2}/.test(after)) dateName = after.trim();
        rest = rest.replace(date, ' '); ahead = ahead.replace(date, ' ');
      }
      /* where the tube is: a well like A01 or B9, the box or plate it sits
         in, and a shelf or rack code. Worth keeping, because a primer you
         cannot find in the freezer is no use. */
      const well = (rest.match(/(^|\s)([A-H](?:0?[1-9]|1[0-2]))(\s|$)/) || [, , ''])[2] || '';
      const plate = (rest.match(/(^|\s)([A-Z]{2,}[A-Z0-9]*\d+)(\s|$)/) || [, , ''])[2] || '';
      /* a shelf is a letter and a digit or two, but the well looks the
         same, so the well is taken out of the running first */
      const forShelf = well ? rest.replace(new RegExp('(^|\\s)' + well + '(\\s|$)'), ' ') : rest;
      const shelf = (forShelf.match(/(^|\s)([A-Z]\d{1,2})(\s|$)/) || [, , ''])[2] || '';
      const tmm = rest.match(/(^|\s)([3-8]\d)(\.\d)?(\s|$)/);
      const tm = tmm ? parseFloat(tmm[2] + (tmm[3] || '')) : null;
      const nameFrom = txt => {
        let n = '';
        splitRow(txt).map(c => c.trim()).filter(Boolean).forEach(c => {
          if (/^\d+(\.\d+)?$/.test(c)) return;                   /* a bare number */
          if (/^[A-H](0?[1-9]|1[0-2])$/.test(c)) return;         /* a well */
          if (/^[A-Z]\d{1,2}$/.test(c)) return;                  /* a shelf */
          if (/^[A-Z]{2,}[A-Z0-9]*\d{2,}$/.test(c) && !/[a-z]/.test(c)) return;  /* a box code */
          if (/^(TE|pH|uM|nmol|H2O|water)\b/i.test(c) && c.length < 14) return;
          const letters = (c.match(/[A-Za-z]/g) || []).length;
          if (letters >= 2 && c.length > n.length) n = c;
        });
        return n;
      };
      /* ahead of the sequence first, and only then anywhere in the row */
      let name = tidyName(dateName || nameFrom(ahead) || nameFrom(rest) || nameFrom(behind));
      name = name.replace(/\s{2,}/g, ' ').trim();
      /* strip a leading date or order code that rode along with the name */
      name = name.replace(/^\s*(\d{4}[-/.]\d{1,2}[-/.]\d{1,2})\s*/, '')
                 .replace(/^\s*[A-Z]{2,4}-\d{2,4}\s*/, '').trim();
      const loc = [shelf, plate, well].filter(Boolean).join(' / ');
      push({ name: name || id, seq: seq, id: id, tm: tm, conc: null, loc: loc,
             note: '', stock: true, date: date, method: '' }, gl[0]);
    });
    return { primers: out, skipped: skipped, columns: null, mode: mode };
  };

  /* =======================================================================
     WHAT ONE PRIMER CAN READ

     All coordinates are 0-based half-open on the top strand, and on a
     circle a span may run across base 1, exactly as a feature does.
  ======================================================================= */

  /* the window a primer reads, as spans on the top strand */
  E.sangerWindow = function (site, L, opts) {
    opts = opts || {};
    const dead = opts.dead == null ? E.SANGER.dead : opts.dead;
    const good = opts.good == null ? E.SANGER.good : opts.good;
    const max = opts.max == null ? E.SANGER.max : opts.max;
    const circ = !!opts.circular;
    const tp = site.threePrime;
    let gs, ge, ms, me;
    if (site.strand === 1) {
      /* synthesis runs toward higher coordinates from the base after the
         primer's last annealed base */
      gs = tp + 1 + dead; ge = tp + 1 + good;
      ms = tp + 1 + dead; me = tp + 1 + max;
    } else {
      /* and toward lower coordinates for a primer on the bottom strand */
      ge = tp - dead; gs = tp - good;
      me = tp - dead; ms = tp - max;
    }
    const fix = (a, b) => {
      if (circ) {
        const len = Math.min(b - a, L);
        let s = ((a % L) + L) % L;
        return { start: s, end: s + len > L ? s + len - L : s + len, wraps: s + len > L, length: len };
      }
      const s = Math.max(0, a), e = Math.min(L, b);
      return { start: s, end: Math.max(s, e), wraps: false, length: Math.max(0, e - s) };
    };
    return { good: fix(gs, ge), reach: fix(ms, me), dead: fix(site.strand === 1 ? tp + 1 : tp - dead,
             site.strand === 1 ? tp + 1 + dead : tp) };
  };

  /* how much of [a,b) a span covers, on a line or around a circle */
  function overlapLen(span, a, b, L, circ) {
    if (!circ) return Math.max(0, Math.min(span.end, b) - Math.max(span.start, a));
    /* unroll both onto [0, 2L) starting at a */
    const rel = p => ((p - a) % L + L) % L;
    const want = ((b - a) % L + L) % L || L;
    const s = rel(span.start);
    const e = s + span.length;
    let n = 0;
    /* count the bases of the span that land in [0, want) */
    const lo = Math.max(0, s), hi = Math.min(want, e);
    if (hi > lo) n += hi - lo;
    if (e > L) { const lo2 = 0, hi2 = Math.min(want, e - L); if (hi2 > lo2) n += hi2 - lo2; }
    return Math.min(n, want);
  }
  E.spanOverlap = overlapLen;

  /* Score a primer as a *sequencing* primer. Returns a penalty in degrees
     and the reasons, so a caller can rank and also explain. */
  E.gradeSeqPrimer = function (annealSeq, opts) {
    opts = opts || {};
    const s = clean(annealSeq);
    const S = E.SANGER;
    const tm = E.tmNN(s, opts.conditions || E.PCR_CONDITIONS);
    const gc = s.length ? (s.match(/[GC]/g) || []).length / s.length : 0;
    const notes = [];
    let pen = 0;
    if (tm < S.tmLo) { pen += (S.tmLo - tm); notes.push('Tm ' + tm.toFixed(1) + ' C is below the ' + S.tmLo + '-' + S.tmHi + ' C range sequencing services ask for'); }
    if (tm > S.tmHi) { pen += (tm - S.tmHi) * 0.5; notes.push('Tm ' + tm.toFixed(1) + ' C is above the usual ' + S.tmLo + '-' + S.tmHi + ' C range'); }
    if (s.length < S.minLen) { pen += (S.minLen - s.length) * 1.5; notes.push('only ' + s.length + ' bases anneal; 18-24 is the usual ask'); }
    if (gc < S.gcLo) { pen += (S.gcLo - gc) * 40; notes.push('GC ' + Math.round(gc * 100) + '% is low'); }
    if (gc > S.gcHi) { pen += (gc - S.gcHi) * 40; notes.push('GC ' + Math.round(gc * 100) + '% is high'); }
    const run = (s.match(/A{5,}|C{5,}|G{5,}|T{5,}/g) || [])[0];
    if (run) { pen += run.length; notes.push('a run of ' + run.length + ' ' + run[0] + 's inside the primer'); }
    if (/[GC]{5,}$/.test(s)) { pen += 2; notes.push('a GC-rich 3′ end can mis-prime'); }
    return { tm: tm, gc: gc, length: s.length, penalty: pen, notes: notes };
  };

  /* =======================================================================
     THE PLAN

     sangerPlan(template, region, primers, opts) ->
       { candidates, plan, coverage, gaps, both, unusable }

     `candidates` is every stock primer that reads some part of the region,
     best first. `plan` is the smallest set that covers the region on each
     strand. `gaps` says what no tube in the freezer reaches.
  ======================================================================= */
  E.sangerPlan = function (template, region, primers, opts) {
    opts = opts || {};
    const T = clean(template), L = T.length;
    const circ = !!opts.circular;
    const dead = opts.dead == null ? E.SANGER.dead : opts.dead;
    const good = opts.good == null ? E.SANGER.good : opts.good;
    const max = opts.max == null ? E.SANGER.max : opts.max;
    const cond = opts.conditions || E.PCR_CONDITIONS;
    const allowRepeat = !!opts.allowRepeat;
    const site0 = { minAnneal: opts.minAnneal || 15, seed: opts.seed || 8,
                    maxMismatchFrac: opts.maxMismatchFrac == null ? 0.10 : opts.maxMismatchFrac,
                    circular: circ };
    if (!L) return { candidates: [], plan: { forward: [], reverse: [] }, coverage: null, gaps: [], unusable: [] };

    /* the region, as a start and a length going clockwise */
    let rs = region.start, re = region.end;
    if (circ) { rs = ((rs % L) + L) % L; re = ((re % L) + L) % L; }
    let rlen = circ ? (((re - rs) % L) + L) % L : Math.max(0, re - rs);
    if (circ && rlen === 0) rlen = L;
    if (!rlen) return { candidates: [], plan: { forward: [], reverse: [] }, coverage: null, gaps: [], unusable: [] };

    const list = (primers || []).filter(p => p && p.seq && clean(p.seq).length >= 12);
    const cands = [], unusable = [];

    list.forEach(p => {
      let sites;
      try { sites = E.primerSites(T, p.seq, site0); } catch (e) { sites = []; }
      if (!sites.length) return;              /* not in this construct at all */
      sites.forEach(st => {
        const win = E.sangerWindow(st, L, { dead: dead, good: good, max: max, circular: circ });
        const covers = overlapLen(win.good, rs, rs + rlen, L, circ);
        const reaches = overlapLen(win.reach, rs, rs + rlen, L, circ);
        if (!reaches) return;                 /* reads somewhere else entirely */
        /* the annealed part is what primes, so that is what gets graded */
        const annealSeq = st.strand === 1
          ? clean(p.seq).slice(clean(p.seq).length - st.annealed)
          : clean(p.seq).slice(clean(p.seq).length - st.annealed);
        const grade = E.gradeSeqPrimer(annealSeq, { conditions: cond });
        /* How far the primer sits from the near edge of the region, in the
           direction it reads. A primer that anneals *inside* the region
           has no meaningful run-up to it - going forward from there, the
           region's first base is nearly a whole plasmid away - so that
           case is marked rather than measured, or the modular distance
           would look like a disqualifying 3,000-base gap. */
        const tpRel = ((st.threePrime - rs) % L + L) % L;
        const inside = tpRel < rlen;
        const toStart = inside ? 0 : (st.strand === 1
          ? ((rs - (st.threePrime + 1)) % L + L) % L
          : ((st.threePrime - (rs + rlen)) % L + L) % L);
        const rec = {
          primer: p, name: p.name, seq: p.seq, id: p.id || '', loc: p.loc || '',
          site: st, strand: st.strand, window: win,
          covers: covers, reach: reaches, fraction: covers / rlen,
          gap: toStart,                       /* bases between primer and region */
          inside: inside,                     /* the primer anneals within the region */
          annealed: st.annealed, tail: st.tail, mismatches: st.mismatches.slice(),
          tm: grade.tm, gc: grade.gc, grade: grade,
          sites: sites.length, warnings: []
        };
        if (sites.length > 1 && !allowRepeat)
          rec.warnings.push('this primer anneals in ' + sites.length + ' places on the template, which gives a mixed, unreadable trace');
        if (st.mismatches.length)
          rec.warnings.push(st.mismatches.length + ' mismatch' + (st.mismatches.length === 1 ? '' : 'es') +
            ' against this construct (primer base ' + st.mismatches.join(', ') + ')');
        if (covers < reaches)
          rec.warnings.push('part of the region falls past the confident ' + good + '-base read, out in the ' + good + '-' + max + ' tail');
        const deadIn = overlapLen(win.dead, rs, rs + rlen, L, circ);
        if (deadIn > 0)
          rec.warnings.push(inside
            ? 'this primer anneals inside the region, so the ' + deadIn + ' bases just after it are in the unreadable window and need reading from somewhere else'
            : 'the primer sits only ' + toStart + ' bases from the region, so its first ' + deadIn + ' bases fall in the unreadable window');
        grade.notes.forEach(n => rec.warnings.push(n));
        /* slippage: a long homopolymer between the primer and the region
           can throw the whole read out of register */
        const leadS = st.strand === 1 ? st.threePrime + 1 : rs + rlen;
        const leadN = inside ? 0 : Math.min(Math.max(0, toStart), 400);
        if (leadN > 0) {
          let lead = '';
          for (let i = 0; i < leadN; i++) lead += T[(((leadS + i) % L) + L) % L];
          const hp = (lead.match(/A{8,}|C{8,}|G{8,}|T{8,}/g) || [])[0];
          if (hp) rec.warnings.push('a run of ' + hp.length + ' ' + hp[0] + 's lies between the primer and the region, which can make the polymerase slip and put the rest of the read out of register');
        }
        if (sites.length > 1 && !allowRepeat) unusable.push(rec); else cands.push(rec);
      });
    });

    /* best first: how much of the region it reads, then how good a
       sequencing primer it is, then how comfortably it clears the dead zone */
    const ideal = dead + 60;
    cands.forEach(c => {
      c.score = c.covers
        - c.grade.penalty * 3
        - (c.inside ? 0 : Math.abs(c.gap - ideal) * 0.05)
        - c.mismatches.length * 25;
    });
    cands.sort((a, b) => b.score - a.score);

    /* --- the smallest set that covers the region on one strand --- */
    function cover(strandCands) {
      const picked = [];
      const done = new Uint8Array(rlen);
      let left = rlen, guard = 0;
      const relOf = c => {
        /* the covered run, as offsets into the region */
        const w = c.window.good;
        const startRel = ((w.start - rs) % L + L) % L;
        return { from: startRel, to: startRel + w.length };
      };
      while (left > 0 && guard++ < 200) {
        let best = null, bestGain = 0;
        /* the first base still missing */
        let first = -1;
        for (let i = 0; i < rlen; i++) if (!done[i]) { first = i; break; }
        if (first < 0) break;
        strandCands.forEach(c => {
          if (picked.indexOf(c) >= 0) return;
          const r = relOf(c);
          let gain = 0, touchesFirst = false;
          for (let i = Math.max(0, r.from); i < Math.min(rlen, r.to); i++) {
            if (!done[i]) { gain++; if (i === first) touchesFirst = true; }
          }
          /* wrapped window: the part that came round the circle */
          if (r.to > L) for (let i = 0; i < Math.min(rlen, r.to - L); i++) {
            if (!done[i]) { gain++; if (i === first) touchesFirst = true; }
          }
          if (!gain) return;
          /* prefer one that starts at the hole and runs furthest, and
             among those the better primer */
          const w = (touchesFirst ? 1e6 : 0) + gain * 10 - c.grade.penalty;
          if (w > bestGain) { bestGain = w; best = c; }
        });
        if (!best) break;
        const r = relOf(best);
        for (let i = Math.max(0, r.from); i < Math.min(rlen, r.to); i++) done[i] = 1;
        if (r.to > L) for (let i = 0; i < Math.min(rlen, r.to - L); i++) done[i] = 1;
        picked.push(best);
        left = 0; for (let i = 0; i < rlen; i++) if (!done[i]) left++;
      }
      /* what is still missing, as spans in region offsets */
      const holes = [];
      let i = 0;
      while (i < rlen) {
        if (done[i]) { i++; continue; }
        let j = i; while (j < rlen && !done[j]) j++;
        holes.push({ from: i, to: j, length: j - i,
                     start: (rs + i) % L, end: (rs + j) % L });
        i = j;
      }
      return { picked: picked, covered: rlen - left, holes: holes, mask: done };
    }

    const fwd = cover(cands.filter(c => c.strand === 1));
    const rev = cover(cands.filter(c => c.strand === -1));

    /* bases read on both strands: the standard for calling a sequence
       confirmed rather than merely seen */
    let both = 0;
    for (let i = 0; i < rlen; i++) if (fwd.mask[i] && rev.mask[i]) both++;
    let either = 0;
    for (let i = 0; i < rlen; i++) if (fwd.mask[i] || rev.mask[i]) either++;

    return {
      region: { start: rs, end: (rs + rlen) % L || (circ ? rs : rlen), length: rlen },
      settings: { dead: dead, good: good, max: max },
      candidates: cands,
      unusable: unusable,
      plan: { forward: fwd.picked, reverse: rev.picked },
      coverage: {
        forward: fwd.covered, reverse: rev.covered,
        both: both, either: either, length: rlen,
        fractionBoth: both / rlen, fractionEither: either / rlen
      },
      gaps: { forward: fwd.holes, reverse: rev.holes }
    };
  };

  /* =======================================================================
     WHEN THE FREEZER CANNOT DO IT

     Design a sequencing primer that reads a span the stock does not
     reach: scan the window upstream of it for the best unique primer.
  ======================================================================= */
  E.designSeqPrimer = function (template, target, strand, opts) {
    opts = opts || {};
    const T = clean(template), L = T.length;
    const circ = !!opts.circular;
    const dead = opts.dead == null ? E.SANGER.dead : opts.dead;
    const cond = opts.conditions || E.PCR_CONDITIONS;
    const S = E.SANGER;
    /* put the 3' end this far ahead of the target: clear of the dead zone
       with room to spare, so the first bases of interest are clean */
    const want = opts.lead == null ? dead + 60 : opts.lead;
    const at = i => T[(((i % L) + L) % L)];
    const slice = (a, n) => { let s = ''; for (let k = 0; k < n; k++) s += at(a + k); return s; };
    const best = [];
    /* try 3' ends from a little before the ideal lead to a good way back */
    for (let off = -40; off <= 240; off++) {
      const tp = strand === 1 ? target - want - off : target + want + off;
      if (!circ && (tp < 0 || tp >= L)) continue;
      for (let len = S.minLen; len <= S.maxLen + 2; len++) {
        let seq;
        if (strand === 1) {
          const a = tp - len + 1;
          if (!circ && a < 0) continue;
          seq = slice(a, len);
        } else {
          if (!circ && tp + len > L) continue;
          seq = E.revcomp(slice(tp, len));
        }
        if (!/^[ACGT]+$/.test(seq)) continue;
        const g = E.gradeSeqPrimer(seq, { conditions: cond });
        let n;
        try { n = E.primerSites(T, seq, { circular: circ }).length; } catch (e) { n = 2; }
        if (n !== 1) continue;                      /* must be the only site */
        const reach = strand === 1 ? ((target - (tp + 1)) % L + L) % L
                                   : ((tp - target) % L + L) % L;
        if (reach < dead) continue;                 /* target would be unreadable */
        best.push({ seq: seq, strand: strand, threePrime: ((tp % L) + L) % L,
                    tm: g.tm, gc: g.gc, length: seq.length, lead: reach,
                    penalty: g.penalty + Math.abs(reach - want) * 0.03, notes: g.notes });
      }
    }
    best.sort((a, b) => a.penalty - b.penalty);
    return best.length ? best[0] : null;
  };

  /* =======================================================================
     SEQUENCING A PCR PRODUCT INSTEAD

     Pairs in the stock that amplify the region, for when the plan is to
     sequence an amplicon rather than the plasmid.
  ======================================================================= */
  E.sangerAmplicons = function (template, region, primers, opts) {
    opts = opts || {};
    const T = clean(template), L = T.length;
    const circ = !!opts.circular;
    const minProd = opts.minProduct || 150, maxProd = opts.maxProduct || 5000;
    const margin = opts.margin == null ? (opts.dead == null ? E.SANGER.dead : opts.dead) + 20 : opts.margin;
    const cond = opts.conditions || E.PCR_CONDITIONS;
    let rs = region.start, re = region.end;
    if (circ) { rs = ((rs % L) + L) % L; re = ((re % L) + L) % L; }
    let rlen = circ ? (((re - rs) % L) + L) % L : Math.max(0, re - rs);
    if (circ && rlen === 0) rlen = L;

    const F = [], R = [];
    (primers || []).forEach(p => {
      if (!p || !p.seq || clean(p.seq).length < 12) return;
      let sites; try { sites = E.primerSites(T, p.seq, { circular: circ }); } catch (e) { return; }
      if (sites.length !== 1) return;              /* one site only, both ways */
      const st = sites[0];
      const annealSeq = clean(p.seq).slice(clean(p.seq).length - st.annealed);
      const entry = { primer: p, name: p.name, id: p.id || '', loc: p.loc || '', site: st,
                      tm: E.tmNN(annealSeq, cond), annealed: st.annealed, tail: st.tail,
                      mismatches: st.mismatches.slice() };
      (st.strand === 1 ? F : R).push(entry);
    });

    const out = [];
    F.forEach(f => R.forEach(r => {
      /* product runs from the forward primer's 5'-most annealed base to
         the reverse primer's 5'-most annealed base, clockwise */
      const a = f.site.strand === 1 ? f.site.start : f.site.end;
      const b = r.site.end;
      const span = circ ? (((b - a) % L) + L) % L : b - a;
      /* The product is longer than the stretch of template between the two
         annealed footprints: a 5' tail does not anneal, but it is copied
         into the product and shows up on the gel. Reporting the template
         span as the product length would disagree with the PCR itself for
         any primer carrying an arm, which most cloning primers do. */
      const len = span + f.tail + r.tail;
      if (span < minProd || span > maxProd) return;
      /* the region must sit inside, clear of both ends */
      const relR = ((rs - a) % L + L) % L;
      if (relR < margin) return;
      if (relR + rlen > span - margin) return;
      const dTm = Math.abs(f.tm - r.tm);
      out.push({
        forward: f, reverse: r, length: len, span: span, dTm: dTm,
        start: a, end: b,
        inset: { fromStart: relR, fromEnd: span - (relR + rlen) },
        mismatches: f.mismatches.length + r.mismatches.length,
        score: -dTm * 2 - Math.abs(span - (rlen + 300)) * 0.01 - (f.mismatches.length + r.mismatches.length) * 10,
        warnings: [].concat(
          dTm > 5 ? ['the two Tms differ by ' + dTm.toFixed(1) + ' C, so one anneal temperature will not suit both'] : [],
          f.mismatches.length || r.mismatches.length
            ? ['mismatches against this construct: ' + (f.mismatches.length + r.mismatches.length)] : [],
          f.tail > 8 || r.tail > 8
            ? ['one of these carries a ' + Math.max(f.tail, r.tail) + '-base 5′ tail that will end up in the product'] : []
        )
      });
    }));
    out.sort((a, b) => b.score - a.score);
    return out;
  };

  /* =======================================================================
     REPLACING ONE STRETCH WITH ANOTHER, INSIDE A CODING SEQUENCE

     Swapping a solubility tag is not the same kind of edit as swapping a
     vector. The new piece lands in the middle of a reading frame, and
     everything behind it - the protease site, the protein itself - is only
     still protein if the join lands on a codon boundary and the length
     change is a whole number of codons. One base out and the rest of the
     gene is read in the wrong frame: the protease site is gone, the
     protein is nonsense, and nothing about the plasmid map looks wrong.

     So rather than trusting arithmetic, this translates the coding
     sequence before and after the edit and compares them.
  ======================================================================= */

  /* sites people put between a tag and a protein, as they read in protein */
  E.PROTEASE_SITES = [
    { name: 'TEV', motif: 'ENLYFQ', note: 'cleaves after Q, leaving a G or S' },
    { name: 'PreScission / HRV 3C', motif: 'LEVLFQ', note: 'cleaves after Q, leaving GP' },
    { name: 'Thrombin', motif: 'LVPR', note: 'cleaves after R' },
    { name: 'Enterokinase', motif: 'DDDDK', note: 'cleaves after K' },
    { name: 'Factor Xa', motif: 'IEGR', note: 'cleaves after R' }
  ];

  /* Replace [from, to) with `insert` and say what it does to the protein.
     opts: { circular, features, code, cds }  */
  E.swapFrameCheck = function (template, from, to, insert, opts) {
    opts = opts || {};
    const circ = !!opts.circular;
    const T = clean(template), L = T.length;
    const ins = clean(insert || '');
    if (!L) return { applies: false, why: 'There is no sequence.' };
    const a = ((from % L) + L) % L, b = ((to % L) + L) % L;
    const removed = circ ? (((b - a) % L) + L) % L : Math.max(0, b - a);
    const delta = ins.length - removed;

    /* the edited sequence, plus a way to carry a coordinate across the edit */
    let neu, map;
    if (!circ || a <= b) {
      neu = T.slice(0, a) + ins + T.slice(b);
      map = p => (p <= a ? p : p >= b ? p + delta : null);
    } else {
      /* the removed span runs across the origin: rotate so it sits at the
         end, which makes the edit a plain truncation plus the insert */
      const rot = T.slice(b) + T.slice(0, b);
      neu = rot.slice(0, L - removed) + ins;
      map = p => {
        const q = ((p - b) % L + L) % L;
        return q <= L - removed ? q : null;
      };
    }

    /* --- what reading frame is this edit sitting in? --- */
    const inSpan = (s, e) => {               /* does [s,e) contain the whole edit */
      if (!circ) return s <= a && e >= b;
      const rel = p => ((p - s) % L + L) % L;
      const len = ((e - s) % L + L) % L || L;
      return rel(a) < len && rel(a) + removed <= len;
    };
    let cds = opts.cds || null;
    if (!cds) {
      const feats = (opts.features || []).filter(f =>
        /^(CDS|gene|ORF)$/i.test(f.type || '') && f.start != null && f.end != null);
      cds = feats.filter(f => inSpan(f.start, f.end))
        .sort((x, y) => (((y.end - y.start) % L + L) % L) - (((x.end - x.start) % L + L) % L))[0] || null;
    }
    if (!cds) {
      /* no annotation to go on: take the longest open reading frame that
         contains the edit, which is what a tag-and-protein cassette is */
      let best = null;
      try {
        const orfs = E.findORFs(T, { circular: circ, minLen: 150 }) || [];
        orfs.forEach(o => {
          if ((o.strand || 1) !== 1) return;
          if (!inSpan(o.start, o.end)) return;
          const len = ((o.end - o.start) % L + L) % L;
          if (!best || len > (((best.end - best.start) % L + L) % L)) best = o;
        });
      } catch (e) { best = null; }
      cds = best;
    }
    if (!cds) {
      return { applies: false, removed: removed, inserted: ins.length, delta: delta,
               why: 'Nothing coding was found around this edit, so there is no reading frame to keep. ' +
                    'If this sits inside a gene, annotate it as a CDS and run the check again.' };
    }
    if ((cds.strand || 1) === -1) {
      /* a tag on the bottom strand: say so rather than guess */
      return { applies: false, removed: removed, inserted: ins.length, delta: delta,
               cds: cds,
               why: 'The coding sequence here runs on the bottom strand; this check only reads forward frames.' };
    }

    /* --- translate it, before and after --- */
    const cdsLen = ((cds.end - cds.start) % L + L) % L || L;
    const grab = (seq, start, len) => {
      const n = seq.length;
      let out = '';
      for (let i = 0; i < len; i++) out += seq[((start + i) % n + n) % n];
      return out;
    };
    const oldCds = grab(T, cds.start, cdsLen);
    const newStart = map(cds.start);
    if (newStart == null) {
      return { applies: false, why: 'The edit removes the start of the coding sequence itself.' };
    }
    /* read on past the old end in case the stop has moved */
    const newCds = grab(neu, newStart, Math.min(neu.length, cdsLen + delta + 300));
    const code = opts.code;
    const P0 = E.translate(oldCds, code ? { code: code } : undefined);
    const P1 = E.translate(newCds, code ? { code: code } : undefined);

    /* where the edit sits, in codons from the start of the gene */
    const offset = ((a - cds.start) % L + L) % L;
    const startsOnCodon = offset % 3 === 0;
    const endsOnCodon = (offset + removed) % 3 === 0;
    const inFrame = delta % 3 === 0;

    /* the protein behind the join: the part that must survive */
    const tailOld = P0.slice(Math.ceil((offset + removed) / 3));
    const tailNew = P1.slice(Math.ceil((offset + ins.length) / 3));
    const tailKept = tailOld.length > 0 && tailNew.indexOf(tailOld.replace(/\*.*$/, '')) === 0;

    /* a stop where there was none */
    const bodyNew = P1.slice(0, Math.floor((cdsLen + delta) / 3));
    const stopAt = bodyNew.slice(0, -1).indexOf('*');

    /* a protease site that used to sit behind the tag */
    /** @type {any} */
    let protease = null;
    E.PROTEASE_SITES.forEach(ps => {
      if (protease) return;
      const near = P0.slice(Math.max(0, Math.floor(offset / 3) - 2),
                            Math.ceil((offset + removed) / 3) + 12);
      if (near.indexOf(ps.motif) < 0) return;
      const stillThere = P1.indexOf(ps.motif) >= 0;
      protease = { name: ps.name, motif: ps.motif, note: ps.note, kept: stillThere };
    });

    /** @type {any[]} */
    const notes = [];
    if (!startsOnCodon) notes.push('The edit starts ' + (offset % 3) + ' base' + (offset % 3 === 1 ? '' : 's') +
      ' into a codon, so the residue at the join is a hybrid of the old and the new.');
    if (!endsOnCodon) notes.push('The edit ends part-way through a codon.');
    if (!inFrame) notes.push('The piece going in is ' + Math.abs(delta) + ' base' + (Math.abs(delta) === 1 ? '' : 's') +
      ' ' + (delta > 0 ? 'longer' : 'shorter') + ' than the piece coming out, which is not a whole number of codons: ' +
      'everything downstream will be read in the wrong frame.');
    if (inFrame && !tailKept) notes.push('The protein behind the join does not read as it did before.');
    if (stopAt >= 0) notes.push('A stop codon now falls at residue ' + (stopAt + 1) +
      ', which truncates the product.');
    if (protease && !protease.kept) notes.push('The ' + protease.name + ' site (' + protease.motif +
      ') that sat behind the tag is gone.');

    const good = inFrame && tailKept && stopAt < 0 && (!protease || protease.kept);
    return {
      applies: true, ok: good,
      cds: { start: cds.start, end: cds.end, label: cds.label || cds.type || 'CDS', length: cdsLen },
      removed: removed, inserted: ins.length, delta: delta,
      offset: offset, codon: Math.floor(offset / 3) + 1,
      startsOnCodon: startsOnCodon, endsOnCodon: endsOnCodon,
      inFrame: inFrame, tailKept: tailKept,
      stopAt: stopAt >= 0 ? stopAt + 1 : null,
      protease: protease,
      proteinBefore: P0, proteinAfter: P1,
      /* A look at each junction rather than at the whole insert: what the
         protein reads where the new piece starts, and where it rejoins
         what was already there. Printing the entire tag in between would
         bury the only part worth reading. */
      join: (function () {
        const win = 8;
        const cut = (P, at) => ({
          left: P.slice(Math.max(0, at - win), at),
          right: P.slice(at, at + win)
        });
        const o5 = Math.floor(offset / 3);
        return {
          before: { start: cut(P0, o5), end: cut(P0, Math.ceil((offset + removed) / 3)) },
          after: { start: cut(P1, o5), end: cut(P1, Math.ceil((offset + ins.length) / 3)) },
          insertResidues: Math.floor(ins.length / 3)
        };
      })(),
      residuesBefore: Math.floor(cdsLen / 3), residuesAfter: Math.floor((cdsLen + delta) / 3),
      notes: notes
    };
  };

  if (typeof module !== 'undefined' && module.exports) module.exports = E;
  root.Bio = E;
})(typeof window !== 'undefined' ? window : globalThis);
