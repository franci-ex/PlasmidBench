
/* ============================================================================
   PlasmidBench — cloning operations, assembly and analysis layer.
   Extends the Bio namespace from engine.js.
   ============================================================================ */
(function (root) {
  'use strict';
  const E = root.Bio;

  /* ======================= END / LIGATION MODEL =======================
     A fragment end is { type: "5'"|"3'"|"blunt", seq: <overhang 5'->3'> }.
     Two ends ligate iff both blunt, or same type and one overhang is the
     reverse complement of the other.
  ==================================================================== */
  E.endOf = function (seq, cutTop, cutBottom, side, L) {
    // side: 'left' (end of a fragment that starts at the cut) or 'right'
    const d = cutBottom - cutTop;
    if (d === 0) return { type: 'blunt', seq: '' };
    const mod = i => ((i % L) + L) % L;
    if (d > 0) {
      // 5' overhang; top-strand bases cutTop..cutBottom-1 are single stranded
      let oh = '';
      for (let i = 0; i < d; i++) oh += seq[mod(cutTop + i)];
      return side === 'left' ? { type: "5'", seq: oh } : { type: "5'", seq: E.revcomp(oh) };
    } else {
      const n = -d;
      let oh = '';
      for (let i = 0; i < n; i++) oh += seq[mod(cutBottom + i)];
      // 3' overhang; bottom-strand bases are single stranded at the left end
      return side === 'left' ? { type: "3'", seq: E.revcomp(oh) } : { type: "3'", seq: oh };
    }
  };

  E.endsCompatible = function (a, b) {
    if (!a || !b) return false;
    if (a.type === 'blunt' && b.type === 'blunt') return true;
    if (a.type !== b.type) return false;
    if (!a.seq.length || a.seq.length !== b.seq.length) return false;
    return a.seq === E.revcomp(b.seq);
  };

  // Digest that also annotates each fragment with structured ends.
  E.digestWithEnds = function (seq, enzymeNames, opts) {
    opts = opts || {};
    seq = E.cleanSeq(seq);
    const L = seq.length;
    const d = E.digest(seq, enzymeNames, opts);
    const cutAt = {};
    for (const c of d.cuts) cutAt[c.top] = c;
    for (const f of d.fragments) {
      const cl = cutAt[f.start], cr = cutAt[f.end % L] || cutAt[f.end];
      f.left = cl ? E.endOf(seq, cl.top, cl.bottom, 'left', L) : { type: 'blunt', seq: '' };
      f.right = cr ? E.endOf(seq, cr.top, cr.bottom, 'right', L) : { type: 'blunt', seq: '' };
    }
    return d;
  };

  // Which enzymes produce ends compatible with a given end?
  E.compatibleEnzymes = function (end) {
    const out = [];
    for (const enz of E.ALL_ENZYMES) {
      if (enz.typeIIS) continue;
      const probe = 'ATCGATCG' + concrete(enz.site) + 'ATCGATCG';
      const sites = E.findSites(probe, enz, false);
      if (!sites.length) continue;
      const s = sites[0];
      const left = E.endOf(probe, s.cutTop, s.cutBottom, 'left', probe.length);
      if (E.endsCompatible(end, left)) out.push({ name: enz.name, site: enz.site, overhang: left });
    }
    return out;
  };
  function concrete(site) {
    const pick = { R:'A',Y:'C',S:'G',W:'A',K:'G',M:'A',B:'C',D:'A',H:'A',V:'A',N:'A' };
    return site.split('').map(c => 'ACGT'.includes(c) ? c : (pick[c] || 'A')).join('');
  }

  /* ======================= RESTRICTION CLONING =======================
     Digest vector and insert, pick the fragments the user wants, verify the
     two junctions, and build the circular product.
  ==================================================================== */
  E.restrictionClone = function (opts) {
    const vec = E.cleanSeq(opts.vector);
    const ins = E.cleanSeq(opts.insert);
    const vd = E.digestWithEnds(vec, opts.vectorEnzymes, { circular: opts.vectorCircular !== false, methylation: opts.vectorMethylation });
    const idg = E.digestWithEnds(ins, opts.insertEnzymes, { circular: !!opts.insertCircular, methylation: opts.insertMethylation });

    // choose vector backbone: largest fragment unless told otherwise
    let vfrag = opts.vectorFragmentIndex != null
      ? vd.fragments[opts.vectorFragmentIndex]
      : vd.fragments.slice().sort((a, b) => b.length - a.length)[0];
    let ifrag = opts.insertFragmentIndex != null
      ? idg.fragments[opts.insertFragmentIndex]
      : idg.fragments.slice().sort((a, b) => b.length - a.length)[0];
    if (!vfrag || !ifrag) return { ok: false, error: 'No fragments produced. Check that both enzymes cut.' };

    const tries = [
      { orientation: 'forward', seq: ifrag.seq, left: ifrag.left, right: ifrag.right },
      { orientation: 'reverse', seq: E.revcomp(ifrag.seq), left: flip(ifrag.right), right: flip(ifrag.left) }
    ];
    const results = [];
    for (const t of tries) {
      const j1 = E.endsCompatible(vfrag.right, t.left);
      const j2 = E.endsCompatible(t.right, vfrag.left);
      if (j1 && j2) {
        const product = vfrag.seq + t.seq;
        results.push({
          orientation: t.orientation,
          seq: product, length: product.length, circular: true,
          junction5: describeJunction(vfrag.right, t.left),
          junction3: describeJunction(t.right, vfrag.left)
        });
      }
    }
    if (!results.length) {
      return {
        ok: false,
        error: 'Ends are not compatible.',
        detail: {
          vectorLeft: vfrag.left, vectorRight: vfrag.right,
          insertLeft: ifrag.left, insertRight: ifrag.right
        },
        vectorFragment: vfrag, insertFragment: ifrag
      };
    }
    return { ok: true, products: results, vectorFragment: vfrag, insertFragment: ifrag };
  };
  function flip(end) {
    if (!end || end.type === 'blunt') return { type: 'blunt', seq: '' };
    return { type: end.type, seq: E.revcomp(end.seq) };
  }
  function describeJunction(a, b) {
    if (a.type === 'blunt') return 'blunt ligation';
    return a.type + ' ' + a.seq + ' / ' + b.seq + ' (' + a.seq.length + ' nt)';
  }

  /* ======================= GOLDEN GATE =======================
     Excise each part with a Type IIS enzyme, read the 4-nt overhangs, and
     order the parts by matching overhangs into a circle.
  ==================================================================== */
  E.goldenGate = function (parts, enzymeName, opts) {
    opts = opts || {};
    const enz = E.enzymeByName(enzymeName);
    if (!enz || !enz.typeIIS) return { ok: false, error: enzymeName + ' is not a Type IIS enzyme.' };
    const pieces = [];
    for (const p of parts) {
      const d = E.digestWithEnds(p.seq, [enzymeName], { circular: !!p.circular });
      // the released insert is the fragment flanked by two enzyme cuts
      const internal = d.fragments.filter(f => f.leftEnzyme === enzymeName && f.rightEnzyme === enzymeName);
      const chosen = internal.length
        ? internal.sort((a, b) => b.length - a.length)[0]
        : d.fragments.sort((a, b) => b.length - a.length)[0];
      if (!chosen) return { ok: false, error: 'No fragment released from ' + p.name };
      pieces.push({
        name: p.name, seq: chosen.seq, length: chosen.length,
        left: chosen.left, right: chosen.right, nCuts: d.cuts.length
      });
    }
    // overhang uniqueness check
    const seen = {}, clashes = [];
    for (const pc of pieces) {
      const k = pc.left.seq;
      if (k) { if (seen[k]) clashes.push(k + ' (' + seen[k] + ' and ' + pc.name + ')'); else seen[k] = pc.name; }
    }
    // order parts into a circle by matching right -> left
    const order = [];
    const used = new Set();
    let cur = pieces[0];
    order.push(cur); used.add(0);
    let guard = 0;
    while (order.length < pieces.length && guard++ < pieces.length * pieces.length) {
      let found = -1;
      for (let i = 0; i < pieces.length; i++) {
        if (used.has(i)) continue;
        if (E.endsCompatible(cur.right, pieces[i].left)) { found = i; break; }
      }
      if (found === -1) break;
      used.add(found); cur = pieces[found]; order.push(cur);
    }
    const closes = order.length === pieces.length &&
      E.endsCompatible(order[order.length - 1].right, order[0].left);
    const product = order.map(p => p.seq).join('');
    return {
      ok: order.length === pieces.length && closes,
      pieces, order: order.map(p => p.name),
      overhangClashes: clashes,
      closes,
      product: { seq: product, length: product.length, circular: true },
      error: order.length !== pieces.length
        ? 'Could not chain all parts: overhangs do not match up.'
        : (!closes ? 'Parts chain but the assembly does not close into a circle.' : null)
    };
  };

  /* ======================= GIBSON / HiFi / In-Fusion =======================
     Join fragments that share terminal homology. Finds the overlap between
     consecutive fragments and fuses them, removing the duplicate.
  ==================================================================== */
  E.findTerminalOverlap = function (a, b, minLen, maxLen) {
    minLen = minLen || 15; maxLen = maxLen || 60;
    const lim = Math.min(maxLen, a.length, b.length);
    for (let n = lim; n >= minLen; n--) {
      if (a.slice(-n) === b.slice(0, n)) return n;
    }
    return 0;
  };

  E.gibsonAssemble = function (fragments, opts) {
    opts = opts || {};
    const minOv = opts.minOverlap || 15;
    const maxOv = opts.maxOverlap || 60;
    const circular = opts.circular !== false;
    let frags = fragments.map(f => ({ name: f.name, seq: E.cleanSeq(f.seq) }));
    if (frags.length === 1) {
      const n = E.findTerminalOverlap(frags[0].seq, frags[0].seq, minOv, maxOv);
      if (circular && n) {
        return { ok: true, order: [frags[0].name], overlaps: [n],
                 product: { seq: frags[0].seq.slice(0, frags[0].seq.length - n), circular: true,
                            length: frags[0].seq.length - n } };
      }
      return { ok: false, error: 'Single fragment has no terminal self-overlap to circularise.' };
    }
    // greedy chain: start from fragment 0, repeatedly find the next with overlap
    const order = [frags[0]];
    const used = new Set([0]);
    const overlaps = [];
    let guard = 0;
    while (order.length < frags.length && guard++ < frags.length * frags.length) {
      const cur = order[order.length - 1];
      let best = -1, bestN = 0, rev = false;
      for (let i = 0; i < frags.length; i++) {
        if (used.has(i)) continue;
        const n = E.findTerminalOverlap(cur.seq, frags[i].seq, minOv, maxOv);
        if (n > bestN) { bestN = n; best = i; rev = false; }
        const rc = E.revcomp(frags[i].seq);
        const n2 = E.findTerminalOverlap(cur.seq, rc, minOv, maxOv);
        if (n2 > bestN) { bestN = n2; best = i; rev = true; }
      }
      if (best === -1) break;
      used.add(best);
      order.push({ name: frags[best].name + (rev ? ' (flipped)' : ''),
                   seq: rev ? E.revcomp(frags[best].seq) : frags[best].seq });
      overlaps.push(bestN);
    }
    if (order.length < frags.length) {
      return { ok: false, error: 'Fragments do not share enough terminal homology (need ' +
               minOv + ' bp or more). Chained ' + order.length + ' of ' + frags.length + '.',
               order: order.map(o => o.name), overlaps };
    }
    // fuse
    let product = order[0].seq;
    for (let i = 1; i < order.length; i++) {
      const n = overlaps[i - 1];
      product += order[i].seq.slice(n);
    }
    let closing = 0;
    if (circular) {
      closing = E.findTerminalOverlap(product, product, minOv, maxOv);
      if (!closing) {
        return { ok: false, error: 'Assembly chains but the last fragment does not overlap the first, so it cannot circularise.',
                 order: order.map(o => o.name), overlaps,
                 product: { seq: product, length: product.length, circular: false } };
      }
      product = product.slice(0, product.length - closing);
    }
    return {
      ok: true, order: order.map(o => o.name), overlaps, closingOverlap: closing,
      product: { seq: product, length: product.length, circular: circular }
    };
  };

  // Design the primer pairs for a Gibson/HiFi/In-Fusion assembly of N inserts
  // into a linearised vector. Returns primers per fragment.
  /* The annealing part of an assembly primer is sized the way a primer is
     sized anywhere else: the shortest length that reaches the target melting
     temperature, within sensible bounds. A fixed length instead gives a
     45 degree primer on an AT-rich end and a 70 degree one on a GC-rich end,
     which is how one primer of a pair ends up mis-priming before the other
     has annealed. Pass opts.annealTm to size them; opts.anneal keeps a fixed
     length. */
  /* The two annealing lengths of one insert, chosen together.

     Only the 3' annealing part is moved. The homology arm is untouched, so
     the junction, the overlap and the finished construct are all exactly
     what they were — the primer simply reaches a few bases further into the
     insert, or a few bases less far. Nothing about the assembly depends on
     where the annealing part stops, which is why this is free to optimise
     and the arm is not.

     Chosen apart, the two lengths each land as near the target as their own
     end allows, and a GC-rich end that bottoms out at 66 is then paired
     with an AT-rich end sitting at 60 — six degrees apart in one tube.
     Chosen together, the end with room moves towards the end without it.
     The match is weighted above the nominal target because once the target
     is out of reach for one side, a pair that agrees runs cleanly at one
     temperature and a pair that disagrees runs well at none. */
  E.annealPair = function (seq, opts) {
    opts = opts || {};
    const target = opts.annealTm || 60;
    const lo = opts.minAnneal || 16, hi = opts.maxAnneal || 36;
    const cond = opts.conditions || E.PCR_CONDITIONS;
    const s = E.cleanSeq(seq);
    const ends = fromStart => {
      const list = [];
      for (let n = lo; n <= Math.min(hi, s.length); n++) {
        const part = fromStart ? s.slice(0, n) : s.slice(s.length - n);
        const tm = E.tmNN(part, cond);
        if (isFinite(tm)) list.push({ n: n, tm: tm });
      }
      return list;
    };
    const fs = ends(true), rs = ends(false);
    if (!fs.length || !rs.length) return null;
    /* The gain here is in moving the end that HAS room towards the end that
       does not. Moving the stuck end further out as well would be a loss
       dressed up as a match, so the worst primer of the pair is never
       allowed to sit further from the target than it would have on its own.
       Within that ceiling, closing the gap is what is optimised. */
    const nearest = list => list.reduce((a, b) =>
      Math.abs(b.tm - target) < Math.abs(a.tm - target) - 1e-9 ? b : a);
    const ceiling = Math.max(Math.abs(nearest(fs).tm - target),
                             Math.abs(nearest(rs).tm - target)) + 0.5;
    let best = null;
    fs.forEach(x => rs.forEach(y => {
      if (Math.max(Math.abs(x.tm - target), Math.abs(y.tm - target)) > ceiling) return;
      const c = Math.abs(x.tm - target) + Math.abs(y.tm - target) + 2 * Math.abs(x.tm - y.tm);
      /* the shorter pair breaks a tie: no reason to buy bases that buy nothing */
      if (!best || c < best.c - 1e-9 || (Math.abs(c - best.c) < 1e-9 && x.n + y.n < best.f.n + best.r.n))
        best = { c: c, f: x, r: y };
    }));
    return best ? { f: best.f.n, r: best.r.n, tmF: best.f.tm, tmR: best.r.tm } : null;
  };

  E.annealLength = function (seq, fromStart, opts) {
    opts = opts || {};
    const target = opts.annealTm || 60;
    const lo = opts.minAnneal || 16, hi = opts.maxAnneal || 36;
    /* the same salt model the caller will report the Tm under, or the
       length is chosen on one scale and shown on another */
    const cond = opts.conditions || E.PCR_CONDITIONS;
    const s = E.cleanSeq(seq);
    /* the length whose melting temperature lands closest to the target, and
       the shortest of those when two are equally close. Taking the first
       length to pass the target instead would stop at the 16-base floor on a
       GC-rich end (70 degrees) while an AT-rich end grew to the ceiling
       (56 degrees) — the mismatched pair this exists to avoid. */
    let best = null;
    for (let n = lo; n <= Math.min(hi, s.length); n++) {
      const part = fromStart ? s.slice(0, n) : s.slice(s.length - n);
      const tm = E.tmNN(part, cond);
      if (!isFinite(tm)) continue;
      const d = Math.abs(tm - target);
      if (!best || d < best.d - 1e-9) best = { n: n, d: d, tm: tm };
    }
    return best ? best.n : Math.min(20, s.length);
  };

  E.designAssemblyPrimers = function (vectorLinear, inserts, opts) {
    opts = opts || {};
    const ov = opts.overlap || (opts.method === 'infusion' ? 15 : 20);
    const anneal = opts.anneal || 20;
    const auto = !!opts.annealTm;
    const vec = E.cleanSeq(vectorLinear);
    // order: vector -> insert1 -> insert2 ... -> back to vector
    const chain = [{ name: 'vector', seq: vec }].concat(
      inserts.map((s, i) => ({ name: s.name || ('insert ' + (i + 1)), seq: E.cleanSeq(s.seq) })));
    const primers = [];
    for (let i = 1; i < chain.length; i++) {
      const prev = chain[i - 1], cur = chain[i];
      const next = chain[(i + 1) % chain.length];
      const fTail = prev.seq.slice(-ov);
      const cond = opts.conditions || E.PCR_CONDITIONS;
      const forceA = (opts.forceAnneal && i === 1) ? opts.forceAnneal : {};
      /* both ends together, unless the caller has fixed one of them */
      const pair = (auto && !forceA.f && !forceA.r)
        ? E.annealPair(cur.seq, Object.assign({}, opts, { conditions: cond })) : null;
      const fLen = forceA.f || (pair ? pair.f
        : (auto ? E.annealLength(cur.seq, true, Object.assign({}, opts, { conditions: cond })) : anneal));
      const rLen = forceA.r || (pair ? pair.r
        : (auto ? E.annealLength(cur.seq, false, Object.assign({}, opts, { conditions: cond })) : anneal));
      const fAnneal = cur.seq.slice(0, fLen);
      const rTailSrc = next.seq.slice(0, ov);
      const rTail = E.revcomp(rTailSrc);
      const rAnneal = E.revcomp(cur.seq.slice(cur.seq.length - rLen));
      primers.push({
        fragment: cur.name,
        forward: {
          name: cur.name + '-F', seq: fTail + fAnneal,
          tail: fTail, anneal: fAnneal, length: fTail.length + fAnneal.length,
          tmAnneal: E.tmNN(fAnneal, cond),
          tmFull: E.tmNN(fTail + fAnneal, cond)
        },
        reverse: {
          name: cur.name + '-R', seq: rTail + rAnneal,
          tail: rTail, anneal: rAnneal, length: rTail.length + rAnneal.length,
          tmAnneal: E.tmNN(E.revcomp(rAnneal), cond),
          tmFull: E.tmNN(rTail + rAnneal, cond)
        }
      });
    }
    return { overlap: ov, primers };
  };

  /* ======================= OPENING A CIRCULAR VECTOR =======================
     Gibson, HiFi and In-Fusion take their homology arms from the ENDS of the
     linear vector, so a circular plasmid has to be opened at a definite place
     first. Opening it anywhere else than where you will really cut gives arms
     that match nothing in the tube, so this is never done implicitly.
  ==================================================================== */
  // Blunt opening: the circle rewritten to start at `at` (0-based).
  E.openAt = function (seq, at) {
    seq = E.cleanSeq(seq);
    const L = seq.length;
    if (!L) return '';
    const i = ((at % L) + L) % L;
    return seq.slice(i) + seq.slice(0, i);
  };

  /* Inverse PCR: amplify the backbone outward, removing [a, b) (0-based,
     half-open; a === b removes nothing and simply opens the circle there).
     The product runs from b round to a, which is the linear vector. */
  E.inversePCR = function (seq, a, b, opts) {
    opts = opts || {};
    const anneal = opts.anneal || 20;
    seq = E.cleanSeq(seq);
    const L = seq.length;
    if (!L) return { ok: false, error: 'No sequence.' };
    a = ((a % L) + L) % L; b = ((b % L) + L) % L;
    const span = (from, to) => { const n = ((to - from) % L + L) % L || (from === to ? L : 0); return (seq + seq).substr(from, n); };
    const backbone = span(b, a);
    if (backbone.length < anneal * 2 + 10) return { ok: false, error: 'The backbone left after that deletion is too short to amplify.' };
    const removed = L - backbone.length;
    const cond = opts.conditions || E.PCR_CONDITIONS;
    /* Primer length is chosen to hit a target Tm, not fixed: a backbone
       primer landing in an AT-rich stretch needs to be longer, and a fixed
       20-mer there melts far too low to amplify several kb reliably. */
    const target = opts.targetTm || 60;
    const minLen = opts.minLen || 18, maxLen = Math.min(opts.maxLen || 36, Math.floor(backbone.length / 2));
    /* A primer that matches the plasmid twice gives two products, and on a
       construct carrying a repeated element — two TEV sites, two linkers,
       a duplicated tag — the stretch beside the opening is exactly where
       that happens. Both ends of a backbone primer grow away from the
       opening, so a few more bases reach into sequence the copies do not
       share: length is therefore chosen for a unique site first and for the
       target melting temperature second. Searching a little past the usual
       ceiling is allowed when that is what uniqueness costs. */
    const hardMax = Math.min(Math.max(maxLen, opts.uniqueMaxLen || 45), Math.floor(backbone.length / 2));
    const dbl = seq + seq;
    const sitesOf = pr => {
      try { return E.primerSites(seq, pr, { circular: true }).length; } catch (e) { return 1; }
    };
    function best(makeBind, toPrimer) {
      let pick = null, uniq = null;
      for (let n = minLen; n <= hardMax; n++) {
        const bind = makeBind(n);
        if (bind.length < n) break;
        const tm = E.tmNN(bind, cond);
        const score = Math.abs(tm - target);
        const cand = { bind: bind, tm: tm, score: score, length: n };
        if (n <= maxLen && (!pick || score < pick.score)) pick = cand;
        if ((!uniq || score < uniq.score) && sitesOf(toPrimer(bind)) === 1) uniq = cand;
      }
      if (uniq) return Object.assign({}, uniq, { unique: true, lengthened: !!(pick && uniq.length > pick.length) });
      return pick ? Object.assign({}, pick, { unique: false, lengthened: false }) : null;
    }
    const force = opts.forceLen || {};
    const fixed = (n, make, toPrimer) => {
      const bind = make(n);
      return { bind: bind, tm: E.tmNN(bind, cond), score: 0, length: n,
               unique: (function () { try { return E.primerSites(seq, toPrimer(bind), { circular: true }).length === 1; }
                                      catch (e) { return true; } })(), lengthened: false, forced: true };
    };
    const mkF = n => dbl.substr(b, n), mkR = n => dbl.substr(((a - n) % L + L) % L, n);
    let f = force.f ? fixed(force.f, mkF, x => x) : best(mkF, x => x);
    let r = force.r ? fixed(force.r, mkR, x => E.revcomp(x)) : best(mkR, x => E.revcomp(x));
    if (!f || !r) return { ok: false, error: 'The backbone is too short to design primers on.' };

    /* The two primers go in one tube and share one annealing temperature,
       so choosing each one's length on its own distance from the target is
       the wrong objective: it can leave a pair six degrees apart, where no
       single temperature suits both. Where the target is reachable this
       changes nothing — both sit on it and the spread is already zero. It
       matters when the target is NOT reachable, which is the common case at
       a GC-rich opening, because the shortest length allowed already melts
       above it. Then a matched pair a few degrees high runs cleanly at one
       temperature, and a pair split between the target and the floor runs
       at neither. Matching is therefore weighted above the nominal target,
       and uniqueness still outranks both. */
    if (!force.f && !force.r) {
      const cands = (make, toPrimer) => {
        const list = [];
        for (let n = minLen; n <= hardMax; n++) {
          const bind = make(n);
          if (bind.length < n) break;
          list.push({ bind: bind, tm: E.tmNN(bind, cond), length: n,
                      unique: sitesOf(toPrimer(bind)) === 1 });
        }
        return list;
      };
      const fs = cands(mkF, x => x), rs = cands(mkR, x => E.revcomp(x));
      /* Uniqueness is judged per primer, not across the pair: one side of
         an opening can sit in a repeat where no length binds in one place
         only, and that is no reason to stop insisting on it for the other
         side, which is free to be unique. */
      const only = list => list.some(x => x.unique) ? list.filter(x => x.unique) : list;
      const fsOk = only(fs), rsOk = only(rs);
      const cost = (x, y) => Math.abs(x.tm - target) + Math.abs(y.tm - target) +
                             2 * Math.abs(x.tm - y.tm);
      let bestPair = null;
      fsOk.forEach(x => {
        rsOk.forEach(y => {
          const c = cost(x, y);
          if (!bestPair || c < bestPair.c) bestPair = { c: c, f: x, y: y };
        });
      });
      if (bestPair) {
        const keep = (old, pickd) => Object.assign({}, pickd, {
          score: Math.abs(pickd.tm - target),
          lengthened: !!(old && pickd.length > old.length), balanced: true });
        f = keep(f, bestPair.f);
        r = keep(r, bestPair.y);
      }
    }
    return {
      ok: true,
      backbone: { seq: backbone, length: backbone.length },
      removed: removed,
      primers: {
        forward: { name: 'bb-F', seq: f.bind, length: f.length, tm: f.tm, unique: f.unique, lengthened: f.lengthened },
        reverse: { name: 'bb-R', seq: E.revcomp(r.bind), binds: r.bind, length: r.length, tm: r.tm,
                   unique: r.unique, lengthened: r.lengthened }
      }
    };
  };

  /* ======================= GATEWAY (att-site recombination) =======================
     Integrase recombination happens inside a short region that the two
     partners share base for base (the 7-bp overlap core and the identical
     bases immediately around it). Because both molecules read the same there,
     the exact crossover point inside that block cannot change the product —
     so the product is cut and joined from the USER'S OWN two sequences and
     nothing is pasted in from a stored att-site table. The motifs below are
     used only to FIND candidate sites when a record carries no att features
     of its own; if a motif were wrong the site would simply not be found,
     which is visible, instead of a wrong product being built silently.
     attB x attP -> attL + attR (BP).   attL x attR -> attB + attP (LR).
  ==================================================================== */
  // search motifs only — every att1 site (B1/P1/L1/R1) contains the first,
  // every att2 site the second
  E.ATT_PROBES = { 1: 'TTTGTACAAAAAAGC', 2: 'TTGTACAAAGT' };

  // longest block of exact identity between two strings (>= minLen)
  E.longestCommonBlock = function (x, y, minLen) {
    minLen = minLen || 7;
    let prev = new Int32Array(y.length + 1), cur = new Int32Array(y.length + 1);
    let best = null;
    for (let i = 1; i <= x.length; i++) {
      for (let j = 1; j <= y.length; j++) {
        cur[j] = x[i - 1] === y[j - 1] ? prev[j - 1] + 1 : 0;
        if (cur[j] > (best ? best.len : minLen - 1)) best = { len: cur[j], aEnd: i, bEnd: j };
      }
      const t = prev; prev = cur; cur = t; cur.fill(0);
    }
    return best ? { len: best.len, a: best.aEnd - best.len, b: best.bEnd - best.len } : null;
  };

  // Find att sites of family 1 and 2 in a sequence. Prefers the record's own
  // annotation (features labelled attB1/attP2/...), falls back to the motifs.
  E.findAttSites = function (seq, opts) {
    opts = opts || {};
    seq = E.cleanSeq(seq);
    const L = seq.length, out = [];
    (opts.features || []).forEach(f => {
      const m = /att([BPLR])([12])/i.exec(f.label || '');
      if (!m) return;
      out.push({ family: +m[2], kind: m[1].toUpperCase(), start: f.start, end: f.end > f.start ? f.end : L,
                 strand: f.strand === -1 ? -1 : 1, source: 'annotation', label: f.label });
    });
    if (out.length) return out;
    [1, 2].forEach(fam => {
      const probe = E.ATT_PROBES[fam], rc = E.revcomp(probe);
      const hay = opts.circular ? seq + seq.slice(0, probe.length) : seq;
      for (let i = 0; i + probe.length <= hay.length; i++) {
        if (i >= L) break;
        if (hay.startsWith(probe, i)) out.push({ family: fam, start: i, end: i + probe.length, strand: 1, source: 'motif' });
        else if (hay.startsWith(rc, i)) out.push({ family: fam, start: i, end: i + probe.length, strand: -1, source: 'motif' });
      }
    });
    return out;
  };

  function rotate(seq, at) { at = ((at % seq.length) + seq.length) % seq.length; return seq.slice(at) + seq.slice(0, at); }

  /* Two-point recombination between two circular molecules.
     a1/a2 and b1/b2 are the starts, in each sequence, of the shared block at
     site 1 and site 2. Returns both reciprocal circular products, each built
     only from A and B. */
  E.recombineTwoPoint = function (A, B, a1, a2, b1, b2) {
    const la = A.length, lb = B.length;
    const span = (s, from, to) => { const n = ((to - from) % s.length + s.length) % s.length; return (s + s).substr(((from % s.length) + s.length) % s.length, n); };
    return {
      product: span(A, a1, a2) + span(B, b2, b1),       // A's site1..site2 piece joined to B's complement piece
      byproduct: span(A, a2, a1) + span(B, b1, b2)
    };
  };

  /* opts: {aFeatures, bFeatures, aCircular, bCircular, window}
     A is the molecule contributing the cassette (entry clone for LR, the attB
     insert for BP); B is the vector contributing the backbone. */
  E.gatewayRecombine = function (Aseq, Bseq, opts) {
    opts = opts || {};
    const A0 = E.cleanSeq(Aseq);
    let B0 = E.cleanSeq(Bseq);
    if (!A0.length || !B0.length) return { ok: false, error: 'Both molecules need a sequence.' };
    const win = opts.window || 260;
    const aSites = E.findAttSites(A0, { features: opts.aFeatures, circular: opts.aCircular !== false });
    let bSites = E.findAttSites(B0, { features: opts.bFeatures, circular: opts.bCircular !== false });
    const pick = (list, fam) => list.filter(s => s.family === fam);
    const notes = [];
    for (const [who, list] of [['first molecule', aSites], ['second molecule', bSites]]) {
      if (!pick(list, 1).length || !pick(list, 2).length) {
        return { ok: false, error: 'No att1 and att2 sites found in the ' + who + '. Annotate the att sites (a feature whose name contains attL1, attR2, …) or check the sequence.' };
      }
    }
    if (aSites.some(s => s.source === 'motif') || bSites.some(s => s.source === 'motif')) {
      notes.push('att sites were located by sequence motif because the records carry no att features; annotate them (attL1, attR2, …) to have the sites taken from your own annotation instead.');
    }
    // crossover: the block the two partners share around each site
    function crossover(fam, flip) {
      const as = pick(aSites, fam)[0];
      const bl = flip ? E.revcomp(B0) : B0;
      const bf = (opts.bFeatures && flip)
        ? opts.bFeatures.map(f => ({ label: f.label, strand: f.strand === -1 ? 1 : -1,
            start: B0.length - (f.end > f.start ? f.end : B0.length), end: B0.length - f.start }))
        : opts.bFeatures;
      const bs2 = E.findAttSites(bl, { features: bf, circular: opts.bCircular !== false });
      const bs = pick(bs2, fam)[0];
      if (!as || !bs) return null;
      const aw = Math.max(0, as.start - win), bw = Math.max(0, bs.start - win);
      const ax = A0.slice(aw, Math.min(A0.length, as.end + win));
      const bx = bl.slice(bw, Math.min(bl.length, bs.end + win));
      const blk = E.longestCommonBlock(ax, bx, 7);
      if (!blk) return null;
      return { a: aw + blk.a, b: bw + blk.b, len: blk.len, bSeq: bl };
    }
    let flip = false, c1 = crossover(1, false), c2 = crossover(2, false);
    // the vector may be annotated on the other strand; try it flipped and keep
    // whichever orientation gives the longer shared blocks at both sites
    const f1 = crossover(1, true), f2 = crossover(2, true);
    const score = (x, y) => (x ? x.len : 0) + (y ? y.len : 0);
    if (score(f1, f2) > score(c1, c2)) { flip = true; c1 = f1; c2 = f2; }
    if (!c1 || !c2) return { ok: false, error: 'The att sites do not share a recombination core — check that these are matching partners (attL × attR, or attB × attP).' };
    const Bm = flip ? E.revcomp(B0) : B0;
    if (flip) notes.push('The second molecule was read on the other strand so its att sites match the first.');
    const aLinear = opts.aCircular === false;
    if (aLinear && c1.a > c2.a) return { ok: false, error: 'On a linear insert the att1 site must come before the att2 site; reverse-complement it first.' };
    const r = E.recombineTwoPoint(A0, Bm, c1.a, c2.a, c1.b, c2.b);
    // Start the expression clone at the vector's own base 1 when that base
    // survives, so the numbering a user already knows is kept; likewise the
    // by-product keeps the first molecule's origin.
    const lenA1 = ((c2.a - c1.a) % A0.length + A0.length) % A0.length;
    const lenB1 = ((c1.b - c2.b) % Bm.length + Bm.length) % Bm.length;
    let prod = r.product;
    const bZero = ((0 - c2.b) % Bm.length + Bm.length) % Bm.length;
    if (!aLinear && bZero < lenB1) prod = rotate(prod, lenA1 + bZero);
    let byp = r.byproduct;
    const lenA2 = ((c1.a - c2.a) % A0.length + A0.length) % A0.length;
    const aZero = ((0 - c2.a) % A0.length + A0.length) % A0.length;
    if (aZero < lenA2) byp = rotate(byp, aZero);
    return {
      ok: true,
      product: { seq: prod, length: prod.length, circular: true },
      byproduct: aLinear ? null : { seq: byp, length: byp.length, circular: true },
      cores: [{ length: c1.len, atA: c1.a, atB: c1.b }, { length: c2.len, atA: c2.a, atB: c2.b }],
      flipped: flip, notes: notes,
      aSites: aSites, bSites: bSites
    };
  };

  // kept for callers that name the reaction; both are the same swap
  E.gatewayLR = function (entrySeq, destSeq, opts) {
    const r = E.gatewayRecombine(entrySeq, destSeq, opts);
    if (r.ok) r.note = 'attL × attR: the cassette between the entry clone’s att sites replaces the region between the destination vector’s att sites. Every base of the product comes from the two sequences you supplied.';
    return r;
  };
  E.gatewayBP = function (attBInsert, donorSeq, opts) {
    opts = Object.assign({ aCircular: false }, opts || {});
    const r = E.gatewayRecombine(attBInsert, donorSeq, opts);
    if (r.ok) r.note = 'attB × attP: the insert between the attB sites moves into the donor vector, giving the entry clone. Every base of the product comes from the two sequences you supplied.';
    return r;
  };

  /* ======================= TOPO / TA / GC CLONING ======================= */
  E.topoClone = function (insert, opts) {
    opts = opts || {};
    const mode = opts.mode || 'blunt';   // 'blunt' | 'ta' | 'directional'
    const vec = E.cleanSeq(opts.vectorSeq || '');
    let ins = E.cleanSeq(insert);
    const notes = [];
    if (mode === 'ta') {
      notes.push("TA cloning: Taq leaves a single untemplated 3'-A on each end of the PCR product, which pairs with the 3'-T overhangs of the T-tailed vector.");
      if (!/^[ACGT]/.test(ins)) notes.push('Insert looks empty.');
    } else if (mode === 'directional') {
      notes.push("Directional TOPO: the forward primer must start with CACC so the GTGG overhang of the vector captures the insert in one orientation.");
      if (!ins.startsWith('CACC')) {
        return { ok: false, error: "Directional TOPO needs the insert to begin with CACC (add it to the 5' end of the forward primer).", notes };
      }
    } else {
      notes.push('Blunt TOPO: the PCR product must be blunt-ended (use a proofreading polymerase, not Taq).');
    }
    const product = vec ? vec + ins : ins;
    return {
      ok: true, mode,
      product: { seq: product, length: product.length, circular: !!vec },
      orientation: mode === 'directional' ? 'fixed (CACC/GTGG)' : 'both orientations expected (~50:50)',
      notes
    };
  };

  /* ======================= SITE-DIRECTED MUTAGENESIS =======================
     Design back-to-back or overlapping primers carrying a point change,
     insertion or deletion, and report the mutated template.
  ==================================================================== */
  E.sdm = function (template, opts) {
    const seq = E.cleanSeq(template);
    const pos = opts.position;            // 0-based position to change
    const kind = opts.kind || 'substitute'; // substitute | insert | delete
    const newBases = E.cleanSeq(opts.newBases || '');
    const flank = opts.flank || 15;
    if (pos < 0 || pos > seq.length) return { ok: false, error: 'Position is outside the sequence.' };
    if (kind !== 'delete' && !newBases.length) return { ok: false, error: 'Enter the new bases.' };
    const span = kind === 'substitute' ? newBases.length : (kind === 'delete' ? (opts.deleteLength || newBases.length || 1) : 0);
    if (pos + span > seq.length) return { ok: false, error: 'The change runs past the end of the sequence.' };
    let mutated, changeDesc;
    if (kind === 'substitute') {
      const n = newBases.length || 1;
      mutated = seq.slice(0, pos) + newBases + seq.slice(pos + n);
      changeDesc = seq.slice(pos, pos + n) + ' -> ' + newBases + ' at ' + (pos + 1);
    } else if (kind === 'insert') {
      mutated = seq.slice(0, pos) + newBases + seq.slice(pos);
      changeDesc = 'insert ' + newBases + ' after ' + pos;
    } else {
      const n = opts.deleteLength || newBases.length || 1;
      mutated = seq.slice(0, pos) + seq.slice(pos + n);
      changeDesc = 'delete ' + n + ' nt at ' + (pos + 1);
    }
    // primers on the mutated sequence, centred on the change
    const changeLen = kind === 'delete' ? 0 : (newBases.length || 1);
    let fwd;
    if (opts.circular && mutated.length > changeLen + 2 * flank) {
      // on a plasmid the primer can run across the origin
      const M = mutated.length, start = ((pos - flank) % M + M) % M, n = changeLen + 2 * flank;
      fwd = (mutated + mutated).slice(start, start + n);
    } else {
      const left = Math.max(0, pos - flank);
      const right = Math.min(mutated.length, pos + changeLen + flank);
      fwd = mutated.slice(left, right);
    }
    const rev = E.revcomp(fwd);
    return {
      ok: true, kind, change: changeDesc,
      mutated: { seq: mutated, length: mutated.length },
      primers: {
        forward: { name: 'SDM-F', seq: fwd, length: fwd.length, tm: E.tmNN(fwd, E.PCR_CONDITIONS) },
        reverse: { name: 'SDM-R', seq: rev, length: rev.length, tm: E.tmNN(rev, E.PCR_CONDITIONS) }
      },
      mismatchCentred: true
    };
  };

  /* ======================= SILENT MUTATION =======================
     Add or remove a restriction site inside a CDS without changing protein.
  ==================================================================== */
  /* A recognition site can straddle up to three codons, so a silent edit may
     need more than one codon changed at once (e.g. GAG TTT -> GAA TTC creates
     GAATTC while keeping Glu-Phe). We therefore search windows of 1, 2 and 3
     consecutive codons, testing each candidate against a local context slice
     rather than rescanning the whole CDS.                                    */
  function enumerateWindow(seq, table, ci, w, cb) {
    const synLists = [];
    let combos = 1;
    for (let k = 0; k < w; k++) {
      const codon = seq.substr((ci + k) * 3, 3);
      const aa = E.translate(codon, { table });
      const syn = E.synonymousCodons(aa, table);
      if (!syn.length) return;
      synLists.push(syn);
      combos *= syn.length;
      if (combos > 500) return;          // guard against pathological blow-up
    }
    const idx = new Array(w).fill(0);
    for (;;) {
      let variant = '', changedCodons = 0;
      for (let k = 0; k < w; k++) {
        const c = synLists[k][idx[k]];
        variant += c;
        if (c !== seq.substr((ci + k) * 3, 3)) changedCodons++;
      }
      if (changedCodons > 0 && cb(variant, changedCodons) === false) return;
      let k = w - 1;
      while (k >= 0) { idx[k]++; if (idx[k] < synLists[k].length) break; idx[k] = 0; k--; }
      if (k < 0) return;
    }
  }

  E.silentMutations = function (cds, enzymeName, opts) {
    opts = opts || {};
    const table = opts.table || 1;
    const mode = opts.mode || 'remove';   // 'remove' | 'add'
    const maxResults = opts.maxResults || 40;
    const maxWindow = opts.maxWindow || 3;
    const enz = E.enzymeByName(enzymeName);
    if (!enz) return { ok: false, error: 'Unknown enzyme ' + enzymeName };
    const seq = E.cleanSeq(cds);
    const protein = E.translate(seq, { table });
    const nCodons = Math.floor(seq.length / 3);
    const pad = enz.site.length - 1 + (enz.typeIIS ? 14 : 0);
    const baseSites = E.findSites(seq, enz, false);
    const results = [];

    // Search one window; returns the best candidate (fewest bases changed).
    function searchWindow(ci, w, want) {
      const wStart = ci * 3, wEnd = (ci + w) * 3;
      if (wEnd > seq.length) return null;
      const ctxStart = Math.max(0, wStart - pad);
      const ctxEnd = Math.min(seq.length, wEnd + pad);
      const prefix = seq.slice(ctxStart, wStart);
      const suffix = seq.slice(wEnd, ctxEnd);
      const origLocal = prefix + seq.slice(wStart, wEnd) + suffix;
      const origCount = E.findSites(origLocal, enz, false).length;
      let best = null;
      enumerateWindow(seq, table, ci, w, (variant) => {
        const local = prefix + variant + suffix;
        const n = E.findSites(local, enz, false).length;
        const improved = want === 'add' ? n > origCount : n < origCount;
        if (!improved) return;
        const cand = seq.slice(0, wStart) + variant + seq.slice(wEnd);
        const g = E.findSites(cand, enz, false).length;
        const globalOk = want === 'add' ? g > baseSites.length : g < baseSites.length;
        if (!globalOk) return;
        if (E.translate(cand, { table }) !== protein) return;
        const changed = diffPositions(seq, cand);
        if (!best || changed.length < best.changed.length) {
          best = {
            codonIndex: ci, aaPosition: ci + 1, codonSpan: w,
            from: seq.slice(wStart, wEnd), to: variant,
            aa: E.translate(seq.slice(wStart, wEnd), { table }),
            seq: cand, changed, siteCount: g
          };
        }
      });
      return best;
    }

    if (mode === 'remove') {
      for (const s of baseSites) {
        const first = Math.floor(s.pos / 3);
        const last = Math.floor((s.pos + enz.site.length - 1) / 3);
        let fixed = null;
        for (let w = 1; w <= maxWindow && !fixed; w++) {
          for (let ci = Math.max(0, first - w + 1); ci <= last && !fixed; ci++) {
            const r = searchWindow(ci, w, 'remove');
            if (r) fixed = Object.assign({ sitePos: s.pos }, r);
          }
        }
        results.push(fixed || { sitePos: s.pos, error: 'No synonymous change removes this site.' });
      }
      return { ok: true, mode, enzyme: enzymeName, siteCount: baseSites.length, results };
    }

    // mode 'add'
    for (let w = 1; w <= maxWindow; w++) {
      for (let ci = 0; ci + w <= nCodons; ci++) {
        const r = searchWindow(ci, w, 'add');
        if (r) {
          // skip if we already have a hit at this position with fewer changes
          const dup = results.find(x => x.codonIndex === ci);
          if (!dup) results.push(r);
          else if (r.changed.length < dup.changed.length) Object.assign(dup, r);
        }
        if (results.length >= maxResults) break;
      }
      if (results.length >= maxResults) break;
    }
    results.sort((a, b) => a.changed.length - b.changed.length || a.codonIndex - b.codonIndex);
    return { ok: true, mode, enzyme: enzymeName, results };
  };
  function diffPositions(a, b) {
    const out = [];
    for (let i = 0; i < Math.min(a.length, b.length); i++) if (a[i] !== b[i]) out.push({ pos: i, from: a[i], to: b[i] });
    return out;
  }

  /* ======================= PRIMER DESIGN =======================
     Pick primers flanking a target region, optimising length for a target Tm,
     then score hairpin / self-dimer / 3'-end stability.
  ==================================================================== */
  E.designPrimers = function (template, opts) {
    opts = opts || {};
    const seq = E.cleanSeq(template);
    const start = opts.start || 0;
    const end = opts.end != null ? opts.end : seq.length;
    const targetTm = opts.targetTm || 60;
    const minLen = opts.minLen || 18, maxLen = opts.maxLen || 30;
    const cond = opts.conditions || E.PCR_CONDITIONS;
    const fTail = E.cleanSeq(opts.forwardTail || '');
    const rTail = E.cleanSeq(opts.reverseTail || '');

    let bestF = null, bestR = null;
    for (let len = minLen; len <= maxLen; len++) {
      const f = seq.substr(start, len);
      if (f.length < len) break;
      const tm = E.tmNN(f, cond);
      const score = Math.abs(tm - targetTm);
      if (!bestF || score < bestF.score) bestF = { seq: f, tm, score, length: len };
    }
    for (let len = minLen; len <= maxLen; len++) {
      const s = seq.substr(Math.max(0, end - len), Math.min(len, end));
      const r = E.revcomp(s);
      const tm = E.tmNN(r, cond);
      const score = Math.abs(tm - targetTm);
      if (!bestR || score < bestR.score) bestR = { seq: r, tm, score, length: r.length, binds: s };
    }
    if (!bestF || !bestR) return { ok: false, error: 'Region too short to design primers.' };
    const fFull = fTail + bestF.seq, rFull = rTail + bestR.seq;
    const product = E.pcr(seq, fFull, rFull, { anchor: Math.min(18, bestF.length) });
    return {
      ok: true,
      forward: primerReport('F', fFull, bestF.seq, fTail, cond),
      reverse: primerReport('R', rFull, bestR.seq, rTail, cond),
      tmDifference: Math.abs(bestF.tm - bestR.tm),
      heterodimer: E.heteroDimer(fFull, rFull),
      amplicon: product.length ? { length: product[0].length, seq: product[0].seq } : null,
      annealingSuggestion: Math.round(Math.min(bestF.tm, bestR.tm) - 3)
    };
  };
  /* ============================================================
     THE MELTING TEMPERATURE OF A PRIMER THAT HAS AN OVERHANG

     A primer with a 5' tail — an assembly arm, a restriction site, a tag —
     does not anneal along its whole length in the first cycles. Only the 3'
     part that matches the template binds, and that part alone sets the
     temperature the reaction has to be run at. The whole primer's Tm is
     nearly always the higher of the two, often by ten degrees or more, and
     annealing at it means annealing above what actually holds: the primer
     does not bind and the reaction gives nothing. (Nearly always, not
     always: an AT-rich arm in front of a short GC-rich annealing stretch
     can come out a few tenths lower. It happens in well under one primer in
     a hundred, and it is a reason to report both numbers rather than to
     assume which is larger.)

     Both numbers matter, at different times. The annealing part governs the
     first cycles, when the template is the original. Once the tail has been
     copied into the product, later cycles anneal along the whole primer and
     the full-length Tm governs those. Protocols that use two annealing
     temperatures do exactly this; one that uses a single temperature uses
     the annealing part's.

     Which part anneals is a fact about a primer AND a template, not about
     the primer alone, so this takes the annealing portion from whoever
     knows it — the designer that built the primer, or a binding site found
     on a given sequence — and says in `basis` where the split came from
     rather than inventing one.
  ============================================================ */
  E.primerTm = function (primer, opts) {
    opts = opts || {};
    const full = E.cleanSeq(primer || '').replace(/U/g, 'T');
    const cond = opts.conditions;
    const out = { seq: full, length: full.length, tmFull: full ? E.tmNN(full, cond) : NaN,
                  basis: 'full' };
    if (!full) { out.tm = NaN; out.anneal = ''; out.annealLen = 0; out.tailLen = 0; out.tailed = false; return out; }
    /* the annealing portion: given outright, or as the number of 3' bases
       of it that bind */
    let anneal = null;
    if (typeof opts.anneal === 'string' && opts.anneal) anneal = E.cleanSeq(opts.anneal).replace(/U/g, 'T');
    else if (typeof opts.annealLen === 'number' && opts.annealLen > 0 && opts.annealLen <= full.length)
      anneal = full.slice(full.length - opts.annealLen);
    /* An annealing portion that is not the 3' end of this primer describes
       some other oligo; rather than quietly measure the wrong thing, fall
       back to the whole primer and say so. */
    if (anneal && full.slice(full.length - anneal.length) !== anneal) anneal = null;
    if (!anneal || anneal.length === full.length) {
      out.tm = out.tmFull; out.anneal = full; out.annealLen = full.length;
      out.tailLen = 0; out.tailed = false;
      if (anneal) out.basis = opts.basis || 'full';
      return out;
    }
    out.anneal = anneal;
    out.annealLen = anneal.length;
    out.tailLen = full.length - anneal.length;
    out.tail = full.slice(0, out.tailLen);
    out.tailed = true;
    out.tm = E.tmNN(anneal, cond);
    out.basis = opts.basis || 'given';
    return out;
  };

  /* ============================================================
     WORKING OUT WHERE A PRIMER'S OVERHANG ENDS

     A primer saved before its design was recorded is just a string of
     bases, and nothing in the string says which part of it anneals. Two
     kinds of evidence can settle it, and both are statements about this
     primer against real sequence rather than guesses about its shape:

       1. Its 3' end anneals to a construct and its 5' end does not. That
          is the ordinary case, and primerSites already measures it.

       2. Its 5' end is found in a construct and its 3' end is not. That is
          what an assembly primer looks like against the VECTOR it was
          designed to go into: the homology arm is vector sequence, so the
          arm matches and the insert-binding part does not. It is the case
          that matters here, because the vector is usually what is open.

     What is deliberately not used is the shape of the library itself.
     Assembly primers made for one site share a 5' arm, so a common prefix
     looks like a free answer — but on a real set it runs past the arm into
     whatever the inserts happen to share. On the set this was written for,
     the arms are 15 bases and the common prefix is 18, because every insert
     also begins ATG. Three bases too many means an annealing temperature a
     degree or two low, from evidence that was never about the arm at all.
  ============================================================ */
  E.findOverhang = function (primer, records, opts) {
    opts = opts || {};
    const minAnneal = opts.minAnneal || 14;
    /* 10 bases of exact match turn up by chance about once in 200 kb of
       sequence, so a prefix at least that long is the arm and not noise */
    const minTail = opts.minTail || 10;
    const p = E.cleanSeq(primer || '').replace(/U/g, 'T');
    if (!p || p.length < minAnneal + minTail) return null;
    const recs = (records || []).filter(r => r && r.seq && r.seq.length > 20);
    if (!recs.length) return null;

    /* 1 — where it anneals, if it anneals */
    let byAnneal = null;
    for (const r of recs) {
      let sites = [];
      try { sites = E.primerSites(r.seq, p, { circular: r.topology === 'circular' }); }
      catch (e) { sites = []; }
      sites.forEach(h => {
        if (h.tail <= 0 || h.annealed < minAnneal) return;
        if (!byAnneal || h.annealed > byAnneal.annealLen)
          byAnneal = { annealLen: h.annealed, tailLen: h.tail, how: 'anneals', on: r.name || '', strand: h.strand };
      });
    }

    /* 2 — how much of its 5' end is the construct's own sequence.

       The guard matters more than the search. A backbone primer is vector
       sequence from end to end, so SOME prefix of it is always in the
       vector — all of it is — and a search for the longest prefix that
       matches will happily report the first fourteen bases as a homology
       arm and halve the temperature. An arm is a 5' end that is in the
       construct while the rest of the primer is not; if the whole primer is
       there, it anneals along its whole length and has no arm at all. */
    let byTail = null;
    for (const r of recs) {
      const s = E.cleanSeq(r.seq).replace(/U/g, 'T');
      const circ = r.topology === 'circular';
      const hay = circ ? s + s : s;
      const rc = E.revcomp(hay);
      const inThere = x => hay.indexOf(x) >= 0 || rc.indexOf(x) >= 0;
      if (inThere(p)) return null;                 // a plain primer on this construct
      /* longest first: the arm is whatever length of 5' end is really
         there. If even the longest prefix tried still matches, the boundary
         lies outside the window and where it falls cannot be seen — better
         to say nothing than to put it at the edge of the search. */
      const top = p.length - minAnneal;
      if (inThere(p.slice(0, top + 1))) continue;
      for (let n = top; n >= minTail; n--) {
        if (byTail && byTail.tailLen >= n) break;
        const pre = p.slice(0, n);
        if (pre.indexOf('N') >= 0) continue;
        const fwd = hay.indexOf(pre) >= 0, rev = rc.indexOf(pre) >= 0;
        if (!fwd && !rev) continue;
        byTail = { annealLen: p.length - n, tailLen: n, how: 'tail-matches',
                   on: r.name || '', strand: fwd ? 1 : -1 };
        break;
      }
    }

    if (!byAnneal && !byTail) return null;
    /* Where it anneals is the more direct statement, so it wins; the other
       is carried along so a disagreement can be shown rather than hidden. */
    const best = byAnneal || byTail;
    const other = byAnneal && byTail && byAnneal.annealLen !== byTail.annealLen ? byTail : null;
    return Object.assign({}, best, {
      anneal: p.slice(p.length - best.annealLen),
      tail: p.slice(0, best.tailLen),
      disagrees: other ? { annealLen: other.annealLen, tailLen: other.tailLen, on: other.on } : null
    });
  };

  function primerReport(tag, full, annealPart, tail, cond) {
    return {
      name: 'primer-' + tag,
      seq: full, length: full.length,
      tail: tail || null,
      anneal: annealPart,
      tm: E.tmNN(annealPart, cond),
      tmFull: E.tmNN(full, cond),
      gc: +E.gcContent(full).toFixed(1),
      gc3: /[GC]$/.test(full) ? 'G/C 3\' end' : 'A/T 3\' end',
      hairpin: E.hairpin(full),
      selfDimer: E.selfDimer(full),
      runs: longestRun(full),
      thermo: E.nnThermo(annealPart, cond)
    };
  }
  function longestRun(s) {
    let best = 0, cur = 1;
    for (let i = 1; i < s.length; i++) {
      if (s[i] === s[i - 1]) { cur++; if (cur > best) best = cur; } else cur = 1;
    }
    return Math.max(best, s.length ? 1 : 0);
  }

  // Hairpin: best intramolecular stem with a loop of >= 3 nt
  E.hairpin = function (s, opts) {
    opts = opts || {};
    const minLoop = opts.minLoop || 3, minStem = opts.minStem || 4;
    s = E.cleanSeq(s);
    let best = null;
    for (let i = 0; i < s.length; i++) {
      for (let j = i + minStem + minLoop; j < s.length; j++) {
        let n = 0;
        while (i + n < j - n - minLoop + 1 && j - n >= 0 &&
               pairs(s[i + n], s[j - n])) n++;
        if (n >= minStem) {
          const loop = (j - n + 1) - (i + n);
          if (loop >= minLoop) {
            const stem = s.substr(i, n);
            const dg = stemDG(stem);
            if (!best || dg < best.dG) best = { stem, stemLength: n, loopLength: loop, start: i, end: j, dG: dg };
          }
        }
      }
    }
    return best || { stemLength: 0, dG: 0, note: 'no significant hairpin' };
  };
  function pairs(a, b) {
    return (a === 'A' && b === 'T') || (a === 'T' && b === 'A') ||
           (a === 'G' && b === 'C') || (a === 'C' && b === 'G');
  }
  function stemDG(stem) {
    const th = E.nnThermo(stem);
    return th ? th.dG37 : 0;
  }

  // Self-dimer: best alignment of the oligo against its own reverse complement
  E.selfDimer = function (s) { return E.dimer(s, s); };
  E.heteroDimer = function (a, b) { return E.dimer(a, b); };
  // Slide b's reverse complement along a and score contiguous complementarity,
  // weighting matches that involve either 3' end (those actually extend).
  E.dimer = function (a, b) {
    a = E.cleanSeq(a); b = E.cleanSeq(b);
    const brc = E.revcomp(b);
    let best = { matches: 0, offset: 0, threePrime: false, dG: 0, a: '', b: '' };
    for (let off = -(brc.length - 1); off < a.length; off++) {
      let run = 0, bestRun = 0, runStartA = 0, bestStartA = 0;
      for (let i = 0; i < brc.length; i++) {
        const ai = off + i;
        if (ai < 0 || ai >= a.length) { run = 0; continue; }
        if (a[ai] === brc[i]) {
          if (run === 0) runStartA = ai;
          run++;
          if (run > bestRun) { bestRun = run; bestStartA = runStartA; }
        } else run = 0;
      }
      if (bestRun > best.matches) {
        const seg = a.substr(bestStartA, bestRun);
        const involves3 = (bestStartA + bestRun >= a.length - 2);
        best = {
          matches: bestRun, offset: off, threePrime: involves3,
          dG: bestRun >= 2 ? stemDG(seg) : 0,
          segment: seg
        };
      }
    }
    return best;
  };

  /* ======================= MULTIPLE SEQUENCE ALIGNMENT =======================
     Center-star progressive alignment built on the pairwise DP in engine.js.
     Deterministic and valid; not an implementation of Clustal/MAFFT/MUSCLE.
  ==================================================================== */
  E.msa = function (seqs, opts) {
    opts = opts || {};
    const recs = seqs.map((s, i) => ({ name: s.name || ('seq' + (i + 1)), seq: (s.seq || s).toUpperCase() }));
    if (recs.length < 2) return { ok: false, error: 'Need at least two sequences.' };
    if (recs.length === 2) {
      const a = E.align(recs[0].seq, recs[1].seq, opts);
      return { ok: true, rows: [{ name: recs[0].name, seq: a.a }, { name: recs[1].name, seq: a.b }],
               identity: a.identity, length: a.a.length, consensus: consensusOf([a.a, a.b]) };
    }
    // pick the center: highest summed identity to all others
    let center = 0, bestSum = -1;
    const simCache = {};
    for (let i = 0; i < recs.length; i++) {
      let sum = 0;
      for (let j = 0; j < recs.length; j++) {
        if (i === j) continue;
        const k = Math.min(i, j) + ':' + Math.max(i, j);
        if (simCache[k] == null) simCache[k] = E.align(recs[i].seq, recs[j].seq, opts).identity;
        sum += simCache[k];
      }
      if (sum > bestSum) { bestSum = sum; center = i; }
    }
    // align each sequence to the center, accumulating gaps into the center
    let centerAln = recs[center].seq;
    const pending = [];
    for (let i = 0; i < recs.length; i++) {
      if (i === center) continue;
      const a = E.align(centerAln.replace(/-/g, ''), recs[i].seq, opts);
      pending.push({ idx: i, centerRow: a.a, otherRow: a.b });
    }
    // merge: build a center consensus of gaps, then project each row onto it
    let merged = centerAln;
    for (const p of pending) merged = mergeGapPattern(merged, p.centerRow);
    const rows = [];
    rows[center] = { name: recs[center].name, seq: merged };
    for (const p of pending) {
      rows[p.idx] = { name: recs[p.idx].name, seq: projectRow(merged, p.centerRow, p.otherRow) };
    }
    const seqRows = rows.map(r => r.seq);
    const maxLen = Math.max(...seqRows.map(s => s.length));
    for (const r of rows) while (r.seq.length < maxLen) r.seq += '-';
    return {
      ok: true, rows, length: maxLen,
      consensus: consensusOf(rows.map(r => r.seq)),
      identity: columnIdentity(rows.map(r => r.seq)),
      method: 'center-star progressive (center: ' + recs[center].name + ')'
    };
  };
  // combine two gap patterns of the same underlying sequence
  function mergeGapPattern(a, b) {
    let out = '', i = 0, j = 0;
    while (i < a.length || j < b.length) {
      const ga = i < a.length && a[i] === '-';
      const gb = j < b.length && b[j] === '-';
      if (ga && gb) { out += '-'; i++; j++; }
      else if (ga) { out += '-'; i++; }
      else if (gb) { out += '-'; j++; }
      else { out += (i < a.length ? a[i] : b[j]); i++; j++; }
    }
    return out;
  }
  function projectRow(merged, centerRow, otherRow) {
    // walk merged; consume centerRow/otherRow together
    let out = '', ci = 0;
    for (let m = 0; m < merged.length; m++) {
      if (merged[m] === '-') {
        if (ci < centerRow.length && centerRow[ci] === '-') { out += otherRow[ci]; ci++; }
        else out += '-';
      } else {
        while (ci < centerRow.length && centerRow[ci] === '-') ci++;
        out += (ci < otherRow.length ? otherRow[ci] : '-');
        ci++;
      }
    }
    return out;
  }
  function consensusOf(rows) {
    if (!rows.length) return '';
    const L = Math.max(...rows.map(r => r.length));
    let out = '';
    for (let i = 0; i < L; i++) {
      const counts = {};
      for (const r of rows) { const c = r[i] || '-'; counts[c] = (counts[c] || 0) + 1; }
      let best = '-', bn = -1;
      for (const k in counts) if (counts[k] > bn && k !== '-') { bn = counts[k]; best = k; }
      out += (bn > rows.length / 2) ? best : (bn > 0 ? best.toLowerCase() : '-');
    }
    return out;
  }
  function columnIdentity(rows) {
    if (!rows.length) return 0;
    const L = Math.max(...rows.map(r => r.length));
    let same = 0;
    for (let i = 0; i < L; i++) {
      const first = rows[0][i];
      if (first && first !== '-' && rows.every(r => r[i] === first)) same++;
    }
    return L ? same / L * 100 : 0;
  }

  /* ======================= CONTIG ASSEMBLY =======================
     Overlap-layout-consensus for a handful of Sanger reads. Finds the best
     overlaps (including reverse-complemented reads) and merges greedily.
  ==================================================================== */
  E.assembleContigs = function (reads, opts) {
    opts = opts || {};
    const minOverlap = opts.minOverlap || 20;
    const minIdentity = opts.minIdentity || 95;
    let pool = reads.map((r, i) => ({ name: r.name || ('read' + (i + 1)), seq: E.cleanSeq(r.seq || r), merged: [r.name || ('read' + (i + 1))] }));
    let changed = true;
    while (changed && pool.length > 1) {
      changed = false;
      let best = null;
      for (let i = 0; i < pool.length; i++) {
        for (let j = 0; j < pool.length; j++) {
          if (i === j) continue;
          for (const flip of [false, true]) {
            const bseq = flip ? E.revcomp(pool[j].seq) : pool[j].seq;
            const ov = bestSuffixPrefix(pool[i].seq, bseq, minOverlap, minIdentity);
            if (ov && (!best || ov.length > best.ov.length)) best = { i, j, flip, ov };
          }
        }
      }
      if (best) {
        const a = pool[best.i], b = pool[best.j];
        const bseq = best.flip ? E.revcomp(b.seq) : b.seq;
        const mergedSeq = a.seq + bseq.slice(best.ov.length);
        const merged = { name: a.name + '+' + b.name, seq: mergedSeq,
                         merged: a.merged.concat(b.merged.map(m => m + (best.flip ? '(rc)' : ''))) };
        pool = pool.filter((_, k) => k !== best.i && k !== best.j);
        pool.push(merged);
        changed = true;
      }
    }
    return {
      contigs: pool.map(p => ({ name: p.name, seq: p.seq, length: p.seq.length, reads: p.merged })),
      contigCount: pool.length
    };
  };
  function bestSuffixPrefix(a, b, minOverlap, minIdentity) {
    const lim = Math.min(a.length, b.length);
    for (let n = lim; n >= minOverlap; n--) {
      const s1 = a.slice(-n), s2 = b.slice(0, n);
      let m = 0;
      for (let k = 0; k < n; k++) if (s1[k] === s2[k]) m++;
      const id = m / n * 100;
      if (id >= minIdentity) return { length: n, identity: id };
    }
    return null;
  }

  /* ======================= RNA SECONDARY STRUCTURE =======================
     Nussinov-Jacobson base-pair maximisation with a minimum hairpin loop.
     This maximises pair count, not free energy, so it is a structural sketch
     rather than an MFE prediction.
  ==================================================================== */
  E.rnaFold = function (seq, opts) {
    opts = opts || {};
    const minLoop = opts.minLoop != null ? opts.minLoop : 3;
    const s = E.cleanSeq(seq).replace(/T/g, 'U');
    const n = s.length;
    if (n < 2) return { pairs: [], structure: '.'.repeat(n), pairCount: 0, sequence: s };
    if (n > 1200) return { error: 'Sequence too long for this folding routine (limit 1200 nt).' };
    const dp = [];
    for (let i = 0; i < n; i++) dp.push(new Int16Array(n));
    const canPair = (a, b) =>
      (a === 'A' && b === 'U') || (a === 'U' && b === 'A') ||
      (a === 'G' && b === 'C') || (a === 'C' && b === 'G') ||
      (a === 'G' && b === 'U') || (a === 'U' && b === 'G');
    for (let len = minLoop + 1; len < n; len++) {
      for (let i = 0; i + len < n; i++) {
        const j = i + len;
        let best = dp[i][j - 1];                       // j unpaired
        if (canPair(s[i], s[j]) && (j - i - 1) >= minLoop) {
          const v = (i + 1 <= j - 1 ? dp[i + 1][j - 1] : 0) + 1;
          if (v > best) best = v;
        }
        for (let k = i + 1; k < j; k++) {
          const v = dp[i][k] + dp[k + 1][j];
          if (v > best) best = v;
        }
        dp[i][j] = best;
      }
    }
    const pairsOut = [];
    (function trace(i, j) {
      if (i >= j) return;
      if (dp[i][j] === dp[i][j - 1]) { trace(i, j - 1); return; }
      if (canPair(s[i], s[j]) && (j - i - 1) >= minLoop &&
          dp[i][j] === (i + 1 <= j - 1 ? dp[i + 1][j - 1] : 0) + 1) {
        pairsOut.push([i, j]); trace(i + 1, j - 1); return;
      }
      for (let k = i + 1; k < j; k++) {
        if (dp[i][j] === dp[i][k] + dp[k + 1][j]) { trace(i, k); trace(k + 1, j); return; }
      }
    })(0, n - 1);
    const st = new Array(n).fill('.');
    for (const [i, j] of pairsOut) { st[i] = '('; st[j] = ')'; }
    return {
      sequence: s, pairs: pairsOut.sort((a, b) => a[0] - b[0]),
      structure: st.join(''), pairCount: pairsOut.length,
      method: 'Nussinov base-pair maximisation (not a free-energy minimisation)'
    };
  };

  /* ======================= AGAROSE GEL MODEL =======================
     Migration distance of a fragment is modelled as a linear function of
     log10(size) over the gel's resolving window, which is what an agarose gel
     approximates in its usable range. Resolving window depends on % agarose.
  ==================================================================== */
  E.GEL_RANGES = {
    0.5: [1000, 30000], 0.7: [800, 12000], 1.0: [500, 10000],
    1.2: [400, 7000], 1.5: [200, 3000], 2.0: [50, 2000], 3.0: [25, 1000]
  };
  // Band sizes as published by the manufacturers (NEB product pages N3232,
  // N3200, N3231; Thermo manuals MAN0013004, MAN0013047), retrieved 2026-10.
  // The 500/517 bp doublets are listed as two bands.
  E.LADDERS = {
    'NEB 1 kb Plus': [10002, 8001, 6001, 5001, 4001, 3001, 2017, 1517, 1200, 1000, 900, 800, 700, 600, 517, 500, 400, 300, 200, 100],
    'NEB 1 kb': [10002, 8001, 6001, 5001, 4001, 3001, 2000, 1500, 1000, 517, 500],
    'NEB 100 bp': [1517, 1200, 1000, 900, 800, 700, 600, 517, 500, 400, 300, 200, 100],
    'GeneRuler 1 kb': [10000, 8000, 6000, 5000, 4000, 3500, 3000, 2500, 2000, 1500, 1000, 750, 500, 250],
    'GeneRuler 1 kb Plus': [20000, 10000, 7000, 5000, 4000, 3000, 2000, 1500, 1000, 700, 500, 400, 300, 200, 75],
    'λ DNA HindIII': [23130, 9416, 6557, 4361, 2322, 2027, 564, 125]
  };
  // returns 0 (top of gel, well) .. 1 (bottom)
  E.gelMigration = function (bp, gelPercent) {
    const range = E.GEL_RANGES[gelPercent] || E.GEL_RANGES[1.0];
    const lo = range[0], hi = range[1];
    const x = Math.log10(Math.max(bp, 1));
    const xHi = Math.log10(hi), xLo = Math.log10(lo);
    // large fragments barely move (near 0.05), small fragments run far (near 0.95)
    let t = (xHi - x) / (xHi - xLo);
    t = Math.max(0, Math.min(1, t));
    // compress beyond the resolving window so out-of-range bands pile up
    if (bp > hi) return 0.04;
    if (bp < lo) return 0.96;
    return 0.08 + t * 0.84;
  };


  /* ======================= METHYLATION SENSITIVITY =======================
     Data: NEB "Dam-Dcm and CpG Methylation" selection chart (retrieved
     2026-10). Codes: B blocked, I impaired, Bo/Io blocked/impaired when the
     methylated base overlaps the site, Bs/Is blocked/impaired by SOME
     overlapping combinations (cannot be decided from sequence alone).
     Enzymes not in the chart (non-NEB) are reported as "unknown".
     Methylated bases (top-strand coordinates):
       Dam  G[A]TC   (N6-mA, both strands: A at +1 and the bottom A at +2)
       Dcm  C[C]WGG  (C5-mC, both strands: +1 and the bottom C at +3)
       CpG  [C]G     (M.SssI / mammalian, both strands: +0 and +1)       */
  const METH_RAW = {
    dam: 'AlwI:B BcgI:Io BclI:B BsaBI:Bo BspDI:Bo BspEI:Bo BspHI:Io ClaI:Bo DpnII:B HphI:B Hpy188I:Bo Hpy188III:Bo MboI:B MboII:Bo NruI:Bo Nt.AlwI:B TaqI:Bo XbaI:Bo',
    dcm: 'Acc65I:Bs AlwNI:Bo ApaI:Bo AvaII:Bo BanI:Bs BsaHI:Bs BsaI:Is BslI:Bs BsmFI:Bo BstXI:Bs EaeI:Bo FokI:Io FseI:Is HphI:B MscI:Bo NlaIV:Bo PflMI:Bo PpuMI:Bo PspGI:B PspOMI:Is Sau96I:Bo ScrFI:Bo SexAI:B SfiI:Io SfoI:Bs StuI:Bo StyD4I:Bo',
    cpg: 'AatII:B AccI:Bo Acc65I:Bs AciI:B AclI:B AfeI:B AgeI:B AhdI:Is AleI:Io ApaI:Bo ApaLI:Bo ApeKI:Bo AscI:B AsiSI:B AvaI:B AvaII:Bo BaeI:Bs BanI:Bs BbvCI:Io BceAI:B BcgI:Bs BcoDI:Is BfuAI:Io BglI:Bs BmgBI:B BsaAI:B BsaBI:Bs BsaHI:B BsaI:Bs BseYI:Bo BsiEI:B BsiWI:B BslI:Bs BsmAI:Bs BsmBI:B BsmFI:Bo BspDI:B BspEI:I BsrBI:Bs BsrFI:B BssHII:B BstAPI:Bs BstBI:B BstUI:B BstZ17I:Bs BtgZI:I Cac8I:Bs ClaI:B DraIII:Io DrdI:Bs EaeI:Bo EagI:B EarI:Io EciI:Bs Eco53kI:Bs EcoRI:Bs EcoRV:Is Esp3I:B FauI:B Fnu4HI:Bo FokI:Io FseI:B FspI:B HaeII:B HgaI:B HhaI:B HincII:Bs HinP1I:B HinfI:Bs HpaI:Bs HpaII:B HpyAV:Io HpyCH4IV:B Hpy99I:B Hpy166II:Bo Hpy188III:Bo KasI:B MluI:B MmeI:Bo MspA1I:Bo MwoI:Bs NaeI:B NarI:B NciI:Io NgoMIV:B NheI:Bs NlaIV:Bo NotI:B NruI:B Nt.BbvCI:Bs Nt.BsmAI:B Nt.CviPII:B PaeR7I:B PaqCI:Io PleI:Bs PluTI:B PmeI:Bs PmlI:B PshAI:Bs PspOMI:Bo PspXI:I PvuI:B RsaI:Bs RsrII:B SacI:Bs SacII:B SalI:B Sau3AI:Bo Sau96I:Bo ScrFI:Bo SfaNI:Is SfiI:Bs SfoI:B SgrAI:B SmaI:B SnaBI:B SrfI:B StyD4I:Io TfiI:Bs TseI:Bs TspMI:B XhoI:I XmaI:I ZraI:B'
  };
  E.METHYLATION = {};
  Object.keys(METH_RAW).forEach(kind => METH_RAW[kind].split(' ').forEach(tok => {
    const [n, code] = tok.split(':');
    (E.METHYLATION[n] = E.METHYLATION[n] || {})[kind] = code;
  }));
  E.METHYLATION_SOURCE = 'NEB Dam-Dcm and CpG Methylation chart';

  // Every methylated top-strand coordinate for the chosen systems.
  E.methylatedBases = function (seq, circular, systems) {
    seq = seq.toUpperCase();
    const L = seq.length, hay = circular ? seq + seq.slice(0, 4) : seq;
    const out = { dam: new Set(), dcm: new Set(), cpg: new Set() };
    for (let i = 0; i < L; i++) {
      if (systems.dam && hay.substr(i, 4) === 'GATC') { out.dam.add((i + 1) % L); out.dam.add((i + 2) % L); }
      if (systems.dcm && hay[i] === 'C' && hay[i + 1] === 'C' && (hay[i + 2] === 'A' || hay[i + 2] === 'T') &&
          hay[i + 3] === 'G' && hay[i + 4] === 'G') { out.dcm.add((i + 1) % L); out.dcm.add((i + 3) % L); }
      if (systems.cpg && hay[i] === 'C' && hay[i + 1] === 'G') { out.cpg.add(i % L); out.cpg.add((i + 1) % L); }
    }
    return out;
  };
  // For one recognition site: is a methylated base inside it, and what does
  // the enzyme do about it? Returns null when the site is unaffected.
  E.methylationEffect = function (enz, site, L, methyl) {
    const data = E.METHYLATION[enz.name] || E.METHYLATION[enz.name.replace(/-HF.*$|-v2$/i, '')];
    const n = enz.site.length, hits = [];
    ['dam', 'dcm', 'cpg'].forEach(kind => {
      const set = methyl[kind];
      if (!set || !set.size) return;
      let overlap = false;
      for (let j = 0; j < n && !overlap; j++) if (set.has((site.pos + j) % L)) overlap = true;
      if (!overlap) return;
      const code = data && data[kind];
      // an NEB enzyme absent from the chart's sensitive lists is not sensitive
      if (!data && /N/.test(enz.suppliers || '')) return;
      if (!data) hits.push({ kind, effect: 'unknown' });
      else if (code) hits.push({ kind, effect: code[0] === 'B' ? (code.length > 1 && code[1] === 's' ? 'maybe-blocked' : 'blocked')
                                                         : (code.length > 1 && code[1] === 's' ? 'maybe-impaired' : 'impaired') });
    });
    if (!hits.length) return null;
    const rank = { blocked: 4, 'maybe-blocked': 3, impaired: 2, 'maybe-impaired': 1, unknown: 0 };
    hits.sort((a, b) => rank[b.effect] - rank[a.effect]);
    return { effect: hits[0].effect, by: hits.map(h => h.kind), hits };
  };

  /* ======================= SEQUENCE UTILITIES ======================= */
  // Map a position through an edit so features can follow an insertion/deletion
  /* ======================= FEATURE GEOMETRY UNDER EDITS =======================
     Every transform works on a feature's pieces (top-strand order, never
     wrapping) and then reassembles it, so joined features keep their parts
     and a feature that ends up straddling the origin becomes a wrapped one.
  ============================================================================= */
  function spanOf(parts) {
    if (parts.length === 1) return { start: parts[0].start, end: parts[0].end };
    for (let i = 1; i < parts.length; i++) {
      if (parts[i].start < parts[i - 1].start) {
        return { start: parts[0].start, end: parts[parts.length - 1].end };
      }
    }
    return { start: Math.min.apply(null, parts.map(p => p.start)),
             end: Math.max.apply(null, parts.map(p => p.end)) };
  }
  E.assembleFeature = function (f, pieces, L) {
    const g = Object.assign({}, f);
    delete g.parts; delete g.mixed;
    if (f.mixed) {
      const parts = pieces.filter(p => p.end > p.start);
      if (!parts.length) return null;
      g.parts = parts; g.mixed = true;
      g.start = Math.min.apply(null, parts.map(p => p.start));
      g.end = Math.max.apply(null, parts.map(p => p.end));
      g.strand = parts[0].strand;
      return g;
    }
    const merged = [];
    for (const p of pieces) {
      if (!(p.end > p.start)) continue;
      const last = merged[merged.length - 1];
      if (last && !last.wrapped && last.end === p.start) { last.end = p.end; continue; }
      if (last && !last.wrapped && last.end === L && p.start === 0) { last.end = p.end; last.wrapped = true; continue; }
      merged.push({ start: p.start, end: p.end });
    }
    if (!merged.length) return null;
    const parts = merged.map(p => ({ start: p.start, end: p.end }));
    if (parts.length === 1) { g.start = parts[0].start; g.end = parts[0].end; return g; }
    const sp = spanOf(parts);
    g.start = sp.start; g.end = sp.end; g.parts = parts;
    return g;
  };

  // replace bases [a, b) with n new bases; returns the moved features
  E.editFeatures = function (features, L, a, b, n) {
    const d = n - (b - a);
    const newL = L + d;
    const mapStart = x => x < a ? x : (x >= b ? x + d : a + n);
    const mapEnd = x => x <= a ? x : (x >= b ? x + d : a);
    return (features || []).map(f => {
      if (f.between) {
        const x = mapStart(f.start);
        return x >= 0 && x <= newL ? Object.assign({}, f, { start: x, end: x }) : null;
      }
      const src = f.mixed ? f.parts : E.featurePieces(f, L);
      const moved = src.map(p => Object.assign({}, p, { start: mapStart(p.start), end: mapEnd(p.end) }));
      return E.assembleFeature(f, moved, newL);
    }).filter(Boolean);
  };
  /* ======================= SEARCH =======================
     Finding a pattern in a sequence: IUPAC codes in the query match the
     bases they stand for (N matches anything, R matches A or G), both
     strands are searched, and on a circle a hit may run across base 1.
     Hits are half-open top-strand ranges, with `end` smaller than `start`
     when the hit wraps, the same convention features and primer sites use.
  ==================================================================== */
  E.findPattern = function (seq, pattern, opts) {
    opts = opts || {};
    const s = E.cleanSeq(seq).replace(/U/g, 'T');
    const p = E.cleanSeq(pattern).replace(/U/g, 'T');
    const L = s.length, n = p.length;
    const out = [];
    if (!n || !L || n > L) return out;
    const circ = !!opts.circular;
    const limit = opts.limit || 2000;
    const strands = opts.strand === 1 ? [1] : opts.strand === -1 ? [-1] : [1, -1];
    /* a palindromic query would otherwise be reported twice at one place */
    const rcP = E.revcomp(p);
    const selfRC = rcP === p;
    strands.forEach(strand => {
      if (strand === -1 && selfRC) return;
      const q = strand === 1 ? p : rcP;
      const last = circ ? L - 1 : L - n;
      for (let i = 0; i <= last && out.length < limit; i++) {
        let hit = true;
        for (let k = 0; k < n; k++) {
          const b = s[circ ? (i + k) % L : i + k];
          if (!E.iupacMatch(q[k], b)) { hit = false; break; }
        }
        if (hit) out.push({ start: i, end: (i + n) <= L ? i + n : (i + n) % L,
                            strand: strand, length: n, wraps: i + n > L });
      }
    });
    return out.sort((a, b) => a.start - b.start || a.strand - b.strand);
  };

  /* A peptide anywhere in the six reading frames. The hit is reported as
     the nucleotides that encode it, so it can be shown on the map like any
     other range. Ambiguous bases translate to X and match nothing. */
  E.findPeptide = function (seq, peptide, opts) {
    opts = opts || {};
    const table = opts.table || 1;
    const pep = String(peptide || '').toUpperCase().replace(/[^A-Z*]/g, '');
    const s = E.cleanSeq(seq).replace(/U/g, 'T');
    const L = s.length, out = [];
    if (pep.length < 2 || !L) return out;
    const circ = !!opts.circular;
    const limit = opts.limit || 500;
    [1, -1].forEach(dir => {
      const t = dir === 1 ? s : E.revcomp(s);
      const ext = circ ? t + t.slice(0, Math.min(t.length, pep.length * 3 + 3)) : t;
      for (let f = 0; f < 3; f++) {
        const nt = ext.slice(f, f + Math.floor((ext.length - f) / 3) * 3);
        let aa = '';
        try { aa = E.translate(nt, { table: table }); } catch (e) { continue; }
        let at = aa.indexOf(pep);
        while (at !== -1 && out.length < limit) {
          const a0 = f + at * 3;                       // on t
          if (a0 < L) {
            const a1 = a0 + pep.length * 3;
            /* a range [a0,a1) on the reverse strand is [L-a1, L-a0) on the top one */
            const start = dir === 1 ? a0 % L : ((L - a1) % L + L) % L;
            const end = (start + pep.length * 3) % L;
            out.push({ start: start, end: (start + pep.length * 3) <= L ? start + pep.length * 3 : end,
                       strand: dir, length: pep.length * 3, frame: f,
                       wraps: start + pep.length * 3 > L, peptide: pep });
          }
          at = aa.indexOf(pep, at + 1);
        }
      }
    });
    return out.sort((a, b) => a.start - b.start || a.strand - b.strand);
  };

  /* What a typed query means. A query is read as whatever it can be, and
     every reading that finds something is reported, so "GGATCC" finds both
     the bases and BamHI without the user having to say which they meant. */
  E.classifyQuery = function (q, opts) {
    opts = opts || {};
    const raw = String(q || '').trim();
    const out = { raw: raw, kinds: [] };
    if (!raw) return out;
    const up = raw.toUpperCase();
    const asDNA = E.cleanSeq(up).replace(/U/g, 'T');
    const letters = up.replace(/[\s-]/g, '');
    /* Bases: every character is an IUPAC code, and the query is mostly
       plain ACGT. Without that last test a gene name like kanR reads as
       bases (K, A, N and R are all codes) and a search for it would return
       noise from every matching stretch of the library. */
    const acgt = (asDNA.match(/[ACGT]/g) || []).length;
    if (asDNA.length >= 3 && asDNA.length === letters.length &&
        (acgt === asDNA.length || acgt / asDNA.length >= 0.7)) out.kinds.push('dna');
    /* protein: letters that cannot be a nucleotide code, or told so */
    const aaOnly = up.replace(/[\s-]/g, '');
    if (aaOnly.length >= 3 && /^[ACDEFGHIKLMNPQRSTVWY*]+$/.test(aaOnly) && /[EFILPQJOZX*]/.test(aaOnly))
      out.kinds.push('protein');
    /* an enzyme by name, or by typing its recognition site */
    if (opts.enzymes && opts.enzymes.some(e => e.name.toLowerCase() === raw.toLowerCase() ||
        (asDNA.length >= 4 && e.site === asDNA))) out.kinds.push('enzyme');
    /* any query at all can match a name or a note */
    if (raw.length >= 2) out.kinds.push('text');
    return out;
  };

  /* Search every construct at once. `records` are the constructs,
     `extras.primers` the primer library, `extras.enzymes` the enzymes
     available. Returns a flat list of hits, each naming the construct it
     is in and carrying a range that can be selected and shown. */
  E.searchLibrary = function (records, query, opts) {
    opts = opts || {};
    const cls = E.classifyQuery(query, { enzymes: opts.enzymes });
    const want = opts.kinds && opts.kinds.length ? opts.kinds : cls.kinds;
    const q = cls.raw, qlow = q.toLowerCase();
    const hits = [];
    const perRecord = opts.perRecord || 200;
    const table = opts.table || 1;
    (records || []).forEach(rec => {
      if (!rec || !rec.seq) return;
      const circ = rec.topology === 'circular';
      let n = 0;
      const add = h => { if (n < perRecord) { hits.push(Object.assign({ recId: rec.id, recName: rec.name }, h)); n++; } };
      if (want.indexOf('dna') >= 0) {
        E.findPattern(rec.seq, q, { circular: circ, limit: perRecord }).forEach(h =>
          add({ kind: 'dna', start: h.start, end: h.end, strand: h.strand, wraps: h.wraps,
                what: 'bases ' + E.cleanSeq(q.toUpperCase()) + (h.strand === -1 ? ' (reverse strand)' : '') }));
      }
      if (want.indexOf('protein') >= 0) {
        E.findPeptide(rec.seq, q, { circular: circ, table: table, limit: perRecord }).forEach(h =>
          add({ kind: 'protein', start: h.start, end: h.end, strand: h.strand, wraps: h.wraps,
                what: 'peptide ' + h.peptide + ' (frame ' + (h.strand === 1 ? '+' : '−') + (h.frame + 1) + ')' }));
      }
      if (want.indexOf('enzyme') >= 0 && opts.enzymes) {
        const qdna = E.cleanSeq(q.toUpperCase()).replace(/U/g, 'T');
        const enz = opts.enzymes.find(e => e.name.toLowerCase() === qlow) ||
                    opts.enzymes.find(e => qdna.length >= 4 && e.site === qdna);
        if (enz) {
          let sites = [];
          try { sites = E.findSites(rec.seq, enz, circ); } catch (e) { sites = []; }
          sites.forEach(s => add({ kind: 'enzyme', start: s.pos, end: (s.pos + enz.site.length) % (rec.seq.length || 1) || rec.seq.length,
                                   strand: s.strand || 1, what: enz.name + ' site, cuts at ' + (s.cutTop + 1) }));
        }
      }
      if (want.indexOf('text') >= 0) {
        (rec.features || []).forEach((f, i) => {
          const fields = [f.label, f.type, f.note, f.quals && f.quals.gene, f.quals && f.quals.product];
          if (fields.some(x => x && String(x).toLowerCase().indexOf(qlow) >= 0))
            add({ kind: 'feature', start: f.start, end: f.end, strand: f.strand || 1, featureIndex: i,
                  what: (f.label || f.type) + ' · ' + f.type });
        });
        if (String(rec.name || '').toLowerCase().indexOf(qlow) >= 0)
          add({ kind: 'name', start: 0, end: Math.min(rec.seq.length, 1), what: 'construct name' });
      }
    });
    /* primers are not inside any one construct */
    const primerHits = [];
    if (want.indexOf('text') >= 0 || want.indexOf('dna') >= 0) {
      (opts.primers || []).forEach((p, i) => {
        const byName = want.indexOf('text') >= 0 && String(p.name || '').toLowerCase().indexOf(qlow) >= 0;
        const bySeq = want.indexOf('dna') >= 0 && E.cleanSeq(q.toUpperCase()).length >= 3 &&
          (p.seq.toUpperCase().indexOf(E.cleanSeq(q.toUpperCase())) >= 0 ||
           p.seq.toUpperCase().indexOf(E.revcomp(E.cleanSeq(q.toUpperCase()))) >= 0);
        if (byName || bySeq) primerHits.push({ kind: 'primer', primerIndex: i, name: p.name, seq: p.seq,
          what: byName ? 'primer name' : 'inside the primer sequence' });
      });
    }
    return { query: q, kinds: want, hits: hits, primers: primerHits,
             total: hits.length + primerHits.length };
  };

  /* Where a primer anneals on a template, on either strand.
     A polymerase extends only from a 3' end that is base-paired, so a site
     is seeded on the 3'-most `seed` bases matching exactly (default 8, and
     the last `clamp` bases, default 5, may never mismatch). From there the
     primer is laid along the template toward its 5' end without gaps; the
     annealed part is the longest stretch from the 3' end in which no more
     than `maxMismatchFrac` of bases mismatch (default 10%), ending inside
     a block of at least `minBlock` consecutive matching bases (default 3) —
     one or two bases that happen to match beyond a mismatch do not anneal
     on their own. Whatever lies 5' of that is a tail — a Gibson arm, an
     added restriction site — which is reported as a tail, not as
     mismatches. A site is kept when at least `minAnneal` bases anneal
     (default 15).
     Coordinates are the annealed footprint on the top strand, 0-based and
     half-open; on a circle `end` may be smaller than `start` when the
     footprint runs across base 1, exactly as for a feature. */
  E.primerSites = function (template, primer, opts) {
    opts = opts || {};
    const seed = opts.seed || 8, clamp = Math.min(opts.clamp == null ? 5 : opts.clamp, seed);
    const minAnneal = opts.minAnneal || 15, maxFrac = opts.maxMismatchFrac == null ? 0.10 : opts.maxMismatchFrac;
    const minBlock = opts.minBlock || 3;
    const circ = !!opts.circular;
    const T = E.cleanSeq(template).replace(/U/g, 'T'), P = E.cleanSeq(primer).replace(/U/g, 'T');
    const L = T.length, n = P.length;
    if (n < Math.max(seed, minAnneal) || L < minAnneal) return [];
    const at = i => circ ? T[((i % L) + L) % L] : (i >= 0 && i < L ? T[i] : '');
    const out = [];
    /* the primer as it reads on the top strand: itself for a forward site,
       its reverse complement for a reverse site */
    [1, -1].forEach(strand => {
      const top = strand === 1 ? P : E.revcomp(P);
      /* the seed is the primer's 3' end: the right end of `top` for a
         forward site, the left end for a reverse one */
      const seedSeq = strand === 1 ? top.slice(n - seed) : top.slice(0, seed);
      const hay = circ ? T + T.slice(0, Math.min(L, n)) : T;
      let i = hay.indexOf(seedSeq);
      while (i !== -1) {
        if (i < L || !circ) {
          /* lay the primer along the template: top[k] sits on base off + k */
          const off = strand === 1 ? i - (n - seed) : i;
          /* walk from the 3' end toward the 5' end */
          const order = [];
          for (let k = 0; k < n; k++) order.push(strand === 1 ? n - 1 - k : k);
          let mism = 0, best = 0, run = 0, okClamp = true;
          const mmIdx = [];
          for (let s = 0; s < n; s++) {
            const k = order[s];
            const tb = at(off + k);
            const match = tb !== '' && tb === top[k];
            if (!match) {
              if (s < clamp) { okClamp = false; break; }
              mism++; mmIdx.push(k); run = 0;
            } else run++;
            const len = s + 1;
            if (match && run >= minBlock && mism <= Math.floor(len * maxFrac)) best = len;
            /* once mismatches outrun the allowance there is no way back */
            if (mism > Math.floor(n * maxFrac) + 1 && s >= minAnneal) break;
          }
          if (okClamp && best >= minAnneal) {
            /* the annealed footprint on the top strand */
            let lo = strand === 1 ? off + (n - best) : off;
            let hi = lo + best;
            if (lo < 0) { lo += L; hi += L; }          // only possible on a circle
            if (circ || (lo >= 0 && hi <= L)) {
              /* mismatches inside the annealed part, as 1-based primer
                 positions counted from its 5' end */
              const inAnneal = mmIdx.filter(k => strand === 1 ? k >= n - best : k < best)
                .map(k => strand === 1 ? k + 1 : n - k).sort((a, b) => a - b);
              out.push({ strand: strand, start: lo, end: hi <= L ? hi : hi - L,
                         annealed: best, tail: n - best, mismatches: inAnneal,
                         wraps: hi > L,
                         threePrime: strand === 1 ? (hi - 1) % L : lo,
                         exact: inAnneal.length === 0 && best === n });
            }
          }
        }
        i = hay.indexOf(seedSeq, i + 1);
      }
    });
    /* one report per footprint */
    const seen = new Set();
    return out.filter(h => { const k = h.strand + ':' + h.start + ':' + h.end; if (seen.has(k)) return false; seen.add(k); return true; })
      .sort((a, b) => a.start - b.start || a.strand - b.strand);
  };

  /* Replace bases a..b (0-based, half-open) with `ins`, moving every
     feature with the edit. On a circular sequence b may run past the end
     (b > length), meaning the range carries on across base 1; the edit is
     then made as a deletion of the part after the origin followed by a
     replacement of the part before it, and the new bases close the circle
     at the end of the sequence. Insertion is a = b. */
  /* On a circle, base 1 is not a real boundary: a feature that spans it is
     stored as two pieces split there, and an edit that touches that split
     cannot tell whether new bases belong inside the feature. So a circular
     edit is made in the middle of the sequence — rotate until the edit sits
     strictly inside, away from base 1, edit, and rotate back — where the
     ordinary rule applies unambiguously: inserted bases join a feature
     exactly when the feature runs right through the edited range. */
  E.editSequence = function (seq, features, a, b, ins, circular) {
    const s = E.cleanSeq(seq), L = s.length;
    ins = E.cleanSeq(ins || '');
    const n = ins.length, d = n - (b - a);
    if (!(a >= 0 && a <= L)) throw new Error('Edit starts outside the sequence.');
    if (b < a) throw new Error('Edit ends before it starts.');
    if (!circular) {
      if (b > L) throw new Error('Only a circular sequence can be edited across base 1.');
      return { seq: s.slice(0, a) + ins + s.slice(b),
               features: E.editFeatures(features, L, a, b, n), wrapped: false, shift: 0 };
    }
    if (b - a >= L) throw new Error('That range covers the whole circle.');
    const wraps = b > L, head = wraps ? b - L : 0;
    /* the result, written out directly: the numbering is kept from base 1
       when base 1 survives; when the edit removed it, the circle now starts
       at the first base after the edit */
    const target = wraps ? s.slice(head, a) + ins : s.slice(0, a) + ins + s.slice(b);
    const compLen = L - (b - a);          // bases the edit leaves alone
    if (compLen < 2) {
      /* almost nothing survives: there is no interior to move the edit to,
         and only a feature inside those last bases could remain */
      return { seq: target, features: wraps ? E.editFeatures(E.editFeatures(features, L, 0, head, 0),
               L - head, a - head, L - head, n) : E.editFeatures(features, L, a, b, n),
               wrapped: wraps, shift: head };
    }
    const compStart = b % L;
    const t = Math.max(1, Math.floor(compLen / 2));
    const o = (compStart + t) % L;        // strictly inside the untouched stretch
    const R = E.setOrigin({ seq: s, features: features, topology: 'circular' }, o);
    const ea = (a - o + L) % L, eb = ea + (b - a);
    const Rs = R.seq.slice(0, ea) + ins + R.seq.slice(eb);
    const Rf = E.editFeatures(R.features, L, ea, eb, n);
    const L2 = L + d;
    let k;
    if (wraps) k = ea + n;                                // first base after the edit
    else if (a === 0) k = ea;                             // the new bases are base 1
    else { const p0 = (L - o) % L; k = p0 < ea ? p0 : p0 + d; }   // old base 1
    const back = L2 ? E.setOrigin({ seq: Rs, features: Rf, topology: 'circular' }, k % L2)
                    : { seq: '', features: [] };
    if (back.seq !== target) throw new Error('Internal: circular edit did not land where expected.');
    return { seq: back.seq, features: back.features, wrapped: wraps, shift: head };
  };

  // kept for older callers: insertion (delta > 0) at `at`, or deletion of
  // -delta bases that ended at `at`
  E.shiftFeatures = function (features, at, delta, L) {
    L = L || Math.max.apply(null, [0].concat((features || []).map(f => f.end)));
    return delta >= 0 ? E.editFeatures(features, L, at, at, delta)
                      : E.editFeatures(features, L, at + delta, at, 0);
  };

  /* ======================= ABIF (.ab1) CHROMATOGRAM READER =======================
     Applied Biosystems ABIF container: a 4-byte magic, a root directory entry
     at offset 6, then fixed 28-byte directory entries. Fields of 4 bytes or
     fewer are stored inline in the dataoffset slot instead of being pointed to.
     We pull the base calls (PBAS), quality (PCON), peak locations (PLOC),
     channel order (FWO_) and the four processed traces (DATA 9-12).
  ============================================================================= */
  /* ======================= SNAPGENE (.dna) READER =======================
     SnapGene files are a flat stream of packets: one type byte, a big-endian
     uint32 length, then that many bytes. The packets we need are 0x09 (header,
     which must begin "SnapGene"), 0x00 (the sequence, whose first byte carries
     the circular flag), 0x0A (features, as XML), 0x05 (primers) and 0x06
     (notes). Feature ranges are 1-based inclusive and directionality 2 means
     the reverse strand.
  ==================================================================== */
  function sgText(bytes) {
    if (typeof TextDecoder !== 'undefined') {
      try { return new TextDecoder('utf-8').decode(bytes); } catch (e) { /* fall through */ }
    }
    let s = '';
    for (let i = 0; i < bytes.length; i++) s += String.fromCharCode(bytes[i]);
    return s;
  }
  function xmlDecode(s) {
    return String(s == null ? '' : s)
      .replace(/&lt;/g, '<').replace(/&gt;/g, '>').replace(/&quot;/g, '"')
      .replace(/&apos;/g, "'").replace(/&#x([0-9a-f]+);/gi, (m, h) => String.fromCharCode(parseInt(h, 16)))
      .replace(/&#(\d+);/g, (m, d) => String.fromCharCode(+d))
      .replace(/&amp;/g, '&');
  }
  function xmlAttr(attrs, name) {
    const m = new RegExp('\\b' + name + '\\s*=\\s*"([^"]*)"', 'i').exec(attrs || '');
    return m ? m[1] : null;
  }
  /* SnapGene coordinates. Feature segment ranges are 1-based inclusive
     ("101-961"); a range whose start is at or past its end wraps the origin
     of a circular sequence ("3925-2"). Primer binding-site locations are
     0-based inclusive (verified against real files). Segments typed "gap" are
     spacers inside a feature and are not part of it. Each part is returned as
     a 0-based half-open {start, end}; a wrapped part has end <= start.      */
  function sgRange(range, zeroBased) {
    const parts = String(range).split(/[-–]/).map(Number);
    if (parts.length !== 2 || !isFinite(parts[0]) || !isFinite(parts[1])) return null;
    let start = zeroBased ? parts[0] : parts[0] - 1;
    let end = zeroBased ? parts[1] + 1 : parts[1];
    return { start: start, end: end };
  }
  function parseSGFeatures(xml) {
    const out = [];
    const featRe = /<Feature\b([^>]*?)(\/)?>([\s\S]*?)(?:<\/Feature>|(?=<Feature\b)|$)/gi;
    let m;
    while ((m = featRe.exec(xml))) {
      const attrs = m[1] || '';
      const inner = m[2] ? '' : (m[3] || '');
      const name = xmlAttr(attrs, 'name');
      const type = xmlAttr(attrs, 'type') || 'misc_feature';
      const dir = xmlAttr(attrs, 'directionality');
      const parts = [];
      let color = null;
      const segRe = /<Segment\b([^>]*)\/?>/gi;
      let sm;
      while ((sm = segRe.exec(inner))) {
        if ((xmlAttr(sm[1], 'type') || '') === 'gap') continue;
        const r = sgRange(xmlAttr(sm[1], 'range'), false);
        if (!r) continue;
        if (!color) {
          const c = xmlAttr(sm[1], 'color');
          if (c && /^#[0-9a-f]{3,8}$/i.test(c)) color = c;
        }
        parts.push(r);
      }
      if (!parts.length) continue;
      let label = name;
      const q = /<Q\s+name="label"[^>]*>\s*<V\s+text="([^"]*)"/i.exec(inner);
      if (q) label = q[1];
      // overall span: a wrapped part defines it; otherwise min start..max end
      const wrapped = parts.find(r => r.end <= r.start);
      const start = wrapped ? wrapped.start : Math.min.apply(null, parts.map(r => r.start));
      const end = wrapped ? wrapped.end : Math.max.apply(null, parts.map(r => r.end));
      const f = {
        type: type, start: start, end: end,
        strand: dir === '2' ? -1 : 1,
        label: xmlDecode(label || type),
        color: color, segments: parts.length
      };
      if (parts.length > 1) f.parts = parts.slice().sort((a, b) => a.start - b.start);
      out.push(f);
    }
    return out;
  }
  function parseSGPrimers(xml) {
    const out = [];
    // SnapGene keeps weak off-target binding sites in the file but only shows
    // those meeting the file's hybridisation thresholds; apply the same rule.
    let minLen = 0, minTm = 0;
    const hp = /<HybridizationParams\b([^>]*)\/?>/i.exec(xml);
    if (hp) {
      minLen = parseInt(xmlAttr(hp[1], 'minContinuousMatchLen') || '0', 10) || 0;
      minTm = parseInt(xmlAttr(hp[1], 'minMeltingTemperature') || '0', 10) || 0;
    }
    const re = /<Primer\b([^>]*?)(\/>|>([\s\S]*?)<\/Primer>)/gi;
    let m;
    while ((m = re.exec(xml))) {
      const seq = xmlAttr(m[1], 'sequence');
      if (!seq) continue;
      const sites = [];
      const seen = new Set();
      const bsRe = /<BindingSite\b([^>]*)>/gi;
      let b;
      while ((b = bsRe.exec(m[3] || ''))) {
        const r = sgRange(xmlAttr(b[1], 'location'), true);
        if (!r) continue;
        const strand = (xmlAttr(b[1], 'boundStrand') || '0') === '1' ? -1 : 1;
        const annealed = xmlAttr(b[1], 'annealedBases');
        if (annealed != null && annealed.length < minLen) continue;
        const tmAttr = xmlAttr(b[1], 'meltingTemperature');
        if (tmAttr != null && parseInt(tmAttr, 10) < minTm) continue;
        const key = r.start + ':' + r.end + ':' + strand;
        if (seen.has(key)) continue;        // "simplified" duplicates
        seen.add(key);
        sites.push({ start: r.start, end: r.end, strand: strand });
      }
      out.push({ name: xmlDecode(xmlAttr(m[1], 'name') || 'primer'), seq: E.cleanSeq(seq), sites: sites });
    }
    return out;
  }

  E.parseSnapGene = function (buffer) {
    const dv = new DataView(buffer);
    const u8 = new Uint8Array(buffer);
    if (dv.byteLength < 14) throw new Error('That file is too small to be a SnapGene .dna file.');
    let pos = 0, seq = '', circular = false, features = [], primers = [], notes = '';
    let sawHeader = false, guard = 0;
    while (pos + 5 <= dv.byteLength && guard++ < 10000) {
      const type = dv.getUint8(pos);
      const len = dv.getUint32(pos + 1);
      pos += 5;
      if (len > dv.byteLength - pos) break;
      const slice = u8.subarray(pos, pos + len);
      if (type === 9) {
        if (sgText(slice.subarray(0, 8)) === 'SnapGene') sawHeader = true;
      } else if (type === 0) {
        if (len > 0) { circular = (slice[0] & 1) === 1; seq = sgText(slice.subarray(1)); }
      } else if (type === 10) {
        try { features = parseSGFeatures(sgText(slice)); } catch (e) { features = []; }
      } else if (type === 5) {
        try { primers = parseSGPrimers(sgText(slice)); } catch (e) { primers = []; }
      } else if (type === 6) {
        notes = sgText(slice);
      }
      pos += len;
    }
    if (!sawHeader) throw new Error('This is not a SnapGene .dna file — the SnapGene header is missing.');
    const clean = E.cleanSeq(seq);
    if (!clean.length) throw new Error('No DNA sequence found inside that .dna file.');
    let name = null;
    const nm = /<CustomMapLabel>([^<]*)<\/CustomMapLabel>/i.exec(notes) ||
               /<Description>([^<]*)<\/Description>/i.exec(notes);
    if (nm) name = xmlDecode(nm[1]).trim() || null;
    const primerFeatures = [];
    primers.forEach(pr => (pr.sites || []).forEach(site => {
      if (site.start < 0 || site.start >= clean.length || site.end > clean.length) return;
      primerFeatures.push({ type: 'primer_bind', start: site.start, end: site.end,
                            strand: site.strand, label: pr.name, color: null, segments: 1 });
    }));
    features = features.concat(primerFeatures);
    return {
      seq: clean,
      topology: circular ? 'circular' : 'linear',
      features: features.filter(f => f.start >= 0 && f.start < clean.length &&
                                     f.end > 0 && f.end <= clean.length),
      primers: primers,
      name: name
    };
  };

  const ABIF_ENTRY = 28;
  E.parseAB1 = function (buffer) {
    const dv = new DataView(buffer);
    if (dv.byteLength < 34) throw new Error('File is too small to be an .ab1 chromatogram.');
    const magic = String.fromCharCode(dv.getUint8(0), dv.getUint8(1), dv.getUint8(2), dv.getUint8(3));
    if (magic !== 'ABIF') throw new Error('Not an ABIF (.ab1) file — magic bytes are "' + magic + '".');
    const version = dv.getInt16(4);
    const rootNumElements = dv.getInt32(6 + 12);
    const rootOffset = dv.getInt32(6 + 20);
    const entries = {};
    for (let i = 0; i < rootNumElements; i++) {
      const b = rootOffset + i * ABIF_ENTRY;
      if (b + ABIF_ENTRY > dv.byteLength) break;
      const name = String.fromCharCode(dv.getUint8(b), dv.getUint8(b + 1), dv.getUint8(b + 2), dv.getUint8(b + 3));
      const number = dv.getInt32(b + 4);
      entries[name + '.' + number] = {
        name, number,
        elementType: dv.getInt16(b + 8),
        elementSize: dv.getInt16(b + 10),
        numElements: dv.getInt32(b + 12),
        dataSize: dv.getInt32(b + 16),
        dataOffset: dv.getInt32(b + 20),
        inlineAt: b + 20
      };
    }
    function read(key) {
      const e = entries[key];
      if (!e) return null;
      const base = e.dataSize <= 4 ? e.inlineAt : e.dataOffset;
      if (base < 0 || base + e.dataSize > dv.byteLength) return null;
      const t = e.elementType, n = e.numElements;
      if (t === 2) {   // char / byte
        let s = '';
        for (let i = 0; i < n; i++) s += String.fromCharCode(dv.getUint8(base + i));
        return s;
      }
      if (t === 1) {   // unsigned byte
        const out = new Array(n);
        for (let i = 0; i < n; i++) out[i] = dv.getUint8(base + i);
        return out;
      }
      if (t === 4) {   // short
        const out = new Array(n);
        for (let i = 0; i < n; i++) out[i] = dv.getInt16(base + i * 2);
        return out;
      }
      if (t === 5) {   // long
        const out = new Array(n);
        for (let i = 0; i < n; i++) out[i] = dv.getInt32(base + i * 4);
        return out;
      }
      if (t === 18) {  // pString: leading length byte
        const len = dv.getUint8(base);
        let s = '';
        for (let i = 0; i < len; i++) s += String.fromCharCode(dv.getUint8(base + 1 + i));
        return s;
      }
      if (t === 19) {  // cString: NUL-terminated
        let s = '';
        for (let i = 0; i < n; i++) {
          const c = dv.getUint8(base + i);
          if (c === 0) break;
          s += String.fromCharCode(c);
        }
        return s;
      }
      return null;
    }
    const seq = read('PBAS.2') || read('PBAS.1') || '';
    const qualRaw = read('PCON.2') || read('PCON.1');
    let quality = null;
    if (typeof qualRaw === 'string') {
      quality = [];
      for (let i = 0; i < qualRaw.length; i++) quality.push(qualRaw.charCodeAt(i));
    } else if (Array.isArray(qualRaw)) quality = qualRaw;
    const order = (read('FWO_.1') || 'GATC').replace(/[^ACGT]/g, '') || 'GATC';
    const traces = {};
    const chans = [read('DATA.9'), read('DATA.10'), read('DATA.11'), read('DATA.12')];
    for (let i = 0; i < 4 && i < order.length; i++) {
      if (Array.isArray(chans[i])) traces[order[i]] = chans[i];
    }
    return {
      version,
      sequence: E.cleanSeq(seq.replace(/N/g, 'N')),
      rawSequence: seq,
      quality,
      peakLocations: read('PLOC.2') || read('PLOC.1') || null,
      channelOrder: order,
      traces,
      sampleName: read('SMPL.1') || null,
      runDate: null,
      length: seq.length,
      meanQuality: quality && quality.length
        ? quality.reduce((a, b) => a + b, 0) / quality.length : null,
      hasTraces: Object.keys(traces).length > 0
    };
  };

  /* Trim a read by quality using a sliding window, the way sequencing
     pipelines do: drop leading/trailing bases until the window mean passes. */
  E.trimByQuality = function (seq, quality, opts) {
    opts = opts || {};
    const win = opts.window || 10, min = opts.minMean || 20;
    if (!quality || !quality.length) return { start: 0, end: seq.length, seq, trimmed: 0 };
    const n = Math.min(seq.length, quality.length);
    let start = 0, end = n;
    for (let i = 0; i + win <= n; i++) {
      let s = 0;
      for (let k = 0; k < win; k++) s += quality[i + k];
      if (s / win >= min) { start = i; break; }
    }
    for (let i = n - win; i >= 0; i--) {
      let s = 0;
      for (let k = 0; k < win; k++) s += quality[i + k];
      if (s / win >= min) { end = i + win; break; }
    }
    if (end <= start) { start = 0; end = n; }
    return { start, end, seq: seq.slice(start, end), trimmed: seq.length - (end - start) };
  };

  /* ======================= ALIGN TO REFERENCE =======================
     Place a read (Sanger or any sequence) on a reference and list the
     discrepancies, which is the actual sequencing-verification workflow.
  ==================================================================== */
  E.alignToReference = function (reference, read, opts) {
    opts = opts || {};
    const ref = E.cleanSeq(reference);
    let rd = E.cleanSeq(read);
    let flipped = false;
    // decide orientation by local alignment score
    const fwdScore = E.align(ref, rd, { mode: 'local' }).score;
    const revScore = E.align(ref, E.revcomp(rd), { mode: 'local' }).score;
    if (revScore > fwdScore) { rd = E.revcomp(rd); flipped = true; }
    const aln = E.align(ref, rd, { mode: 'local', match: 2, mismatch: -3, gap: -5 });
    // walk the alignment and record differences with reference coordinates
    const startInRef = findOffset(ref, aln.a.replace(/-/g, ''));
    const diffs = [];
    let refPos = startInRef;
    for (let i = 0; i < aln.a.length; i++) {
      const a = aln.a[i], b = aln.b[i];
      if (a === '-') { diffs.push({ type: 'insertion', refPos, base: b }); continue; }
      if (b === '-') { diffs.push({ type: 'deletion', refPos, base: a }); refPos++; continue; }
      if (a !== b) diffs.push({ type: 'mismatch', refPos, ref: a, read: b });
      refPos++;
    }
    return {
      flipped, identity: aln.identity, alignedLength: aln.length,
      refStart: startInRef, refEnd: refPos,
      coverage: ref.length ? (refPos - startInRef) / ref.length * 100 : 0,
      differences: diffs, alignment: aln,
      clean: diffs.length === 0
    };
  };
  function findOffset(hay, needle) {
    if (!needle) return 0;
    const i = hay.indexOf(needle);
    if (i !== -1) return i;
    // fall back to anchoring on the first 20 bases
    const j = hay.indexOf(needle.slice(0, 20));
    return j === -1 ? 0 : j;
  }

  // Rotate a circular sequence so that `origin` becomes position 0
  /* ======================= CONSENSUS FROM SEVERAL READS =======================
     Align every read to the reference and vote base by base. What comes back
     is the sequence the reads actually support, plus how many reads stand
     behind each position — because a base no read covered is the reference's
     guess, not a result, and a protein built partly from those is not
     verified. Insertions in a read are recorded but not pushed into the
     consensus, which stays in reference coordinates.
  ==================================================================== */
  E.readPileup = function (reference, reads, opts) {
    opts = opts || {};
    const ref = E.cleanSeq(reference), L = ref.length;
    const minIdentity = opts.minIdentity != null ? opts.minIdentity : 80;
    const minLen = opts.minLength || 30;
    const votes = new Array(L), used = [], rejected = [];
    for (const r of (reads || [])) {
      const seq = E.cleanSeq(r.seq || r);
      const name = r.name || 'read';
      if (seq.length < minLen) { rejected.push({ name, why: 'shorter than ' + minLen + ' bases' }); continue; }
      /* An alignment the caller already made against this same reference is
         reused: aligning is by far the expensive part, and the view has
         usually done it once already when the read was loaded. */
      let a = r.aln || null;
      if (!a) { try { a = E.alignToReference(ref, seq, {}); } catch (e) { a = null; } }
      if (!a) { rejected.push({ name, why: 'could not be aligned' }); continue; }
      /* Identity alone means nothing here: a local alignment of an unrelated
         read finds some short chance match and reports 100%. The alignment
         has to account for a real share of the read before it may vote. */
      const need = Math.max(opts.minAligned || 60, Math.round(seq.length * (opts.minFraction || 0.35)));
      if (a.alignedLength < need) {
        rejected.push({ name, why: 'only ' + a.alignedLength + ' of its ' + seq.length +
          ' bases align to this construct — is this the right reference?' });
        continue;
      }
      if (a.identity < minIdentity) { rejected.push({ name, why: 'only ' + a.identity.toFixed(1) + '% identical over the part that aligns' }); continue; }
      used.push({ name, flipped: a.flipped, identity: a.identity, refStart: a.refStart, refEnd: a.refEnd });
      let refPos = a.refStart;
      const A = a.alignment.a, Bq = a.alignment.b;
      for (let i = 0; i < A.length; i++) {
        const x = A[i], y = Bq[i];
        if (x === '-') continue;                       // insertion in the read
        if (refPos >= 0 && refPos < L) {
          const v = votes[refPos] || (votes[refPos] = {});
          const call = (y === '-') ? '-' : y;
          v[call] = (v[call] || 0) + 1;
        }
        refPos++;
      }
    }
    let cons = '';
    const coverage = new Array(L), disagreements = [];
    for (let i = 0; i < L; i++) {
      const v = votes[i];
      if (!v) { cons += ref[i]; coverage[i] = 0; continue; }
      let best = null, total = 0;
      for (const k in v) { total += v[k]; if (best === null || v[k] > v[best]) best = k; }
      coverage[i] = total;
      cons += (best === '-') ? ref[i] : best;
      if (best !== ref[i] || v[best] !== total) {
        disagreements.push({ pos: i, ref: ref[i], call: best, depth: total, supporting: v[best], votes: Object.assign({}, v) });
      }
    }
    return { consensus: cons, coverage: coverage, disagreements: disagreements, reads: used, rejected: rejected };
  };

  /* How much of a stretch the reads actually cover, and whether they agree. */
  E.coverageOver = function (pile, positions) {
    let min = Infinity, once = 0, twice = 0, conflicts = 0;
    const seen = {};
    pile.disagreements.forEach(d => { seen[d.pos] = d; });
    positions.forEach(p => {
      const c = pile.coverage[p] || 0;
      if (c < min) min = c;
      if (c >= 1) once++;
      if (c >= 2) twice++;
      if (seen[p] && seen[p].supporting !== seen[p].depth) conflicts++;
    });
    const n = positions.length || 1;
    return { minDepth: min === Infinity ? 0 : min, covered: once / n * 100, doubleCovered: twice / n * 100, conflicts: conflicts };
  };

  /* ======================= READ TRIAGE =======================
     Deciding which of a read's six frames is the protein, and whether the
     read is worth believing at all. Three independent signals are combined,
     because any one of them alone is wrong on real sequencing data:

       1. Longest stretch with no stop codon AND no ambiguous residue.
          Stops alone are not enough: an N base translates to X, not to a
          stop, so a dead trace reports a deceptively low stop count.

       2. Amino-acid composition against the Swiss-Prot average. A nonsense
          frame is wildly skewed (runs of F, Q, H, C) while a real protein
          sits near the average. This is what separates two frames that are
          the same length, which the stop count cannot do.

       3. Peptide landmarks — affinity tags, protease sites, vector
          junctions, and the coding sequences of the loaded construct. A
          landmark pins the frame outright and outranks the other two.

     Nothing here is a trained model; every threshold is in E.TRIAGE and
     every signal is reported alongside the verdict so a call can be
     checked by hand.
  ============================================================ */

  /* Swiss-Prot average amino-acid composition, in percent.
     UniProtKB/Swiss-Prot release statistics, amino-acid composition table. */
  E.AA_BACKGROUND = { A: 8.25, R: 5.53, N: 4.06, D: 5.45, C: 1.38, Q: 3.93,
                      E: 6.75, G: 7.07, H: 2.27, I: 5.91, L: 9.65, K: 5.81,
                      M: 2.41, F: 3.86, P: 4.74, S: 6.65, T: 5.36, W: 1.10,
                      Y: 2.92, V: 6.86 };
  E.AA_BACKGROUND_SOURCE = 'UniProtKB/Swiss-Prot amino-acid composition';

  const AA_EXP = (function () {
    const t = Object.keys(E.AA_BACKGROUND).reduce((s, k) => s + E.AA_BACKGROUND[k], 0);
    const o = {};
    Object.keys(E.AA_BACKGROUND).forEach(k => { o[k] = E.AA_BACKGROUND[k] / t; });
    return o;
  })();

  /* Peptide landmarks. Each is at least 6 residues, because anything
     shorter turns up by chance in a nonsense frame. */
  E.PEPTIDE_LANDMARKS = [
    { name: 'His6', seq: 'HHHHHH', kind: 'tag' },
    { name: 'His8', seq: 'HHHHHHHH', kind: 'tag' },
    { name: 'His6-TEV (pET)', seq: 'HHHHHHSSGVDLGTENLYFQ', kind: 'tag' },
    { name: 'His6-thrombin (pET-28)', seq: 'SSHHHHHHSSGLVPRGS', kind: 'tag' },
    { name: 'TEV site', seq: 'ENLYFQ', kind: 'protease' },
    { name: 'HRV 3C / PreScission site', seq: 'LEVLFQGP', kind: 'protease' },
    { name: 'thrombin site', seq: 'LVPRGS', kind: 'protease' },
    { name: 'GST (Sj26) C-terminus', seq: 'WPLQGWQATFGGGDHPPK', kind: 'tag' },
    { name: 'MBP linker', seq: 'NNNNNNNNNNLGIEGR', kind: 'tag' },
    { name: 'thioredoxin (TrxA)', seq: 'WAEWCGPCKMIAPILDEIADEY', kind: 'tag' },
    { name: 'FLAG', seq: 'DYKDDDDK', kind: 'tag' },
    { name: 'Strep-tag II', seq: 'WSHPQFEK', kind: 'tag' },
    { name: 'myc', seq: 'EQKLISEEDL', kind: 'tag' },
    { name: 'HA', seq: 'YPYDVPDYA', kind: 'tag' },
    { name: 'V5', seq: 'GKPIPNPLLGLDST', kind: 'tag' },
    { name: 'S-tag', seq: 'KETAAAKFERQHMDS', kind: 'tag' },
    { name: 'T7 tag', seq: 'MASMTGGQQMG', kind: 'tag' },
    { name: 'AviTag', seq: 'GLNDIFEAQKIEWHE', kind: 'tag' },
    { name: 'GFP', seq: 'FSVSGEGEGDAT', kind: 'reporter' },
    { name: 'mCherry', seq: 'MRFKVHMEGSVNGHEFEIEGEGEGRPY', kind: 'reporter' },
    { name: 'GS linker', seq: 'GGGGSGGGGS', kind: 'linker' }
  ];

  /* Thresholds, all in one place so a verdict can be argued with. */
  E.TRIAGE = {
    minLen: 60,          // nt below which a read is not worth translating
    ambigBad: 0.08,      // non-ACGT fraction that condemns a read
    ambigWarn: 0.03,
    cleanMin: 50,        // aa; below this there is no protein to speak of
    cleanGood: 150,      // aa; at or above this a called frame is "clean"
    domMin: 1.5,         // how far the called frame must beat the runner-up
                         // when no landmark pins it
    shiftMin: 40,        // aa of clean peptide in a shifted frame to call an indel
    shiftRatio: 1.5,
    shiftFit: 0.45,     // composition both halves must beat before an indel is claimed
    compWarn: 0.35,     // composition beyond which an unpinned frame is not called protein
    compWeight: 2.2,     // composition penalty strength
    compFloor: 0.25,
    landmarkBonus: 1000, // a landmark outranks any length difference
    pinMinClean: 40,     // aa of clean peptide before a landmark may pin a frame
    pinGap: 8            // aa a landmark may sit away from that stretch
  };

  /* Chi-square per residue against the Swiss-Prot average. Low is
     protein-like. Returns null for peptides too short to judge. */
  E.compositionFit = function (pep) {
    let n = 0; const c = {};
    for (let i = 0; i < pep.length; i++) {
      const a = pep[i];
      if (AA_EXP[a] === undefined) continue;
      c[a] = (c[a] || 0) + 1; n++;
    }
    if (n < 40) return null;
    let chi = 0;
    Object.keys(AA_EXP).forEach(k => {
      const e = n * AA_EXP[k], o = c[k] || 0;
      chi += (o - e) * (o - e) / e;
    });
    return chi / n;
  };

  /* Ambiguity-tolerant landmark search. X matches anything, but a window
     has to carry real residues to count, or an all-X run would match every
     landmark in the list — which is exactly what happens on a dead trace. */
  E.findLandmarks = function (pep, extra, opts) {
    opts = opts || {};
    const list = E.PEPTIDE_LANDMARKS.concat(extra || []);
    const out = [];
    list.forEach(lm => {
      const pat = lm.seq, n = pat.length;
      if (n < 6) return;
      /* Relaxed matching accepts a window that is mostly ambiguous, which
         recovers a six-residue site carrying three N calls. It is only safe
         on a read whose bases are otherwise confident — on a dead trace a
         run of X matches every landmark in the list — so readTriage turns it
         on only below E.TRIAGE.ambigWarn, and it allows no mismatches. */
      const allowed = opts.relaxed ? 0 : Math.floor(n / 8);
      const needSolid = opts.relaxed ? Math.max(3, Math.ceil(n / 2) - 1) : Math.max(4, Math.ceil(n / 2));
      for (let i = 0; i + n <= pep.length; i++) {
        let mm = 0, solid = 0, amb = 0, ok = true;
        for (let k = 0; k < n; k++) {
          const a = pep[i + k];
          if (a === 'X') { amb++; continue; }
          if (a !== pat[k]) { if (++mm > allowed) { ok = false; break; } }
          solid++;
        }
        const ambOK = opts.relaxed ? amb * 3 < n * 2 : amb * 2 < n;
        if (!ok || !ambOK || solid < needSolid) continue;
        out.push({ name: lm.name, kind: lm.kind || 'tag', at: i,
                   mismatches: mm, ambiguous: amb, exact: mm === 0 && amb === 0 });
        break; // one hit per landmark is enough to pin the frame
      }
    });
    return out;
  };

  /* Landmarks taken from a loaded construct: every annotated coding
     sequence, sampled so that a read landing anywhere in it still matches.
     These are long, so they survive a few N calls. */
  E.landmarksFromRecord = function (rec, opts) {
    opts = opts || {};
    const win = opts.window || 14, step = opts.step || 24;
    const out = [];
    (rec && rec.features || []).forEach(f => {
      const t = (f.type || '').toUpperCase();
      if (t !== 'CDS' && t !== 'GENE' && t !== 'ORF') return;
      let aa = f.quals && (f.quals.translation || f.quals.Translation);
      if (!aa) {
        const nt = rec && rec.seq ? E.featureSequence(rec.seq, f) : null;
        if (!nt || nt.length < win * 3) return;
        try { aa = E.translate(nt, { table: opts.table || 1 }); } catch (e) { return; }
      }
      aa = String(aa).replace(/[^A-Z*]/gi, '').toUpperCase().replace(/\*+$/, '');
      if (aa.length < win) return;
      const label = f.label || f.quals && f.quals.gene || 'CDS';
      for (let i = 0; i + win <= aa.length; i += step) {
        out.push({ name: label, seq: aa.substr(i, win), kind: 'construct' });
      }
    });
    return out;
  };

  function openRuns(s, breakers) {
    const out = []; let cur = 0, curAt = 0;
    for (let i = 0; i <= s.length; i++) {
      if (i < s.length && breakers.indexOf(s[i]) < 0) { if (cur === 0) curAt = i; cur++; }
      else { if (cur) out.push({ len: cur, at: curAt }); cur = 0; }
    }
    return out;
  }
  function longestOpenRun(s, breakers) {
    return openRuns(s, breakers).reduce((a, b) => (b.len > a.len ? b : a), { len: 0, at: 0 });
  }

  /* All six frames of a read, with every signal attached. */
  E.frameTranslate = function (seq, opts) {
    opts = opts || {};
    const table = opts.table || 1;
    const extra = opts.landmarks || [];
    const s0 = E.cleanSeq(seq).replace(/U/g, 'T');
    const rows = [];
    [1, -1].forEach(dir => {
      const s = dir === 1 ? s0 : E.revcomp(s0);
      for (let f = 0; f < 3; f++) {
        const nt = s.slice(f, f + Math.floor((s.length - f) / 3) * 3);
        let aa = '';
        try { aa = E.translate(nt, { table: table }); } catch (e) { continue; }
        const pieces = aa.split('*');
        let best = '', bestAt = 0, at = 0;
        pieces.forEach(x => { if (x.length > best.length) { best = x; bestAt = at; } at += x.length + 1; });
        const clean = longestOpenRun(aa, ['*', 'X']);
        const cleanSeq = aa.substr(clean.at, clean.len);
        const fit = E.compositionFit(cleanSeq.length >= 40 ? cleanSeq : best);
        const cw = E.TRIAGE.compWeight;
        const cf = fit === null ? 0.70
          : Math.max(E.TRIAGE.compFloor, Math.min(1, 1 / (1 + cw * fit)));
        /* A landmark only pins the frame when the frame also reads like a
           frame around it. A tag matched inside a stretch of N calls, in a
           frame broken by stops every few residues, is a coincidence, not a
           reading frame — so the match is still reported but it carries no
           weight in the call. */
        const lms = E.findLandmarks(aa, extra, { relaxed: !!opts.relaxed });
        const T = E.TRIAGE;
        const c0 = clean.at, c1 = clean.at + clean.len;
        let pinned = false;
        lms.forEach(l => {
          const n = (extra.concat(E.PEPTIDE_LANDMARKS).find(x => x.name === l.name) || { seq: '' }).seq.length;
          const near = l.at + n >= c0 - T.pinGap && l.at <= c1 + T.pinGap;
          l.pins = clean.len >= T.pinMinClean && near;
          if (l.pins) pinned = true;
        });
        const m = best.indexOf('M');
        let amb = 0;
        for (let i = 0; i < aa.length; i++) if (aa[i] === 'X') amb++;
        /* The bases each peptide was translated from, cut here where the
           strand is in hand. A caller that wanted to recover them itself
           would have to know that a minus frame is read off the reverse
           complement and that the offsets below are measured on THAT
           strand, and getting either wrong yields bases that look
           plausible and translate to something else. Cut here, the
           invariant holds by construction: translate(xxxNt) === xxx, so
           what is copied is the coding strand of that reading frame,
           5'->3', without the stop codon that ended it. */
        const longestAt = f + bestAt * 3, cleanAt = f + clean.at * 3;
        rows.push({
          label: (dir === 1 ? '+' : '−') + (f + 1), dir: dir, frame: f, offset: f,
          full: aa, stops: pieces.length - 1,
          longest: best, longestAt: longestAt,
          longestNt: s.substr(longestAt, best.length * 3),
          clean: cleanSeq, cleanLen: clean.len, cleanAt: cleanAt,
          cleanNt: s.substr(cleanAt, cleanSeq.length * 3),
          ambiguous: amb, ambiguousFrac: aa.length ? amb / aa.length : 0,
          fit: fit, compFactor: cf,
          landmarks: lms, pinned: pinned,
          fromMet: m >= 0 ? best.slice(m) : '',
          fromMetNt: m >= 0 ? s.substr(longestAt + m * 3, (best.length - m) * 3) : '',
          score: clean.len * cf + (pinned ? E.TRIAGE.landmarkBonus : 0)
        });
      }
    });
    return rows;
  };

  /* Does the read carry an indel? If the peptide continues cleanly in a
     different frame of the same strand just past where the called frame
     breaks, one base went missing or got added — the usual Sanger failure
     near the end of a trace, and the reason a real protein can show up
     split across two frames with neither one complete. */
  /* Spans are compared in nucleotides, and either side of the called
     stretch is considered, because a deletion early in a read leaves the
     longer half downstream and the called frame is then the second half. */
  function shiftEvidence(frames, chosen) {
    const T = E.TRIAGE;
    /* Only worth asking once the called frame itself looks like a protein.
       Non-coding DNA has fairly long stop-free runs in several frames at
       once, and calling that a frameshift would be nonsense. */
    if (!chosen.pinned && !(chosen.cleanLen >= T.cleanMin &&
        chosen.fit !== null && chosen.fit < T.shiftFit)) return null;
    const a0 = chosen.cleanAt, a1 = chosen.cleanAt + chosen.cleanLen * 3;
    const outside = (f, r) => {
      const s = f.offset + r.at * 3, e = s + r.len * 3;
      const ov = Math.max(0, Math.min(e, a1) - Math.max(s, a0));
      return ov * 5 < (e - s) * 2 ? { start: s, end: e, len: r.len } : null;
    };
    /* the best stretch the called frame itself has outside its own call —
       a real tail of noise gives almost nothing here */
    let own = 0;
    openRuns(chosen.full, ['*', 'X']).forEach(r => {
      const o = outside(chosen, r); if (o && o.len > own) own = o.len;
    });
    let best = null;
    frames.filter(f => f.dir === chosen.dir && f.frame !== chosen.frame).forEach(f => {
      openRuns(f.full, ['*', 'X']).forEach(r => {
        const o = outside(f, r);
        if (!o || o.len < T.shiftMin) return;
        /* the continuation has to look like protein too */
        const fit = E.compositionFit(f.full.substr(r.at, r.len));
        if (fit !== null && fit >= T.shiftFit) return;
        if (!best || o.len > best.len) best = { len: o.len, label: f.label, start: o.start, end: o.end, fit: fit };
      });
    });
    if (best && best.len >= T.shiftRatio * Math.max(own, 1)) {
      /* Where the frames swap over. When the two stretches overlap, the
         missing or extra base is inside the overlap, because both frames
         still read through it. */
      const lo = Math.max(a0, best.start), hi = Math.min(a1, best.end);
      const from = hi > lo ? lo : (best.start >= a1 ? a1 : best.end);
      const to = hi > lo ? hi : from;
      return { at: Math.round((from + to) / 2), from: from, to: to,
               into: best.label, continues: best.len, sameFrame: own,
               fit: best.fit, before: best.start < a0 };
    }
    return null;
  }

  /* The readable part of a read that came without quality values — a
     pasted sequence or a .seq file. Sanger base calls are N-rich for the
     first few dozen bases and again where the trace dies, and judging the
     whole read on those ends would condemn a good middle. Each end is cut
     back until a 20-base window carries at most one ambiguous base. A read
     with no such window anywhere is left whole, and fails on its own. */
  E.usableSpan = function (seq, opts) {
    opts = opts || {};
    const w = opts.window || 20, maxBad = opts.maxBad == null ? 1 : opts.maxBad;
    const s = E.cleanSeq(seq).replace(/U/g, 'T');
    const L = s.length;
    const isBad = c => 'ACGT'.indexOf(c) < 0;
    if (L < w) return { start: 0, end: L, trimmed: 0 };
    const badIn = new Array(L + 1).fill(0);
    for (let i = 0; i < L; i++) badIn[i + 1] = badIn[i] + (isBad(s[i]) ? 1 : 0);
    const count = (a, b) => badIn[b] - badIn[a];
    let start = -1, end = -1;
    for (let i = 0; i + w <= L; i++) if (count(i, i + w) <= maxBad) { start = i; break; }
    for (let i = L; i - w >= 0; i--) if (count(i - w, i) <= maxBad) { end = i; break; }
    if (start < 0 || end <= start) return { start: 0, end: L, trimmed: 0 };
    /* do not keep a lone ambiguous base at either edge */
    while (start < end && isBad(s[start])) start++;
    while (end > start && isBad(s[end - 1])) end--;
    return { start: start, end: end, trimmed: L - (end - start) };
  };

  /* Reads from pasted text or a text file: FASTA with any number of
     records, a .seq file, or bare bases. A FASTA header becomes the read
     name; without one the file name is used, numbered if there are several. */
  E.parseReadsText = function (text, fallbackName) {
    const t = String(text || '').replace(/\r\n?/g, '\n');
    const out = [];
    if (/^\s*>/.test(t)) {
      let name = null, buf = [];
      const flush = () => {
        if (name === null) return;
        const seq = E.cleanSeq(buf.join('')).replace(/U/g, 'T');
        if (seq.length) out.push({ name: name || (fallbackName || 'read') + '_' + (out.length + 1), seq: seq });
      };
      t.split('\n').forEach(line => {
        if (/^\s*>/.test(line)) {
          flush();
          name = line.replace(/^\s*>\s*/, '').trim().split(/\s+/)[0] || '';
          buf = [];
        } else if (!/^\s*[;#]/.test(line)) buf.push(line);
      });
      flush();
    } else {
      const seq = E.cleanSeq(t).replace(/U/g, 'T');
      if (seq.length) out.push({ name: fallbackName || 'read', seq: seq });
    }
    return out;
  };

  /* The whole verdict on one read. */
  E.readTriage = function (seq, opts) {
    opts = opts || {};
    const T = E.TRIAGE;
    const whole = E.cleanSeq(seq).replace(/U/g, 'T');
    /* judge the readable part, unless told the read is already trimmed */
    const span = opts.trimEnds === false ? { start: 0, end: whole.length, trimmed: 0 }
      : E.usableSpan(whole);
    const s = whole.slice(span.start, span.end);
    let bad = 0;
    for (let i = 0; i < s.length; i++) if ('ACGT'.indexOf(s[i]) < 0) bad++;
    const ambigFrac = s.length ? bad / s.length : 1;
    const base = { length: whole.length, span: span, ambiguous: bad, ambiguousFrac: ambigFrac,
                   frames: [], chosen: null, runnerUp: null, dominance: 0,
                   shift: null, flags: [], verdict: 'too-short' };
    /* a position on either strand of the judged span, as a 1-based base of
       the read as it was given */
    const readPos = (pos, dir) => dir === 1 ? span.start + pos + 1 : span.start + (s.length - pos);
    if (s.length < T.minLen) {
      base.flags.push({ level: 'bad', text: 'Only ' + s.length + ' readable bases' +
        (span.trimmed ? ' once the N-rich ends are set aside' : '') + ' — too short to translate.' });
      return base;
    }
    if (span.trimmed)
      base.flags.push({ level: 'note', text: 'Judged on bases ' + (span.start + 1) + '–' + span.end +
        ' of ' + whole.length + '; the ends outside that are mostly N, as Sanger ends usually are.' });
    const frames = E.frameTranslate(s, Object.assign({}, opts,
      { relaxed: opts.relaxed !== undefined ? opts.relaxed : ambigFrac <= T.ambigWarn }));
    frames.sort((a, b) => b.score - a.score);
    const chosen = frames[0], next = frames[1] || null;
    const pinned = !!chosen.pinned;
    const nextPinned = next && next.pinned;
    const dom = next ? (pinned && !nextPinned ? Infinity : chosen.score / Math.max(next.score, 1)) : Infinity;
    base.frames = frames; base.chosen = chosen; base.runnerUp = next; base.dominance = dom;

    let verdict = null;
    if (ambigFrac > T.ambigBad) verdict = 'low-quality';
    else if (chosen.cleanLen < T.cleanMin) verdict = 'no-protein';
    else if (!pinned && dom < T.domMin) verdict = 'frame-ambiguous';
    else if (!pinned && chosen.fit !== null && chosen.fit >= T.compWarn) verdict = 'open-only';
    /* An indel is only worth claiming on a read that otherwise looks like a
       protein. Non-coding DNA has several half-open frames at once, and
       calling that a frameshift would be noise dressed up as a finding. */
    const proteinLike = pinned || (chosen.fit !== null && chosen.fit < T.compWarn);
    if ((verdict === null || verdict === 'frame-ambiguous') && proteinLike)
      base.shift = shiftEvidence(frames, chosen);
    if (base.shift) {
      /* shiftEvidence works on the strand it was read from; report where it
         is in the read the user actually has in front of them */
      const d = chosen.dir, a = readPos(base.shift.from, d), b = readPos(base.shift.to, d);
      base.shift.readFrom = Math.min(a, b); base.shift.readTo = Math.max(a, b);
      base.shift.readAt = Math.round((base.shift.readFrom + base.shift.readTo) / 2);
    }
    if (verdict === null)
      verdict = base.shift ? 'frameshift' : (chosen.cleanLen < T.cleanGood ? 'partial' : 'clean');
    base.verdict = verdict;

    if (ambigFrac > T.ambigBad)
      base.flags.push({ level: 'bad', text: (ambigFrac * 100).toFixed(0) +
        '% of the bases are ambiguous — the trace failed. A low stop count here means nothing, ' +
        'because an N translates to X rather than to a stop.' });
    else if (ambigFrac > T.ambigWarn)
      base.flags.push({ level: 'warn', text: (ambigFrac * 100).toFixed(1) + '% ambiguous bases.' });
    if (pinned)
      base.flags.push({ level: 'ok', text: 'Frame pinned by ' +
        chosen.landmarks.filter(l => l.pins)
          .map(l => l.name + (l.exact ? '' : ' (partial match)')).join(', ') + '.' });
    const loose = chosen.landmarks.filter(l => !l.pins);
    if (loose.length && !pinned)
      base.flags.push({ level: 'note', text: 'A match to ' + loose.map(l => l.name).join(', ') +
        ' turned up in this frame but sits outside the readable stretch, so it was not used to call the frame.' });
    if (!pinned && dom < T.domMin && next)
      base.flags.push({ level: 'warn', text: 'No frame stands out: ' + chosen.label + ' gives ' +
        chosen.cleanLen + ' aa and ' + next.label + ' gives ' + next.cleanLen + ' aa. ' +
        ((opts.landmarks && opts.landmarks.length)
          ? 'Nothing in the read matches a tag or any coding sequence of the loaded construct, so the frame is a guess.'
          : 'Nothing in the read matches a tag, and no construct is loaded to check it against, so the frame is a guess.') });
    if (base.shift)
      base.flags.push({ level: 'warn', text: 'The peptide carries on in frame ' + base.shift.into +
        ' for another ' + base.shift.continues + ' aa ' +
        (base.shift.readTo > base.shift.readFrom
          ? 'after a switch between bases ' + base.shift.readFrom + ' and ' + base.shift.readTo
          : 'from base ' + base.shift.readAt) +
        ' of the read — one base is probably inserted or deleted there, so neither frame holds ' +
        'the whole protein. Look at the trace at that spot before deciding it is real.' });
    if (chosen.fit !== null && chosen.fit >= T.compWarn)
      base.flags.push({ level: 'warn', text: 'Amino-acid composition is far from any ordinary protein ' +
        '(chi-square per residue ' + chosen.fit.toFixed(2) + ' against the Swiss-Prot average; real frames ' +
        'in a good read sit below 0.3). There is an open stretch here, but it may not be coding at all.' });
    if (verdict === 'clean' && chosen.stops > 3)
      base.flags.push({ level: 'note', text: chosen.stops + ' stops in this frame, all outside the ' +
        chosen.cleanLen + ' aa stretch — normal tail decay past the good sequence.' });
    return base;
  };

  /* short forms, for counting up a set of reads */
  E.VERDICTS_SHORT = {
    'clean': 'clean', 'partial': 'short but protein-like', 'open-only': 'not protein-like',
    'frameshift': 'indel', 'frame-ambiguous': 'frame unclear', 'low-quality': 'trace failed',
    'no-protein': 'no protein', 'too-short': 'too short'
  };
  E.VERDICTS = {
    'clean': 'Clean coding frame',
    'partial': 'Short but protein-like',
    'open-only': 'Open stretch, not protein-like',
    'frameshift': 'Indel splits the frame',
    'frame-ambiguous': 'No frame stands out',
    'low-quality': 'Trace failed',
    'no-protein': 'No protein in any frame',
    'too-short': 'Too short'
  };

  /* Reads that are the same DNA. Nine independent clones cannot give
     byte-identical sequence, so a group here means one template was
     sequenced several times, the wells were mixed up, or the facility
     returned one trace under several names. Canonical k-mers, so a read
     and its reverse complement group together. */
  E.duplicateReads = function (reads, opts) {
    opts = opts || {};
    const k = opts.k || 24, minShare = opts.minShare || 0.8, minLen = opts.minLen || 120;
    const sets = reads.map(r => {
      const s = E.cleanSeq(typeof r === 'string' ? r : (r.seq || '')).replace(/U/g, 'T');
      const core = s.replace(/^[^ACGT]+/, '').replace(/[^ACGT]+$/, '');
      const set = new Set();
      if (core.length >= Math.max(k, minLen)) {
        for (let i = 0; i + k <= core.length; i++) {
          const w = core.substr(i, k);
          if (/[^ACGT]/.test(w)) continue;
          const rc = E.revcomp(w);
          set.add(w < rc ? w : rc);
        }
      }
      return set;
    });
    const n = reads.length, seen = new Array(n).fill(false), groups = [];
    for (let i = 0; i < n; i++) {
      if (seen[i] || !sets[i].size) continue;
      const g = [i];
      for (let j = i + 1; j < n; j++) {
        if (seen[j] || !sets[j].size) continue;
        let shared = 0;
        sets[j].forEach(w => { if (sets[i].has(w)) shared++; });
        const containment = shared / Math.min(sets[i].size, sets[j].size);
        if (containment >= minShare) { g.push(j); seen[j] = true; }
      }
      if (g.length > 1) { g.forEach(x => { seen[x] = true; }); groups.push(g); }
    }
    return groups;
  };

  E.setOrigin = function (rec, origin) {
    const L = rec.seq.length;
    const o = ((origin % L) + L) % L;
    const seq = rec.seq.slice(o) + rec.seq.slice(0, o);
    const mapX = x => ((x - o) % L + L) % L;
    const features = (rec.features || []).map(f => {
      if (f.between) return Object.assign({}, f, { start: mapX(f.start), end: mapX(f.start) });
      const src = f.mixed ? f.parts : E.featurePieces(f, L);
      const out = [];
      for (const p of src) {
        const s0 = mapX(p.start);
        let e0 = mapX(p.end); if (e0 === 0) e0 = L;
        if (e0 > s0) out.push(Object.assign({}, p, { start: s0, end: e0 }));
        else {                                   // this piece now straddles the origin
          const hi = Object.assign({}, p, { start: s0, end: L });
          const lo = Object.assign({}, p, { start: 0, end: e0 });
          if (f.mixed && p.strand === -1) out.push(lo, hi); else out.push(hi, lo);
        }
      }
      return E.assembleFeature(f, out, L);
    }).filter(Boolean);
    return Object.assign({}, rec, { seq, features });
  };

  // Flip a whole record (reverse complement), moving features with it
  E.flipRecord = function (rec) {
    const L = rec.seq.length;
    return Object.assign({}, rec, {
      seq: E.revcomp(rec.seq),
      features: (rec.features || []).map(f => {
        if (f.between) return Object.assign({}, f, { start: L - f.start, end: L - f.start, strand: -(f.strand || 1) });
        if (f.mixed) {
          const parts = f.parts.map(p => ({ start: L - p.end, end: L - p.start, strand: -(p.strand || 1) }));
          return E.assembleFeature(f, parts, L);
        }
        const pieces = E.featurePieces(f, L).map(p => ({ start: L - p.end, end: L - p.start })).reverse();
        return E.assembleFeature(Object.assign({}, f, { strand: -(f.strand || 1) }), pieces, L);
      }).filter(Boolean)
    });
  };


  /* ======================= BENCH ARITHMETIC =======================
     What turns a design into a reaction: how many moles of each fragment,
     how much of each stock, and what to put in the block.

     Molecular weight is taken from the real base composition whenever the
     sequence is known. The usual shortcut of 650 g/mol per base pair is
     out by about 5%, and on an assembly that error goes straight into the
     molar ratio of insert to vector.

     The kit recipes below are the published conditions for each kit. They
     are a starting point for a master mix, not a substitute for the
     manual that came with the box.
  ------------------------------------------------------------------ */
  E.NT_MW = { A: 313.21, T: 304.20, C: 289.18, G: 329.21 };
  E.BP_MW_AVG = 617.96;        // mean of the A:T (617.41) and G:C (618.39) pairs
  E.DS_ENDS_MW = 36.04;        // the water on the two 5'-OH ends of a duplex

  /* grams per mole of a linear duplex. Pass a sequence for the exact
     figure, or a number of base pairs for the average-composition one. */
  E.dnaMW = function (seqOrBp) {
    if (typeof seqOrBp === 'number') {
      return seqOrBp > 0 ? seqOrBp * E.BP_MW_AVG + E.DS_ENDS_MW : 0;
    }
    const s = E.cleanSeq(seqOrBp);
    if (!s.length) return 0;
    const comp = { A: 'T', T: 'A', G: 'C', C: 'G' };
    let mw = 0;
    for (let i = 0; i < s.length; i++) {
      const c = s[i], rc = comp[c];
      mw += rc ? (E.NT_MW[c] + E.NT_MW[rc]) : E.BP_MW_AVG;   // an N contributes the average pair
    }
    return mw + E.DS_ENDS_MW;
  };
  E.ngToPmol = function (ng, seqOrBp) {
    const mw = E.dnaMW(seqOrBp);
    return mw > 0 ? (ng * 1000) / mw : 0;
  };
  E.pmolToNg = function (pmol, seqOrBp) { return pmol * E.dnaMW(seqOrBp) / 1000; };
  E.ngToFmol = function (ng, seqOrBp) { return E.ngToPmol(ng, seqOrBp) * 1000; };

  /* C1V1 = C2V2 */
  E.dilute = function (stock, want, finalVol) {
    if (!(stock > 0) || !(want >= 0) || !(finalVol > 0)) return null;
    if (want > stock) return { ok: false, error: 'The stock is weaker than the concentration asked for.' };
    const v = want * finalVol / stock;
    return { ok: true, stockVol: v, diluentVol: finalVol - v, fold: stock / (want || 1) };
  };
  /* a lyophilised oligo: nmol delivered -> volume for a round stock */
  E.resuspend = function (nmol, uM) {
    if (!(nmol > 0) || !(uM > 0)) return null;
    return { volUl: nmol * 1000 / uM };           // nmol / (pmol/uL) = uL
  };

  /* ---------- polymerases ----------
     Volumes are given as a fraction of the reaction so any reaction size
     scales exactly. Annealing and extension rules are each supplier's own.  */
  E.POLYMERASES = {
    phusion: {
      name: 'Phusion High-Fidelity (NEB M0530)',
      denature: 98, initialDenatureS: 30, denatureS: 10, finalExtendMin: 5,
      extendC: 72, extendSPerKb: [15, 30], cycles: [25, 35],
      /* NEB: Ta = Tm of the lower primer + 3 when both are longer than 20 nt */
      anneal: (lo, lenF, lenR) => (lenF > 20 && lenR > 20) ? lo + 3 : lo,
      annealNote: 'lower primer Tm + 3 °C (both primers over 20 nt); the lower Tm itself if either is 20 nt or shorter',
      parts: [
        { name: '5× Phusion HF buffer', per: 1 / 5, final: '1×' },
        { name: '10 mM dNTP mix', per: 1 / 50, final: '200 µM each' },
        { name: 'Phusion polymerase (2 U/µL)', per: 1 / 100, final: '1 U per 50 µL' }
      ],
      optional: [{ name: 'DMSO', per: 0.03, final: '3%', why: 'GC-rich or structured template' }]
    },
    q5: {
      name: 'Q5 High-Fidelity (NEB M0491)',
      denature: 98, initialDenatureS: 30, denatureS: 10, finalExtendMin: 2,
      extendC: 72, extendSPerKb: [20, 30], cycles: [25, 35],
      anneal: lo => lo,
      annealNote: 'the lower primer Tm, as NEB’s Tm calculator reports it for Q5',
      parts: [
        { name: '5× Q5 reaction buffer', per: 1 / 5, final: '1×' },
        { name: '10 mM dNTP mix', per: 1 / 50, final: '200 µM each' },
        { name: 'Q5 polymerase (2 U/µL)', per: 1 / 100, final: '1 U per 50 µL' }
      ],
      optional: [{ name: '5× Q5 High GC enhancer', per: 1 / 5, final: '1×', why: 'GC-rich template' }]
    },
    mm2x: {
      name: 'A 2× master mix (Q5/Phusion/KAPA)',
      denature: 98, initialDenatureS: 30, denatureS: 10, finalExtendMin: 2,
      extendC: 72, extendSPerKb: [20, 30], cycles: [25, 35],
      anneal: lo => lo,
      annealNote: 'the lower primer Tm — check the rule for the mix you have',
      parts: [{ name: '2× master mix', per: 1 / 2, final: '1×' }],
      optional: []
    },
    taq: {
      name: 'Taq (NEB M0273)',
      denature: 95, initialDenatureS: 30, denatureS: 30, finalExtendMin: 5,
      extendC: 68, extendSPerKb: [60, 60], cycles: [30, 35],
      anneal: lo => lo - 5,
      annealNote: 'the lower primer Tm − 5 °C',
      parts: [
        { name: '10× standard Taq buffer', per: 1 / 10, final: '1×' },
        { name: '10 mM dNTP mix', per: 1 / 50, final: '200 µM each' },
        { name: 'Taq polymerase (5 U/µL)', per: 1 / 200, final: '1.25 U per 50 µL' }
      ],
      optional: []
    }
  };

  /* Everything a PCR needs, for n reactions with an overage. Primer and
     template volumes come from the stocks actually on the bench. */
  E.pcrMix = function (o) {
    o = o || {};
    const pol = E.POLYMERASES[o.polymerase || 'phusion'] || E.POLYMERASES.phusion;
    const vol = o.volume || 50;
    const n = Math.max(1, o.reactions || 1);
    const over = o.overage == null ? 0.1 : o.overage;
    const primerStock = o.primerStockUM || 10;      // µM
    const primerFinal = o.primerFinalUM || (pol === E.POLYMERASES.taq ? 0.2 : 0.5);
    const perPrimer = vol * primerFinal / primerStock;
    const tplStock = o.templateStockNgUl || 0;
    const tplNg = o.templateNg == null ? 1 : o.templateNg;
    const tplVol = tplStock > 0 ? tplNg / tplStock : null;

    const rows = [];
    pol.parts.forEach(p => rows.push({ name: p.name, perRxn: vol * p.per, final: p.final, mix: true }));
    (o.extras || []).forEach(k => {
      const x = (pol.optional || []).find(y => y.name === k);
      if (x) rows.push({ name: x.name, perRxn: vol * x.per, final: x.final, mix: true });
    });
    rows.push({ name: 'Forward primer (' + f2(primerStock) + ' µM)', perRxn: perPrimer,
                final: f2(primerFinal) + ' µM', mix: !!o.primersInMix });
    rows.push({ name: 'Reverse primer (' + f2(primerStock) + ' µM)', perRxn: perPrimer,
                final: f2(primerFinal) + ' µM', mix: !!o.primersInMix });
    rows.push({ name: 'Template' + (tplStock > 0 ? ' (' + f2(tplStock) + ' ng/µL)' : ''),
                perRxn: tplVol, final: f2(tplNg) + ' ng', mix: false, unknown: tplVol == null });

    const used = rows.reduce((a, r) => a + (r.perRxn || 0), 0);
    const water = vol - used;
    rows.splice(pol.parts.length, 0, { name: 'Nuclease-free water', perRxn: water, final: 'to ' + vol + ' µL',
                                       mix: true, water: true });
    const mixRows = rows.filter(r => r.mix);
    const mixPer = mixRows.reduce((a, r) => a + (r.perRxn || 0), 0);
    const factor = n * (1 + over);
    return {
      ok: water >= 0,
      error: water < 0 ? 'The volumes already come to more than ' + vol + ' µL. Use a stronger template stock or a bigger reaction.' : '',
      polymerase: pol, volume: vol, reactions: n, overage: over, factor: factor,
      rows: rows.map(r => Object.assign({}, r, { total: r.perRxn == null ? null : r.perRxn * factor })),
      masterMixPerRxn: mixPer,
      masterMixTotal: mixPer * factor,
      addPerTube: vol - mixPer,
      water: water
    };
  };

  /* Annealing temperature and extension time for a designed pair. */
  E.pcrProgram = function (o) {
    o = o || {};
    const pol = E.POLYMERASES[o.polymerase || 'phusion'] || E.POLYMERASES.phusion;
    const lo = Math.min(o.tmF, o.tmR);
    const ta = pol.anneal(lo, o.lenF || 25, o.lenR || 25);
    const kb = (o.productBp || 1000) / 1000;
    const ext = pol.extendSPerKb.map(s => Math.max(15, Math.round(s * kb)));
    return {
      polymerase: pol, anneal: ta, annealNote: pol.annealNote,
      lower: lo, spread: Math.abs(o.tmF - o.tmR),
      extendS: ext, extendC: pol.extendC,
      denature: pol.denature, initialDenatureS: pol.initialDenatureS, denatureS: pol.denatureS,
      cycles: pol.cycles, finalExtendMin: pol.finalExtendMin
    };
  };

  /* ---------- assembly kits ---------- */
  E.ASSEMBLY_KITS = {
    nebuilder: {
      name: 'NEBuilder HiFi (NEB E2621)', mixName: '2× NEBuilder HiFi master mix',
      mixFraction: 0.5, totalVol: 20, tempC: 50, minutes: [15, 60],
      pmolRange: [0.03, 0.2], pmolNote: '0.03–0.2 pmol total for two or three fragments',
      vectorNg: 100
    },
    gibson: {
      name: 'Gibson Assembly (NEB E5510)', mixName: '2× Gibson Assembly master mix',
      mixFraction: 0.5, totalVol: 20, tempC: 50, minutes: [15, 60],
      pmolRange: [0.02, 0.5], pmolNote: '0.02–0.5 pmol total',
      vectorNg: 100
    },
    infusion: {
      name: 'In-Fusion Snap Assembly (Takara)', mixName: '5× In-Fusion Snap Assembly master mix',
      mixFraction: 0.2, totalVol: 10, tempC: 50, minutes: [15, 15],
      pmolRange: [0.02, 0.3], pmolNote: '50–200 ng of linearised vector, insert at about 2:1',
      vectorNg: 100
    },
    kld: {
      name: 'KLD blunt closure (NEB M0554)', mixName: '2× KLD reaction buffer',
      mixFraction: 0.5, totalVol: 10, tempC: 22, minutes: [5, 5],
      pmolRange: null, pmolNote: '1–5 ng/µL of PCR product, 1 µL into the reaction',
      vectorNg: null, enzyme: { name: 'KLD enzyme mix', vol: 1 }
    }
  };

  /* How much of each fragment, by moles, and how many µL of each stock.
     fragments: [{ name, seq | bp, stockNgUl, role: 'vector'|'insert' }]      */
  E.assemblySetup = function (o) {
    o = o || {};
    const kit = E.ASSEMBLY_KITS[o.kit || 'nebuilder'] || E.ASSEMBLY_KITS.nebuilder;
    const frags = (o.fragments || []).map(f => {
      const size = f.seq && f.seq.length ? f.seq : (f.bp || 0);
      return { name: f.name || 'fragment', role: f.role || 'insert',
               bp: f.seq && f.seq.length ? E.cleanSeq(f.seq).length : (f.bp || 0),
               mw: E.dnaMW(size), stockNgUl: f.stockNgUl || 0 };
    }).filter(f => f.bp > 0);
    if (!frags.length) return { ok: false, error: 'Nothing to assemble.' };
    const vec = frags.find(f => f.role === 'vector') || frags[0];
    const ratio = o.ratio == null ? 2 : o.ratio;
    const vecNg = o.vectorNg == null ? (kit.vectorNg || 100) : o.vectorNg;
    const vecPmol = (vecNg * 1000) / vec.mw;

    const total = o.totalVol || kit.totalVol;
    const rows = frags.map(f => {
      const pmol = f === vec ? vecPmol : vecPmol * ratio;
      const ng = pmol * f.mw / 1000;
      return { name: f.name, role: f.role, bp: f.bp, mw: f.mw, pmol: pmol, ng: ng,
               stockNgUl: f.stockNgUl,
               vol: f.stockNgUl > 0 ? ng / f.stockNgUl : null };
    });
    const dnaVol = rows.reduce((a, r) => a + (r.vol || 0), 0);
    const mixVol = total * kit.mixFraction;
    const enzVol = kit.enzyme ? kit.enzyme.vol : 0;
    const water = total - mixVol - enzVol - dnaVol;
    const totalPmol = rows.reduce((a, r) => a + r.pmol, 0);
    const warn = [];
    if (kit.pmolRange && totalPmol > kit.pmolRange[1])
      warn.push('that comes to ' + f2(totalPmol) + ' pmol of DNA, above the ' + kit.pmolNote +
                ' this kit asks for — scale the vector down');
    if (kit.pmolRange && totalPmol < kit.pmolRange[0])
      warn.push('that comes to ' + f2(totalPmol) + ' pmol of DNA, below the ' + kit.pmolNote +
                ' this kit asks for — use more vector');
    if (rows.some(r => r.vol == null)) warn.push('enter the concentration of each stock to get volumes');
    else if (water < 0) warn.push('the DNA alone is more than the reaction holds — concentrate the stocks or scale the reaction up');
    return {
      ok: true, kit: kit, ratio: ratio, rows: rows, totalPmol: totalPmol,
      mixVol: mixVol, enzVol: enzVol, waterVol: water, totalVol: total,
      vectorNg: vecNg, warnings: warn
    };
  };

  function f2(x) { return (Math.round(x * 100) / 100).toString(); }


  /* ======================= PRIMER QUALITY =======================
     A primer is not good because its melting temperature is right. It is
     good because it anneals where it should, once, and does not spend the
     reaction folded on itself or stuck to its partner. Length is usually
     the only free variable — the 5' end is pinned by the junction — but
     within a few bases either way the structural behaviour changes a lot,
     so it is worth searching rather than taking the first length that
     reaches the target.

     Penalties are expressed in degrees-Celsius-equivalent so they can be
     weighed against the melting temperature on one scale.
  ------------------------------------------------------------------ */
  E.QUALITY = {
    tmWeight: 1,
    hairpinPerBp: 1.6, hairpin3Mult: 2.2,
    dimerPerBp: 1.1, dimer3Mult: 2.6,
    crossPerBp: 1.0, cross3Mult: 2.4,
    noClamp: 2.0, tooSticky: 1.4,
    runPenalty: 2.0, gcOut: 0.35,
    extraSite: 12, lengthOut: 0.4,
    minStemCount: 4, minDimerCount: 5
  };

  E.primerQuality = function (seq, opts) {
    opts = opts || {};
    const Q = E.QUALITY;
    const s = E.cleanSeq(seq);
    const cond = opts.conditions || E.PCR_CONDITIONS;
    const target = opts.targetTm || 60;
    const tm = E.tmNN(s, cond);
    const gc = E.gcContent(s);
    const flags = [];
    let p = 0;

    const dTm = isFinite(tm) ? Math.abs(tm - target) : 99;
    p += dTm * Q.tmWeight;

    const hp = E.hairpin(s);
    if (hp.stemLength >= Q.minStemCount) {
      const near3 = hp.end >= s.length - 5;
      const add = (hp.stemLength - (Q.minStemCount - 1)) * Q.hairpinPerBp * (near3 ? Q.hairpin3Mult : 1);
      p += add;
      flags.push({ level: near3 ? 'bad' : 'warn',
                   text: hp.stemLength + ' bp hairpin stem' + (near3 ? ' at the 3′ end' : '') });
    }
    const sd = E.selfDimer(s);
    if (sd.matches >= Q.minDimerCount) {
      const add = (sd.matches - (Q.minDimerCount - 1)) * Q.dimerPerBp * (sd.threePrime ? Q.dimer3Mult : 1);
      p += add;
      flags.push({ level: sd.threePrime ? 'bad' : 'warn',
                   text: sd.matches + ' bp self-dimer' + (sd.threePrime ? ' involving the 3′ end' : '') });
    }
    const tail5 = s.slice(-5);
    const clamp = (tail5.match(/[GC]/g) || []).length;
    if (clamp === 0) { p += Q.noClamp; flags.push({ level: 'warn', text: 'no G or C in the last 5 bases' }); }
    else if (clamp >= 4) { p += Q.tooSticky; flags.push({ level: 'warn', text: 'four or more G/C in the last 5 bases' }); }
    const runs = s.match(/([ACGT])\1{3,}/g) || [];
    runs.forEach(r => {
      p += Q.runPenalty * (r.length >= 5 ? 1.5 : 1);
      flags.push({ level: r.length >= 5 ? 'bad' : 'warn', text: r.length + ' identical bases in a row (' + r[0] + ')' });
    });
    if (gc < 40 || gc > 60) {
      p += Math.abs(gc < 40 ? 40 - gc : gc - 60) * Q.gcOut;
      flags.push({ level: 'warn', text: 'GC ' + Math.round(gc) + '%' });
    }
    if (s.length < 18 || s.length > 32) p += Math.abs(s.length < 18 ? 18 - s.length : s.length - 32) * Q.lengthOut;

    let sites = 1;
    if (opts.template) {
      try { sites = E.primerSites(opts.template, s, { circular: !!opts.circular }).length; } catch (e) { sites = 1; }
      if (sites !== 1) {
        p += Q.extraSite * Math.max(1, Math.abs(sites - 1));
        flags.push({ level: 'bad', text: sites === 0 ? 'does not anneal to the template' : 'anneals in ' + sites + ' places' });
      }
    }
    return { seq: s, length: s.length, tm: tm, gc: gc, sites: sites,
             hairpin: hp, selfDimer: sd, penalty: p, flags: flags };
  };

  /* Candidates for one primer: the 5' end is fixed by the junction, so
     length is varied, and the whole set is scored. `shift` optionally
     allows the 5' end to move, which moves the junction with it — the
     caller decides whether that is allowed and the result says by how
     much. */
  E.tunePrimer = function (template, o) {
    o = o || {};
    const L = E.cleanSeq(template).length;
    const dbl = E.cleanSeq(template) + E.cleanSeq(template);
    const minLen = o.minLen || 18, maxLen = o.maxLen || 36;
    const shift = o.shift || 0;
    const out = [];
    for (let d = -shift; d <= shift; d++) {
      for (let n = minLen; n <= maxLen; n++) {
        let bind;
        if (o.direction === 'forward') bind = dbl.substr(((o.at + d) % L + L) % L, n);
        else bind = dbl.substr(((o.at + d - n) % L + L) % L, n);
        if (bind.length < n) continue;
        const seq = o.direction === 'forward' ? bind : E.revcomp(bind);
        const q = E.primerQuality(seq, { targetTm: o.targetTm, conditions: o.conditions,
                                         template: o.template || template, circular: o.circular });
        out.push(Object.assign(q, { shift: d, bindStart: ((o.at + d - (o.direction === 'forward' ? 0 : n)) % L + L) % L }));
      }
    }
    out.sort((a, b) => a.penalty - b.penalty);
    return out;
  };

  /* The pair together: the two have to match each other as well as the
     template, so the cross-dimer and the gap between their melting
     temperatures are scored on top of the two individual scores. */
  E.tunePair = function (template, o) {
    o = o || {};
    const Q = E.QUALITY;
    const fs = E.tunePrimer(template, Object.assign({}, o, { at: o.fAt, direction: 'forward' })).slice(0, o.depth || 14);
    const rs = E.tunePrimer(template, Object.assign({}, o, { at: o.rAt, direction: 'reverse' })).slice(0, o.depth || 14);
    if (!fs.length || !rs.length) return { ok: false, error: 'No primer fits there.' };
    /* A structural problem is worth a few degrees, not ten. Candidates
       whose melting temperature strays outside the band are kept and
       labelled so the trade can be seen, but they are not chosen
       automatically: an 8 °C overshoot costs more than the dimer it
       avoids. */
    const band = o.tmBand == null ? 5 : o.tmBand;
    const target = o.targetTm || 60;
    const pairs = [];
    fs.forEach(f => rs.forEach(r => {
      if (f.shift !== 0 && r.shift !== 0 && o.shiftOne) return;
      const x = E.dimer(f.seq, r.seq);
      let p = f.penalty + r.penalty;
      const spread = Math.abs(f.tm - r.tm);
      p += spread * 1.2;
      const flags = [];
      if (x.matches >= Q.minDimerCount) {
        p += (x.matches - (Q.minDimerCount - 1)) * Q.crossPerBp * (x.threePrime ? Q.cross3Mult : 1);
        flags.push({ level: x.threePrime ? 'bad' : 'warn',
                     text: x.matches + ' bp between the two primers' + (x.threePrime ? ', 3′ end involved' : '') });
      }
      if (spread > 5) flags.push({ level: 'warn', text: f1(spread) + ' °C apart' });
      const out = Math.max(Math.abs(f.tm - target), Math.abs(r.tm - target)) > band;
      pairs.push({ forward: f, reverse: r, cross: x, spread: spread, penalty: p,
                   flags: flags, shift: { f: f.shift, r: r.shift }, outOfBand: out });
    }));
    pairs.sort((a, b) => (a.outOfBand ? 1 : 0) - (b.outOfBand ? 1 : 0) || a.penalty - b.penalty);
    /* keep one entry per distinct pair of lengths-and-shifts, best first */
    const seen = {}, uniq = [];
    pairs.forEach(p => {
      const k = p.forward.seq + '|' + p.reverse.seq;
      if (seen[k]) return; seen[k] = 1; uniq.push(p);
    });
    return { ok: true, best: uniq[0], options: uniq.slice(0, o.keep || 24) };
  };
  /* An assembly primer is a tail that must match the neighbour exactly and
     an annealing part that must sit on its own template. Both lengths are
     free — a longer arm is still a correct arm — so both are searched. The
     finished construct is the same whichever is chosen. */
  E.tuneAssemblyPrimer = function (o) {
    o = o || {};
    const insert = E.cleanSeq(o.insert || '');
    const neighbour = E.cleanSeq(o.neighbour || '');     // the fragment this end joins to
    const atStart = o.end !== '3';
    const armLo = o.armMin || 15, armHi = Math.min(o.armMax || 40, neighbour.length);
    const anLo = o.annealMin || 16, anHi = Math.min(o.annealMax || 34, insert.length);
    const target = o.targetTm || 60;
    const band = o.tmBand == null ? 5 : o.tmBand;
    const out = [];
    for (let arm = armLo; arm <= armHi; arm++) {
      const tail = atStart ? neighbour.slice(-arm) : E.revcomp(neighbour.slice(0, arm));
      for (let an = anLo; an <= anHi; an++) {
        const anneal = atStart ? insert.slice(0, an) : E.revcomp(insert.slice(insert.length - an));
        const seq = tail + anneal;
        /* the annealing part alone is what binds in the first cycles, so
           that is what the melting temperature is judged on */
        const qa = E.primerQuality(anneal, { targetTm: target, conditions: o.conditions,
                                             template: o.template, circular: false });
        const qf = E.primerQuality(seq, { targetTm: target, conditions: o.conditions });
        /* the whole primer's structure matters, the whole primer's Tm does not */
        const penalty = qa.penalty + (qf.penalty - Math.abs(qf.tm - target)) * 0.6 + (arm - armLo) * 0.05;
        out.push({ seq: seq, tail: tail, anneal: anneal, arm: arm, annealLen: an,
                   length: seq.length, tmAnneal: qa.tm, tmFull: qf.tm,
                   flags: qa.flags.concat(qf.flags.filter(x => !/GC |no G or C|identical bases/.test(x.text))),
                   penalty: penalty, outOfBand: Math.abs(qa.tm - target) > band });
      }
    }
    out.sort((a, b) => (a.outOfBand ? 1 : 0) - (b.outOfBand ? 1 : 0) || a.penalty - b.penalty);
    return { ok: out.length > 0, best: out[0], options: out.slice(0, o.keep || 24) };
  };

  function f1(x) { return (Math.round(x * 10) / 10).toString(); }


  /* ======================= WHAT IS ON THE ENDS OF A FRAGMENT =======================
     A synthesised fragment rarely arrives as bare coding sequence. Whoever
     ordered it put something on the ends: a protease site, a tag, a
     restriction site, an overhang for whatever vector they had in mind.
     When the vector already supplies the same thing, including it again
     puts two of them in the protein — and nothing about the DNA looks
     wrong, so it is only caught by reading the translation.
  ------------------------------------------------------------------ */
  E.PEPTIDE_ADDONS = [
    { name: 'TEV protease site', kind: 'protease site', aa: 'ENLYFQ', cleaves: 'after the Q' },
    { name: 'TEV protease site (partial)', kind: 'protease site', aa: 'LYFQG', partialOf: 'TEV protease site' },
    { name: 'TEV protease site (partial)', kind: 'protease site', aa: 'LYFQS', partialOf: 'TEV protease site' },
    { name: 'HRV 3C / PreScission site', kind: 'protease site', aa: 'LEVLFQGP', cleaves: 'between Q and G' },
    { name: 'Thrombin site', kind: 'protease site', aa: 'LVPRGS', cleaves: 'between R and G' },
    { name: 'Factor Xa site', kind: 'protease site', aa: 'IEGR', cleaves: 'after the R' },
    { name: 'Enterokinase site', kind: 'protease site', aa: 'DDDDK', cleaves: 'after the K' },
    { name: 'His tag', kind: 'tag', aa: 'HHHHHH' },
    { name: 'FLAG tag', kind: 'tag', aa: 'DYKDDDDK' },
    { name: 'Strep-II tag', kind: 'tag', aa: 'WSHPQFEK' },
    { name: 'Myc tag', kind: 'tag', aa: 'EQKLISEEDL' },
    { name: 'HA tag', kind: 'tag', aa: 'YPYDVPDYA' },
    { name: 'S-tag', kind: 'tag', aa: 'KETAAAKFERQHMDS' },
    { name: 'flexible linker', kind: 'linker', aa: 'GGGGS' },
    { name: 'flexible linker', kind: 'linker', aa: 'GSGSGS' }
  ];

  /* What sits on the ends of a fragment, and whether the vector already
     provides it. `o.vectorUpstream` is the backbone's last bases, read in
     the frame the fusion will have; `o.vectorDownstream` its first bases. */
  E.insertAddons = function (ins, o) {
    o = o || {};
    const s = E.cleanSeq(ins);
    const table = o.table || 11;
    if (s.length < 30) return { ok: false, leading: [], trailing: [] };
    const prot = E.translate(s, { table: table });
    const stopAt = prot.replace(/\*+$/, '').indexOf('*');
    const coding = stopAt >= 0 ? prot.slice(0, stopAt) : prot.replace(/\*+$/, '');
    const upProt = o.vectorUpstream ? E.translate(o.vectorUpstream, { table: table }) : '';

    const window = coding.slice(0, 25);
    const leading = [];
    E.PEPTIDE_ADDONS.forEach(m => {
      const i = window.indexOf(m.aa);
      if (i < 0 || i > 12) return;
      /* a longer motif covering the same bases wins */
      const from = i * 3, to = (i + m.aa.length) * 3;
      if (leading.some(x => x.from <= from && x.to >= to)) return;
      leading.push({ kind: m.kind, name: m.name, aa: m.aa, from: from, to: to,
                     aaFrom: i + 1, aaTo: i + m.aa.length,
                     inVector: !!(upProt && upProt.indexOf(m.aa) >= 0) ||
                               !!(upProt && m.partialOf === 'TEV protease site' && /ENLYFQ|LYFQ/.test(upProt)),
                     cleaves: m.cleaves || '' });
    });
    /* drop a partial that sits inside a full one already found */
    const kept = leading.filter(x => !leading.some(y => y !== x && y.from <= x.from && y.to >= x.to && (y.to - y.from) > (x.to - x.from)));
    kept.sort((a, b) => a.from - b.from);

    /* restriction sites riding on either end */
    const sites = [];
    if (o.enzymes && o.enzymes.length) {
      const head = s.slice(0, Math.min(45, s.length)), tail = s.slice(Math.max(0, s.length - 45));
      o.enzymes.forEach(e => {
        try {
          E.findSites(head, e, false).forEach(h => sites.push({ where: 'start', name: e.name, at: h.pos }));
          E.findSites(tail, e, false).forEach(h => sites.push({ where: 'end', name: e.name,
            at: (s.length - tail.length) + h.pos }));
        } catch (err) { /* an enzyme that will not search is simply skipped */ }
      });
    }

    /* where the coding sequence really begins: after the last add-on, at
       the next ATG in frame if there is one close by */
    let startBase = kept.length ? kept[kept.length - 1].to : 0;
    let atg = -1;
    for (let i = startBase; i <= Math.min(s.length - 3, startBase + 18); i += 3) {
      if (s.substr(i, 3) === 'ATG') { atg = i; break; }
    }
    if (atg >= 0) startBase = atg;

    /* trailing: everything after the stop codon */
    let trailing = null;
    if (stopAt >= 0) {
      const afterStop = (stopAt + 1) * 3;
      if (s.length - afterStop > 0)
        trailing = { from: afterStop, to: s.length, bp: s.length - afterStop,
                     sites: sites.filter(x => x.where === 'end' && x.at >= afterStop - 3) };
    }
    return {
      ok: true, protein: coding, leading: kept, sites: sites,
      suggestFrom: startBase, suggestTo: stopAt >= 0 ? (stopAt + 1) * 3 : s.length,
      stopAtResidue: stopAt >= 0 ? stopAt + 1 : null, trailing: trailing
    };
  };


  /* ======================= READING A FRAGMENT =======================
     A tube of DNA from a previous student, a gBlock ordered two years ago,
     a sequence pasted out of a thesis: what is actually in it? The coding
     part, and everything somebody added around it. Each base is assigned
     to one segment, so the parts add up to the whole and nothing is
     quietly left out of the account.
  ------------------------------------------------------------------ */
  E.ANATOMY_EXTRA = [
    { name: 'Avi tag (biotinylation)', kind: 'tag', aa: 'GLNDIFEAQKIEWHE' },
    { name: 'His tag (8x)', kind: 'tag', aa: 'HHHHHHHH' },
    { name: 'Twin-Strep tag', kind: 'tag', aa: 'SAWSHPQFEK' },
    { name: 'V5 tag', kind: 'tag', aa: 'GKPIPNPLLGLDST' },
    { name: 'SUMO C-terminus', kind: 'fusion partner', aa: 'LRLRGG' },
    { name: 'flexible linker', kind: 'linker', aa: 'GGSGG' },
    { name: 'rigid linker', kind: 'linker', aa: 'EAAAK' },
    { name: 'TEV protease site', kind: 'protease site', aa: 'ENLYFQG' },
    { name: 'TEV protease site', kind: 'protease site', aa: 'ENLYFQS' }
  ];

  /* Pick the frame to read it in: the one whose longest run without a stop
     is longest, preferring a run that begins at an ATG. The reverse strand
     is tried too and reported, because a fragment is as often pasted the
     wrong way round as not. */
  E.bestFrame = function (seq, opts) {
    opts = opts || {};
    const table = opts.table || 11;
    const s = E.cleanSeq(seq);
    const tries = [];
    [0, 1, 2].forEach(f => {
      [1, -1].forEach(str => {
        if (str === -1 && opts.forward) return;
        const src = str === 1 ? s : E.revcomp(s);
        const body = src.slice(f, f + Math.floor((src.length - f) / 3) * 3);
        const prot = E.translate(body, { table: table });
        let best = { len: 0, at: 0 };
        let run = 0, start = 0;
        for (let i = 0; i <= prot.length; i++) {
          const stop = i === prot.length || prot[i] === '*';
          if (stop) {
            if (run > best.len) best = { len: run, at: start };
            run = 0; start = i + 1;
          } else { if (run === 0) start = i; run++; }
        }
        /* the real discriminator is the longest ATG-to-stop reading frame,
           not the longest stretch without a stop: a frame can run a long
           way without stopping and still code for nothing */
        let orf = { len: 0, at: best.at };
        for (let i = 0; i < prot.length; i++) {
          if (prot[i] !== 'M') continue;
          const stop = prot.indexOf('*', i);
          const len = (stop < 0 ? prot.length : stop) - i;
          if (len > orf.len) orf = { len: len, at: i, closed: stop >= 0 };
          if (stop >= 0) i = i;        // keep scanning; a later ATG may open a longer one
        }
        const hasMet = orf.len > 0;
        tries.push({ frame: f, strand: str, protein: prot,
                     runLength: hasMet ? orf.len : best.len, runAt: hasMet ? orf.at : best.at,
                     hasMet: hasMet, closed: !!orf.closed,
                     score: (hasMet ? orf.len * 2 + (orf.closed ? 20 : 0) : best.len * 0.5) +
                            (str === 1 ? 0.5 : 0) });
      });
    });
    tries.sort((a, b) => b.score - a.score);
    return tries[0];
  };


  /* A solubility or purification partner is a whole protein domain, not a
     six-residue motif: MBP is 370 residues, NusA nearly 500. Carrying full
     reference sequences for each would be a database; what is used instead
     is a signature from each one, long enough to be unmistakable, matched
     with a few mismatches allowed so that a variant or a different species'
     version is still recognised. The extent is then taken from the typical
     length of that partner, trimmed back to a protease site when one sits
     where the join should be.

     This is a starter set of the partners met most often, not a complete
     list. Anything not in it is reported as an unrecognised stretch rather
     than guessed at. */
  E.FUSION_PARTNERS = [
    { name: 'MBP (maltose-binding protein)', len: 370, why: 'solubility and affinity',
      sig: 'KIEEGKLVIWINGDKGYNGLAEVGKKFEKDTGIKV' },
    { name: 'GST (glutathione S-transferase)', len: 220, why: 'solubility and affinity',
      sig: 'MSPILGYWKIKGLVQPTRLLLEYLEEKYEE' },
    { name: 'NusA', len: 495, why: 'solubility',
      sig: 'MNKEILAVVEAVSNEKALPREKIFEALESALATATKKKY' },
    { name: 'thioredoxin (TrxA)', len: 109, why: 'solubility',
      sig: 'DFWAEWCGPCKMIAPILDEIA' },
    { name: 'SUMO (Smt3)', len: 98, why: 'solubility, cleaved by Ulp1',
      sig: 'LYDGIRIQADQTPEDLDMEDNDIIEAHREQIGG' },
    { name: 'GB1 (protein G B1 domain)', len: 56, why: 'solubility and expression',
      sig: 'MTYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE' }
  ];

  /* the best place a signature sits in a protein, allowing mismatches */
  E.findPeptideFuzzy = function (prot, sig, maxMismatchFrac) {
    const lim = Math.floor(sig.length * (maxMismatchFrac == null ? 0.2 : maxMismatchFrac));
    let best = null;
    for (let i = 0; i + sig.length <= prot.length; i++) {
      let miss = 0;
      for (let j = 0; j < sig.length; j++) {
        if (prot[i + j] !== sig[j]) { miss++; if (miss > lim) break; }
      }
      if (miss <= lim && (!best || miss < best.miss)) best = { at: i, miss: miss };
      if (best && best.miss === 0) break;
    }
    return best;
  };

  E.fragmentAnatomy = function (seq, opts) {
    opts = opts || {};
    const table = opts.table || 11;
    const s = E.cleanSeq(seq);
    if (s.length < 24) return { ok: false, error: 'Too short to read.' };
    const fr = opts.frame != null
      ? { frame: opts.frame, strand: opts.strand || 1,
          protein: E.translate(E.cleanSeq(opts.strand === -1 ? E.revcomp(s) : s).slice(opts.frame), { table: table }),
          runLength: 0, runAt: 0 }
      : E.bestFrame(s, { table: table, forward: opts.forward });
    const work = fr.strand === 1 ? s : E.revcomp(s);
    const off = fr.frame;
    const prot = E.translate(work.slice(off, off + Math.floor((work.length - off) / 3) * 3), { table: table });

    /* the open reading frame proper: the longest stop-free run, trimmed to
       begin at its first methionine when it has one */
    let runAt = fr.runAt, runLen = fr.runLength;
    if (!runLen) {
      let best = { len: 0, at: 0 }, run = 0, start = 0;
      for (let i = 0; i <= prot.length; i++) {
        const stop = i === prot.length || prot[i] === '*';
        if (stop) { if (run > best.len) best = { len: run, at: start }; run = 0; start = i + 1; }
        else { if (run === 0) start = i; run++; }
      }
      runAt = best.at; runLen = best.len;
    }
    const motifs = E.PEPTIDE_ADDONS.concat(E.ANATOMY_EXTRA);
    const found = [];
    motifs.forEach(m => {
      let i = prot.indexOf(m.aa);
      while (i >= 0) {
        found.push({ kind: m.kind, name: m.name, aa: m.aa, aaFrom: i, aaTo: i + m.aa.length });
        i = prot.indexOf(m.aa, i + 1);
      }
    });
    /* a shorter motif inside a longer one is the same thing seen twice */
    const motifsKept = found.filter(x => !found.some(y => y !== x && y.aaFrom <= x.aaFrom && y.aaTo >= x.aaTo &&
      (y.aaTo - y.aaFrom) > (x.aaTo - x.aaFrom)));
    motifsKept.sort((a, b) => a.aaFrom - b.aaFrom);

    /* where the protein that matters starts and ends: inside the run,
       after any add-on that sits at its front */
    let coreFrom = runAt, coreTo = runAt + runLen;
    motifsKept.forEach(m => {
      if (m.aaFrom <= coreFrom + 2 && m.aaTo > coreFrom && m.aaTo < coreTo) coreFrom = m.aaTo;
      if (m.aaTo >= coreTo - 2 && m.aaFrom < coreTo && m.aaFrom > coreFrom) coreTo = m.aaFrom;
    });
    const met = prot.indexOf('M', coreFrom);
    if (met >= 0 && met < coreFrom + 6) coreFrom = met;
    /* a coding sequence conventionally takes in its own stop codon, and
       leaving it out would report three stray bases as something added */
    const stopIncluded = prot[coreTo] === '*';
    const codeTo = stopIncluded ? coreTo + 1 : coreTo;

    /* ---- fusion partners carried on the front ----
       A solubility tag is a domain, so it is found by its signature and its
       extent taken from its known length, trimmed back to a protease site
       if one sits near where the join should be. */
    const aaToBase = i => off + i * 3;
    const partners = [];
    if (!opts.noPartners) {
      E.FUSION_PARTNERS.forEach(fp => {
        const hit = E.findPeptideFuzzy(prot.slice(coreFrom, coreTo), fp.sig, 0.2);
        if (!hit) return;
        let at = coreFrom + hit.at;
        /* the initiator methionine in front of the signature belongs to the
           partner, not to a one-residue protein of its own */
        if (at - coreFrom <= 2) at = coreFrom;
        let end = Math.min(coreTo, at + fp.len);
        /* the first protease site after it is where the join really is, and
           it is a better boundary than any nominal length */
        let cut = null;
        motifsKept.forEach(m => {
          if (m.kind !== 'protease site' || m.aaFrom <= at) return;
          if (m.aaFrom < end && (cut == null || m.aaFrom < cut)) cut = m.aaFrom;
        });
        if (cut != null) end = cut;
        partners.push({ name: fp.name, why: fp.why, aaFrom: at, aaTo: end,
                        mismatches: hit.miss, residues: end - at });
      });
      partners.sort((a, b) => a.aaFrom - b.aaFrom);
    }

    /* every base assigned to exactly one segment */
    const segs = [];
    const push = (from, to, kind, name, extra) => {
      if (to <= from) return;
      /* Every segment carries its own bases. Without them the caller has to
         slice the fragment itself and get the strand and the offset right to
         do it, which is the one place a reader of this table cannot check
         the answer. `dna` is on the working strand, like from and to. */
      segs.push(Object.assign({ from: from, to: to, bp: to - from, kind: kind, name: name,
        dna: work.slice(from, to) }, extra || {}));
    };
    const coreA = aaToBase(coreFrom), coreB = aaToBase(codeTo);
    /* Which motifs get a segment of their own. A motif outside the reading
       frame always does. Inside it, only one that genuinely divides the
       protein: a protease site, a recognised partner, or a tag of six
       residues or more. A five-residue linker occurring by chance in the
       middle of a real protein must not split it in two. */
    const inner = [], outerM = [];
    motifsKept.forEach(m => {
      const a = aaToBase(m.aaFrom), b = aaToBase(m.aaTo);
      const mk = { from: a, to: b, kind: m.kind, name: m.name, aa: m.aa };
      if (a >= coreA && b <= coreB) {
        if (m.kind === 'protease site' || m.aa.length >= 6) inner.push(mk);
      } else outerM.push(mk);
    });
    /* the partner's whole peptide, not a preview of it: a caller that wants
       twelve residues can take twelve, but one that wants to copy the domain
       out cannot put back what was elided */
    partners.forEach(p => inner.push({ from: aaToBase(p.aaFrom), to: aaToBase(p.aaTo),
      kind: 'fusion partner', name: p.name, aa: prot.slice(p.aaFrom, p.aaTo),
      partner: p }));
    inner.sort((a, b) => a.from - b.from || b.to - a.to);
    /* drop an inner mark swallowed by a bigger one, so the parts still tile */
    const innerKept = inner.filter((m, i) => !inner.some((n, j) => j !== i && n.from <= m.from && n.to >= m.to &&
      (n.to - n.from) > (m.to - m.from)));

    const marks = outerM.concat(innerKept).sort((a, b) => a.from - b.from);
    let cur = 0;
    let codingParts = 0;
    marks.forEach(m => {
      if (m.from > cur) {
        const inCore = cur >= coreA && m.from <= coreB;
        if (inCore) { codingParts++; push(cur, m.from, 'coding', 'protein', { aa: prot.slice(cur === coreA ? coreFrom : (cur - off) / 3, (m.from - off) / 3) }); }
        else {
          /* every coordinate in this result is on the working strand, so the
             bases have to be read from there too — off the pasted strand a
             fragment given in reverse would be named from the wrong letters */
          const piece = work.slice(cur, m.from);
          push(cur, m.from, 'other', cur === 0 ? (piece === 'ATG' ? 'start codon' : 'untranslated start') : 'spacer');
        }
      }
      push(Math.max(cur, m.from), m.to, m.kind, m.name, { aa: m.aa, partner: m.partner });
      cur = Math.max(cur, m.to);
    });
    if (cur < coreB) {
      codingParts++;
      push(Math.max(cur, coreA), coreB, 'coding', 'protein',
           { aa: prot.slice(Math.max(coreFrom, (Math.max(cur, coreA) - off) / 3), coreTo) });
      cur = coreB;
    }
    if (cur < s.length) {
      const stopHere = prot[coreTo] === '*';
      push(cur, s.length, 'other', stopHere ? 'after the stop codon' : 'untranslated end');
    }

    /* Restriction sites, on the working strand like everything else here.
       Searching the pasted strand instead put every site of a fragment given
       in reverse at length-minus-its-position, which then read as inside the
       protein when it was in the flank, and sent "how this was cloned before"
       to the wrong end of the fragment. A site is double-stranded, so the
       same cuts are found either way round; only where they are reported
       changes. */
    const sites = [];
    (opts.enzymes || []).forEach(e => {
      try {
        E.findSites(work, e, false).forEach(h => {
          const inCore = h.pos >= aaToBase(coreFrom) && h.pos < aaToBase(coreTo);
          sites.push({ name: e.name, at: h.pos, inCore: inCore });
        });
      } catch (err) { /* skip */ }
    });
    sites.sort((a, b) => a.at - b.at);

    return {
      ok: true, length: s.length, strand: fr.strand, frame: off,
      /* the strand everything here is measured on: the fragment as pasted
         when it reads forward, its reverse complement when it does not */
      work: work,
      partners: partners,
      proteinParts: segs.filter(g => g.kind === 'coding').length,
      /* what to clone when someone says "just the protein": the largest
         coding stretch, which on a fusion is the passenger and not the tag */
      main: (function () {
        const cs = segs.filter(g => g.kind === 'coding');
        if (!cs.length) return null;
        let b = cs[0];
        cs.forEach(g => { if (g.bp > b.bp) b = g; });
        return { from: b.from, to: b.to, bp: b.bp, aa: b.aa || '',
                 residues: Math.floor(b.bp / 3) - ((prot[coreTo] === '*' && b.to === coreB) ? 1 : 0) };
      })(),
      protein: prot, core: { aaFrom: coreFrom, aaTo: coreTo, from: aaToBase(coreFrom), to: aaToBase(codeTo),
                             aa: prot.slice(coreFrom, coreTo), residues: coreTo - coreFrom,
                             includesStop: stopIncluded },
      startsWithMet: prot[coreFrom] === 'M',
      endsWithStop: prot[coreTo] === '*',
      motifs: motifsKept.map(m => Object.assign({}, m, { from: aaToBase(m.aaFrom), to: aaToBase(m.aaTo),
        inside: aaToBase(m.aaFrom) >= aaToBase(coreFrom) && aaToBase(m.aaTo) <= aaToBase(codeTo) })),
      segments: segs, sites: sites
    };
  };

  /* ============================================================
     WHAT OF THIS READ IS ALREADY IN THE PLASMID

     A read off a sequencing primer starts in the vector and runs into the
     insert, so most of its front is backbone: promoter, operator, ribosome
     binding site, the vector's own ATG and tag. Deciding what to clone from
     it means knowing where the vector stops and the new sequence starts —
     and that is a question about two sequences, not about the read alone.

     The test is k-mer containment: a base of the read counts as already in
     the plasmid when some window of `k` bases covering it occurs in the
     plasmid, on either strand, with the plasmid read round the origin when
     it is circular. At k = 18 a window hits by chance about once in 7e10
     bases, so a match is the sequence and not a coincidence; and because
     the window must match exactly, a single real difference opens a gap of
     up to k bases rather than being passed over. Those gaps are reported as
     their own kind rather than as new sequence.

     Ambiguous bases cannot match anything, so a run of N is called
     unreadable and never new. Reporting an N-rich Sanger tail as "sequence
     to add" would be the worst answer this function could give: it reads
     like a finding and it would be cloned.
  ============================================================ */
  E.readInVector = function (read, rec, opts) {
    opts = opts || {};
    const k = opts.k || 18;
    const r = E.cleanSeq(read).replace(/U/g, 'T');
    const vec = rec && rec.seq ? E.cleanSeq(rec.seq).replace(/U/g, 'T') : '';
    if (r.length < k || vec.length < k) return null;
    const circular = rec.topology === 'circular';
    const L = vec.length;
    /* the plasmid read round the origin, so a window spanning base 1 is
       found rather than lost at the join */
    const hay = circular ? vec + vec.slice(0, k - 1) : vec;
    const rcHay = E.revcomp(hay);

    const index = new Map();
    const add = (text, strand) => {
      for (let i = 0; i + k <= text.length; i++) {
        const key = text.substr(i, k);
        if (key.indexOf('N') >= 0) continue;
        if (!index.has(key)) {
          /* a position on the plasmid itself, whichever strand it was found
             on, so the caller can say where in the plasmid this came from */
          const at = strand === 1 ? i % L : (text.length - i - k) % L;
          index.set(key, { at: at, strand: strand });
        }
      }
    };
    add(hay, 1);
    add(rcHay, -1);

    const hit = new Array(r.length).fill(null);
    for (let i = 0; i + k <= r.length; i++) {
      const w = r.substr(i, k);
      if (w.indexOf('N') >= 0) continue;
      const h = index.get(w);
      if (!h) continue;
      for (let j = 0; j < k; j++) {
        if (hit[i + j]) continue;
        /* the plasmid base this read base sits on */
        const vp = h.strand === 1 ? (h.at + j) % L : ((h.at + k - 1 - j) % L + L) % L;
        hit[i + j] = { at: vp, strand: h.strand };
      }
    }

    /* Classify every base, then run-length encode — but a matched run is
       only broken where the plasmid positions stop running on. Containment
       alone cannot see a join: read the last 150 bases of a plasmid
       followed by its first 150 and every base is found in the plasmid, yet
       that order exists nowhere in a linear one. Requiring consecutive read
       bases to sit on consecutive plasmid bases, on the same strand, turns
       the answer into "this stretch of your read is this stretch of the
       plasmid" and shows the join as the break it is. */
    const kindAt = i => hit[i] ? 'vector' : (r[i] === 'N' ? 'unreadable' : 'new');
    const runsOn = (a, b) => a && b && a.strand === b.strand &&
      b.at === (circular ? ((a.at + a.strand) % L + L) % L : a.at + a.strand);
    const blocks = [];
    let i = 0;
    while (i < r.length) {
      const kd = kindAt(i);
      let j = i + 1;
      while (j < r.length && kindAt(j) === kd &&
             (kd !== 'vector' || runsOn(hit[j - 1], hit[j]))) j++;
      const a = hit[i], z = hit[j - 1];
      blocks.push({ from: i, to: j, bp: j - i, kind: kd,
        /* half-open on the plasmid, in the direction the match runs */
        vecFrom: a ? (a.strand === 1 ? a.at : z.at) : null,
        vecTo: a ? ((a.strand === 1 ? z.at : a.at) + 1) : null,
        strand: a ? a.strand : 0 });
      i = j;
    }
    /* Nothing shorter than k can be resolved by a containment test at all,
       so a short unmatched run is never reported as sequence to clone.
       Between two matched runs it is a difference — an exact window cannot
       see past a changed base, so one substitution drops out a gap of up to
       k, and calling that new sequence would invite someone to clone a
       mismatch. At the very start or end of the read, next to a matched
       run, it is the ragged edge every Sanger trace has. Neighbours are
       taken past any ambiguous run, since an N tells you nothing either. */
    const near = (n, step) => {
      for (let m = n + step; m >= 0 && m < blocks.length; m += step) {
        if (blocks[m].kind !== 'unreadable') return blocks[m];
      }
      return null;
    };
    const firstReal = blocks.findIndex(b => b.kind !== 'unreadable');
    let lastReal = -1;
    blocks.forEach((b, n) => { if (b.kind !== 'unreadable') lastReal = n; });
    blocks.forEach((b, n) => {
      if (b.kind !== 'new' || b.bp >= k) return;
      const before = near(n, -1), after = near(n, 1);
      if (before && before.kind === 'vector' && after && after.kind === 'vector') b.kind = 'differs';
      else if ((n === firstReal && after && after.kind === 'vector') ||
               (n === lastReal && before && before.kind === 'vector')) b.kind = 'edge';
    });
    /* an unreadable run with plasmid on both sides is plasmid you could not
       read, not sequence you do not have */
    blocks.forEach((b, n) => {
      if (b.kind !== 'unreadable') return;
      const before = blocks[n - 1], after = blocks[n + 1];
      if (before && before.kind === 'vector' && after && after.kind === 'vector') b.inVectorRun = true;
    });

    const sum = kd => blocks.filter(b => b.kind === kd).reduce((t, b) => t + b.bp, 0);
    const edge = sum('edge');
    /* what you would be adding: the longest stretch the plasmid does not
       already have. Ambiguous runs inside it are carried along, because a
       bad patch in the middle of an insert does not split the insert — but
       it does mean the insert is not fully read, which is said plainly. */
    let span = null;
    blocks.forEach((b, n) => {
      if (b.kind !== 'new') return;
      let from = b.from, to = b.to, m = n + 1;
      while (m < blocks.length && (blocks[m].kind === 'new' ||
             (blocks[m].kind === 'unreadable' && blocks[m + 1] && blocks[m + 1].kind === 'new'))) {
        to = blocks[m].to; m++;
      }
      if (!span || to - from > span.bp) span = { from: from, to: to, bp: to - from };
    });
    if (span) {
      let amb = 0;
      for (let x = span.from; x < span.to; x++) if (r[x] === 'N') amb++;
      span.ambiguous = amb;
      span.seq = r.slice(span.from, span.to);
      /* Where the ambiguity sits matters more than how much of it there is.
         Sixty N scattered through a kilobase is a mediocre read; the same
         sixty packed into the last hundred and seventy bases means the end
         of the insert is not established at all, and a stop codon could be
         missing without anything looking wrong. */
      const firstN = span.seq.indexOf('N');
      span.cleanBp = firstN < 0 ? span.bp : firstN;
      /* What decides whether the far end is established is not where the
         ambiguity starts but whether the read recovers from it: a patch in
         the middle with clean sequence after it still fixes the end, while
         a read that never comes back does not, however little of it is N. */
      const lastN = span.seq.lastIndexOf('N');
      span.finalCleanRun = lastN < 0 ? span.bp : span.bp - lastN - 1;
      span.tailAmbiguous = amb > 0 && span.finalCleanRun < Math.max(30, span.bp * 0.15);
    }

    /* which of the plasmid's own features the matched parts land on, so the
       answer is "your vector already gives you the T7 promoter and the His
       tag" rather than a pair of coordinates */
    const names = [];
    (rec.features || []).forEach(f => {
      if (f.virtual) return;
      const pieces = E.featurePieces(f, L);
      const touched = blocks.some(b => b.kind === 'vector' &&
        pieces.some(p => {
          const a = Math.min(b.vecFrom, b.vecTo == null ? b.vecFrom : b.vecTo);
          const z = Math.max(b.vecFrom, b.vecTo == null ? b.vecFrom : b.vecTo);
          return a < p.end && z > p.start;
        }));
      if (touched && f.label && names.indexOf(f.label) < 0) names.push(f.label);
    });

    return {
      k: k, length: r.length, blocks: blocks,
      inVector: sum('vector'), isNew: sum('new'), unreadable: sum('unreadable'),
      differs: sum('differs'), edge: edge,
      newSpan: span, features: names,
      vectorName: rec.name || 'the construct', vectorLength: L, circular: circular
    };
  };

  if (typeof module !== 'undefined' && module.exports) module.exports = E;
  root.Bio = E;
})(typeof window !== 'undefined' ? window : globalThis);

