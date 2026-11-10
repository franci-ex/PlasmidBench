/* ============================================================================
   PlasmidBench — tolerant annotation.

   THE PROBLEM
   The matcher this project started with does exact substring search, on both
   strands, origin-aware, with a coincidence estimate. It is a good exact
   matcher. But a CMV promoter one base from the stored sequence is not
   found, and that is the common case, not the edge: the CMV promoter in any
   given plasmid is some lab's variant, a cloning scar has eaten three bases
   off a terminator, somebody's GFP is EGFP with two silent changes.
   pLannotate's advantage over exact matchers is not its library, it is
   BLAST — homology search, which tolerates that.

   BLAST does not run in a browser offline. So this is seed-and-extend,
   written to work on the scale a plasmid annotator actually needs: a few
   thousand library entries against a sequence of a few tens of kilobases.

   HOW IT WORKS

   1. INDEX. Every library entry is cut into overlapping k-mers (k=12 by
      default) and each is recorded against (entry, offset). Built once for
      a library and reused across sequences, because the library does not
      change between calls and rebuilding it per sequence is the difference
      between a tenth of a second and ten.

   2. SEED. The query is walked k-mer by k-mer. A hit gives a candidate
      (entry, diagonal): the offset between where the k-mer sits in the
      query and where it sits in the entry. Hits that agree on a diagonal
      are the same alignment seen several times, so they are grouped.

   3. EXTEND. Each diagonal group is extended ungapped in both directions
      while the running identity holds up, which settles the common case —
      point substitutions — with no alignment at all.

   4. ALIGN, in a band. If the ungapped extension stops short of the entry's
      ends, the remainder is aligned with Needleman-Wunsch restricted to a
      band around the diagonal. A full alignment of a 2 kb entry against a
      10 kb query is 20 million cells; in a band of 2w+1 it is 2 kb × 2w+1,
      which is thousands. The band is what makes this possible in a browser,
      and its width is the largest indel the alignment can contain — stated,
      not hidden, and reported when a match is cut off by it.

   5. REPORT. Each match carries where it is, which strand, the identity
      over the aligned part, and how much of the library entry was covered.
      Partial coverage is a finding in its own right rather than something
      to discard: a terminator present at 60% of its length is a cloning
      remnant, which is exactly what you want to be told about a plasmid
      whose history you are trying to work out.

   WHY NOT ONE BIG ALIGNMENT
   Because the question is not "how do these two sequences relate", it is
   "which of these few thousand short things is in this one long thing".
   Seeding answers the second cheaply and only pays for alignment where
   there is already evidence.

   WHAT THIS IS NOT
   It is not BLAST, and does not claim BLAST's statistics: there is no
   E-value here, because an E-value needs a database model this does not
   have. What it reports instead is identity, coverage, and the chance of a
   seed of that length occurring by coincidence in a sequence of this
   length — which is the quantity that actually matters when deciding
   whether a 12-base agreement means anything.
   ============================================================================ */
(function (root) {
  'use strict';
  const E = (typeof require !== 'undefined' && typeof module !== 'undefined')
    ? require('./engine.js') : root.Bio;
  if (typeof require !== 'undefined' && typeof module !== 'undefined') require('./engine7.js');

  const DEFAULTS = {
    k: 12,              /* seed length */
    minIdentity: 0.85,  /* over the aligned part */
    /* Of the library entry. Half rather than more, because a fragment is a
       finding in itself: twenty bases of a terminator left behind by a
       cloning step is what someone reconstructing a plasmid's history wants
       to see, and a floor of 0.6 threw away a feature trimmed ten bases at
       each end. */
    minCoverage: 0.50,
    band: 12,           /* the largest indel an alignment may contain */
    minLen: 15,         /* an entry shorter than this is left to the exact matcher */
    maxPerEntry: 8      /* how many places one entry may be reported in */
  };

  /**
   * The mismatch penalty that makes a local extension stop at a given
   * identity. Extension pays while the identity of what is added is above
   * m/(1+m), so m = id/(1-id): at 0.85, 5.67; at 0.95, 19. Clamped, because
   * at 0.99 the arithmetic asks for 99 and nothing is gained past about 20.
   */
  function mismatchFor(minIdentity) {
    const id = Math.min(0.98, Math.max(0.5, minIdentity || 0.85));
    return Math.max(1, Math.min(20, Math.round(id / (1 - id))));
  }
  E._mismatchFor = mismatchFor;

  /* --------------------------------------------------------------- index -- */

  /**
   * Build the k-mer index for a library. Keep it and pass it back in: it
   * does not depend on the sequence being annotated.
   *
   * @param {LibraryEntry[]} lib
   * @param {{k?:number, minLen?:number}} [opts]
   */
  E.buildFeatureIndex = function (lib, opts) {
    opts = Object.assign({}, DEFAULTS, opts || {});
    const k = opts.k;
    /** @type {Map<string, number[]>} kmer -> flat list of entry, offset, strand */
    const map = new Map();
    const entries = [];
    (lib || []).forEach(f => {
      const seq = E.cleanSeq(f.seq || '');
      if (seq.length < Math.max(k, opts.minLen)) return;
      const id = entries.length;
      const rc = E.revcomp(seq);
      /* Both strands are kept on the entry rather than recomputed for each
         candidate: with a few thousand entries and several candidate
         diagonals each, reverse-complementing in the inner loop was most of
         the running time. */
      entries.push({ name: f.name, type: f.type || 'misc_feature', seq: seq, rc: rc,
                     len: seq.length, source: f.source || null, note: f.note || null,
                     accession: f.accession || null });
      /* both strands are indexed rather than the query being searched twice:
         one pass over a long sequence instead of two */
      for (const [s, strand] of [[seq, 1], [rc, -1]]) {
        for (let i = 0; i + k <= s.length; i++) {
          const key = s.substr(i, k);
          if (/[^ACGT]/.test(key)) continue;      /* an ambiguous seed matches too much */
          let arr = map.get(key);
          if (!arr) { arr = []; map.set(key, arr); }
          arr.push(id, i, strand);
        }
      }
    });
    return { k: k, map: map, entries: entries, size: entries.length };
  };

  /* ------------------------------------------------------- banded aligner -- */

  /**
   * Smith-Waterman confined to a band around a known diagonal.
   *
   * Local at both ends on purpose. A library entry is often present in a
   * plasmid only in part — twenty bases of a terminator left behind by a
   * cloning step — and a global alignment of the whole entry would score
   * that as 40% identity and throw it away, when the finding is "most of
   * this terminator is missing", which is exactly what someone working out
   * a plasmid's history wants to be told. So the alignment is free to start
   * and end anywhere, and how much of the entry it covered is reported
   * separately from how well it matched.
   *
   * Confining it to a band is what makes it affordable: a full matrix for a
   * 2 kb entry against a 10 kb query is 20 million cells, and in a band of
   * 2w+1 it is 2 kb x 2w+1, a few tens of thousands. The band width is
   * therefore the largest indel an alignment can contain — stated here
   * rather than hidden, and reported on any hit that reached the edge of it.
   *
   * Scores are the plain ones that suit nucleotide matching at this
   * identity: +1 match, -1 mismatch, -2 to open a gap, -1 to extend. They
   * are stated rather than tuned, because tuning them against the test set
   * would be fitting the scoring to the answer.
   *
   * @param {string} entry      the library sequence, already on the right strand
   * @param {string} query      the window of the sequence being annotated
   * @param {number} diag       where entry[0] is expected to sit in query
   * @param {number} band       half-width: the largest indel allowed
   */
  function bandedAlign(entry, query, diag, band, mismatch) {
    const m = entry.length, n = query.length;
    /* The mismatch penalty is derived from the identity being asked for,
       not chosen. A local alignment with +1 match and -1 mismatch will
       always extend into the flanking plasmid, because two matches pay for
       one mismatch: a terminator present only in its middle 28 bases came
       back as 39 bases at 92% identity, which is both numbers wrong from
       one cause. Extension is profitable while the identity of what is
       being added exceeds m/(1+m), so setting m from the threshold makes
       the alignment stop exactly where the caller said it should: at 85%,
       m = 5.67, so 6. Gap costs follow the mismatch, or a gap would be
       cheaper than a mismatch and the aligner would invent indels. */
    const MATCH = 1;
    const MISMATCH = -(mismatch || 6);
    const GAP_OPEN = MISMATCH - 1, GAP_EXT = Math.min(-1, Math.round(MISMATCH / 2));
    const w = band, width = 2 * w + 1;
    /* j = i + diag + (c - w), so c = j - i - diag + w */
    const colOf = (i, j) => j - i - diag + w;
    const jOf = (i, c) => i + diag + (c - w);

    let prev = new Float64Array(width);          /* row i-1 */
    let prevUp = new Float64Array(width);        /* best ending in a gap in the query */
    const dir = new Uint8Array((m + 1) * width);
    const len = new Uint16Array((m + 1) * width);     /* aligned columns, for identity */
    const mat = new Uint16Array((m + 1) * width);     /* of which matches */
    const gap = new Uint16Array((m + 1) * width);     /* of which gaps */
    const cov = new Uint16Array((m + 1) * width);     /* entry bases consumed */

    let bestScore = 0, bestI = 0, bestC = 0;
    for (let i = 1; i <= m; i++) {
      const cur = new Float64Array(width);
      const curUp = new Float64Array(width);
      let left = 0, leftRun = 0;                 /* gap in the entry, within this row */
      for (let c = 0; c < width; c++) {
        const j = jOf(i, c);
        if (j < 1 || j > n) { cur[c] = 0; curUp[c] = 0; left = 0; continue; }

        /* diagonal: (i-1, j-1) is the same column c */
        const sub = entry[i - 1] === query[j - 1] ? MATCH : MISMATCH;
        const diagScore = prev[c] + sub;

        /* gap in the query, entry advancing alone: (i-1, j) is column c+1 */
        const up = c + 1 < width
          ? Math.max(prev[c + 1] + GAP_OPEN, prevUp[c + 1] + GAP_EXT)
          : -1e9;
        curUp[c] = Math.max(up, -1e9);

        /* gap in the entry, query advancing alone: (i, j-1) is column c-1 */
        const leftScore = c >= 1 ? Math.max(cur[c - 1] + GAP_OPEN, leftRun + GAP_EXT) : -1e9;

        let best = 0, from = 0;                  /* 0 = start afresh here */
        if (diagScore > best) { best = diagScore; from = 1; }
        if (up > best) { best = up; from = 2; }
        if (leftScore > best) { best = leftScore; from = 3; }
        cur[c] = best;
        leftRun = leftScore;
        left = leftScore;
        const at = i * width + c;
        dir[at] = from;
        if (from === 1) {
          const p = (i - 1) * width + c;
          len[at] = len[p] + 1; mat[at] = mat[p] + (sub === MATCH ? 1 : 0);
          gap[at] = gap[p]; cov[at] = cov[p] + 1;
        } else if (from === 2) {
          const p = (i - 1) * width + (c + 1);
          len[at] = len[p] + 1; mat[at] = mat[p]; gap[at] = gap[p] + 1; cov[at] = cov[p] + 1;
        } else if (from === 3) {
          const p = i * width + (c - 1);
          len[at] = len[p] + 1; mat[at] = mat[p]; gap[at] = gap[p] + 1; cov[at] = cov[p];
        } else { len[at] = 0; mat[at] = 0; gap[at] = 0; cov[at] = 0; }

        if (best > bestScore) { bestScore = best; bestI = i; bestC = c; }
      }
      prev = cur; prevUp = curUp;
      void left;
    }
    if (!bestScore) return null;

    const at = bestI * width + bestC;
    const aligned = len[at];
    if (!aligned) return null;
    /* walk back only far enough to find where it started in the query */
    let i = bestI, c = bestC;
    while (i > 0) {
      const f = dir[i * width + c];
      if (f === 0) break;
      if (f === 1) { i--; }
      else if (f === 2) { i--; c = c + 1; }
      else { c = c - 1; }
      if (c < 0 || c >= width) break;
    }
    return {
      score: bestScore,
      identity: mat[at] / aligned,
      matches: mat[at], aligned: aligned, gaps: gap[at],
      entryCovered: cov[at],
      entryFrom: i, entryTo: bestI,
      qStart: Math.max(0, jOf(i, c)),
      qEnd: jOf(bestI, bestC),
      clippedByBand: bestC === 0 || bestC === width - 1
    };
  }

  E._bandedAlign = bandedAlign;

  /* ------------------------------------------------------------ annotate -- */

  /**
   * Find library entries in a sequence, allowing for substitutions and
   * small indels.
   *
   * @param {string} seq
   * @param {{k?:number, map?:Map, entries?:Array}|Array} indexOrLib
   * @param {{circular?:boolean, minIdentity?:number, minCoverage?:number, band?:number, maxPerEntry?:number}} [opts]
   * @returns {{hits:Array, searched:number, seeds:number, ms:number}}
   */
  E.annotateApprox = function (seq, indexOrLib, opts) {
    opts = Object.assign({}, DEFAULTS, opts || {});
    const t0 = Date.now();
    const index = (indexOrLib && indexOrLib.map) ? indexOrLib
      : E.buildFeatureIndex(indexOrLib || E.FEATURE_LIB, opts);
    const k = index.k;
    const clean = E.cleanSeq(seq);
    const L = clean.length;
    if (!L || !index.size) return { hits: [], searched: index.size || 0, seeds: 0, ms: 0 };

    /* A circular sequence is searched with enough of its start appended that
       a feature crossing the origin is a contiguous string; the overhang is
       the longest entry, so nothing can be cut in half by the join. */
    const longest = index.entries.reduce((a, e) => Math.max(a, e.len), 0);
    const circular = !!opts.circular && L > longest;
    const hay = circular ? clean + clean.slice(0, Math.min(L, longest + opts.band)) : clean;

    /* 1-2: seed, and group hits by (entry, strand, diagonal) */
    /** @type {Map<string, {id:number, strand:number, diag:number, hits:number, qMin:number, qMax:number}>} */
    const groups = new Map();
    let seeds = 0;
    for (let q = 0; q + k <= hay.length; q++) {
      const key = hay.substr(q, k);
      const arr = index.map.get(key);
      if (!arr) continue;
      for (let a = 0; a < arr.length; a += 3) {
        const id = arr[a], off = arr[a + 1], strand = arr[a + 2];
        seeds++;
        const diag = q - off;
        const gk = id + ':' + strand + ':' + diag;
        let g = groups.get(gk);
        if (!g) {
          g = { id: id, strand: strand, diag: diag, hits: 0, qMin: q, qMax: q + k };
          groups.set(gk, g);
        }
        g.hits++;
        if (q < g.qMin) g.qMin = q;
        if (q + k > g.qMax) g.qMax = q + k;
      }
    }

    /* 3-4: extend and align each group */
    const perEntry = new Map();
    const hits = [];
    for (const g of groups.values()) {
      const e = index.entries[g.id];
      /* BLAST's two-hit rule, for the same reason: a single k-mer agreement
         is expected by chance several times in a plasmid, so extending
         every one of them is almost all wasted work. Two seeds on one
         diagonal is not a coincidence. An entry too short to produce two
         non-overlapping seeds is exempt, or an 18 bp promoter could never
         be found at all. */
      if (g.hits < 2 && e.len >= k * 2) continue;
      const entrySeq = g.strand === 1 ? e.seq : e.rc;
      /* where the entry would start in the query if the diagonal is right */
      const qStart = g.diag;
      if (qStart + e.len < 0 || qStart > hay.length) continue;

      /* UNGAPPED EXTENSION, from the seed outwards.
         Not a comparison of the whole entry against the whole window: that
         is what this did first, and it gets a truncated feature wrong in
         both directions at once. A terminator present only in its middle
         thirty bases was compared across all forty-eight, so the eighteen
         bases of surrounding plasmid went into the identity (pulling it
         down to 92%) and into the coverage (pushing it up to 79%) — one
         number too low and the other too high, from the same mistake.

         So the match is grown out from the seed instead, with BLAST's
         X-drop rule: keep a running score, remember where it was highest,
         and stop when it has fallen X below that. Where it was highest is
         where the feature ends, which is the quantity both identity and
         coverage have to be measured over. */
      /* The extension's penalty is not the aligner's, and the reason is
         worth stating because setting them equal made things worse in a way
         that looked like an improvement.

         The aligner's job is to not run off the end of a real match into
         the surrounding plasmid, so its penalty is tuned to the identity
         being asked for. The extension's job is the opposite: to cross a
         rough patch inside a real match and find where it truly ends. With
         the aligner's penalty (-6 at 85%) it stopped inside a uniformly
         10%-diverged feature, so starts and ends came back wrong and
         recall fell. With -3, a run of random flanking sequence costs about
         2 a base and the extension stops within four or five bases of the
         truth, while a 10%-diverged feature still pays +0.5 a base and is
         crossed. Identity is then measured over the extent found, and the
         caller's threshold does the filtering it is for. */
      const MM = 3;
      const X = 12;
      const seedI = g.qMin - qStart;              /* entry offset of the first seed */
      let i0 = Math.max(0, seedI), i1 = Math.min(e.len, g.qMax - qStart);
      {
        /* right */
        let score = 0, bestScore = 0, bestI = i1;
        for (let i = i1; i < e.len; i++) {
          const j = qStart + i;
          if (j < 0 || j >= hay.length) break;
          score += entrySeq[i] === hay[j] ? 1 : -MM;
          if (score > bestScore) { bestScore = score; bestI = i + 1; }
          else if (score < bestScore - X) break;
        }
        i1 = bestI;
        /* left */
        score = 0; bestScore = 0;
        let bestLeft = i0;
        for (let i = i0 - 1; i >= 0; i--) {
          const j = qStart + i;
          if (j < 0 || j >= hay.length) break;
          score += entrySeq[i] === hay[j] ? 1 : -MM;
          if (score > bestScore) { bestScore = score; bestLeft = i; }
          else if (score < bestScore - X) break;
        }
        i0 = bestLeft;
      }
      let matches = 0, compared = 0;
      for (let i = i0; i < i1; i++) {
        const j = qStart + i;
        if (j < 0 || j >= hay.length) continue;
        compared++;
        if (entrySeq[i] === hay[j]) matches++;
      }
      const ungappedIdentity = compared ? matches / compared : 0;
      const ungappedCoverage = compared / e.len;

      let best = { identity: ungappedIdentity, coverage: ungappedCoverage,
                   qStart: Math.max(0, qStart + i0), qEnd: Math.min(hay.length, qStart + i1),
                   gaps: 0, how: 'substitutions only', clippedByBand: false };

      /* The aligner is paid for whenever the ungapped extension stopped
         short of the entry's ends. That is what an insertion or deletion
         looks like from here — the extension runs into the indel, the
         X-drop stops it, and the rest of the entry is on a different
         diagonal that only a gapped alignment can reach.

         It is also what a genuine fragment looks like, and there the
         aligner finds the same fragment again and nothing is gained. Trying
         it in both cases rather than guessing between them is the right
         trade: the band bounds the cost to a few tens of thousands of
         cells, and guessing wrongly means reporting a 3 bp deletion as a
         42%-coverage fragment, which is a different and wrong finding. The
         better of the two results is kept, so the attempt can only help. */
      if (ungappedCoverage < 0.97) {
        const pad = opts.band + 2;
        const from = Math.max(0, qStart - pad);
        const to = Math.min(hay.length, qStart + e.len + pad);
        const al = bandedAlign(entrySeq, hay.slice(from, to), Math.max(0, qStart - from),
                               opts.band, mismatchFor(opts.minIdentity));
        if (al && (al.entryCovered / e.len) > best.coverage && al.identity >= opts.minIdentity) {
          best = { identity: al.identity,
                   coverage: al.entryCovered / e.len,
                   qStart: from + al.qStart, qEnd: from + al.qEnd,
                   gaps: al.gaps, how: al.gaps ? al.gaps + ' base indel' : 'substitutions only',
                   clippedByBand: al.clippedByBand };
        }
      }

      if (best.identity < opts.minIdentity) continue;
      if (best.coverage < opts.minCoverage) continue;

      /* Kept in the coordinates of the searched string, where every span is
         a plain interval. Unwrapping to the sequence's own numbering turns a
         span that crosses the origin into one whose end is before its
         start, and comparing those as intervals — which the overlap check
         below does — gives nonsense: a wrapped hit and the half of it that
         also matched on its own were both reported, because the arithmetic
         that should have collapsed them was working with a negative width.
         So unwrapping happens last, after everything that compares spans. */
      const hayStart = best.qStart, hayEnd = best.qEnd;
      if (hayStart >= L) continue;               /* entirely in the appended copy */

      const exact = best.identity === 1 && best.coverage === 1;
      const hit = {
        name: e.name, type: e.type, strand: g.strand,
        hayStart: hayStart, hayEnd: hayEnd,
        start: hayStart, end: hayEnd, wrapped: false,
        identity: best.identity, coverage: best.coverage,
        gaps: best.gaps, seeds: g.hits, exact: exact,
        partial: best.coverage < 0.98,
        how: exact ? 'exact' : best.how,
        clippedByBand: best.clippedByBand,
        source: e.source, accession: e.accession,
        /* The chance of the seed agreement alone being coincidence, which is
           the honest substitute for an E-value: with k=12 in a 5 kb
           sequence that is about 0.6, so a single seed means nothing and
           several on one diagonal mean a great deal. */
        seedChance: E.chanceHits('N'.repeat(0) + e.seq.substr(0, k), L),
        libLength: e.len
      };

      /* one entry may be in a plasmid several times, but not hundreds */
      const seenN = perEntry.get(e.name) || 0;
      if (seenN >= opts.maxPerEntry) continue;
      perEntry.set(e.name, seenN + 1);
      hits.push(hit);
    }

    /* Overlapping reports of the same thing: keep the better. Two entries
       for GFP variants will both match; reporting both at the same place
       with the same type is noise, and the one with the higher identity is
       the answer. */
    hits.sort((a, b) => (b.identity - a.identity) || (b.coverage - a.coverage) ||
                        (b.end - b.start) - (a.end - a.start));
    /* Overlap has to be judged on the circle, not on the searched string.
       A feature split across base 1 is found twice: once as its second half
       at the start of the sequence, and once whole, in the appended copy
       that crosses the join. Those two spans do not overlap as intervals —
       0..24 and 924..972 — while on the molecule the first is entirely
       inside the second. Reducing both to the circle before comparing is
       what collapses them, and leaves the whole match rather than a half of
       it. */
    const circSpans = h => {
      const s = ((h.hayStart % L) + L) % L;
      const n = h.hayEnd - h.hayStart;
      return s + n <= L ? [[s, s + n]] : [[s, L], [0, s + n - L]];
    };
    const circOverlap = (a, b) => {
      let n = 0;
      for (const [as, ae] of circSpans(a)) for (const [bs, be] of circSpans(b))
        n += Math.max(0, Math.min(ae, be) - Math.max(as, bs));
      return n;
    };
    const kept = [];
    for (const h of hits) {
      const clash = kept.find(x => x.strand === h.strand &&
        circOverlap(x, h) > 0.8 * Math.min(x.hayEnd - x.hayStart, h.hayEnd - h.hayStart));
      if (clash) {
        /* the same thing seen twice, or two library entries for the same
           feature; the better match is the answer and the other is a note
           on it. Types are not compared: the half of a wrapped match found
           on its own has the same type, and so does a GFP variant. */
        if (clash.name !== h.name) (clash.alsoMatched = clash.alsoMatched || []).push(h.name);
        continue;
      }
      kept.push(h);
    }
    /* now, and only now, into the sequence's own numbering */
    kept.forEach(h => {
      h.start = h.hayStart;
      h.end = h.hayEnd;
      if (h.end > L) { h.end -= L; h.wrapped = true; }
      delete h.hayStart; delete h.hayEnd;
    });
    kept.sort((a, b) => a.start - b.start);
    return { hits: kept, searched: index.size, seeds: seeds, ms: Date.now() - t0 };
  };

  /* ------------------------------------------------- library importing -- */

  /**
   * Read a feature library from a file.
   *
   * Section 10 of the project evaluation is firm that the thing to do is
   * import a curated library rather than build one, and firm about which:
   * ApE ships a plain-text, user-editable feature list in a format close to
   * this project's own, FPbase has an open API for fluorescent proteins,
   * and GenoLIB is SBOL. Crawling Addgene would circumvent an access
   * control; their Developers Portal is free and takes about a week.
   *
   * So three shapes are read here, and each entry keeps where it came from.
   * Provenance is not decoration: a library assembled from several sources
   * has different licence terms and different reliability per entry, and a
   * feature called in a plasmid map should be traceable to the thing that
   * called it.
   *
   *   ApE           tab-separated: name, sequence, type, colour, ...
   *   FPbase        JSON array of proteins with a dna sequence
   *   PlasmidBench  JSON {format, source, entries:[{name,type,seq,...}]}
   *
   * @param {string} text
   * @param {{source?:string, licence?:string}} [meta] what to record as the origin
   */
  E.parseFeatureLibrary = function (text, meta) {
    meta = meta || {};
    const out = [];
    const note = { source: meta.source || 'imported', licence: meta.licence || null };
    const t = String(text || '').trim();
    if (!t) return { ok: false, error: 'That file is empty.' };

    /* --- JSON, either shape --- */
    if (t[0] === '{' || t[0] === '[') {
      let j;
      try { j = JSON.parse(t); }
      catch (e) { return { ok: false, error: 'That looks like JSON but does not parse: ' + e.message }; }

      /* this project's own */
      if (j && j.format === 'plasmidbench-feature-library' && Array.isArray(j.entries)) {
        j.entries.forEach(f => {
          const seq = E.cleanSeq(f.seq || f.sequence || '');
          if (seq.length < 8) return;
          out.push({ name: f.name || 'unnamed', type: f.type || 'misc_feature', seq: seq,
                     source: f.source || j.source || note.source,
                     licence: f.licence || j.licence || note.licence,
                     accession: f.accession || null, note: f.note || null });
        });
        return finish(out, j.source || note.source, j.licence || note.licence);
      }

      /* FPbase: an array of proteins, each with a DNA sequence somewhere */
      const arr = Array.isArray(j) ? j : (j.results || j.proteins || j.data || []);
      if (Array.isArray(arr) && arr.length) {
        arr.forEach(p => {
          const seq = E.cleanSeq(p.dna_sequence || p.dna || p.seq ||
            (p.default_state && p.default_state.dna) || '');
          if (seq.length < 8) return;
          out.push({ name: p.name || p.slug || 'protein', type: 'CDS', seq: seq,
                     source: note.source || 'FPbase', licence: note.licence,
                     accession: p.genbank || p.uuid || p.slug || null,
                     note: p.agg || p.switch_type ? ('FPbase: ' + (p.agg || p.switch_type)) : null });
        });
        if (out.length) return finish(out, note.source || 'FPbase', note.licence);
      }
      return { ok: false, error: 'That JSON has no feature entries this reads. It understands ' +
               'a PlasmidBench feature library and an FPbase protein list.' };
    }

    /* --- ApE's feature file: tab-separated, one feature a line --- */
    const lines = t.split(/\r?\n/);
    let bad = 0;
    lines.forEach(line => {
      if (!line.trim() || line[0] === '#' || line[0] === ';') return;
      const col = line.split('\t');
      if (col.length < 2) { bad++; return; }
      /* ApE writes name, sequence, then type and colour; some copies have
         the first two the other way round, so whichever column looks like
         DNA is the sequence */
      let name = col[0].trim(), seq = E.cleanSeq(col[1]);
      if (!seq.length || seq.length < col[1].replace(/\s/g, '').length * 0.8) {
        const alt = E.cleanSeq(col[0]);
        if (alt.length >= col[0].replace(/\s/g, '').length * 0.8 && alt.length >= 8) {
          seq = alt; name = col[1].trim();
        }
      }
      if (seq.length < 8) { bad++; return; }
      out.push({ name: name || 'unnamed', type: (col[2] || '').trim() || 'misc_feature',
                 seq: seq, source: note.source, licence: note.licence,
                 accession: null, note: null });
    });
    if (!out.length) {
      return { ok: false, error: 'No features could be read from that. A tab-separated file ' +
               'wants a name and a sequence in its first two columns.' };
    }
    return finish(out, note.source, note.licence, bad);

    function finish(entries, source, licence, skipped) {
      /* Duplicates are the norm in a library assembled from several places,
         and the same sequence under two names is one feature with two names
         rather than two features. */
      const bySeq = new Map();
      entries.forEach(f => {
        const key = f.seq;
        const have = bySeq.get(key);
        if (!have) { bySeq.set(key, f); return; }
        if (have.name !== f.name) (have.aliases = have.aliases || []).push(f.name);
      });
      const uniq = Array.from(bySeq.values());
      const lens = uniq.map(f => f.seq.length).sort((a, b) => a - b);
      return {
        ok: true, entries: uniq,
        source: source, licence: licence,
        count: uniq.length, duplicates: entries.length - uniq.length,
        skipped: skipped || 0,
        shortest: lens[0] || 0, longest: lens[lens.length - 1] || 0,
        median: lens[Math.floor(lens.length / 2)] || 0,
        /* How many of them are short enough to turn up by chance in a
           plasmid, which is what a person needs to know before trusting a
           library: an 8-base "feature" is found in everything. */
        coincidental: uniq.filter(f => E.chanceHits(f.seq, 5000) > 0.05).length
      };
    }
  };

  /** Write a library back out, provenance and all. */
  E.featureLibraryToJSON = function (entries, meta) {
    meta = meta || {};
    return JSON.stringify({
      format: 'plasmidbench-feature-library',
      version: 1,
      source: meta.source || null,
      licence: meta.licence || null,
      exported: new Date().toISOString().slice(0, 10),
      note: 'Each entry keeps the source it came from. Check the licence of each ' +
            'source before redistributing a library assembled from several.',
      entries: (entries || []).map(f => ({
        name: f.name, type: f.type, seq: f.seq,
        source: f.source || null, licence: f.licence || null,
        accession: f.accession || null, note: f.note || null,
        aliases: f.aliases || undefined
      }))
    }, null, 1);
  };

  /* --------------------------------------------- Aho-Corasick exact pass -- */

  /**
   * Exact multi-pattern search in one pass over the sequence.
   *
   * The per-entry scan the exact matcher does is fine for 27 entries and is
   * 27 passes over the sequence; at a few thousand it is a few thousand
   * passes, and section 10 of the evaluation says so. This builds one
   * automaton and makes one pass, whatever the library's size.
   *
   * Both strands: the reverse complement of each pattern is a pattern too,
   * which keeps it to the single pass.
   */
  E.buildExactAutomaton = function (lib) {
    const goto_ = [new Map()], out = [[]], fail = [0];
    const add = (pat, meta) => {
      let s = 0;
      for (const ch of pat) {
        let nxt = goto_[s].get(ch);
        if (nxt == null) {
          nxt = goto_.length;
          goto_.push(new Map()); out.push([]); fail.push(0);
          goto_[s].set(ch, nxt);
        }
        s = nxt;
      }
      out[s].push(meta);
    };
    (lib || []).forEach((f, i) => {
      const seq = E.cleanSeq(f.seq || '');
      if (seq.length < 6 || /[^ACGT]/.test(seq)) return;
      add(seq, { i: i, len: seq.length, strand: 1, f: f });
      const rc = E.revcomp(seq);
      if (rc !== seq) add(rc, { i: i, len: rc.length, strand: -1, f: f });
    });
    /* failure links, breadth first */
    const queue = [];
    for (const [, s] of goto_[0]) { fail[s] = 0; queue.push(s); }
    while (queue.length) {
      const s = queue.shift();
      for (const [ch, t] of goto_[s]) {
        queue.push(t);
        let f = fail[s];
        while (f && !goto_[f].has(ch)) f = fail[f];
        fail[t] = goto_[f].has(ch) ? goto_[f].get(ch) : 0;
        if (fail[t] === t) fail[t] = 0;
        out[t] = out[t].concat(out[fail[t]]);
      }
    }
    return { goto_: goto_, out: out, fail: fail, states: goto_.length };
  };

  /** Run the automaton. Returns [{ f, start, end, strand }]. */
  E.runExactAutomaton = function (auto, seq) {
    const found = [];
    let s = 0;
    for (let i = 0; i < seq.length; i++) {
      const ch = seq[i];
      while (s && !auto.goto_[s].has(ch)) s = auto.fail[s];
      s = auto.goto_[s].has(ch) ? auto.goto_[s].get(ch) : 0;
      const hits = auto.out[s];
      for (let h = 0; h < hits.length; h++) {
        const m = hits[h];
        found.push({ f: m.f, start: i - m.len + 1, end: i + 1, strand: m.strand });
      }
    }
    return found;
  };

  if (typeof module !== 'undefined' && module.exports) module.exports = E;
  else root.Bio = E;
})(typeof window !== 'undefined' ? window : globalThis);
