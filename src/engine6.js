/* ===========================================================================
   ALIGNMENT FIGURES, IN THE ESPRIPT IDIOM

   An alignment you can read on screen and an alignment you can put in a
   paper are not the same drawing. The second one has conventions, and they
   have been the same since ESPript (Robert & Gouet, Nucleic Acids Research
   2014, 42:W320) made them standard:

     - a strictly conserved column: white letter on a red box
     - a column similar above a threshold: red letter, and a blue frame
       drawn around the whole run of such columns
     - everything else: black letter, no box
     - a numbering ruler every ten residues on a reference sequence
     - secondary structure drawn above the alignment: a squiggle for a
       helix, an arrow for a strand, TT for a turn
     - a consensus line underneath, in MultAlin's notation

   WHAT IS IMPLEMENTED EXACTLY, AND WHAT IS NOT

   ESPript's own user guide specifies three column scores completely, and
   those three are implemented here to the letter:

     strict     the fraction of the column holding the same residue
     multalin   MultAlin similarity over the groups IV, LM, FY, NDQEBZ
     physchem   equivalence by physico-chemical property, over the groups
                HKR, DE, STNQ, AVLIM, FYW, PGC

   Its default mode is a "global score" over the Risler matrix, and the
   guide does not publish the formula for it. So it is NOT reproduced here
   and nothing claims to be it. In its place is a fourth mode, `blosum62`,
   defined in the open: the fraction of residue pairs in the column that
   score above zero in BLOSUM62 — the same notion of "similar" that a
   BLAST alignment marks with a + between its two rows. The default mode
   is `physchem`, because it is both fully specified and the one most
   figures actually want.

   The matrix below was dumped from Biopython's own BLOSUM62 rather than
   typed out, and the test suite compares all 576 cells back against it.

   N-ESPript's nucleotide conventions are the same drawing with one extra
   option: purines red, pyrimidines blue.
   ========================================================================= */
(function (root) {
  'use strict';
  const E = (typeof require !== 'undefined' && typeof module !== 'undefined')
    ? require('./engine.js') : root.Bio;
  if (typeof require !== 'undefined' && typeof module !== 'undefined') require('./engine5.js');

  /* ------------------------------------------------------------------
     BLOSUM62, as Biopython ships it. Upper triangle, row-major, over the
     alphabet below; the matrix is symmetric so the lower half is read
     back off the same numbers.
     ------------------------------------------------------------------ */
  const B62_ALPHA = 'ARNDCQEGHILKMFPSTWYVBZX*';
  const B62_PACKED = '4,-1,-2,-2,0,-1,-1,0,-2,-1,-1,-1,-1,-2,-1,1,0,-3,-2,0,-2,-1,0,-4,5,0,-2,-3,1,0,-2,0,-3,-2,2,-1,-3,-2,-1,-1,-3,-2,-3,-1,0,-1,-4,6,1,-3,0,0,0,1,-3,-3,0,-2,-3,-2,1,0,-4,-2,-3,3,0,-1,-4,6,-3,0,2,-1,-1,-3,-4,-1,-3,-3,-1,0,-1,-4,-3,-3,4,1,-1,-4,9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1,-3,-3,-2,-4,5,2,-2,0,-3,-2,1,0,-3,-1,0,-1,-2,-1,-2,0,3,-1,-4,5,-2,0,-3,-3,1,-2,-3,-1,0,-1,-3,-2,-2,1,4,-1,-4,6,-2,-4,-4,-2,-3,-3,-2,0,-2,-2,-3,-3,-1,-2,-1,-4,8,-3,-3,-1,-2,-1,-2,-1,-2,-2,2,-3,0,0,-1,-4,4,2,-3,1,0,-3,-2,-1,-3,-1,3,-3,-3,-1,-4,4,-2,2,0,-3,-2,-1,-2,-1,1,-4,-3,-1,-4,5,-1,-3,-1,0,-1,-3,-2,-2,0,1,-1,-4,5,0,-2,-1,-1,-1,-1,1,-3,-1,-1,-4,6,-4,-2,-2,1,3,-1,-3,-3,-1,-4,7,-1,-1,-4,-3,-2,-2,-1,-2,-4,4,1,-3,-2,-2,0,0,0,-4,5,-2,-2,0,-1,-1,0,-4,11,2,-3,-4,-3,-2,-4,7,-1,-3,-2,-1,-4,4,-3,-2,-1,-4,4,1,-1,-4,4,-1,-4,-1,-4,1';
  const B62 = (function () {
    const n = B62_ALPHA.length, v = B62_PACKED.split(',').map(Number), m = {};
    let k = 0;
    for (let i = 0; i < n; i++) {
      for (let j = i; j < n; j++) {
        const s = v[k++];
        m[B62_ALPHA[i] + B62_ALPHA[j]] = s;
        m[B62_ALPHA[j] + B62_ALPHA[i]] = s;
      }
    }
    return m;
  })();
  E.BLOSUM62 = B62;
  E.blosum62 = function (a, b) {
    const s = B62[String(a).toUpperCase() + String(b).toUpperCase()];
    return s === undefined ? null : s;
  };

  /* A BLOSUM62 scorer for the pairwise dynamic program in engine.js, and
     the protein defaults that go with it. The DP uses a linear gap
     penalty, so -8 rather than BLAST's affine -11/-1: a linear -11 would
     refuse to open the gaps a real alignment needs.                      */
  E.proteinScorer = function () {
    return function (a, b) {
      const s = B62[String(a).toUpperCase() + String(b).toUpperCase()];
      return s === undefined ? -4 : s;
    };
  };
  E.PROTEIN_ALIGN = { gap: -8, scorer: null };
  E.alignProtein = function (a, b, opts) {
    return E.align(a, b, Object.assign({ gap: -8 }, opts || {}, { scorer: E.proteinScorer() }));
  };
  E.msaProtein = function (seqs, opts) {
    return E.msa(seqs, Object.assign({ gap: -8 }, opts || {}, { scorer: E.proteinScorer() }));
  };

  /* Is this protein or nucleotide? Letters that are not IUPAC nucleotide
     codes settle it; a sequence written only in A, C, G, T, U and N is
     read as nucleotide, which is the right way round — a peptide of only
     those residues is vanishingly rare and a DNA sequence of them is the
     normal case. */
  E.looksProteinSeq = function (s) {
    const up = String(s || '').toUpperCase().replace(/[^A-Z*]/g, '');
    if (!up.length) return false;
    if (/[EFILPQJOZX*]/.test(up)) return true;
    const acgtun = (up.match(/[ACGTUN]/g) || []).length;
    return acgtun / up.length < 0.85;
  };

  /* ------------------------------------------------------------------
     The two published equivalence groupings. A residue not in any group
     is only ever equivalent to itself.
     ------------------------------------------------------------------ */
  const GROUPS = {
    /* MultAlin's similarity groups, as ESPript's guide lists them */
    multalin: ['IV', 'LM', 'FY', 'NDQEBZ'],
    /* equivalence by physico-chemical property */
    physchem: ['HKR', 'DE', 'STNQ', 'AVLIM', 'FYW', 'PGC'],
    /* The nucleotide analogue, and the one N-ESPript's "colour by
       chemical properties" expresses: purines together, pyrimidines
       together. A column of A and G is a transition, which is the
       substitution a real alignment sees most of; a column of A and C is
       a transversion and is not shaded. */
    purine: ['AGR', 'CTUY']
  };
  function groupMap(list) {
    const m = {};
    list.forEach((g, i) => { for (const ch of g) m[ch] = i; });
    return m;
  }
  const GMAP = {
    multalin: groupMap(GROUPS.multalin),
    physchem: groupMap(GROUPS.physchem),
    purine: groupMap(GROUPS.purine)
  };
  E.ALN_GROUPS = GROUPS;
  /* which scoring modes make sense for which kind of sequence: the
     protein groupings say nothing about a nucleotide column, and the
     purine grouping says nothing about a residue */
  E.ALN_MODES = {
    protein: ['physchem', 'multalin', 'strict', 'blosum62'],
    nucleotide: ['strict', 'purine']
  };

  /* ------------------------------------------------------------------
     One column's score, on [0, 1].

     A gap is never equivalent to a residue and never to another gap: a
     column of six residues and two gaps is six-eighths conserved at best,
     which is what you want a figure to say. `strictOnly` counts gaps as
     breaking strict identity outright, which is also what ESPript does —
     a column is strictly conserved only when every sequence has the same
     residue in it.
     ------------------------------------------------------------------ */
  function columnScore(col, mode) {
    const n = col.length;
    if (!n) return 0;
    const res = col.filter(c => c && c !== '-' && c !== '.');
    if (res.length < 2) return 0;
    if (mode === 'strict') {
      const first = res[0];
      let same = 0;
      for (const c of res) if (c === first) same++;
      return same / n;
    }
    if (mode === 'blosum62') {
      /* The fraction of residue pairs in the column scoring above zero in
         BLOSUM62 — the same notion of "similar" that a BLAST alignment
         marks with a + between the two rows. A column of one residue
         scores 1; a column of chemically unrelated residues scores 0.
         Pairs involving a gap count against it, as in every other mode.   */
      let good = 0, pairs = 0;
      for (let i = 0; i < n; i++) {
        for (let j = i + 1; j < n; j++) {
          pairs++;
          const a = col[i], b = col[j];
          if (!a || !b || a === '-' || a === '.' || b === '-' || b === '.') continue;
          const s = B62[a + b];
          if (s !== undefined && s > 0) good++;
        }
      }
      return pairs ? good / pairs : 0;
    }
    /* multalin and physchem: the largest set of residues in one group */
    const gm = GMAP[mode] || GMAP.physchem;
    const counts = {};
    let best = 0;
    for (const c of res) {
      const g = gm[c] !== undefined ? 'g' + gm[c] : 'r' + c;
      counts[g] = (counts[g] || 0) + 1;
      if (counts[g] > best) best = counts[g];
    }
    return best / n;
  }
  E.alnColumnScore = columnScore;

  function strictlyConserved(col) {
    const first = col[0];
    if (!first || first === '-' || first === '.') return false;
    for (const c of col) if (c !== first) return false;
    return true;
  }

  /* ------------------------------------------------------------------
     The consensus line, in MultAlin's notation as ESPript's guide gives
     it: upper case for strict identity, lower case where one residue
     holds more than half the column, and the four symbols for the
     classic equivalence pairs.
     ------------------------------------------------------------------ */
  const CONSENSUS_SYMBOL = [['IV', '!'], ['LM', '$'], ['FY', '%'], ['NDQEBZ', '#']];
  function consensusLine(cols) {
    let out = '';
    for (const col of cols) {
      if (strictlyConserved(col)) { out += col[0]; continue; }
      const res = col.filter(c => c && c !== '-' && c !== '.');
      if (!res.length) { out += ' '; continue; }
      /* all of the column inside one of the four classic sets */
      let sym = '';
      for (const [set, s] of CONSENSUS_SYMBOL) {
        if (res.length === col.length && res.every(c => set.indexOf(c) !== -1)) { sym = s; break; }
      }
      if (sym) { out += sym; continue; }
      const counts = {};
      let bestC = '', bestN = 0;
      for (const c of res) { counts[c] = (counts[c] || 0) + 1; if (counts[c] > bestN) { bestN = counts[c]; bestC = c; } }
      out += bestN > col.length / 2 ? bestC.toLowerCase() : ' ';
    }
    return out;
  }
  E.alnConsensus = consensusLine;

  /* ------------------------------------------------------------------
     Secondary structure, given as one character per column of the
     alignment (or per residue of the reference row, with gaps inserted
     for us). DSSP's own letters are understood, which means a DSSP file,
     a PDB-derived string or a hand-typed one all work:

       H  alpha helix        G  3-10 helix       I  pi helix
       E  beta strand        B  beta bridge
       T  turn               S  bend             - / space  none
     ------------------------------------------------------------------ */
  const SS_KIND = { H: 'helix', G: 'helix310', I: 'helixPi', E: 'strand', B: 'strand', T: 'turn', S: 'bend' };
  /* N-ESPript draws a nucleic-acid structure in dot-bracket notation, so
     that is understood too: a paired base either way round is drawn as a
     stem, an unpaired one in a loop region is left blank. Detected by the
     characters themselves, because no DSSP string contains a bracket. */
  const DOT_BRACKET = { '(': 'stem5', ')': 'stem3', '[': 'stem5', ']': 'stem3',
                        '{': 'stem5', '}': 'stem3', '<': 'stem5', '>': 'stem3' };
  function parseSS(text, width, row) {
    const raw = String(text || '').replace(/[\r\n]/g, '');
    if (!raw.trim()) return null;
    let s = raw;
    /* given per residue of a row, the gaps of that row are opened up */
    if (row && raw.replace(/[-\s.]/g, '').length && raw.length !== width) {
      let out = '', k = 0;
      for (let i = 0; i < row.length; i++) {
        if (row[i] === '-' || row[i] === '.') out += ' ';
        else out += (k < raw.length ? raw[k++] : ' ');
      }
      s = out;
    }
    while (s.length < width) s += ' ';
    const dotBracket = /[()\[\]{}<>]/.test(raw);
    const kinds = [];
    for (let i = 0; i < width; i++) {
      kinds.push(dotBracket ? (DOT_BRACKET[s[i]] || null) : (SS_KIND[s[i]] || null));
    }
    /* contiguous runs, which is what gets drawn */
    const runs = [];
    let i = 0;
    while (i < width) {
      if (!kinds[i]) { i++; continue; }
      let j = i;
      while (j + 1 < width && kinds[j + 1] === kinds[i]) j++;
      runs.push({ kind: kinds[i], start: i, end: j + 1 });
      i = j + 1;
    }
    return { kinds: kinds, runs: runs, dotBracket: dotBracket };
  }
  E.parseSecondaryStructure = parseSS;

  /* ======================================================================
     WHY THERE IS NO SECONDARY-STRUCTURE PREDICTOR HERE

     Chou-Fasman (1974) is the only prediction that can be made from
     composition alone, with no database search, and it was implemented
     here and then taken out again. On myoglobin, which is all helix and
     has no beta strand at all, it predicted 36% strand. On GFP, a beta
     barrel with long loops, it predicted 36% helix and 2% coil. Those
     are not marginal errors; a track like that drawn above an alignment
     in a paper is worse than no track.

     A predictor worth the name (PSIPRED, JPred, the structure-prediction
     networks) works from a position-specific scoring matrix built by
     searching a sequence database, and a page with no network cannot
     build one. So nothing is predicted here. What the figure does
     instead is take a real assignment in whatever form it arrives:

       - a DSSP string, one letter a residue (H G I E B T S, space)
       - the H/E/C output of PSIPRED or JPred, pasted as it comes,
         including PSIPRED's horizontal blocks with their Conf/Pred/AA
         lines
       - dot-bracket notation, for nucleic acids

     Parsing those exactly is something that can be tested. Guessing the
     structure is not.
     ====================================================================== */
  E.SS_SOURCES = ['DSSP', 'PSIPRED', 'JPred', 'dot-bracket'];

  /* Pull a per-residue string out of whatever was pasted. Returns
     { text, source } so the caller can say where it came from, or null. */
  E.parseSSInput = function (raw) {
    const t = String(raw == null ? '' : raw).replace(/\r\n?/g, '\n');
    if (!t.trim()) return null;
    /* PSIPRED's horizontal format: repeated Conf/Pred/AA blocks */
    const pred = t.match(/^\s*Pred:\s*(\S+)\s*$/gm);
    if (pred && pred.length) {
      const seq = pred.map(l => (/^\s*Pred:\s*(\S+)/.exec(l) || [, ''])[1]).join('');
      if (seq.length) return { text: seq.replace(/C/g, ' '), source: 'PSIPRED' };
    }
    /* JPred's jnetpred line, comma-separated */
    const jnet = /^\s*jnetpred\s*[:,]?\s*(.+)$/mi.exec(t);
    if (jnet) {
      const seq = jnet[1].replace(/[\s,]/g, '');
      if (/^[HEC\-]+$/i.test(seq) && seq.length > 3) {
        return { text: seq.toUpperCase().replace(/[C\-]/g, ' '), source: 'JPred' };
      }
    }
    const flat = t.replace(/\n/g, '');
    if (/[()\[\]{}<>]/.test(flat)) return { text: flat, source: 'dot-bracket' };
    /* a bare per-residue string. C means coil in every predictor's output
       and is not a DSSP code, so it is read as coil rather than as a
       structure element. */
    if (/^[HGIEBTSC\-. ]+$/i.test(flat)) {
      const up = flat.toUpperCase();
      const src = /C/.test(up) && !/[GIBS]/.test(up) ? 'H/E/C' : 'DSSP';
      return { text: up.replace(/[C\-.]/g, ' '), source: src };
    }
    return null;
  };

  /* ------------------------------------------------------------------
     Everything the figure needs, computed once.
     ------------------------------------------------------------------ */
  E.alnAnalyse = function (rows, opts) {
    opts = opts || {};
    const mode = opts.mode || 'physchem';
    const threshold = opts.threshold != null ? opts.threshold : 0.7;
    const width = rows.reduce((w, r) => Math.max(w, r.seq.length), 0);
    const seqs = rows.map(r => {
      let s = String(r.seq || '').toUpperCase();
      while (s.length < width) s += '-';
      return s;
    });
    const cols = [];
    for (let i = 0; i < width; i++) cols.push(seqs.map(s => s[i]));
    const scores = cols.map(c => columnScore(c, mode));
    const strict = cols.map(strictlyConserved);
    /* a frame is drawn around a run of columns at or above the threshold,
       strictly conserved ones included — that is what makes an ESPript
       figure readable as blocks rather than as speckle */
    const over = scores.map((s, i) => strict[i] || s >= threshold);
    const frames = [];
    let i = 0;
    while (i < width) {
      if (!over[i]) { i++; continue; }
      let j = i;
      while (j + 1 < width && over[j + 1]) j++;
      frames.push({ start: i, end: j + 1 });
      i = j + 1;
    }
    return {
      rows: rows.map((r, k) => ({ name: r.name || ('seq' + (k + 1)), seq: seqs[k] })),
      width: width, mode: mode, threshold: threshold,
      scores: scores, strict: strict, over: over, frames: frames,
      consensus: consensusLine(cols),
      identity: width ? strict.filter(Boolean).length / width * 100 : 0
    };
  };

  /* ------------------------------------------------------------------
     The drawing itself, as standalone SVG: real text, so every label is
     still editable in Illustrator or Inkscape.
     ------------------------------------------------------------------ */
  const PURINE = 'AG', PYRIMIDINE = 'CTU';
  function esc(s) {
    return String(s == null ? '' : s).replace(/&/g, '&amp;').replace(/</g, '&lt;')
      .replace(/>/g, '&gt;').replace(/"/g, '&quot;');
  }

  E.ALN_FIGURE_DEFAULTS = {
    mode: 'physchem', threshold: 0.7, columns: 60, cell: 11.5, lineH: 15,
    consensus: true, ruler: true, mono: false, nucleotideColour: false,
    molecule: 'auto', fontSize: 11, nameFont: 11, title: '', legend: true
  };
  E.MODE_LEGEND = {
    strict: 'the same residue in every sequence',
    multalin: 'one MultAlin group (IV · LM · FY · NDQEBZ)',
    physchem: 'one physico-chemical group (HKR · DE · STNQ · AVLIM · FYW · PGC)',
    blosum62: 'pairs scoring above zero in BLOSUM62',
    purine: 'both purines, or both pyrimidines'
  };

  E.alnFigureSVG = function (rows, opts) {
    opts = Object.assign({}, E.ALN_FIGURE_DEFAULTS, opts || {});
    if (!rows || rows.length < 2) throw new Error('An alignment figure needs at least two sequences.');
    /* protein, DNA or RNA. RNA is DNA with the thymines shown as uracils,
       the way N-ESPript does it; the scoring never sees the difference
       because U and T are in the same group. */
    let molecule = opts.molecule;
    if (!molecule || molecule === 'auto') {
      molecule = rows.some(r => E.looksProteinSeq(String(r.seq).replace(/[-.]/g, ''))) ? 'protein' : 'dna';
    }
    if (molecule === 'rna') {
      rows = rows.map(r => ({ name: r.name, seq: String(r.seq).toUpperCase().replace(/T/g, 'U') }));
    }
    /* a scoring mode that means nothing for this kind of sequence falls
       back to strict identity rather than shading the figure wrongly */
    const allowed = E.ALN_MODES[molecule === 'protein' ? 'protein' : 'nucleotide'];
    if (allowed.indexOf(opts.mode) === -1) opts = Object.assign({}, opts, { mode: 'strict' });
    const A = E.alnAnalyse(rows, opts);
    A.molecule = molecule;
    const per = Math.max(10, Math.round(opts.columns));
    const cw = opts.cell, lh = opts.lineH;
    const ssIn = E.parseSSInput(opts.ss);
    const ss = ssIn ? parseSS(ssIn.text, A.width, opts.ssRow != null ? A.rows[opts.ssRow].seq : A.rows[0].seq) : null;

    /* colours. Monochrome keeps the same drawing and takes the hue out,
       which is what a journal that charges for colour wants. */
    const C = opts.mono
      ? { strictBg: '#1c1c1c', strictFg: '#ffffff', simFg: '#1c1c1c', frame: '#8a8a8a',
          text: '#000000', dim: '#6a6a6a', ss: '#1c1c1c', bg: '#ffffff' }
      : { strictBg: '#e02020', strictFg: '#ffffff', simFg: '#e02020', frame: '#2b4fd8',
          text: '#111111', dim: '#777777', ss: '#111111', bg: '#ffffff' };

    /* the name gutter has to hold the longest row name AND the word
       "consensus", or the consensus label runs into the first column */
    const nameLens = A.rows.map(r => r.name.length).concat([4]);
    if (opts.consensus) nameLens.push(9);
    const nameW = Math.min(220, Math.max.apply(null, nameLens) * opts.nameFont * 0.62 + 10);
    const left = 14 + nameW;
    const numW = 52;                                   /* the residue number on the right */
    /* a block is only as wide as the columns actually in it: an alignment
       of 28 columns shown 60 to a line must not reserve 60 columns of
       empty paper and push its residue numbers off to the right */
    const drawnW = Math.min(per, A.width) * cw;
    let W = left + drawnW + numW + 14;

    /* how tall one block is: an optional structure track, one line per
       sequence, an optional consensus line, and air between blocks */
    const ssH = ss ? lh + 6 : 0;
    const rulerH = opts.ruler ? 13 : 0;
    const blockH = ssH + rulerH + A.rows.length * lh + (opts.consensus ? lh + 2 : 0) + 16;
    const nBlocks = Math.ceil(A.width / per) || 1;
    const titleH = opts.title ? 26 : 0;
    /* a key, so the figure explains its own symbols rather than needing a
       caption that says what # means */
    /** @type {Array<Array<{swatch?:string, text:string}>>} each line is a row of chips */
    const legendLines = [];
    if (opts.legend) {
      legendLines.push([
        { swatch: 'strict', text: 'identical in all ' + A.rows.length },
        { swatch: 'sim', text: 'conserved — ' + (E.MODE_LEGEND[A.mode] || A.mode) +
            ' — at or above ' + A.threshold }
      ]);
      if (opts.consensus) legendLines.push([
        { text: 'consensus:  UPPER identical · lower in more than half · ! I/V · $ L/M · % F/Y · # N/D/Q/E/B/Z · blank no majority' }
      ]);
      if (ss) legendLines.push([
        { text: (ss.dotBracket
            ? 'structure:  arrow = base-paired stem, pointing 5′→3′'
            : 'structure:  squiggle = α-helix, small squiggle = 3₁₀, large = π, arrow = β-strand, T = turn') +
            (ssIn ? '   (read as ' + ssIn.source + ')' : '') }
      ]);
    }
    const legendH = legendLines.length ? 10 + legendLines.length * 14 : 0;
    /* a key longer than the alignment widens the page rather than running
       off the edge of it */
    legendLines.forEach(line => {
      const w = 20 + line.reduce((t, it) => t + (it.swatch ? 15 : 0) + it.text.length * 4.7 + 18, 0);
      if (w > W) W = w;
    });
    const H = 16 + titleH + nBlocks * blockH + legendH + 10;

    /* residue counts along each row, so the number at the end of a block
       is the residue number and not the column number */
    const seen = A.rows.map(() => 0);

    let s = '';
    s += '<rect x="0" y="0" width="' + W + '" height="' + H + '" fill="' + C.bg + '"/>';
    if (opts.title) {
      s += '<text x="14" y="' + (20 + 0) + '" font-size="' + (opts.fontSize + 2) +
        '" font-weight="600" fill="' + C.text + '">' + esc(opts.title) + '</text>';
    }

    let y = 16 + titleH;
    for (let b = 0; b < nBlocks; b++) {
      const from = b * per, to = Math.min(A.width, from + per);
      let ty = y;

      /* ---- secondary structure above the block ---- */
      if (ss) {
        ss.runs.forEach(run => {
          const a = Math.max(run.start, from), z = Math.min(run.end, to);
          if (z <= a) return;
          const x0 = left + (a - from) * cw, x1 = left + (z - from) * cw;
          const mid = ty + lh * 0.55;
          if (run.kind === 'stem5' || run.kind === 'stem3') {
            /* a base-paired stem, drawn as a bar with the arrow pointing
               the way the strand runs into the pair */
            const bh = 3.2, head = Math.min(8, (x1 - x0) * 0.4);
            const d5 = run.kind === 'stem5';
            s += '<path d="M' + (d5 ? x0 : x0 + head).toFixed(1) + ' ' + (mid - bh).toFixed(1) +
              ' L' + (d5 ? x1 - head : x1).toFixed(1) + ' ' + (mid - bh).toFixed(1) +
              (d5 ? ' L' + x1.toFixed(1) + ' ' + mid.toFixed(1) + ' L' + (x1 - head).toFixed(1) + ' ' + (mid + bh).toFixed(1)
                  : ' L' + x1.toFixed(1) + ' ' + (mid + bh).toFixed(1)) +
              ' L' + (d5 ? x0 : x0 + head).toFixed(1) + ' ' + (mid + bh).toFixed(1) +
              (d5 ? '' : ' L' + x0.toFixed(1) + ' ' + mid.toFixed(1)) +
              ' Z" fill="' + C.ss + '"/>';
          } else if (run.kind === 'strand') {
            const head = Math.min(9, (x1 - x0) * 0.4);
            const bh = 4.6;
            s += '<path d="M' + x0.toFixed(1) + ' ' + (mid - bh).toFixed(1) +
              ' L' + (x1 - head).toFixed(1) + ' ' + (mid - bh).toFixed(1) +
              ' L' + (x1 - head).toFixed(1) + ' ' + (mid - bh - 3).toFixed(1) +
              ' L' + x1.toFixed(1) + ' ' + mid.toFixed(1) +
              ' L' + (x1 - head).toFixed(1) + ' ' + (mid + bh + 3).toFixed(1) +
              ' L' + (x1 - head).toFixed(1) + ' ' + (mid + bh).toFixed(1) +
              ' L' + x0.toFixed(1) + ' ' + (mid + bh).toFixed(1) + ' Z" fill="' + C.ss + '"/>';
          } else if (run.kind === 'turn' || run.kind === 'bend') {
            for (let k = a; k < z; k++) {
              s += '<text x="' + (left + (k - from) * cw + cw / 2).toFixed(1) + '" y="' + (mid + 4).toFixed(1) +
                '" font-size="' + (opts.fontSize - 1) + '" text-anchor="middle" fill="' + C.ss +
                '" font-weight="600">T</text>';
            }
          } else {
            /* a helix as a squiggle: amplitude says which helix it is */
            const amp = run.kind === 'helix310' ? 3.2 : (run.kind === 'helixPi' ? 7 : 5);
            const period = run.kind === 'helix310' ? cw * 1.1 : cw * 1.6;
            let d = 'M' + x0.toFixed(1) + ' ' + mid.toFixed(1);
            const steps = Math.max(2, Math.round((x1 - x0) / (period / 2)));
            for (let k = 1; k <= steps; k++) {
              const px = x0 + (x1 - x0) * (k / steps);
              const pmx = x0 + (x1 - x0) * ((k - 0.5) / steps);
              d += ' Q' + pmx.toFixed(1) + ' ' + (mid + (k % 2 ? -amp * 2 : amp * 2)).toFixed(1) +
                   ' ' + px.toFixed(1) + ' ' + mid.toFixed(1);
            }
            s += '<path d="' + d + '" fill="none" stroke="' + C.ss + '" stroke-width="1.6"/>';
          }
        });
        ty += ssH;
      }

      /* ---- the numbering ruler ---- */
      if (opts.ruler) {
        const refRow = A.rows[opts.numberRow != null ? opts.numberRow : 0].seq;
        let count = 0;
        for (let i = 0; i < from; i++) if (refRow[i] !== '-' && refRow[i] !== '.') count++;
        for (let i = from; i < to; i++) {
          if (refRow[i] === '-' || refRow[i] === '.') continue;
          count++;
          if (count % 10 === 0) {
            const x = left + (i - from) * cw + cw / 2;
            s += '<text x="' + x.toFixed(1) + '" y="' + (ty + 9) + '" font-size="8.5" text-anchor="middle" fill="' +
              C.dim + '">' + count + '</text>';
            s += '<line x1="' + x.toFixed(1) + '" y1="' + (ty + 11) + '" x2="' + x.toFixed(1) + '" y2="' + (ty + 14) +
              '" stroke="' + C.dim + '" stroke-width="0.7"/>';
          }
        }
        ty += rulerH;
      }

      /* ---- the blue frames round runs of conserved columns ---- */
      const boxTop = ty + 2;
      const boxBot = ty + A.rows.length * lh + 1;
      A.frames.forEach(fr => {
        const a = Math.max(fr.start, from), z = Math.min(fr.end, to);
        if (z <= a) return;
        s += '<rect x="' + (left + (a - from) * cw).toFixed(1) + '" y="' + boxTop.toFixed(1) +
          '" width="' + ((z - a) * cw).toFixed(1) + '" height="' + (boxBot - boxTop).toFixed(1) +
          '" fill="none" stroke="' + C.frame + '" stroke-width="0.9"/>';
      });

      /* ---- the rows ---- */
      A.rows.forEach((row, ri) => {
        const baseY = ty + ri * lh + lh - 3;
        s += '<text x="10" y="' + baseY.toFixed(1) + '" font-size="' + opts.nameFont +
          '" fill="' + C.text + '" font-family="inherit">' + esc(row.name) + '</text>';
        for (let i = from; i < to; i++) {
          const ch = row.seq[i] || '-';
          const x = left + (i - from) * cw;
          if (ch !== '-' && ch !== '.') seen[ri]++;
          if (ch === '-' || ch === '.') {
            s += '<text x="' + (x + cw / 2).toFixed(1) + '" y="' + baseY.toFixed(1) + '" font-size="' + opts.fontSize +
              '" text-anchor="middle" fill="' + C.dim + '">-</text>';
            continue;
          }
          let fill = C.text;
          if (A.strict[i]) {
            s += '<rect x="' + x.toFixed(1) + '" y="' + (baseY - lh + 4).toFixed(1) + '" width="' + cw.toFixed(1) +
              '" height="' + lh.toFixed(1) + '" fill="' + C.strictBg + '"/>';
            fill = C.strictFg;
          } else if (A.over[i]) {
            /* N-ESPript's chemical-property colouring replaces the flat
               red of a conserved column rather than sitting on top of it:
               purines red, pyrimidines blue, and only where the column is
               conserved enough to be worth colouring at all */
            fill = (opts.nucleotideColour && !opts.mono)
              ? (PURINE.indexOf(ch) !== -1 ? '#c02626' : (PYRIMIDINE.indexOf(ch) !== -1 ? '#2b4fd8' : C.simFg))
              : C.simFg;
          }
          s += '<text x="' + (x + cw / 2).toFixed(1) + '" y="' + baseY.toFixed(1) + '" font-size="' + opts.fontSize +
            '" text-anchor="middle" fill="' + fill + '" font-weight="' + (A.strict[i] ? '600' : '400') + '">' +
            esc(ch) + '</text>';
        }
        /* the residue number this row has reached, set just past the last
           column of THIS block rather than at the nominal block width */
        if (seen[ri] > 0) {
          s += '<text x="' + (left + (to - from) * cw + 8).toFixed(1) + '" y="' + baseY.toFixed(1) +
            '" font-size="8.5" fill="' + C.dim + '">' + seen[ri] + '</text>';
        }
      });
      ty += A.rows.length * lh;

      /* ---- the consensus line ---- */
      if (opts.consensus) {
        const cy = ty + lh;
        s += '<text x="10" y="' + cy.toFixed(1) + '" font-size="' + opts.nameFont + '" fill="' + C.dim + '">consensus</text>';
        for (let i = from; i < to; i++) {
          const ch = A.consensus[i] || ' ';
          if (ch === ' ') continue;
          s += '<text x="' + (left + (i - from) * cw + cw / 2).toFixed(1) + '" y="' + cy.toFixed(1) +
            '" font-size="' + opts.fontSize + '" text-anchor="middle" fill="' + C.dim + '">' + esc(ch) + '</text>';
        }
        ty += lh + 2;
      }
      y += blockH;
    }

    /* ---- the key ---- */
    if (legendLines.length) {
      let ly = y + 4;
      legendLines.forEach(line => {
        let lx = 10;
        line.forEach(item => {
          if (item.swatch === 'strict') {
            s += '<rect x="' + lx + '" y="' + (ly - 8).toFixed(1) + '" width="11" height="11" fill="' + C.strictBg + '"/>' +
              '<text x="' + (lx + 5.5) + '" y="' + (ly + 0.5).toFixed(1) + '" font-size="8.5" text-anchor="middle" fill="' +
              C.strictFg + '" font-weight="600">A</text>';
            lx += 15;
          } else if (item.swatch === 'sim') {
            s += '<rect x="' + lx + '" y="' + (ly - 8).toFixed(1) + '" width="11" height="11" fill="none" stroke="' +
              C.frame + '" stroke-width="0.9"/>' +
              '<text x="' + (lx + 5.5) + '" y="' + (ly + 0.5).toFixed(1) + '" font-size="8.5" text-anchor="middle" fill="' +
              C.simFg + '">A</text>';
            lx += 15;
          }
          s += '<text x="' + lx + '" y="' + ly.toFixed(1) + '" font-size="8.5" fill="' + C.dim + '">' +
            esc(item.text) + '</text>';
          lx += item.text.length * 4.7 + 18;
        });
        ly += 14;
      });
    }

    const head = '<svg xmlns="http://www.w3.org/2000/svg" version="1.1" width="' + Math.round(W) +
      '" height="' + Math.round(H) + '" viewBox="0 0 ' + Math.round(W) + ' ' + Math.round(H) +
      '" font-family="' + (opts.fontFamily || "'IBM Plex Mono',ui-monospace,Menlo,Consolas,monospace") + '">';
    return '<?xml version="1.0" encoding="UTF-8"?>\n' + head +
      '<title>' + esc(opts.title || 'Sequence alignment') + '</title>' + s + '</svg>';
  };

  if (typeof module !== 'undefined' && module.exports) module.exports = E;
  root.Bio = E;
})(typeof window !== 'undefined' ? window : globalThis);
