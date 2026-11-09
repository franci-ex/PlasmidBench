
/* ===========================================================================
   WHAT A READ ACTUALLY TELLS YOU ABOUT THE PROTEIN

   A Sanger read lands somewhere on a plasmid. The question that matters at
   the bench is never "how many bases differ" — it is "did my protein come
   out right, and which part of it did this read actually see?".

   Answering that needs the map from a genomic position to a codon, and
   that map is where naive implementations go wrong. A coding feature can
   be:

     - on the reverse strand, so the first codon is at its 3' end on the
       top strand and every base is complemented;
     - spliced, written join(a..b,c..d), so codon 40 may straddle an
       intron and the two halves are hundreds of bases apart;
     - wrapped across the origin of a circular plasmid, written
       join(5249..5467,1..438), where slicing start..end returns nothing
       at all;
     - all three at once, which is the CmR gene of pGro7 read backwards.

   So the map is built once, explicitly, as an ordered list of the genomic
   positions the coding sequence visits — index i in that list is base i of
   the CDS — and everything else is a lookup into it. A position is then
   codon Math.floor(i/3) + 1, at offset i % 3.

   Built on top of that:

     cdsIndex       the map, both ways
     readOverlap    which features a read's aligned span touches, and for
                    coding ones which residues it reached
     diffEffect     silent / missense / nonsense / frameshift for one
                    difference, with the codon and the residue number
     readProtein    the peptide the read itself supports over the stretch
                    it covers, against the reference peptide

   And, separately, the BLAST side: the page cannot reach NCBI, so the
   results come back by paste, and parseBlast reads the four shapes NCBI
   actually hands out.
   ========================================================================= */
(function (root) {
  'use strict';
  const E = root.Bio;

  const CODING = /^(CDS|gene|ORF|exon|mat_peptide|sig_peptide|transit_peptide|proprotein)$/i;
  E.isCodingFeature = f => !!f && CODING.test(String(f.type || ''));

  /* ------------------------------------------------------------------
     The map. `pos[i]` is the genomic index of base i of the coding
     sequence, in transcription order. For a reverse-strand feature that
     runs from the 3' end backwards, and the base there must be
     complemented to read the codon.
     ------------------------------------------------------------------ */
  const COMP = { A: 'T', C: 'G', G: 'C', T: 'A', U: 'A', N: 'N' };
  E.cdsIndex = function (rec, f) {
    if (!rec || !f) return null;
    const L = rec.seq.length;
    let pos = [];
    if (f.mixed && f.parts && f.parts.length) {
      /* a trans-spliced feature: each piece carries its own strand */
      f.parts.forEach(p => {
        const run = [];
        for (let i = p.start; i < p.end; i++) run.push(i);
        if (p.strand === -1) run.reverse();
        pos = pos.concat(run);
      });
    } else {
      E.featurePieces(f, L).forEach(p => { for (let i = p.start; i < p.end; i++) pos.push(i); });
      if (f.strand === -1) pos.reverse();
    }
    /* the reverse lookup, so a reference coordinate finds its codon in
       one step rather than by scanning the whole gene for every diff */
    const at = {};
    for (let i = 0; i < pos.length; i++) if (at[pos[i]] === undefined) at[pos[i]] = i;
    return {
      pos: pos, at: at, length: pos.length,
      strand: f.strand === -1 ? -1 : 1,
      /* base i of the coding sequence, read in transcription orientation */
      base: i => {
        const g = pos[i];
        if (g === undefined) return 'N';
        const c = rec.seq[g];
        return (f.mixed ? c : (f.strand === -1 ? (COMP[c] || 'N') : c));
      },
      codonOf: g => {
        const i = at[g];
        if (i === undefined) return null;
        return { cdsPos: i, codon: Math.floor(i / 3), within: i % 3, residue: Math.floor(i / 3) + 1 };
      }
    };
  };

  /* the codon at index `k`, as three characters in transcription order */
  function codonAt(idx, k) {
    return idx.base(k * 3) + idx.base(k * 3 + 1) + idx.base(k * 3 + 2);
  }
  E.cdsCodon = codonAt;

  /* ------------------------------------------------------------------
     One difference, turned into its consequence.

     `d` is what alignToReference produces: {type, refPos, ref, read} for
     a mismatch, {type, refPos, base} for an indel.
     ------------------------------------------------------------------ */
  E.diffEffect = function (rec, f, d, opts) {
    opts = opts || {};
    const idx = E.cdsIndex(rec, f);
    const hit = idx && idx.codonOf(d.refPos);
    if (!hit) return null;
    const out = {
      feature: f.label, residue: hit.residue, cdsPos: hit.cdsPos + 1,
      strand: idx.strand, kind: null, from: null, to: null, codon: null, mutCodon: null
    };
    if (d.type !== 'mismatch') {
      /* an indel inside a coding sequence shifts the frame unless it is a
         multiple of three, and a single-base call from a Sanger trace is
         never a multiple of three */
      out.kind = 'frameshift';
      out.text = 'frameshift from residue ' + hit.residue;
      return out;
    }
    const k = hit.codon;
    if ((k + 1) * 3 > idx.length) { out.kind = 'partial'; out.text = 'partial codon'; return out; }
    const codon = codonAt(idx, k);
    const readBase = idx.strand === -1 && !f.mixed ? (COMP[d.read] || 'N') : d.read;
    const arr = codon.split('');
    arr[hit.within] = readBase;
    const mut = arr.join('');
    const aa1 = E.translate(codon, { table: opts.table });
    const aa2 = E.translate(mut, { table: opts.table });
    out.codon = codon; out.mutCodon = mut; out.from = aa1; out.to = aa2;
    if (aa1 === aa2) { out.kind = 'silent'; out.text = 'silent · ' + aa1 + hit.residue; }
    else if (aa2 === '*') { out.kind = 'nonsense'; out.text = 'nonsense · ' + aa1 + hit.residue + '→stop'; }
    else if (aa1 === '*') { out.kind = 'readthrough'; out.text = 'stop lost · ' + hit.residue + '→' + aa2; }
    else { out.kind = 'missense'; out.text = 'missense · ' + aa1 + hit.residue + aa2; }
    return out;
  };

  /* ------------------------------------------------------------------
     What this one read covers.

     The aligned span is [refStart, refEnd) in reference coordinates. For
     every feature it touches, how much of it the read reached; for coding
     features, which residues, and whether the read agrees with the
     reference over them.
     ------------------------------------------------------------------ */
  E.readOverlap = function (rec, aln, opts) {
    opts = opts || {};
    if (!rec || !aln) return [];
    const L = rec.seq.length;
    const s = aln.refStart, e = aln.refEnd;
    const out = [];
    rec.features.forEach((f, fi) => {
      if (f.virtual && !opts.includeVirtual) return;
      const pieces = E.featurePieces(f, L);
      let covered = 0;
      pieces.forEach(p => { covered += Math.max(0, Math.min(p.end, e) - Math.max(p.start, s)); });
      if (!covered) return;
      const total = pieces.reduce((a, p) => a + (p.end - p.start), 0) || 1;
      const row = {
        index: fi, label: f.label, type: f.type, strand: f.strand,
        bases: covered, total: total, fraction: covered / total,
        coding: E.isCodingFeature(f), residues: null, aaFrom: null, aaTo: null, diffs: []
      };
      if (row.coding) {
        const idx = E.cdsIndex(rec, f);
        /* which coding positions the read actually reached — a spliced
           gene can be covered in two separate residue stretches, and
           reporting only the span between them would overstate it */
        let lo = Infinity, hi = -1, n = 0;
        for (let g = s; g < e && g < L; g++) {
          const i = idx.at[g];
          if (i === undefined) continue;
          n++; if (i < lo) lo = i; if (i > hi) hi = i;
        }
        if (n) {
          row.residues = { from: Math.floor(lo / 3) + 1, to: Math.floor(hi / 3) + 1, bases: n };
          row.aaTotal = Math.floor(idx.length / 3);
          row.aaCovered = Math.floor(n / 3);
        }
      }
      out.push(row);
    });
    /* the differences this read found, each placed in its feature */
    (aln.differences || []).forEach(d => {
      out.forEach(row => {
        if (!row.coding) return;
        const f = rec.features[row.index];
        const eff = E.diffEffect(rec, f, d, opts);
        if (eff) row.diffs.push(Object.assign({ refPos: d.refPos, type: d.type }, eff));
      });
    });
    out.sort((a, b) => b.bases - a.bases);
    return out;
  };

  /* The peptide this read supports over the stretch it covers, next to
     the reference peptide for the same stretch. Positions the read did
     not reach are not filled in from the reference — they come back as
     '.', because a residue nobody sequenced is not a result. */
  E.readProtein = function (rec, aln, f, opts) {
    opts = opts || {};
    const idx = E.cdsIndex(rec, f);
    if (!idx || !idx.length) return null;
    const s = aln.refStart, e = aln.refEnd;
    /* the read's own bases, in reference coordinates */
    const readAt = {};
    if (aln.alignment) {
      let g = s;
      const A = aln.alignment.a, Bq = aln.alignment.b;
      for (let i = 0; i < A.length; i++) {
        if (A[i] === '-') continue;              /* insertion in the read */
        if (Bq[i] !== '-') readAt[g] = Bq[i];
        g++;
      }
    }
    const nAA = Math.floor(idx.length / 3);
    let ref = '', obs = '';
    const changes = [];
    for (let k = 0; k < nAA; k++) {
      const gs = [idx.pos[k * 3], idx.pos[k * 3 + 1], idx.pos[k * 3 + 2]];
      const inRead = gs.every(g => g !== undefined && g >= s && g < e);
      const refCodon = codonAt(idx, k);
      const aaRef = E.translate(refCodon, { table: opts.table });
      ref += aaRef;
      if (!inRead) { obs += '.'; continue; }
      let mut = '';
      for (let j = 0; j < 3; j++) {
        const g = gs[j];
        let b = readAt[g] === undefined ? idx.base(k * 3 + j) : readAt[g];
        if (readAt[g] !== undefined && idx.strand === -1 && !f.mixed) b = COMP[b] || 'N';
        mut += b;
      }
      const aaObs = E.translate(mut, { table: opts.table });
      obs += aaObs;
      if (aaObs !== aaRef) changes.push({ residue: k + 1, from: aaRef, to: aaObs, codon: refCodon, mutCodon: mut });
    }
    const coveredAA = obs.split('').filter(c => c !== '.').length;
    return {
      name: f.label, reference: ref, observed: obs, changes: changes,
      length: nAA, covered: coveredAA,
      percent: nAA ? coveredAA / nAA * 100 : 0,
      identical: changes.length === 0 && coveredAA > 0
    };
  };

  /* ======================================================================
     BLAST RESULTS, PASTED BACK

     A published page cannot reach NCBI: the sandbox blocks outside
     requests, and NCBI's URL API needs a polling loop against a Request
     ID anyway. So the round trip is: the tool opens BLAST with the
     sequence, you run it, and you paste the output back. These are the
     four shapes NCBI's "Download" menu actually produces, plus the one
     people copy off the screen:

       Hit Table (text) / (csv)   one tabular row per HSP, 12 columns
       Descriptions (CSV)          the summary table above the alignments
       Single-file JSON (JSON2)    the full report
       Pairwise text               what the results page shows

     Returns a list of hits, best first, each with whatever the format
     carried. Nothing is invented: a field the format does not hold comes
     back null rather than guessed.
     ====================================================================== */
  function num(v) { const n = parseFloat(v); return isFinite(n) ? n : null; }

  E.parseBlast = function (text) {
    const t = String(text == null ? '' : text).replace(/\r\n?/g, '\n');
    if (!t.trim()) return { hits: [], format: null };
    /* ---- JSON2 ---- */
    if (/^\s*[{[]/.test(t)) {
      let data;
      /* only the parse is guarded: wrapping the reader too would relabel
         its own "no hits in here" as a syntax error */
      try { data = JSON.parse(t); }
      catch (e) { const err = new Error('That looks like JSON but does not parse: ' + e.message); err.code = 'badjson'; throw err; }
      return blastJson(data);
    }
    /* ---- Hit table: 12 tab- or comma-separated columns, outfmt 6/7/10 ---- */
    const tabRows = t.split('\n').filter(l => l.trim() && l[0] !== '#');
    const split = l => (l.indexOf('\t') !== -1 ? l.split('\t') : l.split(','));
    const looksHitTable = tabRows.length && tabRows.filter(l => {
      const c = split(l);
      return c.length >= 12 && num(c[2]) !== null && num(c[3]) !== null && num(c[10]) !== null;
    }).length >= Math.max(1, Math.ceil(tabRows.length * 0.6));
    if (looksHitTable) {
      const by = {};
      tabRows.forEach(l => {
        const c = split(l).map(x => x.trim());
        if (c.length < 12) return;
        const id = c[1];
        const hsp = {
          identity: num(c[2]), length: num(c[3]), mismatches: num(c[4]), gaps: num(c[5]),
          queryFrom: num(c[6]), queryTo: num(c[7]), hitFrom: num(c[8]), hitTo: num(c[9]),
          evalue: num(c[10]), bits: num(c[11])
        };
        if (!by[id]) by[id] = { id: id, title: null, accession: id, hsps: [], evalue: hsp.evalue, bits: hsp.bits, identity: hsp.identity, cover: null };
        by[id].hsps.push(hsp);
        if (hsp.bits != null && (by[id].bits == null || hsp.bits > by[id].bits)) {
          by[id].bits = hsp.bits; by[id].evalue = hsp.evalue; by[id].identity = hsp.identity;
        }
      });
      const hits = Object.keys(by).map(k => by[k]).sort((a, b) => (b.bits || 0) - (a.bits || 0));
      return { hits: hits, format: 'Hit table' };
    }
    /* ---- Descriptions CSV: a header row naming the columns ---- */
    const head = tabRows[0] || '';
    if (/description/i.test(head) && /(e[\s_-]?value|accession)/i.test(head)) {
      const cols = splitCsv(head).map(h => h.trim().toLowerCase());
      const col = (...names) => { for (const n of names) { const i = cols.findIndex(c => c.indexOf(n) !== -1); if (i !== -1) return i; } return -1; };
      const ci = { desc: col('description', 'scientific'), acc: col('accession'), ev: col('e value', 'evalue', 'e-value'),
                   bits: col('total score', 'max score', 'score'), ident: col('per. ident', 'ident'), cover: col('query cover') };
      const hits = tabRows.slice(1).map(l => {
        const c = splitCsv(l);
        if (c.length < 2) return null;
        const acc = ci.acc >= 0 ? (c[ci.acc] || '').trim() : null;
        return {
          id: acc || (c[ci.desc] || '').trim(), accession: acc,
          title: ci.desc >= 0 ? (c[ci.desc] || '').trim() : null,
          evalue: ci.ev >= 0 ? num(c[ci.ev]) : null,
          bits: ci.bits >= 0 ? num(c[ci.bits]) : null,
          identity: ci.ident >= 0 ? num(c[ci.ident]) : null,
          cover: ci.cover >= 0 ? num(String(c[ci.cover]).replace('%', '')) : null,
          hsps: []
        };
      }).filter(Boolean);
      if (hits.length) return { hits: hits, format: 'Descriptions table' };
    }
    /* ---- Pairwise text ---- */
    if (/Sequences producing significant alignments|^>\s*\S+.*\n\s*Length=/m.test(t) || /\bScore\s*=\s*[\d.]+\s*bits/.test(t)) {
      const hits = [];
      /* the per-subject blocks start at a line beginning with > */
      const blocks = t.split(/\n(?=>)/).slice(1);
      blocks.forEach(b => {
        const firstLine = (b.split('\n')[0] || '').replace(/^>\s*/, '').trim();
        const acc = (firstLine.match(/^(\S+)/) || [, null])[1];
        const sc = /Score\s*=\s*([\d.]+)\s*bits[^\n]*Expect(?:\(\d+\))?\s*=\s*([\deE.+-]+)/.exec(b);
        const id = /Identities\s*=\s*(\d+)\/(\d+)\s*\((\d+)%\)/.exec(b);
        const len = /Length\s*=\s*(\d+)/.exec(b);
        hits.push({
          id: acc, accession: acc,
          title: firstLine.replace(/^\S+\s*/, '') || firstLine,
          bits: sc ? num(sc[1]) : null,
          evalue: sc ? num(sc[2]) : null,
          identity: id ? num(id[3]) : null,
          alignedLength: id ? num(id[2]) : null,
          subjectLength: len ? num(len[1]) : null,
          cover: null, hsps: []
        });
      });
      if (hits.length) return { hits: hits, format: 'Pairwise text' };
    }
    const err = new Error('That is not a BLAST result PlasmidBench recognises. From the BLAST page, ' +
      'Download → Hit Table (text), Descriptions (CSV) or Single-file JSON all work, and so does the ' +
      'pairwise text copied straight off the screen.');
    err.code = 'unknown';
    throw err;
  };

  function blastJson(data) {
    /* NCBI's JSON2 nests the report; a few tools flatten it */
    const reports = [];
    const walk = o => {
      if (!o || typeof o !== 'object') return;
      if (o.report && o.report.results) reports.push(o.report);
      if (Array.isArray(o.BlastOutput2)) o.BlastOutput2.forEach(walk);
      else if (o.BlastOutput2) walk(o.BlastOutput2);
      if (Array.isArray(o)) o.forEach(walk);
    };
    walk(data);
    if (!reports.length && data && data.results) reports.push(data);
    const hits = [];
    reports.forEach(rep => {
      const search = rep.results && (rep.results.search || (rep.results.bl2seq && rep.results.bl2seq[0]));
      if (!search || !Array.isArray(search.hits)) return;
      search.hits.forEach(h => {
        const d = (h.description && h.description[0]) || {};
        const best = (h.hsps || []).reduce((a, x) => (!a || (x.bit_score || 0) > (a.bit_score || 0) ? x : a), null);
        hits.push({
          id: d.accession || d.id || String(h.num),
          accession: d.accession || null,
          title: d.title || d.sciname || null,
          bits: best ? num(best.bit_score) : null,
          evalue: best ? num(best.evalue) : null,
          identity: best && best.align_len ? (best.identity / best.align_len * 100) : null,
          alignedLength: best ? num(best.align_len) : null,
          subjectLength: num(h.len),
          queryFrom: best ? num(best.query_from) : null,
          queryTo: best ? num(best.query_to) : null,
          hitSeq: best ? (best.hseq || null) : null,
          querySeq: best ? (best.qseq || null) : null,
          cover: null,
          hsps: (h.hsps || []).map(x => ({
            bits: num(x.bit_score), evalue: num(x.evalue),
            queryFrom: num(x.query_from), queryTo: num(x.query_to),
            hitFrom: num(x.hit_from), hitTo: num(x.hit_to),
            identity: x.align_len ? (x.identity / x.align_len * 100) : null, length: num(x.align_len)
          }))
        });
      });
    });
    if (!hits.length) {
      const err = new Error('That JSON holds no BLAST hits. NCBI’s Download → Single-file JSON is the one that works.');
      err.code = 'nohits';
      throw err;
    }
    hits.sort((a, b) => (b.bits || 0) - (a.bits || 0));
    return { hits: hits, format: 'Single-file JSON' };
  }

  /* a quote-aware CSV splitter — a BLAST description routinely contains
     commas, and splitting on them naively shifts every later column */
  function splitCsv(line) {
    const out = []; let cur = '', q = false;
    for (let i = 0; i < line.length; i++) {
      const c = line[i];
      if (q) {
        if (c === '"' && line[i + 1] === '"') { cur += '"'; i++; }
        else if (c === '"') q = false;
        else cur += c;
      } else if (c === '"') q = true;
      else if (c === ',') { out.push(cur); cur = ''; }
      else cur += c;
    }
    out.push(cur);
    return out;
  }
  E.splitCsvLine = splitCsv;

  if (typeof module !== 'undefined' && module.exports) module.exports = E;
  root.Bio = E;
})(typeof window !== 'undefined' ? window : globalThis);

