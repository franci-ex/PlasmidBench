/* ============================================================================
   PlasmidBench — core molecular biology engine
   Pure, dependency-free functions. Written to be unit-testable in Node and
   inlined verbatim into the browser app. Accuracy is the priority:
   algorithms follow published references (SantaLucia 1998 NN thermodynamics,
   NCBI genetic code tables, REBASE recognition sites).
   ============================================================================ */
/* The types these functions are annotated with — Feature, SeqRecord, Enzyme,
   PcrResult, FeatureList, SeqTextReport, CodedError — are declared in
   types.d.ts, which holds no code. See the header of that file for the two
   conventions worth knowing: positions are 0-based and end-exclusive
   everywhere, and a few functions return an array that also carries named
   properties. Checked with: npx tsc -p tsconfig.json */

(function (root) {
  'use strict';
  const E = {};

  /* ---------- Alphabet / IUPAC ---------- */
  const COMP = {
    A:'T',T:'A',G:'C',C:'G',U:'A',
    R:'Y',Y:'R',S:'S',W:'W',K:'M',M:'K',
    B:'V',V:'B',D:'H',H:'D',N:'N','-':'-','.':'.'
  };
  const IUPAC = {
    A:['A'],C:['C'],G:['G'],T:['T'],U:['T'],
    R:['A','G'],Y:['C','T'],S:['G','C'],W:['A','T'],K:['G','T'],M:['A','C'],
    B:['C','G','T'],D:['A','G','T'],H:['A','C','T'],V:['A','C','G'],N:['A','C','G','T']
  };

  E.clean = function (s) {
    return (s || '').toUpperCase().replace(/[^A-Z\-\.]/g, '');
  };
  E.cleanSeq = function (s) {
    // keep only valid nucleotide/IUPAC letters
    return (s || '').toUpperCase().replace(/[^ACGTURYSWKMBDHVN]/g, '');
  };
  /* What cleanSeq would silently throw away from a block of sequence text.
     Whitespace, digits and alignment gaps are layout; anything else dropped
     is a real character (X, *, E, F, I, L, P, Q, ...) and usually means the
     text is protein or corrupted.                                          */
  /** @param {string} raw @returns {SeqTextReport} */
  E.seqTextReport = function (raw) {
    const letters = String(raw || '').replace(/[\s\d\/\-.]/g, '');
    const up = letters.toUpperCase();
    const bad = up.replace(/[ACGTURYSWKMBDHVN]/g, '');
    const acgtn = up.replace(/[^ACGTUN]/g, '').length;
    const frac = up.length ? acgtn / up.length : 1;
    const proteinOnly = /[EFILPQ*JOZX]/.test(up);
    return {
      seq: E.cleanSeq(raw).replace(/U/g, 'T'),
      /* The same letters with nothing thrown away. E, F, I, L, P and Q are
         not nucleotide codes, so DNA cleaning silently deletes two residues
         in five from a protein record and leaves a string that still reads
         as DNA. A reader that knows it has protein takes this instead, so
         the residues survive intact and are never counted as damage. */
      protein: up.replace(/[^A-Z*]/g, ''),
      dropped: bad.length,
      droppedChars: Array.from(new Set(bad.split(''))).slice(0, 12).join(''),
      hadU: /U/.test(up),
      looksProtein: up.length >= 8 && ((proteinOnly && frac < 0.9) || frac < 0.35)
    };
  };
  E.complement = function (s) {
    let out = '';
    for (const ch of s.toUpperCase()) out += (COMP[ch] || 'N');
    return out;
  };
  E.reverse = function (s) { return s.split('').reverse().join(''); };
  E.revcomp = function (s) {
    let out = '';
    for (let i = s.length - 1; i >= 0; i--) out += (COMP[s[i].toUpperCase()] || 'N');
    return out;
  };
  // does base `b` (concrete ACGT) match IUPAC code `code`?
  E.iupacMatch = function (code, b) {
    const set = IUPAC[code];
    return set ? set.indexOf(b) !== -1 : false;
  };

  /* ---------- Composition ---------- */
  E.gcContent = function (s) {
    s = E.cleanSeq(s);
    if (!s.length) return 0;
    let gc = 0;
    for (const ch of s) if (ch === 'G' || ch === 'C') gc++;
    return gc / s.length * 100;
  };
  E.baseCounts = function (s) {
    s = E.cleanSeq(s);
    const c = { A:0,C:0,G:0,T:0,N:0 };
    for (const ch of s) { if (ch in c) c[ch]++; else if (ch==='U') c.T++; else c.N++; }
    return c;
  };
  // ssDNA / dsDNA molecular weight (g/mol). Anhydrous MW of dNMP residues
  // minus one water-equiv per phosphodiester; add 79.0 (HPO3) + 18.02 for 5'-phosphate?
  // We report the standard ssDNA formula used by common calculators.
  const DNA_MONO = { A:313.21, C:289.18, G:329.21, T:304.2 }; // monophosphate residues
  E.ssDNAWeight = function (s) {
    s = E.cleanSeq(s);
    if (!s.length) return 0;
    let w = 0;
    for (const ch of s) w += (DNA_MONO[ch] || DNA_MONO[ch === 'U' ? 'T' : 'A'] || 0);
    // -61.96 removes HPO2 from the 5' end, giving a 5'-OH oligo (IDT formula)
    return w - 61.96;
  };
  E.dsDNAWeight = function (s) {
    return E.ssDNAWeight(s) + E.ssDNAWeight(E.revcomp(s));
  };

  /* ---------- Melting temperature ----------
     Two methods:
     1. Wallace rule (2(A+T)+4(G+C)) for short oligos < 14 nt.
     2. SantaLucia (1998) unified nearest-neighbor with salt correction.
        Returns salt-adjusted Tm in Celsius.
     Parameters follow SantaLucia & Hicks 2004 defaults.
  ------------------------------------------------------------------ */
  // Unified NN parameters (SantaLucia 1998): dH (kcal/mol), dS (cal/K/mol)
  const NN_H = {
    AA:-7.9, AT:-7.2, TA:-7.2, CA:-8.5, GT:-8.4, CT:-7.8, GA:-8.2,
    CG:-10.6, GC:-9.8, GG:-8.0,
    TT:-7.9, AC:-8.4, TG:-8.5, AG:-7.8, TC:-8.2, CC:-8.0
  };
  const NN_S = {
    AA:-22.2, AT:-20.4, TA:-21.3, CA:-22.7, GT:-22.4, CT:-21.0, GA:-22.2,
    CG:-27.2, GC:-24.4, GG:-19.9,
    TT:-22.2, AC:-22.4, TG:-22.7, AG:-21.0, TC:-22.2, CC:-19.9
  };
  // initiation terms
  const INIT_H_GC = 0.1, INIT_S_GC = -2.8;   // with terminal G·C
  const INIT_H_AT = 2.3, INIT_S_AT = 4.1;    // with terminal A·T
  const R = 1.987; // cal/(K·mol)

  E.tmWallace = function (s) {
    s = E.cleanSeq(s);
    let at = 0, gc = 0;
    for (const ch of s) { if (ch==='G'||ch==='C') gc++; else at++; }
    return 2 * at + 4 * gc;
  };

  // Nearest-neighbour thermodynamics for a duplex oligo.
  // Returns { dH (kcal/mol), dS (cal/mol/K), dG37 (kcal/mol), tm1M (C) }
  E.nnThermo = function (s, opts) {
    opts = opts || {};
    s = E.cleanSeq(s).replace(/U/g, 'T');
    if (s.length < 2) return null;
    // Concentration term follows Bio.SeqUtils.MeltingTemp.Tm_NN exactly:
    // dnac1 is the higher-concentration strand (the primer), dnac2 the lower
    // (the template, ignored in PCR). k = dnac1 - dnac2/2.
    const dnac1 = (opts.oligo != null ? opts.oligo : 0.25e-6);
    const dnac2 = (opts.oligo2 != null ? opts.oligo2 : 0);
    let dH = 0, dS = 0;
    for (let i = 0; i < s.length - 1; i++) {
      const nn = s[i] + s[i + 1];
      if (!(nn in NN_H)) return null; // ambiguous base: NN model not defined
      dH += NN_H[nn];
      dS += NN_S[nn];
    }
    const first = s[0], last = s[s.length - 1];
    dH += (first === 'G' || first === 'C') ? INIT_H_GC : INIT_H_AT;
    dS += (first === 'G' || first === 'C') ? INIT_S_GC : INIT_S_AT;
    dH += (last === 'G' || last === 'C') ? INIT_H_GC : INIT_H_AT;
    dS += (last === 'G' || last === 'C') ? INIT_S_GC : INIT_S_AT;
    const selfComp = (s === E.revcomp(s));
    let k;
    if (selfComp) { dS += -1.4; k = dnac1; }   // symmetry correction
    else k = dnac1 - dnac2 / 2;
    const tm1M = (dH * 1000) / (dS + R * Math.log(k)) - 273.15;
    const dG37 = dH - (310.15 * dS / 1000);
    return { dH, dS, dG37, tm1M, selfComp, length: s.length };
  };

  /* Salt-corrected Tm.
     opts: { oligo (M, default 0.25e-6), Na (M monovalent, default 0.05),
             Mg (mM, default 0), dNTP (mM, default 0) }
     Monovalent-only: Owczarzy et al. 2004, Biochemistry 43:3537, eq. 22.
     With Mg2+: Owczarzy et al. 2008, Biochemistry 47:5336, eq. 16 (+ the
     monovalent-dependent a/d/g terms of eq. 19 in the mixed-ion regime).
     This is the same treatment IDT's OligoAnalyzer uses.                     */
  E.tmNN = function (s, opts) {
    opts = opts || {};
    const th = E.nnThermo(s, opts);
    if (!th) return NaN;
    s = E.cleanSeq(s).replace(/U/g, 'T');
    const Na = (opts.Na != null ? opts.Na : 0.05);          // M
    const MgTot = (opts.Mg != null ? opts.Mg : 0) / 1000;   // mM -> M
    const dNTP = (opts.dNTP != null ? opts.dNTP : 0) / 1000;
    // Free Mg2+ after dNTP chelation, via the Mg:dNTP dissociation quadratic
    // used by Bio.SeqUtils.MeltingTemp (ka = 3e4).
    let Mg = MgTot;
    if (dNTP > 0 && MgTot > 0) {
      const ka = 3e4;
      Mg = (-(ka * dNTP - ka * MgTot + 1.0) +
            Math.sqrt(Math.pow(ka * dNTP - ka * MgTot + 1.0, 2) + 4.0 * ka * MgTot)) / (2.0 * ka);
    }
    if (Mg < 0) Mg = 0;
    const fgc = E.gcContent(s) / 100;
    const Nbp = s.length;
    const invTm1M = 1 / (th.tm1M + 273.15);
    let inv;

    const ratio = Na > 0 ? Math.sqrt(Mg) / Na : (Mg > 0 ? Infinity : 0);

    if (Mg <= 0 || ratio < 0.22) {
      // monovalent dominates
      const lnNa = Math.log(Na <= 0 ? 1e-9 : Na);
      inv = invTm1M + (4.29 * fgc - 3.95) * 1e-5 * lnNa + 9.40e-6 * lnNa * lnNa;
    } else {
      const lnMg = Math.log(Mg);
      let a = 3.92e-5, b = -9.11e-6, c = 6.26e-5,
          d = 1.42e-5, e = -4.82e-4, f = 5.25e-4, g = 8.31e-5;
      if (ratio < 6.0) {
        // mixed monovalent / divalent: a, d, g become monovalent-dependent
        const lnNa = Math.log(Na);
        a = 3.92e-5 * (0.843 - 0.352 * Math.sqrt(Na) * lnNa);
        d = 1.42e-5 * (1.279 - 4.03e-3 * lnNa - 8.03e-3 * lnNa * lnNa);
        g = 8.31e-5 * (0.486 - 0.258 * lnNa + 5.25e-3 * lnNa * lnNa * lnNa);
      }
      // parameters are signed, so every term adds directly
      inv = invTm1M
          + a
          + b * lnMg
          + fgc * (c + d * lnMg)
          + (1 / (2 * (Nbp - 1))) * (e + f * lnMg + g * lnMg * lnMg);
    }
    return 1 / inv - 273.15;
  };

  // Typical PCR conditions preset (what a thermocycler reaction actually sees)
  E.PCR_CONDITIONS = { oligo: 0.25e-6, Na: 0.05, Mg: 1.5, dNTP: 0.2 };

  /* ---------- Genetic codes (NCBI translation tables) ----------
     Codon order in the standard NCBI string is TTT,TTC,...; we build a map
     from the AA string + start string. Base order: T,C,A,G for all three pos.
  ------------------------------------------------------------------ */
  const BASES = ['T', 'C', 'A', 'G'];
  function buildCodons() {
    const list = [];
    for (const b1 of BASES) for (const b2 of BASES) for (const b3 of BASES) list.push(b1 + b2 + b3);
    return list;
  }
  const CODON_ORDER = buildCodons();
  // NCBI aa / starts strings
  const CODE_TABLES = {
    1:  { name:'Standard',
          aa:'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
          starts:'---M------**--*----M---------------M----------------------------' },
    2:  { name:'Vertebrate Mitochondrial',
          aa:'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG',
          starts:'----------**--------------------MMMM----------**---M------------' },
    3:  { name:'Yeast Mitochondrial',
          aa:'FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
          starts:'----------**----------------------MM---------------M------------' },
    4:  { name:'Mold/Protozoan Mitochondrial',
          aa:'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
          starts:'--MM------**-------M---------------M------------MMMM---------M---' },
    5:  { name:'Invertebrate Mitochondrial',
          aa:'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG',
          starts:'---M------**--------------------MMMM---------------M------------' },
    11: { name:'Bacterial/Archaeal/Plastid',
          aa:'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
          starts:'---M------**--*----M------------MMMM---------------M------------' },
    6:  { name:'Ciliate/Dasycladacean',
          aa:'FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
          starts:'---M------**--------------------------------------------------' + '--' },
    13: { name:'Ascidian Mitochondrial',
          aa:'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG',
          starts:'---M------**----------------------MM---------------M------------' }
  };
  function buildAAmap(tableId) {
    const t = CODE_TABLES[tableId] || CODE_TABLES[1];
    const aaMap = {}, startMap = {};
    for (let i = 0; i < CODON_ORDER.length; i++) {
      aaMap[CODON_ORDER[i]] = t.aa[i];
      startMap[CODON_ORDER[i]] = (t.starts[i] === 'M');
    }
    return { aaMap, startMap, name: t.name };
  }
  E.geneticCodes = function () {
    return Object.keys(CODE_TABLES).map(k => ({ id: +k, name: CODE_TABLES[k].name }));
  };

  // translate; opts {table, toStop, frame(0..2)}
  E.translate = function (s, opts) {
    opts = opts || {};
    const { aaMap } = buildAAmap(opts.table || 1);
    s = E.cleanSeq(s).replace(/U/g, 'T');
    const frame = opts.frame || 0;
    let prot = '';
    for (let i = frame; i + 3 <= s.length; i += 3) {
      const codon = s.substr(i, 3);
      let aa = aaMap[codon];
      if (aa === undefined) {
        // ambiguous codon: resolve if unambiguous in product
        aa = resolveAmbiguous(codon, aaMap);
      }
      if (aa === '*' && opts.toStop) break;
      prot += aa;
    }
    return prot;
  };
  function resolveAmbiguous(codon, aaMap) {
    const opts0 = IUPAC[codon[0]] || ['N'];
    const opts1 = IUPAC[codon[1]] || ['N'];
    const opts2 = IUPAC[codon[2]] || ['N'];
    if (opts0[0] === 'N' || opts1[0] === 'N' || opts2[0] === 'N') return 'X';
    let aa = null;
    for (const a of opts0) for (const b of opts1) for (const c of opts2) {
      const x = aaMap[a + b + c];
      if (x === undefined) return 'X';
      if (aa === null) aa = x; else if (aa !== x) return 'X';
    }
    return aa || 'X';
  }

  // three-letter + full names
  E.AA3 = {A:'Ala',R:'Arg',N:'Asn',D:'Asp',C:'Cys',E:'Glu',Q:'Gln',G:'Gly',
    H:'His',I:'Ile',L:'Leu',K:'Lys',M:'Met',F:'Phe',P:'Pro',S:'Ser',T:'Thr',
    W:'Trp',Y:'Tyr',V:'Val','*':'***','X':'Xaa'};

  /* ---------- ORF finder ---------- */
  // Find ORFs in all 6 frames. opts {table, minAA, circular, requireStart}
  E.findORFs = function (seq, opts) {
    opts = opts || {};
    const table = opts.table || 1;
    const minAA = opts.minAA != null ? opts.minAA : 30;
    const requireStart = opts.requireStart !== false;
    const { aaMap, startMap } = buildAAmap(table);
    seq = E.cleanSeq(seq).replace(/U/g, 'T');
    const L = seq.length;
    const orfs = [];
    /* On a circle a reading frame runs past base 1, so the search is done on
       the sequence written twice and anything starting in the second copy, or
       longer than one turn, is dropped. Coordinates come back in the feature
       convention: end <= start means the ORF crosses the origin. */
    const circular = !!opts.circular;
    const scanSeq = circular ? seq + seq : seq;
    const strands = [{ s: scanSeq, dir: 1 }, { s: E.revcomp(scanSeq), dir: -1 }];
    for (const { s, dir } of strands) {
      for (let frame = 0; frame < 3; frame++) {
        let i = frame;
        let startIdx = -1;
        while (i + 3 <= s.length) {
          const codon = s.substr(i, 3);
          const aa = aaMap[codon];
          const isStart = startMap[codon];
          if (startIdx === -1) {
            if (requireStart ? isStart : aa !== '*') startIdx = i;
          }
          if (startIdx !== -1 && aa === '*') {
            const nt = i + 3 - startIdx;
            const aaLen = nt / 3 - 1;
            if (aaLen >= minAA && (!circular || nt <= L)) {
              const o = makeOrf(startIdx, i + 3, dir, circular ? L * 2 : L, aaLen, s, startIdx, i + 3, table);
              if (!circular) orfs.push(o);
              else if (o.start < L) {                       // starts in the first copy
                o.end = o.end % L === 0 && o.end !== o.start ? L : o.end % L;
                if (o.end === o.start && nt < L) o.end = o.start;   // full turn only when it is one
                o.wrapped = o.end <= o.start;
                o.frame = dir === 1 ? (o.start % 3) + 1 : -((o.start % 3) + 1);
                orfs.push(o);
              }
            }
            startIdx = -1;
          }
          i += 3;
        }
      }
    }
    orfs.sort((a, b) => a.start - b.start || b.length - a.length);
    if (!circular) return orfs;
    const seen = new Set(), out = [];
    for (const o of orfs) {
      const k = o.start + ':' + o.end + ':' + o.strand;
      if (seen.has(k)) continue;
      seen.add(k); out.push(o);
    }
    return out;
  };
  function makeOrf(sStart, sEnd, dir, L, aaLen, s, a, b, table) {
    let start, end;
    if (dir === 1) { start = sStart; end = sEnd; }
    else { start = L - sEnd; end = L - sStart; }
    return {
      start, end, strand: dir, length: aaLen,
      ntLength: sEnd - sStart,
      frame: (dir === 1 ? (sStart % 3) + 1 : -((sStart % 3) + 1)),
      protein: E.translate(s.substring(a, b), { table, toStop: true })
    };
  }

  /* ---------- Restriction enzyme catalogue ----------
     The full set of commercially available enzymes from REBASE, as bundled
     with Biopython (Bio.Restriction). Cut coordinates were converted from
     Biopython's (fst5, fst3) characterisation and then cross-validated against
     Bio.Restriction.search() on hundreds of thousands of enzyme x sequence x
     topology comparisons (see test-rebase.js). Methylation-dependent enzymes
     (which cut only methylated DNA and cannot be modelled on plain sequence)
     are excluded. Each record: {name, site (IUPAC, top strand 5'->3'),
     top/bot cut positions from the site start for ordinary enzymes, or o1/o2
     offsets beyond the site end for Type IIS enzymes, typeIIS, suppliers}.
  ------------------------------------------------------------------ */
var REBASE_PACKED = "Eco72I|CACGTG|3|3|0|B\nEco53kI|GAGCTC|3|3|0|N\nMseI|TTAA|1|3|0|N\nAcc65I|GGTACC|1|5|0|BINV\nBclI|TGATCA|1|5|0|BCJMNORS\nBlpI|GCTNAGC|2|5|0|N\nBbvI|GCAGC|8|12|1|N\nMaeII|ACGT|1|3|0|M\nFspEI|CC|12|16|1|N\nBmrI|ACTGGG|5|4|1|N\nVspI|ATTAAT|2|4|0|BIV\nSspDI|GGCGCC|1|5|0|B\nEcoO65I|GGTNACC|1|6|0|K\nBslI|CCNNNNNNNGG|7|4|0|N\nHpyCH4III|ACNGT|3|2|0|N\nBbrPI|CACGTG|3|3|0|M\nEco32I|GATATC|3|3|0|B\nEco147I|AGGCCT|3|3|0|B\nMluI|ACGCGT|1|5|0|BIJKMNOQRSV\nBspT104I|TTCGAA|2|4|0|K\nBanI|GGYRCC|1|5|0|N\nKzo9I|GATC|0|4|0|I\nTaqI|TCGA|1|3|0|BCIJKMNQRVX\nFaqI|GGGAC|10|14|1|B\nSchI|GAGTC|5|5|1|B\nZraI|GACGTC|3|3|0|INV\nAvaII|GGWCC|1|4|0|JNX\nHapII|CCGG|1|3|0|BK\nSgsI|GGCGCGCC|2|6|0|B\nTsp45I|GTSAC|0|5|0|N\nEco91I|GGTNACC|1|6|0|B\nSbfI|CCTGCAGG|6|2|0|INV\nBglI|GCCNNNNNGGC|7|4|0|BCIJKNOQRVX\nEco24I|GRGCYC|5|1|0|B\nAor51HI|AGCGCT|3|3|0|K\nDinI|GGCGCC|3|3|0|V\nMboI|GATC|0|4|0|BCKNQRXY\nBse118I|RCCGGY|1|5|0|IV\nBspEI|TCCGGA|1|5|0|N\nKroI|GCCGGC|1|5|0|I\nXagI|CCTNNNNNAGG|5|6|0|B\nBtrI|CACGTC|3|3|0|I\nCseI|GACGC|5|10|1|B\nSseBI|AGGCCT|3|3|0|C\nSpeI|ACTAGT|1|5|0|BJKMNOQRSX\nAspS9I|GGNCC|1|4|0|IV\nEcoRI|GAATTC|1|5|0|BCIJKMNOQRSVXY\nBstACI|GRCGYC|2|4|0|I\nSfuI|TTCGAA|2|4|0|MS\nTfiI|GAWTC|1|4|0|N\nEco81I|CCTNAGG|2|5|0|BK\nSacI|GAGCTC|5|1|0|BJKMNOQRSX\nBaeGI|GKGCMC|5|1|0|N\nDseDI|GACNNNNNNGTC|7|5|0|IV\nBsp1407I|TGTACA|1|5|0|BK\nKasI|GGCGCC|1|5|0|N\nBsuRI|GGCC|2|2|0|BI\nAfaI|GTAC|2|2|0|BK\nBpu14I|TTCGAA|2|4|0|IV\nBstSFI|CTRYAG|1|5|0|I\nTseFI|GTSAC|0|5|0|I\nSfr303I|CCGCGG|4|2|0|IV\nAlwNI|CAGNNNCTG|6|3|0|N\nHsp92I|GRCGYC|2|4|0|R\nSphI|GCATGC|5|1|0|BCIJKMNOQRVX\nSau3AI|GATC|0|4|0|CJKMNX\nAjiI|CACGTC|3|3|0|B\nBstV2I|GAAGAC|2|6|1|IV\nBstYI|RGATCY|1|5|0|N\nSmiMI|CAYNNNNRTG|5|5|0|IV\nAjnI|CCWGG|0|5|0|I\nSaqAI|TTAA|1|3|0|B\nStyI|CCWWGG|1|5|0|CJN\nPluTI|GGCGCC|5|1|0|N\nCspI|CGGWCCG|2|5|0|R\nEaeI|YGGCCR|1|5|0|KN\nBssMI|GATC|0|4|0|V\nCaiI|CAGNNNCTG|6|3|0|B\nBauI|CACGAG|1|5|0|B\nAcc16I|TGCGCA|3|3|0|IV\nHpyCH4IV|ACGT|1|3|0|N\nAsuII|TTCGAA|2|4|0|C\nBsp120I|GGGCCC|1|5|0|B\nBstPI|GGTNACC|1|6|0|K\nSgeI|CNNG|9|13|1|B\nSdaI|CCTGCAGG|6|2|0|B\nAfiI|CCNNNNNNNGG|7|4|0|V\nBstMAI|GTCTC|1|5|1|IV\nRseI|CAYNNNNRTG|5|5|0|B\nPspOMI|GGGCCC|1|5|0|INV\nBspHI|TCATGA|1|5|0|N\nHin6I|GCGC|1|3|0|B\nAccI|GTMKAC|2|4|0|BJKMNQX\nSgfI|GCGATCGC|5|3|0|R\nCviQI|GTAC|1|3|0|N\nBssAI|RCCGGY|1|5|0|C\nPteI|GCGCGC|1|5|0|B\nSinI|GGWCC|1|4|0|X\nCpoI|CGGWCCG|2|5|0|BK\nNsiI|ATGCAT|5|1|0|JMNQSX\nXceI|RCATGY|5|1|0|B\nBstMWI|GCNNNNNNNGC|7|4|0|I\nHindIII|AAGCTT|1|5|0|BCIJKMNOQRSVXY\nBssSI|CACGAG|1|5|0|N\nAsp718I|GGTACC|1|5|0|MS\nBsp19I|CCATGG|1|5|0|IV\nBsp13I|TCCGGA|1|5|0|IV\nBstDSI|CCRYGG|1|5|0|IV\nSatI|GCNGC|2|3|0|B\nRgaI|GCGATCGC|5|3|0|I\nAccB7I|CCANNNNNTGG|7|4|0|IV\nNotI|GCGGCCGC|2|6|0|BCJKMNOQRSX\nCfr9I|CCCGGG|1|5|0|B\nBsaI|GGTCTC|1|5|1|N\nBsmFI|GGGAC|10|14|1|N\nBspTNI|GGTCTC|1|5|1|X\nBsp68I|TCGCGA|3|3|0|B\nPceI|AGGCCT|3|3|0|IV\nPspN4I|GGNNCC|3|3|0|I\nFspBI|CTAG|1|3|0|B\nXspI|CTAG|1|3|0|K\nGsuI|CTGGAG|16|14|1|B\nClaI|ATCGAT|2|4|0|BKMNQRSX\nPspXI|VCTCGAGB|2|6|0|IN\nScrFI|CCNGG|2|3|0|JN\nBst2UI|CCWGG|2|3|0|IV\nKpnI|GGTACC|5|1|0|BCIJKMNOQRSVXY\nSstI|GAGCTC|5|1|0|C\nBst4CI|ACNGT|3|2|0|IV\nBbvCI|CCTCAGC|2|5|0|N\nBst2BI|CACGAG|1|5|0|I\nBsrI|ACTGG|1|-1|1|N\nSwaI|ATTTAAAT|4|4|0|JMNS\nBshVI|ATCGAT|2|4|0|V\nBshNI|GGYRCC|1|5|0|B\nBssT1I|CCWWGG|1|5|0|IV\nPsyI|GACNNNGTC|4|5|0|B\nPsp124BI|GAGCTC|5|1|0|IV\nNsbI|TGCGCA|3|3|0|B\nPpu21I|YACGTR|3|3|0|B\nBtuMI|TCGCGA|3|3|0|V\nBsrBI|CCGCTC|3|3|0|N\nBfuAI|ACCTGC|4|8|1|N\nBsmAI|GTCTC|1|5|1|N\nBspPI|GGATC|4|5|1|B\nBshFI|GGCC|2|2|0|C\nFauNDI|CATATG|2|4|0|IV\nBsmI|GAATGC|1|-1|1|JMNS\nPsp1406I|AACGTT|2|4|0|BK\nRsrII|CGGWCCG|2|5|0|NQX\nHaeII|RGCGCY|5|1|0|JKN\nBsiEI|CGRYCG|4|2|0|N\nEarI|CTCTTC|1|4|1|N\nAclWI|GGATC|4|5|1|I\nSnaBI|TACGTA|3|3|0|CKMN\nBseX3I|CGGCCG|1|5|0|IV\nBseBI|CCWGG|2|3|0|C\nBsp1720I|GCTNAGC|2|5|0|IV\nPspEI|GGTNACC|1|6|0|IV\nPaeI|GCATGC|5|1|0|B\nMssI|GTTTAAAC|4|4|0|B\nPmaCI|CACGTG|3|3|0|K\nBslFI|GGGAC|10|14|1|I\nPshAI|GACNNNNGTC|5|5|0|KN\nBseJI|GATNNNNATC|5|5|0|B\nEclXI|CGGCCG|1|5|0|S\nBseMI|GCAATG|2|0|1|B\nPscI|ACATGT|1|5|0|B\nPspGI|CCWGG|0|5|0|N\nCfr42I|CCGCGG|4|2|0|B\nBseSI|GKGCMC|5|1|0|B\nAarI|CACCTGC|4|8|1|B\nSfoI|GGCGCC|3|3|0|N\nBspFNI|CGCG|2|2|0|I\nAoxI|GGCC|0|4|0|I\nBpuMI|CCSGG|2|3|0|V\nPsp5II|RGGWCCY|2|5|0|B\nMph1103I|ATGCAT|5|1|0|B\nBpiI|GAAGAC|2|6|1|B\nNruI|TCGCGA|3|3|0|BCIJKMNQSX\nBsaBI|GATNNNNATC|5|5|0|N\nMslI|CAYNNNNRTG|5|5|0|N\nCsp6I|GTAC|1|3|0|B\nBse3DI|GCAATG|2|0|1|IV\nPfl23II|CGTACG|1|5|0|B\nNciI|CCSGG|2|3|0|JN\nAsiSI|GCGATCGC|5|3|0|IN\nBmeRI|GACNNNNNGTC|6|5|0|V\nHindII|GTYRAC|3|3|0|IMSV\nRsaI|GTAC|2|2|0|BCIJMNQRVXY\nAlw44I|GTGCAC|1|5|0|BJ\nBmrFI|CCNGG|2|3|0|V\nPflFI|GACNNNGTC|4|5|0|N\nHsp92II|CATG|4|0|0|R\nNaeI|GCCGGC|3|3|0|CKN\nPvuII|CAGCTG|3|3|0|BCIJKMNOQSVX\nBoxI|GACNNNNGTC|5|5|0|B\nMox20I|TGGCCA|3|3|0|I\nCciI|TCATGA|1|5|0|I\nBmuI|ACTGGG|5|4|1|I\nPauI|GCGCGC|1|5|0|B\nHinfI|GANTC|1|4|0|BCIJKMNOQRVXY\nAatII|GACGTC|5|1|0|BIKMNV\nHaeIII|GGCC|2|2|0|BIJKMNOQRSX\nBstZ17I|GTATAC|3|3|0|N\nAcyI|GRCGYC|2|4|0|J\nBmeT110I|CYCGRG|1|5|0|BK\nPasI|CCCWGGG|2|5|0|B\nFaeI|CATG|4|0|0|I\nXcmI|CCANNNNNNNNNTGG|8|7|0|N\nBalI|TGGCCA|3|3|0|BJKQX\nBstUI|CGCG|2|2|0|N\nMscI|TGGCCA|3|3|0|NO\nMlsI|TGGCCA|3|3|0|B\nBstZI|CGGCCG|1|5|0|R\nFnu4HI|GCNGC|2|3|0|N\nBciVI|GTATCC|6|5|1|N\nPagI|TCATGA|1|5|0|B\nEcoRV|GATATC|3|3|0|BCIJKMNOQRSVX\nBstPAI|GACNNNNGTC|5|5|0|IV\nAcoI|YGGCCR|1|5|0|I\nBme18I|GGWCC|1|4|0|IV\nMteI|GCGCNGCGC|4|5|0|I\nCfoI|GCGC|3|1|0|MRS\nSfiI|GGCCNNNNNGGCC|8|5|0|BCIJKMNOQVX\nVan91I|CCANNNNNTGG|7|4|0|BK\nBstMBI|GATC|0|4|0|IV\nAluI|AGCT|2|2|0|BCIJKMNOQRSVXY\nKspAI|GTTAAC|3|3|0|B\nTspGWI|ACGGA|11|9|1|QX\nAbaSI|C|11|9|1|N\nNdeII|GATC|0|4|0|JM\nEcoO109I|RGGNCCY|2|5|0|BJKN\nBtsIMutI|CAGTG|2|0|1|N\nCviJI|RGCY|2|2|0|QX\nBstBAI|YACGTR|3|3|0|IV\nTauI|GCSGC|4|1|0|B\nXmaJI|CCTAGG|1|5|0|B\nBfmI|CTRYAG|1|5|0|B\nMaeIII|GTNAC|0|5|0|MS\nBstKTI|GATC|3|1|0|I\nMwoI|GCNNNNNNNGC|7|4|0|N\nHpyCH4V|TGCA|2|2|0|N\nAseI|ATTAAT|2|4|0|JN\nAccII|CGCG|2|2|0|JK\nHpaI|GTTAAC|3|3|0|BCIJKMNQSVX\nBstAUI|TGTACA|1|5|0|IV\nDdeI|CTNAG|1|4|0|KMNOQRSX\nBsgI|GTGCAG|16|14|1|N\nNmeAIII|GCCGAG|21|19|1|N\nMspCI|CTTAAG|1|5|0|C\nBsuI|GTATCC|6|5|1|I\nBssNAI|GTATAC|3|3|0|IV\nBstH2I|RGCGCY|5|1|0|IV\nHpy99I|CGWCG|5|0|0|N\nTspDTI|ATGAA|11|9|1|QX\nVneI|GTGCAC|1|5|0|IV\nAxyI|CCTNAGG|2|5|0|J\nKflI|GGGWCCC|2|5|0|B\nSetI|ASST|4|0|0|I\nBstNI|CCWGG|2|3|0|N\nNgoMIV|GCCGGC|1|5|0|N\nGlaI|GCGC|2|2|0|I\nBfrI|CTTAAG|1|5|0|MS\nApoI|RAATTY|1|5|0|N\nMmeI|TCCRAC|20|18|1|NX\nMroI|TCCGGA|1|5|0|MO\nFokI|GGATG|9|13|1|BIJKMNVXY\nPciSI|GCTCTTC|1|4|1|I\nDrdI|GACNNNNNNGTC|7|5|0|N\nPstNI|CAGNNNCTG|6|3|0|I\nSapI|GCTCTTC|1|4|1|N\nMva1269I|GAATGC|1|-1|1|B\nTru9I|TTAA|1|3|0|IMV\nAma87I|CYCGRG|1|5|0|IV\nHpy188III|TCNNGA|2|4|0|N\nBspMAI|CTGCAG|5|1|0|X\nFaiI|YATR|2|2|0|I\nAgeI|ACCGGT|1|5|0|JNR\nBcuI|ACTAGT|1|5|0|B\nBsaWI|WCCGGW|1|5|0|N\nBssECI|CCNNGG|1|5|0|I\nNcoI|CCATGG|1|5|0|BCJKMNOQRSXY\nEco31I|GGTCTC|1|5|1|B\nHpyAV|CCTTC|6|5|1|N\nMly113I|GGCGCC|2|4|0|I\nXmaI|CCCGGG|1|5|0|INV\nBseMII|CTCAG|10|8|1|B\nAciI|CCGC|1|3|0|N\nLweI|GCATC|5|9|1|B\nBstXI|CCANNNNNNTGG|8|4|0|BIJKMNQVXY\nPcsI|WCGNNNNNNNCGW|7|6|0|I\nPaqCI|CACCTGC|4|8|1|N\nTatI|WGTACW|1|5|0|B\nFblI|GTMKAC|2|4|0|IV\nAspLEI|GCGC|3|1|0|IV\nAsuNHI|GCTAGC|1|5|0|IV\nEgeI|GGCGCC|3|3|0|I\nMspI|CCGG|1|3|0|BIJKNQRVX\nAclI|AACGTT|2|4|0|INV\nBglII|AGATCT|1|5|0|BCIJKMNOQRSVX\nBsaJI|CCNNGG|1|5|0|N\nOliI|CACNNNNGTG|5|5|0|B\nAlwI|GGATC|4|5|1|N\nBccI|CCATC|4|5|1|N\nEco57I|CTGAAG|16|14|1|B\nMauBI|CGCGCGCG|2|6|0|B\nPspFI|CCCAGC|1|5|0|B\nLguI|GCTCTTC|1|4|1|B\nBpmI|CTGGAG|16|14|1|IN\nXbaI|TCTAGA|1|5|0|BCIJKMNQRSVX\nSspMI|CTAG|1|3|0|I\nEcoT14I|CCWWGG|1|5|0|K\nBsp1286I|GDGCHC|5|1|0|JKN\nHpyF10VI|GCNNNNNNNGC|7|4|0|B\nMluCI|AATT|0|4|0|N\nBmcAI|AGTACT|3|3|0|V\nEco47III|AGCGCT|3|3|0|BMS\nEcoICRI|GAGCTC|3|3|0|IV\nAccIII|TCCGGA|1|5|0|JK\nBcnI|CCSGG|2|3|0|BK\nMaeI|CTAG|1|3|0|M\nFauI|CCCGC|4|6|1|IN\nAcuI|CTGAAG|16|14|1|IN\nZrmI|AGTACT|3|3|0|I\nTspMI|CCCGGG|1|5|0|N\nAgsI|TTSAA|3|2|0|I\nSlaI|CTCGAG|1|5|0|CY\nTth111I|GACNNNGTC|4|5|0|IKNQVX\nEco130I|CCWWGG|1|5|0|B\nBsc4I|CCNNNNNNNGG|7|4|0|I\nFriOI|GRGCYC|5|1|0|IV\nAsp700I|GAANNNNTTC|5|5|0|M\nDpnI|GATC|2|2|0|BEKMNOQRSX\nEcl136II|GAGCTC|3|3|0|B\nEco105I|TACGTA|3|3|0|B\nMfeI|CAATTG|1|5|0|IN\nBseAI|TCCGGA|1|5|0|C\nBspTI|CTTAAG|1|5|0|B\nKsp22I|TGATCA|1|5|0|IV\nXmiI|GTMKAC|2|4|0|B\nStuI|AGGCCT|3|3|0|BJKMNQX\nMbiI|CCGCTC|3|3|0|B\nEam1104I|CTCTTC|1|4|1|B\nGsaI|CCCAGC|5|1|0|I\nSse9I|AATT|0|4|0|IV\nAvaI|CYCGRG|1|5|0|JNQX\nFatI|CATG|0|4|0|IN\nBstAFI|CTTAAG|1|5|0|I\nSgrDI|CGTCGACG|2|6|0|B\nTseI|GCWGC|1|4|0|N\nEco88I|CYCGRG|1|5|0|B\nSacII|CCGCGG|4|2|0|BJKNOQRX\nBanII|GRGCYC|5|1|0|KNX\nEam1105I|GACNNNNNGTC|6|5|0|B\nAfeI|AGCGCT|3|3|0|IN\nKpn2I|TCCGGA|1|5|0|B\nBsa29I|ATCGAT|2|4|0|I\nBspDI|ATCGAT|2|4|0|N\nVpaK11BI|GGWCC|1|4|0|K\nBbv12I|GWGCWC|5|1|0|IV\nBmgBI|CACGTC|3|3|0|N\nBveI|ACCTGC|4|8|1|B\nSrfI|GCCCGGGC|4|4|0|N\nSgrAI|CRCCGGYG|2|6|0|N\nApeKI|GCWGC|1|4|0|N\nSse8387I|CCTGCAGG|6|2|0|K\nEagI|CGGCCG|1|5|0|N\nBssNI|GRCGYC|2|4|0|V\nSfr274I|CTCGAG|1|5|0|IV\nStyD4I|CCNGG|0|5|0|N\nEco47I|GGWCC|1|4|0|B\nPstI|CTGCAG|5|1|0|BCIJKMNOQRSVX\nAhdI|GACNNNNNGTC|6|5|0|N\nDriI|GACNNNNNGTC|6|5|0|I\nHspAI|GCGC|1|3|0|IV\nBspMI|ACCTGC|4|8|1|N\nAcvI|CACGTG|3|3|0|QX\nBlnI|CCTAGG|1|5|0|KMS\nBsp143I|GATC|0|4|0|B\nBstSCI|CCNGG|0|5|0|I\nSmoI|CTYRAG|1|5|0|B\nSfaAI|GCGATCGC|5|3|0|B\nAlw21I|GWGCWC|5|1|0|B\nSmiI|ATTTAAAT|4|4|0|BIKV\nSalI|GTCGAC|1|5|0|BCIJKMNOQRVX\nAccBSI|CCGCTC|3|3|0|IV\nBstV1I|GCAGC|8|12|1|I\nBsrFI|RCCGGY|1|5|0|N\nHpySE526I|ACGT|1|3|0|I\nAflIII|ACRYGT|1|5|0|MNS\nSgrBI|CCGCGG|4|2|0|C\nDpnII|GATC|0|4|0|N\nBssHII|GCGCGC|1|5|0|JKMNQX\nRsaNI|GTAC|1|3|0|I\nSmlI|CTYRAG|1|5|0|N\nCsiI|ACCWGGT|1|6|0|B\nNspI|RCATGY|5|1|0|N\nZsp2I|ATGCAT|5|1|0|IV\nBstSLI|GKGCMC|5|1|0|I\nAanI|TTATAA|3|3|0|B\nHpaII|CCGG|1|3|0|BINQRVX\nAspA2I|CCTAGG|1|5|0|IV\nBsp119I|TTCGAA|2|4|0|B\nBstENI|CCTNNNNNAGG|5|6|0|IV\nSfcI|CTRYAG|1|5|0|N\nRigI|GGCCGGCC|6|2|0|I\nAdeI|CACNNNGTG|6|3|0|B\nBsiWI|CGTACG|1|5|0|N\nCfr10I|RCCGGY|1|5|0|BK\nBbsI|GAAGAC|2|6|1|N\nBst6I|CTCTTC|1|4|1|IV\nSfaNI|GCATC|5|9|1|INV\nPdiI|GCCGGC|3|3|0|B\nRruI|TCGCGA|3|3|0|B\nPaeR7I|CTCGAG|1|5|0|N\nHin1I|GRCGYC|2|4|0|BK\nLmnI|GCTCC|1|-1|1|I\nCviAII|CATG|1|3|0|N\nBsrGI|TGTACA|1|5|0|N\nPsuI|RGATCY|1|5|0|B\nSexAI|ACCWGGT|1|6|0|MN\nBtgI|CCRYGG|1|5|0|N\nNlaIII|CATG|4|0|0|N\nTaiI|ACGT|4|0|0|B\nBstMCI|CGRYCG|4|2|0|IV\nHinP1I|GCGC|1|3|0|N\nBspACI|CCGC|1|3|0|I\nXmnI|GAANNNNTTC|5|5|0|N\nBsiSI|CCGG|1|3|0|CY\nBsiHKCI|CYCGRG|1|5|0|QX\nBstDEI|CTNAG|1|4|0|IV\nRsr2I|CGGWCCG|2|5|0|IV\nPvuI|CGATCG|4|2|0|BKMNOQRSX\nAasI|GACNNNNNNGTC|7|5|0|B\nFbaI|TGATCA|1|5|0|K\nPspCI|CACGTG|3|3|0|IV\nBseCI|ATCGAT|2|4|0|C\nCviKI_1|RGCY|2|2|0|N\nBsnI|GGCC|2|2|0|V\nMlyI|GAGTC|5|5|1|N\nBpu10I|CCTNAGC|2|5|0|BINV\nBsmBI|CGTCTC|1|5|1|N\nBspQI|GCTCTTC|1|4|1|N\nBsu15I|ATCGAT|2|4|0|B\nBsrDI|GCAATG|2|0|1|N\nEciI|GGCGGA|11|9|1|N\nPspLI|CGTACG|1|5|0|I\nSau96I|GGNCC|1|4|0|JN\nHhaI|GCGC|3|1|0|BJKNQRX\nBsiHKAI|GWGCWC|5|1|0|N\nBseYI|CCCAGC|1|5|0|N\nBspCNI|CTCAG|9|7|1|N\nSspI|AATATT|3|3|0|BCIJKNQVX\nBshTI|ACCGGT|1|5|0|B\nBseDI|CCNNGG|1|5|0|B\nBspT107I|GGYRCC|1|5|0|K\nPspPPI|RGGWCCY|2|5|0|I\nPle19I|CGATCG|4|2|0|I\nMvnI|CGCG|2|2|0|M\nPmeI|GTTTAAAC|4|4|0|N\nBcoDI|GTCTC|1|5|1|N\nBso31I|GGTCTC|1|5|1|IV\nPsiI|TTATAA|3|3|0|IN\nBsh1236I|CGCG|2|2|0|B\nEco52I|CGGCCG|1|5|0|BK\nBseNI|ACTGG|1|-1|1|B\nPshBI|ATTAAT|2|4|0|K\nPspPI|GGNCC|1|4|0|C\nFseI|GGCCGGCC|6|2|0|N\nBsh1285I|CGRYCG|4|2|0|B\nAsiGI|ACCGGT|1|5|0|IV\nBtgZI|GCGATG|10|14|1|N\nAcc36I|ACCTGC|4|8|1|I\nSmaI|CCCGGG|3|3|0|BCIJKMNOQRSVXY\nBsePI|GCGCGC|1|5|0|IV\nBse21I|CCTNAGG|2|5|0|IV\nPsp6I|CCWGG|0|5|0|I\nPacI|TTAATTAA|5|3|0|BNO\nBseXI|GCAGC|8|12|1|B\nPmlI|CACGTG|3|3|0|N\nBse8I|GATNNNNATC|5|5|0|IV\nPdmI|GAANNNNTTC|5|5|0|B\nMsp20I|TGGCCA|3|3|0|V\nCspAI|ACCGGT|1|5|0|C\nBseGI|GGATG|2|0|1|B\nPinAI|ACCGGT|1|5|0|QX\nPpuMI|RGGWCCY|2|5|0|N\nBmtI|GCTAGC|5|1|0|INV\nBseLI|CCNNNNNNNGG|7|4|0|B\nScaI|AGTACT|3|3|0|BCJKMNOQRX\nBspANI|GGCC|2|2|0|X\nAor13HI|TCCGGA|1|5|0|K\nBpu1102I|GCTNAGC|2|5|0|BK\nPfoI|TCCNGGA|1|6|0|B\nKspI|CCGCGG|4|2|0|MS\nBmsI|GCATC|5|9|1|B\nNlaIV|GGNNCC|3|3|0|N\nMroXI|GAANNNNTTC|5|5|0|IV\nCciNI|GCGGCCGC|2|6|0|IV\nBse1I|ACTGG|1|-1|1|IV\nPciI|ACATGT|1|5|0|IN\nMvaI|CCWGG|2|3|0|BM\nApaI|GGGCCC|5|1|0|BIJKMNQRVX\nBlsI|GCNGC|3|2|0|I\nAhlI|ACTAGT|1|5|0|IV\nHincII|GTYRAC|3|3|0|BJKNOQX\nBmgT120I|GGNCC|1|4|0|K\nPfeI|GAWTC|1|4|0|B\nHin1II|CATG|4|0|0|B\nBsaAI|YACGTR|3|3|0|N\nMspA1I|CMGCKG|3|3|0|INV\nBmiI|GGNNCC|3|3|0|V\nMluNI|TGGCCA|3|3|0|M\nBsuTUI|ATCGAT|2|4|0|X\nBfuI|GTATCC|6|5|1|B\nPalAI|GGCGCGCC|2|6|0|I\nFsp4HI|GCNGC|2|3|0|I\nFspI|TGCGCA|3|3|0|JN\nBstSNI|TACGTA|3|3|0|IV\nAcsI|RAATTY|1|5|0|IV\nBme1390I|CCNGG|2|3|0|B\nNmuCI|GTSAC|0|5|0|B\nEcoT22I|ATGCAT|5|1|0|BK\nTspRI|CASTG|2|-7|1|N\nAluBI|AGCT|2|2|0|I\nBstC8I|GCNNGC|3|3|0|I\nHpy8I|GTNNAC|3|3|0|B\nMalI|GATC|2|2|0|I\nBstX2I|RGATCY|1|5|0|IV\nAsuHPI|GGTGA|8|7|1|IV\nNspV|TTCGAA|2|4|0|J\nEcoRII|CCWGG|0|5|0|BJ\nAbsI|CCTCGAGG|2|6|0|I\nBtsCI|GGATG|2|0|1|N\nDraI|TTTAAA|3|3|0|BIJKMNQSVX\nBstFNI|CGCG|2|2|0|IV\nBisI|GCNGC|2|3|0|I\nMspR9I|CCNGG|2|3|0|I\nBstNSI|RCATGY|5|1|0|IV\nPflMI|CCANNNNNTGG|7|4|0|N\nTscAI|CASTG|2|-7|1|B\nBamHI|GGATCC|1|5|0|BCIJKMNOQRSVXY\nAleI|CACNNNNGTG|5|5|0|N\nKroNI|GCCGGC|3|3|0|I\nBstBI|TTCGAA|2|4|0|N\nTaqII|GACCGA|11|9|1|QX\nNarI|GGCGCC|2|4|0|JMNQX\nEcoNI|CCTNNNNNAGG|5|6|0|N\nBtsI|GCAGTG|2|0|1|N\nCac8I|GCNNGC|3|3|0|N\nBst1107I|GTATAC|3|3|0|BK\nXapI|RAATTY|1|5|0|B\nBciT130I|CCWGG|2|3|0|K\nMabI|ACCWGGT|1|6|0|I\nBstHHI|GCGC|3|1|0|IV\nHpy188I|TCNGA|3|2|0|N\nTaaI|ACNGT|3|2|0|B\nNheI|GCTAGC|1|5|0|BCJKMNOQRSX\nHpy166II|GTNNAC|3|3|0|N\nAscI|GGCGCGCC|2|6|0|N\nBsu36I|CCTNAGG|2|5|0|N\nBseRI|GAGGAG|10|8|1|N\nMnlI|CCTC|7|6|1|BINQVX\nMroNI|GCCGGC|1|5|0|IV\nCfr13I|GGNCC|1|4|0|B\nPpsI|GAGTC|4|5|1|I\nBspLI|GGNNCC|3|3|0|B\nPctI|GAATGC|1|-1|1|IV\nVha464I|CTTAAG|1|5|0|V\nAsuC2I|CCSGG|2|3|0|I\nHpyF3I|CTNAG|1|4|0|B\nBspOI|GCTAGC|5|1|0|B\nEcoT38I|GRGCYC|5|1|0|J\nSduI|GDGCHC|5|1|0|B\nApaLI|GTGCAC|1|5|0|CKN\nNdeI|CATATG|2|4|0|BJKMNQRSX\nFspAI|RTGCGCAY|4|4|0|B\nBfaI|CTAG|1|3|0|N\nBstEII|GGTNACC|1|6|0|CJN\nEsp3I|CGTCTC|1|5|1|BN\nMboII|GAAGA|8|7|1|BIJKNQVX\nMreI|CGCCGGCG|2|6|0|B\nHgaI|GACGC|5|10|1|IN\nMspJI|CNNR|9|13|1|N\nPleI|GAGTC|4|5|1|N\nTru1I|TTAA|1|3|0|B\nAccB1I|GGYRCC|1|5|0|IV\nGluI|GCNGC|2|3|0|I\nBfoI|RGCGCY|5|1|0|B\nDraIII|CACNNNGTG|6|3|0|IMNV\nPkrI|GCNGC|3|2|0|I\nMunI|CAATTG|1|5|0|BKMS\nEheI|GGCGCC|3|3|0|B\nAflII|CTTAAG|1|5|0|JKN\nBsaHI|GRCGYC|2|4|0|N\nAvrII|CCTAGG|1|5|0|N\nBsoBI|CYCGRG|1|5|0|N\nAlw26I|GTCTC|1|5|1|B\nBceAI|ACGGC|12|14|1|N\nHphI|GGTGA|8|7|1|BN\nMflI|RGATCY|1|5|0|K\nXhoI|CTCGAG|1|5|0|BJKMNOQRX\nBpuEI|CTTGAG|16|14|1|N\nSsiI|CCGC|1|3|0|B\nLpnPI|CCDG|10|14|1|N\nBstF5I|GGATG|2|0|1|IV\nLsp1109I|GCAGC|8|12|1|B\nTasI|AATT|0|4|0|B\nErhI|CCWWGG|1|5|0|I\nBstAPI|GCANNNNNTGC|7|4|0|IN\nMhlI|GDGCHC|5|1|0|IV";
  const METH_DEP = ["MspJI", "LpnPI", "FspEI", "AbaSI", "SgeI"];
  const _ENZ = [];
  REBASE_PACKED.split("\n").forEach(function (line) {
    const p = line.split("|");
    if (p.length < 5) return;
    const name = p[0];
    if (METH_DEP.indexOf(name) !== -1) return;
    const site = p[1], a = +p[2], b = +p[3], iis = p[4] === "1", sup = p[5] || "";
    const rec = iis
      ? { name: name, site: site, o1: a, o2: b, typeIIS: true, suppliers: sup }
      : { name: name, site: site, top: a, bot: b, typeIIS: false, suppliers: sup };
    _ENZ.push(rec);
  });
  // common aliases / isoschizomers people type that REBASE lists under one name
  const ALIASES = {
    "BsmBI": "Esp3I", "PaqCI": "AarI", "Eco31I": "BsaI",
    "BbsI-HF": "BbsI", "LguI": "SapI", "BspQI": "SapI", "BsaI-HFv2": "BsaI",
    "Esp3I": null, "AarI": null
  };
  Object.keys(ALIASES).forEach(function (alias) {
    if (_ENZ.some(function (e) { return e.name === alias; })) return;
    const target = ALIASES[alias];
    if (!target) return;
    const base = _ENZ.find(function (e) { return e.name === target; });
    if (base) { const c = Object.assign({}, base, { name: alias, alias: target }); _ENZ.push(c); }
  });
  E.ALL_ENZYMES = _ENZ;
  // REBASE supplier codes and names as distributed with Biopython 1.88
  E.SUPPLIER_NAMES = {"B": "Thermo Fisher Scientific", "C": "Minotech Biotechnology", "E": "Agilent Technologies", "I": "SibEnzyme Ltd.", "J": "Nippon Gene Co., Ltd.", "K": "Takara Bio Inc.", "M": "Roche Applied Science", "N": "New England Biolabs", "O": "Toyobo Biochemicals", "Q": "Molecular Biology Resources - CHIMERx", "R": "Promega Corporation", "S": "Sigma Chemical Corporation", "V": "Vivantis Technologies", "X": "EURx Ltd.", "Y": "SinaClon BioScience Co."};

  E.enzymeByName = function (name) {
    return E.ALL_ENZYMES.find(e => e.name.toLowerCase() === name.toLowerCase());
  };

  // find all sites for one enzyme. circular wraps. returns list of
  // {pos (0-based start of recognition on top), cutTop, cutBottom, strand}
  /** @param {string} seq @param {Enzyme} enz @param {boolean} [circular] @returns {Site[]} */
  E.findSites = function (seq, enz, circular) {
    seq = E.cleanSeq(seq);
    const L = seq.length;
    const search = circular ? seq + seq.substring(0, Math.max(enz.site.length, 20)) : seq;
    const sites = [];
    const site = enz.site;
    const rc = E.revcomp(site);
    const palindrome = site === rc;

    function scan(pattern, strand) {
      for (let i = 0; i + pattern.length <= search.length; i++) {
        let ok = true;
        for (let j = 0; j < pattern.length; j++) {
          if (!matchCode(pattern[j], search[i + j])) { ok = false; break; }
        }
        if (!ok) continue;
        if (i >= L) continue; // only report starts within sequence
        const cut = computeCut(enz, i, strand, L, circular);
        if (cut) sites.push(cut);
      }
    }
    scan(site, 1);
    if (!palindrome) scan(rc, -1);
    // dedupe by cutTop
    const seen = new Set();
    const uniq = [];
    for (const s of sites.sort((a, b) => a.pos - b.pos)) {
      const key = s.strand + ':' + s.pos;
      if (seen.has(key)) continue;
      seen.add(key); uniq.push(s);
    }
    return uniq;
  };
  function matchCode(code, base) {
    if (code === base) return true;
    const set = IUPAC[code];
    return set ? set.indexOf(base) !== -1 : false;
  }
  function computeCut(enz, pos, strand, L, circular) {
    let cutTop, cutBottom;
    if (!enz.typeIIS) {
      if (strand === 1) { cutTop = pos + enz.top; cutBottom = pos + enz.bot; }
      else {
        // site found on bottom strand; recognition occupies pos..pos+len-1
        const len = enz.site.length;
        cutTop = pos + len - enz.bot;
        cutBottom = pos + len - enz.top;
      }
    } else {
      const len = enz.site.length;
      if (strand === 1) { cutTop = pos + len + enz.o1; cutBottom = pos + len + enz.o2; }
      else { cutTop = pos - enz.o2; cutBottom = pos - enz.o1; }
    }
    if (circular) {
      cutTop = ((cutTop % L) + L) % L;
      cutBottom = ((cutBottom % L) + L) % L;
    } else {
      // A linear molecule has scissile bonds only between adjacent bases, i.e.
      // at cut positions 1..L-1 on each strand. A recognition site whose cut
      // would fall at or beyond either end (position <=0 or >=L) does not
      // produce a cut. Biopython drops these; so do we.
      if (cutTop < 1 || cutTop > L - 1 || cutBottom < 1 || cutBottom > L - 1) return null;
    }
    return { pos, strand, cutTop, cutBottom, enzyme: enz.name };
  }

  // Digest: returns fragments with sizes and overhang descriptions.
  // opts {circular}
  E.digest = function (seq, enzymeNames, opts) {
    opts = opts || {};
    seq = E.cleanSeq(seq);
    const L = seq.length;
    const enzymes = enzymeNames.map(n => typeof n === 'string' ? E.enzymeByName(n) : n).filter(Boolean);
    let cuts = [];
    const blocked = [], caution = [];
    // opts.methylation: {dam, dcm, cpg} — sites the chart says are blocked by
    // methylation of this DNA are left uncut; uncertain ones are reported
    const meth = opts.methylation && E.methylatedBases
      ? E.methylatedBases(seq, opts.circular, opts.methylation) : null;
    for (const enz of enzymes) {
      for (const s of E.findSites(seq, enz, opts.circular)) {
        const m = meth ? E.methylationEffect(enz, s, L, meth) : null;
        if (m && m.effect === 'blocked') { blocked.push({ enzyme: enz.name, pos: s.pos, by: m.by }); continue; }
        if (m) caution.push({ enzyme: enz.name, pos: s.pos, by: m.by, effect: m.effect });
        cuts.push({ top: s.cutTop, bottom: s.cutBottom, enzyme: enz.name });
      }
    }
    // sort by top cut
    cuts.sort((a, b) => a.top - b.top);
    if (cuts.length === 0) {
      return { fragments: [{ seq, start: 0, end: L, length: L, leftEnd: 'blunt', rightEnd: 'blunt', circular: !!opts.circular, note: 'uncut' }], cuts: [], blocked, caution };
    }
    const fragments = [];
    if (opts.circular) {
      for (let i = 0; i < cuts.length; i++) {
        const c1 = cuts[i];
        const c2 = cuts[(i + 1) % cuts.length];
        const fseq = sliceCircular(seq, c1.top, c2.top, L);
        fragments.push({
          start: c1.top, end: c2.top, length: fseq.length, seq: fseq,
          leftEnd: endDesc(c1, seq, L), rightEnd: endDesc(c2, seq, L),
          leftEnzyme: c1.enzyme, rightEnzyme: c2.enzyme
        });
      }
    } else {
      let prev = { top: 0, bottom: 0, enzyme: null };
      const bounded = cuts.filter(c => c.top > 0 && c.top < L);
      let cursor = 0;
      const points = [{ top: 0, bottom: 0, enzyme: 'end' }].concat(bounded).concat([{ top: L, bottom: L, enzyme: 'end' }]);
      for (let i = 0; i < points.length - 1; i++) {
        const c1 = points[i], c2 = points[i + 1];
        const fseq = seq.substring(c1.top, c2.top);
        fragments.push({
          start: c1.top, end: c2.top, length: fseq.length, seq: fseq,
          leftEnd: c1.enzyme === 'end' ? 'blunt(end)' : endDesc(c1, seq, L),
          rightEnd: c2.enzyme === 'end' ? 'blunt(end)' : endDesc(c2, seq, L),
          leftEnzyme: c1.enzyme, rightEnzyme: c2.enzyme
        });
      }
    }
    return { fragments, cuts, blocked, caution };
  };
  function sliceCircular(seq, from, to, L) {
    if (to > from) return seq.substring(from, to);
    return seq.substring(from) + seq.substring(0, to);
  }
  function endDesc(cut, seq, L) {
    const d = cut.bottom - cut.top;
    if (d === 0) return 'blunt';
    const type = d > 0 ? "5'" : "3'";
    const len = Math.abs(d);
    let oh;
    if (d > 0) oh = seq.substr(((cut.top % L) + L) % L, len);
    else oh = seq.substr(((cut.bottom % L) + L) % L, len);
    return type + ' ' + oh + ' (' + len + ')';
  }

  // fragment list with just sizes (for gel), sorted desc
  E.digestSizes = function (seq, enzymeNames, opts) {
    return E.digest(seq, enzymeNames, opts).fragments.map(f => f.length).sort((a, b) => b - a);
  };

  // Enzyme survey: cutters count across sequence. returns map name->count
  E.surveyEnzymes = function (seq, opts) {
    opts = opts || {};
    const list = opts.enzymes || E.ALL_ENZYMES;
    const out = [];
    for (const enz of list) {
      const n = E.findSites(seq, enz, opts.circular).length;
      out.push({ name: enz.name, site: enz.site, count: n, typeIIS: enz.typeIIS });
    }
    return out;
  };

  /* ---------- Pairwise alignment ----------
     Needleman-Wunsch (global) and Smith-Waterman (local), affine-ish simple
     linear gap. For DNA default match/mismatch; returns aligned strings +
     identity. This is a real DP implementation (O(nm)).
  ------------------------------------------------------------------ */
  E.align = function (a, b, opts) {
    opts = opts || {};
    const match = opts.match != null ? opts.match : 2;
    const mismatch = opts.mismatch != null ? opts.mismatch : -1;
    const gap = opts.gap != null ? opts.gap : -2;
    const local = opts.mode === 'local';
    a = (a || '').toUpperCase(); b = (b || '').toUpperCase();
    const n = a.length, m = b.length;
    const H = new Int32Array((n + 1) * (m + 1));
    const T = new Uint8Array((n + 1) * (m + 1)); // 0 diag,1 up,2 left,3 stop
    const W = m + 1;
    if (!local) {
      for (let i = 1; i <= n; i++) { H[i * W] = gap * i; T[i * W] = 1; }
      for (let j = 1; j <= m; j++) { H[j] = gap * j; T[j] = 2; }
    }
    /* A caller may supply its own substitution score — a protein alignment
       scored +2/-1 is a poor alignment, because it cannot tell a
       conservative substitution from a disruptive one. engine6 provides a
       BLOSUM62 scorer; without one, the flat match/mismatch is used. */
    const scorer = typeof opts.scorer === 'function' ? opts.scorer : null;
    let best = 0, bi = 0, bj = 0;
    for (let i = 1; i <= n; i++) {
      for (let j = 1; j <= m; j++) {
        const sc = scorer ? scorer(a[i - 1], b[j - 1])
                          : ((a[i - 1] === b[j - 1]) ? match : mismatch);
        const diag = H[(i - 1) * W + (j - 1)] + sc;
        const up = H[(i - 1) * W + j] + gap;
        const left = H[i * W + (j - 1)] + gap;
        let v = diag, t = 0;
        if (up > v) { v = up; t = 1; }
        if (left > v) { v = left; t = 2; }
        if (local && v < 0) { v = 0; t = 3; }
        H[i * W + j] = v; T[i * W + j] = t;
        if (local && v > best) { best = v; bi = i; bj = j; }
      }
    }
    let i, j, score;
    if (local) { i = bi; j = bj; score = best; }
    else { i = n; j = m; score = H[n * W + m]; }
    let A = '', B = '';
    while (i > 0 || j > 0) {
      if (local && (H[i * W + j] === 0)) break;
      const t = T[i * W + j];
      if (i > 0 && j > 0 && t === 0) { A = a[i - 1] + A; B = b[j - 1] + B; i--; j--; }
      else if (i > 0 && t === 1) { A = a[i - 1] + A; B = '-' + B; i--; }
      else if (j > 0 && t === 2) { A = '-' + A; B = b[j - 1] + B; j--; }
      else break;
    }
    let ident = 0, alnLen = A.length;
    let midline = '';
    for (let k = 0; k < A.length; k++) {
      if (A[k] === B[k] && A[k] !== '-') { ident++; midline += '|'; }
      else midline += ' ';
    }
    return { a: A, b: B, midline, score, identity: alnLen ? ident / alnLen * 100 : 0, length: alnLen, matches: ident };
  };

  /* ---------- PCR simulation ----------
     Given template (linear or circular), fwd and rev primers (5'->3'),
     find binding sites (allow 3' exact match of a 3' anchor of given length)
     and produce the amplicon including primer tails.
  ------------------------------------------------------------------ */
  /**
   * @param {string} template @param {string} fwd @param {string} rev
   * @param {{anchor?:number, circular?:boolean, maxPairs?:number, maxLen?:number, maxProducts?:number}} [opts]
   * @returns {PcrResult} the products, with `found`, `sites` and `truncated` on the array
   */
  E.pcr = function (template, fwd, rev, opts) {
    opts = opts || {};
    const anchor = opts.anchor != null ? opts.anchor : 15; // 3' bases that must match
    const circular = !!opts.circular;
    template = E.cleanSeq(template);
    fwd = E.cleanSeq(fwd); rev = E.cleanSeq(rev);
    const L = template.length;
    const search = circular ? template + template : template;

    // fwd primer: its 3' end anchors on top strand
    const fAnchorSeq = fwd.slice(-anchor);
    const fIdx = indexOfAll(search, fAnchorSeq).filter(x => x < L);
    // rev primer: anchors on bottom strand => its revcomp appears on top strand
    const rAnchorRC = E.revcomp(rev.slice(-anchor));
    const rIdx = indexOfAll(search, rAnchorRC).filter(x => x < L + (circular ? L : 0));

    /* Pairs are collected as coordinates first and the amplicon sequence is
       built only for the ones that are kept. On an ordinary template there
       are a handful of sites and this makes no difference, but a repetitive
       or low-complexity template can give a primer a thousand sites, and a
       thousand by a thousand amplicons built as strings is hundreds of
       megabytes - enough to take the tab down. */
    const pairCap = opts.maxPairs || 200000;
    const products = [];
    let overflowed = false;
    for (const f of fIdx) {
      // fwd 3' end is at position f + anchor - 1 (on top). primer extends rightward.
      const fStart = f + anchor - fwd.length; // where full primer (with tail) would align (tail may hang off)
      if (overflowed) break;
      for (const r of rIdx) {
        // rev anchor rc starts at r; the amplicon end (top strand) is r+anchor
        const ampEndTop = r + anchor;
        const primeStartTop = f; // 3'-anchor start
        if (ampEndTop <= primeStartTop) continue;
        let coreStart = f;            // first matched base of fwd anchor
        let coreEnd = r + anchor;     // last matched base (exclusive) of rev anchor rc
        /* A circular template is searched twice over so that a product may
           cross the origin, but a product can never be longer than the
           plasmid itself: once the polymerase reaches the reverse site it
           stops, so a second lap is not a PCR product. Without this the
           same pair is reported again at length + L. */
        if (coreEnd - coreStart > L) continue;
        if (coreEnd - coreStart > (opts.maxLen || 50000)) continue;
        /* the length without building the string: the primers plus what
           lies between their matched anchors */
        const internalLen = Math.max(0, r - (f + anchor));
        products.push({
          start: coreStart % L, end: coreEnd % L,
          length: fwd.length + internalLen + rev.length,
          from: f + anchor, to: r,
          fwdSite: f % L, revSite: r % L
        });
        if (products.length >= pairCap) { overflowed = true; break; }
      }
    }
    products.sort((a, b) => a.length - b.length);
    /* the same pair of sites found once per lap is one product */
    const seenPair = {}, uniq = [];
    products.forEach(p => {
      const k = p.fwdSite + ':' + p.revSite + ':' + p.length;
      if (seenPair[k]) return;
      seenPair[k] = 1; uniq.push(p);
    });
    /* the shortest products are the ones worth having; anything past the
       cap is reported as a count rather than built */
    const cap = opts.maxProducts || 500;
    const keep = uniq.slice(0, cap);
    keep.forEach(p => {
      // build amplicon: fwd full primer + internal template + rev full primer(revcomp)
      const internal = sliceLinear(search, p.from, p.to, circular, L);
      p.seq = fwd + internal + E.revcomp(rev);
      p.length = p.seq.length;
      delete p.from; delete p.to;
    });
    keep.sort((a, b) => a.length - b.length);
    /* The products, with what the search found attached to the array; see
       the PcrResult typedef at the top of this file. */
    const result = /** @type {PcrResult} */ (keep);
    result.truncated = overflowed || uniq.length > keep.length;
    result.found = uniq.length;
    /* Distinct sites on the plasmid, not matches in the doubled sequence
       the circular search runs over: a primer with one site would
       otherwise be reported as having two. */
    const distinct = idx => { const seen = {}; let n = 0; idx.forEach(i => { const k = i % L; if (!seen[k]) { seen[k] = 1; n++; } }); return n; };
    result.sites = { forward: distinct(fIdx), reverse: distinct(rIdx) };
    return result;
  };
  function indexOfAll(hay, needle) {
    const out = []; if (!needle) return out;
    let i = hay.indexOf(needle);
    while (i !== -1) { out.push(i); i = hay.indexOf(needle, i + 1); }
    return out;
  }
  function sliceLinear(search, from, to, circular, L) {
    if (to >= from) return search.substring(from, to);
    if (circular) return search.substring(from, to + L);
    return '';
  }

  /* ---------- Gibson / overlap assembly primer design ----------
     Design primers to add homology arms for seamless assembly of an insert
     into a linearized vector (or chaining fragments).
     overlapLen default 20-40 bp. Returns primers with tails.
  ------------------------------------------------------------------ */
  E.gibsonJunction = function (upstreamSeq, downstreamSeq, overlap) {
    overlap = overlap || 20;
    // homology arm = last `overlap` of upstream + first `overlap` of downstream
    const up = upstreamSeq.slice(-overlap);
    const down = downstreamSeq.slice(0, overlap);
    return { overlapRegion: up + down, upTail: up, downTail: down };
  };

  // Design the 4 primers to clone insert into vector (linearized ends given)
  // vectorEndLeft = sequence at the 5' side of the linearized vector (the end
  // that will join insert 5'); vectorEndRight = 3' side.
  E.gibsonPrimers = function (insert, vecUpstream, vecDownstream, opts) {
    opts = opts || {};
    const overlap = opts.overlap || 20;
    const anneal = opts.anneal || 20;
    insert = E.cleanSeq(insert);
    const fwdTail = vecUpstream.slice(-overlap);
    const fwdAnneal = insert.slice(0, anneal);
    const revTail = E.revcomp(vecDownstream.slice(0, overlap));
    const revAnneal = E.revcomp(insert.slice(-anneal));
    return {
      insertF: { seq: fwdTail + fwdAnneal, tail: fwdTail, anneal: fwdAnneal,
                 tmAnneal: E.tmNN(fwdAnneal), length: fwdTail.length + fwdAnneal.length },
      insertR: { seq: revTail + revAnneal, tail: revTail, anneal: revAnneal,
                 tmAnneal: E.tmNN(E.revcomp(revAnneal)), length: revTail.length + revAnneal.length }
    };
  };

  /* ---------- Golden Gate ----------
     For a Type IIS enzyme, excise the fragment between the two cut sites and
     report the 4-nt fusion overhangs. Check overhang uniqueness across parts.
  ------------------------------------------------------------------ */
  E.goldenGateOverhangs = function (parts, enzymeName) {
    const enz = E.enzymeByName(enzymeName);
    if (!enz || !enz.typeIIS) throw new Error('Golden Gate needs a Type IIS enzyme (BsaI, BsmBI, BbsI...)');
    const results = [];
    for (const p of parts) {
      const d = E.digest(p.seq, [enzymeName], { circular: !!p.circular });
      results.push({ name: p.name, fragments: d.fragments, cuts: d.cuts });
    }
    return results;
  };

  /* ---------- Silent mutation / restriction site removal ----------
     For a coding region, find synonymous codon changes that remove or add a
     target restriction site without changing the protein.
  ------------------------------------------------------------------ */
  E.synonymousCodons = function (aa, table) {
    const { aaMap } = buildAAmap(table || 1);
    const out = [];
    for (const codon of CODON_ORDER) if (aaMap[codon] === aa) out.push(codon);
    return out;
  };

  /* ---------- Codon optimization ----------
     Re-code a protein (or CDS) to a host's preferred codons. Uses provided
     usage table (codon->fraction within its AA). Greedy "best codon".
  ------------------------------------------------------------------ */
  // opts.keepFirst leaves the start codon alone (alternative starts such as
  // GTG/TTG are read as Met by the ribosome). opts.avoid lists enzymes whose
  // sites must not be created: any new site is removed by trying the other
  // synonymous codons in it, best usage first, and if that fails by putting
  // the original codons back. Bases after the last whole codon are kept.
  E.codonOptimize = function (cds, usageTable, table, opts) {
    opts = opts || {};
    const { aaMap } = buildAAmap(table || 1);
    cds = E.cleanSeq(cds).replace(/U/g, 'T');
    const nC = Math.floor(cds.length / 3);
    const orig = [], codons = [], choices = [];
    for (let k = 0; k < nC; k++) {
      const codon = cds.substr(k * 3, 3);
      orig.push(codon);
      const aa = aaMap[codon];
      if (aa === undefined || aa === '*' || (k === 0 && opts.keepFirst)) { codons.push(codon); choices.push([codon]); continue; }
      const syn = E.synonymousCodons(aa, table).slice()
        .sort((a, b) => (usageTable[b] || 0) - (usageTable[a] || 0));
      codons.push(syn[0] || codon); choices.push(syn.length ? syn : [codon]);
    }
    const tail = cds.slice(nC * 3);
    const join = () => codons.join('') + tail;
    const avoid = (opts.avoid || []).map(n => typeof n === 'string' ? E.enzymeByName(n) : n).filter(Boolean);
    const siteStarts = (s, enz) => {
      const out = [], pats = [enz.site, E.revcomp(enz.site)];
      for (const pat of pats) for (let i = 0; i + pat.length <= s.length; i++) {
        let ok = true;
        for (let j = 0; j < pat.length; j++) if (!matchCode(pat[j], s[i + j])) { ok = false; break; }
        if (ok) out.push({ i, n: pat.length });
      }
      return out;
    };
    const before = {};
    avoid.forEach(e => { before[e.name] = siteStarts(cds, e).length; });
    for (let round = 0; round < 6; round++) {
      let changed = false;
      for (const enz of avoid) {
        let hits = siteStarts(join(), enz);
        if (hits.length <= before[enz.name]) continue;
        for (const h of hits) {
          const k0 = Math.floor(h.i / 3), k1 = Math.min(nC - 1, Math.floor((h.i + h.n - 1) / 3));
          const inside = s => siteStarts(s.slice(k0 * 3, (k1 + 1) * 3), enz).length > 0;
          if (!inside(join())) continue;
          // was this site already in the input at the same place? then leave it
          if (orig.slice(k0, k1 + 1).join('') === codons.slice(k0, k1 + 1).join('')) continue;
          let fixed = false;
          for (let k = k0; k <= k1 && !fixed; k++) {
            const keep = codons[k];
            for (const alt of choices[k]) {
              if (alt === keep) continue;
              codons[k] = alt;
              if (!inside(join())) { fixed = true; break; }
            }
            if (!fixed) codons[k] = keep;
          }
          if (!fixed) for (let k = k0; k <= k1; k++) codons[k] = orig[k];
          changed = true;
        }
      }
      if (!changed) break;
    }
    return join();
  };

  /* ---------- CRISPR guide finder ----------
     Find all NGG (SpCas9) PAM sites and extract 20-nt protospacers on both
     strands. Report GC, and simple self-complementarity flag.
  ------------------------------------------------------------------ */
  E.findGuides = function (seq, opts) {
    opts = opts || {};
    const pam = opts.pam || 'NGG';
    const gLen = opts.length || 20;
    seq = E.cleanSeq(seq);
    const L = seq.length;
    const guides = [];
    // on a circle a guide may straddle base 1, so search a sequence carrying
    // enough of the start again, then fold the positions back
    const pad = gLen + pam.length;
    const scan = (opts.circular && L > pad) ? seq + seq.slice(0, pad) : seq;
    const strands = [{ s: scan, dir: 1 }, { s: E.revcomp(scan), dir: -1 }];
    const seen = new Set();
    for (const { s, dir } of strands) {
      for (let i = gLen; i + pam.length <= s.length; i++) {
        const pamSeq = s.substr(i, pam.length);
        if (!matchPam(pam, pamSeq)) continue;
        const proto = s.substr(i - gLen, gLen);
        if (proto.length < gLen) continue;
        const gc = E.gcContent(proto);
        let st = dir === 1 ? i - gLen : s.length - i;
        if (st >= L) continue;                       // a copy of one already found
        st = ((st % L) + L) % L;
        const key = st + ':' + dir + ':' + proto;
        if (seen.has(key)) continue;
        seen.add(key);
        guides.push({
          protospacer: proto, pam: pamSeq, strand: dir,
          start: st,
          wrapped: st + gLen + pam.length > L,
          gc: +gc.toFixed(1),
          hasBsaI: /GGTCTC|GAGACC/.test(proto)
        });
      }
    }
    return guides;
  };
  function matchPam(pam, s) {
    if (pam.length !== s.length) return false;
    for (let k = 0; k < pam.length; k++) if (!matchCode(pam[k], s[k])) return false;
    return true;
  }
  function pamToRegexTail() { return null; }

  /* ---------- Intron splicing (GT-AG) ----------
     Given a genomic sequence and explicit intron ranges, remove them to make
     the mRNA/cDNA. Also a GT..AG scanner to suggest canonical introns.
  ------------------------------------------------------------------ */
  E.spliceIntrons = function (seq, introns) {
    seq = E.cleanSeq(seq);
    const sorted = introns.slice().sort((a, b) => b.start - a.start);
    let out = seq;
    for (const it of sorted) out = out.substring(0, it.start) + out.substring(it.end);
    return out;
  };

  /* ---------- Parsers ---------- */
  // every text reader normalises Windows (\r\n) and old Mac (\r) line endings
  function nl(text) { return String(text == null ? '' : text).replace(/\r\n?/g, '\n'); }
  E.normalizeNewlines = nl;

  /** @param {string} text @returns {SeqRecord[]} */
  E.parseFASTA = function (text) {
    text = nl(text);
    const recs = [];
    const blocks = text.split(/^>/m).filter(b => b.trim());
    for (const b of blocks) {
      const i = b.indexOf('\n');
      const header = (i === -1 ? b : b.slice(0, i)).trim();
      const body = i === -1 ? '' : b.slice(i + 1);
      const rep = E.seqTextReport(body);
      recs.push(/** @type {SeqRecord} */ ({ name: header.split(/\s+/)[0] || 'seq', description: header,
                  /* a protein record keeps its residues; see seqTextReport */
                  seq: rep.looksProtein ? rep.protein : rep.seq, features: [], topology: 'linear',
                  molecule: rep.looksProtein ? 'protein' : 'DNA',
                  dropped: rep.looksProtein ? 0 : rep.dropped,
                  droppedChars: rep.looksProtein ? '' : rep.droppedChars, hadU: rep.hadU }));
    }
    return recs;
  };

  /* ---------- Feature geometry shared by every view and exporter ----------
     A feature is {start, end, strand, parts?, between?} in 0-based half-open
     coordinates. end <= start means it wraps the origin of a circular
     sequence (start == end with no parts covers the whole circle). `parts`
     holds the pieces of a joined feature in top-strand order. `between`
     marks a zero-width site between two bases (GenBank a^b).            */
  /* A feature that crosses the origin can be written two ways: with its end
     before its start (4057..46), or with its end past the end of the
     sequence (4057..4888 on a 4842 bp circle). Both are produced by this
     program's own code and by the files it reads, and both mean the same
     thing. Only the first used to be split here; the second was handed back
     whole, and slice() then quietly returned the 785 bases up to the end of
     the string instead of the 831 the feature covers. Nothing raised an
     error — the length and the location still printed correctly, because
     those are arithmetic on the bounds — so the truncation only showed up
     as a feature that silently failed to carry across into an assembly
     product, where the bases of the copy no longer matched the original. */
  /**
   * The pieces of a feature as plain ranges on the top strand, an origin
   * crossing split in two.
   * @param {Feature} f @param {number} L
   * @returns {Array<{start:number,end:number}>}
   */
  E.featurePieces = function (f, L) {
    if (!f || f.between) return [];
    const src = (f.parts && f.parts.length) ? f.parts : [{ start: f.start, end: f.end }];
    const out = [];
    for (const p of src) {
      let a = p.start, b = p.end;
      if (L) {
        if (a >= L) { a -= L; b -= L; }                 // written entirely past the end
        if (b > L) {                                     // runs off the end: split at the origin
          out.push({ start: a, end: L });
          const tail = Math.min(b - L, a);               // never past its own start
          if (tail > 0) out.push({ start: 0, end: tail });
          continue;
        }
      }
      if (b > a) out.push({ start: a, end: b });
      else if (L) {
        out.push({ start: a, end: L });
        if (b > 0) out.push({ start: 0, end: b });
      }
    }
    return out;
  };
  E.featureLength = function (f, L) {
    return E.featurePieces(f, L).reduce((a, p) => a + (p.end - p.start), 0);
  };
  /** @param {string} seq @param {Feature} f @returns {string} */
  E.featureSequence = function (seq, f) {
    if (f && f.mixed && f.parts) {
      /* each part on its own strand, and each able to cross the origin */
      const L = seq.length;
      return f.parts.map(p => {
        let a = p.start, b = p.end;
        if (L && a >= L) { a -= L; b -= L; }
        const piece = (L && b > L) ? seq.slice(a) + seq.slice(0, Math.min(b - L, a))
                                   : seq.slice(a, b);
        return p.strand === -1 ? E.revcomp(piece) : piece;
      }).join('');
    }
    const s = E.featurePieces(f, seq.length).map(p => seq.slice(p.start, p.end)).join('');
    return f.strand === -1 ? E.revcomp(s) : s;
  };

  /* ---------- GenBank location grammar ----------
     complement(), join(), order(), bond(), one-of(), fuzzy ends (<, >),
     between-base sites (a^b), and ranges a..b. Remote references
     (ACC:1..10) cannot be resolved locally and drop the feature, as does any
     coordinate outside the sequence. A range written backwards (a..b, a>b)
     on a circular record is read as wrapping the origin, which is how
     SnapGene and Addgene often write such features.                    */
  function splitTop(str) {
    const out = []; let depth = 0, cur = '';
    for (const ch of str) {
      if (ch === '(') depth++;
      if (ch === ')') depth--;
      if (ch === ',' && depth === 0) { out.push(cur); cur = ''; } else cur += ch;
    }
    if (cur) out.push(cur);
    return out;
  }
  function parseLocExpr(expr, comp) {
    expr = expr.trim();
    let m;
    if ((m = /^complement\((.*)\)$/i.exec(expr))) {
      const inner = parseLocExpr(m[1], !comp);
      return inner ? inner.reverse() : null;
    }
    if ((m = /^(join|order|bond)\((.*)\)$/i.exec(expr))) {
      const out = [];
      for (const part of splitTop(m[2])) {
        const r = parseLocExpr(part, comp);
        if (!r) return null;
        out.push.apply(out, r);
      }
      return out;
    }
    if (/:/.test(expr)) return null;   // remote reference
    const pos = t => {
      t = t.trim();
      let om = /^one-of\(([^)]*)\)$/i.exec(t);
      if (om) t = om[1].split(',')[0];
      t = t.replace(/[<>]/g, '');
      if (/^\d+$/.test(t)) return +t;
      return NaN;
    };
    if ((m = /^(.+?)\^(.+)$/.exec(expr))) {
      const a = pos(m[1]);
      if (!isFinite(a)) return null;
      return [{ a: a, b: a, between: true, strand: comp ? -1 : 1 }];
    }
    if ((m = /^(.+?)\.\.(.+)$/.exec(expr))) {
      const a = pos(m[1]), b = pos(m[2]);
      if (!isFinite(a) || !isFinite(b)) return null;
      return [{ a: a, b: b, strand: comp ? -1 : 1 }];
    }
    if ((m = /^(\S+?)\.(\d+)$/.exec(expr)) && isFinite(pos(m[1]))) {   // old "a.b": one base in a..b
      const a = pos(m[1]);
      return [{ a: a, b: a, strand: comp ? -1 : 1 }];
    }
    const single = pos(expr);
    if (isFinite(single)) return [{ a: single, b: single, strand: comp ? -1 : 1 }];
    return null;
  }
  E.parseLocation = function (loc, L, circular) {
    const raw = parseLocExpr(String(loc).replace(/\s+/g, ''), false);
    if (!raw || !raw.length) return null;
    const mixedStrand = raw.some(r => r.strand === -1) && raw.some(r => r.strand !== -1);
    if (mixedStrand) {
      // e.g. trans-spliced genes: join(complement(a..b),c..d). Keep every piece
      // with its own strand, in biological order.
      const parts = [];
      for (const r of raw) {
        if (r.between) continue;
        if (r.a < 1 || r.b < 1 || (L > 0 && (r.a > L || r.b > L)) || r.a > r.b) return null;
        parts.push({ start: r.a - 1, end: r.b, strand: r.strand === -1 ? -1 : 1 });
      }
      if (!parts.length) return null;
      return { start: Math.min.apply(null, parts.map(p => p.start)),
               end: Math.max.apply(null, parts.map(p => p.end)),
               strand: parts[0].strand, parts: parts, mixed: true };
    }
    const strand = raw.every(r => r.strand === -1) ? -1 : 1;
    if (raw.length === 1 && raw[0].between) {
      const a = raw[0].a;
      if (a < 0 || (L > 0 && a > L)) return null;
      return { start: a, end: a, strand: raw[0].strand, between: true };
    }
    // pieces in biological order -> 0-based half-open, splitting origin wraps
    const bio = [];
    for (const r of raw) {
      if (r.between) continue;
      /* L is 0 when the record carries no ORIGIN at all — a CONTIG entry, or
         a features-only file. There is then nothing to be out of range of,
         and dropping every feature would throw away the whole annotation
         set, which is the only thing such a file holds. */
      if (r.a < 1 || r.b < 1 || (L > 0 && (r.a > L || r.b > L))) return null;
      if (r.a <= r.b) bio.push({ start: r.a - 1, end: r.b });
      else if (circular && L > 0) {              // written backwards: wraps the origin
        if (r.strand === -1) { bio.push({ start: 0, end: r.b }); bio.push({ start: r.a - 1, end: L }); }
        else { bio.push({ start: r.a - 1, end: L }); bio.push({ start: 0, end: r.b }); }
      } else if (L === 0) {
        /* no bases to wrap around: keep the two ends as given, so the
           annotation can still be applied to a sequence from elsewhere */
        bio.push({ start: r.a - 1, end: r.a });
        bio.push({ start: 0, end: r.b });
      } else return null;
    }
    if (!bio.length) return null;
    // top-strand order: minus-strand pieces are listed 3'->5' on the top strand
    const top = strand === -1 ? bio.slice().reverse() : bio;
    // merge an origin split that the file wrote as two adjacent pieces
    if (top.length === 2 && circular && top[0].end === L && top[1].start === 0) {
      return { start: top[0].start, end: top[1].end, strand: strand };
    }
    if (top.length === 1) return { start: top[0].start, end: top[0].end, strand: strand };
    let wraps = false;
    for (let i = 1; i < top.length; i++) if (top[i].start < top[i - 1].start) { wraps = true; break; }
    const start = wraps ? top[0].start : Math.min.apply(null, top.map(p => p.start));
    const end = wraps ? top[top.length - 1].end : Math.max.apply(null, top.map(p => p.end));
    return { start: start, end: end, strand: strand, parts: top };
  };

  E.formatLocation = function (f, L) {
    if (f.between) {
      const t = f.start + '^' + (f.start + 1);
      return f.strand === -1 ? 'complement(' + t + ')' : t;
    }
    const txt = p => (p.end - p.start === 1) ? String(p.start + 1) : (p.start + 1) + '..' + p.end;
    if (f.mixed && f.parts) {
      return 'join(' + f.parts.map(p => p.strand === -1 ? 'complement(' + txt(p) + ')' : txt(p)).join(',') + ')';
    }
    const pieces = E.featurePieces(f, L);
    if (!pieces.length) return null;
    let inner = pieces.length === 1 ? txt(pieces[0]) : 'join(' + pieces.map(txt).join(',') + ')';
    return f.strand === -1 ? 'complement(' + inner + ')' : inner;
  };

  /* GenBank is a column format, and a tab in it destroys the columns the
     parser reads. Editors and home-made export scripts put them there all
     the same, so they are expanded to 8-column tab stops first — the width
     every terminal and editor assumes — which puts the feature keys and
     qualifiers back where the format says they are. */
  function untab(s) {
    if (s.indexOf('\t') === -1) return s;
    return s.split('\n').map(line => {
      let out = '';
      for (const ch of line) {
        if (ch === '\t') { do { out += ' '; } while (out.length % 8); }
        else out += ch;
      }
      return out;
    }).join('\n');
  }

  /** @param {string} text @returns {SeqRecord[]} */
  E.parseGenBank = function (text) {
    text = untab(nl(text));
    const records = [];
    const chunks = text.split(/^\/\/[ \t]*$/m).filter(c => /\S/.test(c) && /^LOCUS/mi.test(c));
    for (const chunk of chunks) {
      const rec = /** @type {SeqRecord} */ ({ name: 'Untitled', description: '', seq: '', features: [], topology: 'linear', molecule: 'DNA' });
      const locusLine = (chunk.match(/^LOCUS[^\n]*/mi) || [''])[0];
      const locus = locusLine.match(/^LOCUS\s+(\S+)/i);
      if (locus) rec.name = locus[1];
      if (/\bcircular\b/i.test(locusLine)) rec.topology = 'circular';
      /* "5467 bp" and "412 aa" both appear in this column; ApE writes
         "ds-DNA" where NCBI writes "DNA", and neither says protein */
      if (/\s\d+\s+aa\b/i.test(locusLine)) rec.molecule = 'protein';
      const defm = chunk.match(/^DEFINITION {2}([^\n]*(?:\n {12}[^\n]*)*)/mi);
      if (defm) rec.description = defm[1].replace(/\s+/g, ' ').trim();
      const originAt = chunk.search(/^ORIGIN/mi);
      if (originAt !== -1) {
        const rep = E.seqTextReport(chunk.substring(originAt).replace(/^ORIGIN[^\n]*/mi, ''));
        if (rep.looksProtein) rec.molecule = 'protein';
        /* LOCUS said "aa", or the letters themselves say so: keep the
           residues rather than the DNA-cleaned wreck of them */
        const prot = rec.molecule === 'protein';
        rec.seq = prot ? rep.protein : rep.seq;
        rec.dropped = prot ? 0 : rep.dropped;
        rec.droppedChars = prot ? '' : rep.droppedChars;
      } else if (/^CONTIG/mi.test(chunk)) {
        /* a CONTIG record points at pieces held elsewhere and carries no
           bases of its own; saying so beats importing an empty construct */
        rec.contigOnly = true;
      }
      const featStart = chunk.search(/^FEATURES/mi);
      if (featStart !== -1) {
        const after = chunk.substring(featStart);
        // the block runs until the next line that starts in column 1
        const stop = after.search(/^(?!FEATURES)\S/m);
        const block = stop > 0 ? after.substring(0, stop) : after;
        /* Addgene's exports, and SnapGene's, write "linear" on the LOCUS line
           for plasmids and then annotate a feature that runs off the end and
           back to base 1 — pSV40's enhancer at 8569..276 on an 8,596 bp
           circle, for instance. Taken at its word the molecule is linear, the
           location is impossible, and the feature is dropped: on the common
           case of an Addgene plasmid that silently loses an annotation. A
           location written backwards cannot exist on a linear molecule, so it
           is itself the better evidence of the two, and the correction is
           recorded rather than made quietly. */
        if (rec.topology !== 'circular' && rec.seq.length > 0 && backwardsLocation(block, rec.seq.length)) {
          rec.topology = 'circular';
          rec.topologyFromFeature = true;
        }
        rec.features = parseGBFeatures(block, rec.seq.length, rec.topology === 'circular');
        /* the refused locations travel on the array; they are moved onto
           the record and taken off it, so a caller sees one list of
           features and one list of what was declined */
        const fl = /** @type {FeatureList} */ (rec.features);
        if (fl.skipped) { rec.skippedFeatures = fl.skipped; delete fl.skipped; }
      }
      records.push(rec);
    }
    return records;
  };

  function unquote(v, key) {
    if (v === true) return true;
    if (v.length >= 2 && v[0] === '"' && v[v.length - 1] === '"') v = v.slice(1, -1).replace(/""/g, '"');
    if (key === 'translation') v = v.replace(/\s+/g, '');
    return v;
  }
  /* The feature table is laid out by indent: the key in column 6, the
     qualifiers in column 22. Files that came through an editor, a script,
     or a tool with its own idea of the format use other depths, so the
     depth a feature starts at is measured from the block rather than
     assumed to be 5. Anything deeper than that is a continuation; anything
     at column 1 is the next GenBank keyword and ends the block. */
  function featureIndent(lines) {
    let min = Infinity;
    for (const line of lines) {
      if (!/\S/.test(line)) continue;
      const ind = line.length - line.replace(/^ +/, '').length;
      if (ind === 0) continue;
      const rest = line.slice(ind);
      if (rest[0] === '/') continue;         /* a qualifier is never the shallowest */
      if (ind < min) min = ind;
    }
    return min === Infinity ? 5 : min;
  }

  /* Is there a plain "a..b" with a > b, both inside the sequence, in this
     feature block? That is a span running off the end and round to base 1,
     which only a circular molecule can hold. Joins and cross-file
     references are left out of it: a join written out of order says nothing
     about topology, and an accession-qualified part is not on this sequence
     at all. */
  function backwardsLocation(block, L) {
    const re = /^ {2,}\S+ +(\d+)\.\.(\d+) *$/gm;
    let m;
    while ((m = re.exec(block))) {
      const a = +m[1], b = +m[2];
      if (a > b && a <= L && b >= 1) return true;
    }
    return false;
  }

  /**
   * @param {string} block @param {number} L sequence length, 0 when there is none
   * @param {boolean} circular
   * @returns {FeatureList}
   */
  function parseGBFeatures(block, L, circular) {
    const lines = block.split('\n').slice(1);
    const base = featureIndent(lines);
    const raw = [];
    let cur = null;
    for (const line of lines) {
      if (!/\S/.test(line)) continue;
      const ind = line.length - line.replace(/^ +/, '').length;
      if (ind === 0) continue;
      if (ind === base) {
        const m = line.slice(ind).match(/^(\S+)\s*(.*)$/);
        if (!m) continue;
        cur = { type: m[1], loc: m[2].trim(), qlines: [], inQual: false };
        raw.push(cur);
      } else if (cur && ind > base) {
        const t = line.replace(/^\s+/, '').replace(/\s+$/, '');
        if (!cur.inQual && t[0] !== '/') cur.loc += t;
        else { cur.inQual = true; cur.qlines.push(t); }
      }
    }
    const feats = /** @type {FeatureList} */ ([]);
    for (const r of raw) {
      // qualifiers: a new one starts with "/" only when no quote is open
      const quals = {};
      let key = null, val = null;
      const open = v => typeof v === 'string' && v[0] === '"' && ((v.match(/"/g) || []).length % 2 === 1);
      const flush = () => {
        if (key == null) return;
        const u = unquote(val, key);
        if (quals[key] === undefined) quals[key] = u;
        key = null; val = null;
      };
      for (const t of r.qlines) {
        if (t[0] === '/' && !open(val)) {
          flush();
          const m = /^\/([^=]+)(?:=(.*))?$/.exec(t);
          if (!m) continue;
          key = m[1]; val = m[2] != null ? m[2] : true;
        } else if (key != null && typeof val === 'string') {
          val += (key === 'translation' ? '' : ' ') + t;
        }
      }
      flush();
      const loc = E.parseLocation(r.loc, L, circular);
      if (!loc) { feats.skipped = (feats.skipped || []).concat([r.type + ' ' + r.loc]); continue; }
      const pick = k => (typeof quals[k] === 'string' && quals[k]) ? quals[k] : null;
      const f = {
        type: r.type, start: loc.start, end: loc.end, strand: loc.strand,
        label: pick('label') || pick('gene') || pick('product') || pick('locus_tag') ||
               pick('standard_name') || pick('note') || r.type,
        quals: quals
      };
      /* the colour the file was drawn in, where it carries one: ApE writes
         /ApEinfo_fwdcolor, SnapGene's GenBank export writes it inside the
         note as "color: #rrggbb", and some tools use /color or /colour.
         Keeping it means a map imported from elsewhere still looks like the
         map the person knows. */
      const col = pick('ApEinfo_fwdcolor') || pick('ApEinfo_revcolor') || pick('color') || pick('colour') ||
        ((pick('note') || '').match(/colou?r:\s*(#[0-9a-f]{3,8})/i) || [])[1];
      if (col && /^#[0-9a-f]{3}(?:[0-9a-f]{3}(?:[0-9a-f]{2})?)?$/i.test(col.trim())) f.color = col.trim();
      if (loc.parts) f.parts = loc.parts;
      if (loc.between) f.between = true;
      if (loc.mixed) f.mixed = true;
      /* The label is kept whole. It used to be cut to 80 characters with
         an ellipsis here, in the reader, which is a drawing decision made
         in the wrong place: the shortened string then went into the feature
         table, into exports, and back out to GenBank, so a long /note used
         as a label — the fallback when there is no /label, /gene or
         /product — was permanently truncated by the act of reading the
         file. Shortening for display happens where the label is drawn. */
      feats.push(f);
    }
    return feats;
  }

  /* EMBL flat files share GenBank's feature-table columns behind an "FT"
     prefix, so the same verified location and qualifier parser is reused. */
  /** @param {string} text @returns {SeqRecord[]} */
  E.parseEMBL = function (text) {
    text = nl(text);
    const records = [];
    text.split(/^\/\/[ \t]*$/m).forEach(chunk => {
      if (!/^ID {3}/m.test(chunk)) return;
      const idLine = (chunk.match(/^ID {3}[^\n]*/m) || [''])[0];
      const rec = /** @type {SeqRecord} */ ({ name: 'EMBL record', description: '', seq: '', features: [], topology: 'linear', molecule: 'DNA' });
      const idm = /^ID {3}([^;\s]+)/.exec(idLine);
      if (idm) rec.name = idm[1];
      if (/;\s*circular\s*;/i.test(idLine)) rec.topology = 'circular';
      if (/;\s*(protein|PRT)\s*;/i.test(idLine) || /\bAA\.?\s*$/.test(idLine)) rec.molecule = 'protein';
      const acc = /^AC {3}([^;\s]+)/m.exec(chunk);
      if (acc && (!idm || /^XXX/.test(rec.name))) rec.name = acc[1];
      const de = chunk.match(/^DE {3}.*(?:\nDE {3}.*)*/m);
      if (de) rec.description = de[0].replace(/^DE {3}/gm, '').replace(/\s+/g, ' ').trim();
      const sq = chunk.search(/^SQ {3}/m);
      if (sq !== -1) {
        const rep = E.seqTextReport(chunk.substring(sq).replace(/^SQ[^\n]*/m, ''));
        if (rep.looksProtein) rec.molecule = 'protein';
        const prot = rec.molecule === 'protein';
        rec.seq = prot ? rep.protein : rep.seq;
        rec.dropped = prot ? 0 : rep.dropped;
        rec.droppedChars = prot ? '' : rep.droppedChars;
      }
      const ft = (chunk.match(/^FT[^\n]*/gm) || []).map(l => '  ' + l.slice(2));
      if (ft.length) {
        rec.features = parseGBFeatures('FEATURES\n' + ft.join('\n'), rec.seq.length, rec.topology === 'circular');
        /* the refused locations travel on the array; they are moved onto
           the record and taken off it, so a caller sees one list of
           features and one list of what was declined */
        const fl = /** @type {FeatureList} */ (rec.features);
        if (fl.skipped) { rec.skippedFeatures = fl.skipped; delete fl.skipped; }
      }
      records.push(rec);
    });
    return records;
  };

  E.toFASTA = function (name, seq, width) {
    width = width || 70;
    let out = '>' + name + '\n';
    for (let i = 0; i < seq.length; i += width) out += seq.substr(i, width) + '\n';
    return out;
  };

  function wrapQual(key, value) {
    // /key="value" wrapped at column 79 with a 21-space hanging indent
    const ind = ' '.repeat(21);
    const text = '/' + key + '="' + String(value).replace(/"/g, '""') + '"';
    const lines = [];
    let rest = text;
    while (rest.length > 58) {
      /* Only at a space. A reader has no way to tell a line break that
         replaced a space from one that fell inside a word, so it rejoins
         with a space either way — which means breaking inside a word puts
         a space into the text. A /note holding a URL or a long identifier
         (ApE and SnapGene both write them) came back altered. GenBank
         allows the line to run long instead, which is what Biopython's
         writer does, so a word with no space in it is left whole. */
      const cut = rest.lastIndexOf(' ', 58);
      if (cut < 20) {
        const next = rest.indexOf(' ', 58);
        if (next === -1) break;                 // nothing left to break at
        lines.push(ind + rest.slice(0, next));
        rest = rest.slice(next).replace(/^\s+/, '');
        continue;
      }
      lines.push(ind + rest.slice(0, cut).replace(/\s+$/, ''));
      rest = rest.slice(cut).replace(/^\s+/, '');
    }
    lines.push(ind + rest);
    return lines.join('\n') + '\n';
  }
  function wrapLocation(loc) {
    if (loc.length <= 58) return loc;
    const out = [];
    let rest = loc;
    while (rest.length > 58) {
      let cut = rest.lastIndexOf(',', 58);
      if (cut < 1) cut = 58; else cut += 1;
      out.push(rest.slice(0, cut));
      rest = rest.slice(cut);
    }
    out.push(rest);
    return out.join('\n' + ' '.repeat(21));
  }
  /** @param {SeqRecord} rec @returns {string} */
  E.toGenBank = function (rec) {
    const L = rec.seq.length;
    const d = new Date();
    const MON = ['JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC'];
    const date = String(d.getDate()).padStart(2, '0') + '-' + MON[d.getMonth()] + '-' + d.getFullYear();
    const nm = String(rec.name || 'construct').replace(/\s+/g, '_').slice(0, 16);
    let out = 'LOCUS       ' + pad(nm, 16) + ' ' + String(L).padStart(11) + ' bp    DNA     ' +
      pad(rec.topology === 'circular' ? 'circular' : 'linear', 8) + ' SYN ' + date + '\n';
    out += 'DEFINITION  ' + (rec.description ? String(rec.description).replace(/\s+/g, ' ').slice(0, 900) : '.') + '\n';
    out += 'ACCESSION   .\nVERSION     .\nKEYWORDS    .\nSOURCE      synthetic DNA construct\n' +
           '  ORGANISM  synthetic DNA construct\n';
    out += 'FEATURES             Location/Qualifiers\n';
    out += '     source          1..' + L + '\n' + wrapQual('mol_type', 'other DNA');
    for (const f of (rec.features || [])) {
      if (f.type === 'source') continue;           // written once, above
      const loc = E.formatLocation(f, L);
      if (!loc) continue;
      const type = String(f.type || 'misc_feature').replace(/\s+/g, '_').slice(0, 15);
      out += '     ' + pad(type, 15) + ' ' + wrapLocation(loc) + '\n';
      if (f.label) out += wrapQual('label', f.label);
      const q = f.quals || {};
      Object.keys(q).forEach(k => {
        if (k === 'label' || k === 'translation' || k[0] === '_') return;
        const v = q[k];
        if (v === true) out += ' '.repeat(21) + '/' + k + '\n';
        else if (typeof v === 'string' && v.length) out += wrapQual(k, v);
      });
    }
    out += 'ORIGIN\n';
    const sq = rec.seq.toLowerCase();
    for (let i = 0; i < sq.length; i += 60) {
      out += String(i + 1).padStart(9) + ' ' + (sq.substr(i, 60).match(/.{1,10}/g) || []).join(' ') + '\n';
    }
    out += '//\n';
    return out;
  };
  function pad(s, n) { s = s || ''; return (s + ' '.repeat(n)).slice(0, n); }

  /* ---------- Common feature library (for auto-annotation) ----------
     Short signature sequences for frequent plasmid elements. Matching is by
     substring (both strands). This is a curated starter set.
  ------------------------------------------------------------------ */
  E.FEATURE_LIB = [
    { name:'T7 promoter', type:'promoter', seq:'TAATACGACTCACTATAG' },
    { name:'T7 terminator', type:'terminator', seq:'CTAGCATAACCCCTTGGGGCCTCTAAACGGGTCTTGAGGGGTTTTTTG' },
    { name:'lac promoter', type:'promoter', seq:'TTTACACTTTATGCTTCCGGCTCGTATGTTG' },
    { name:'CMV promoter', type:'promoter', seq:'GTTGACATTGATTATTGACTAGTTATTAATAGTAATCAATTACGGGGTCATTAGTTCAT' },
    { name:'SV40 polyA', type:'polyA_signal', seq:'AACTTGTTTATTGCAGCTTATAATGGTTACAAATAAAGCAATAGCATCACAAATTTCACAAATAAAGCATTTTTTTCACTGC' },
    { name:'AmpR (bla) CDS start', type:'CDS', seq:'ATGAGTATTCAACATTTCCGTGTCGCCCTTATTCCCTTTTTTGCGGCATTTTGCCTTCCTGTTTTTGCTCACCCAGAAACGCTGGTG' },
    { name:'KanR start', type:'CDS', seq:'ATGAGCCATATTCAACGGGAAACGTCTTGCTCGAGGCCGCGATTAAATTCCAACATGG' },
    { name:'pUC ori', type:'rep_origin', seq:'TTGAGATCCTTTTTTTCTGCGCGTAATCTGCTGCTTGCAAACAAAAAAACCACCGCTACCAGCGG' },
    { name:'f1 ori', type:'rep_origin', seq:'ACGCGCCCTGTAGCGGCGCATTAAGCGCGGCGGGTGTGGTGGTTACGCGCAGCGTGACCGCTACACTTGCCAGCG' },
    { name:'His6 tag', type:'CDS', seq:'CATCATCATCATCATCAT' },
    { name:'FLAG tag', type:'CDS', seq:'GATTACAAGGATGACGACGATAAG' },
    { name:'HA tag', type:'CDS', seq:'TACCCATACGATGTTCCAGATTACGCT' },
    { name:'Myc tag', type:'CDS', seq:'GAACAAAAACTCATCTCAGAAGAGGATCTG' },
    { name:'EGFP start', type:'CDS', seq:'ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGG' },
    { name:'mCherry start', type:'CDS', seq:'ATGGTGAGCAAGGGCGAGGAGGATAACATGGCC' },
    { name:'SV40 NLS', type:'CDS', seq:'CCAAAGAAGAAGCGGAAGGTC' },
    { name:'P2A', type:'CDS', seq:'GGAAGCGGAGCTACTAACTTCAGCCTGCTGAAGCAGGCTGGAGACGTGGAGGAGAACCCTGGACCT' },
    { name:'WPRE', type:'regulatory', seq:'AATCAACCTCTGGATTACAAAATTTGTGAAAGATTGACTGGTATTCTTAACTATGTTGCTCCTTTTACGCTATGTGGATACGCTGCTTTAATGCCTTTGTATCATGCTATTGCTTCCCGTATGGCTTTCATTTTCTCCTCCTTGTATAAATCCTGGTTGCTGTCTCTTTATGAGG' },
    { name:'lac operator', type:'protein_bind', seq:'AATTGTGAGCGGATAACAATT' },
    { name:'RBS (consensus)', type:'RBS', seq:'AGGAGG', needsStartCodon: [4, 14] },
    { name:'SP6 promoter', type:'promoter', seq:'ATTTAGGTGACACTATAG' },
    { name:'T3 promoter', type:'promoter', seq:'AATTAACCCTCACTAAAGG' },
    { name:'BGH polyA', type:'polyA_signal', seq:'CTGTGCCTTCTAGTTGCCAGCCATCTGTTGTTTGCCCCTCCCCCGTGCCTTCCTTGACCCTGGAAGGTGCCACTCCCACTGTCCTTTCCTAAT' },
    { name:'attB1', type:'recomb', seq:'ACAAGTTTGTACAAAAAAGCAGGCT' },
    { name:'attB2', type:'recomb', seq:'ACCCAGCTTTCTTGTACAAAGTGG' },
    { name:'loxP', type:'recomb', seq:'ATAACTTCGTATAATGTATGCTATACGAAGTTAT' },
    { name:'FRT', type:'recomb', seq:'GAAGTTCCTATTCTCTAGAAAGTATAGGAACTTC' }
  ];
  /* opts.circular lets an element that crosses base 1 be found: the search
     runs over the sequence plus enough of its own start, and a hit that runs
     past the end comes back wrapped (end <= start), as features are written
     everywhere else. A second argument that is not a library is read as opts,
     so detectFeatures(seq, {circular:true}) does the expected thing. */
  /* How often a signature of this length would turn up in a sequence of this
     length by chance alone, counting both strands and any ambiguity codes.
     A match whose expectation is not small is not evidence of anything. */
  E.chanceHits = function (sig, seqLen) {
    const AMB = { R:2, Y:2, S:2, W:2, K:2, M:2, B:3, D:3, H:3, V:3, N:4 };
    let fold = 1;
    for (const c of sig.toUpperCase()) fold *= (AMB[c] || 1);
    const places = Math.max(0, seqLen - sig.length + 1) * 2;
    return places * fold / Math.pow(4, sig.length);
  };

  E.detectFeatures = function (seq, lib, opts) {
    if (lib && !Array.isArray(lib)) { opts = lib; lib = null; }
    opts = opts || {};
    lib = lib || E.FEATURE_LIB;
    seq = E.cleanSeq(seq);
    const L = seq.length;
    const longest = lib.reduce((a, f) => Math.max(a, (f.seq || '').length), 0);
    const circular = !!opts.circular && L > longest;
    const hay = circular ? seq + seq.slice(0, Math.max(0, longest - 1)) : seq;
    const found = [], seen = new Set();
    const at = i => hay[((i % hay.length) + hay.length) % hay.length];
    /* A short motif only means something in context: a ribosome binding site
       has to sit a few bases in front of a start codon, so one in the middle
       of a gene is a coincidence, not a finding. */
    function hasStartCodon(idx, len, strand, win) {
      for (let d = win[0]; d <= win[1]; d++) {
        if (strand === 1) {
          const j = idx + len + d - 1;
          if (at(j) === 'A' && at(j + 1) === 'T' && at(j + 2) === 'G') return true;
        } else {
          const j = idx - d - 2;                 // reading leftwards on the bottom strand
          if (at(j) === 'C' && at(j + 1) === 'A' && at(j + 2) === 'T') return true;
        }
      }
      return false;
    }
    const add = (f, idx, len, strand) => {
      if (idx >= L) return;                       // a repeat of one already found
      if (f.needsStartCodon && !hasStartCodon(idx, len, strand, f.needsStartCodon)) return;
      const end = idx + len;
      const g = Object.assign({}, f, { start: idx, end: end > L ? end - L : end, strand: strand });
      if (end > L) g.wrapped = true;
      g.expected = E.chanceHits(E.cleanSeq(f.seq), L);
      g.confident = g.expected < 0.05;
      const k = g.start + ':' + g.end + ':' + strand + ':' + (f.name || '');
      if (seen.has(k)) return;
      seen.add(k); found.push(g);
    };
    for (const f of lib) {
      const sig = E.cleanSeq(f.seq);
      if (sig.length < 6) continue;
      let idx = hay.indexOf(sig);
      while (idx !== -1) { add(f, idx, sig.length, 1); idx = hay.indexOf(sig, idx + 1); }
      const rc = E.revcomp(sig);
      idx = hay.indexOf(rc);
      while (idx !== -1) { add(f, idx, rc.length, -1); idx = hay.indexOf(rc, idx + 1); }
    }
    return found.sort((a, b) => a.start - b.start);
  };

  /* ---------- E. coli codon usage (relative fractions within AA) ----------
     From the Kazusa E. coli K-12 table (fraction per AA). Used by codon opt.
  ------------------------------------------------------------------ */
  E.ECOLI_USAGE = {
    TTT:0.57,TTC:0.43,TTA:0.13,TTG:0.13,CTT:0.10,CTC:0.10,CTA:0.04,CTG:0.50,
    ATT:0.51,ATC:0.42,ATA:0.07,ATG:1.00,GTT:0.26,GTC:0.22,GTA:0.15,GTG:0.37,
    TCT:0.15,TCC:0.15,TCA:0.12,TCG:0.15,AGT:0.15,AGC:0.28,CCT:0.16,CCC:0.12,
    CCA:0.19,CCG:0.53,ACT:0.17,ACC:0.44,ACA:0.13,ACG:0.27,GCT:0.16,GCC:0.27,
    GCA:0.21,GCG:0.36,TAT:0.57,TAC:0.43,TAA:0.64,TAG:0.07,TGA:0.29,CAT:0.57,
    CAC:0.43,CAA:0.34,CAG:0.66,AAT:0.49,AAC:0.51,AAA:0.74,AAG:0.26,GAT:0.63,
    GAC:0.37,GAA:0.68,GAG:0.32,TGT:0.46,TGC:0.54,TGG:1.00,CGT:0.36,CGC:0.36,
    CGA:0.07,CGG:0.11,AGA:0.07,AGG:0.04,GGT:0.35,GGC:0.37,GGA:0.13,GGG:0.15
  };
  E.HUMAN_USAGE = {
    TTT:0.46,TTC:0.54,TTA:0.08,TTG:0.13,CTT:0.13,CTC:0.20,CTA:0.07,CTG:0.40,
    ATT:0.36,ATC:0.47,ATA:0.17,ATG:1.00,GTT:0.18,GTC:0.24,GTA:0.12,GTG:0.46,
    TCT:0.19,TCC:0.22,TCA:0.15,TCG:0.05,AGT:0.15,AGC:0.24,CCT:0.29,CCC:0.32,
    CCA:0.28,CCG:0.11,ACT:0.25,ACC:0.36,ACA:0.28,ACG:0.11,GCT:0.27,GCC:0.40,
    GCA:0.23,GCG:0.11,TAT:0.44,TAC:0.56,TAA:0.30,TAG:0.24,TGA:0.47,CAT:0.42,
    CAC:0.58,CAA:0.27,CAG:0.73,AAT:0.47,AAC:0.53,AAA:0.43,AAG:0.57,GAT:0.46,
    GAC:0.54,GAA:0.42,GAG:0.58,TGT:0.46,TGC:0.54,TGG:1.00,CGT:0.08,CGC:0.18,
    CGA:0.11,CGG:0.20,AGA:0.21,AGG:0.21,GGT:0.16,GGC:0.34,GGA:0.25,GGG:0.25
  };

  if (typeof module !== 'undefined' && module.exports) module.exports = E;
  root.Bio = E;
})(typeof window !== 'undefined' ? window : globalThis);
