
/* ===========================================================================
   EXPRESSION STRATEGY — will this thing be soluble, and in what?

   Two separate jobs live here, and they are kept separate on purpose.

   The first is reading a protein and saying how hard it is likely to be.
   That is done with published, checkable descriptors - hydropathy, charge,
   the instability index, predicted membrane helices and so on - each
   reported with its own value and what it means, never rolled into a
   single black-box number. Every one of them is cross-checked against
   Biopython to the last decimal, so what the tool says about a sequence is
   at least arithmetically beyond doubt. Whether a protein folds is not
   something any composition-based method can settle, and the wording
   throughout says so.

   The second is knowing what a lab actually has on its shelf. A vector
   name like "pACYCDuet-1" carries a replication origin, a copy number, a
   selection marker, a promoter and an induction route, and whether two of
   them can live in the same cell is a hard constraint, not a preference.
   That part is knowledge, not prediction, and it is where the tool can be
   held to account: Cm chaperones on a p15A backbone cannot go in with a
   p15A vector, a Rosetta host already carries one, and two ColE1 plasmids
   will not be kept together whatever the markers say.
   ========================================================================= */
(function (root) {
  'use strict';
  const E = root.Bio;
  
  /* ---------------------------------------------------------------------
     PUBLISHED SCALES
     Kyte & Doolittle hydropathy (1982); Guruprasad dipeptide instability
     weights (1990); the pK set Biopython uses for the isoelectric point.
     These are data, not choices, and they are here verbatim so that the
     numbers agree with the reference implementation exactly.
  --------------------------------------------------------------------- */
  E.KD = {"A":1.8,"R":-4.5,"N":-3.5,"D":-3.5,"C":2.5,"Q":-3.5,"E":-3.5,"G":-0.4,"H":-3.2,"I":4.5,"L":3.8,"K":-3.9,"M":1.9,"F":2.8,"P":-1.6,"S":-0.8,"T":-0.7,"W":-0.9,"Y":-1.3,"V":4.2};
  /* Guruprasad dipeptide instability weights (1990) */
  E.DIWV = {"AA":1.0,"AC":44.94,"AD":-7.49,"AE":1.0,"AF":1.0,"AG":1.0,"AH":-7.49,"AI":1.0,"AK":1.0,"AL":1.0,"AM":1.0,"AN":1.0,"AP":20.26,"AQ":1.0,"AR":1.0,"AS":1.0,"AT":1.0,"AV":1.0,"AW":1.0,"AY":1.0,"CA":1.0,"CC":1.0,"CD":20.26,"CE":1.0,"CF":1.0,"CG":1.0,"CH":33.6,"CI":1.0,"CK":1.0,"CL":20.26,"CM":33.6,"CN":1.0,"CP":20.26,"CQ":-6.54,"CR":1.0,"CS":1.0,"CT":33.6,"CV":-6.54,"CW":24.68,"CY":1.0,"DA":1.0,"DC":1.0,"DD":1.0,"DE":1.0,"DF":-6.54,"DG":1.0,"DH":1.0,"DI":1.0,"DK":-7.49,"DL":1.0,"DM":1.0,"DN":1.0,"DP":1.0,"DQ":1.0,"DR":-6.54,"DS":20.26,"DT":-14.03,"DV":1.0,"DW":1.0,"DY":1.0,"EA":1.0,"EC":44.94,"ED":20.26,"EE":33.6,"EF":1.0,"EG":1.0,"EH":-6.54,"EI":20.26,"EK":1.0,"EL":1.0,"EM":1.0,"EN":1.0,"EP":20.26,"EQ":20.26,"ER":1.0,"ES":20.26,"ET":1.0,"EV":1.0,"EW":-14.03,"EY":1.0,"FA":1.0,"FC":1.0,"FD":13.34,"FE":1.0,"FF":1.0,"FG":1.0,"FH":1.0,"FI":1.0,"FK":-14.03,"FL":1.0,"FM":1.0,"FN":1.0,"FP":20.26,"FQ":1.0,"FR":1.0,"FS":1.0,"FT":1.0,"FV":1.0,"FW":1.0,"FY":33.601,"GA":-7.49,"GC":1.0,"GD":1.0,"GE":-6.54,"GF":1.0,"GG":13.34,"GH":1.0,"GI":-7.49,"GK":-7.49,"GL":1.0,"GM":1.0,"GN":-7.49,"GP":1.0,"GQ":1.0,"GR":1.0,"GS":1.0,"GT":-7.49,"GV":1.0,"GW":13.34,"GY":-7.49,"HA":1.0,"HC":1.0,"HD":1.0,"HE":1.0,"HF":-9.37,"HG":-9.37,"HH":1.0,"HI":44.94,"HK":24.68,"HL":1.0,"HM":1.0,"HN":24.68,"HP":-1.88,"HQ":1.0,"HR":1.0,"HS":1.0,"HT":-6.54,"HV":1.0,"HW":-1.88,"HY":44.94,"IA":1.0,"IC":1.0,"ID":1.0,"IE":44.94,"IF":1.0,"IG":1.0,"IH":13.34,"II":1.0,"IK":-7.49,"IL":20.26,"IM":1.0,"IN":1.0,"IP":-1.88,"IQ":1.0,"IR":1.0,"IS":1.0,"IT":1.0,"IV":-7.49,"IW":1.0,"IY":1.0,"KA":1.0,"KC":1.0,"KD":1.0,"KE":1.0,"KF":1.0,"KG":-7.49,"KH":1.0,"KI":-7.49,"KK":1.0,"KL":-7.49,"KM":33.6,"KN":1.0,"KP":-6.54,"KQ":24.64,"KR":33.6,"KS":1.0,"KT":1.0,"KV":-7.49,"KW":1.0,"KY":1.0,"LA":1.0,"LC":1.0,"LD":1.0,"LE":1.0,"LF":1.0,"LG":1.0,"LH":1.0,"LI":1.0,"LK":-7.49,"LL":1.0,"LM":1.0,"LN":1.0,"LP":20.26,"LQ":33.6,"LR":20.26,"LS":1.0,"LT":1.0,"LV":1.0,"LW":24.68,"LY":1.0,"MA":13.34,"MC":1.0,"MD":1.0,"ME":1.0,"MF":1.0,"MG":1.0,"MH":58.28,"MI":1.0,"MK":1.0,"ML":1.0,"MM":-1.88,"MN":1.0,"MP":44.94,"MQ":-6.54,"MR":-6.54,"MS":44.94,"MT":-1.88,"MV":1.0,"MW":1.0,"MY":24.68,"NA":1.0,"NC":-1.88,"ND":1.0,"NE":1.0,"NF":-14.03,"NG":-14.03,"NH":1.0,"NI":44.94,"NK":24.68,"NL":1.0,"NM":1.0,"NN":1.0,"NP":-1.88,"NQ":-6.54,"NR":1.0,"NS":1.0,"NT":-7.49,"NV":1.0,"NW":-9.37,"NY":1.0,"PA":20.26,"PC":-6.54,"PD":-6.54,"PE":18.38,"PF":20.26,"PG":1.0,"PH":1.0,"PI":1.0,"PK":1.0,"PL":1.0,"PM":-6.54,"PN":1.0,"PP":20.26,"PQ":20.26,"PR":-6.54,"PS":20.26,"PT":1.0,"PV":20.26,"PW":-1.88,"PY":1.0,"QA":1.0,"QC":-6.54,"QD":20.26,"QE":20.26,"QF":-6.54,"QG":1.0,"QH":1.0,"QI":1.0,"QK":1.0,"QL":1.0,"QM":1.0,"QN":1.0,"QP":20.26,"QQ":20.26,"QR":1.0,"QS":44.94,"QT":1.0,"QV":-6.54,"QW":1.0,"QY":-6.54,"RA":1.0,"RC":1.0,"RD":1.0,"RE":1.0,"RF":1.0,"RG":-7.49,"RH":20.26,"RI":1.0,"RK":1.0,"RL":1.0,"RM":1.0,"RN":13.34,"RP":20.26,"RQ":20.26,"RR":58.28,"RS":44.94,"RT":1.0,"RV":1.0,"RW":58.28,"RY":-6.54,"SA":1.0,"SC":33.6,"SD":1.0,"SE":20.26,"SF":1.0,"SG":1.0,"SH":1.0,"SI":1.0,"SK":1.0,"SL":1.0,"SM":1.0,"SN":1.0,"SP":44.94,"SQ":20.26,"SR":20.26,"SS":20.26,"ST":1.0,"SV":1.0,"SW":1.0,"SY":1.0,"TA":1.0,"TC":1.0,"TD":1.0,"TE":20.26,"TF":13.34,"TG":-7.49,"TH":1.0,"TI":1.0,"TK":1.0,"TL":1.0,"TM":1.0,"TN":-14.03,"TP":1.0,"TQ":-6.54,"TR":1.0,"TS":1.0,"TT":1.0,"TV":1.0,"TW":-14.03,"TY":1.0,"VA":1.0,"VC":1.0,"VD":-14.03,"VE":1.0,"VF":1.0,"VG":-7.49,"VH":1.0,"VI":1.0,"VK":-1.88,"VL":1.0,"VM":1.0,"VN":1.0,"VP":20.26,"VQ":1.0,"VR":1.0,"VS":1.0,"VT":-7.49,"VV":1.0,"VW":1.0,"VY":-6.54,"WA":-14.03,"WC":1.0,"WD":1.0,"WE":1.0,"WF":1.0,"WG":-9.37,"WH":24.68,"WI":1.0,"WK":1.0,"WL":13.34,"WM":24.68,"WN":13.34,"WP":1.0,"WQ":1.0,"WR":1.0,"WS":1.0,"WT":-14.03,"WV":-7.49,"WW":1.0,"WY":1.0,"YA":24.68,"YC":1.0,"YD":24.68,"YE":-6.54,"YF":1.0,"YG":-7.49,"YH":13.34,"YI":1.0,"YK":1.0,"YL":1.0,"YM":44.94,"YN":1.0,"YP":13.34,"YQ":1.0,"YR":-15.91,"YS":1.0,"YT":-7.49,"YV":1.0,"YW":-9.37,"YY":13.34};
  /* amino-acid average residue masses, and water, for the molecular weight */
  /* Monomer masses of the free amino acids. A chain weighs their sum less
     one water for every bond formed, which is how the reference computes
     it; taking a water off each residue first and adding one back rounds
     differently and drifts by a tenth of a dalton over a long protein. */
  E.AA_MONOMER = {"A":89.0932,"C":121.1582,"D":133.1027,"E":147.1293,"F":165.1891,
    "G":75.0666,"H":155.1546,"I":131.1729,"K":146.1876,"L":131.1729,"M":149.2113,
    "N":132.1179,"O":255.3134,"P":115.1305,"Q":146.1445,"R":174.201,"S":105.0926,
    "T":119.1192,"U":168.0532,"V":117.1463,"W":204.2252,"Y":181.1885};
  E.WATER = 18.0153;

  const PK_POS = { Nterm: 7.5, K: 10.0, R: 12.0, H: 5.98 };
  const PK_NEG = { Cterm: 3.55, D: 4.05, E: 4.45, C: 9.0, Y: 10.0 };
  const PK_CTERM = { D: 4.55, E: 4.75 };
  const PK_NTERM = { A: 7.59, M: 7.0, S: 6.93, P: 8.36, T: 6.82, V: 7.44, E: 7.7 };

  const AAs = 'ACDEFGHIKLMNPQRSTVWY';
  function cleanAA(s) { return String(s || '').toUpperCase().replace(/[^A-Z*]/g, ''); }

  /* charge at a given pH, by the same method the reference uses */
  function chargeAt(counts, pH, nterm, cterm) {
    const pos = Object.assign({}, PK_POS), neg = Object.assign({}, PK_NEG);
    if (PK_NTERM[nterm] != null) pos.Nterm = PK_NTERM[nterm];
    if (PK_CTERM[cterm] != null) neg.Cterm = PK_CTERM[cterm];
    let q = 0;
    q += 1 / (1 + Math.pow(10, pH - pos.Nterm));
    ['K', 'R', 'H'].forEach(a => { q += (counts[a] || 0) / (1 + Math.pow(10, pH - pos[a])); });
    q -= 1 / (1 + Math.pow(10, neg.Cterm - pH));
    ['D', 'E', 'C', 'Y'].forEach(a => { q -= (counts[a] || 0) / (1 + Math.pow(10, neg[a] - pH)); });
    return q;
  }

  /* Everything a composition can tell you about a protein, with each number
     reported on its own rather than mixed into a score. */
  E.proteinProps = function (seq) {
    const s = cleanAA(seq).replace(/\*+$/, '');
    const aa = s.replace(/[^A-Z]/g, '');
    const n = aa.length;
    const counts = {};
    for (let i = 0; i < n; i++) counts[aa[i]] = (counts[aa[i]] || 0) + 1;
    const known = AAs.split('').reduce((t, a) => t + (counts[a] || 0), 0);
    const f = a => (counts[a] || 0) / (n || 1);

    let mw = 0;
    for (let i = 0; i < n; i++) mw += E.AA_MONOMER[aa[i]] || 0;
    if (n > 1) mw -= (n - 1) * E.WATER;

    let gravy = 0;
    for (let i = 0; i < n; i++) gravy += E.KD[aa[i]] || 0;
    gravy = known ? gravy / known : 0;

    /* Guruprasad instability index: 10/L times the sum of the dipeptide
       weights along the chain */
    let inst = 0;
    for (let i = 0; i + 1 < n; i++) {
      const w = E.DIWV[aa[i] + aa[i + 1]];
      if (w != null) inst += w;
    }
    inst = n > 1 ? (10 / n) * inst : 0;

    /* Isoelectric point by bisection, bracketed and terminated exactly as
       the reference implementation does, so the two agree to the digit. */
    const nt = aa[0], ct = aa[n - 1];
    let pI = 7.775, lo = 4.05, hi = 12;
    while (hi - lo > 0.0001) {
      if (chargeAt(counts, pI, nt, ct) > 0) lo = pI; else hi = pI;
      pI = (lo + hi) / 2;
    }

    const aliphatic = 100 * (f('A') + 2.9 * f('V') + 3.9 * (f('I') + f('L')));
    const aromatic = f('F') + f('W') + f('Y');
    const turnFrac = f('N') + f('G') + f('P') + f('S');
    const chargeAvg = ((counts.K || 0) + (counts.R || 0) - (counts.D || 0) - (counts.E || 0)) / (n || 1);

    /* Wilkinson-Harrison: the two composition terms that their 1991 work
       found to track solubility. Reported as the two terms and their
       combination, not as a probability - published constants for the
       probability curve differ between sources, and a number carried to
       three decimals would imply a confidence this method does not have. */
    const wh = 15.43 * turnFrac - 29.56 * Math.abs(chargeAvg);

    return {
      length: n, unknown: n - known,
      mw: mw, mwKDa: mw / 1000,
      pI: pI, chargeAt7: chargeAt(counts, 7.0, nt, ct),
      gravy: gravy, instability: inst, aliphatic: aliphatic,
      aromaticity: aromatic,
      counts: counts, fraction: a => f(a),
      cys: counts.C || 0, met: counts.M || 0,
      turnFraction: turnFrac, chargeAverage: chargeAvg, whIndex: wh,
      lowComplexity: E.lowComplexityFraction(aa)
    };
  };

  /* How much of the chain sits in windows of very limited alphabet - a
     rough stand-in for disordered or repetitive stretches, which tend to
     stay soluble themselves but drag structure prediction around. */
  E.lowComplexityFraction = function (aa, win) {
    const w = win || 20, n = aa.length;
    if (n < w) return 0;
    let flagged = 0;
    const mark = new Uint8Array(n);
    for (let i = 0; i + w <= n; i++) {
      const seen = {};
      let distinct = 0;
      for (let k = 0; k < w; k++) { const c = aa[i + k]; if (!seen[c]) { seen[c] = 1; distinct++; } }
      if (distinct <= 7) for (let k = 0; k < w; k++) mark[i + k] = 1;
    }
    for (let i = 0; i < n; i++) if (mark[i]) flagged++;
    return flagged / n;
  };

  /* Stretches hydrophobic enough, and long enough, to cross a membrane.
     A soluble-expression plan for a protein with several of these is the
     wrong plan, whatever else the numbers say. */
  E.tmSegments = function (seq, opts) {
    opts = opts || {};
    const aa = cleanAA(seq).replace(/[^A-Z]/g, '');
    const win = opts.window || 19, cut = opts.threshold == null ? 1.6 : opts.threshold;
    const out = [];
    if (aa.length < win) return out;
    let run = null;
    for (let i = 0; i + win <= aa.length; i++) {
      let h = 0;
      for (let k = 0; k < win; k++) h += E.KD[aa[i + k]] || 0;
      h /= win;
      if (h >= cut) {
        if (run && i <= run.end) { run.end = i + win; run.peak = Math.max(run.peak, h); }
        else { run = { start: i, end: i + win, peak: h }; out.push(run); }
      }
    }
    return out.map(r => ({ start: r.start, end: r.end, length: r.end - r.start,
                           peak: r.peak, seq: aa.slice(r.start, r.end) }));
  };

  /* Short stretches that are both hydrophobic and beta-prone: the kind of
     patch that drives inclusion bodies. Deliberately conservative. */
  E.aggregationRuns = function (seq) {
    const aa = cleanAA(seq).replace(/[^A-Z]/g, '');
    const beta = 'VILFYWMC';
    const out = [];
    let i = 0;
    while (i < aa.length) {
      if (beta.indexOf(aa[i]) < 0) { i++; continue; }
      let j = i;
      while (j < aa.length && beta.indexOf(aa[j]) >= 0) j++;
      if (j - i >= 5) out.push({ start: i, end: j, length: j - i, seq: aa.slice(i, j) });
      i = j;
    }
    return out;
  };

  /* =======================================================================
     WHAT A VECTOR IS

     Two plasmids share a cell only if they replicate from different
     origins and are held by different drugs. Origin families are the hard
     part: everything ColE1-derived - pUC, pBR322, pET, pGEX, pMAL - counts
     as one family however different the names look, so a pET and a pGEX
     will not be maintained together no matter which antibiotics are used.
     Copy numbers are the manufacturers' nominal figures and move with host
     and temperature; they are here to rank, not to quote.
  ======================================================================= */
  E.ORI_FAMILIES = {
    ColE1:   { label: 'ColE1 / pMB1 (pBR322, pUC, pET, pGEX, pMAL)', copy: 40,
               note: 'the commonest backbone family; two of them cannot be co-maintained' },
    p15A:    { label: 'p15A (pACYC, pRARE, pLysS, pGro7)', copy: 12 },
    RSF1030: { label: 'RSF1030 (pRSFDuet)', copy: 100 },
    CloDF13: { label: 'CloDF13 (pCDFDuet)', copy: 25 },
    pSC101:  { label: 'pSC101', copy: 5 },
    pUChigh: { label: 'pUC (high-copy ColE1 derivative)', copy: 500, sameAs: 'ColE1' },
    mammal:  { label: 'mammalian episome / integration', copy: 0 },
    baculo:  { label: 'baculovirus transfer', copy: 0 }
  };
  /* what actually clashes: pUC and ColE1 are one family for this purpose */
  E.oriGroup = o => (E.ORI_FAMILIES[o] && E.ORI_FAMILIES[o].sameAs) || o;

  E.MARKERS = {
    Amp:  { label: 'ampicillin / carbenicillin', note: 'beta-lactamase is secreted, so selection decays in dense culture; carbenicillin holds better' },
    Kan:  { label: 'kanamycin' },
    Cm:   { label: 'chloramphenicol' },
    Spec: { label: 'spectinomycin / streptomycin' },
    Tet:  { label: 'tetracycline' },
    Gent: { label: 'gentamicin' },
    Neo:  { label: 'G418 / neomycin (mammalian)' },
    Hygro:{ label: 'hygromycin (mammalian)' },
    Puro: { label: 'puromycin (mammalian)' }
  };

  /* Fusion partners, with what they are actually for. Sizes are the mature
     tag in residues. */
  E.FUSION_TAGS = {
    His:   { label: 'His6/His8', size: 1, solubility: 0, purify: 'IMAC', note: 'purification only; no help with folding' },
    GST:   { label: 'GST', size: 220, solubility: 2, purify: 'glutathione', note: 'dimerises, which can be a problem for a protein that must stay monomeric' },
    MBP:   { label: 'MBP', size: 370, solubility: 3, purify: 'amylose (or His, if the vector adds one)', note: 'the strongest general solubiliser; large, so cleave it off before any structural work' },
    SUMO:  { label: 'SUMO', size: 100, solubility: 2, purify: 'IMAC', note: 'cleaved by Ulp1 to leave a native N-terminus, with no scar' },
    Trx:   { label: 'thioredoxin', size: 109, solubility: 2, purify: 'IMAC (if His is present)', note: 'helps most with proteins that need a reducing, folding-friendly environment' },
    NusA:  { label: 'NusA', size: 495, solubility: 3, purify: 'IMAC (if His is present)', note: 'very large and very effective; the fusion dominates the mass' },
    GB1:   { label: 'GB1', size: 56, solubility: 2, purify: 'IMAC (if His is present)', note: 'small, highly soluble, well behaved in NMR' },
    Halo:  { label: 'HaloTag', size: 297, solubility: 2, purify: 'HaloLink' },
    Intein:{ label: 'intein-CBD', size: 200, solubility: 1, purify: 'chitin', note: 'self-cleaving, so no protease is needed, but cleavage is finicky' },
    Strep: { label: 'Strep-tag II', size: 1, solubility: 0, purify: 'Strep-Tactin' },
    FLAG:  { label: 'FLAG', size: 1, solubility: 0, purify: 'anti-FLAG' },
    Bio:   { label: 'biotin acceptor peptide', size: 15, solubility: 0, purify: 'streptavidin' }
  };

  /* The vector families a lab of this kind actually holds. Matched loosely
     by name, because inventories are typed by hand. First match wins, so
     the specific patterns come before the general ones. */
  E.VECTOR_KB = [
    /* ---- the Duet set: built to be used together ---- */
    { re: /^pacyc-?duet/i, name: 'pACYCDuet-1', system: 'E. coli', promoter: 'T7lac', induce: 'IPTG',
      ori: 'p15A', marker: 'Cm', tags: ['His'], mcs: 2, needsT7: true,
      note: 'two cassettes; designed to sit alongside pET, pRSF and pCDF' },
    { re: /^prsf-?duet/i, name: 'pRSFDuet-1', system: 'E. coli', promoter: 'T7lac', induce: 'IPTG',
      ori: 'RSF1030', marker: 'Kan', tags: ['His'], mcs: 2, needsT7: true,
      note: 'highest copy of the Duet set, so it gives the most of whatever goes in it' },
    { re: /^pcdf/i, name: 'pCDFDuet-1', system: 'E. coli', promoter: 'T7lac', induce: 'IPTG',
      ori: 'CloDF13', marker: 'Spec', tags: ['His'], mcs: 2, needsT7: true },
    { re: /^petduet/i, name: 'pETDuet-1', system: 'E. coli', promoter: 'T7lac', induce: 'IPTG',
      ori: 'ColE1', marker: 'Amp', tags: ['His'], mcs: 2, needsT7: true },
    /* ---- helper plasmids that take up a slot of their own ---- */
    { re: /^pgro7|^pg-kje8|^pkje7|^pg-tf2|^ptf16/i, name: 'chaperone plasmid (Takara set)', system: 'E. coli',
      promoter: 'araB / Pzt-1', induce: 'arabinose (and tetracycline for some)', ori: 'p15A', marker: 'Cm',
      helper: 'chaperones', tags: [],
      note: 'GroEL/GroES and friends; occupies the p15A slot and the Cm marker, so it cannot be used with pACYC or a Rosetta host' },
    { re: /^prare|rosetta/i, name: 'pRARE (rare-codon tRNAs)', system: 'E. coli', promoter: 'native tRNA',
      ori: 'p15A', marker: 'Cm', helper: 'tRNA', tags: [],
      note: 'supplies tRNAs for codons E. coli uses rarely; carried by Rosetta hosts' },
    { re: /^plys[se]/i, name: 'pLysS / pLysE', system: 'E. coli', ori: 'p15A', marker: 'Cm',
      helper: 'T7 lysozyme', tags: [], note: 'damps leaky T7 expression of a toxic product' },
    /* ---- fusion vectors, most specific first ---- */
    { re: /nusa/i, name: 'NusA fusion vector', system: 'E. coli', promoter: 'T7lac', induce: 'IPTG',
      ori: 'ColE1', marker: 'Amp', tags: ['His', 'NusA'], cleave: 'TEV', needsT7: true },
    { re: /sumo/i, name: 'SUMO fusion vector', system: 'E. coli', promoter: 'T7lac', induce: 'IPTG',
      ori: 'ColE1', marker: 'Kan', tags: ['His', 'SUMO'], cleave: 'Ulp1 / SUMO protease', needsT7: true },
    { re: /trx/i, name: 'thioredoxin fusion vector', system: 'E. coli', promoter: 'T7lac', induce: 'IPTG',
      ori: 'ColE1', marker: 'Amp', tags: ['His', 'Trx'], cleave: 'TEV', needsT7: true },
    { re: /pmal|6hismal|mbp|mal-?[pc0-9]/i, name: 'MBP fusion vector', system: 'E. coli', promoter: 'tac', induce: 'IPTG',
      ori: 'ColE1', marker: 'Amp', tags: ['His', 'MBP'], cleave: 'TEV',
      note: 'a pMAL-p backbone sends the product to the periplasm; pMAL-c keeps it in the cytoplasm' },
    { re: /gb1|ngb/i, name: 'GB1 fusion vector', system: 'E. coli', promoter: 'T7lac', induce: 'IPTG',
      ori: 'ColE1', marker: 'Kan', tags: ['His', 'GB1'], cleave: 'TEV', needsT7: true },
    { re: /^pgex|gst-parallel|^pacyc-gst|^popink/i, name: 'GST fusion vector', system: 'E. coli',
      promoter: 'tac', induce: 'IPTG', ori: 'ColE1', marker: 'Amp', tags: ['GST'],
      cleave: 'thrombin, PreScission or TEV depending on the variant' },
    { re: /^ptxb|^ptyb|^ptby/i, name: 'IMPACT intein vector', system: 'E. coli', promoter: 'T7lac',
      induce: 'IPTG', ori: 'ColE1', marker: 'Amp', tags: ['Intein'], needsT7: true,
      note: 'self-cleaving; no protease needed, but the cleavage step is temperamental' },
    { re: /halotag|^phtc|^phtn/i, name: 'HaloTag vector', system: 'mammalian', promoter: 'CMV',
      ori: 'mammal', marker: 'Neo', tags: ['Halo'] },
    /* ---- plain expression backbones ---- */
    { re: /^pnic|^p28bioh|bio3/i, name: 'pNIC (SGC pET derivative)', system: 'E. coli', promoter: 'T7lac',
      induce: 'IPTG', ori: 'ColE1', marker: 'Kan', tags: ['His'], cleave: 'TEV', needsT7: true },
    { re: /^pet-?(28|24|26|29|30|31|41|42|45)/i, name: 'pET (kanamycin)', system: 'E. coli',
      promoter: 'T7lac', induce: 'IPTG', ori: 'ColE1', marker: 'Kan', tags: ['His'], needsT7: true },
    { re: /^pet/i, name: 'pET (ampicillin)', system: 'E. coli', promoter: 'T7lac', induce: 'IPTG',
      ori: 'ColE1', marker: 'Amp', tags: ['His'], needsT7: true },
    { re: /^pq(e|iq)/i, name: 'pQE (Qiagen T5)', system: 'E. coli', promoter: 'T5-lac', induce: 'IPTG',
      ori: 'ColE1', marker: 'Amp', tags: ['His'],
      note: 'read by E. coli’s own polymerase, so no DE3 lysogen is needed - but the host must supply lacIq' },
    { re: /^pcw/i, name: 'pCW tet-regulated', system: 'E. coli', promoter: 'tac / tet', induce: 'IPTG or tetracycline',
      ori: 'ColE1', marker: 'Amp', tags: ['His'] },
    { re: /^pst39/i, name: 'pST39 polycistronic', system: 'E. coli', promoter: 'T7', induce: 'IPTG',
      ori: 'ColE1', marker: 'Amp', tags: [], needsT7: true, mcs: 4,
      note: 'four cassettes on one plasmid: a complex from a single selection' },
    { re: /parallel|^pmcgs|^pgtvl/i, name: 'T7 expression vector', system: 'E. coli', promoter: 'T7',
      induce: 'IPTG', ori: 'ColE1', marker: 'Amp', tags: ['His'], needsT7: true },
    /* ---- other hosts ---- */
    { re: /^pacyc/i, name: 'pACYC derivative', system: 'E. coli', promoter: 'varies',
      ori: 'p15A', marker: 'Cm', tags: [],
      note: 'a p15A backbone: compatible with ColE1 vectors, but it clashes with any other p15A or chloramphenicol plasmid' },
    { re: /^pfastbac|^pfb|^pfbnbg|hy_pmt/i, name: 'baculovirus transfer vector', system: 'insect',
      promoter: 'polyhedrin', ori: 'baculo', marker: 'Amp', tags: ['His'],
      note: 'insect cells fold difficult eukaryotic proteins that E. coli will not' },
    { re: /^pcdna|^pcmv|^rccmv|^rc\/cmv|^pci|^pcag|^pece|^pcs3|myc-pcs3|^pebg|^pegfp|^plp-cmv|^pnlf|^pbs-cmv|^pcmx|^pflag-cmv/i,
      name: 'mammalian expression vector', system: 'mammalian', promoter: 'CMV', ori: 'mammal',
      marker: 'Amp', tags: [], note: 'for proteins that need mammalian folding and modification' },
    { re: /^pdonr|^pentr/i, name: 'Gateway entry clone', system: 'cloning only', ori: 'ColE1', marker: 'Kan',
      tags: [], note: 'a holding vector: the insert is recombined onward into a destination vector' },
    { re: /^puc|^pbluescript|^pcr-blunt|topo|^pot[b2]|^pflc|^pkw|^poat/i, name: 'cloning vector',
      system: 'cloning only', ori: 'pUChigh', marker: 'Amp', tags: [],
      note: 'holds and propagates a sequence; it does not express it' }
  ];

  /* Work out what a vector is from how it was written down. */
  E.classifyVector = function (name) {
    const s = String(name || '').trim();
    if (!s) return null;
    for (let i = 0; i < E.VECTOR_KB.length; i++) {
      if (E.VECTOR_KB[i].re.test(s)) {
        const k = E.VECTOR_KB[i];
        const out = Object.assign({}, k);
        delete out.re;
        out.matched = s;
        out.copy = (E.ORI_FAMILIES[k.ori] || {}).copy || 0;
        out.group = E.oriGroup(k.ori);
        /* a pMAL-p backbone exports to the periplasm, which changes the plan */
        if (/pmal-?p/i.test(s)) { out.periplasm = true; out.note = (out.note || '') + ' This one is the periplasmic version.'; }
        return out;
      }
    }
    return null;
  };

  /* Hosts, and what they already carry. A strain that comes with a Cm
     plasmid has spent that slot before you start. */
  E.HOSTS = [
    { re: /lobstr.*ril|ril.*lobstr/i, name: 'LOBSTR-RIL', t7: true, carries: { ori: 'p15A', marker: 'Cm', what: 'rare-codon tRNAs' },
      note: 'LOBSTR with the RIL tRNAs: fewer His-tag contaminants, and rare codons covered' },
    { re: /lobstr/i, name: 'LOBSTR', t7: true,
      note: 'ArnA and SlyD altered, so far less of the usual background binds the IMAC column' },
    { re: /rosetta|prare/i, name: 'Rosetta / pRARE', t7: true, carries: { ori: 'p15A', marker: 'Cm', what: 'rare-codon tRNAs' } },
    { re: /shuffle/i, name: 'SHuffle T7', t7: true, oxidising: true,
      note: 'cytoplasm is oxidising and carries DsbC, so disulfides can form and isomerise inside the cell' },
    { re: /origami/i, name: 'Origami', t7: true, oxidising: true },
    { re: /arctic/i, name: 'ArcticExpress', t7: true, cold: true,
      note: 'carries cold-adapted chaperonins for growth and folding at 10-13 C' },
    { re: /c4[13]/i, name: 'C41 / C43 (Walker)', t7: true, toxic: true,
      note: 'tolerates products that kill ordinary BL21, including membrane proteins' },
    { re: /plyss/i, name: 'BL21(DE3) pLysS', t7: true, carries: { ori: 'p15A', marker: 'Cm', what: 'T7 lysozyme' } },
    { re: /bl21/i, name: 'BL21(DE3)', t7: true },
    { re: /hek|293|expi/i, name: 'HEK293', system: 'mammalian' },
    { re: /sf9|sf21|high ?five|insect/i, name: 'insect cells', system: 'insect' }
  ];
  /* Everything a strain already carries, whether that came from the lab's
     own sheet or from what is known about the strain by name. */
  E.hostCarries = function (host) {
    if (!host) return [];
    const out = [];
    const list = Array.isArray(host.carries) ? host.carries
               : (host.carries && (host.carries.ori || host.carries.marker)) ? [host.carries] : [];
    list.forEach(c => { if (c && (c.ori || c.marker)) out.push(c); });
    /* Every drug the strain is held by is spent, whichever of its plasmids
       carries which. A sheet that reads "Chl, Str" for two resident
       plasmids does not say which is which, and pairing them off in the
       order written gets it wrong as often as right - so both are taken. */
    (host.markers || []).forEach(m => {
      if (!out.some(c => c.marker === m)) out.push({ marker: m, what: 'its own selection (' + m + ')' });
    });
    return out;
  };

  E.classifyHost = function (name) {
    const s = String(name || '');
    for (let i = 0; i < E.HOSTS.length; i++) {
      if (E.HOSTS[i].re.test(s)) { const o = Object.assign({}, E.HOSTS[i]); delete o.re; return o; }
    }
    return null;
  };

  /* =======================================================================
     WHAT THIS PROTEIN ACTUALLY IS

     Composition is blind to identity. A kinase domain and a random coil of
     the same amino acids read the same way on hydropathy and charge, and
     only one of them is a famously bad actor in E. coli. So before any of
     the composition arithmetic is believed, the sequence is checked for
     the handful of conserved motifs that carry hard expression
     consequences - a catalytic kinase fold, a cyclin-dependent kinase's
     signature helix, a zinc finger that needs zinc, a signal peptide that
     should not be there at all.

     Each is a published, checkable motif, and each is reported with the
     position it was found at, so a claim can be argued with rather than
     taken on faith. What they contribute is not a score but a fact about
     the protein that outranks anything the composition suggests.
  ======================================================================= */
  E.FAMILY_RULES = [
    {
      id: 'kinase',
      name: 'protein kinase catalytic domain',
      /* the two most conserved landmarks of the eukaryotic kinase fold:
         the catalytic loop HRD..N and the Mg-binding DFG, which sit a
         characteristic distance apart */
      find: function (aa) {
        /* The catalytic loop is conserved as [HY]-[RC]-D-x-K-x-x-N rather
           than literally HRD: PKA and a good many others carry YRD, and
           the two residues before the asparagine are unrestricted. */
        const hrd = /[HY][RC]D[LIVMFA][KRHA].{2}[NGS]/.exec(aa);
        if (!hrd) return null;
        const after = aa.slice(hrd.index + 8, hrd.index + 8 + 80);
        const dfg = /D[FLWY]G/.exec(after);
        if (!dfg) return null;
        const gloop = /G.G..G.V/.exec(aa.slice(0, Math.max(0, hrd.index)));
        return { at: hrd.index, evidence:
          'catalytic loop ' + hrd[0] + ' at residue ' + (hrd.index + 1) +
          ', ' + dfg[0] + ' motif ' + (dfg.index + 8) + ' residues later' +
          (gloop ? ', glycine-rich loop ' + gloop[0] + ' at ' + (gloop.index + 1) : '') };
      },
      risk: 'high',
      says: 'Eukaryotic protein kinases are one of the hardest classes to get out of E. coli in soluble, folded form. The catalytic domain folds slowly, the activation loop is disordered until it is phosphorylated, and the usual result of a bare tag at 37 \u00b0C is inclusion bodies.',
      does: ['a solubilising fusion rather than a plain His tag',
             'induction at 16-18 \u00b0C, never at 37',
             'co-expression with lambda phosphatase if the protein is toxic or heterogeneously phosphorylated',
             'insect cells as the fallback, which is where most commercial kinase is made']
    },
    {
      id: 'cdk',
      name: 'cyclin-dependent kinase',
      needs: 'kinase',
      find: function (aa) {
        const m = /P[SICFLNV][ST][AT][ILV][RLK]E/.exec(aa);
        if (!m) return null;
        return { at: m.index, evidence: 'the ' + m[0] + ' helix at residue ' + (m.index + 1) +
                 ', the signature of the cyclin-binding helix in this family' };
      },
      risk: 'high',
      partner: 'its cyclin',
      says: 'A CDK on its own is not a folded protein in the way a single-domain enzyme is. The PSTAIRE-type helix and the activation loop only take up their ordered positions when a cyclin binds, so the monomer expressed alone is partly disordered and aggregates readily, whatever its amino-acid composition suggests.',
      does: ['co-express the cyclin \u2014 this is the single thing that most changes the outcome',
             'failing that, a GST, MBP or SUMO fusion, which is how commercially sold CDK is made',
             'a Cks/Suc1-family partner also helps for several CDKs',
             '16-18 \u00b0C induction']
    },
    {
      id: 'ploop',
      name: 'P-loop NTPase (Walker A and B)',
      find: function (aa) {
        /* Walker A, then a Walker B downstream: the pair together is
           specific where the A motif alone is not - GxxxxGK[ST] turns up
           by chance in long sequences often enough to be useless alone */
        const a = /[AG].{4}GK[ST]/.exec(aa);
        if (!a) return null;
        const tail = aa.slice(a.index + 8, a.index + 260);
        const b = /[ILVFYM]{4}D[ED]/.exec(tail);
        if (!b) return null;
        return { at: a.index, evidence: 'Walker A ' + a[0] + ' at residue ' + (a.index + 1) +
                 ', Walker B ' + b[0] + ' ' + (b.index + 8) + ' residues later' };
      },
      risk: 'medium',
      says: 'This binds and hydrolyses a nucleotide. Many such proteins are only properly folded with it bound, and the empty form is the one that aggregates.',
      does: ['keep ATP or ADP with magnesium in the lysis and gel-filtration buffers',
             'a non-hydrolysable analogue if the protein chews through ATP during purification',
             'expect the apo protein to behave worse than the loaded one, not better']
    },
    {
      id: 'wd40',
      name: 'WD40 beta-propeller',
      find: function (aa) {
        const at = [];
        for (let i = 0; i + 2 <= aa.length; i++) if (aa[i] === 'W' && aa[i + 1] === 'D') at.push(i);
        /* the blades repeat every 40 residues or so; four of them spaced
           like that is a propeller, two is a coincidence */
        let runs = 0;
        for (let i = 1; i < at.length; i++) { const d = at[i] - at[i - 1]; if (d >= 30 && d <= 55) runs++; }
        if (runs < 3) return null;
        return { at: at[0], count: runs + 1,
                 evidence: (runs + 1) + ' WD motifs spaced 30-55 residues apart, first at residue ' + (at[0] + 1) };
      },
      risk: 'medium',
      says: 'A seven-bladed propeller closes only when every blade is present, so a construct that cuts through the repeat region gives a protein that cannot fold at all. Propellers also fold slowly, and many of them are built to sit on a partner rather than alone.',
      does: ['take the domain boundaries from a structure or an alignment, never from a guess \u2014 one blade short is fatal',
             'induce cold and give it time',
             'if it still will not fold, co-express the partner it binds in the cell']
    },
    {
      id: 'protease',
      name: 'chymotrypsin-family serine protease',
      find: function (aa) {
        const m = /GDSGG[PAST]/.exec(aa);
        if (!m) return null;
        return { at: m.index, evidence: 'the GDSGG active-site motif at residue ' + (m.index + 1) };
      },
      risk: 'high',
      says: 'An active protease digests itself and everything around it. Expressed in a cytoplasm it is either toxic, degraded, or both.',
      does: ['express the zymogen, not the mature enzyme, and activate it afterwards',
             'or use inclusion bodies deliberately and refold \u2014 for this class that is the normal route, not a failure',
             'inhibitors in every buffer from the moment of lysis']
    },
    {
      id: 'collagen',
      name: 'collagen-like Gly-X-Y repeat',
      find: function (aa) {
        const m = /(G..){12,}/.exec(aa);
        if (!m) return null;
        const seg = m[0];
        const pro = (seg.match(/P/g) || []).length / seg.length;
        if (pro < 0.15) return null;
        return { at: m.index, evidence: Math.floor(seg.length / 3) + ' Gly-X-Y triplets from residue ' +
                 (m.index + 1) + ', ' + Math.round(pro * 100) + '% proline' };
      },
      risk: 'high',
      says: 'A collagen triple helix is held together by hydroxyproline, and E. coli has no prolyl hydroxylase. What comes out of a bacterial culture is the unmodified chain, which does not form a stable helix at body temperature.',
      does: ['use an insect or mammalian host, where the modification happens',
             'or co-express prolyl-4-hydroxylase, which is the established bacterial workaround',
             'expect a lower melting temperature than the native protein whatever you do']
    },
    {
      id: 'ankyrin',
      name: 'ankyrin repeat',
      find: function (aa) {
        const re = /[TSV]PLH[LMVIA]A/g;
        const hits = []; let m;
        while ((m = re.exec(aa))) { hits.push(m.index + 1); re.lastIndex = m.index + 1; }
        if (hits.length < 3) return null;
        return { at: hits[0] - 1, count: hits.length,
                 evidence: hits.length + ' ankyrin repeats, first at residue ' + hits[0] };
      },
      risk: 'good',
      says: 'Ankyrin repeats stack into a solenoid that folds locally and robustly. This is one of the better-behaved folds in E. coli, which is why designed binders are built on it.',
      does: ['a plain tag is a fair first attempt here',
             'keep whole repeats in the construct: cutting through one leaves a frayed, sticky end']
    },
    {
      id: 'rrm',
      name: 'RNA recognition motif',
      find: function (aa) {
        const rnp1 = /[KR]G[FY][GA][FY][ILVA].[FY]/.exec(aa);
        if (!rnp1) return null;
        const before = aa.slice(Math.max(0, rnp1.index - 90), rnp1.index);
        const rnp2 = /[ILV][FY][ILV].[NGK]L/.exec(before);
        if (!rnp2) return null;
        return { at: rnp1.index, evidence: 'RNP2 ' + rnp2[0] + ' followed by RNP1 ' + rnp1[0] +
                 ' at residue ' + (rnp1.index + 1) };
      },
      risk: 'medium',
      says: 'This binds RNA, and in a bacterial lysate it will bind whatever bacterial RNA it can find. The usual symptom is a protein that purifies well, runs as a smear on gel filtration and has an A260/A280 well above 0.6.',
      does: ['wash with high salt, or add a nuclease step, before the protein ever reaches the column',
             'check the A260/A280 of the eluate rather than assuming it is clean']
    },
    {
      id: 'c4zinc',
      name: 'C4-type zinc module',
      find: function (aa) {
        const re = /C.{2}C.{9,20}C.{2}C/g;
        const hits = []; let m;
        while ((m = re.exec(aa))) { hits.push(m.index + 1); re.lastIndex = m.index + 1; }
        if (hits.length < 2) return null;
        return { at: hits[0] - 1, count: hits.length,
                 evidence: hits.length + ' C4 zinc-binding modules, first at residue ' + hits[0] };
      },
      risk: 'medium',
      says: 'Paired cysteine quartets of this kind hold a structural zinc. Without the metal the module is unfolded, and with EDTA anywhere in the workflow it will lose it.',
      does: ['supplement the culture with zinc at induction',
             'no EDTA in any buffer, including the one the protease comes in']
    },
    {
      id: 'zf',
      name: 'C2H2 zinc finger',
      find: function (aa) {
        const re = /C.{2,4}C.{10,14}H.{3,5}H/g;
        const hits = []; let m;
        while ((m = re.exec(aa))) { hits.push(m.index + 1); re.lastIndex = m.index + 1; }
        if (!hits.length) return null;
        return { at: hits[0] - 1, count: hits.length,
                 evidence: hits.length + ' C2H2 finger' + (hits.length === 1 ? '' : 's') +
                   ', first at residue ' + hits[0] };
      },
      risk: 'medium',
      says: 'Zinc fingers do not fold without zinc. In a cytoplasm that has none to spare, the apo form is a sticky, partly unfolded peptide.',
      does: ['add 50-100 \u00b5M ZnCl2 or ZnSO4 to the medium at induction',
             'keep a little zinc in the lysis and purification buffers',
             'avoid EDTA throughout, which will strip the metal straight back out']
    },
    {
      id: 'signal',
      name: 'N-terminal signal peptide',
      find: function (aa) {
        const head = aa.slice(0, 40);
        if (head.length < 20) return null;
        /* A real signal peptide has three parts, and all three are needed:
           a short positively charged n-region, a genuinely hydrophobic
           h-region, and a small-residue cleavage site. Without the charge
           requirement an alanine-rich N-terminus reads as a signal
           peptide, which it is not. */
        const nRegion = /[KR]/.test(head.slice(0, 6));
        if (!nRegion) return null;
        const core = /[AILMFVWC]{8,}/.exec(head.slice(3, 30));
        if (!core) return null;
        let h = 0;
        for (let i = 0; i < core[0].length; i++) h += E.KD[core[0][i]] || 0;
        if (h / core[0].length < 1.8) return null;      /* alanine alone is not hydrophobic enough */
        const ala = /[AGS][^DEKRHP][AGS]/.test(head.slice(10, 35));
        if (!ala) return null;
        return { at: 0, evidence: 'a basic n-region, a hydrophobic core at residues ' + (core.index + 4) +
                 '-' + (core.index + 3 + core[0].length) + ', and a small-residue cleavage site behind it' };
      },
      risk: 'medium',
      says: 'The N-terminus looks like a secretion signal. If this came from a eukaryotic protein, E. coli will not process it properly, and leaving it on is a common cause of a construct that expresses nothing useful.',
      does: ['express the mature protein, starting after the predicted cleavage site',
             'or keep it and use a periplasmic vector such as a pMAL-p backbone, where the export is the point']
    },
    {
      id: 'glyco',
      name: 'secreted glycoprotein',
      needs: 'signal',
      find: function (aa) {
        const re = /N[^P][ST]/g;
        const hits = []; let m;
        while ((m = re.exec(aa))) { hits.push(m.index + 1); re.lastIndex = m.index + 1; }
        if (hits.length < 3) return null;
        return { at: hits[0] - 1, count: hits.length,
                 evidence: hits.length + ' N-x-S/T glycosylation sequons behind a signal peptide' };
      },
      risk: 'high',
      says: 'A signal peptide and a run of glycosylation sequons together say this protein is normally made in the secretory pathway and carries sugars. E. coli adds none, and for many such proteins the sugar is part of what keeps them folded and soluble.',
      does: ['an insect or mammalian host, which is what this class is usually made in',
             'if it must be E. coli, express the mature chain without the signal peptide and expect to work at it',
             'a periplasmic route at least gives you the disulfides, if not the sugars']
    },
    {
      id: 'ef',
      name: 'EF-hand calcium-binding loop',
      find: function (aa) {
        const re = /D.[DNS].[DNSG]G.[ILVFYW].[DENSTG]/g;
        const hits = []; let m;
        while ((m = re.exec(aa))) { hits.push(m.index + 1); re.lastIndex = m.index + 1; }
        if (hits.length < 2) return null;
        return { at: hits[0] - 1, count: hits.length,
                 evidence: hits.length + ' EF-hand-like loops, first at residue ' + hits[0] };
      },
      risk: 'low',
      says: 'Paired EF-hand loops bind calcium, and many such proteins are only stable with it bound.',
      does: ['keep a little CaCl2 in the buffers rather than EDTA, unless the apo form is what you want']
    }
  ];

  /* What the sequence is, as opposed to what it is made of. */
  E.recogniseFamily = function (seq) {
    const aa = String(seq || '').toUpperCase().replace(/[^A-Z]/g, '');
    const found = {}, out = [];
    E.FAMILY_RULES.forEach(rule => {
      if (rule.needs && !found[rule.needs]) return;
      let hit = null;
      try { hit = rule.find(aa); } catch (e) { hit = null; }
      if (!hit) return;
      found[rule.id] = true;
      out.push({ id: rule.id, name: rule.name, risk: rule.risk, partner: rule.partner,
                 says: rule.says, does: rule.does, at: hit.at, count: hit.count || 1,
                 evidence: hit.evidence });
    });
    return out;
  };

  /* =======================================================================
     HOW HARD IS THIS GOING TO BE

     Each descriptor is reported with its own value and what it implies.
     No composition-based method can tell you whether a protein will fold;
     what these can do is point at the specific thing likely to go wrong,
     which is what decides the vector, the tag and the temperature.
  ======================================================================= */
  E.solubilityOutlook = function (seq, opts) {
    opts = opts || {};
    const isDNA = /^[ACGTUNacgtun\s]+$/.test(String(seq || '')) && String(seq).replace(/\s/g, '').length > 60;
    let aa = seq, dna = null;
    if (isDNA) { dna = E.cleanSeq(seq); aa = E.translate(dna).replace(/\*+$/, ''); }
    aa = String(aa || '').toUpperCase().replace(/[^A-Z]/g, '');
    if (!aa) return null;
    const p = E.proteinProps(aa);
    const tm = E.tmSegments(aa);
    const agg = E.aggregationRuns(aa);
    const families = E.recogniseFamily(aa);
    const f = [];
    const add = (name, value, verdict, why) => f.push({ name: name, value: value, verdict: verdict, why: why });

    /* size */
    const kda = p.mwKDa;
    add('Size', kda.toFixed(1) + ' kDa, ' + p.length + ' residues',
        kda > 100 ? 'bad' : kda > 60 ? 'watch' : 'good',
        kda > 100 ? 'Large single chains fold slowly and give aggregation more chance to win. Expressing a domain at a time usually works where the whole thing does not.'
        : kda > 60 ? 'On the large side for E. coli; worth having a domain boundary in mind as a fallback.'
        : 'A size E. coli handles comfortably.');

    /* membrane helices */
    if (tm.length >= 3)
      add('Predicted membrane helices', tm.length + ' hydrophobic stretches of ≥ 19 residues', 'bad',
          'This looks like an integral membrane protein. Soluble cytoplasmic expression is the wrong target for it — it needs a detergent or nanodisc route, or a construct that leaves the membrane region out.');
    else if (tm.length)
      add('Predicted membrane helices', tm.length + ' hydrophobic stretch' + (tm.length === 1 ? '' : 'es'), 'watch',
          'One strongly hydrophobic stretch may be a signal peptide, a buried helix, or a genuine anchor. Check where it falls before deciding.');
    else add('Predicted membrane helices', 'none', 'good', 'Nothing that reads as a membrane anchor.');

    /* hydropathy */
    add('Hydropathy (GRAVY)', p.gravy.toFixed(3), p.gravy > 0.0 ? 'bad' : p.gravy > -0.2 ? 'watch' : 'good',
        p.gravy > 0 ? 'Net hydrophobic over the whole chain, which tracks strongly with inclusion bodies. A solubilising fusion is close to essential here.'
        : p.gravy > -0.2 ? 'Borderline; a solubilising fusion is worth having.'
        : 'Comfortably hydrophilic overall.');

    /* cysteines */
    const cys = p.cys;
    add('Cysteines', cys + (cys ? (cys % 2 ? ' (odd)' : ' (even)') : ''),
        cys >= 6 ? 'watch' : 'good',
        cys === 0 ? 'No cysteines, so disulfides are not a consideration.'
        : cys >= 6 ? 'Enough cysteines that disulfides are likely. The ordinary E. coli cytoplasm is reducing and will not form them: an oxidising strain such as SHuffle, or export to the periplasm, is the usual answer.'
        : 'Few enough that free cysteines are more likely than structural disulfides; keep reductant in every buffer.');

    /* pI against the buffer */
    const buf = opts.bufferPH == null ? 7.5 : opts.bufferPH;
    const near = Math.abs(p.pI - buf) < 1.0;
    add('Isoelectric point', p.pI.toFixed(2), near ? 'watch' : 'good',
        near ? 'The pI sits within a unit of your buffer pH (' + buf + '), where a protein carries least charge and is most likely to come out of solution. Moving the buffer a point either way often fixes a protein that looks hopeless.'
             : 'Far enough from the working pH to carry a useful net charge.');

    /* instability index */
    add('Instability index', p.instability.toFixed(1), p.instability > 40 ? 'watch' : 'good',
        p.instability > 40 ? 'Above 40, which Guruprasad’s index calls unstable. It predicts a short half-life in vivo rather than insolubility, so expect yield to suffer more than folding.'
                           : 'Below 40, the stable side of the index.');

    /* aggregation-prone runs */
    const longRuns = agg.filter(r => r.length >= 7);
    add('Hydrophobic / beta-prone runs', agg.length + ' of ≥ 5 residues' + (longRuns.length ? ', ' + longRuns.length + ' of ≥ 7' : ''),
        longRuns.length >= 3 ? 'watch' : 'good',
        longRuns.length >= 3 ? 'Several long hydrophobic stretches, the kind of patch that nucleates inclusion bodies. Slower expression at lower temperature gives folding time to win.'
                             : 'No great excess of sticky stretches.');

    /* low complexity / disorder */
    add('Low-complexity content', Math.round(p.lowComplexity * 100) + '%',
        p.lowComplexity > 0.3 ? 'watch' : 'good',
        p.lowComplexity > 0.3 ? 'A large disordered or repetitive fraction. Such regions usually stay soluble but resist crystallisation, and they invite proteolysis — trimming them often improves both yield and behaviour.'
                              : 'Mostly ordered sequence.');

    /* Wilkinson-Harrison */
    add('Wilkinson–Harrison index', p.whIndex.toFixed(2) +
        '  (turn-forming ' + (p.turnFraction * 100).toFixed(1) + '%, charge average ' + p.chargeAverage.toFixed(3) + ')',
        p.whIndex < -1.5 ? 'watch' : 'good',
        'A composition rule from 1991 that weighs turn-forming residues against charge imbalance. Treat it as a weak hint, not a verdict — it is right rather more often than chance and no more than that.');

    /* How open the start of the message is. Translation begins by the
       ribosome landing on the region around the start codon, and a stem
       sitting over it suppresses initiation however good the promoter is
       - which shows up as a construct that is perfectly correct and makes
       nothing. Only checkable when DNA was given, since it is a property
       of the message rather than of the protein. */
    if (dna && dna.length >= 60) {
      const open = E.startAccessibility(dna, opts.utr);
      if (open) {
        add('Start-codon accessibility', Math.round(open.paired * 100) + '% of the first ' +
            open.window + ' bases are paired' + (open.overStart ? ', including the start codon itself' : ''),
            open.paired > 0.6 || open.overStart ? 'watch' : 'good',
            open.paired > 0.6 || open.overStart
              ? 'The message folds over its own start. A stem across the ribosome binding site and the first codons is one of the commonest reasons a construct that is entirely correct still makes nothing: silent changes in the first few codons, or a different 5\u2032 region, usually open it up.'
              : 'The start of the message is largely unpaired, so the ribosome can get on to it.');
      }
    }

    /* rare codons, when DNA was given */
    if (dna) {
      const rare = E.rareCodonLoad(dna);
      add('Rare codons for E. coli', rare.count + ' of ' + rare.total + ' (' + (rare.fraction * 100).toFixed(1) + '%)' +
          (rare.worstRun > 1 ? ', longest run ' + rare.worstRun : ''),
          rare.fraction > 0.08 || rare.worstRun >= 3 ? 'watch' : 'good',
          rare.fraction > 0.08 || rare.worstRun >= 3
            ? 'Enough rare codons to stall translation, which both lowers yield and gives partly made chains time to aggregate. A tRNA-supplementing host such as Rosetta or a RIL strain, or a synthetic gene recoded for E. coli, removes the problem.'
            : 'Codon usage is not going to be the limiting factor.');
    }

    /* What the protein is, put in front of what it is made of. A reading
       taken from composition alone would call a cyclin-dependent kinase
       unremarkable - the amino acids are unremarkable - and send you off
       with a bare His tag. Knowing it is a CDK settles the question
       before the arithmetic gets a say. */
    families.forEach(fam => {
      f.unshift({ name: 'Recognised as: ' + fam.name, value: fam.evidence,
                  verdict: fam.risk === 'high' ? 'bad' : fam.risk === 'good' ? 'good' : 'watch',
                  family: fam.id, why: fam.says });
    });

    /* the overall reading */
    const bad = f.filter(x => x.verdict === 'bad' && !x.family).length;
    const watch = f.filter(x => x.verdict === 'watch' && !x.family).length;
    const membrane = tm.length >= 3;
    const famHigh = families.filter(x => x.risk === 'high');
    const famWatch = families.filter(x => x.risk === 'medium');
    let band = membrane ? 'wrong target'
             : bad >= 1 || watch >= 4 ? 'hard' : watch >= 2 ? 'tricky' : 'straightforward';
    /* A known-difficult family is a floor on the verdict, never a ceiling.
       A family that is good news - a solenoid repeat, say - says nothing
       against the protein and must not drag the reading down with it. */
    if (famHigh.length && band !== 'wrong target') band = 'hard';
    else if (famWatch.length && band === 'straightforward') band = 'tricky';

    let summary;
    if (membrane) summary = 'This reads as a membrane protein, and no amount of fusion tag makes a membrane protein soluble.';
    else if (famHigh.length) {
      summary = famHigh[0].says + ' Nothing in the amino-acid composition would have told you that, and on composition alone this looks ordinary \u2014 which is exactly why it is worth saying.';
    } else if (band === 'hard') summary = 'Several things here point towards inclusion bodies. Plan for a strong solubilising fusion and a slow, cold induction from the start rather than after a failure.';
    else if (band === 'tricky') summary = 'Nothing disqualifying, but enough to be worth a solubilising fusion rather than a bare His tag.';
    else summary = 'Nothing in the composition argues against straightforward expression. A plain His tag is a reasonable first attempt.';

    /* Silence is not the same as a clean bill of health, and the two are
       easy to confuse. If nothing was recognised, the reading rests on
       composition alone - which is the weakest thing here, and the reader
       should be told so rather than left to infer it. */
    if (!families.length)
      f.push({ name: 'Recognised as', value: 'nothing the tool knows',
               verdict: 'watch', family: 'none',
               why: 'No fold or motif the tool carries was found, so everything above rests on composition alone. ' +
                    'That is the weakest kind of evidence there is: a cyclin-dependent kinase reads as perfectly ' +
                    'ordinary on composition. If you already know this protein needs a partner, a cofactor, a ' +
                    'modification E. coli cannot make, or has a reputation for going into inclusion bodies, trust that over anything here.' });

    return { protein: aa, dna: dna, props: p, tm: tm, aggregation: agg, families: families,
             factors: f, band: band, summary: summary, concerns: { bad: bad, watch: watch } };
  };

  /* Is the start of the message free enough for a ribosome to land on it?
     Folded over the first stretch of the coding sequence, with whatever
     5' untranslated region was supplied in front of it. */
  E.startAccessibility = function (dna, utr) {
    const s = (E.cleanSeq(utr || '') + E.cleanSeq(dna)).replace(/U/g, 'T');
    const win = Math.min(90, s.length);
    const head = s.slice(0, win);
    let fold;
    try { fold = E.rnaFold(head); } catch (e) { return null; }
    if (!fold || fold.error || !fold.structure) return null;
    const startAt = E.cleanSeq(utr || '').length;
    const paired = (fold.structure.match(/[()]/g) || []).length / win;
    /* is the start codon itself inside a stem */
    const around = fold.structure.slice(Math.max(0, startAt - 4), startAt + 8);
    const overStart = (around.match(/[()]/g) || []).length >= 8;
    return { window: win, paired: paired, overStart: overStart, structure: fold.structure };
  };

  /* Codons E. coli uses rarely, and whether they come in clusters. */
  E.RARE_CODONS = ['AGG', 'AGA', 'CGA', 'CGG', 'AUA', 'ATA', 'CTA', 'CUA', 'CCC', 'GGA'];
  E.rareCodonLoad = function (dna) {
    const s = E.cleanSeq(dna).replace(/U/g, 'T');
    const rare = { AGG: 1, AGA: 1, CGA: 1, CGG: 1, ATA: 1, CTA: 1, CCC: 1, GGA: 1 };
    let count = 0, total = 0, run = 0, worst = 0;
    for (let i = 0; i + 3 <= s.length; i += 3) {
      const c = s.substr(i, 3);
      total++;
      if (rare[c]) { count++; run++; if (run > worst) worst = run; }
      else run = 0;
    }
    return { count: count, total: total, fraction: total ? count / total : 0, worstRun: worst };
  };

  /* =======================================================================
     WILL THESE LIVE IN THE SAME CELL
  ======================================================================= */
  E.compatibility = function (items, opts) {
    opts = opts || {};
    const host = opts.host ? (typeof opts.host === 'string' ? E.classifyHost(opts.host) : opts.host) : null;
    const slots = items.slice();
    /* A strain is not an empty vessel. Several of the useful ones arrive
       carrying a plasmid of their own - rare-codon tRNAs, T7 lysozyme,
       chaperones - and that plasmid has spent an origin and a drug before
       anything of yours goes in. Some carry two. */
    E.hostCarries(host).forEach(c => {
      slots.push({ name: host.name + ' (carries ' + (c.what || c.plasmid || 'a plasmid') + ')',
                   ori: c.ori, marker: c.marker, group: E.oriGroup(c.ori), fromHost: true });
    });
    const problems = [], notes = [];
    /* replication origins */
    const byGroup = {};
    slots.forEach(s => {
      const g = s.group || E.oriGroup(s.ori);
      if (!g || g === 'mammal' || g === 'baculo') return;
      (byGroup[g] = byGroup[g] || []).push(s);
    });
    Object.keys(byGroup).forEach(g => {
      if (byGroup[g].length > 1)
        problems.push({ kind: 'origin', severity: 'blocking',
          what: byGroup[g].map(s => s.name).join(' and ') + ' all replicate from ' +
                ((E.ORI_FAMILIES[g] || {}).label || g) + '.',
          why: 'Plasmids of one origin family are not kept together: within a few generations most cells hold only one of them, whatever the antibiotics say.' });
    });
    /* selection markers */
    const byMarker = {};
    slots.forEach(s => { if (s.marker) (byMarker[s.marker] = byMarker[s.marker] || []).push(s); });
    Object.keys(byMarker).forEach(m => {
      if (byMarker[m].length > 1)
        problems.push({ kind: 'marker', severity: 'blocking',
          what: byMarker[m].map(s => s.name).join(' and ') + ' are all selected with ' +
                ((E.MARKERS[m] || {}).label || m) + '.',
          why: 'One drug cannot hold two plasmids: the cell keeps whichever it likes and the other is lost.' });
    });
    /* T7 polymerase */
    const needsT7 = slots.filter(s => s.needsT7);
    if (needsT7.length && host && host.t7 === false)
      problems.push({ kind: 'host', severity: 'blocking',
        what: needsT7.map(s => s.name).join(', ') + ' need T7 polymerase, and ' + host.name + ' has none.',
        why: 'A T7 promoter is silent without a DE3 lysogen or an equivalent source of the polymerase.' });
    if (needsT7.length > 1)
      notes.push('All ' + needsT7.length + ' of these run from T7 promoters, so one dose of IPTG starts them together and they compete for the same polymerase. The higher-copy plasmid tends to win; if the ratio matters, put the partner you need more of on the higher-copy backbone, or move one onto a promoter with a different inducer.');
    /* ampicillin in long cultures */
    if (byMarker.Amp && byMarker.Amp.length)
      notes.push('Ampicillin selection decays as beta-lactamase accumulates in the medium, which matters most in the long cold inductions these plans usually call for. Carbenicillin holds better.');
    /* copy balance */
    const copies = slots.filter(s => s.copy).map(s => s.name + ' ~' + s.copy);
    if (copies.length > 1) notes.push('Nominal copy numbers: ' + copies.join(', ') + '. They set the rough ratio you will get before induction is taken into account.');
    return { ok: problems.length === 0, problems: problems, notes: notes, slots: slots };
  };

  /* =======================================================================
     WHAT TO BUILD

     Ranked strategies, drawn from what the lab actually holds rather than
     from a catalogue. The reasons are the point: each option says what it
     is answering.
  ======================================================================= */
  E.strategyFor = function (outlook, inventory, opts) {
    opts = opts || {};
    const inv = inventory || [];
    const want = [];
    const o = outlook;
    const hard = o.band === 'hard', tricky = o.band === 'tricky';
    const membrane = o.band === 'wrong target';
    const cys = o.props.cys;
    const bigCys = cys >= 6;

    /* which fusion partners are worth having, and why, for this protein */
    const fams = o.families || [];
    const cdk = fams.filter(x => x.id === 'cdk')[0];
    const kinase = fams.filter(x => x.id === 'kinase')[0];
    const wishlist = [];
    if (membrane) wishlist.push({ tag: null, why: 'A membrane protein does not become soluble by fusion. Either remove the membrane region and express the soluble domain, or move to a detergent route.' });
    if (membrane) { /* nothing else to offer: the question itself is wrong */ }
    else if (cdk || kinase) {
      /* a kinase is not a composition problem, and the ranking says so:
         GST, MBP and SUMO are what commercially sold kinase is made with,
         and a bare tag is the one thing not to start with */
      wishlist.push({ tag: 'GST', why: 'how commercially sold ' + (cdk ? 'CDK' : 'kinase') +
        ' is actually made \u2014 a large, well-behaved partner that keeps the catalytic domain in solution while it folds' });
      wishlist.push({ tag: 'MBP', why: 'the strongest general solubiliser, and the usual alternative to GST for a kinase domain' });
      wishlist.push({ tag: 'SUMO', why: 'strong, and cleaves to leave a native N-terminus, which matters if the first residues are part of the fold' });
      wishlist.push({ tag: 'GB1', why: 'small, if you would rather not carry several hundred residues of fusion through purification' });
      wishlist.push({ tag: 'His', why: 'worth one plate as a control, but a bare tag on a kinase domain usually ends in the pellet' });
    }
    else if (hard) {
      wishlist.push({ tag: 'MBP', why: 'the strongest general solubiliser, and the usual first answer to a protein that goes into inclusion bodies' });
      wishlist.push({ tag: 'NusA', why: 'comparable to MBP in effect, worth trying when MBP has not worked' });
      wishlist.push({ tag: 'SUMO', why: 'strong solubilising effect and cleaves to leave a native N-terminus' });
      wishlist.push({ tag: 'GB1', why: 'small and very soluble, so it adds little mass to carry through purification' });
    } else if (tricky) {
      wishlist.push({ tag: 'SUMO', why: 'a modest solubility gain with a clean, scarless cut' });
      wishlist.push({ tag: 'GB1', why: 'small, soluble, and barely changes the protein it carries' });
      wishlist.push({ tag: 'MBP', why: 'the fallback if the smaller tags are not enough' });
      wishlist.push({ tag: 'His', why: 'worth one attempt first — nothing here says it must fail' });
    } else {
      wishlist.push({ tag: 'His', why: 'nothing in the composition argues against a plain tag' });
      wishlist.push({ tag: 'GB1', why: 'cheap insurance if the first attempt disappoints' });
      wishlist.push({ tag: 'SUMO', why: 'the next step up if needed' });
    }
    if (bigCys && !membrane) wishlist.push({ tag: 'Trx', why: 'suits proteins with several cysteines, and pairs with an oxidising strain' });

    /* match the wishlist against the shelf */
    const seen = {};
    wishlist.forEach((w, rank) => {
      if (!w.tag) { want.push({ tag: null, why: w.why, have: [], rank: rank }); return; }
      const have = inv.filter(v => v.kb && v.kb.tags && v.kb.tags.indexOf(w.tag) >= 0 &&
                                   v.kb.system === 'E. coli' && !v.kb.helper);
      const key = w.tag;
      if (seen[key]) return; seen[key] = 1;
      want.push({ tag: w.tag, info: E.FUSION_TAGS[w.tag], why: w.why, have: have, rank: rank });
    });

    /* conditions, which matter as much as the vector */
    const conditions = [];
    /* a partner is not a condition you tune, it is part of the construct */
    if (cdk && !(opts.partners || []).some(n => /cyclin|ccn/i.test(n)))
      conditions.push({ what: 'Co-express the cyclin',
        detail: 'This is the one change that matters most for a CDK. The monomer does not take up its folded form until a cyclin binds, so expressing it alone is asking a half-folded protein to stay in solution. Put the cyclin on a compatible second plasmid and induce them together.' });
    fams.forEach(fam => {
      (fam.does || []).forEach((d, i) => {
        if (/co-express the cyclin/i.test(d) && cdk) return;    /* already said, above */
        conditions.push({ what: i === 0 ? 'Because it is a ' + fam.name : '\u2026also', detail: d });
      });
    });
    if (hard || tricky) conditions.push({ what: 'Induce cold', detail: 'Drop to 16-18 °C before adding IPTG and leave it overnight. Slower synthesis is the single most reliable way to move protein from the pellet to the supernatant, and it costs nothing to try.' });
    if (hard) conditions.push({ what: 'Induce lightly', detail: 'Use 0.1-0.25 mM IPTG rather than 1 mM. Less polymerase pressure gives each chain more time to fold before the next one is made.' });
    if (bigCys) conditions.push({ what: 'Consider an oxidising cytoplasm', detail: 'With ' + cys + ' cysteines, disulfides are likely. SHuffle or Origami will form them inside the cell; the periplasm via a pMAL-p backbone is the other route.' });
    if (o.props.mwKDa > 100) conditions.push({ what: 'Split it', detail: 'At ' + o.props.mwKDa.toFixed(0) + ' kDa, expressing defined domains separately usually succeeds where the whole chain does not.' });
    const rareFactor = o.factors.filter(x => /Rare codons/.test(x.name))[0];
    if (rareFactor && rareFactor.verdict !== 'good')
      conditions.push({ what: 'Cover the rare codons', detail: 'A RIL or Rosetta host supplies the missing tRNAs. Note that those tRNAs ride on a chloramphenicol p15A plasmid, which then rules out any other p15A or Cm plasmid in the same cell.' });
    const pIf = o.factors.filter(x => /Isoelectric/.test(x.name))[0];
    if (pIf && pIf.verdict === 'watch')
      conditions.push({ what: 'Move the buffer pH', detail: 'With a pI of ' + o.props.pI.toFixed(1) + ', work at least a unit away from it. A protein that looks insoluble at its pI is often perfectly soluble one pH unit off.' });
    if (!conditions.length) conditions.push({ what: 'Standard conditions are a fair first try', detail: '37 °C to mid-log, then 0.5 mM IPTG at 18 °C overnight is a sensible default for a first attempt.' });

    return { options: want, conditions: conditions, band: o.band };
  };

  /* =======================================================================
     READING A PLASMID INVENTORY
  ======================================================================= */
  const PHEAD = [
    ['id', /^(plasmid ?id|id|stock ?id|clone ?id)$/i],
    ['name', /^(construct ?id|construct|name|description|insert)$/i],
    ['vector', /^(vector|backbone|plasmid)$/i],
    ['resistance', /^(resistance|marker|selection|antibiotic|res)$/i],
    ['tag', /^(tag|tags|fusion)$/i],
    ['target', /^(target|gene|protein)$/i],
    ['organism', /^(organism|species|source)$/i],
    ['type', /^(type|system|use)$/i],
    ['domain', /^(domain\/design|domain|design|boundaries)$/i],
    ['conc', /^(ng\/ul|ng\/µl|conc\.?|concentration)$/i],
    ['expression', /^(expression.*|expresses\??|yield)$/i],
    ['notes', /^(notes.*|note|comment.*)$/i],
    ['shelf', /^(shelf|level)$/i],
    ['box', /^(box|plate|rack)$/i],
    ['well', /^(position|well|pos|slot)$/i],
    ['room', /^(location|rm|room)$/i],
    ['seq', /^(dna ?sequence|sequence)$/i]
  ];

  /* What a recorded expression note actually says. These are free text in
     every lab, so the reading is deliberately cautious: anything unclear
     comes back as unknown rather than as a guess. */
  E.expressionOutcome = function (text) {
    const s = String(text || '').trim();
    if (!s) return null;
    const low = s.toLowerCase();
    let soluble = null;
    if (/\binsoluble\b|inclusion bod/.test(low)) soluble = false;
    else if (/\bsoluble\b/.test(low)) soluble = true;
    if (/insoluble\?|soluble\?/.test(low)) soluble = null;      /* they were unsure */
    let level = null;
    if (/\bnot seen|no expression|\bno\b|low\/no/.test(low)) level = 'none';
    else if (/very high|\bhigh\b/.test(low)) level = 'high';
    else if (/\bgood\b|expresses|\bmedium\b|\bmoderate\b/.test(low)) level = 'moderate';
    else if (/\blow\b/.test(low)) level = 'low';
    const host = E.classifyHost(s);
    const t = s.match(/(\d{2})\s*°?\s*c\b/i);
    return { text: s, soluble: soluble, level: level,
             host: host ? host.name : null, temperature: t ? +t[1] : null };
  };

  E.parsePlasmidStock = function (text) {
    const rows = [];
    const lines = String(text || '').replace(/\r\n?/g, '\n').split('\n').filter(l => l.trim());
    if (!lines.length) return { plasmids: [], skipped: [], columns: null };
    /* Tabs if there are any, otherwise proper comma-separated fields with
       quoting, and only as a last resort the column-aligned split. Treating
       a run of spaces as a separator while also splitting on commas
       quietly cuts quoted fields in half - "N-6his, C-strep tag II" became
       two columns and shifted every field after it along by one. */
    const splitCSV = l => {
      const out = [];
      let cur = '', q = false;
      for (let i = 0; i < l.length; i++) {
        const c = l[i];
        if (q) {
          if (c === '"') { if (l[i + 1] === '"') { cur += '"'; i++; } else q = false; }
          else cur += c;
        } else if (c === '"') q = true;
        else if (c === ',') { out.push(cur); cur = ''; }
        else cur += c;
      }
      out.push(cur);
      return out;
    };
    const split = l => l.indexOf('\t') >= 0 ? l.split('\t')
                     : l.indexOf(',') >= 0 ? splitCSV(l)
                     : l.split(/\s{2,}/);
    let cols = null, start = 0;
    for (let i = 0; i < Math.min(6, lines.length); i++) {
      const cells = split(lines[i]).map(c => c.replace(/^"|"$/g, '').trim());
      const map = {};
      let hits = 0;
      cells.forEach((c, j) => {
        for (const [field, re] of PHEAD) {
          if (re.test(c) && map[field] == null) { map[field] = j; hits++; break; }
        }
      });
      if (map.vector != null && hits >= 3) { cols = map; start = i + 1; break; }
    }
    if (!cols) return { plasmids: [], skipped: lines.slice(0, 3).map(l => ({ line: l, why: 'no header row recognised' })), columns: null };
    const skipped = [];
    for (let i = start; i < lines.length; i++) {
      const cells = split(lines[i]).map(c => c.replace(/^"|"$/g, '').trim());
      const g = k => (cols[k] != null && cells[cols[k]] != null) ? cells[cols[k]].trim() : '';
      const vector = g('vector');
      if (!vector || /^\?+$/.test(vector)) { skipped.push({ line: lines[i].slice(0, 60), why: 'no vector named' }); continue; }
      const loc = [g('room'), g('shelf'), g('box'), g('well')].filter(Boolean).join(' / ');
      const rec = {
        id: g('id'), name: g('name') || g('target') || g('id'), vector: vector,
        resistance: g('resistance'), tag: g('tag'), target: g('target'),
        organism: g('organism'), type: g('type'), domain: g('domain'),
        conc: parseFloat(g('conc')) || null, loc: loc,
        expressionText: g('expression'), notes: g('notes'), seq: g('seq')
      };
      rec.kb = E.classifyVector(vector);
      rec.outcome = E.expressionOutcome(rec.expressionText);
      rows.push(rec);
    }
    return { plasmids: rows, skipped: skipped, columns: cols };
  };

  /* Every distinct vector in an inventory, with what is known about it and
     what the lab's own records say it did. */
  E.vectorSummary = function (plasmids) {
    const by = {};
    (plasmids || []).forEach(p => {
      const key = (p.kb && p.kb.name) || p.vector;
      const e = by[key] = by[key] || { name: key, kb: p.kb, spellings: {}, count: 0,
                                       soluble: 0, insoluble: 0, noExpression: 0, outcomes: [], examples: [] };
      e.count++;
      e.spellings[p.vector] = (e.spellings[p.vector] || 0) + 1;
      if (e.examples.length < 4 && p.loc) e.examples.push(p);
      if (p.outcome) {
        e.outcomes.push(p.outcome);
        if (p.outcome.soluble === true) e.soluble++;
        else if (p.outcome.soluble === false) e.insoluble++;
        if (p.outcome.level === 'none') e.noExpression++;
      }
    });
    return Object.keys(by).map(k => by[k]).sort((a, b) => b.count - a.count);
  };

  /* =======================================================================
     THE COMPETENT CELLS ON THE SHELF

     Which strain to use is half the answer, and a lab's freezer settles it
     more than any catalogue does. The sheet is authoritative for what each
     strain carries and what holds it; what the tool adds is the rest -
     whether the strain makes T7 polymerase, whether its cytoplasm is
     oxidising, whether it tolerates a product that kills BL21.
  ======================================================================= */
  const CHEAD = [
    ['name', /^(strain|cells?|name)$/i],
    ['carries', /^(existing plasmid|plasmid|carries|resident)/i],
    ['use', /^(use ?case|use|purpose)/i],
    ['resistance', /^(resistance|marker|antibiotic|selection)/i],
    ['stock', /^(stock\??|in ?stock|available|have)/i],
    ['competency', /^(competency|competence|method)/i],
    ['note', /^(note|notes|comment)/i],
    ['loc', /^(location|where|freezer|box)/i]
  ];
  const MARKER_WORDS = [
    [/chl|cam\b|chloramph/i, 'Cm'], [/kan/i, 'Kan'], [/amp|carb/i, 'Amp'],
    [/spec|\bstr\b|strep/i, 'Spec'], [/tet/i, 'Tet'], [/gent|genta/i, 'Gent']
  ];
  /* An origin named outright, in a plasmid name or in a note. */
  E.oriFromName = function (text) {
    const t = String(text || '');
    if (/psc101/i.test(t)) return 'pSC101';
    if (/pacyc|p15a/i.test(t)) return 'p15A';
    if (/rsf ?1030|^prsf/i.test(t)) return 'RSF1030';
    if (/clodf|^pcdf/i.test(t)) return 'CloDF13';
    if (/cole1|pmb1|pbr322/i.test(t)) return 'ColE1';
    return null;
  };

  E.markersIn = function (text) {
    const out = [];
    String(text || '').split(/[,;/]/).forEach(part => {
      MARKER_WORDS.forEach(([re, m]) => { if (re.test(part) && out.indexOf(m) < 0) out.push(m); });
    });
    return out;
  };

  /* Read a competent-cell sheet. The resident plasmid and the drug come
     from the sheet, because that is what is in the freezer; the origin
     each plasmid replicates from comes from the vector knowledge. */
  E.parseCellStock = function (text) {
    const lines = String(text || '').replace(/\r\n?/g, '\n').split('\n').filter(l => l.trim());
    if (!lines.length) return { strains: [], skipped: [] };
    const splitCSV = l => {
      const out = []; let cur = '', q = false;
      for (let i = 0; i < l.length; i++) {
        const c = l[i];
        if (q) { if (c === '"') { if (l[i + 1] === '"') { cur += '"'; i++; } else q = false; } else cur += c; }
        else if (c === '"') q = true;
        else if (c === ',') { out.push(cur); cur = ''; }
        else cur += c;
      }
      out.push(cur); return out;
    };
    const split = l => l.indexOf('\t') >= 0 ? l.split('\t') : l.indexOf(',') >= 0 ? splitCSV(l) : l.split(/\s{2,}/);
    let cols = null, start = 0;
    for (let i = 0; i < Math.min(5, lines.length); i++) {
      const cells = split(lines[i]).map(c => c.trim());
      const map = {}; let hits = 0;
      cells.forEach((c, j) => {
        for (const [f, re] of CHEAD) if (re.test(c) && map[f] == null) { map[f] = j; hits++; break; }
      });
      if (map.name != null && hits >= 2) { cols = map; start = i + 1; break; }
    }
    if (!cols) return { strains: [], skipped: [{ line: lines[0].slice(0, 60), why: 'no header row with a Strain column' }] };
    const strains = [], skipped = [];
    for (let i = start; i < lines.length; i++) {
      const cells = split(lines[i]).map(c => c.trim());
      const g = k => (cols[k] != null && cells[cols[k]] != null) ? cells[cols[k]].trim() : '';
      const name = g('name');
      if (!name) { skipped.push({ line: lines[i].slice(0, 60), why: 'no strain named' }); continue; }
      /* The resistance cell often carries an instruction in brackets -
         "(only use Str for comp cell prep)" - which names a drug that is
         not a resistance the strain holds. The instruction is dropped
         before the drugs are read, or every such strain would appear to
         block a marker it does not use. */
      const note = g('note');
      const markers = E.markersIn(g('resistance').replace(/\([^)]*\)/g, ' '));
      /* each resident plasmid, with the origin it replicates from */
      const carries = g('carries').split(/[;,]/).map(x => x.trim()).filter(Boolean).map(pl => {
        const kb = E.classifyVector(pl);
        /* a resident plasmid is often written as its origin rather than as
           a vector name - "pSC101" is the replicon, not a catalogue item */
        const bare = E.oriFromName(pl) || E.oriFromName(note);
        return { plasmid: pl, what: pl, ori: (kb && kb.ori) || bare || null };
      });
      /* An extra origin named only in the notes still occupies a slot:
         "pACYC184 derivative" means p15A whatever the plasmid is called. */
      const noteOri = E.oriFromName(note);
      if (noteOri && !carries.some(c => E.oriGroup(c.ori) === E.oriGroup(noteOri)))
        carries.push({ plasmid: note.slice(0, 40), what: 'the backbone named in its note', ori: noteOri, fromNote: true });
      const known = E.classifyHost(name) || {};
      const use = g('use');
      strains.push({
        name: name, carries: carries,
        markers: markers,
        use: use, expression: /express/i.test(use), cloning: /replicat|clon/i.test(use),
        inStock: /^y/i.test(g('stock')),
        competency: g('competency'), note: g('note'), loc: g('loc'),
        /* Whether the strain makes T7 polymerase decides whether a pET is
           silent in it, so it is read from the note as well as the name:
           a sheet says "BL21 (DE3) derivative" for strains whose names
           give nothing away, and names get misspelled. */
        t7: known.t7 != null ? known.t7
            : /\(?de3\)?|bl21|lobstr|ros+et+a|origami|arctic|shuffle|solubl21|nico|c4[13]/i.test(name + ' ' + note),
        oxidising: !!known.oxidising || /origami|shuffle/i.test(name),
        cold: !!known.cold || /arctic/i.test(name),
        toxic: !!known.toxic || /c4[13]/i.test(name),
        ccdB: /db3\.?1/i.test(name),
        chaperones: carries.some(c => /gro|kje|tf\d|chaperon/i.test(c.plasmid || '')) || /arctic/i.test(name),
        rareCodons: carries.some(c => /rare|ril|rare2/i.test(c.plasmid || '')) || /ril|ros+et+a|codon ?\+/i.test(name)
      });
    }
    return { strains: strains, skipped: skipped };
  };

  /* Which strains suit this protein, given what the reading said. */
  E.pickStrains = function (outlook, strains, opts) {
    opts = opts || {};
    const list = (strains || []).filter(s => s.expression !== false);
    const o = outlook;
    const needOx = o && o.props.cys >= 6;
    const needRare = o && o.factors.some(f => /Rare codons/.test(f.name) && f.verdict !== 'good');
    const hard = o && (o.band === 'hard' || o.band === 'wrong target');
    const out = list.map(s => {
      const why = [];
      let score = 0;
      if (s.inStock) { score += 3; } else { why.push('not made up as competent cells yet'); }
      if (!s.t7) { score -= 6; why.push('no T7 polymerase, so a pET-style vector would stay silent in it'); }
      if (s.cloning && !s.expression) { score -= 10; why.push('a cloning strain, not an expression one'); }
      if (needOx && s.oxidising) { score += 5; why.push('oxidising cytoplasm, which suits the cysteines in this one'); }
      if (needRare && s.rareCodons) { score += 4; why.push('supplies the rare tRNAs this sequence needs'); }
      if (hard && s.chaperones) { score += 3; why.push('carries chaperones, which is worth having for a protein this awkward'); }
      if (hard && s.cold) { score += 2; why.push('built for folding at low temperature'); }
      if (/lobstr|nico/i.test(s.name)) { score += 2; why.push('less background on the IMAC column'); }
      if (s.toxic) { score += 1; why.push('tolerates products that kill ordinary BL21'); }
      if (s.carries.length) why.push('already carries ' + s.carries.map(c => c.plasmid).join(' and ') +
        ', which takes up ' + (s.carries.length === 1 ? 'an origin and a drug' : 'origins and drugs') + ' you then cannot use');
      return { strain: s, score: score, why: why };
    }).sort((a, b) => b.score - a.score);
    return out;
  };

  /* =======================================================================
     CO-EXPRESSION

     Several proteins in one cell. The constraint is not preference: two
     plasmids of one origin family will not stay together, and one drug
     cannot hold two plasmids. The Duet set exists precisely because its
     four members take four different origins and four different markers.
  ======================================================================= */
  E.COEX_SLOTS = [
    { ori: 'ColE1',   marker: 'Amp',  vector: 'pETDuet-1 / any pET-Amp',  copy: 40 },
    { ori: 'RSF1030', marker: 'Kan',  vector: 'pRSFDuet-1',               copy: 100 },
    { ori: 'p15A',    marker: 'Cm',   vector: 'pACYCDuet-1',              copy: 12 },
    { ori: 'CloDF13', marker: 'Spec', vector: 'pCDFDuet-1',               copy: 25 }
  ];

  E.coexpressionPlan = function (partners, inventory, opts) {
    opts = opts || {};
    /* the host may arrive as a name, or as a strain already read off the
       lab's own sheet - which is the better of the two, since it carries
       what that strain actually holds */
    const host = !opts.host ? null
               : (typeof opts.host === 'string' ? E.classifyHost(opts.host) : opts.host);
    const inv = inventory || [];
    const taken = {};
    const clashes = [];
    E.hostCarries(host).forEach(c => {
      if (c.ori) taken[E.oriGroup(c.ori)] = true;
      if (c.marker) taken['marker:' + c.marker] = true;
    });
    if (opts.chaperones) {
      /* the chaperone plasmids are themselves p15A and chloramphenicol, so
         a host that already carries one of those has no room for them */
      const clash = E.hostCarries(host).filter(c => E.oriGroup(c.ori) === 'p15A' || c.marker === 'Cm')[0];
      if (clash)
        clashes.push({ kind: 'helper', severity: 'blocking',
          what: 'The chaperone plasmid and ' + host.name + '\u2019s own ' + (clash.what || clash.plasmid) +
                ' both want p15A and chloramphenicol.',
          why: 'They cannot be held in the same cell. Use a host without a carried plasmid \u2014 plain BL21(DE3) or LOBSTR \u2014 if you want the chaperones, or drop the chaperones and keep the tRNAs.' });
      taken.p15A = true; taken['marker:Cm'] = true;
    }

    /* most of what you need first, so it lands on the highest copy */
    const order = partners.map((p, i) => Object.assign({ index: i }, p))
      .sort((a, b) => (b.want || 1) - (a.want || 1));
    let slots = E.COEX_SLOTS.filter(s => !taken[E.oriGroup(s.ori)] && !taken['marker:' + s.marker]);
    /* When one partner is wanted in greater amount than another, the
       backbones are ranked by copy number so the one you need most of sits
       on the highest. With nothing to choose between them the conventional
       order stands, since that is what most shelves are stocked for. */
    const wants = partners.map(p => p.want || 1);
    if (wants.some(w => w !== wants[0])) slots = slots.slice().sort((a, b) => b.copy - a.copy);
    const assigned = [], unplaced = [];
    order.forEach((p, i) => {
      if (i < slots.length) {
        const s = slots[i];
        /* anything on the shelf that already sits on this backbone */
        const have = inv.filter(v => v.kb && v.kb.ori === s.ori && v.kb.marker === s.marker &&
                                     v.kb.system === 'E. coli');
        assigned.push({ partner: p, slot: s, have: have.slice(0, 6), haveCount: have.length });
      } else unplaced.push(p);
    });

    const items = assigned.map(a => ({ name: a.partner.name || ('partner ' + (a.partner.index + 1)),
                                       ori: a.slot.ori, marker: a.slot.marker, group: E.oriGroup(a.slot.ori),
                                       copy: a.slot.copy, needsT7: true }));
    const compat = E.compatibility(items, { host: host });
    compat.problems = clashes.concat(compat.problems);
    if (clashes.length) compat.ok = false;

    const advice = [];
    if (unplaced.length)
      advice.push({ kind: 'blocking', text: unplaced.length + ' partner' + (unplaced.length === 1 ? '' : 's') +
        ' could not be given a slot: only ' + slots.length + ' compatible origin/marker combinations are left' +
        (E.hostCarries(host).length ? ', because ' + host.name + ' already occupies ' +
          (E.hostCarries(host).length === 1 ? 'one' : String(E.hostCarries(host).length)) + ' with its own plasmid' + (E.hostCarries(host).length === 1 ? '' : 's') : '') +
        (opts.chaperones ? ', and the chaperone plasmid occupies another' : '') +
        '. Put two of them on one Duet plasmid — each carries two cassettes — or use a polycistronic vector such as pST39, which holds four.' });
    if (assigned.length > 1)
      advice.push({ kind: 'note', text: 'Transform them together and select on all ' + assigned.length +
        ' drugs at once. Keeping several plasmids is a burden on the cell, so grow from a fresh co-transformation rather than from a glycerol stock that has been passaged.' });
    if (assigned.length > 2)
      advice.push({ kind: 'note', text: 'Three or more plasmids in one cell is where yields usually fall away. If any two partners can go on a single Duet backbone, that is almost always the better construct.' });
    if (opts.chaperones)
      advice.push({ kind: 'note', text: 'The chaperone plasmid is induced with arabinose, separately from IPTG: give the chaperones a head start so they are already present when the target appears.' });
    return { assigned: assigned, unplaced: unplaced, compatibility: compat, advice: advice,
             host: host, slotsAvailable: slots.length };
  };

  /* =======================================================================
     EVERYTHING ON THE SHELF, AND EVERY COMBINATION OF IT

     Not "here is a tag, here are vectors that carry it", but the whole
     cross-product: every expression backbone the lab holds against every
     strain it can be put in, each one scored and each one either valid or
     struck out with a reason. The useful output of that is short - half a
     dozen combinations worth trying - but it is short because everything
     else was considered and rejected, not because it was never looked at.
  ======================================================================= */

  /* What the collection contains, counted by element rather than by tube. */
  E.elementInventory = function (plasmids, strains) {
    const inv = { origins: {}, promoters: {}, markers: {}, tags: {}, systems: {},
                  hosts: [], vectors: {} };
    (plasmids || []).forEach(p => {
      const k = p.kb;
      if (!k) return;
      const bump = (m, key) => { if (key) m[key] = (m[key] || 0) + 1; };
      bump(inv.origins, k.ori);
      bump(inv.promoters, k.promoter);
      bump(inv.markers, k.marker);
      bump(inv.systems, k.system);
      (k.tags || []).forEach(t => bump(inv.tags, t));
      const vk = k.name;
      const v = inv.vectors[vk] = inv.vectors[vk] || { kb: k, n: 0, spellings: {}, where: null,
                                                       soluble: 0, insoluble: 0, none: 0 };
      v.n++;
      v.spellings[p.vector] = (v.spellings[p.vector] || 0) + 1;
      if (!v.where && p.loc) v.where = p.loc;
      if (p.outcome) {
        if (p.outcome.soluble === true) v.soluble++;
        else if (p.outcome.soluble === false) v.insoluble++;
        if (p.outcome.level === 'none') v.none++;
      }
    });
    inv.hosts = (strains || []).filter(s => s.expression !== false);
    return inv;
  };

  /* Every backbone-and-strain pairing, scored. */
  E.rankCombinations = function (outlook, plasmids, strains, opts) {
    opts = opts || {};
    const inv = E.elementInventory(plasmids, strains);
    const strat = E.strategyFor(outlook, plasmids || [], opts);
    /* the order the fusions came out in, as a lookup */
    const tagRank = {};
    strat.options.forEach((o, i) => { if (o.tag) tagRank[o.tag] = i; });
    const nTags = Math.max(1, Object.keys(tagRank).length);

    const hosts = (strains || []).filter(s => s.expression !== false);
    const fallbackHost = { name: 'BL21(DE3)', t7: true, carries: [], markers: [],
                           inStock: true, expression: true };
    const hostList = hosts.length ? hosts : [fallbackHost];
    const hostScore = {};
    E.pickStrains(outlook, hostList).forEach(x => { hostScore[x.strain.name] = x; });

    const out = [], rejected = [];
    Object.keys(inv.vectors).forEach(vname => {
      const v = inv.vectors[vname];
      const k = v.kb;
      if (k.system !== 'E. coli' || k.helper) return;      /* not an E. coli expression backbone */
      hostList.forEach(h => {
        const item = { name: vname, ori: k.ori, marker: k.marker, group: k.group,
                       needsT7: k.needsT7, copy: k.copy };
        const compat = E.compatibility([item], { host: h });
        const why = [], against = [];
        let score = 0;

        /* does it carry a fusion this protein actually wants */
        const tags = (k.tags || []).filter(t => t !== 'His');
        const best = tags.map(t => tagRank[t]).filter(r => r != null).sort((a, b) => a - b)[0];
        if (best != null) {
          score += (nTags - best) * 6;
          why.push('carries ' + ((E.FUSION_TAGS[tags[tags.map(t => tagRank[t]).indexOf(best)]] || {}).label || tags[0]) +
                   ', which is ' + (best === 0 ? 'the first fusion this protein calls for' : 'number ' + (best + 1) + ' on the list for this protein'));
        } else if (tagRank.His != null && (k.tags || []).indexOf('His') >= 0) {
          score += (nTags - tagRank.His) * 2;
          why.push('a plain affinity tag only');
        }

        /* the lab's own record for this backbone */
        if (v.soluble || v.insoluble || v.none) {
          score += v.soluble * 4 - v.insoluble * 4 - v.none * 2;
          why.push('your own records on this backbone: ' +
            [v.soluble ? v.soluble + ' soluble' : null,
             v.insoluble ? v.insoluble + ' insoluble' : null,
             v.none ? v.none + ' no expression' : null].filter(Boolean).join(', '));
        }

        /* how many of them are in the freezer */
        score += Math.min(6, Math.log(1 + v.n) * 2);

        /* the strain */
        const hs = hostScore[h.name];
        if (hs) { score += hs.score; hs.why.slice(0, 2).forEach(w => why.push(w)); }

        if (!compat.ok) {
          compat.problems.forEach(pr => against.push(pr.what + ' ' + pr.why));
        }
        if (k.needsT7 && h.t7 === false)
          against.push('this backbone runs from a T7 promoter and ' + h.name + ' makes no T7 polymerase, so it would stay silent.');
        if (!h.inStock) why.push('these cells are not made up yet');

        const common = Object.keys(v.spellings).sort((a, b) => v.spellings[b] - v.spellings[a])[0];
      const row = { vector: vname, spelling: common, kb: k, have: v.n,
                      where: v.where, host: h, score: score, why: why, against: against,
                      records: { soluble: v.soluble, insoluble: v.insoluble, none: v.none },
                      ok: against.length === 0 };
        (row.ok ? out : rejected).push(row);
      });
    });
    out.sort((a, b) => b.score - a.score);
    rejected.sort((a, b) => b.score - a.score);
    /* One line per backbone, paired with the strain that suits it best.
       Six rows of the same vector against six strains is not six options,
       it is one option and a strain ranking, and the two are easier to
       read apart than together. */
    const seenVec = {}, best = [];
    out.forEach(r => {
      if (seenVec[r.vector]) { seenVec[r.vector].alternatives.push(r); return; }
      r.alternatives = [];
      seenVec[r.vector] = r;
      best.push(r);
    });
    return { inventory: inv, combinations: best, everyPairing: out, rejected: rejected,
             strategy: strat,
             counts: { considered: out.length + rejected.length, valid: out.length,
                       backbones: best.length } };
  };

  if (typeof module !== 'undefined' && module.exports) module.exports = E;
  root.Bio = E;
})(typeof window !== 'undefined' ? window : globalThis);

