/* ============================================================================
   Build the three stock fixtures.

   WHY THESE ARE GENERATED AND NOT REAL

   The suites that read them — uitest-sanger.js and uitest-express.js — test
   the importer: that it finds the sequence column in a sheet laid out the way
   a lab's own spreadsheet is laid out, that it recognises vector families,
   that it counts what it read and says what it could not parse. None of that
   depends on which oligos a particular lab owns.

   The files that used to be here were the Tong lab's actual inventory: 874
   real oligo sequences with their freezer positions, and 854 real constructs
   with their targets, designs, concentrations and boxes. That is unpublished
   research data and a map of what the lab is working on, and a public
   repository is not the place for it — whatever the fixtures folder is
   nominally for. Replacing it with sequence of the same shape loses nothing
   the tests were measuring.

   WHAT IS KEPT FROM THE REAL FILES

   Only the shape, because the shape is the thing being tested:

     - the exact header lines, including "Sequence 5'->3'" with the arrow,
       which is the header the reader failed to recognise once and imported
       874 oligos named after their freezer positions;
     - the "L3 / WTOL001 / A01" style of a three-part location, in one cell;
     - one row that cannot be parsed, so "874 read, 1 skipped" is still
       exercised;
     - real vector NAMES in the plasmid sheet (pET28a, pGEX-6P-1, …), which
       are published products and are what the family recogniser is for;
     - real strain names (DH5a, BL21(DE3), Rosetta 2, …), likewise.

   Run:  node tests/fixtures/make-stock-fixtures.js
   ============================================================================ */
const fs = require('fs');
const path = require('path');

/* mulberry32: deterministic, so the fixtures are reproducible, and well
   mixed, so two seeds do not share a cycle the way a plain LCG's do */
function stream(seed) {
  let a = seed >>> 0;
  return () => {
    a = (a + 0x6D2B79F5) >>> 0;
    let t = a;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}
const r = stream(20261110);
const pick = a => a[Math.floor(r() * a.length)];
const dna = n => { let s = ''; for (let i = 0; i < n; i++) s += 'ACGT'[Math.floor(r() * 4)]; return s; };
const gcOf = s => Math.round(100 * (s.match(/[GC]/g) || []).length / s.length);

/* invented gene-like names, so nothing points at anyone's project */
const GENES = ['ARTX', 'BLNK2', 'CRDM', 'DPLX', 'EVRN', 'FKBX3', 'GLSN', 'HRTA',
  'IVSD', 'JNKL', 'KPRB', 'LMNX', 'MRTF2', 'NDLS', 'ORVC', 'PLTK', 'QRSM',
  'RVDX', 'STLN4', 'TQRP', 'USHN', 'VRLM', 'WDLK', 'XNTP', 'YRBC', 'ZPLQ'];
const KINDS = ['F', 'R', 'seq', 'mut', 'bb_F', 'bb_R', 'chk'];
const METHODS = ['Infusion', 'Gibson', 'restriction', 'SDM', 'Gateway'];

/* Some of the oligos have to actually anneal to the vector the suites load,
   because the point of the sequencing planner is to find tubes that can read
   a chosen region — and an inventory of sequence unrelated to the construct
   gives it nothing to find. A real inventory is full of primers made against
   the lab's own vectors, so the fixture takes a few dozen stretches out of
   pRSFDuet-1, which is a published Novagen product and in this folder
   already. Everything else is invented. */
function vectorPrimers(n, span) {
  const out = [];
  try {
    const B = require('../../src/engine8.js');
    const buf = fs.readFileSync(path.join(__dirname, 'pRSFDuet-1.dna'));
    const rec = B.parseSnapGene(buf.buffer.slice(buf.byteOffset, buf.byteOffset + buf.byteLength));
    const seq = (Array.isArray(rec) ? rec[0] : rec).seq;
    const reach = Math.floor(seq.length * (span || 1));
    const step = Math.floor(reach / n);
    for (let i = 0; i < n; i++) {
      const at = 40 + i * step;
      const len = 19 + Math.floor(r() * 5);
      const fwd = seq.substr(at, len);
      if (fwd.length < 18) continue;
      out.push(r() < 0.5 ? fwd : B.revcomp(fwd));
    }
  } catch (e) {
    /* no vector file: the inventory is then entirely invented, and the
       suites that need a match will say so rather than this failing quietly */
    console.error('note: pRSFDuet-1.dna not readable, so no vector-derived oligos');
  }
  return out;
}

/* ---------------------------------------------------------- oligo sheet -- */
function primerStock() {
  const lines = ['Stock ID\tName\tWhere it is\tTm (C)\tLength (nt)\tGC (%)\tSequence 5′→′3′'
    .replace('5′→′3′', '5′→3′')];
  /* 874 rows, all of them parseable: uitest-sanger.js asserts that nothing
     in this inventory is skipped, which is the point of that assertion —
     the header "Sequence 5′3′" with the arrow was once not recognised and
     every oligo came in named after its freezer position. */
  const ROWS = 874;
  /* Vector-derived oligos are clustered in the first two thirds, so that the
     far end of the plasmid has no coverage and the planner's "nothing reads
     this stretch" path is exercised. A real inventory is lopsided the same
     way: nobody has primers spread evenly round a vector. */
  const fromVector = vectorPrimers(30, 0.66);
  let vi = 0;
  /* The universal primers every lab has a tube of. They are published
     sequences, and the T7 promoter one earns its place: a Duet vector has
     two T7 promoters, so the planner must refuse it as ambiguous rather
     than plan a read with it. That behaviour has its own assertion. */
  const UNIVERSAL = [
    ['T7 promoter', 'TAATACGACTCACTATAGGG'],
    ['T7 terminator', 'GCTAGTTATTGCTCAGCGG'],
    ['M13 forward (-20)', 'GTAAAACGACGGCCAGT'],
    ['M13 reverse', 'CAGGAAACAGCTATGAC'],
    ['pBAD forward', 'ATGCCATAGCATTTTTATCC']
  ];
  UNIVERSAL.forEach(([nm, sq], k) => {
    lines.push(['OL-U' + String(k + 1).padStart(2, '0'), nm,
      'L1 / WTOL000 / A' + String(k + 1).padStart(2, '0'),
      50 + Math.round(r() * 10), sq.length, gcOf(sq), sq].join('\t'));
  });

  for (let i = 1; i <= ROWS - UNIVERSAL.length; i++) {
    const id = 'OL-' + String(i).padStart(3, '0');
    const gene = pick(GENES);
    const kind = pick(KINDS);
    const name = gene + ':' + (r() < 0.4 ? pick(['S', 'E', 'K', 'D', 'R']) +
      (20 + Math.floor(r() * 400)) + pick(['A', 'F', 'C']) + '_' : '') + kind;
    /* a tail on roughly a third of them, which is what a real inventory of
       assembly primers looks like */
    /* roughly one in eighteen reads the vector, which is what gives the
       sequencing planner something to find */
    const anneal = (vi < fromVector.length && i % 18 === 5)
      ? fromVector[vi++]
      : dna(18 + Math.floor(r() * 7));
    const seq = (r() < 0.35 && i % 18 !== 5
      ? pick(['TTGTATTTCCAGGGC', 'ATTATGCGGCCGCAATT', 'CAGGGCGCCGGC']) : '') + anneal;
    const loc = 'L' + (1 + Math.floor(r() * 4)) + ' / WTOL' +
      String(1 + Math.floor(i / 96)).padStart(3, '0') + ' / ' +
      'ABCDEFGH'[Math.floor(r() * 8)] + String(1 + (i % 12)).padStart(2, '0');
    const tm = 52 + Math.round(r() * 14);
    lines.push([id, name, loc, tm, seq.length, gcOf(seq), seq].join('\t'));
  }
  return lines.join('\n') + '\n';
}

/* -------------------------------------------------------- plasmid sheet -- */
/* Published vector products. The family recogniser is meant to know these,
   so they are the one thing that has to be real. */
const VECTORS = ['pET28a', 'pET21a', 'pET53', 'pETDuet-1', 'pRSFDuet-1',
  'pACYCDuet-1', 'pCDFDuet-1', 'pGEX-6P-1', 'pMAL-c5X', 'pSUMO', 'pTXB1',
  'pQE30', 'pcDNA3.1', 'pFastBac1', 'pDONR221', 'pLVX', 'pUC19', 'pGro7',
  'pRARE2', 'pGB1', 'pET15b', 'pBAD33'];
const RES = ['Amp', 'Kan', 'Cm', 'Spec', 'Str', 'Tet'];
const TAGS = ['N-6his', 'C-6his', 'N-6his, TEV', 'N-GST', 'N-MBP', 'N-SUMO',
  'C-strep tag II', 'N-6his, C-strep tag II', 'none'];
const ORGS = ['Synthetic', 'H. sapiens', 'M. musculus', 'E. coli', 'S. cerevisiae'];
const TYPES = ['bacteria expression', 'mammalian expression', 'insect expression',
  'cloning intermediate', 'entry clone'];
const EXPR = ['High, soluble, LOBSTR', 'Low, inclusion bodies', 'Moderate, soluble',
  'Not tested', 'High, soluble, BL21(DE3)', 'Insoluble at 37, soluble at 18'];

function plasmidStock() {
  const head = ['Plasmid ID', 'Construct ID', 'Vector', 'Resistance', 'Tag', 'Target',
    'Domain/Design', 'Organism', 'Type', 'ng/ul', 'Location', 'Shelf', 'Box',
    'Position', 'Expression?', 'Notes'].join('\t');
  const lines = [head];
  const ROWS = 854;
  for (let i = 1; i <= ROWS; i++) {
    const gene = pick(GENES);
    const dom = r() < 0.5 ? pick(['WD40', 'kinase', 'SH3', 'BRCT', 'PH', 'full length', 'CTD', 'NTD']) : '';
    lines.push([
      'FIX-' + String(i).padStart(4, '0'),
      (r() < 0.2 ? 'UbV ' : '') + gene + (dom ? '-' + dom : '') + '.' + (1 + Math.floor(r() * 4)),
      pick(VECTORS), pick(RES), pick(TAGS), gene.toLowerCase(), dom,
      pick(ORGS), pick(TYPES),
      String(15 + Math.floor(r() * 400)),
      'LAB-' + (100 + Math.floor(r() * 150)),
      'L' + (1 + Math.floor(r() * 4)),
      'WFIX-' + String(1 + Math.floor(i / 81)).padStart(3, '0'),
      'ABCDEFGHI'[Math.floor(r() * 9)] + (1 + Math.floor(r() * 9)),
      pick(EXPR),
      r() < 0.3 ? pick(['W', 'sequenced', 'verify before use', 'glycerol stock too']) : ''
    ].join('\t'));
  }
  return lines.join('\n') + '\n';
}

/* ----------------------------------------------------------- cell sheet -- */
/* Standard, published strains: what they are is the point of the sheet. */
function cellStock() {
  const head = ['Location', 'Strain', 'Existing plasmid', 'Use Case', 'Resistance',
    'stock?', 'Competency', 'Note'].join('\t');
  const rows = [
    ['DH5a', '', 'Plasmid Replication', '-', 'Yes', 'Chemical', 'endA and recA mutation, cloning, blue/white screening'],
    ['Mach1', '', 'Plasmid Replication', '-', 'Yes', 'Chemical', 'fast growing, high copy'],
    ['Stbl3', '', 'Plasmid Replication', '-', 'Yes', 'Chemical', 'for repeat-containing plasmids'],
    ['TOP10', '', 'Plasmid Replication', '-', 'Yes', 'Chemical', 'general cloning'],
    ['BL21(DE3)', '', 'Expression', '-', 'Yes', 'Chemical', 'T7 polymerase, lon and ompT deficient'],
    ['BL21(DE3) pLysS', 'pLysS', 'Expression', 'Cm', 'Yes', 'Chemical', 'tighter control of T7 expression'],
    ['Rosetta 2(DE3)', 'pRARE2', 'Expression', 'Cm', 'Yes', 'Chemical', 'rare tRNAs for eukaryotic codons'],
    ['Rosetta 2(DE3) pLysS', 'pRARE2', 'Expression', 'Cm', 'Yes', 'Chemical', 'rare tRNAs plus T7 lysozyme'],
    ['LOBSTR RIL', 'pSC101', 'Expression', 'Str', 'Yes', 'Chemical', 'reduced background binding on Ni resin'],
    ['LOBSTR BL21', 'pSC101', 'Expression', 'Str', 'Yes', 'Chemical', 'ArnA and SlyD knocked down'],
    ['Origami 2(DE3)', '', 'Expression', 'Str', 'No', 'Chemical', 'trxB gor mutant, for disulphide bonds'],
    ['SHuffle T7', '', 'Expression', 'Spec', 'Yes', 'Chemical', 'cytoplasmic disulphide bond formation'],
    ['ArcticExpress(DE3)', '', 'Expression', 'Gent', 'No', 'Chemical', 'cold-adapted chaperonins, 10-13 C'],
    ['C41(DE3)', '', 'Expression', '-', 'Yes', 'Chemical', 'tolerates toxic and membrane proteins'],
    ['C43(DE3)', '', 'Expression', '-', 'No', 'Chemical', 'tolerates toxic and membrane proteins'],
    ['BL21-AI', '', 'Expression', '-', 'No', 'Chemical', 'arabinose-inducible T7, tight control'],
    ['Tuner(DE3)', '', 'Expression', '-', 'No', 'Chemical', 'lacY mutant, IPTG titratable'],
    ['BL21(DE3) pGro7', 'pGro7', 'Expression', 'Cm', 'Yes', 'Chemical', 'GroEL/GroES co-expression'],
    ['BL21(DE3) pKJE7', 'pKJE7', 'Expression', 'Cm', 'No', 'Chemical', 'DnaK/DnaJ/GrpE co-expression'],
    ['BL21 Star(DE3)', '', 'Expression', '-', 'No', 'Chemical', 'rne131, more stable mRNA'],
    ['NiCo21(DE3)', '', 'Expression', '-', 'No', 'Chemical', 'tagged contaminants removable on Ni'],
    ['ER2566', '', 'Expression', '-', 'No', 'Chemical', 'for intein fusions'],
    ['DH10Bac', '', 'Baculovirus', 'Kan, Tet', 'Yes', 'Chemical', 'bacmid generation'],
    ['CcdB Survival', '', 'Gateway', 'Cm', 'No', 'Chemical', 'propagating ccdB-containing donor vectors'],
    ['BW25113', '', 'Other', '-', 'No', 'Electro', 'Keio parent strain'],
    ['XL10-Gold', '', 'Plasmid Replication', 'Tet', 'No', 'Chemical', 'for difficult ligations'],
    ['NEB 10-beta', '', 'Plasmid Replication', '-', 'Yes', 'Electro', 'for large plasmids and assemblies'],
    ['SoluBL21', '', 'Expression', '-', 'No', 'Chemical', 'selected for improved solubility']
  ];
  const lines = [head];
  rows.forEach((row, i) => {
    const loc = (200 + Math.floor(r() * 30)) + ', -80 L' + (1 + (i % 4));
    lines.push([loc].concat(row).join('\t'));
  });
  return lines.join('\n') + '\n';
}

const here = __dirname;
fs.writeFileSync(path.join(here, 'Primer_stock.tsv'), primerStock());
fs.writeFileSync(path.join(here, 'Plasmid_stock.tsv'), plasmidStock());
fs.writeFileSync(path.join(here, 'Competent_cells.tsv'), cellStock());
console.log('wrote Primer_stock.tsv (874 oligos, 5 of them universal), ' +
            'Plasmid_stock.tsv (854 constructs), Competent_cells.tsv (28 strains)');
