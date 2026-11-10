/* ===========================================================================
   Dump PlasmidBench's protein descriptors for every protein in Biopython's
   own test corpus, for cmpprops.py to set beside Bio.SeqUtils.ProtParam.

   The descriptors compared are the ones with one published definition and so
   one right answer: molecular weight (average isotopic), isoelectric point
   (Bjellqvist), GRAVY (Kyte & Doolittle), aromaticity (Lobry), the
   instability index (Guruprasad) and the aliphatic index (Ikai). Anything
   this project defines for itself is not in the list, because an oracle can
   only check a shared definition.
   ========================================================================= */
const B = require('../../src/engine7.js');
const fs = require('fs');
const path = require('path');

const outPath = process.argv[2] || '/dev/stdout';
const out = {};
for (const file of process.argv.slice(3)) {
  let recs = [];
  const text = fs.readFileSync(file, 'utf8');
  try {
    recs = /^LOCUS/m.test(text) ? B.parseGenBank(text)
         : /^ID {3}/m.test(text) ? B.parseEMBL(text)
         : B.parseFASTA(text);
  } catch (e) { continue; }
  for (const r of recs) {
    if (r.molecule !== 'protein') continue;
    const aa = (r.seq || '').toUpperCase().replace(/[^A-Z]/g, '');
    if (aa.length < 2) continue;
    const p = B.proteinProps(aa);
    out[path.basename(file) + '|' + r.name] = {
      seq: aa,
      length: p.length,
      mw: p.mw,
      pI: p.pI,
      gravy: p.gravy,
      aromaticity: p.aromaticity,
      instability: p.instability,
      aliphatic: p.aliphatic
    };
  }
}

/* A handful of sequences chosen to exercise the edges rather than to be
   typical: one residue, every residue once, an unknown X, a selenocysteine,
   runs of the charged residues the pI depends on. */
const EDGES = {
  'edge|single-A': 'A',
  'edge|one-of-each': 'ACDEFGHIKLMNPQRSTVWY',
  'edge|all-acidic': 'DDDDEEEEDDDDEEEE',
  'edge|all-basic': 'KKKKRRRRHHHHKKKK',
  'edge|no-charge': 'GGGGAAAAVVVVLLLL',
  'edge|cys-rich': 'CCCCCWWWWYYYY',
  'edge|with-X': 'MKVXLAGGXW',
  'edge|insulin-A': 'GIVEQCCTSICSLYQLENYCN'
};
for (const k in EDGES) {
  const aa = EDGES[k];
  const p = B.proteinProps(aa);
  out[k] = { seq: aa, length: p.length, mw: p.mw, pI: p.pI, gravy: p.gravy,
             aromaticity: p.aromaticity, instability: p.instability, aliphatic: p.aliphatic };
}

fs.writeFileSync(outPath, JSON.stringify(out));
