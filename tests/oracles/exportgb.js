/* ===========================================================================
   Write every record in a set of files back out as GenBank, for rtcheck.py
   to read with Biopython and compare against what went in.

   An export is only worth anything if something else can read it. This half
   writes; the Python half reads the written file with Bio.SeqIO and checks
   it against the original file, also read with Bio.SeqIO — so the engine's
   reader is not standing behind its own writer.

   Usage: node exportgb.js <output dir> <file> [<file> ...]
   ========================================================================= */
const B = require('../../src/engine7.js');
const fs = require('fs');
const path = require('path');

const outDir = process.argv[2];
fs.mkdirSync(outDir, { recursive: true });

let written = 0, skipped = 0;
for (const file of process.argv.slice(3)) {
  const name = path.basename(file);
  let recs = [];
  try {
    if (/\.dna$/i.test(file)) {
      const nb = fs.readFileSync(file);
      const r = B.parseSnapGene(nb.buffer.slice(nb.byteOffset, nb.byteOffset + nb.byteLength));
      recs = Array.isArray(r) ? r : [r];
    } else {
      const text = fs.readFileSync(file, 'utf8');
      recs = /^LOCUS/m.test(text) ? B.parseGenBank(text)
           : /^ID {3}/m.test(text) ? B.parseEMBL(text)
           : B.parseFASTA(text);
    }
  } catch (e) { skipped++; continue; }

  recs.forEach((r, i) => {
    /* a record with no bases has nothing to round-trip, and a protein one is
       not DNA GenBank; both are counted rather than quietly dropped */
    if (!r.seq || !r.seq.length || r.molecule === 'protein') { skipped++; return; }
    const gb = B.toGenBank(r);
    fs.writeFileSync(path.join(outDir, name + '.' + i + '.gb'), gb);
    written++;
  });
}
console.log(written + ' written, ' + skipped + ' had nothing to write');
