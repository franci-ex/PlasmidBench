/* ===========================================================================
   Dump what PlasmidBench makes of real binary files — SnapGene .dna and ABI
   .ab1 traces — for cmpreal.py to set beside Biopython's answer.

   Binary formats are where a reader is most likely to be quietly wrong: a
   wrong offset gives a sequence that still looks like DNA. Only an
   independent implementation of the same format settles it.
   ========================================================================= */
const B = require('../../src/engine7.js');
const fs = require('fs');
const path = require('path');

const out = {};
for (const file of process.argv.slice(2)) {
  const name = path.basename(file);
  /* the browser hands these readers an ArrayBuffer, so a Node Buffer is
     turned into exactly that rather than into a view over the pool */
  const nb = fs.readFileSync(file);
  const buf = nb.buffer.slice(nb.byteOffset, nb.byteOffset + nb.byteLength);
  try {
    if (/\.dna$/i.test(file)) {
      const r = B.parseSnapGene(buf);
      const rec = Array.isArray(r) ? r[0] : r;
      out[name] = {
        kind: 'dna',
        seq: (rec.seq || '').toUpperCase(),
        length: (rec.seq || '').length,
        topology: rec.topology,
        features: (rec.features || []).map(f => ({
          type: f.type, start: f.start, end: f.end,
          strand: f.strand, label: f.label
        }))
      };
    } else {
      const r = B.parseAB1(buf);
      const rec = Array.isArray(r) ? r[0] : r;
      const seq = rec.sequence || rec.seq || '';
      out[name] = {
        kind: 'ab1',
        seq: seq.toUpperCase(),
        length: seq.length,
        /* the per-base confidences, which is what trimming decisions rest on */
        qual: rec.quality || rec.qual || null,
        sampleName: rec.sampleName || null
      };
    }
  } catch (e) {
    out[name] = { error: String(e && e.message || e) };
  }
}
process.stdout.write(JSON.stringify(out));
