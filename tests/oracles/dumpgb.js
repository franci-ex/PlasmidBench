/* ===========================================================================
   Dump what PlasmidBench's GenBank reader makes of a set of files, in a form
   cmpgb.py can set beside Biopython's answer for the same files.

   The point of the pair is that neither half knows what the other expects.
   This side only reports; the comparison and every expected value live in
   the Python half, which reads the same files with Bio.SeqIO and never looks
   at this code.

   Positions are reported exactly as the engine holds them: 0-based, end
   exclusive, strand +1 or -1, a joined feature as its parts in order. That
   is also how Biopython holds them, so the comparison needs no arithmetic
   on either side — arithmetic in a test is somewhere for an error to hide.
   ========================================================================= */
const B = require('../../src/engine7.js');
const fs = require('fs');
const path = require('path');
const crypto = require('crypto');

const sha = s => crypto.createHash('sha1').update(String(s).toUpperCase()).digest('hex').slice(0, 16);

/** @param {{start:number,end:number,strand:number,parts?:Array<{start:number,end:number,strand?:number}>}} f */
function loc(f) {
  const parts = (f.parts && f.parts.length) ? f.parts : [{ start: f.start, end: f.end, strand: f.strand }];
  return parts.map(p => p.start + '..' + p.end + ':' + ((p.strand == null ? f.strand : p.strand) < 0 ? '-' : '+')).join(',');
}

const out = {};
for (const file of process.argv.slice(2)) {
  let recs;
  try {
    recs = B.parseGenBank(fs.readFileSync(file, 'utf8'));
  } catch (e) {
    out[path.basename(file)] = { error: String(e && e.message || e) };
    continue;
  }
  if (!Array.isArray(recs)) recs = [recs];
  out[path.basename(file)] = recs.map(r => ({
    name: r.name,
    length: (r.seq || '').length,
    seq: (r.seq || '').toUpperCase(),
    topology: r.topology,
    molecule: r.molecule,
    topologyFromFeature: !!r.topologyFromFeature,
    skipped: r.skippedFeatures || [],
    features: (r.features || []).map(f => ({
      type: f.type,
      loc: loc(f),
      /* The bases this feature actually covers, hashed. Location order is a
         convention — this library lists the pieces of a minus-strand join in
         ascending top-strand order, Biopython in transcription order — so
         comparing the locations alone turns a convention into a failure and
         leaves the thing that matters untested. What matters is the sequence
         the feature yields, and that has one right answer. */
      extract: sha(f.between ? '' : B.featureSequence(r.seq || '', f)),
      length: f.between ? 0 : B.featureLength(f, (r.seq || '').length),
      label: f.label,
      quals: f.quals || {}
    }))
  }));
}
process.stdout.write(JSON.stringify(out));
