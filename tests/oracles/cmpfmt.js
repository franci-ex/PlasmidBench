/* ===========================================================================
   The new readers against an independent oracle.

   mkfmt.py takes real GenBank files, parses them with Biopython, and writes
   each record out four more times: as GFF3, as a Benchling export, as an
   openVectorEditor file and as INSDSeq XML. Alongside each it writes the
   bases Biopython's own extract() gives for every feature.

   This file reads those four back with PlasmidBench and compares, feature by
   feature, the bases PlasmidBench would hand to a primer designer against
   the bases Biopython extracted from the original. Nothing is compared
   against PlasmidBench's own GenBank reader, so an error shared by both
   would still show up here.

   A reader that is one base out anywhere fails. Which is the point: an
   off-by-one in a feature boundary is invisible on a map and fatal in a
   fusion junction.

     node cmpfmt.js <dir written by mkfmt.py>
   ========================================================================= */
const B = require('../../src/engine5.js');
const fs = require('fs'), path = require('path');

const dir = process.argv[2];
if (!dir) { console.error('usage: node cmpfmt.js <dir>'); process.exit(2); }

let files = 0, checked = 0, bad = 0;
const fails = [];

function featSeq(rec, f) {
  return B.featureSequence(rec.seq, Object.assign({}, f, { strand: f.strand }));
}

/* match up by label, because the formats do not promise an order */
function compare(tag, file, recs, truth, which) {
  if (!recs || !recs.length) { fails.push(tag + ' ' + file + ': nothing parsed'); bad++; return; }
  const rec = recs[0];
  if (rec.seq !== truth.seq) { fails.push(tag + ' ' + file + ': sequence differs'); bad++; return; }
  const want = truth[which];
  const got = {};
  rec.features.forEach(f => { (got[f.label] = got[f.label] || []).push(featSeq(rec, f)); });
  want.forEach(w => {
    checked++;
    const list = got[w.label];
    if (!list || list.indexOf(w.seq) === -1) {
      bad++;
      fails.push(tag + ' ' + file + ' [' + w.label + ' ' + w.type + ']: ' +
        (list ? 'got ' + (list[0] || '').slice(0, 30) + '… want ' + w.seq.slice(0, 30) + '…'
              : 'feature missing'));
    }
  });
  /* circularity must survive the trip, or a wrap-round primer search is wrong */
  const wantTopo = truth.circular ? 'circular' : 'linear';
  if (rec.topology !== wantTopo) { bad++; fails.push(tag + ' ' + file + ': topology ' + rec.topology + ' want ' + wantTopo); }
}

fs.readdirSync(dir).filter(f => f.endsWith('.truth.json')).sort().forEach(tf => {
  const base = tf.replace(/\.truth\.json$/, '');
  const truth = JSON.parse(fs.readFileSync(path.join(dir, tf), 'utf8'));
  files++;
  const read = (ext, fn) => {
    const p = path.join(dir, base + ext);
    if (!fs.existsSync(p)) return null;
    try { return fn(fs.readFileSync(p, 'utf8')); }
    catch (e) { bad++; fails.push(ext + ' ' + base + ': threw ' + e.message); return null; }
  };
  compare('gff3', base, read('.gff3', t => B.parseGFF(t)), truth, 'simple');
  compare('benchling', base, read('.benchling.json', t => B.parseSeqJSON(t)), truth, 'simple');
  compare('ove', base, read('.ove.json', t => B.parseSeqJSON(t)), truth, 'simple');
  compare('insdxml', base, read('.insd.xml', t => B.parseSeqXML(t)), truth, 'all');

  /* the sniffer has to pick the right door with no help from the name */
  [['.gff3', 'gff'], ['.benchling.json', 'json'], ['.ove.json', 'json'], ['.insd.xml', 'xml']].forEach(([ext, want]) => {
    const p = path.join(dir, base + ext);
    if (!fs.existsSync(p)) return;
    checked++;
    const got = B.sniffFormat(fs.readFileSync(p, 'utf8'), 'unnamed');
    if (got !== want) { bad++; fails.push('sniff ' + base + ext + ': ' + got + ' want ' + want); }
  });
});

fails.slice(0, 20).forEach(f => console.log('  FAIL ' + f));
console.log('  ' + files + ' record(s) x 4 formats: ' + checked + ' checks, ' + bad + ' differ from Biopython');
process.exit(bad ? 1 : 0);
