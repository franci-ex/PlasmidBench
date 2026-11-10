const B = require('../src/engine.js');
require('../src/engine2.js');
let pass = 0, fail = 0;
function ok(name, cond, extra) {
  if (cond) { pass++; console.log('  PASS  ' + name); }
  else { fail++; console.log('  FAIL  ' + name + (extra !== undefined ? '   -> ' + extra : '')); }
}
function prng(seed) { let x = seed; return () => (x = (x * 1103515245 + 12345) & 0x7fffffff) / 0x7fffffff; }
function randDNA(n, seed) { const r = prng(seed); let s = ''; for (let i = 0; i < n; i++) s += 'ACGT'[Math.floor(r() * 4)]; return s; }

/* Build a structurally valid SnapGene .dna file (packet stream). */
function packet(type, payload) {
  const head = Buffer.alloc(5);
  head.writeUInt8(type, 0);
  head.writeUInt32BE(payload.length, 1);
  return Buffer.concat([head, payload]);
}
function buildDna(seq, circular, featuresXml, primersXml, notesXml) {
  const header = Buffer.alloc(14);
  header.write('SnapGene', 0, 'latin1');
  header.writeUInt16BE(1, 8);    // f_type = DNA
  header.writeUInt16BE(14, 10);  // export version
  header.writeUInt16BE(14, 12);  // import version
  const dna = Buffer.concat([Buffer.from([circular ? 1 : 0]), Buffer.from(seq, 'latin1')]);
  const parts = [packet(9, header), packet(0, dna)];
  if (notesXml) parts.push(packet(6, Buffer.from(notesXml, 'utf8')));
  if (primersXml) parts.push(packet(5, Buffer.from(primersXml, 'utf8')));
  if (featuresXml) parts.push(packet(10, Buffer.from(featuresXml, 'utf8')));
  return Buffer.concat(parts);
}
const toAB = b => b.buffer.slice(b.byteOffset, b.byteOffset + b.byteLength);

console.log('\n== SnapGene .dna reader ==');
const seq = randDNA(3000, 901);
const featuresXml =
  '<?xml version="1.0"?>\n<Features nextValidID="5">' +
  '<Feature recentID="0" name="AmpR" type="CDS" directionality="1">' +
    '<Segment range="101-961" color="#ccffcc" type="standard"/>' +
    '<Q name="codon_start"><V int="1"/></Q>' +
    '<Q name="label"><V text="AmpR (bla)"/></Q>' +
  '</Feature>' +
  '<Feature recentID="1" name="T7 promoter" type="promoter" directionality="1">' +
    '<Segment range="20-37" color="#ffffff" type="standard"/>' +
  '</Feature>' +
  '<Feature recentID="2" name="rev gene" type="CDS" directionality="2">' +
    '<Segment range="1200-1700" color="#ff9ccd" type="standard"/>' +
  '</Feature>' +
  '<Feature recentID="3" name="amp &amp; kan &lt;dual&gt;" type="misc_feature" directionality="1">' +
    '<Segment range="2000-2100" color="#b4abac" type="standard"/>' +
  '</Feature>' +
  '<Feature recentID="4" name="spliced" type="CDS" directionality="1">' +
    '<Segment range="2400-2500" color="#cc99b2" type="standard"/>' +
    '<Segment range="2600-2700" color="#cc99b2" type="standard"/>' +
  '</Feature>' +
  '</Features>';
const primersXml = '<Primers><Primer name="M13 fwd" sequence="GTAAAACGACGGCCAGT"><BindingSite location="5-21"/></Primer>' +
  '<Primer name="seq_rev" sequence="CAGGAAACAGCTATGAC"/></Primers>';
const notesXml = '<Notes><Type>Synthetic</Type><CustomMapLabel>pTest-SG</CustomMapLabel></Notes>';

const dna = buildDna(seq, true, featuresXml, primersXml, notesXml);
let p;
try { p = B.parseSnapGene(toAB(dna)); } catch (e) { p = { error: e.message }; }
ok('parses without error', !p.error, p.error);
if (!p.error) {
  ok('sequence recovered exactly', p.seq === seq, p.seq.slice(0, 30));
  ok('length correct', p.seq.length === 3000, p.seq.length);
  ok('circular flag read', p.topology === 'circular', p.topology);
  /* Five <Feature> elements, plus the one <Primer> that carries a
     <BindingSite> — SnapGene shows those as primer_bind features and so does
     this reader, which is why the count is six, not five. */
  ok('all 5 features parsed, plus the primer binding site', p.features.length === 6, p.features.length);
  ok('and the binding site is the primer_bind one',
    p.features.filter(f => f.type === 'primer_bind').map(f => f.label).join(',') === 'M13 fwd',
    p.features.filter(f => f.type === 'primer_bind').map(f => f.label).join(','));
  const amp = p.features.find(f => /AmpR/.test(f.label));
  ok('label qualifier overrides the name attribute', amp && amp.label === 'AmpR (bla)', amp && amp.label);
  ok('feature start converted to 0-based', amp && amp.start === 100, amp && amp.start);
  ok('feature end is inclusive-to-exclusive', amp && amp.end === 961, amp && amp.end);
  ok('feature length is 861 bp', amp && (amp.end - amp.start) === 861, amp && (amp.end - amp.start));
  ok('feature type preserved', amp && amp.type === 'CDS', amp && amp.type);
  ok('SnapGene colour captured', amp && amp.color === '#ccffcc', amp && amp.color);
  const rev = p.features.find(f => f.label === 'rev gene');
  ok('directionality 2 becomes the reverse strand', rev && rev.strand === -1, rev && rev.strand);
  const fwd = p.features.find(f => f.label === 'T7 promoter');
  ok('directionality 1 stays forward', fwd && fwd.strand === 1, fwd && fwd.strand);
  ok('promoter coords correct', fwd && fwd.start === 19 && fwd.end === 37, fwd && (fwd.start + '-' + fwd.end));
  const ent = p.features.find(f => /dual/.test(f.label));
  ok('XML entities decoded in labels', ent && ent.label === 'amp & kan <dual>', ent && ent.label);
  const spl = p.features.find(f => f.label === 'spliced');
  ok('multi-segment feature spans min..max', spl && spl.start === 2399 && spl.end === 2700,
     spl && (spl.start + '-' + spl.end));
  ok('multi-segment count reported', spl && spl.segments === 2, spl && spl.segments);
  ok('primers extracted', p.primers.length === 2, p.primers.length);
  ok('primer sequence correct', p.primers[0].seq === 'GTAAAACGACGGCCAGT', p.primers[0].seq);
  ok('primer name correct', p.primers[0].name === 'M13 fwd', p.primers[0].name);
  ok('map label used as the name', p.name === 'pTest-SG', p.name);
  // features must land on the right bases
  ok('AmpR slice matches the source sequence',
     seq.slice(amp.start, amp.end) === seq.substr(100, 861));
}

// linear file
const lin = B.parseSnapGene(toAB(buildDna(randDNA(500, 902), false, '<Features></Features>', '', '')));
ok('linear flag read', lin.topology === 'linear', lin.topology);
ok('empty feature set is fine', lin.features.length === 0, lin.features.length);

// rejection paths
let bad = false, msg = '';
try { B.parseSnapGene(toAB(Buffer.from('>not a dna file\nACGTACGT'))); } catch (e) { bad = true; msg = e.message; }
ok('a FASTA file is rejected with a clear message', bad && /SnapGene/.test(msg), msg);
let noSeq = false;
try {
  const h = Buffer.alloc(14); h.write('SnapGene', 0, 'latin1');
  B.parseSnapGene(toAB(packet(9, h)));
} catch (e) { noSeq = /No DNA sequence/.test(e.message); }
ok('a header with no sequence is rejected', noSeq);
let tiny = false;
try { B.parseSnapGene(toAB(Buffer.from([1, 2, 3]))); } catch (e) { tiny = true; }
ok('a truncated file is rejected', tiny);
// a truncated packet length must not hang or throw unexpectedly
const trunc = Buffer.concat([dna.subarray(0, dna.length - 400)]);
let okTrunc = true, tp = null;
try { tp = B.parseSnapGene(toAB(trunc)); } catch (e) { okTrunc = /SnapGene|No DNA/.test(e.message); }
ok('a truncated file degrades gracefully', okTrunc, tp ? ('recovered ' + tp.seq.length + ' bp') : 'threw cleanly');

console.log('\n== round trip into a record ==');
const rec = {
  name: p.name, seq: p.seq, topology: p.topology,
  features: p.features
};
ok('every feature is inside the sequence',
   rec.features.every(f => f.start >= 0 && f.end <= rec.seq.length && f.end > f.start));
const dig = B.digest(rec.seq, ['EcoRI'], { circular: rec.topology === 'circular' });
ok('an imported record digests normally', dig.fragments.length >= 1, dig.fragments.length);
ok('digest fragments sum to the plasmid',
   dig.fragments.reduce((a, f) => a + f.length, 0) === rec.seq.length);

console.log('\n==========================================');
console.log('   ' + pass + ' passed, ' + fail + ' failed');
console.log('==========================================\n');
process.exit(fail ? 1 : 0);
