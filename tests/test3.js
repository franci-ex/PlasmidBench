const B = require('../src/engine.js');
require('../src/engine2.js');
let pass = 0, fail = 0;
function ok(name, cond, extra) {
  if (cond) { pass++; console.log('  PASS  ' + name); }
  else { fail++; console.log('  FAIL  ' + name + (extra !== undefined ? '   -> ' + extra : '')); }
}
function prng(seed) { let x = seed; return () => (x = (x * 1103515245 + 12345) & 0x7fffffff) / 0x7fffffff; }
function randDNA(n, seed) {
  const r = prng(seed); let s = '';
  for (let i = 0; i < n; i++) s += 'ACGT'[Math.floor(r() * 4)];
  return s;
}

/* ---- Build a synthetic but structurally valid ABIF (.ab1) file so the
       parser can be verified without a proprietary sample file. ---- */
function buildAB1(seq, quality, traces, order) {
  const fields = [];
  function addField(name, number, etype, esize, nelem, bytes) {
    fields.push({ name, number, etype, esize, nelem, bytes });
  }
  const enc = s => Buffer.from(s, 'latin1');
  addField('PBAS', 1, 2, 1, seq.length, enc(seq));
  addField('PCON', 1, 2, 1, quality.length, Buffer.from(quality.map(q => q & 0xff)));
  addField('FWO_', 1, 2, 1, 4, enc(order));              // 4 bytes -> inline
  const shorts = arr => {
    const b = Buffer.alloc(arr.length * 2);
    arr.forEach((v, i) => b.writeInt16BE(v, i * 2));
    return b;
  };
  for (let i = 0; i < 4; i++) {
    const ch = traces[order[i]];
    addField('DATA', 9 + i, 4, 2, ch.length, shorts(ch));
  }
  addField('SMPL', 1, 18, 1, 1, Buffer.concat([Buffer.from([4]), enc('cl12')]));

  const headerLen = 6 + 28;
  let dataPos = headerLen;
  for (const f of fields) {
    f.dataSize = f.bytes.length;
    if (f.dataSize > 4) { f.dataOffset = dataPos; dataPos += f.dataSize; }
    else f.dataOffset = null;   // inline
  }
  const dirOffset = dataPos;
  const total = dirOffset + fields.length * 28;
  const buf = Buffer.alloc(total, 0);
  buf.write('ABIF', 0, 'latin1');
  buf.writeInt16BE(101, 4);
  // root directory entry
  buf.write('tdir', 6, 'latin1');
  buf.writeInt32BE(1, 10);
  buf.writeInt16BE(1023, 14);
  buf.writeInt16BE(28, 16);
  buf.writeInt32BE(fields.length, 18);
  buf.writeInt32BE(fields.length * 28, 22);
  buf.writeInt32BE(dirOffset, 26);
  buf.writeInt32BE(0, 30);
  // data blobs
  for (const f of fields) if (f.dataOffset !== null) f.bytes.copy(buf, f.dataOffset);
  // directory entries
  fields.forEach((f, i) => {
    const b = dirOffset + i * 28;
    buf.write(f.name, b, 'latin1');
    buf.writeInt32BE(f.number, b + 4);
    buf.writeInt16BE(f.etype, b + 8);
    buf.writeInt16BE(f.esize, b + 10);
    buf.writeInt32BE(f.nelem, b + 12);
    buf.writeInt32BE(f.dataSize, b + 16);
    if (f.dataOffset !== null) buf.writeInt32BE(f.dataOffset, b + 20);
    else f.bytes.copy(buf, b + 20);   // inline payload
    buf.writeInt32BE(0, b + 24);
  });
  return buf;
}

console.log('\n== ABIF (.ab1) parser ==');
const readSeq = randDNA(600, 301);
const qual = [];
for (let i = 0; i < 600; i++) qual.push(i < 20 ? 8 : (i > 560 ? 9 : 55));  // poor ends
const order = 'GATC';
const traces = {};
for (const c of order) traces[c] = Array.from({ length: 1200 }, (_, i) => (i * 7 + c.charCodeAt(0)) % 900);
const ab1 = buildAB1(readSeq, qual, traces, order);
let parsed;
try { parsed = B.parseAB1(ab1.buffer.slice(ab1.byteOffset, ab1.byteOffset + ab1.byteLength)); }
catch (e) { parsed = { error: e.message }; }
ok('ABIF parses without error', !parsed.error, parsed.error);
if (!parsed.error) {
  ok('base calls recovered exactly', parsed.sequence === readSeq, parsed.sequence.slice(0, 30));
  ok('read length correct', parsed.length === 600, parsed.length);
  ok('quality array recovered', parsed.quality && parsed.quality.length === 600,
     parsed.quality && parsed.quality.length);
  ok('quality values correct', parsed.quality[100] === 55 && parsed.quality[0] === 8,
     parsed.quality && parsed.quality.slice(0, 3).join(','));
  ok('mean quality computed', parsed.meanQuality > 40 && parsed.meanQuality < 56, parsed.meanQuality);
  ok('channel order read from FWO_ (inline field)', parsed.channelOrder === 'GATC', parsed.channelOrder);
  ok('all four traces present', parsed.hasTraces && Object.keys(parsed.traces).length === 4,
     Object.keys(parsed.traces).join(','));
  ok('trace length correct', parsed.traces.G.length === 1200, parsed.traces.G.length);
  ok('trace values correct', parsed.traces.G[3] === (3 * 7 + 71) % 900, parsed.traces.G[3]);
  ok('pString sample name read', parsed.sampleName === 'cl12', parsed.sampleName);
  ok('version read', parsed.version === 101, parsed.version);
}
// a non-ABIF file must be rejected clearly, not parsed into nonsense
let rejected = false, msg = '';
try { B.parseAB1(Buffer.from('>not an ab1 file at all, this is fasta\nACGT').buffer); }
catch (e) { rejected = true; msg = e.message; }
ok('a non-ABIF file is rejected with a clear message', rejected && /ABIF/.test(msg), msg);
let tooSmall = false;
try { B.parseAB1(Buffer.from('AB').buffer); } catch (e) { tooSmall = true; }
ok('a truncated file is rejected', tooSmall);

console.log('\n== quality trimming ==');
const trim = B.trimByQuality(readSeq, qual, { window: 10, minMean: 20 });
ok('poor 5-prime end is trimmed', trim.start >= 10, trim.start);
ok('poor 3-prime end is trimmed', trim.end <= 570, trim.end);
ok('trimmed read is a substring of the original', readSeq.includes(trim.seq));
ok('trimmed length is consistent', trim.seq.length === trim.end - trim.start);
const noTrim = B.trimByQuality(readSeq, readSeq.split('').map(() => 60), { minMean: 20 });
ok('a uniformly good read is not trimmed', noTrim.start === 0 && noTrim.end === readSeq.length,
   noTrim.start + '-' + noTrim.end);
ok('missing quality returns the whole read', B.trimByQuality('ACGT', null).seq === 'ACGT');

console.log('\n== align to reference (sequencing verification) ==');
const ref = randDNA(1500, 311);
// a perfect read from the middle
const perfect = ref.slice(400, 1000);
const a1 = B.alignToReference(ref, perfect);
ok('perfect read is clean', a1.clean, JSON.stringify(a1.differences.slice(0, 3)));
ok('perfect read is 100% identity', Math.abs(a1.identity - 100) < 1e-9, a1.identity);
ok('read is placed at the right reference position', a1.refStart === 400, a1.refStart);
ok('read end is correct', a1.refEnd === 1000, a1.refEnd);
ok('not reported as flipped', a1.flipped === false);
// a reverse-complement read must be detected and flipped
const a2 = B.alignToReference(ref, B.revcomp(perfect));
ok('reverse read is detected as flipped', a2.flipped === true);
ok('flipped read is clean', a2.clean, JSON.stringify(a2.differences.slice(0, 3)));
ok('flipped read placed correctly', a2.refStart === 400, a2.refStart);
// a single point mutation
const mut = perfect.slice(0, 300) + (perfect[300] === 'A' ? 'G' : 'A') + perfect.slice(301);
const a3 = B.alignToReference(ref, mut);
ok('point mutation detected', a3.differences.length === 1, JSON.stringify(a3.differences));
ok('mutation typed as mismatch', a3.differences[0].type === 'mismatch', a3.differences[0].type);
ok('mutation at the right reference coordinate', a3.differences[0].refPos === 700,
   a3.differences[0].refPos);
ok('mutation reports both bases',
   a3.differences[0].ref === perfect[300] && a3.differences[0].read !== perfect[300],
   JSON.stringify(a3.differences[0]));
ok('mutated read is not clean', !a3.clean);
// a deletion
const delRead = perfect.slice(0, 300) + perfect.slice(303);
const a4 = B.alignToReference(ref, delRead);
ok('deletion produces deletion differences', a4.differences.some(d => d.type === 'deletion'),
   JSON.stringify(a4.differences.slice(0, 5)));
ok('deletion read is not clean', !a4.clean);
// coverage reporting
ok('coverage is reported as a percentage', a1.coverage > 39 && a1.coverage < 41, a1.coverage);

console.log('\n== end-to-end: clone, then verify by simulated sequencing ==');
// build a construct, then "sequence" it with a noisy read and confirm the
// verification workflow flags exactly the planted error.
const vector = 'G'.repeat(1500) + 'GAATTC' + 'TTTT' + 'AAGCTT' + 'C'.repeat(400);
const gene = randDNA(450, 321);
const insert = 'AA' + 'GAATTC' + gene + 'AAGCTT' + 'TT';
const built = B.restrictionClone({
  vector, insert, vectorEnzymes: ['EcoRI', 'HindIII'], insertEnzymes: ['EcoRI', 'HindIII'],
  vectorCircular: true
});
ok('construct built', built.ok, built.error);
const product = built.ok ? built.products[0].seq : '';
ok('construct contains the gene', product.includes(gene));
// a good read over the insert
const goodRead = product.slice(0, 500);
const v1 = B.alignToReference(product, goodRead);
ok('good read verifies clean', v1.clean, JSON.stringify(v1.differences.slice(0, 3)));
// a read carrying a single planted substitution
const badRead = goodRead.slice(0, 250) + (goodRead[250] === 'T' ? 'C' : 'T') + goodRead.slice(251);
const v2 = B.alignToReference(product, badRead);
ok('planted error is found', v2.differences.length === 1, JSON.stringify(v2.differences));
ok('planted error is at the right position', v2.differences[0].refPos === 250, v2.differences[0].refPos);
// and the construct's own features survive a digest check
const checkDigest = B.digest(product, ['EcoRI', 'HindIII'], { circular: true });
ok('diagnostic digest gives 2 fragments', checkDigest.fragments.length === 2,
   checkDigest.fragments.length);
ok('digest fragment sizes sum to the construct',
   checkDigest.fragments.reduce((a, f) => a + f.length, 0) === product.length,
   checkDigest.fragments.map(f => f.length).join(','));
console.log('    construct', product.length, 'bp; diagnostic digest:',
            checkDigest.fragments.map(f => f.length).sort((a, b) => b - a).join(' + '));

console.log('\n==========================================');
console.log('   ' + pass + ' passed, ' + fail + ' failed');
console.log('==========================================\n');
process.exit(fail ? 1 : 0);
