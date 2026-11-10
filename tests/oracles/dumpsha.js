/* PlasmidBench's SHA-256 over a set of inputs, for cmpsha.py to check
   against Python's hashlib. The function is in the engine because
   crypto.subtle does not exist on a file:// page, and a checksum that is
   missing on a lab machine is no checksum at all — so it needs an oracle. */
const B = require('../../src/engine7.js');
const fs = require('fs');

const inputs = [
  '', 'a', 'abc', 'message digest',
  'abcdefghijklmnopqrstuvwxyz',
  'The quick brown fox jumps over the lazy dog',
  'ATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGTAACGGTGCGGGCTGA',
  'x'.repeat(55),   // one byte short of a padding block
  'x'.repeat(56),   // exactly the boundary
  'x'.repeat(57),   // one over, so a second block is needed
  'x'.repeat(63), 'x'.repeat(64), 'x'.repeat(65),
  'x'.repeat(1000), 'x'.repeat(100000),
  'héllo wörld 🦬',          // multi-byte and a surrogate pair
  '\u0000\u0001\u007f\u0080߿ࠀ￿'   // the UTF-8 length boundaries
];

const rows = inputs.map(s => ({
  /* long inputs are described rather than carried, and rebuilt the same way
     on the Python side */
  repeat: s.length > 200 ? { ch: s[0], n: s.length } : null,
  text: s.length > 200 ? null : s,
  sha: B.sha256(s)
}));

/* and the sequence checksum, which normalises before hashing */
rows.push({ label: 'seqChecksum linear', text: 'atg aaa\n12 cgc', normalisedFrom: true,
            sha: B.seqChecksum('atg aaa\n12 cgc', { topology: 'linear' }) });
rows.push({ label: 'seqChecksum circular', text: 'atg aaa\n12 cgc', normalisedFrom: true,
            sha: B.seqChecksum('atg aaa\n12 cgc', { topology: 'circular' }) });
rows.push({ label: 'seqChecksum U as T', text: 'AUGAAACGC', normalisedFrom: true,
            sha: B.seqChecksum('AUGAAACGC', { topology: 'linear' }) });

fs.writeFileSync(process.argv[2] || '/dev/stdout', JSON.stringify(rows));
