/* Dump the embedded BLOSUM62 so Biopython can be asked whether it is right. */
const B = require('../../src/engine6.js');
const alpha = 'ARNDCQEGHILKMFPSTWYVBZX*';
const out = {};
for (const a of alpha) for (const b of alpha) out[a + b] = B.blosum62(a, b);
process.stdout.write(JSON.stringify({ alphabet: alpha, matrix: out }));
