/* ===========================================================================
   Dump PlasmidBench's nucleic-acid molecular weights for cmpmw.py to set
   beside Bio.SeqUtils.molecular_weight.

   Weight is where a convention hides: whether the ends carry a phosphate,
   whether a circle loses two waters, whether U is weighed as U or as T. The
   comparison names the convention rather than assuming the two agree.
   ========================================================================= */
const B = require('../../src/engine7.js');

const SEQS = [
  'A', 'C', 'G', 'T',
  'ACGT',
  'ACGTACGTACGTAA',
  'GAATTCGGGAAGCTTCCCGGGTACCAAGCTT',
  'ATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGTAACGGTGCGGGCTGA',
  'AAAAAAAAAAAAAAAAAAAA',
  'GCGCGCGCGCGCGCGCGCGC',
  'N', 'ACGTN',
  /* a real plasmid-sized string, built deterministically so the value is
     reproducible from this file alone */
  (function () { let x = 7, s = ''; for (let i = 0; i < 3000; i++) { x = (x * 1103515245 + 12345) & 0x7fffffff; s += 'ACGT'[x % 4]; } return s; })()
];

const out = SEQS.map(s => ({
  seq: s.length > 60 ? null : s,
  length: s.length,
  sha: require('crypto').createHash('sha1').update(s).digest('hex').slice(0, 16),
  full: s,
  ss: B.ssDNAWeight(s),
  ds: B.dsDNAWeight(s)
}));
process.stdout.write(JSON.stringify(out));
