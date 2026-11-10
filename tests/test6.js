// Recoding, input hygiene and mutagenesis edge cases.
const B = require('../src/engine.js'); require('../src/engine2.js');
let pass = 0, fail = 0;
function ok(name, cond, extra) {
  if (cond) { pass++; console.log('  PASS  ' + name); }
  else { fail++; console.log('  FAIL  ' + name + (extra !== undefined ? '   -> ' + extra : '')); }
}
const EGFP = 'ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAA';
const prot = B.translate(EGFP, { table: 1 });

console.log('\n== codon optimisation ==');
for (const [host, U] of [['E. coli', B.ECOLI_USAGE], ['human', B.HUMAN_USAGE]]) {
  const o = B.codonOptimize(EGFP, U, 1);
  ok(host + ': protein identical', B.translate(o, { table: 1 }) === prot);
  ok(host + ': same length', o.length === EGFP.length);
}
const avoid = ['BsaI', 'BsmBI', 'EcoRI', 'BamHI', 'XhoI', 'NotI', 'HindIII', 'NcoI', 'NdeI'];
const count = (s, n) => B.findSites(s, B.enzymeByName(n), false).length;
for (const [host, U] of [['E. coli', B.ECOLI_USAGE], ['human', B.HUMAN_USAGE]]) {
  const o = B.codonOptimize(EGFP, U, 1, { avoid });
  ok(host + ' with site avoidance: protein identical', B.translate(o, { table: 1 }) === prot);
  const created = avoid.filter(n => count(o, n) > count(EGFP, n));
  ok(host + ' with site avoidance: no new site for ' + avoid.length + ' enzymes', !created.length, created.join(','));
}
/* Glu-Phe used to stand here, on the assumption that greedy E. coli recoding
   spells GAA TTC and so makes an EcoRI site. It does not: E. coli prefers TTT
   over TTC for Phe, so the greedy answer is GAA TTT and no site appears. The
   dipeptide that really does trip over itself is Asp-Pro — greedy gives
   GAT CCG, and the junction CCG|GAT reads GGATCC, a BamHI site. The greedy
   codons this engine picks for E. coli (GCG A, GAA E, TTT F, CTG L, AAA K,
   CGT R …) are the standard high-usage ones. */
const ef = 'GACCCAGACCCAGACCCA';          // (Asp-Pro) x3
const plain = B.codonOptimize(ef, B.ECOLI_USAGE, 1);
const safe = B.codonOptimize(ef, B.ECOLI_USAGE, 1, { avoid: ['BamHI'] });
ok('greedy recoding of Asp-Pro makes BamHI sites at the codon junctions',
   count(plain, 'BamHI') > 0, plain);
ok('…and avoidance removes them all', count(safe, 'BamHI') === 0, safe);
ok('…keeping the protein', B.translate(safe, { table: 1 }) === B.translate(ef, { table: 1 }));
ok('…and Glu-Phe is left alone, because greedy E. coli Phe is TTT not TTC',
   count(B.codonOptimize('GAGTTTGAGTTTGAGTTT', B.ECOLI_USAGE, 1), 'EcoRI') === 0,
   B.codonOptimize('GAGTTTGAGTTTGAGTTT', B.ECOLI_USAGE, 1));
ok('trailing partial codon is kept', B.codonOptimize('ATGAAAGG', B.ECOLI_USAGE, 1).length === 8);
ok('keepFirst leaves an alternative start alone',
   B.codonOptimize('GTGAAACTG', B.ECOLI_USAGE, 11, { keepFirst: true }).slice(0, 3) === 'GTG');
// pre-existing sites are not the optimiser's to remove
const withSite = 'ATGGAATTCAAA';
ok('an input site is left as it was', count(B.codonOptimize(withSite, B.ECOLI_USAGE, 1, { avoid: ['EcoRI'] }), 'EcoRI') <= 1);

console.log('\n== sequence text hygiene ==');
let r = B.seqTextReport('ACGT ACGT 1234\nacgtn');
ok('layout characters are not reported as dropped', r.dropped === 0 && r.seq === 'ACGTACGTACGTN');
r = B.seqTextReport('MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKL');
ok('a protein is recognised', r.looksProtein);
r = B.seqTextReport('ACGTACGTXXACGT*');
ok('stray characters are counted', r.dropped === 3 && r.droppedChars === 'X*', JSON.stringify(r));
ok('RNA converted to DNA', B.seqTextReport('AUGGCCUAA').seq === 'ATGGCCTAA');
ok('mostly-N DNA is still DNA', !B.seqTextReport('NNNNNNNNNNNNACGT').looksProtein);
ok('IUPAC-rich primer is still DNA', !B.seqTextReport('ACGTRYKMSWBDHVNACGT').looksProtein);
let recs = B.parseFASTA('>p1 protein\nMSKGEELFTGVVPILVELDGDVNG\n>d1\nATGCATGCATGC\n');
ok('FASTA: protein record flagged, DNA record not', recs[0].molecule === 'protein' && recs[1].molecule === 'DNA');
recs = B.parseGenBank('LOCUS       P1    20 aa    linear   PRT 01-JAN-2000\nORIGIN\n        1 mskgeelftg vvpilvelpq\n//\n');
ok('GenBank protein record flagged', recs[0].molecule === 'protein');

console.log('\n== site-directed mutagenesis ==');
const T = 'GATCGGATCCAAGCTTGCATGCCTGCAGGTCGACTCTAGAGGATCCCCGGGTACCGAGCTCGAATTCACTGGCCG';
ok('substitution past the end is refused', !B.sdm(T, { position: T.length - 1, kind: 'substitute', newBases: 'AAA' }).ok);
ok('deletion past the end is refused', !B.sdm(T, { position: T.length - 2, kind: 'delete', deleteLength: 5 }).ok);
ok('empty substitution is refused', !B.sdm(T, { position: 5, kind: 'substitute', newBases: '' }).ok);
const c = B.sdm(T, { position: 2, kind: 'substitute', newBases: 'A', flank: 15, circular: true });
ok('on a plasmid the primer runs across the origin at full length', c.ok && c.primers.forward.length === 31, c.primers && c.primers.forward.length);
ok('…and is the mutated sequence read round the origin',
   (c.mutated.seq + c.mutated.seq).indexOf(c.primers.forward.seq) !== -1);
const l = B.sdm(T, { position: 2, kind: 'substitute', newBases: 'A', flank: 15 });
ok('on a linear template it is clipped at the end', l.ok && l.primers.forward.length === 18);
ok('reverse primer is the reverse complement', c.primers.reverse.seq === B.revcomp(c.primers.forward.seq));

console.log('\n==========================================');
console.log('   ' + pass + ' passed, ' + fail + ' failed');
console.log('==========================================\n');
process.exit(fail ? 1 : 0);
