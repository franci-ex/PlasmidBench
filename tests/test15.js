/* ===========================================================================
   The file formats.

   cmpfmt.js checks the readers against Biopython on real records. This file
   checks the things an oracle cannot: the dialects a real export uses that
   the oracle does not write, the edge cases that make a reader quietly wrong
   (a wrap, a zero as an end marker, a tab where a column should be), and the
   refusals — because a format read half-wrongly is worse than one refused.
   ========================================================================= */
const B = require('../src/engine5.js');
let pass = 0, fail = 0;
function ok(name, cond, extra) {
  if (cond) { pass++; }
  else { fail++; console.log('  FAIL  ' + name + (extra !== undefined ? '   -> ' + extra : '')); }
}
function section(s) { console.log('\n== ' + s + ' =='); }

const SEQ = 'ATGGCTAGCAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGGCACAAATTTTCTGTC';
const L = SEQ.length;   /* 90 */

/* ------------------------------------------------------------------ sniff */
section('the format is read from the content, not the name');
ok('GenBank recognised with a .txt name',
  B.sniffFormat('LOCUS       x   10 bp DNA circular\nORIGIN\n  1 acgtacgtac\n//\n', 'emailed.txt') === 'genbank');
ok('GenBank recognised with an ApE .ape name',
  B.sniffFormat('LOCUS       x   10 bp ds-DNA circular 1-JAN-2020\nORIGIN\n  1 acgtacgtac\n//\n', 'p.ape') === 'genbank');
ok('FASTA recognised with no name at all',
  B.sniffFormat('>s\n' + SEQ + '\n', '') === 'fasta');
ok('GFF3 recognised from its header', B.sniffFormat('##gff-version 3\n', 'x') === 'gff');
ok('GFF3 recognised from its columns alone when the header is missing',
  B.sniffFormat('c1\tsrc\tCDS\t1\t9\t.\t+\t0\tID=a\nc1\tsrc\tgene\t2\t8\t.\t-\t.\tID=b\n', 'x') === 'gff');
ok('JSON recognised from its first character', B.sniffFormat('{"bases":"ACGT"}', 'x') === 'json');
ok('XML recognised from its declaration', B.sniffFormat('<?xml version="1.0"?><TSeqSet/>', 'x') === 'xml');
ok('bare bases recognised as raw sequence', B.sniffFormat(SEQ, 'x') === 'raw');
ok('an empty file is nothing at all', B.sniffFormat('   \n\n', 'x.gb') === null);
ok('EMBL recognised', B.sniffFormat('ID   X; SV 1; circular; DNA; STD; SYN; 10 BP.\nSQ   Sequence 10 BP;\n     acgtacgtac 10\n//\n', 'x') === 'embl');

/* --------------------------------------------------------------- GenBank */
section('GenBank dialects');
const gbTabs = 'LOCUS\tpTab\t90 bp\tDNA\tcircular\tSYN\t01-JAN-2024\n' +
  'FEATURES\tLocation/Qualifiers\n' +
  '\tCDS\t1..30\n' +
  '\t\t/label="tagged"\n' +
  'ORIGIN\n' +
  '        1 ' + SEQ.toLowerCase().slice(0, 60) + '\n' +
  '       61 ' + SEQ.toLowerCase().slice(60) + '\n//\n';
const tabbed = B.parseGenBank(gbTabs);
ok('a file written with tabs instead of columns still reads', tabbed.length === 1 && tabbed[0].seq === SEQ,
  tabbed.length && tabbed[0].seq.length);
ok('and its features land in the right place',
  tabbed[0].features.length === 1 && tabbed[0].features[0].start === 0 && tabbed[0].features[0].end === 30,
  JSON.stringify(tabbed[0].features));
ok('and its topology survives', tabbed[0].topology === 'circular');

const gbDeep = 'LOCUS       pDeep                     90 bp    DNA     linear   SYN 01-JAN-2024\n' +
  'FEATURES             Location/Qualifiers\n' +
  '        CDS             4..30\n' +
  '                        /label=deep\n' +
  'ORIGIN\n        1 ' + SEQ.toLowerCase() + '\n//\n';
const deep = B.parseGenBank(gbDeep);
ok('a feature table indented to 8 rather than 5 still reads',
  deep[0].features.length === 1 && deep[0].features[0].label === 'deep', JSON.stringify(deep[0].features));

const gbColor = 'LOCUS       pCol                      90 bp    DNA     linear   SYN 01-JAN-2024\n' +
  'FEATURES             Location/Qualifiers\n' +
  '     CDS             1..30\n' +
  '                     /label=ape\n' +
  '                     /ApEinfo_fwdcolor=#ff6600\n' +
  '     misc_feature    31..60\n' +
  '                     /label=sg\n' +
  '                     /note="color: #0066ff"\n' +
  '     misc_feature    61..90\n' +
  '                     /label=plain\n' +
  'ORIGIN\n        1 ' + SEQ.toLowerCase() + '\n//\n';
const col = B.parseGenBank(gbColor)[0];
ok("ApE's feature colour is kept", col.features[0].color === '#ff6600', col.features[0].color);
ok("SnapGene's colour-in-the-note is kept", col.features[1].color === '#0066ff', col.features[1].color);
ok('a feature with no colour does not get one', !col.features[2].color, col.features[2].color);

const gbContig = 'LOCUS       pCon                   90 bp    DNA     linear   CON 01-JAN-2024\n' +
  'CONTIG      join(AB000001.1:1..90)\n//\n';
ok('a CONTIG record is flagged rather than imported as empty', B.parseGenBank(gbContig)[0].contigOnly === true);

const gbNoEnd = 'LOCUS       pNoEnd                    90 bp    DNA     linear   SYN 01-JAN-2024\n' +
  'ORIGIN\n        1 ' + SEQ.toLowerCase() + '\n';
ok('a file with no // terminator still reads', B.parseGenBank(gbNoEnd)[0].seq === SEQ);

/* ------------------------------------------------------------------- GFF */
section('GFF3');
const gff = ['##gff-version 3',
  '##sequence-region p1 1 90',
  'p1\tt\tregion\t1\t90\t.\t+\t.\tID=p1;Is_circular=true',
  'p1\tt\tCDS\t4\t30\t.\t+\t0\tID=f1;Name=GFP%2C N-term;product=green',
  'p1\tt\tterminator\t61\t90\t.\t-\t.\tID=f2;locus_tag=T7term',
  '##FASTA', '>p1 a test', SEQ].join('\n');
const g = B.parseGFF(gff);
ok('one record comes out', g.length === 1 && g[0].name === 'p1', g.length);
ok('the bases come from the ##FASTA block', g[0].seq === SEQ);
ok('Is_circular is honoured', g[0].topology === 'circular');
ok('column 4 is 1-based and becomes 0-based', g[0].features[0].start === 3, g[0].features[0].start);
ok('column 5 is inclusive and becomes half-open', g[0].features[0].end === 30, g[0].features[0].end);
ok('the percent-encoding in an attribute is decoded', g[0].features[0].label === 'GFP, N-term', g[0].features[0].label);
ok('a minus in column 7 is the reverse strand', g[0].features[1].strand === -1);
ok('locus_tag is used when there is no Name', g[0].features[1].label === 'T7term', g[0].features[1].label);
ok('the contig row itself is not a feature', g[0].features.length === 2, g[0].features.length);
ok('the extracted bases match the sequence',
  B.featureSequence(g[0].seq, g[0].features[0]) === SEQ.slice(3, 30));

const gffOnly = ['##gff-version 3', 'p1\tt\tCDS\t4\t30\t.\t+\t0\tID=f1;Name=x'].join('\n');
const go = B.parseGFF(gffOnly);
ok('a GFF with no bases is flagged as annotations only', go.annotationsOnly === true);
ok('and its features are still read', go[0].features.length === 1, go[0].features.length);

const gffBad = ['##gff-version 3', '##sequence-region p1 1 90',
  'p1\tt\tCDS\t4\t30\t.\t+\t0\tID=ok;Name=ok',
  'p1\tt\tCDS\t80\t400\t.\t+\t0\tID=bad;Name=bad',
  '##FASTA', '>p1', SEQ].join('\n');
const gb2 = B.parseGFF(gffBad);
ok('a row running past the end of the sequence is not clamped to fit',
  gb2[0].features.length === 1 && gb2[0].features[0].label === 'ok', JSON.stringify(gb2[0].features.map(f => f.label)));
ok('and is reported rather than dropped silently',
  (gb2[0].skippedFeatures || []).length === 1, JSON.stringify(gb2[0].skippedFeatures));

/* ------------------------------------------------------------------ JSON */
section('JSON dialects');
const bench = JSON.stringify({
  name: 'pBench', bases: SEQ, isCircular: true,
  annotations: [
    { name: 'gfp', type: 'CDS', start: 3, end: 30, strand: 1 },
    { name: 'toEnd', type: 'misc_feature', start: 60, end: 0, strand: 1 },
    { name: 'wrapped', type: 'CDS', start: 80, end: 10, strand: 1 }
  ]
});
const bj = B.parseSeqJSON(bench);
ok('a Benchling export is recognised as one', bj.dialect === 'benchling', bj.dialect);
ok('its end is exclusive, so no base is lost', bj[0].features[0].end === 30, bj[0].features[0].end);
ok('and the bases match', B.featureSequence(bj[0].seq, bj[0].features[0]) === SEQ.slice(3, 30));
ok('an end of 0 means the end of the sequence, not a zero-length feature',
  bj[0].features[1].end === L, bj[0].features[1].end);
ok('an annotation that wraps the origin stays wrapped',
  bj[0].features[2].start === 80 && bj[0].features[2].end === 10, JSON.stringify(bj[0].features[2]));
ok('and its bases are the two pieces in order',
  B.featureSequence(bj[0].seq, bj[0].features[2]) === SEQ.slice(80) + SEQ.slice(0, 10),
  B.featureSequence(bj[0].seq, bj[0].features[2]));
ok('isCircular is honoured', bj[0].topology === 'circular');

const ove = JSON.stringify({
  name: 'pOve', sequence: SEQ, circular: false,
  features: [{ name: 'gfp', type: 'CDS', start: 3, end: 29, forward: true },
             { name: 'rev', type: 'CDS', start: 60, end: 89, forward: false }]
});
const oj = B.parseSeqJSON(ove);
ok('an openVectorEditor file is recognised as one', oj.dialect === 'ove', oj.dialect);
ok('its ends are inclusive, so one is added', oj[0].features[0].end === 30, oj[0].features[0].end);
ok('the same bases come out as the Benchling file gave',
  B.featureSequence(oj[0].seq, oj[0].features[0]) === SEQ.slice(3, 30));
ok('forward:false is the reverse strand', oj[0].features[1].strand === -1);

const ranged = JSON.stringify({
  name: 'pRange', sequence: SEQ,
  features: [{ name: 'a', type: 'CDS', rangeStart: 3, rangeEnd: 29, strand: 1 }]
});
const rj = B.parseSeqJSON(ranged);
ok('rangeStart/rangeEnd is read as an inclusive range', rj[0].features[0].end === 30, rj[0].features[0].end);

const joined = JSON.stringify({
  name: 'pJoin', sequence: SEQ,
  features: [{ name: 'spliced', type: 'CDS', start: 3, end: 50,
               locations: [{ start: 3, end: 11 }, { start: 30, end: 50 }], forward: true }]
});
const jj = B.parseSeqJSON(joined);
ok('a feature given as pieces keeps its pieces', (jj[0].features[0].parts || []).length === 2,
  JSON.stringify(jj[0].features[0].parts));
ok('and its bases are the pieces joined, not the span',
  B.featureSequence(jj[0].seq, jj[0].features[0]) === SEQ.slice(3, 12) + SEQ.slice(30, 51),
  B.featureSequence(jj[0].seq, jj[0].features[0]));

ok('a nested sequence object is not read as the string "[object Object]"', (() => {
  const r = B.parseSeqJSON(JSON.stringify({ name: 'n', sequence: { sequence: SEQ, name: 'inner' }, features: [] }));
  return r[0].seq === SEQ;
})());

let threw = '';
try { B.parseSeqJSON('{"notASequence":true}'); } catch (e) { threw = e.code; }
ok('JSON with no sequence in it is refused, not guessed at', threw === 'noseq', threw);
threw = '';
try { B.parseSeqJSON('{oops'); } catch (e) { threw = e.code; }
ok('JSON that does not parse says so', threw === 'badjson', threw);

/* ------------------------------------------------------------------- XML */
section('XML');
const insd = '<?xml version="1.0"?><INSDSet><INSDSeq>' +
  '<INSDSeq_locus>pXml</INSDSeq_locus>' +
  '<INSDSeq_definition>a &amp; b</INSDSeq_definition>' +
  '<INSDSeq_topology>circular</INSDSeq_topology>' +
  '<INSDSeq_feature-table>' +
  '<INSDFeature><INSDFeature_key>CDS</INSDFeature_key>' +
  '<INSDFeature_location>complement(join(4..12,31..51))</INSDFeature_location>' +
  '<INSDFeature_quals><INSDQualifier><INSDQualifier_name>gene</INSDQualifier_name>' +
  '<INSDQualifier_value>bla</INSDQualifier_value></INSDQualifier></INSDFeature_quals>' +
  '</INSDFeature></INSDSeq_feature-table>' +
  '<INSDSeq_sequence>' + SEQ.toLowerCase() + '</INSDSeq_sequence></INSDSeq></INSDSet>';
const xi = B.parseSeqXML(insd);
ok('an INSDSeq record reads', xi.length === 1 && xi[0].seq === SEQ, xi.length);
ok('its topology reads', xi[0].topology === 'circular');
ok('its entities are decoded', xi[0].description === 'a & b', xi[0].description);
ok('a joined complement location reads through the same parser GenBank uses',
  xi[0].features[0].strand === -1 && (xi[0].features[0].parts || []).length === 2,
  JSON.stringify(xi[0].features[0]));
ok('and gives the same bases as the GenBank spelling would', (() => {
  const gbv = 'LOCUS       x   90 bp DNA circular SYN 01-JAN-2024\nFEATURES             Location/Qualifiers\n' +
    '     CDS             complement(join(4..12,31..51))\n                     /gene=bla\n' +
    'ORIGIN\n        1 ' + SEQ.toLowerCase() + '\n//\n';
  const g2 = B.parseGenBank(gbv)[0];
  return B.featureSequence(g2.seq, g2.features[0]) === B.featureSequence(xi[0].seq, xi[0].features[0]);
})());
ok('the qualifier becomes the label', xi[0].features[0].label === 'bla', xi[0].features[0].label);

const tseq = '<?xml version="1.0"?><TSeqSet><TSeq><TSeq_seqtype value="nucleotide"/>' +
  '<TSeq_accver>NM_000.1</TSeq_accver><TSeq_defline>a gene</TSeq_defline>' +
  '<TSeq_sequence>' + SEQ + '</TSeq_sequence></TSeq></TSeqSet>';
const xt = B.parseSeqXML(tseq);
ok('a TinySeq record reads', xt.length === 1 && xt[0].seq === SEQ && xt[0].name === 'NM_000.1', JSON.stringify(xt[0] && xt[0].name));
threw = '';
try { B.parseSeqXML('<?xml version="1.0"?><nothing/>'); } catch (e) { threw = e.code; }
ok('XML with no sequence in it is refused', threw === 'noseq', threw);

/* ------------------------------------------------- one door for all of it */
section('parseAnySequence');
[['genbank', gbColor, 'GenBank'], ['gff', gff, 'GFF3'], ['xml', insd, 'XML'],
 ['fasta', '>x\n' + SEQ + '\n', 'FASTA']].forEach(([kind, text, want]) => {
  const r = B.parseAnySequence(text, 'file');
  ok(kind + ' goes through the right reader and says which', r.format === want, r.format);
});
const anyJson = B.parseAnySequence(bench, 'x.json');
ok('a JSON import names the dialect it was read as', /Benchling/.test(anyJson.format), anyJson.format);
ok('raw bases still import', B.parseAnySequence(SEQ, 'x')[0].seq === SEQ);
ok('nothing in, nothing out', B.parseAnySequence('   ', 'x').length === 0);

section('the formats that cannot be read say so');
['clc', 'sbd', 'ma4', 'xdna', 'svg', 'pdf', 'zip', 'scf'].forEach(ext => {
  ok('.' + ext + ' has an explanation naming a way out',
    typeof B.UNREADABLE[ext] === 'string' && B.UNREADABLE[ext].length > 40, B.UNREADABLE[ext]);
});
ok('.gb is not in that list', B.UNREADABLE.gb === undefined);

console.log('\n   ' + pass + ' passed, ' + fail + ' failed\n');
process.exit(fail ? 1 : 0);
