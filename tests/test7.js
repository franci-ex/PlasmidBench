/* ===========================================================================
   PlasmidBench — read triage tests.

   The frame call and the read verdict decide what a user BLASTs, so every
   rule in E.readTriage is pinned here. The peptides marked REAL / NONSENSE
   are verbatim frames from a real Sanger run of CDK and Speedy constructs
   in which the old "longest stop-free stretch" rule picked the wrong frame
   on five of forty-nine reads.
   =========================================================================== */
const E = require('../src/engine2.js');
const path = require('path');
const fs = require('fs');

let pass = 0, fail = 0;
function ok(cond, what, extra) {
  if (cond) { pass++; } else { fail++; console.log('FAIL: ' + what + (extra ? '  [' + extra + ']' : '')); }
}
function eq(a, b, what) { ok(a === b, what, 'got ' + JSON.stringify(a) + ' want ' + JSON.stringify(b)); }
function near(a, b, tol, what) { ok(Math.abs(a - b) <= tol, what, 'got ' + a + ' want ' + b + '±' + tol); }

/* ------------------------------------------------------------------ data */
/* Real coding frames from the run. */
const REAL = {
  polA: 'TAEHLYEGLLKVTRPDSHKVHTRFNQALTQTGRLSSTDPNLQNIPIRLEEGRKIRQAFVPSEKDWLIFAADYSQIELRVLAHISKDENLIEAFTNDMDIHTKTAMDVFHVAKDEVTSAMRRQAKAVNFGIVYGISDYGLSQNLGITRKEAGAFIDRYLESFQGVKAYMEDSVQEAKQKGYVTTLMHRRRYIPELTSRNFNIRSFAERTAMNTPIQGSAADIIKKAMIDMAAKLKEKQLKARLLLQVHDELIFEAPKEEIEILEKLVPEVMEHALALDVPXKXX',
  recA: 'EGDLNMSDLKSRLIKASTSKLTAELTASKFFNEKDVVRTKIPMMNIALSGEITGGMQSGLLILAGPSKSFKSNFGLTMVSSYMRQYPDAVCLFYDSEFGITPAYLRSMGVDPERVIHTPVQSLEQLRIDMVNQLDAIERGEKVVVFIDSLGNLASKKETEDALNEKVVSDMTRAKTMKSLFRIVTPYFSTKNIPCIAINHTYETQEMFSKTVMGGGTGPMYSADTVFIIGKRQIKDGSDLQGYQFVLNVEKSRTVKEKSKFFIDVKFDGGIDPYSXCWIWHWXWVLLXNRKMVGMPVXXG',
  fubi: 'EGDIHMHHHHHHSSGVDLGTENLYFQSMQLFVRAQELHTFEVTGQETVAQIKAHVASLEGIAPEDQVVLLAGAPLEDEATLGQCGVEALTTLEVAGRMLGGKVHGSLARAGKVRGQTPKVAKQEKKKKKTGRAKRRMQYNRRFVNVVPTFGKKKGPNANSSKGGYGLNDIFEAQKIEWHE',
  cdk1: 'XNNNNXXLGENLYFQGMEDYTKIEKIGEGTYGVVYKGRHKTTGQVVAMKKIRLESEEEGVPSTAIREISLLKELRHPNIVSLQDVLMQDSRLYLIFEFLSMDLKKYLDSIPPGQYMDSSLVKSYLYQILQGIVFCHSRRVLHRDLKPQNLLIDDKGTIKLADFGLARAFGIPIRVYTHEVVTLWYRSPEVLLGSARYSTPVDIWSIGTIFQNWRRKNHCFTGIRRSISCSVFSGRSALQTTRYGLKSRVFKTTKIRFRSGSQARLRPM'
};
/* Nonsense frames of the SAME reads, which the old rule preferred. */
const NONSENSE = {
  polAother: 'XXFXRNIQCQCMFHNFRNQLFQDLDFFFRCFKDQFIMHLQQQTCFQLFFFQLRCHVNHCFFDDICCATLNRGIHCGTFCKATNIEITAGQFRNITTTMHQCGNITFLFCLLNAIFHVCFNALKTFQIAINKCACFFTGNTQILAQAIIADAIHNTKINCFGLTTHCAGNFIFRNMKHIHGCFGMNIHIVGEGFDQIFIFANMCQNTQFNLAIICCKNQPIFFARHKCLTNFTTFFQTNRNILQIRIGAAQTTGLGQCLIKTGMNFMAIRTGNFQQTFI',
  recAother: 'PAFXXTDHFXXXSSXXHGHXNHFSVXQQNPVPVPYPAXTIRINATIEFHVDKEFAFLFHGAAFFNIQNKLITLQIATIFNLAFANDKNGICAIHRTGATTHNGFAKHLLGFIGMVDGNTRNVFGAEIRGNNTKQAFHGFGTGHIANNFFIQGIFGFFFACQITQAINKHHHFFTTFNCIQLINHINTQLFQALNRCMNNTFRINTHATQICRCDTKFAIIKQTNCIRILTHIAAHHGQTKIAFKAFARTSQNQQTALHTTGNFTAQCNIHHRNFRAHHIFFVEKFACCQFCGQFAGACFNQTAFQIAHI',
  fubiother: 'CAATGRVPFANPDIVPPFSKKPLKTRLEAPRGYASYCSAVAAANSASFRALLAAGSQWWWWWCSSAAASLSTELEFGSLSFMPFNLLRFKDIIQAVSTFTRVGIGALLLAKGGHNVDKAPVVLHPPLSPTCLLLLLLLFGHLRSLTSHFSSTGQGTMDFTSKHAACYFQGSQGLHPALAQSGLILQGRACQEHDLIFRGNALQ',
  vectorjunk: 'CAATGRVPFANPDIVPPFSKKPLKTRLEAPRGYASYCSAVAAANSASFRALLAAGSQWWWWWCSSAAASFVIMALEIQVFSAAAVMMMMMHGISPS'
};

/* One codon per amino acid, for building a read whose true frame is known. */
const BT = { A: 'GCT', R: 'CGT', N: 'AAT', D: 'GAT', C: 'TGT', Q: 'CAA', E: 'GAA',
             G: 'GGT', H: 'CAT', I: 'ATT', L: 'CTG', K: 'AAA', M: 'ATG', F: 'TTT',
             P: 'CCT', S: 'TCT', T: 'ACT', W: 'TGG', Y: 'TAT', V: 'GTT', X: 'NNN' };
const back = p => p.split('').map(c => BT[c] || 'NNN').join('');

/* ------------------------------------- 1. composition against Swiss-Prot */
(function () {
  Object.keys(REAL).forEach(k => {
    const f = E.compositionFit(REAL[k]);
    ok(f !== null && f < 0.40, 'real frame ' + k + ' has protein-like composition', 'fit=' + f);
  });
  Object.keys(NONSENSE).forEach(k => {
    const f = E.compositionFit(NONSENSE[k]);
    ok(f !== null, 'nonsense frame ' + k + ' is long enough to score');
  });
  ok(E.compositionFit(NONSENSE.polAother) > 1.0, 'F/Q/C-rich nonsense scores far from a protein',
     'fit=' + E.compositionFit(NONSENSE.polAother));
  ok(E.compositionFit(NONSENSE.recAother) > 0.8, 'second nonsense frame scores far from a protein');
  eq(E.compositionFit('MVLSPADK'), null, 'peptide under 40 residues gets no composition score');
  /* The statistic is what breaks a tie the stop count cannot: both frames of
     this read are ~280 aa with 2 and 3 stops. */
  ok(E.compositionFit(REAL.polA) * 10 < E.compositionFit(NONSENSE.polAother),
     'composition separates two same-length frames by more than tenfold');
})();

/* ------------------------------------------------- 2. peptide landmarks */
(function () {
  const h = E.findLandmarks(REAL.fubi).map(x => x.name);
  ok(h.indexOf('His6') >= 0, 'His6 found in a pET leader');
  ok(h.indexOf('TEV site') >= 0, 'TEV site found in a pET leader');
  ok(h.indexOf('AviTag') >= 0, 'AviTag found at the C-terminus');
  ok(h.indexOf('His6-TEV (pET)') >= 0, 'the whole His6-TEV cassette is found as one landmark');

  /* Ambiguous residues inside a landmark are tolerated: this is the real
     reason the old code missed the TEV site on seven reads of this run. */
  const two = E.findLandmarks('GXNLYFXGMRRVSLVLPEHHEAFNRLLEDPVIKRFLAWDKDLRVSDKYLLAMVIAYFSRAG');
  ok(two.some(x => x.name === 'TEV site'), 'TEV site found with two N calls inside it');
  ok(two.some(x => x.name === 'TEV site' && !x.exact), 'the partial match is reported as partial');

  /* Three of six residues ambiguous is past what strict matching will take,
     and is only accepted on a read whose bases are otherwise confident.
     Measured on 294 real frames: strict pins 19 reads, relaxed pins 23, and
     neither ever pins two frames of the same read. */
  const three = 'GXXLYFXGMRRVSLVLPEHHEAFNRLLEDPVIKRFLAWDKDLRVSDKYLLAMVIAYFSRAGLPSWQYQRIHFF';
  eq(E.findLandmarks(three).filter(x => x.name === 'TEV site').length, 0,
     'strict matching refuses a six-residue site with three N calls');
  ok(E.findLandmarks(three, [], { relaxed: true }).some(x => x.name === 'TEV site'),
     'relaxed matching accepts it on an otherwise clean read');

  /* But an ambiguity run must not match everything. */
  eq(E.findLandmarks('AAXXXXXXXXXXXXXXXXXXXXXXXXAA').length, 0,
     'a run of X matches no landmark');
  eq(E.findLandmarks('AAXXXXXXXXXXXXXXXXXXXXXXXXAA', [], { relaxed: true }).length, 0,
     'nor under relaxed matching');
  eq(E.findLandmarks('HVPXQKVXAXKSXADQLGASVVVSXVGGNGVQASKAAVHNLLAPRVSGLII')
     .filter(x => x.name === 'myc').length, 0,
     'a window that is half ambiguous does not match the myc tag');
  eq(E.findLandmarks('RLQFRWHFSGKCARNPYLFIFLXHSNMYPLMXXXXKNSSSINETXXXXSGXSXPXFXK', [], { relaxed: true })
     .filter(x => x.pins).length, 0,
     'a tag matched inside a run of N calls does not pin a broken frame');

  /* No landmark may fire in two frames of the same read — that would make
     the signal useless for pinning a frame. */
  Object.keys(REAL).forEach(k => {
    Object.keys(NONSENSE).forEach(j => {
      if (E.findLandmarks(REAL[k]).length && E.findLandmarks(NONSENSE[j]).length)
        ok(false, 'landmark fired in both a real and a nonsense frame: ' + k + ' / ' + j);
    });
  });
  pass++; // the loop above only fails
  ok(E.PEPTIDE_LANDMARKS.every(l => l.seq.length >= 6), 'no landmark is shorter than 6 residues');
})();

/* ------------------------------------------------ 3. the frame call itself */
(function () {
  /* Scored against each other, the real frame must win even when the
     nonsense frame is longer. These are the five reads the old rule got
     wrong, in peptide form. */
  const pairs = [['polA', 'polAother'], ['recA', 'recAother'], ['fubi', 'fubiother']];
  pairs.forEach(([r, n]) => {
    const rf = E.compositionFit(REAL[r]), nf = E.compositionFit(NONSENSE[n]);
    const cw = E.TRIAGE.compWeight;
    const cf = f => Math.max(E.TRIAGE.compFloor, Math.min(1, 1 / (1 + cw * f)));
    const rs = REAL[r].replace(/X/g, ' ').split(' ').reduce((a, b) => Math.max(a, b.length), 0) * cf(rf)
             + (E.findLandmarks(REAL[r]).length ? E.TRIAGE.landmarkBonus : 0);
    const ns = NONSENSE[n].replace(/X/g, ' ').split(' ').reduce((a, b) => Math.max(a, b.length), 0) * cf(nf)
             + (E.findLandmarks(NONSENSE[n]).length ? E.TRIAGE.landmarkBonus : 0);
    ok(rs > ns, 'the real frame outscores the nonsense frame for ' + r,
       'real=' + rs.toFixed(0) + ' nonsense=' + ns.toFixed(0));
  });
})();

/* ---------------------------------------- 4. six frames of a built read */
const LEADER = 'ATGGGCAGCAGCCATCATCATCATCATCACAGCAGCGGCGTGGATCTGGGTACCGAGAACCTGTACTTCCAGAGC';
const CDS = LEADER + back('MEDYTKIEKIGEGTYGVVYKGRHKTTGQVVAMKKIRLESEEEGVPSTAIREISLLKELRHPNIVSLQDVLMQDSRLYLIFEFLSMDLKKYLDSIPPGQYMDSSLVKSYLYQILQGIVFCHSRRVLHRDLKPQNLLIDDKGTIKLADFGLARAFGIPIRVYTHEVVTLWYRSPEVLLGSARYSTPVDIWSIGTIFAELATKKPLFHGDSEIDQLFRIFRALGTPNNEVWPEVESLQDYKNTFPKWKPGSLASHVKNLDENGLDLLSKMLIYDPAKRISGKMALNHPYFNDLDNQIKKM') + 'TAA';
const FLANK5 = 'GATCCGGCGTAGAGGATCGAGATCTCGATCCCGCGAAATTAATACGACTCACTATAGGG';
/* a real T7 terminator: non-coding, and stop-riddled in every frame, so
   it cannot accidentally supply a tag to a shifted frame. */
const FLANK3 = 'CTAGCATAACCCCTTGGGGCCTCTAAACGGGTCTTGAGGGGTTTTTTGCTGAAAGGAGG';
const READ = FLANK5 + CDS + FLANK3;

(function () {
  const fr = E.frameTranslate(READ);
  eq(fr.length, 6, 'six frames are returned');
  const labels = fr.map(f => f.label).sort().join(',');
  eq(labels, '+1,+2,+3,−1,−2,−3', 'frames are labelled +1..+3 and −1..−3');
  /* every frame agrees with a plain translate of the same offset */
  fr.forEach(f => {
    const s = f.dir === 1 ? READ : E.revcomp(READ);
    const nt = s.slice(f.frame, f.frame + Math.floor((s.length - f.frame) / 3) * 3);
    eq(f.full, E.translate(nt), 'frame ' + f.label + ' matches E.translate');
  });
  fr.forEach(f => {
    ok(f.cleanLen <= f.longest.length, 'frame ' + f.label + ': clean run never exceeds the stop-free run');
  });
})();

/* ------------------------------------------- 5. verdicts on a built read */
(function () {
  const t = E.readTriage(READ);
  eq(t.verdict, 'clean', 'a tagged full-length ORF reads as clean');
  eq(t.chosen.label, '+3', 'the coding frame is called');
  ok(t.chosen.landmarks.length > 0, 'the call is pinned by a landmark');
  ok(t.chosen.landmarks.some(l => l.name === 'TEV site'), 'the TEV site is the landmark');
  ok(t.dominance > 1.5, 'the called frame dominates');
  eq(t.shift, null, 'no indel is reported for an intact ORF');
  ok(t.flags.some(f => f.level === 'ok'), 'a pinned call is flagged as pinned');
  ok(t.chosen.landmarks.some(l => l.pins), 'the pinning landmark is marked as pinning');
  ok(t.chosen.clean.indexOf('MEDYTKIEKIGEGTYGVVYKG') > 0, 'the peptide contains the insert');

  /* The same read on the other strand must give the same peptide. */
  const r = E.readTriage(E.revcomp(READ));
  eq(r.verdict, 'clean', 'the reverse complement also reads as clean');
  ok(r.chosen.label[0] === '−', 'the call lands on a minus frame');
  eq(r.chosen.clean, t.chosen.clean, 'the peptide is identical on either strand');
})();

/* -------------------------------------------------- 6. a failed trace */
(function () {
  /* One base in eight replaced by N: the stop count drops because N
     translates to X, which is exactly the trap this verdict exists for. */
  let s = '';
  for (let i = 0; i < READ.length; i++) s += (i % 8 === 0 ? 'N' : READ[i]);
  const t = E.readTriage(s);
  eq(t.verdict, 'low-quality', 'an N-riddled read is condemned');
  ok(t.ambiguousFrac > E.TRIAGE.ambigBad, 'the ambiguous fraction is reported');
  ok(t.flags.some(f => f.level === 'bad' && /N translates to X/.test(f.text)),
     'the verdict explains why a low stop count means nothing here');
  /* and it is not rescued by having few stops */
  ok(t.chosen.stops < 20, 'the dead trace does indeed show few stops', 'stops=' + t.chosen.stops);
})();

/* ----------------------------------------------- 7. a single-base indel */
(function () {
  const at = FLANK5.length + LEADER.length + 300;           // inside the CDS
  const s = READ.slice(0, at) + READ.slice(at + 1);          // delete one base
  const t = E.readTriage(s);
  ok(t.shift !== null, 'a one-base deletion is detected');
  if (t.shift) {
    near(t.shift.at, at, 60, 'the indel is located near the deleted base');
    ok(t.shift.from <= at + 1 && t.shift.to >= at - 1,
       'the reported window brackets the deleted base',
       'window ' + t.shift.from + '-' + t.shift.to + ' vs ' + at);
    ok(t.shift.continues >= E.TRIAGE.shiftMin, 'the shifted frame carries a long peptide');
    ok(t.flags.some(f => /inserted or deleted/.test(f.text)), 'the flag explains the split');
  }
  ok(t.verdict === 'frameshift' || t.chosen.landmarks.length > 0,
     'the read is either called a frameshift or still pinned by its tag');
})();

/* --------------------------------------- 8. non-coding vector sequence */
(function () {
  const VEC = path.resolve(__dirname, 'fixtures/pRSFDuet-1.dna');
  if (!fs.existsSync(VEC)) { console.log('SKIP: vector file not present'); return; }
  const buf = fs.readFileSync(VEC);
  const rec = E.parseSnapGene(buf.buffer.slice(buf.byteOffset, buf.byteOffset + buf.byteLength));
  ok(rec && rec.seq.length > 3000, 'the real vector parses', 'len=' + (rec && rec.seq.length));

  /* A stretch of the backbone that carries no coding sequence must not be
     offered as a protein. */
  const cds = (rec.features || []).filter(f => (f.type || '').toUpperCase() === 'CDS');
  const covered = new Set();
  cds.forEach(f => E.featurePieces(f, rec.seq.length).forEach(p => {
    for (let i = p.start; i < p.end; i++) covered.add(i);
  }));
  let bestStart = -1, bestLen = 0, run = 0;
  for (let i = 0; i < rec.seq.length; i++) {
    if (!covered.has(i)) { run++; if (run > bestLen) { bestLen = run; bestStart = i - run + 1; } }
    else run = 0;
  }
  ok(bestLen >= 300, 'the vector has a long non-coding stretch to test with', 'len=' + bestLen);
  const nc = rec.seq.substr(bestStart, Math.min(bestLen, 700));
  const t = E.readTriage(nc);
  ok(['no-protein', 'frame-ambiguous', 'open-only'].indexOf(t.verdict) >= 0,
     'a non-coding backbone read is not presented as a protein', 'verdict=' + t.verdict);
  ok(t.verdict !== 'clean' && t.verdict !== 'partial', 'and is never called clean or real');
  ok(t.flags.some(f => f.level === 'warn'), 'and it carries a warning rather than nothing at all');

  /* Landmarks taken from the construct pin a read that carries no tag. */
  const lms = E.landmarksFromRecord(rec);
  ok(lms.length > 0, 'landmarks are derived from the vector CDS features', lms.length + ' probes');
  ok(lms.every(l => l.seq.length >= 14), 'construct landmarks are long enough to survive N calls');
  const big = cds.slice().sort((a, b) =>
    E.featureLength(b, rec.seq.length) - E.featureLength(a, rec.seq.length))[0];
  const inside = E.featureSequence(rec.seq, big).substr(120, 600);
  const withLm = E.readTriage(inside, { landmarks: lms });
  ok(withLm.chosen.landmarks.length > 0, 'a read from a vector gene is pinned by construct landmarks');
  eq(withLm.chosen.label, '+1', 'and the pinned frame is the gene frame');
})();

/* ---------------------------------------------------- 9. short reads */
(function () {
  eq(E.readTriage('ACGTACGTACGT').verdict, 'too-short', 'a 12-base read is refused');
  const t = E.readTriage('ACGT'.repeat(20));
  ok(t.verdict === 'no-protein' || t.verdict === 'frame-ambiguous',
     'a repetitive 80-base read yields no protein', 'verdict=' + t.verdict);
})();

/* ------------------------------------------------- 10. duplicate reads */
(function () {
  const a = READ, b = READ.slice(20, READ.length - 15), c = E.revcomp(READ.slice(5));
  const other = (function () { let s = ''; const r = (function () { let x = 7; return () => (x = (x * 1103515245 + 12345) % 2147483648) / 2147483648; })();
    for (let i = 0; i < 700; i++) s += 'ACGT'[Math.floor(r() * 4)]; return s; })();
  const g = E.duplicateReads([{ seq: a }, { seq: other }, { seq: b }, { seq: c }]);
  eq(g.length, 1, 'one duplicate group is found');
  eq(g[0].slice().sort().join(','), '0,2,3', 'the three copies group together and the unrelated read does not');

  eq(E.duplicateReads([{ seq: a }, { seq: other }]).length, 0, 'two unrelated reads form no group');
  /* strings are accepted as well as read objects */
  eq(E.duplicateReads([a, a]).length, 1, 'plain strings are accepted');
  /* a read that is nearly all N has no k-mers and cannot be grouped */
  eq(E.duplicateReads([{ seq: 'N'.repeat(400) }, { seq: 'N'.repeat(400) }]).length, 0,
     'ambiguous reads are not grouped with each other');
})();

/* ------------------------------- 12. N-rich ends on a read with no quality */
(function () {
  /* a typical .seq file: a ragged N-heavy start, the good read, a dead tail */
  /* The ends have to be long enough to actually matter: E.TRIAGE.ambigBad is
     0.08, so on a 1,087 base read, 83 N-rich bases of flank come to 6.2% and
     are correctly called clean. The point being made here is that untrimmed
     ends CAN condemn a read, which needs the fraction over the threshold —
     62 + 92 bases of flank put it at 11.3%. */
  const head = 'NNNNNNNNNNNNGNNANTNNNCNNNNGN' + 'NNNNNNNNTNNNNNNNNNANNNNNNNNNNCNNNN';
  const tail = 'NNANNNTNNNNGNNNNNNCNNNNANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN' +
               'NNNNNNNNNNNNNNNNNNNNNNNNNNGNNNNNNNNNN';
  const raw = head + READ + tail;
  const sp = E.usableSpan(raw);
  ok(sp.start >= head.length - 8 && sp.start <= head.length + 2, 'the ragged start is cut back',
     'start=' + sp.start + ' head=' + head.length);
  ok(sp.end >= head.length + READ.length - 2 && sp.end <= head.length + READ.length + 8,
     'the dead tail is cut back', 'end=' + sp.end);
  ok('ACGT'.indexOf(raw[sp.start]) >= 0 && 'ACGT'.indexOf(raw[sp.end - 1]) >= 0,
     'the span starts and ends on a real base');

  const whole = E.readTriage(raw, { trimEnds: false });
  const judged = E.readTriage(raw);
  eq(whole.verdict, 'low-quality', 'judged whole, the N-rich ends would condemn the read');
  eq(judged.verdict, 'clean', 'judged on its readable span, the same read is clean');
  ok(judged.flags.some(f => /Judged on bases/.test(f.text)), 'the trimming is reported');
  eq(judged.chosen.clean, E.readTriage(READ).chosen.clean, 'and the peptide is the same as for the clean read');

  /* a read that is dead throughout is not rescued by trimming */
  let dead = '';
  for (let i = 0; i < READ.length; i++) dead += (i % 6 === 0 ? 'N' : READ[i]);
  eq(E.readTriage(dead).verdict, 'low-quality', 'trimming does not rescue a read that is dead throughout');
  const sd = E.usableSpan(dead);
  eq(sd.trimmed, 0, 'and such a read is left whole');

  eq(E.usableSpan('ACGTACGT').trimmed, 0, 'a read shorter than the window is left alone');
  eq(E.usableSpan(READ).trimmed, 0, 'a read with no ambiguous ends loses nothing');
})();

/* -------------------------- 13. indel positions are in the read's own bases */
(function () {
  const at = FLANK5.length + LEADER.length + 300;   // 0-based index of the deleted base
  const fwd = READ.slice(0, at) + READ.slice(at + 1);
  const t1 = E.readTriage(fwd);
  ok(t1.shift && t1.shift.readFrom <= at + 1 && t1.shift.readTo >= at - 1,
     'forward read: the reported window brackets the deletion',
     t1.shift ? t1.shift.readFrom + '-' + t1.shift.readTo + ' vs ' + (at + 1) : 'no shift');

  /* the same read sequenced from the other end: the deletion sits at a
     different index of the read, and the report has to follow it there */
  const rev = E.revcomp(fwd);
  const atRev = fwd.length - at;                    // where the deletion lands, 1-based, in rev
  const t2 = E.readTriage(rev);
  ok(t2.chosen.dir === -1, 'the reverse read is called on a minus frame');
  ok(t2.shift && t2.shift.readFrom <= atRev + 2 && t2.shift.readTo >= atRev - 2,
     'reverse read: the window is given in the reverse read\'s own bases',
     t2.shift ? t2.shift.readFrom + '-' + t2.shift.readTo + ' vs ' + atRev : 'no shift');

  /* and with a trimmed ragged start, positions still count from the
     first base of the read as given */
  const pad = 'NNNNNNNNNNNNNNNNNNNNNNNNN';
  const t3 = E.readTriage(pad + fwd);
  ok(t3.shift && t3.shift.readFrom <= at + pad.length + 1 && t3.shift.readTo >= at + pad.length - 1,
     'after trimming, positions still count from the start of the read',
     t3.shift ? t3.shift.readFrom + '-' + t3.shift.readTo + ' vs ' + (at + pad.length + 1) : 'no shift');
  ok(t1.flags.some(f => /of the read/.test(f.text) && /trace/.test(f.text)),
     'the flag gives the position in the read and says to check the trace');
})();

/* ---------------------------------------------- 14. reads from pasted text */
(function () {
  const many = '>JM02-T7 some description\nACGTACGTAC\nGTACGT\n>JM03-T7\nTTTTGGGGCCCC\n\n>JM04-T7\r\nAAAACCCC\r\n';
  const r = E.parseReadsText(many, 'paste');
  eq(r.length, 3, 'a multi-FASTA paste gives one read per record, not one long read');
  eq(r.map(x => x.name).join(','), 'JM02-T7,JM03-T7,JM04-T7', 'each read is named from its header');
  eq(r[0].seq, 'ACGTACGTACGTACGT', 'wrapped lines are joined within a record');
  eq(r[2].seq, 'AAAACCCC', 'Windows line endings are handled');

  const plain = E.parseReadsText('acgt acgt\nACGTNNNN\n', 'clone7.seq');
  eq(plain.length, 1, 'bare bases give one read');
  eq(plain[0].name, 'clone7.seq', 'named after the file');
  eq(plain[0].seq, 'ACGTACGTACGTNNNN', 'case and spacing are normalised, N kept');

  const noname = E.parseReadsText('>\nACGT\n>\nGGGG\n', 'batch');
  eq(noname.map(x => x.name).join(','), 'batch_1,batch_2', 'records without names are numbered');
  eq(E.parseReadsText('>empty\n\n>real\nACGT\n', 'x').length, 1, 'an empty record is dropped');
  eq(E.parseReadsText('', 'x').length, 0, 'empty text gives no reads');
  eq(E.parseReadsText('>a\nACGU\n', 'x')[0].seq, 'ACGT', 'U is read as T');
})();

/* -------------------------------------- 11. every verdict is documented */
(function () {
  const seen = ['clean', 'partial', 'open-only', 'frameshift', 'frame-ambiguous',
                'low-quality', 'no-protein', 'too-short'];
  seen.forEach(v => ok(!!E.VERDICTS[v], 'verdict ' + v + ' has a label'));
  eq(Object.keys(E.VERDICTS).length, seen.length, 'no undocumented verdicts');
  seen.forEach(v => ok(!!E.VERDICTS_SHORT[v], 'verdict ' + v + ' has a short label'));
  eq(Object.keys(E.VERDICTS_SHORT).length, seen.length, 'short labels match the verdict set');
  ['minLen', 'ambigBad', 'cleanMin', 'cleanGood', 'domMin', 'shiftMin', 'landmarkBonus',
    'compWarn', 'shiftFit', 'pinMinClean', 'pinGap']
    .forEach(k => ok(typeof E.TRIAGE[k] === 'number', 'threshold ' + k + ' is exposed'));
  ok(E.AA_BACKGROUND_SOURCE.length > 10, 'the composition table names its source');
  const tot = Object.keys(E.AA_BACKGROUND).reduce((s, k) => s + E.AA_BACKGROUND[k], 0);
  near(tot, 100, 0.5, 'the background composition sums to 100%');
  eq(Object.keys(E.AA_BACKGROUND).length, 20, 'all twenty residues have a background frequency');
})();

console.log('\ntest7 (read triage): ' + pass + ' passed, ' + fail + ' failed');
process.exit(fail ? 1 : 0);
