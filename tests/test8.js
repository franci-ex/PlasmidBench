/* ===========================================================================
   PlasmidBench — in-place sequence editing.

   An edit that crosses base 1 of a circle is checked against an independent
   route to the same answer: rotate the circle so the range no longer
   crosses the origin, make an ordinary edit, and compare. The two results
   must be the same circle, and every surviving feature must cover exactly
   the same bases in both.
   =========================================================================== */
const E = require('../src/engine2.js');
let pass = 0, fail = 0;
function ok(cond, what, extra) {
  if (cond) pass++; else { fail++; console.log('FAIL: ' + what + (extra ? '  [' + extra + ']' : '')); }
}
function prng(seed) { let x = seed; return () => (x = (x * 1103515245 + 12345) & 0x7fffffff) / 0x7fffffff; }
const r = prng(4242);
const rnd = (n) => Math.floor(r() * n);
const dna = n => { let s = ''; for (let i = 0; i < n; i++) s += 'ACGT'[rnd(4)]; return s; };
const sameCircle = (x, y) => x.length === y.length && (x + x).indexOf(y) >= 0;

/* ---------------------------------------------- 1. ordinary edits */
(function () {
  const seq = 'AAAAACCCCCGGGGGTTTTT';
  const f = [{ start: 5, end: 15, strand: 1, label: 'mid', type: 'misc_feature' }];
  let e = E.editSequence(seq, f, 10, 10, 'NNN', false);
  ok(e.seq === 'AAAAACCCCCNNNGGGGGTTTTT', 'insertion puts the bases before the given position', e.seq);
  ok(e.features[0].start === 5 && e.features[0].end === 18, 'a feature the insertion falls inside grows');
  e = E.editSequence(seq, f, 0, 5, '', false);
  ok(e.seq === 'CCCCCGGGGGTTTTT' && e.features[0].start === 0 && e.features[0].end === 10,
     'deleting before a feature moves it down');
  e = E.editSequence(seq, f, 7, 12, 'TT', false);
  ok(e.seq === 'AAAAACCTTGGGTTTTT', 'replacement swaps the range for the new bases', e.seq);
  ok(E.featureSequence(e.seq, e.features[0]) === 'CCTTGGG', 'the feature follows the replacement',
     E.featureSequence(e.seq, e.features[0]));
  e = E.editSequence(seq, f, 4, 16, '', false);
  ok(e.features.length === 0, 'a feature deleted entirely is dropped');
  let threw = false; try { E.editSequence(seq, f, 18, 23, '', false); } catch (x) { threw = true; }
  ok(threw, 'a linear sequence refuses an edit past its end');
  threw = false; try { E.editSequence(seq, f, 5, 3, '', false); } catch (x) { threw = true; }
  ok(threw, 'a range that ends before it starts is refused');
  threw = false; try { E.editSequence(seq, f, 5, 26, '', true); } catch (x) { threw = true; }
  ok(threw, 'a range covering the whole circle is refused');
})();

/* --------------------- 2. across base 1, against the rotate-then-edit route */
(function () {
  let checked = 0, featuresChecked = 0;
  for (let trial = 0; trial < 400; trial++) {
    const L = 60 + rnd(400);
    const seq = dna(L);
    const feats = [];
    for (let k = 0; k < 6; k++) {
      const st = rnd(L), len = 3 + rnd(Math.floor(L / 3));
      const en = (st + len) % L;            // end <= start means the feature wraps
      feats.push({ start: st, end: en === st ? (en + 1) % L : en, strand: r() < 0.5 ? 1 : -1,
                   label: 'f' + k, type: 'misc_feature' });
    }
    /* a range that crosses the origin */
    const a = L - 1 - rnd(Math.floor(L / 4));
    const b = L + 1 + rnd(Math.floor(L / 4));
    if (b - L > a) continue;
    const ins = r() < 0.3 ? '' : dna(rnd(30));

    const got = E.editSequence(seq, feats, a, b, ins, true);

    /* independent route: rotate so base a is base 1, edit 0..(b-a) */
    const rot = E.setOrigin({ seq: seq, features: feats, topology: 'circular' }, a);
    /* The rotated molecule is still a circle, and features still wrap in it,
       so the reference edit has to be a circular edit too. Passing false here
       used to make this route disagree with the direct one on 48 of 2,159
       features: a feature whose two surviving fragments sit on either side of
       the replacement was seen as ending at the right-hand edge rather than
       continuing through the origin, so the replacement bases fell outside it
       on one route and inside it on the other. With the flag right, the two
       routes agree on every feature of every trial. */
    const ref = E.editSequence(rot.seq, rot.features, 0, b - a, ins, true);

    ok(sameCircle(got.seq, ref.seq), 'trial ' + trial + ': the same circle either way',
       got.seq.length + ' vs ' + ref.seq.length);
    ok(got.seq.length === L - (b - a) + ins.length, 'trial ' + trial + ': length changes by exactly the edit');
    ok(got.wrapped === true, 'trial ' + trial + ': the edit is reported as crossing base 1');
    const byLabel = {};
    ref.features.forEach(f => { byLabel[f.label] = f; });
    ok(got.features.length === ref.features.length, 'trial ' + trial + ': the same features survive',
       got.features.map(f => f.label).join(',') + ' vs ' + ref.features.map(f => f.label).join(','));
    got.features.forEach(f => {
      const g = byLabel[f.label];
      if (!g) { ok(false, 'trial ' + trial + ': feature ' + f.label + ' survives only one route'); return; }
      ok(E.featureSequence(got.seq, f) === E.featureSequence(ref.seq, g),
         'trial ' + trial + ': feature ' + f.label + ' covers the same bases either way');
      featuresChecked++;
    });
    checked++;
  }
  ok(checked > 300, 'enough origin-crossing edits were tried', checked);
  console.log('  ' + checked + ' origin-crossing edits, ' + featuresChecked + ' features compared');
})();

/* -------------------------- 3. ordinary edits on circles, same comparison */
(function () {
  let checked = 0;
  for (let trial = 0; trial < 300; trial++) {
    const L = 60 + rnd(300);
    const seq = dna(L);
    const feats = [];
    for (let k = 0; k < 5; k++) {
      const st = rnd(L), len = 3 + rnd(Math.floor(L / 3));
      const en = (st + len) % L;
      feats.push({ start: st, end: en === st ? (en + 1) % L : en, strand: 1, label: 'f' + k, type: 'misc_feature' });
    }
    const a = rnd(L - 5), b = a + rnd(Math.min(20, L - a));
    const ins = dna(rnd(12));
    const got = E.editSequence(seq, feats, a, b, ins, true);
    /* the plain splice, written out independently */
    ok(got.seq === seq.slice(0, a) + ins + seq.slice(b), 'trial ' + trial + ': splice is exact');
    got.features.forEach(f => {
      const before = feats.find(x => x.label === f.label);
      const outside = before.start <= before.end
        ? (before.end <= a || before.start >= b)
        : (b <= before.start && a >= before.end);
      if (outside) ok(E.featureSequence(got.seq, f) === E.featureSequence(seq, before),
        'trial ' + trial + ': a feature clear of the edit keeps its exact bases (' + f.label + ')');
    });
    checked++;
  }
  ok(checked === 300, 'all ordinary circle edits ran');
})();

console.log('\ntest8 (editing): ' + pass + ' passed, ' + fail + ' failed');
process.exit(fail ? 1 : 0);
