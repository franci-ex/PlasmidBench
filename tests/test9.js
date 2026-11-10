/* ===========================================================================
   PlasmidBench — where primers anneal.

   E.primerSites is checked three ways: hand-built cases with a known answer;
   a brute-force oracle that tries every offset on both strands with the rule
   written out plainly; and agreement with E.pcr, the PCR simulator that has
   been in use (and cross-checked) since the start.
   =========================================================================== */
const E = require('../src/engine2.js');
let pass = 0, fail = 0;
function ok(cond, what, extra) {
  if (cond) pass++; else { fail++; console.log('FAIL: ' + what + (extra !== undefined ? '  [' + extra + ']' : '')); }
}
function prng(seed) { let x = seed; return () => (x = (x * 1103515245 + 12345) & 0x7fffffff) / 0x7fffffff; }
const r = prng(9001), rnd = n => Math.floor(r() * n);
const dna = n => { let s = ''; for (let i = 0; i < n; i++) s += 'ACGT'[rnd(4)]; return s; };
const mut = b => 'ACGT'.replace(b, '')[rnd(3)];
const footprint = (T, h) => h.end > h.start ? T.slice(h.start, h.end) : T.slice(h.start) + T.slice(0, h.end);

/* ------------------------------------------------- 1. known answers */
(function () {
  const T = dna(2000);
  const P = T.slice(500, 524);
  let h = E.primerSites(T, P, {});
  ok(h.length === 1, 'an exact forward primer is found once', h.length);
  ok(h[0] && h[0].strand === 1 && h[0].start === 500 && h[0].end === 524, 'at its exact footprint', JSON.stringify(h[0]));
  ok(h[0] && h[0].exact && h[0].annealed === 24 && h[0].tail === 0, 'fully annealed, no tail');
  ok(h[0] && h[0].threePrime === 523, 'its 3\' end is the last base of the footprint');

  const R = E.revcomp(T.slice(1200, 1226));
  h = E.primerSites(T, R, {});
  ok(h.length === 1 && h[0].strand === -1 && h[0].start === 1200 && h[0].end === 1226,
     'a reverse primer is found on the bottom strand at its footprint', JSON.stringify(h[0]));
  ok(h[0] && h[0].threePrime === 1200, 'a reverse primer\'s 3\' end is the left edge of its footprint');

  /* a 5' tail that does not match: a Gibson arm */
  const tail = 'GATCGGATCCAAGCTTGCAT';
  h = E.primerSites(T, tail + T.slice(800, 822), {});
  ok(h.length === 1 && h[0].annealed === 22 && h[0].tail === 20 && h[0].start === 800 && h[0].end === 822,
     'a 5\' tail is reported as a tail, and the footprint is only the annealed part', JSON.stringify(h[0]));
  ok(h[0] && h[0].mismatches.length === 0, 'the tail is not counted as mismatches');

  /* one internal mismatch */
  let p = T.slice(300, 325).split(''); p[9] = mut(p[9]); p = p.join('');
  h = E.primerSites(T, p, {});
  ok(h.length === 1 && h[0].annealed === 25 && h[0].mismatches.join() === '10',
     'an internal mismatch is tolerated and reported at primer base 10', JSON.stringify(h[0]));
  ok(h[0] && !h[0].exact, 'and the site is not called exact');

  /* a 3' mismatch: no extension, no site */
  for (const k of [24, 22, 20]) {
    let q = T.slice(300, 325).split(''); q[k] = mut(q[k]); q = q.join('');
    ok(E.primerSites(T, q, {}).length === 0, 'a mismatch ' + (25 - k) + ' base(s) from the 3\' end leaves no site');
  }

  /* across base 1 of a circle */
  const C = dna(1500);
  const P2 = C.slice(1490) + C.slice(0, 14);
  h = E.primerSites(C, P2, { circular: true });
  ok(h.length === 1 && h[0].start === 1490 && h[0].end === 14 && h[0].wraps,
     'on a circle a footprint across base 1 runs from 1490 to 14', JSON.stringify(h[0]));
  ok(E.primerSites(C, P2, { circular: false }).length === 0, 'and on the same sequence read as linear there is no site');
  const R2 = E.revcomp(C.slice(1495) + C.slice(0, 15));
  h = E.primerSites(C, R2, { circular: true });
  ok(h.length === 1 && h[0].strand === -1 && h[0].start === 1495 && h[0].end === 15,
     'a reverse primer across base 1 is found too', JSON.stringify(h[0]));
  /* a footprint that ends exactly on the last base */
  h = E.primerSites(C, C.slice(1480), { circular: true });
  ok(h.length === 1 && h[0].start === 1480 && h[0].end === 1500 && !h[0].wraps,
     'a footprint ending on the last base ends at the length, not at 0', JSON.stringify(h[0]));

  /* two sites */
  const D = dna(600) + T.slice(100, 125) + dna(600) + T.slice(100, 125) + dna(300);
  ok(E.primerSites(D, T.slice(100, 125), {}).length === 2, 'a primer with two sites is reported twice');

  ok(E.primerSites(T, 'ACGTACGT', {}).length === 0, 'a primer shorter than the minimum is not placed');
})();

/* ------------------------------ 2. the brute-force oracle, every offset */
function oracle(T, P, circ, o) {
  /* The rule being checked is the one E.primerSites documents: seed on the
     3'-most `seed` bases, never a mismatch inside the last `clamp`, extend
     while no more than `maxMismatchFrac` of the laid-down bases mismatch,
     and end inside a run of at least `minBlock` consecutive matches — a
     lone base that happens to pair beyond a mismatch does not anneal by
     itself. This oracle used to leave out clamp and minBlock, and so read
     56 sites longer at the 5' end than the engine did; those were the
     oracle's error, not the engine's. */
  o = Object.assign({ seed: 8, clamp: 5, minAnneal: 15, maxMismatchFrac: 0.10, minBlock: 3 }, o || {});
  const clamp = Math.min(o.clamp, o.seed);
  const L = T.length, n = P.length, out = [];
  const at = i => circ ? T[((i % L) + L) % L] : (i >= 0 && i < L ? T[i] : '');
  [1, -1].forEach(strand => {
    const top = strand === 1 ? P : E.revcomp(P);
    for (let off = -n; off < L + (circ ? 0 : n); off++) {
      /* the 3' end, base by base from the end */
      const idx = k => strand === 1 ? n - 1 - k : k;
      let seedOK = true;
      for (let s = 0; s < o.seed; s++) if (at(off + idx(s)) !== top[idx(s)]) { seedOK = false; break; }
      if (!seedOK) continue;
      let mism = 0, best = 0, run = 0, okClamp = true;
      for (let s = 0; s < n; s++) {
        const m = at(off + idx(s)) === top[idx(s)];
        if (!m) {
          if (s < clamp) { okClamp = false; break; }
          mism++; run = 0;
        } else run++;
        if (m && run >= o.minBlock && mism <= Math.floor((s + 1) * o.maxMismatchFrac)) best = s + 1;
      }
      if (!okClamp) continue;
      if (best < o.minAnneal) continue;
      let lo = strand === 1 ? off + (n - best) : off, hi = lo + best;
      if (!circ && (lo < 0 || hi > L)) continue;
      if (lo < 0) { lo += L; hi += L; }
      if (circ && lo >= L) continue;          // the same site, one turn later
      out.push(strand + ':' + lo + ':' + (hi <= L ? hi : hi - L) + ':' + best);
    }
  });
  return Array.from(new Set(out)).sort();
}
(function () {
  let cases = 0, sites = 0;
  for (let t = 0; t < 400; t++) {
    const circ = t % 2 === 0;
    const L = 200 + rnd(800), T = dna(L);
    /* primers taken from the template, either strand, sometimes across the
       origin, with a tail and internal mismatches added at random */
    const n = 18 + rnd(14), st = rnd(circ ? L : L - n);
    let core = (T + T).slice(st, st + n);
    if (r() < 0.5) core = E.revcomp(core);
    let p = core.split('');
    const nm = rnd(4);
    for (let k = 0; k < nm; k++) { const at = rnd(n); p[at] = mut(p[at]); }
    p = (r() < 0.3 ? dna(5 + rnd(20)) : '') + p.join('');
    /* make some templates repetitive, so several sites compete */
    const T2 = r() < 0.2 ? T.slice(0, L >> 1) + core + T.slice(L >> 1) : T;
    const got = E.primerSites(T2, p, { circular: circ }).map(h => h.strand + ':' + h.start + ':' + h.end + ':' + h.annealed).sort();
    const want = oracle(T2, p, circ);
    ok(JSON.stringify(got) === JSON.stringify(want), 'oracle case ' + t + (circ ? ' (circular)' : ' (linear)'),
       JSON.stringify(got) + ' vs ' + JSON.stringify(want));
    cases++; sites += want.length;
  }
  console.log('  ' + cases + ' templates against the oracle, ' + sites + ' sites');
})();

/* ------------------------- 3. agreement with the PCR simulator, E.pcr */
(function () {
  let checked = 0;
  for (let t = 0; t < 200; t++) {
    const circ = t % 2 === 1;
    const L = 800 + rnd(2000), T = dna(L);
    const a = rnd(L - 400), len = 150 + rnd(300);
    const f = T.slice(a, a + 20 + rnd(6));
    const b = a + len;
    const rv = E.revcomp(T.slice(b - 22, b));
    const prods = E.pcr(T, f, rv, { circular: circ });
    if (prods.length !== 1) continue;
    const fs = E.primerSites(T, f, { circular: circ }).filter(h => h.strand === 1);
    const rs = E.primerSites(T, rv, { circular: circ }).filter(h => h.strand === -1);
    ok(fs.length === 1 && rs.length === 1, 'pcr case ' + t + ': one site for each primer');
    if (fs.length !== 1 || rs.length !== 1) continue;
    /* E.pcr reports where the 15-base 3' anchors sit */
    ok(fs[0].end - 15 === prods[0].fwdSite, 'pcr case ' + t + ': forward site agrees with E.pcr',
       (fs[0].end - 15) + ' vs ' + prods[0].fwdSite);
    ok(rs[0].start === prods[0].revSite, 'pcr case ' + t + ': reverse site agrees with E.pcr',
       rs[0].start + ' vs ' + prods[0].revSite);
    ok(footprint(T, fs[0]) === f && E.revcomp(footprint(T, rs[0])) === rv, 'pcr case ' + t + ': footprints are the primers');
    checked++;
  }
  ok(checked > 150, 'enough PCR cases compared', checked);
  console.log('  ' + checked + ' primer pairs agree with E.pcr');
})();

/* ----------------------------------------------------------- 4. speed */
(function () {
  const T = dna(10000), ps = [];
  for (let i = 0; i < 100; i++) { const a = rnd(9970); ps.push(T.slice(a, a + 22)); }
  const t0 = Date.now();
  ps.forEach(p => E.primerSites(T, p, { circular: true }));
  const ms = Date.now() - t0;
  ok(ms < 1500, '100 primers on a 10 kb circle in well under two seconds', ms + ' ms');
  console.log('  100 primers x 10 kb: ' + ms + ' ms');
})();

console.log('\ntest9 (primer sites): ' + pass + ' passed, ' + fail + ' failed');
process.exit(fail ? 1 : 0);
