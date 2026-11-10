/* ===========================================================================
   Sequencing-from-stock: tests.

   The approach is to build the truth first and then see whether the engine
   recovers it. A primer is planted at a known position on a random
   template, so the window it reads and the part of a region it covers are
   known by arithmetic before the engine is asked. Nothing here recomputes
   the engine's rules in the engine's own style, which would only prove it
   agrees with itself.

   The plan is then checked against brute force: every base the plan claims
   is re-derived by marking bases one at a time, and every base the plan
   gives up on is checked against every candidate to be sure none of them
   could have read it.
   ========================================================================= */
const E = require('../src/engine3.js');

let pass = 0, fail = 0;
function ok(name, cond, extra) {
  if (cond) { pass++; }
  else { fail++; console.log('  FAIL  ' + name + (extra !== undefined ? '   -> ' + extra : '')); }
}
function section(s) { console.log('\n== ' + s + ' =='); }

/* a deterministic pseudo-random generator, so a failure is reproducible */
let seed = 20261029;
function rnd() { seed = (seed * 1103515245 + 12345) & 0x7fffffff; return seed / 0x7fffffff; }
function randSeq(n) { let s = ''; for (let i = 0; i < n; i++) s += 'ACGT'[Math.floor(rnd() * 4)]; return s; }
const mod = (a, n) => ((a % n) + n) % n;

/* ======================================================================= */
section('the read window, measured from the 3′ end');

{
  const L = 5000;
  const w = E.sangerWindow({ strand: 1, threePrime: 1000, start: 980, end: 1001 }, L,
    { dead: 35, good: 700, max: 900 });
  ok('a forward read begins dead bases after the primer', w.good.start === 1036, w.good.start);
  ok('and ends at the confident read length', w.good.end === 1701, w.good.end);
  ok('its length is good minus dead', w.good.length === 665, w.good.length);
  ok('the furthest reach is max from the primer', w.reach.end === 1901, w.reach.end);
  ok('the unreadable window sits right after the primer', w.dead.start === 1001 && w.dead.end === 1036,
     JSON.stringify(w.dead));

  const r = E.sangerWindow({ strand: -1, threePrime: 3000, start: 3000, end: 3020 }, L,
    { dead: 35, good: 700, max: 900 });
  ok('a reverse read runs the other way', r.good.end === 2965 && r.good.start === 2300,
     JSON.stringify(r.good));
  ok('and reaches max in that direction', r.reach.start === 2100, r.reach.start);
  ok('its unreadable window is just below the primer', r.dead.start === 2965 && r.dead.end === 3000,
     JSON.stringify(r.dead));

  /* a read running off the end of a linear template is cut, not wrapped */
  const lin = E.sangerWindow({ strand: 1, threePrime: 4800, start: 4780, end: 4801 }, L,
    { dead: 35, good: 700, max: 900, circular: false });
  ok('on a linear template the read stops at the end', lin.good.end === 5000 && !lin.good.wraps,
     JSON.stringify(lin.good));
  /* but on a circle it carries on past base 1 */
  const cir = E.sangerWindow({ strand: 1, threePrime: 4800, start: 4780, end: 4801 }, L,
    { dead: 35, good: 700, max: 900, circular: true });
  ok('on a circle it carries on across the origin', cir.good.wraps && cir.good.length === 665,
     JSON.stringify(cir.good));
  ok('and its end has come round to a low coordinate', cir.good.end === mod(4800 + 1 + 700, L),
     cir.good.end);
}

/* ======================================================================= */
section('a planted primer is found where it was planted');

{
  let tested = 0;
  for (let trial = 0; trial < 400; trial++) {
    const L = 2000 + Math.floor(rnd() * 4000);
    let T = randSeq(L);
    const circ = rnd() < 0.5;
    const strand = rnd() < 0.5 ? 1 : -1;
    const len = 18 + Math.floor(rnd() * 7);
    /* plant the primer far from the ends so a linear template is not a
       special case here */
    const at = 300 + Math.floor(rnd() * (L - 1200));
    const footprint = T.slice(at, at + len);
    if (!/^[ACGT]+$/.test(footprint)) continue;
    const primer = { name: 'planted', seq: strand === 1 ? footprint : E.revcomp(footprint) };
    /* its 3' end, by definition of where it was planted */
    const tp = strand === 1 ? at + len - 1 : at;
    /* a region placed a known distance downstream of the read direction */
    const dead = 35, good = 700;
    const lead = 50 + Math.floor(rnd() * 300);
    const rlen = 100 + Math.floor(rnd() * 400);
    const rs = strand === 1 ? tp + 1 + lead : tp - lead - rlen;
    if (rs < 0 || rs + rlen > L) continue;
    /* by arithmetic: the read covers [tp+1+dead, tp+1+good) going forward */
    const rdS = strand === 1 ? tp + 1 + dead : tp - good;
    const rdE = strand === 1 ? tp + 1 + good : tp - dead;
    const expect = Math.max(0, Math.min(rdE, rs + rlen) - Math.max(rdS, rs));

    const plan = E.sangerPlan(T, { start: rs, end: rs + rlen }, [primer],
      { circular: circ, dead: dead, good: good });
    const c = plan.candidates[0];
    tested++;
    if (!c) { ok('the planted primer is a candidate', false, JSON.stringify({ L, at, strand, lead })); continue; }
    if (c.covers !== expect) {
      ok('it covers exactly what the arithmetic says', false,
         JSON.stringify({ L, at, strand, lead, rlen, got: c.covers, expect: expect }));
      continue;
    }
    if (c.site.threePrime !== tp) {
      ok('its 3′ end is where it was planted', false, JSON.stringify({ got: c.site.threePrime, tp }));
      continue;
    }
    if (c.strand !== strand) { ok('on the strand it was planted on', false, strand); continue; }
  }
  ok('400 planted primers are each recovered with the right 3′ end, strand and coverage', fail === 0,
     'tested ' + tested);
}

/* ======================================================================= */
section('a region that runs across the origin');

{
  const L = 4000;
  const T = randSeq(L);
  /* a primer just before the origin reading forward across it */
  const at = L - 120, len = 20;
  const primer = { name: 'origin reader', seq: T.slice(at, at + len) };
  const tp = at + len - 1;
  /* region from 3980 round to 400: 420 bases across base 1 */
  const plan = E.sangerPlan(T, { start: 3980, end: 400 }, [primer], { circular: true, dead: 35, good: 700 });
  ok('the region length is measured the long way round the origin', plan.region.length === 420,
     plan.region.length);
  const c = plan.candidates[0];
  /* the read covers [tp+1+35, tp+1+700) = [3936, 4601) -> 3936..3999 and 0..600
     intersected with 3980..3999 and 0..399 = 20 + 400 = 420 */
  ok('a primer reading across the origin covers the whole of it', c && c.covers === 420,
     c ? c.covers : 'no candidate');
  ok('and the plan says the forward strand is complete', plan.coverage.forward === 420,
     plan.coverage.forward);
  ok('while nothing reads the other strand', plan.coverage.reverse === 0, plan.coverage.reverse);
  ok('so no base is confirmed on both strands', plan.coverage.both === 0, plan.coverage.both);
}

/* ======================================================================= */
section('the plan is honest: brute force agrees with it');

{
  /* a template with many planted primers, some of them repeats */
  let bad = 0, checked = 0;
  for (let trial = 0; trial < 60; trial++) {
    const L = 3000 + Math.floor(rnd() * 3000);
    const T = randSeq(L);
    const circ = rnd() < 0.5;
    const prims = [];
    const nP = 6 + Math.floor(rnd() * 14);
    for (let i = 0; i < nP; i++) {
      const len = 18 + Math.floor(rnd() * 6);
      const at = Math.floor(rnd() * (L - len));
      const fp = T.slice(at, at + len);
      if (!/^[ACGT]+$/.test(fp) || fp.length < len) continue;
      prims.push({ name: 'p' + i, seq: rnd() < 0.5 ? fp : E.revcomp(fp) });
    }
    const rs = Math.floor(rnd() * L);
    const rlen = 200 + Math.floor(rnd() * 1500);
    if (!circ && rs + rlen > L) continue;
    const plan = E.sangerPlan(T, { start: rs, end: circ ? mod(rs + rlen, L) : rs + rlen }, prims,
      { circular: circ });
    checked++;

    /* mark, base by base, what the chosen primers really read */
    const mark = strandWanted => {
      const m = new Uint8Array(plan.region.length);
      plan.plan[strandWanted === 1 ? 'forward' : 'reverse'].forEach(c => {
        const w = c.window.good;
        for (let k = 0; k < w.length; k++) {
          const p = circ ? mod(w.start + k, L) : w.start + k;
          const rel = circ ? mod(p - plan.region.start, L) : p - plan.region.start;
          if (rel >= 0 && rel < plan.region.length) m[rel] = 1;
        }
      });
      return m;
    };
    const mf = mark(1), mr = mark(-1);
    let cf = 0, cr = 0, cb = 0;
    for (let i = 0; i < plan.region.length; i++) { if (mf[i]) cf++; if (mr[i]) cr++; if (mf[i] && mr[i]) cb++; }
    if (cf !== plan.coverage.forward) { bad++; console.log('  forward coverage disagrees', cf, plan.coverage.forward); continue; }
    if (cr !== plan.coverage.reverse) { bad++; console.log('  reverse coverage disagrees', cr, plan.coverage.reverse); continue; }
    if (cb !== plan.coverage.both) { bad++; console.log('  both-strand coverage disagrees', cb, plan.coverage.both); continue; }

    /* a gap must be a real gap: no candidate on that strand may read it */
    ['forward', 'reverse'].forEach(k => {
      const strandWanted = k === 'forward' ? 1 : -1;
      plan.gaps[k].forEach(g => {
        for (let rel = g.from; rel < g.to; rel++) {
          const p = circ ? mod(plan.region.start + rel, L) : plan.region.start + rel;
          const reachable = plan.candidates.some(c => {
            if (c.strand !== strandWanted) return false;
            const w = c.window.good;
            const off = circ ? mod(p - w.start, L) : p - w.start;
            return off >= 0 && off < w.length;
          });
          if (reachable) { bad++; console.log('  a base called a gap is in fact readable', k, p); return; }
        }
      });
    });

    /* nothing with more than one site may be in the plan */
    ['forward', 'reverse'].forEach(k => plan.plan[k].forEach(c => {
      if (c.sites > 1) { bad++; console.log('  a repeat-site primer got into the plan', c.name); }
    }));
  }
  ok('over ' + checked + ' random templates the plan’s coverage survives brute force, its gaps are real, and no repeat-site primer is ever planned', bad === 0, bad + ' problems');
}

/* ======================================================================= */
section('primers that should be refused');

{
  const L = 3000;
  const T = randSeq(L);
  /* the same footprint planted twice gives a mixed trace */
  const fp = T.slice(500, 522);
  const T2 = T.slice(0, 1500) + fp + T.slice(1500 + fp.length);
  const plan = E.sangerPlan(T2, { start: 600, end: 1000 }, [{ name: 'twice', seq: fp }], { circular: false });
  ok('a primer with two sites is kept out of the candidates', plan.candidates.length === 0,
     plan.candidates.length);
  ok('and is reported as unusable with the reason', plan.unusable.length >= 1 &&
     /mixed, unreadable/.test(plan.unusable[0].warnings[0] || ''),
     JSON.stringify(plan.unusable.map(u => u.warnings[0])));

  /* a primer whose 3' end is wrong does not prime at all */
  const good = T.slice(800, 822);
  const broken = good.slice(0, 19) + 'AAA';
  const p2 = E.sangerPlan(T, { start: 900, end: 1200 }, [{ name: 'bad 3′ end', seq: broken }], { circular: false });
  ok('a primer with a mismatched 3′ end is not used', p2.candidates.length === 0, p2.candidates.length);

  /* and one that reads away from the region is not offered */
  const away = { name: 'wrong way', seq: E.revcomp(T.slice(800, 822)) };
  const p3 = E.sangerPlan(T, { start: 1500, end: 1800 }, [away], { circular: false });
  ok('a primer pointing away from the region is not offered', p3.candidates.length === 0, p3.candidates.length);

  /* a primer so close that the region starts inside the dead zone is
     offered, but the problem is stated */
  const close = { name: 'too close', seq: T.slice(1000, 1022) };
  const p4 = E.sangerPlan(T, { start: 1030, end: 1400 }, [close], { circular: false, dead: 35 });
  ok('a primer whose dead zone eats into the region is flagged', p4.candidates.length === 1 &&
     p4.candidates[0].warnings.some(w => /unreadable window/.test(w)),
     JSON.stringify(p4.candidates.map(c => c.warnings)));
}

/* ======================================================================= */
section('reading past the confident length');

{
  const L = 4000, T = randSeq(L);
  const p = { name: 'far', seq: T.slice(1000, 1022) };
  const plan = E.sangerPlan(T, { start: 1800, end: 1900 }, [p], { circular: false, dead: 35, good: 700, max: 900 });
  /* 3' end at 1021: confident read ends at 1722, so 1800-1900 is only in
     the ragged tail and must not count as covered */
  ok('a region beyond the confident read is not counted as covered',
     plan.candidates.length === 1 && plan.candidates[0].covers === 0, JSON.stringify(plan.candidates.map(c => c.covers)));
  ok('but it is still reported as within reach', plan.candidates[0].reach === 100,
     plan.candidates[0].reach);
  ok('with the tail explained', plan.candidates[0].warnings.some(w => /tail/.test(w)),
     JSON.stringify(plan.candidates[0].warnings));
  ok('and the plan admits the region is not covered', plan.coverage.forward === 0, plan.coverage.forward);
}

/* ======================================================================= */
section('designing a primer for what the freezer cannot read');

{
  let bad = 0;
  for (let trial = 0; trial < 30; trial++) {
    const L = 3000 + Math.floor(rnd() * 2000);
    const T = randSeq(L);
    const circ = rnd() < 0.5;
    const strand = rnd() < 0.5 ? 1 : -1;
    const target = 600 + Math.floor(rnd() * (L - 1400));
    const d = E.designSeqPrimer(T, target, strand, { circular: circ });
    if (!d) { bad++; console.log('  no primer designed', { L, target, strand }); continue; }
    /* it must be the only place it anneals */
    const sites = E.primerSites(T, d.seq, { circular: circ });
    if (sites.length !== 1) { bad++; console.log('  designed primer is not unique', sites.length); continue; }
    if (sites[0].strand !== strand) { bad++; console.log('  designed primer is on the wrong strand'); continue; }
    /* and it must clear the dead zone before the target */
    if (d.lead < E.SANGER.dead) { bad++; console.log('  designed primer too close', d.lead); continue; }
    /* and actually read the target, by the planner's own account */
    const region = strand === 1 ? { start: target, end: target + 50 } : { start: target - 50, end: target };
    const plan = E.sangerPlan(T, region, [{ name: 'designed', seq: d.seq }], { circular: circ });
    if (!plan.candidates.length || plan.candidates[0].covers !== 50) {
      bad++; console.log('  designed primer does not read its target', plan.candidates.length ? plan.candidates[0].covers : 'none');
    }
    if (d.tm < 45 || d.tm > 70) { bad++; console.log('  designed primer Tm out of range', d.tm); }
  }
  ok('30 designed primers are each unique, on the right strand, clear of the dead zone and actually read the target', bad === 0, bad + ' problems');
}

/* ======================================================================= */
section('pairs that amplify the region');

{
  const L = 4000, T = randSeq(L);
  const f = { name: 'amp f', seq: T.slice(500, 522) };
  const r = { name: 'amp r', seq: E.revcomp(T.slice(1500, 1522)) };
  const near = { name: 'too near', seq: T.slice(980, 1002) };
  const amps = E.sangerAmplicons(T, { start: 1000, end: 1100 }, [f, r, near], { circular: false });
  ok('a pair flanking the region is found', amps.length >= 1, amps.length);
  const a = amps[0];
  ok('the product length is the distance between the outer ends', a.length === 1022, a.length);
  ok('the region sits inside it, clear of both ends',
     a.inset.fromStart >= 55 && a.inset.fromEnd >= 55, JSON.stringify(a.inset));
  /* The engine's own PCR must agree that this pair gives that product.
     E.pcr returns the list of products directly — an array carrying
     `found`, `sites` and `truncated` as extra properties — so there is no
     .products to read. This test asked for one until 2026-11; the engine
     was checked against test.js and test2.js, which have always used the
     array, and against E.pcr itself at engine.js:653. The test was wrong,
     not the engine. */
  const prod = E.pcr(T, f.seq, r.seq, { circular: false });
  ok('the PCR engine agrees on the product', prod.length === 1 &&
     prod[0].length === a.length,
     JSON.stringify(prod.map(p => p.length)));
  /* a primer that would put the region in the dead zone is not paired */
  ok('a pair that leaves the region against the primer is refused',
     !amps.some(x => x.forward.name === 'too near'),
     JSON.stringify(amps.map(x => x.forward.name + '/' + x.reverse.name)));
}

/* ======================================================================= */
section('reading a stock list');

{
  /* a tab-separated export with a header */
  const tsv = ['Name\tSequence\tTm (C)\tWell\tPlate\tNotes',
               'T7 fwd\tTAATACGACTCACTATAGGG\t54\tA01\tBOX1\tstandard',
               'M13 rev\tCAGGAAACAGCTATGAC\t51\tA02\tBOX1\t',
               'junk row\t\t\t\t\t'].join('\n');
  const a = E.parseStockTable(tsv);
  ok('a tab-separated list with a header is read', a.primers.length === 2 && a.mode === 'header',
     JSON.stringify({ n: a.primers.length, mode: a.mode }));
  ok('its names, Tms and wells come through',
     a.primers[0].name === 'T7 fwd' && a.primers[0].tm === 54 && /A01/.test(a.primers[0].loc),
     JSON.stringify(a.primers[0]));
  ok('a row with no sequence is reported, not dropped in silence',
     a.skipped.length === 1 && /no sequence/.test(a.skipped[0].why), JSON.stringify(a.skipped));

  /* a comma-separated one */
  const csv = 'id,description,seq\nOL-1,my primer,ACGTACGTACGTACGTACGT\nOL-2,another,TTTTGGGGCCCCAAAATTTT';
  const b = E.parseStockTable(csv);
  ok('a comma-separated list is read', b.primers.length === 2, b.primers.length);
  ok('and the id column is kept', b.primers[0].id === 'OL-1', JSON.stringify(b.primers[0]));

  /* sequences typed in groups, with a lowercase tail, no header at all */
  const loose = ['OL-101  L3  BOX2  A01  2024-01-02 my gene F  infusion  ttgtatttccagggc ACG TAC GTA CGT ACG TAC  58  100  TE pH 8.0  note here',
                 'OL-102  L3  BOX2  A02  2024-01-03 my gene R  infusion  CAA GCT TCG TCA TCA GGG TTT AAA CCC  60  100  TE pH 8.0'].join('\n');
  const c = E.parseStockTable(loose);
  ok('a list with no header is read by shape', c.primers.length === 2 && c.mode === 'shape',
     JSON.stringify({ n: c.primers.length, mode: c.mode }));
  ok('the spaces people type inside a sequence are closed up',
     c.primers[0].seq === 'TTGTATTTCCAGGGCACGTACGTACGTACGTAC', c.primers[0].seq);
  ok('the description is taken, not the notes column', c.primers[0].name === 'my gene F',
     c.primers[0].name);
  ok('the tube location is assembled from shelf, box and well',
     /L3/.test(c.primers[0].loc) && /BOX2/.test(c.primers[0].loc) && /A01/.test(c.primers[0].loc),
     c.primers[0].loc);
  ok('the stock Tm is read', c.primers[0].tm === 58, c.primers[0].tm);

  /* the field before a sequence must never bleed into it */
  const bleed = 'OL-200  L3  BOX1  B02  2024-02-02 tagged thing  infusion pHTN     tgggcccaaatctagttaAGGCCTCGCCCCACGAAG  65  100  TE';
  const d = E.parseStockTable(bleed);
  ok('the tail of the column before is not swallowed into the sequence',
     d.primers[0].seq === 'TGGGCCCAAATCTAGTTAAGGCCTCGCCCCACGAAG', d.primers[0].seq);

  /* a header that does not line up with its rows must be abandoned */
  const ragged = ['#   Rm   ShelfPlate   PositionOrderID Date   Description   Clone method   Sequence   Tm (C)',
                  'OL-1  224  L3 BOX1  A01  WTL-1 2019-09-25RFWD3:S425F  Infusion  TTG TAT TTC CAG GGC TCA CCC  62  100',
                  'OL-2  224  L3 BOX1  A02  WTL-2 2019-09-26SIAH1:M32F   Infusion  TTG TAT TTC CAG GGC ATG AGC  59  100'].join('\n');
  const e = E.parseStockTable(ragged);
  ok('a header whose columns do not match its rows is abandoned for shape reading',
     e.primers.length === 2 && e.mode === 'shape', JSON.stringify({ n: e.primers.length, mode: e.mode }));
  ok('and the real descriptions are recovered',
     e.primers[0].name === 'RFWD3:S425F' && e.primers[1].name === 'SIAH1:M32F',
     JSON.stringify(e.primers.map(p => p.name)));
}

/* ======================================================================= */
section('grading a primer as a sequencing primer');

{
  const g = E.gradeSeqPrimer('TAATACGACTCACTATAGGG');
  ok('a standard T7 primer grades with no complaint about length', !/18-24/.test(g.notes.join(' ')),
     JSON.stringify(g.notes));
  const short = E.gradeSeqPrimer('ACGTACGTACGTACG');
  ok('a short primer is called short', short.notes.some(n => /18-24/.test(n)), JSON.stringify(short.notes));
  const runny = E.gradeSeqPrimer('ACGTAAAAAAAACGTACGTA');
  ok('a homopolymer run inside the primer is noticed', runny.notes.some(n => /run of/.test(n)),
     JSON.stringify(runny.notes));
  const gc = E.gradeSeqPrimer('GCGCGCGCGCGCGCGCGCGC');
  ok('a GC-heavy primer is penalised', gc.penalty > 0 && gc.notes.length > 0, JSON.stringify(gc.notes));
  ok('and an AT-heavy one too', E.gradeSeqPrimer('ATATATATATATATATATAT').penalty > 0);
}

/* ======================================================================= */
section('a slippery homopolymer between primer and target');

{
  const L = 3000;
  const T = randSeq(800) + 'A'.repeat(14) + randSeq(L - 814);
  const p = { name: 'before the run', seq: T.slice(700, 722) };
  const plan = E.sangerPlan(T, { start: 900, end: 1200 }, [p], { circular: false });
  ok('a long homopolymer in the run-up is called out',
     plan.candidates.length === 1 && plan.candidates[0].warnings.some(w => /slip/.test(w)),
     JSON.stringify(plan.candidates.map(c => c.warnings)));
}

/* ======================================================================= */
section('the settings are honoured');

{
  const L = 4000, T = randSeq(L);
  const p = { name: 'x', seq: T.slice(1000, 1022) };
  const tight = E.sangerPlan(T, { start: 1100, end: 2000 }, [p], { circular: false, dead: 50, good: 400 });
  /* Both settings are distances from the primer's 3' end, which here is the
     base at 1021: the read is unusable for `dead` bases and confident out
     to `good`. So the confident window is 1021+1+50 = 1072 to 1021+1+400 =
     1422, and what it covers of 1100..2000 is 1422-1100 = 322 bases. Worked
     out from the definition in E.sangerWindow rather than from the answer:
     this line asked for 1072-1100+329 until 2026-11, which corresponds to
     no definition of either setting and had gone unchecked while the suite
     could not run at all. */
  const tipAt = 1021, deadZ = 50, goodTo = 400;
  const windowEnds = tipAt + 1 + goodTo;
  ok('a shorter confident read covers less',
     tight.candidates[0].covers === windowEnds - 1100,
     tight.candidates[0].covers + ' vs ' + (windowEnds - 1100));
  ok('and the window itself starts after the dead zone, not at the primer',
     tight.candidates[0].window.good.start === tipAt + 1 + deadZ,
     tight.candidates[0].window.good.start);
  const loose = E.sangerPlan(T, { start: 1100, end: 2000 }, [p], { circular: false, dead: 20, good: 1000 });
  ok('a longer one covers more', loose.candidates[0].covers === 900, loose.candidates[0].covers);
  ok('and the settings used are reported back', loose.settings.good === 1000 && loose.settings.dead === 20,
     JSON.stringify(loose.settings));
}

console.log('\n   ' + pass + ' passed, ' + fail + ' failed\n');
process.exit(fail ? 1 : 0);
