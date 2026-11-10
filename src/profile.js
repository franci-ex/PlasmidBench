/* ============================================================================
   PlasmidBench — the lab profile.

   Everything a lab has to set before this tool fits the way they work, in
   one named bundle that can be exported, handed to another lab, and
   imported there. The reason it exists: this page's distinctive claim is
   that it knows what a particular freezer holds and what a particular bench
   runs on. That only works if a second lab can make it theirs without
   editing the source, and if a new student can be handed a file rather than
   a list of settings to re-enter.

   WHAT IS IN IT
     who          lab name, PI, contact — used on exported plans and figures
     bench        salt, magnesium, dNTP and oligo concentrations; the
                  polymerase and kit a plan is written for
     enzymes      which suppliers the lab can actually buy from, and any
                  enzymes added by hand
     sequencing   the read length and dead zone the facility really gives
     display      map style, bases per line, codon table, theme
     records      the eLabFTW address (never the key)
     stock        the column mapping that makes this lab's own spreadsheet
                  import without being reformatted

   WHAT IS NOT
     No credentials. Not the eLabFTW API key, not anything else secret: a
     profile is made to be emailed, and a file made to be emailed is the
     wrong place for a key. The exporter cannot reach one — the key is kept
     outside the settings the profile reads, in elab.js.

     No records. Not plasmids, primers, strains or stock lists. Those are
     the lab's data, they are large, and they belong in eLabFTW or in an
     export of their own. A profile is configuration, and keeping the two
     apart is what makes it safe to share.

   Both exclusions are enforced in code below, not left to the author of
   the next field to remember.
   ============================================================================ */
(function (root) {
  'use strict';
  /** @type {any} */
  const PB = root.PB = root.PB || {};

  /* The settings a profile carries, grouped as a person thinks of them.
     Each entry: the setting key, and a line saying what it is for — which
     is also what the panel shows, so there is one description and not two. */
  /** @type {Array<[string, Array<[string, string]>]>} group, then its [key, what it is for] rows */
  const FIELDS = [
    ['who', [
      ['lab.name', 'The lab, as it should appear on a plan or a figure'],
      ['lab.pi', 'Whose lab it is'],
      ['lab.contact', 'Who to ask about this profile']
    ]],
    ['bench', [
      ['tmCond', 'Salt, magnesium, dNTP and oligo concentrations — every melting temperature in the page is calculated at these'],
      ['bench.polymerase', 'The polymerase a cycling plan is written for'],
      ['bench.kit', 'The assembly kit, which sets the homology arm length'],
      ['bench.armLen', 'Homology arm length, if the kit is not the one deciding']
    ]],
    ['enzymes', [
      ['enzSet', 'Which enzymes the panel offers'],
      ['suppliers', 'The suppliers this lab can buy from — the enzyme list can be held to these'],
      ['customEnzymes', 'Enzymes added by hand, with their sites and cut positions']
    ]],
    ['sequencing', [
      ['sanger.good', 'How far a read from your facility is actually confident, in bases'],
      ['sanger.dead', 'The unusable stretch right after the primer'],
      ['sanger.facility', 'Who does the sequencing']
    ]],
    ['display', [
      ['theme', 'Light or dark'],
      ['mapMode', 'How the map is drawn'],
      ['seqOpts', 'Bases per line, and what is shown under them'],
      ['codeTable', 'Which genetic code is the default']
    ]],
    ['records', [
      ['elab.url', 'The lab’s eLabFTW address. The key is never part of a profile.']
    ]],
    ['stock', [
      ['stock.columns', 'How this lab’s oligo spreadsheet is laid out, so it imports as it is'],
      ['plasmid.columns', 'The same for the plasmid collection'],
      ['cell.columns', 'And for the strain list']
    ]]
  ];

  /* Keys that must never leave, whatever anyone adds to FIELDS later. A
     profile is made to be emailed. */
  const NEVER_EXPORT = [
    /key/i, /token/i, /secret/i, /password/i, /auth/i,
    /^stock$/, /^primers$/, /^plasmidStock$/, /^cellStock$/, /^library/, /^elab$/
  ];
  const exportable = k => !NEVER_EXPORT.some(re => re.test(k));

  function allKeys() {
    const out = [];
    FIELDS.forEach(([, rows]) => rows.forEach(([k]) => out.push(k)));
    return out;
  }

  /** The profile as it stands, ready to be written to a file. */
  function current() {
    const values = {};
    allKeys().forEach(k => {
      if (!exportable(k)) return;                /* belt and braces */
      const v = PB.lsGet(k, undefined);
      if (v !== undefined && v !== null && v !== '') values[k] = v;
    });
    /* the eLabFTW address lives inside the elab settings, which as a whole
       are not exportable because the key may be in there; the address is
       lifted out on its own */
    const el = PB.lsGet('elab', {});
    if (el && el.url) values['elab.url'] = el.url;
    return {
      format: 'plasmidbench-lab-profile',
      version: 1,
      name: PB.lsGet('lab.name', '') || 'Unnamed lab',
      exported: new Date().toISOString().slice(0, 10),
      note: 'Configuration only. No API keys, no plasmids, no primer or strain lists.',
      values: values
    };
  }

  /**
   * Read a profile file. Returns what would change rather than changing it,
   * so a person can see what they are about to take on from another lab —
   * importing settings blind is how someone ends up calculating every
   * melting temperature at the wrong salt.
   */
  function inspect(text) {
    let j;
    try { j = typeof text === 'string' ? JSON.parse(text) : text; }
    catch (e) { return { ok: false, error: 'That is not a profile file — it is not even JSON.' }; }
    if (!j || j.format !== 'plasmidbench-lab-profile') {
      return { ok: false, error: 'That JSON is not a PlasmidBench lab profile. A profile says ' +
               'so in its "format" field.' };
    }
    const known = new Set(allKeys());
    const changes = [], ignored = [], refused = [];
    Object.keys(j.values || {}).forEach(k => {
      if (!exportable(k)) { refused.push(k); return; }
      if (!known.has(k)) { ignored.push(k); return; }
      const now = PB.lsGet(k, undefined);
      const next = j.values[k];
      if (JSON.stringify(now) !== JSON.stringify(next))
        changes.push({ key: k, from: now, to: next, what: describe(k) });
    });
    return { ok: true, name: j.name || 'Unnamed lab', exported: j.exported || '',
             changes: changes, ignored: ignored, refused: refused };
  }

  function describe(key) {
    for (const [, rows] of FIELDS) for (const [k, what] of rows) if (k === key) return what;
    return key;
  }

  /** Apply what inspect() reported. Only the keys it listed. */
  function apply(report) {
    if (!report || !report.ok) return 0;
    let n = 0;
    report.changes.forEach(c => {
      if (c.key === 'elab.url') {
        /* folded back into the elab settings, leaving any key alone */
        if (PB.elab) PB.elab.setConfig({ url: c.to });
      } else {
        PB.lsSet(c.key, c.to);
      }
      n++;
    });
    /* the settings the running page already read have to be re-read */
    if (PB.reloadSettings) PB.reloadSettings();
    if (PB.renderAll) PB.renderAll();
    return n;
  }

  PB.profile = {
    FIELDS: FIELDS,
    current: current,
    inspect: inspect,
    apply: apply,
    describe: describe,
    exportable: exportable,
    /** the file a person downloads */
    toText: () => JSON.stringify(current(), null, 2),
    filename: () => (PB.lsGet('lab.name', 'lab') || 'lab')
      .replace(/[^\w.-]+/g, '_').toLowerCase() + '-plasmidbench-profile.json'
  };

  if (typeof module !== 'undefined' && module.exports) module.exports = PB.profile;
})(typeof window !== 'undefined' ? window : globalThis);
