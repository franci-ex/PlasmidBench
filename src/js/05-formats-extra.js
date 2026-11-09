
/* ===========================================================================
   THE REST OF THE SEQUENCE FILE FORMATS

   GenBank, EMBL, FASTA, SnapGene .dna and .ab1 are read elsewhere. This
   file adds the three families a plasmid can also arrive in:

     GFF3        annotations as tab-separated rows, with the sequence
                 either appended after a ##FASTA line or supplied
                 separately. Used by Prokka, bakta, Geneious and every
                 genome browser.

     JSON        what the web-based tools export. There is no single
                 schema, so four real ones are read: Benchling's
                 (`bases` + `annotations`), the Open Vector Format
                 (`sequence` + `features` with `rangeStart`/`rangeEnd`),
                 the jbei/openVectorEditor shape (`sequence.sequence`),
                 and PlasmidBench's own backup file.

     XML         NCBI's two: INSDSeq (the XML spelling of a GenBank
                 record, feature table and all) and TinySeq (sequence
                 only). Both come straight out of efetch.

   Every reader returns the same record shape the GenBank reader returns,
   so the importer, the map, the exporters and the tests do not have to
   know which door a construct came in through:

     { name, description, seq, features[], topology, molecule,
       dropped, droppedChars, hadU, skippedFeatures? }

   A feature is {type, start, end, strand, label, quals} in 0-based
   half-open coordinates, the same convention the rest of the engine uses.

   What this file deliberately does NOT do is guess. A row whose
   coordinates fall outside the sequence is reported in skippedFeatures
   rather than clamped to fit, because a feature silently moved to the
   nearest legal position is worse than one that is missing: you would
   design primers against it.
   ========================================================================= */
(function (root) {
  'use strict';
  const E = root.Bio;

  const nl = s => String(s == null ? '' : s).replace(/\r\n?/g, '\n');

  /* A feature arriving from GFF, JSON or XML, in whatever coordinate
     convention that format uses, turned into ours and checked. Returns
     null when it cannot be placed, so the caller can report it. */
  function place(o, L, circular) {
    let s = Math.round(o.start), e = Math.round(o.end);
    if (!isFinite(s) || !isFinite(e)) return null;
    if (s < 0 || e < 0) return null;
    if (e <= s) {
      /* end before start is a wrap on a circle, and nothing anywhere else */
      if (!circular || !L) return null;
      if (s >= L || e > L) return null;
      return { start: s, end: e, strand: o.strand === -1 ? -1 : 1 };
    }
    if (L && (s >= L || e > L)) return null;
    return { start: s, end: e, strand: o.strand === -1 ? -1 : 1 };
  }

  function mkFeature(o, L, circular) {
    const p = place(o, L, circular);
    if (!p) return null;
    const label = String(o.label || o.name || o.type || 'feature').trim() || 'feature';
    const f = {
      type: String(o.type || 'misc_feature').trim() || 'misc_feature',
      start: p.start, end: p.end, strand: p.strand,
      label: label.length > 80 ? label.slice(0, 77) + '...' : label,
      quals: o.quals || {}
    };
    if (o.color) f.color = o.color;
    if (o.note) f.note = String(o.note);
    if (Array.isArray(o.parts) && o.parts.length > 1) {
      const parts = [];
      for (const q of o.parts) {
        const pp = place({ start: q.start, end: q.end, strand: q.strand }, L, circular);
        if (!pp) return null;
        parts.push({ start: pp.start, end: pp.end });
      }
      f.parts = parts;
    }
    return f;
  }

  function finish(rec) {
    const rep = E.seqTextReport(rec.seq || '');
    rec.seq = rep.seq;
    rec.dropped = rep.dropped; rec.droppedChars = rep.droppedChars; rec.hadU = rep.hadU;
    if (rep.looksProtein) rec.molecule = 'protein';
    rec.molecule = rec.molecule || 'DNA';
    rec.topology = rec.topology === 'circular' ? 'circular' : 'linear';
    rec.features = rec.features || [];
    rec.description = rec.description || '';
    return rec;
  }

  /* =======================================================================
     GFF3

     Nine tab-separated columns; column 4 and 5 are 1-based inclusive, so
     start-1 is ours. Column 7 is + / - / . — a dot means unstranded, which
     we draw forward. Column 9 holds attributes as key=value;key=value with
     percent-encoding, and that is where the name lives: Name, then gene,
     then locus_tag, then product, then ID.

     The sequence may be in the same file after a `##FASTA` line (Prokka
     and bakta both write it that way), or the file may be annotations
     only, in which case the caller is told so and can apply them to a
     construct already open.

     `##sequence-region` lines give each sequence's extent, which is how a
     multi-contig file is split. Circularity is not part of GFF3, but
     bakta and Prokka both write `;Is_circular=true` on the region feature,
     and that is honoured when present.
     ===================================================================== */
  const GFF_SKIP = { region: 1, source: 1, databank_entry: 1, chromosome: 1, contig: 1, supercontig: 1, scaffold: 1 };

  function gffDecode(s) {
    return String(s).replace(/%([0-9A-Fa-f]{2})/g, (m, h) => String.fromCharCode(parseInt(h, 16)));
  }
  function gffAttrs(col9) {
    const out = {};
    String(col9 || '').split(';').forEach(kv => {
      if (!kv.trim()) return;
      const i = kv.indexOf('=');
      if (i === -1) { if (out._ === undefined) out._ = gffDecode(kv.trim()); return; }
      const k = kv.slice(0, i).trim(), v = gffDecode(kv.slice(i + 1).trim());
      if (out[k] === undefined) out[k] = v;
    });
    return out;
  }

  E.parseGFF = function (text, opts) {
    opts = opts || {};
    text = nl(text);
    /* the sequence, if the file carries one */
    const fastaAt = text.search(/^##FASTA[ \t]*$/m);
    let seqs = {};
    let body = text;
    if (fastaAt !== -1) {
      body = text.slice(0, fastaAt);
      E.parseFASTA(text.slice(fastaAt).replace(/^##FASTA[^\n]*\n?/, '')).forEach(r => {
        seqs[r.name] = r;
      });
    }
    /* ##sequence-region <id> <start> <end> gives the length even when the
       bases are not in the file */
    const regions = {};
    (body.match(/^##sequence-region[^\n]*/gm) || []).forEach(l => {
      const m = /^##sequence-region\s+(\S+)\s+(\d+)\s+(\d+)/.exec(l);
      if (m) regions[m[1]] = +m[3];
    });

    const order = [];
    const bySeq = {};
    const circularOf = {};
    body.split('\n').forEach(line => {
      if (!line || line[0] === '#') return;
      const c = line.split('\t');
      if (c.length < 8) return;                       /* not a feature row */
      const id = c[0];
      const a = +c[3], b = +c[4];
      if (!isFinite(a) || !isFinite(b)) return;
      const type = (c[2] || 'misc_feature').trim();
      const at = gffAttrs(c[8]);
      if (/^(true|1)$/i.test(at.Is_circular || '')) circularOf[id] = true;
      if (GFF_SKIP[type.toLowerCase()]) return;       /* the contig itself, not a feature */
      if (!bySeq[id]) { bySeq[id] = []; order.push(id); }
      bySeq[id].push({
        type: type,
        start: a - 1, end: b,
        strand: c[6] === '-' ? -1 : 1,
        label: at.Name || at.gene || at.locus_tag || at.product || at.ID || at._ || type,
        note: at.product || at.note || '',
        quals: at
      });
    });

    /* a sequence-region or a FASTA block for an id nobody annotated is still
       a sequence worth importing */
    Object.keys(seqs).forEach(id => { if (order.indexOf(id) === -1) order.push(id); });
    Object.keys(regions).forEach(id => { if (order.indexOf(id) === -1) order.push(id); });

    const out = [];
    order.forEach(id => {
      const fa = seqs[id];
      const seq = fa ? fa.seq : (opts.seq && order.length === 1 ? opts.seq : '');
      const L = seq.length || regions[id] || 0;
      const circular = !!circularOf[id] || (opts.topology === 'circular' && !seq);
      const rec = {
        name: id, description: fa ? fa.description : '',
        seq: seq, topology: circular ? 'circular' : 'linear',
        molecule: 'DNA', features: [], skippedFeatures: []
      };
      (bySeq[id] || []).forEach(raw => {
        const f = mkFeature(raw, L, circular);
        if (f) rec.features.push(f);
        else rec.skippedFeatures.push(raw.type + ' ' + (raw.start + 1) + '..' + raw.end);
      });
      if (!rec.skippedFeatures.length) delete rec.skippedFeatures;
      out.push(finish(rec));
    });
    /* annotations with nowhere to put them: say so rather than inventing a
       sequence of Ns, which would look like a construct and is not one */
    out.annotationsOnly = out.every(r => !r.seq.length) && out.some(r => r.features.length);
    return out;
  };

  /* =======================================================================
     JSON

     There is no one JSON sequence format, and the shapes that exist
     disagree about what `end` means. Getting that wrong moves every
     feature by one base, which still looks right on a map and is wrong in
     a primer, so the convention is decided from the file's own shape and
     then reported, rather than assumed:

       Benchling      `bases` + `annotations`. "start: 0-based inclusive
                      start index. end: 0-based exclusive end index. The
                      end of the sequence is always represented as 0."
                      (Benchling DnaSequence schema.) So: half-open, and
                      an end of 0 means the end of the sequence.

       openVectorEditor / Teselagen
                      `sequence` + `features`, "start and end are 0-based
                      inclusive for all annotations", strand as
                      `forward: true|false`. So: inclusive, +1 on the end.

       rangeStart / rangeEnd
                      the naming the Open Vector Format and several
                      exporters use for an inclusive range. Inclusive.

       PlasmidBench   its own backup: `records`, already half-open, with
                      ids and history to restore.

     A file that fits none of them is refused with a message naming what is
     read, instead of being guessed at.
     ===================================================================== */
  function jsonSeqOf(o) {
    if (typeof o === 'string') return o;
    if (!o || typeof o !== 'object') return '';
    /* `sequence` is a string in some shapes and a nested object in others,
       so only a string counts here; the nested case is unwrapped by the
       caller. String(anObject) would otherwise yield "[object Object]",
       which cleanSeq would happily turn into the bases CTBJCT. */
    for (const k of ['sequence', 'bases', 'seq']) if (typeof o[k] === 'string') return o[k];
    return '';
  }
  function jsonCircular(o) {
    if (o == null) return false;
    if (o.isCircular !== undefined) return !!o.isCircular;
    if (o.circular !== undefined) return !!o.circular;
    if (o.isSingleStranded !== undefined && o.circular !== undefined) return !!o.circular;
    if (typeof o.topology === 'string') return /circular/i.test(o.topology);
    return false;
  }
  function jsonStrand(a) {
    if (a.strand !== undefined) {
      if (a.strand === -1 || a.strand === '-' || a.strand === 'reverse') return -1;
      return 1;
    }
    if (a.forward !== undefined) return a.forward ? 1 : -1;
    if (a.isForward !== undefined) return a.isForward ? 1 : -1;
    if (a.direction !== undefined) return /rev|-|left/i.test(String(a.direction)) ? -1 : 1;
    return 1;
  }
  /* Which of the four shapes is this? Decided per record, because a file
     can hold a list of them. */
  function jsonDialect(o, anns) {
    if (o && o.__plasmidbench) return 'plasmidbench';
    if (o && typeof o.bases === 'string') return 'benchling';
    if (o && Array.isArray(o.annotations) && !Array.isArray(o.features)) return 'benchling';
    const probe = (anns || []).find(a => a && typeof a === 'object') || {};
    if (probe.rangeStart !== undefined || probe.rangeEnd !== undefined) return 'inclusive';
    if (probe.forward !== undefined || probe.isForward !== undefined) return 'ove';
    if (Array.isArray(probe.locations)) return 'inclusive';
    return 'ove';            /* `sequence` + `features` with start/end */
  }
  /* The named dialect's range turned into ours (0-based half-open), or null.
     `L` lets Benchling's "end of the sequence is always represented as 0"
     be honoured rather than read as a zero-length feature at the origin. */
  function jsonRange(a, dialect, L) {
    const pick = (...keys) => {
      for (const k of keys) if (a[k] !== undefined && a[k] !== null && a[k] !== '') return +a[k];
      return undefined;
    };
    const s = pick('start', 'rangeStart', 'begin', 'from', 'startIndex', 'start_index');
    let e = pick('end', 'rangeEnd', 'to', 'endIndex', 'end_index');
    if (s === undefined || e === undefined || !isFinite(s) || !isFinite(e)) return null;
    if (dialect === 'plasmidbench') return { start: s, end: e };
    if (dialect === 'benchling') {
      if (e === 0 && L) e = L;            /* their marker for the last base */
      return { start: s, end: e };
    }
    return { start: s, end: e + 1 };      /* ove and the rangeStart/rangeEnd shapes */
  }

  E.parseSeqJSON = function (text) {
    let data;
    try { data = typeof text === 'string' ? JSON.parse(text) : text; }
    catch (e) { const err = new Error('That file is not valid JSON: ' + e.message); err.code = 'badjson'; throw err; }
    const out = [];
    const dialects = {};

    const take = o => {
      if (!o || typeof o !== 'object') return;
      const seq = jsonSeqOf(o) || jsonSeqOf(o.sequence) || jsonSeqOf(o.seq);
      if (!seq) return;
      const circular = jsonCircular(o) || jsonCircular(o.sequence);
      const src = o.sequence && typeof o.sequence === 'object' ? o.sequence : o;
      const rec = {
        name: String(o.name || src.name || o.id || 'imported'),
        description: String(o.description || src.description || o.notes || ''),
        seq: seq, topology: circular ? 'circular' : 'linear',
        molecule: 'DNA', features: [], skippedFeatures: []
      };
      const L = E.cleanSeq(seq).length;
      let anns = o.annotations || src.annotations || o.features || src.features || [];
      if (anns && !Array.isArray(anns) && typeof anns === 'object') anns = Object.keys(anns).map(k => anns[k]);
      const dialect = jsonDialect(o, anns);
      dialects[dialect] = (dialects[dialect] || 0) + 1;
      rec.dialect = dialect;
      (anns || []).forEach(a => {
        if (!a || typeof a !== 'object') return;
        /* several exporters put the pieces of a joined feature in locations[] */
        let parts = null;
        if (Array.isArray(a.locations) && a.locations.length) {
          parts = a.locations.map(l => jsonRange(l, dialect, L));
          if (parts.some(p => !p)) parts = null;
        }
        const r = parts ? { start: parts[0].start, end: parts[parts.length - 1].end }
                        : jsonRange(a, dialect, L);
        if (!r) { rec.skippedFeatures.push(String(a.name || a.label || a.type || 'feature') + ' (no coordinates)'); return; }
        const f = mkFeature({
          type: a.type || a.featureType || a.kind || 'misc_feature',
          start: r.start, end: r.end, strand: jsonStrand(a),
          label: a.name || a.label || a.gene || a.type,
          color: a.color || a.colour || null,
          note: a.notes || a.note || '',
          parts: parts && parts.length > 1 ? parts : null,
          quals: a.qualifiers || a.quals || {}
        }, L, circular);
        if (f) rec.features.push(f);
        else rec.skippedFeatures.push(String(a.name || a.type || 'feature') + ' ' + (r.start + 1) + '..' + r.end);
      });
      if (!rec.skippedFeatures.length) delete rec.skippedFeatures;
      out.push(finish(rec));
    };

    if (Array.isArray(data)) data.forEach(take);
    else if (data && Array.isArray(data.records)) data.records.forEach(r => take(Object.assign({ __plasmidbench: true }, r)));
    else if (data && Array.isArray(data.sequences)) data.sequences.forEach(take);
    else take(data);

    if (!out.length) {
      const err = new Error('No sequence in that JSON. PlasmidBench reads Benchling exports (bases + annotations), ' +
        'openVectorEditor files (sequence + features), files that use rangeStart/rangeEnd, and its own backups.');
      err.code = 'noseq';
      throw err;
    }
    out.dialect = Object.keys(dialects).sort((a, b) => dialects[b] - dialects[a])[0] || null;
    return out;
  };
  E.JSON_DIALECT_NOTE = {
    benchling:    'read as a Benchling export: starts 0-based, ends exclusive',
    ove:          'read as an openVectorEditor file: starts and ends 0-based and inclusive',
    inclusive:    'read with inclusive ends (rangeStart / rangeEnd)',
    plasmidbench: 'read as a PlasmidBench backup'
  };

  /* =======================================================================
     XML

     INSDSeq is the XML form of a GenBank record and keeps the whole feature
     table, including the location strings, which the existing verified
     location parser then reads. TinySeq carries the sequence only.

     DOMParser is used in the browser; in Node the tests run against the
     same code through a tiny tag scanner, because the only thing needed
     here is well-known element text, not a general XML tree.
     ===================================================================== */
  function tagList(xml, tag) {
    const out = [];
    const re = new RegExp('<' + tag + '(?:\\s[^>]*)?>([\\s\\S]*?)<\\/' + tag + '>', 'g');
    let m;
    while ((m = re.exec(xml))) out.push(m[1]);
    return out;
  }
  function tagText(xml, tag) {
    const m = new RegExp('<' + tag + '(?:\\s[^>]*)?>([\\s\\S]*?)<\\/' + tag + '>').exec(xml);
    return m ? unent(m[1]) : '';
  }
  function unent(s) {
    return String(s).replace(/<!\[CDATA\[([\s\S]*?)\]\]>/g, '$1')
      .replace(/&lt;/g, '<').replace(/&gt;/g, '>').replace(/&quot;/g, '"')
      .replace(/&apos;/g, "'").replace(/&#(\d+);/g, (m, d) => String.fromCharCode(+d))
      .replace(/&amp;/g, '&');
  }

  E.parseSeqXML = function (text) {
    const xml = nl(text);
    const out = [];

    /* ---- INSDSeq ---- */
    tagList(xml, 'INSDSeq').forEach(block => {
      const rec = {
        name: tagText(block, 'INSDSeq_locus') || tagText(block, 'INSDSeq_accession-version') || 'INSDSeq',
        description: tagText(block, 'INSDSeq_definition'),
        seq: tagText(block, 'INSDSeq_sequence'),
        topology: /circular/i.test(tagText(block, 'INSDSeq_topology')) ? 'circular' : 'linear',
        molecule: /aa|protein/i.test(tagText(block, 'INSDSeq_moltype')) ? 'protein' : 'DNA',
        features: [], skippedFeatures: []
      };
      const L = E.cleanSeq(rec.seq).length;
      const circular = rec.topology === 'circular';
      tagList(block, 'INSDFeature').forEach(fb => {
        const type = tagText(fb, 'INSDFeature_key') || 'misc_feature';
        const locStr = tagText(fb, 'INSDFeature_location');
        const quals = {};
        tagList(fb, 'INSDQualifier').forEach(qb => {
          const k = tagText(qb, 'INSDQualifier_name');
          if (k && quals[k] === undefined) quals[k] = tagText(qb, 'INSDQualifier_value') || true;
        });
        /* the same verified location parser GenBank uses */
        const loc = locStr ? E.parseLocation(locStr, L, circular) : null;
        if (!loc) { rec.skippedFeatures.push(type + ' ' + locStr); return; }
        const pick = k => (typeof quals[k] === 'string' && quals[k]) ? quals[k] : null;
        const f = {
          type: type, start: loc.start, end: loc.end, strand: loc.strand,
          label: pick('label') || pick('gene') || pick('product') || pick('locus_tag') ||
                 pick('standard_name') || pick('note') || type,
          quals: quals
        };
        if (loc.parts) f.parts = loc.parts;
        if (loc.between) f.between = true;
        if (loc.mixed) f.mixed = true;
        if (f.label.length > 80) f.label = f.label.slice(0, 77) + '...';
        rec.features.push(f);
      });
      if (!rec.skippedFeatures.length) delete rec.skippedFeatures;
      out.push(finish(rec));
    });

    /* ---- TinySeq ---- */
    if (!out.length) tagList(xml, 'TSeq').forEach(block => {
      out.push(finish({
        name: tagText(block, 'TSeq_accver') || tagText(block, 'TSeq_sid') || 'TSeq',
        description: tagText(block, 'TSeq_defline'),
        seq: tagText(block, 'TSeq_sequence'),
        topology: 'linear',
        molecule: /protein/i.test(tagText(block, 'TSeq_seqtype')) ? 'protein' : 'DNA',
        features: []
      }));
    });

    /* ---- a generic <sequence> element, which several older web tools wrote ---- */
    if (!out.length) {
      const bare = tagText(xml, 'sequence') || tagText(xml, 'Sequence') || tagText(xml, 'seq');
      if (bare && E.cleanSeq(bare).length >= 10) {
        out.push(finish({
          name: tagText(xml, 'name') || tagText(xml, 'id') || 'imported',
          description: tagText(xml, 'description'),
          seq: bare,
          topology: /circular/i.test(tagText(xml, 'topology')) ? 'circular' : 'linear',
          molecule: 'DNA', features: []
        }));
      }
    }

    if (!out.length) {
      const err = new Error('No sequence in that XML. PlasmidBench reads NCBI INSDSeq and TinySeq, ' +
        'which are what efetch returns with rettype=gb&retmode=xml and rettype=fasta&retmode=xml.');
      err.code = 'noseq';
      throw err;
    }
    return out;
  };

  /* =======================================================================
     WHICH FORMAT IS THIS?

     Decided from the content, not the extension, because the extension lies
     constantly: ApE writes GenBank into .ape, Vector NTI writes GenBank into
     .gb, people rename .gbk to .txt to email it, and a .json may be any of
     four schemas. The extension is used only to break a tie and to give a
     better message when nothing matches.
     ===================================================================== */
  E.sniffFormat = function (text, fileName) {
    const t = nl(String(text || ''));
    const head = t.slice(0, 4000);
    const ext = (/\.([A-Za-z0-9_]+)$/.exec(String(fileName || '')) || [, ''])[1].toLowerCase();

    if (/^\s*LOCUS\s+\S/mi.test(head) || /^ORIGIN[ \t]*$/mi.test(t)) return 'genbank';
    if (/^ID\s{3}\S/m.test(head) && /^SQ\s{3}/m.test(t)) return 'embl';
    if (/^ID\s+\S+;/m.test(head) && /^SQ\s/m.test(t)) return 'embl';
    if (/^##gff-version/m.test(head)) return 'gff';
    if (/^\s*<\?xml/.test(head) || /<INSDSeq[\s>]/.test(head) || /<TSeq[\s>]/.test(head)) return 'xml';
    if (/^\s*[[{]/.test(head)) return 'json';
    if (/^\s*>/.test(t)) return 'fasta';
    /* a GFF without its header line: nine tab-separated columns, column 4
       and 5 numeric, column 7 a strand character */
    const rows = t.split('\n').filter(l => l && l[0] !== '#').slice(0, 20);
    if (rows.length && rows.filter(l => {
      const c = l.split('\t');
      return c.length >= 8 && /^\d+$/.test(c[3]) && /^\d+$/.test(c[4]) && /^[-+.?]$/.test(c[6]);
    }).length >= Math.max(1, Math.ceil(rows.length * 0.6))) return 'gff';
    if (ext === 'gff' || ext === 'gff3') return 'gff';
    if (ext === 'json') return 'json';
    if (ext === 'xml') return 'xml';
    if (E.cleanSeq(t).length >= 10) return 'raw';
    return null;
  };

  /* One door for text of any of the formats above. Throws with a message
     worth showing when the text is recognisable but unusable; returns an
     empty array when there is simply nothing in it. */
  E.parseAnySequence = function (text, fileName, opts) {
    const kind = E.sniffFormat(text, fileName);
    switch (kind) {
      case 'genbank': return E.tag(E.parseGenBank(text), 'GenBank');
      case 'embl':    return E.tag(E.parseEMBL(text), 'EMBL');
      case 'fasta':   return E.tag(E.parseFASTA(text), 'FASTA');
      case 'gff':     return E.tag(E.parseGFF(text, opts), 'GFF3');
      case 'json': {
        const r = E.parseSeqJSON(text);
        return E.tag(r, 'JSON' + (E.JSON_DIALECT_NOTE[r.dialect] ? ' (' + E.JSON_DIALECT_NOTE[r.dialect] + ')' : ''));
      }
      case 'xml':     return E.tag(E.parseSeqXML(text), 'XML');
      case 'raw': {
        const rep = E.seqTextReport(text);
        if (rep.seq.length < 10) return E.tag([], 'raw sequence');
        return E.tag([{ name: 'pasted sequence', description: '', seq: rep.seq, topology: 'linear',
          features: [], molecule: rep.looksProtein ? 'protein' : 'DNA',
          dropped: rep.dropped, droppedChars: rep.droppedChars, hadU: rep.hadU }], 'raw sequence');
      }
      default: return E.tag([], null);
    }
  };
  E.tag = function (recs, format) { recs.format = format; return recs; };

  /* Binary formats that are somebody else's private file layout. Saying
     plainly that a file cannot be read, and what to export instead, is
     worth more than reading half of it wrongly. */
  E.UNREADABLE = {
    clc:  'CLC Main Workbench files are a private binary format. In CLC, File → Export → GenBank gives a .gb that holds the sequence and every annotation.',
    sbd:  'Sequencher project files are a private binary format. In Sequencher, File → Export → Sequence → GenBank gives a readable file.',
    ma4:  'Vector NTI archive files are a private binary format. In Vector NTI, File → Export → GenBank gives a .gb.',
    xdna: 'Serial Cloner .xdna files are a private binary format. In Serial Cloner, File → Export → GenBank gives a readable file.',
    sqn:  'Sequin submission files hold ASN.1, not a sequence file. Converting with NCBI asn2gb, or exporting GenBank from the tool that made it, gives a file PlasmidBench reads.',
    svg:  'An SVG is a picture of a map, not a sequence — it holds no bases to work on. Export GenBank from whatever drew it.',
    png:  'That is a picture of a map, not a sequence. Export GenBank from whatever drew it.',
    jpg:  'That is a picture of a map, not a sequence. Export GenBank from whatever drew it.',
    jpeg: 'That is a picture of a map, not a sequence. Export GenBank from whatever drew it.',
    pdf:  'A PDF holds a picture of a map, not a sequence. Export GenBank from whatever drew it.',
    docx: 'A Word file is not a sequence file. Paste the bases, or the GenBank text, into the import box instead.',
    xlsx: 'A spreadsheet is not a sequence file. Paste the bases into the import box instead.',
    zip:  'That is an archive. Unzip it and import the sequence files inside.',
    scf:  'SCF traces are not supported. Export the read as .ab1, or as FASTA, from your sequencing provider.'
  };

  if (typeof module !== 'undefined' && module.exports) module.exports = E;
  root.Bio = E;
})(typeof window !== 'undefined' ? window : globalThis);

