"""Write GFF3, Benchling-shaped JSON, openVectorEditor-shaped JSON and
INSDSeq XML from real GenBank files, using Biopython as the source of truth,
plus a truth file holding each feature's extracted sequence.

The point is that nothing here shares a line of code with PlasmidBench: the
coordinates written out come from Biopython's own parse, and the sequence
each feature covers comes from Biopython's own extract(). cmpfmt.js then
reads the four files back with PlasmidBench's readers and must produce the
same bases for every feature. A reader that is one base out, or that reads
an inclusive end as an exclusive one, fails immediately.

  python3 -I mkfmt.py <genbank-dir> <out-dir>
"""
import json
import os
import sys
import xml.sax.saxutils as sx

from Bio import SeqIO
from Bio.SeqFeature import SimpleLocation


def simple(loc):
    """A single contiguous forward or reverse range, which is what GFF3 and
    the flat JSON shapes can represent without inventing a convention."""
    return isinstance(loc, SimpleLocation) and int(loc.start) < int(loc.end)


def gff_escape(s):
    out = []
    for ch in str(s):
        if ch in ';=&,\t\n%':
            out.append('%%%02X' % ord(ch))
        else:
            out.append(ch)
    return ''.join(out)


def main(src, dst):
    os.makedirs(dst, exist_ok=True)
    made = 0
    for fn in sorted(os.listdir(src)):
        if not fn.endswith(('.gb', '.gbk')):
            continue
        path = os.path.join(src, fn)
        try:
            recs = list(SeqIO.parse(path, 'genbank'))
        except Exception:
            continue
        if len(recs) != 1:
            continue
        r = recs[0]
        try:
            seq = str(r.seq).upper()
        except Exception:
            continue          # CONTIG records carry no bases
        if not seq or set(seq) - set('ACGTURYSWKMBDHVN'):
            continue
        circular = r.annotations.get('topology', '') == 'circular'
        base = os.path.splitext(fn)[0]

        feats = []          # simple ranges: every format can carry these
        allf = []           # every feature, for the XML (location strings)
        for f in r.features:
            if f.type == 'source' or f.location is None:
                continue
            label = (f.qualifiers.get('label') or f.qualifiers.get('gene')
                     or f.qualifiers.get('product') or f.qualifiers.get('locus_tag')
                     or [f.type])[0]
            try:
                fs = str(f.extract(r.seq)).upper()
            except Exception:
                continue
            allf.append({'type': f.type, 'label': label, 'seq': fs,
                         'loc': str(f.location) if hasattr(f, 'location') else ''})
            if simple(f.location):
                feats.append({'type': f.type, 'label': label, 'seq': fs,
                              'start': int(f.location.start), 'end': int(f.location.end),
                              'strand': -1 if f.location.strand == -1 else 1})

        if not feats:
            continue

        # ---- GFF3, with the sequence appended ------------------------------
        rows = ['##gff-version 3',
                '##sequence-region %s 1 %d' % (base, len(seq))]
        if circular:
            rows.append('%s\toracle\tregion\t1\t%d\t.\t+\t.\tID=%s;Is_circular=true'
                        % (base, len(seq), base))
        for i, f in enumerate(feats):
            rows.append('\t'.join([
                base, 'oracle', f['type'], str(f['start'] + 1), str(f['end']), '.',
                '+' if f['strand'] == 1 else '-', '.',
                'ID=f%d;Name=%s' % (i, gff_escape(f['label']))]))
        rows.append('##FASTA')
        rows.append('>' + base)
        for i in range(0, len(seq), 60):
            rows.append(seq[i:i + 60])
        open(os.path.join(dst, base + '.gff3'), 'w').write('\n'.join(rows) + '\n')

        # ---- Benchling: 0-based start inclusive, end exclusive -------------
        bench = {'name': base, 'bases': seq, 'isCircular': circular,
                 'annotations': [{'name': f['label'], 'type': f['type'],
                                  'start': f['start'], 'end': f['end'],
                                  'strand': f['strand']} for f in feats]}
        open(os.path.join(dst, base + '.benchling.json'), 'w').write(json.dumps(bench))

        # ---- openVectorEditor: 0-based, both ends inclusive ---------------
        ove = {'name': base, 'sequence': seq, 'circular': circular,
               'features': [{'name': f['label'], 'type': f['type'],
                             'start': f['start'], 'end': f['end'] - 1,
                             'forward': f['strand'] == 1} for f in feats]}
        open(os.path.join(dst, base + '.ove.json'), 'w').write(json.dumps(ove))

        # ---- INSDSeq XML, with Biopython's own location strings ------------
        x = ['<?xml version="1.0"?>', '<INSDSet>', '<INSDSeq>',
             '<INSDSeq_locus>%s</INSDSeq_locus>' % sx.escape(base),
             '<INSDSeq_topology>%s</INSDSeq_topology>' % ('circular' if circular else 'linear'),
             '<INSDSeq_feature-table>']
        for f in r.features:
            if f.type == 'source' or f.location is None:
                continue
            label = (f.qualifiers.get('label') or f.qualifiers.get('gene')
                     or f.qualifiers.get('product') or f.qualifiers.get('locus_tag')
                     or [f.type])[0]
            x.append('<INSDFeature>')
            x.append('<INSDFeature_key>%s</INSDFeature_key>' % sx.escape(f.type))
            x.append('<INSDFeature_location>%s</INSDFeature_location>' % sx.escape(gb_loc(f)))
            x.append('<INSDFeature_quals><INSDQualifier>'
                     '<INSDQualifier_name>label</INSDQualifier_name>'
                     '<INSDQualifier_value>%s</INSDQualifier_value>'
                     '</INSDQualifier></INSDFeature_quals>' % sx.escape(str(label)))
            x.append('</INSDFeature>')
        x += ['</INSDSeq_feature-table>',
              '<INSDSeq_sequence>%s</INSDSeq_sequence>' % seq.lower(),
              '</INSDSeq>', '</INSDSet>']
        open(os.path.join(dst, base + '.insd.xml'), 'w').write('\n'.join(x))

        truth = {'name': base, 'seq': seq, 'circular': circular,
                 'simple': feats, 'all': allf}
        open(os.path.join(dst, base + '.truth.json'), 'w').write(json.dumps(truth))
        made += 1
    print('  wrote %d record(s) in four formats each' % made)


def gb_loc(f):
    """GenBank location string for a Biopython feature."""
    loc = f.location
    parts = list(loc.parts)
    def one(p):
        # a^b: a zero-width site between two bases, which Biopython stores
        # as start == end and GenBank writes with a caret
        if int(p.start) == int(p.end):
            return '%d^%d' % (int(p.start), int(p.start) + 1)
        return '%d..%d' % (int(p.start) + 1, int(p.end))
    if len(parts) == 1:
        s = one(parts[0])
    else:
        order = parts if (parts[0].strand != -1) else list(reversed(parts))
        s = 'join(%s)' % ','.join(one(p) for p in order)
    return 'complement(%s)' % s if loc.strand == -1 else s


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
