
/* ======================= A SPREADSHEET, WRITTEN HERE =======================
   An .xlsx file is a zip of XML parts. Everything in this tool runs from
   one file with nothing fetched, so the zip is written by hand: stored
   entries, no compression, which needs only a CRC and keeps the writer
   small enough to read. The result opens in Excel, Numbers, LibreOffice
   and Google Sheets as an ordinary workbook with typed cells, a frozen
   heading row and columns wide enough to read.
============================================================================ */
(function (PB) {
  'use strict';

  /* ---- CRC-32, the one thing a stored zip entry still needs ---- */
  let CRC_TABLE = null;
  function crcTable() {
    if (CRC_TABLE) return CRC_TABLE;
    const t = new Int32Array(256);
    for (let n = 0; n < 256; n++) {
      let c = n;
      for (let k = 0; k < 8; k++) c = (c & 1) ? (0xEDB88320 ^ (c >>> 1)) : (c >>> 1);
      t[n] = c;
    }
    CRC_TABLE = t;
    return t;
  }
  function crc32(bytes) {
    const t = crcTable();
    let c = -1;
    for (let i = 0; i < bytes.length; i++) c = t[(c ^ bytes[i]) & 0xFF] ^ (c >>> 8);
    return (c ^ -1) >>> 0;
  }
  function utf8(str) {
    if (typeof TextEncoder !== 'undefined') return new TextEncoder().encode(str);
    const out = [];
    for (let i = 0; i < str.length; i++) {
      let c = str.charCodeAt(i);
      if (c < 0x80) out.push(c);
      else if (c < 0x800) out.push(0xC0 | (c >> 6), 0x80 | (c & 63));
      else { out.push(0xE0 | (c >> 12), 0x80 | ((c >> 6) & 63), 0x80 | (c & 63)); }
    }
    return new Uint8Array(out);
  }

  /* ---- the smallest zip that is still a zip ---- */
  function zip(files) {
    const enc = files.map(f => ({ name: utf8(f.name), data: utf8(f.data) }));
    let total = 0, central = 0;
    enc.forEach(f => { total += 30 + f.name.length + f.data.length; central += 46 + f.name.length; });
    const out = new Uint8Array(total + central + 22);
    const dv = new DataView(out.buffer);
    let off = 0;
    const dirs = [];
    enc.forEach(f => {
      const crc = crc32(f.data), start = off;
      dv.setUint32(off, 0x04034b50, true);          // local header
      dv.setUint16(off + 4, 20, true);              // version needed
      dv.setUint16(off + 6, 0x0800, true);          // UTF-8 names
      dv.setUint16(off + 8, 0, true);               // stored
      dv.setUint16(off + 10, 0, true);              // time
      dv.setUint16(off + 12, 0x21, true);           // date: 1st Jan 1996, fixed so files compare equal
      dv.setUint32(off + 14, crc, true);
      dv.setUint32(off + 18, f.data.length, true);
      dv.setUint32(off + 22, f.data.length, true);
      dv.setUint16(off + 26, f.name.length, true);
      dv.setUint16(off + 28, 0, true);
      off += 30;
      out.set(f.name, off); off += f.name.length;
      out.set(f.data, off); off += f.data.length;
      dirs.push({ name: f.name, crc: crc, size: f.data.length, start: start });
    });
    const cdStart = off;
    dirs.forEach(d => {
      dv.setUint32(off, 0x02014b50, true);          // central directory header
      dv.setUint16(off + 4, 20, true);
      dv.setUint16(off + 6, 20, true);
      dv.setUint16(off + 8, 0x0800, true);
      dv.setUint16(off + 10, 0, true);
      dv.setUint16(off + 12, 0, true);
      dv.setUint16(off + 14, 0x21, true);
      dv.setUint32(off + 16, d.crc, true);
      dv.setUint32(off + 20, d.size, true);
      dv.setUint32(off + 24, d.size, true);
      dv.setUint16(off + 28, d.name.length, true);
      dv.setUint16(off + 30, 0, true);
      dv.setUint16(off + 32, 0, true);
      dv.setUint16(off + 34, 0, true);
      dv.setUint16(off + 36, 0, true);
      dv.setUint32(off + 38, 0, true);
      dv.setUint32(off + 42, d.start, true);
      off += 46;
      out.set(d.name, off); off += d.name.length;
    });
    dv.setUint32(off, 0x06054b50, true);            // end of central directory
    dv.setUint16(off + 8, dirs.length, true);
    dv.setUint16(off + 10, dirs.length, true);
    dv.setUint32(off + 12, off - cdStart, true);
    dv.setUint32(off + 16, cdStart, true);
    return out;
  }

  function xmlEsc(s) {
    return String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;')
      .replace(/"/g, '&quot;').replace(/\u0000-\u0008\u000b\u000c\u000e-\u001f/g, '');
  }
  function colName(i) {
    let s = '';
    i++;
    while (i > 0) { const m = (i - 1) % 26; s = String.fromCharCode(65 + m) + s; i = Math.floor((i - 1) / 26); }
    return s;
  }
  function isNum(v) {
    return typeof v === 'number' && isFinite(v);
  }

  /* sheets: [{ name, rows: [[cell, …], …], widths: [n, …] }]
     The first row of each sheet is treated as its heading. */
  PB.xlsxBytes = function (sheets) {
    const sheetXml = sheets.map(sh => {
      const rows = sh.rows || [];
      const cols = (sh.widths && sh.widths.length)
        ? '<cols>' + sh.widths.map((w, i) =>
            '<col min="' + (i + 1) + '" max="' + (i + 1) + '" width="' + w + '" customWidth="1"/>').join('') + '</cols>'
        : '';
      const body = rows.map((r, ri) =>
        '<row r="' + (ri + 1) + '">' + r.map((v, ci) => {
          const ref = colName(ci) + (ri + 1);
          if (v == null || v === '') return '';
          if (isNum(v)) return '<c r="' + ref + '"' + (ri === 0 ? ' s="1"' : '') + '><v>' + v + '</v></c>';
          return '<c r="' + ref + '" t="inlineStr"' + (ri === 0 ? ' s="1"' : '') +
                 '><is><t xml:space="preserve">' + xmlEsc(v) + '</t></is></c>';
        }).join('') + '</row>').join('');
      return '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>' +
        '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">' +
        '<sheetViews><sheetView workbookViewId="0">' +
        (rows.length > 1 ? '<pane ySplit="1" topLeftCell="A2" activePane="bottomLeft" state="frozen"/>' : '') +
        '</sheetView></sheetViews>' + cols +
        '<sheetData>' + body + '</sheetData></worksheet>';
    });

    const files = [
      { name: '[Content_Types].xml',
        data: '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>' +
          '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">' +
          '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>' +
          '<Default Extension="xml" ContentType="application/xml"/>' +
          '<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>' +
          '<Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>' +
          sheets.map((s, i) => '<Override PartName="/xl/worksheets/sheet' + (i + 1) +
            '.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>').join('') +
          '</Types>' },
      { name: '_rels/.rels',
        data: '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>' +
          '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">' +
          '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/>' +
          '</Relationships>' },
      { name: 'xl/workbook.xml',
        data: '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>' +
          '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" ' +
          'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"><sheets>' +
          sheets.map((s, i) => '<sheet name="' + xmlEsc((s.name || ('Sheet' + (i + 1))).slice(0, 31)) +
            '" sheetId="' + (i + 1) + '" r:id="rId' + (i + 1) + '"/>').join('') +
          '</sheets></workbook>' },
      { name: 'xl/_rels/workbook.xml.rels',
        data: '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>' +
          '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">' +
          sheets.map((s, i) => '<Relationship Id="rId' + (i + 1) +
            '" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet' +
            (i + 1) + '.xml"/>').join('') +
          '<Relationship Id="rId' + (sheets.length + 1) +
          '" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" Target="styles.xml"/>' +
          '</Relationships>' },
      { name: 'xl/styles.xml',
        data: '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>' +
          '<styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">' +
          '<fonts count="2"><font><sz val="11"/><name val="Calibri"/></font>' +
          '<font><b/><sz val="11"/><name val="Calibri"/></font></fonts>' +
          '<fills count="2"><fill><patternFill patternType="none"/></fill>' +
          '<fill><patternFill patternType="gray125"/></fill></fills>' +
          '<borders count="1"><border/></borders>' +
          '<cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs>' +
          '<cellXfs count="2"><xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/>' +
          '<xf numFmtId="0" fontId="1" fillId="0" borderId="0" xfId="0" applyFont="1"/></cellXfs>' +
          '<cellStyles count="1"><cellStyle name="Normal" xfId="0" builtinId="0"/></cellStyles>' +
          '</styleSheet>' }
    ].concat(sheets.map((s, i) => ({ name: 'xl/worksheets/sheet' + (i + 1) + '.xml', data: sheetXml[i] })));

    return zip(files);
  };
  PB.xlsxBlob = function (sheets) {
    return new Blob([PB.xlsxBytes(sheets)],
      { type: 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet' });
  };
  if (typeof module !== 'undefined' && module.exports) module.exports = PB;
})(typeof window !== 'undefined' ? (window.PB = window.PB || {}) : (globalThis.PB = globalThis.PB || {}));

