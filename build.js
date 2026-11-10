const fs = require('fs');
const path = require('path');
const dir = __dirname;
const read = f => fs.readFileSync(path.join(dir, f), 'utf8');

const parts = [
  read('src/part1-head.html'),
  read('src/part2-body.html'),
  '<script>\n' + read('src/engine.js') + '\n</script>',
  '<script>\n' + read('src/engine2.js') + '\n</script>',
  '<script>\n' + read('src/engine3.js') + '\n</script>',
  '<script>\n' + read('src/engine4.js') + '\n</script>',
  '<script>\n' + read('src/engine5.js') + '\n</script>',
  '<script>\n' + read('src/engine6.js') + '\n</script>',
  '<script>\n' + read('src/engine7.js') + '\n</script>',
  '<script>\n' + read('src/engine8.js') + '\n</script>',
  '<script>\n' + read('src/store.js') + '\n</script>',
  '<script>\n' + read('src/app-core.js') + '\n</script>',
  '<script>\n' + read('src/xlsx.js') + '\n</script>',
  '<script>\n' + read('src/app-views1.js') + '\n</script>',
  '<script>\n' + read('src/app-map.js') + '\n</script>',
  '<script>\n' + read('src/app-views2.js') + '\n</script>',
  '<script>\n' + read('src/app-tools.js') + '\n</script>',
  '<script>\n' + read('src/elab.js') + '\n</script>',
  '<script>\n' + read('src/profile.js') + '\n</script>',
  '<script>\n' + read('src/app-profile.js') + '\n</script>',
  '<script>\n' + read('src/app-annotate.js') + '\n</script>',
  '<script>\n' + read('src/app-elab.js') + '\n</script>',
  '<script>\n' + read('src/app-seq.js') + '\n</script>',
  '<script>\n' + read('src/app-express.js') + '\n</script>',
  '<script>\n' + read('src/app-init.js') + '\n</script>'
];

/* Every engine file after the first uses a CommonJS require() guard so the
   Node test suites can load it on its own. In the browser there is one
   window.Bio and no module loader, so the guards are stripped here. Done by
   pattern over every part rather than by index, so adding an engine file
   does not silently leave a require() in the shipped page — the check at the
   foot of this file would then warn. */
for (let i = 0; i < parts.length; i++) {
  parts[i] = parts[i]
    .replace(/const E = \(typeof require !== 'undefined' && typeof module !== 'undefined'\)\s*\?\s*require\('\.\/engine\.js'\)\s*:\s*root\.Bio;/g, 'const E = root.Bio;')
    .replace(/[ \t]*if \(typeof require !== 'undefined' && typeof module !== 'undefined'\) \{[^}]*\}\n?/g, '')
    .replace(/[ \t]*if \(typeof require !== 'undefined' && typeof module !== 'undefined'\) require\('\.\/engine\d?\.js'\);\n?/g, '');
}

const out = parts.join('\n');
const dest = path.join(dir, 'plasmidbench.html');
fs.writeFileSync(dest, out);
const kb = (Buffer.byteLength(out) / 1024).toFixed(0);
console.log('wrote plasmidbench.html — ' + kb + ' KB');
if (/require\(/.test(out)) {
  const lines = out.split('\n');
  lines.forEach((l, i) => { if (/require\(/.test(l)) console.log('  WARNING require() at line ' + (i + 1) + ': ' + l.trim().slice(0, 90)); });
}
if (/module\.exports/.test(out)) {
  console.log('  note: module.exports present (guarded by typeof check, harmless in browser)');
}
