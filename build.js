#!/usr/bin/env node
// Builds dist/PlasmidBench.html from the files in src/.
// src/index.template.html contains {{path}} markers (paths relative to src/);
// each one is replaced by the contents of that file, unchanged.
// Needs Node, nothing else. Run it from anywhere: node build.js
'use strict';
const fs = require('fs');
const path = require('path');

const SRC = path.join(__dirname, 'src');
const OUT = path.join(__dirname, 'dist', 'PlasmidBench.html');

const template = fs.readFileSync(path.join(SRC, 'index.template.html'), 'utf8');
const used = [];
const html = template.replace(/\{\{([\w./-]+)\}\}/g, (_, rel) => {
  used.push(rel);
  return fs.readFileSync(path.join(SRC, rel), 'utf8');
});

// a source file that the template never includes would silently go missing
const onDisk = ['styles.css'].concat(
  fs.readdirSync(path.join(SRC, 'js')).filter(f => f.endsWith('.js')).map(f => 'js/' + f));
const unused = onDisk.filter(f => used.indexOf(f) < 0);
if (unused.length) {
  console.error('not included by the template: ' + unused.join(', '));
  process.exit(1);
}

fs.mkdirSync(path.dirname(OUT), { recursive: true });
fs.writeFileSync(OUT, html);
console.log('wrote ' + path.relative(process.cwd(), OUT) + ' (' + used.length + ' parts)');
