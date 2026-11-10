# PlasmidBench

A browser-based plasmid editor and cloning workbench. The app is a single HTML file:
no server, no installation.

**Status:** early development (v0.0.1). Features and behaviour may change between releases.

## Features

- Circular and linear maps; sequence, feature, enzyme and ORF views
- Primer design, gel simulation, Sanger read verification and alignments
- Cloning simulation, mutagenesis and protein expression strategy
- Primer stock search for sequencing a region
- Reads GenBank, EMBL, FASTA, SnapGene `.dna`, `.ab1`, GFF3, JSON and NCBI XML; exports `.xlsx`

## Usage

Download `PlasmidBench.html` from the [latest release](https://github.com/franci-ex/PlasmidBench/releases)
(or `dist/` in this repository) and open it in a web browser.

## Building from source

Requires [Node.js](https://nodejs.org/) (tested with v22); no other dependencies.

    node build.js

This combines the files in `src/` into `dist/PlasmidBench.html`.

## Project layout

    build.js                  build script
    src/index.template.html   page markup with {{path}} include markers
    src/styles.css            styles
    src/js/                   scripts, numbered in load order

Files `01`–`07` define the engine (`Bio`); files `08`–`16` define the app (`PB`).
Load order matters, so renumber with care.

## Contributors

- Joseph Francis — creator and maintainer
- Dr. Yufeng Tong — project advisor
