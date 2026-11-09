# PlasmidBench

A plasmid editor and cloning workbench that runs in the browser. The built app is a
single HTML file, so there is no server and nothing to install.

It has a circular and linear map viewer, sequence, feature, enzyme and ORF views, primer
design, a gel simulator, Sanger read verification, alignments, cloning simulators,
mutagenesis, a protein expression strategy tool, and a search of a lab's primer stock for
tubes that can sequence a given region. It reads GenBank, EMBL, FASTA, SnapGene .dna, .ab1,
GFF3, JSON and NCBI XML, and exports spreadsheets as .xlsx.

## Building

    node build.js

This joins the files in `src/` into `dist/PlasmidBench.html`, which you open in a browser.
It uses only Node's standard library. It was run with Node 22.

## Layout

    build.js                    the build script
    src/index.template.html     page markup, with {{path}} markers where files are inlined
    src/styles.css              all CSS
    src/js/                     the scripts, numbered in the order they load

The scripts share two globals, `Bio` (the engine, files 01 to 07) and `PB` (the app,
files 08 to 16), so the load order matters. Renumber with care.

    01-engine.js                  sequence, enzyme, primer and translation functions
    02-cloning.js                 cloning operations and assembly
    03-primer-stock-engine.js     choosing stocked primers to sequence a region
    04-expression-engine.js       protein descriptors and expression strategy
    05-formats-extra.js           GFF3, JSON and NCBI XML readers
    06-alignment-figures.js       alignment figures in the ESPript style
    07-read-interpretation.js     mapping a Sanger read onto codons and features
    08-app-core.js                state, storage, library and enzyme panels, inspector
    09-xlsx-writer.js             .xlsx writer
    10-app-views.js               map, sequence, feature, enzyme and ORF views
    11-app-map.js                 map viewer
    12-app-primers-gel-verify.js  primers, gel, verify, align, history
    13-app-tools.js               cloning simulators, mutagenesis, design, import and export
    14-primer-stock-ui.js         primer stock panel
    15-expression-ui.js           expression strategy panel
    16-app-main.js                wiring and start-up

## Where this came from

These files were produced by splitting an existing built copy of the app back into its
parts, so there is no history before the first commit. The tests are not in this
repository yet.


## Contributors

* Joseph Francis (Creator & Maintainer)
* Dr. Yufeng Tong (Project Advisor / Evaluation)
