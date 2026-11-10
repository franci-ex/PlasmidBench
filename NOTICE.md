# Where the data came from

The code in this repository is MIT (see `LICENSE`). Several of the tables it
ships were compiled from other people's work, and those have their own terms.
This file says which, so that anyone redistributing the built file — or
publishing with it — knows what they are redistributing. Section 9 of the
project evaluation asks for exactly this list.

## Restriction enzymes

**607 enzymes with their recognition sequences, cut positions and supplier
codes, and 15 supplier names.**

From REBASE, the Restriction Enzyme Database (Roberts, R.J. et al., *Nucleic
Acids Research*), by way of the copy Biopython distributes in
`Bio.Restriction.Restriction_Dictionary`.

REBASE asks to be cited and is free for academic use. Its own terms should be
checked before the table is redistributed in a commercial product or
published as a dataset in its own right: academic use is straightforward,
redistribution is the part worth reading.

- Cite: Roberts RJ, Vincze T, Posfai J, Macelis D. *REBASE — a database for
  DNA restriction and modification: enzymes, genes and genomes.* Nucleic Acids
  Research.
- <http://rebase.neb.com>

The validation suite compares every one of those enzymes against Biopython's
copy, so the table is checkable rather than asserted: see
`tests/oracles/test-rebase.py`.

## Genetic codes

**The NCBI translation tables**, from the NCBI Taxonomy genetic code tables.
Public domain as US Government work; cite NCBI.

- <https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi>

## Codon usage

**Codon frequencies per host**, derived from the Kazusa Codon Usage Database
(Nakamura, Y. et al., *Nucleic Acids Research* 2000), which in turn derives
from GenBank. Cite Kazusa and NCBI.

- <https://www.kazusa.or.jp/codon/>

## Amino-acid constants

**Residue masses, hydropathy and the other per-residue values** used by the
protein descriptors come from the published sources named in the engine
beside each table: Kyte & Doolittle 1982 (hydropathy), Lobry & Gautier 1994
(aromaticity), Guruprasad et al. 1990 (instability), Ikai 1980 (aliphatic
index), Bjellqvist et al. 1993/1994 (isoelectric point). They are published
constants, cited where used, and checked against `Bio.SeqUtils.ProtParam`.

## Melting temperature

**Nearest-neighbour thermodynamics** from SantaLucia 1998 (PNAS), with the
salt corrections of Owczarzy et al. 2004 and 2008. Published constants, cited
in `src/engine.js` beside the table.

## Oligonucleotide weight

**The per-residue masses** are IDT's published values for a synthetic
oligonucleotide with a free 5'-OH. That convention is the one a lab needs,
and `tests/oracles/cmpmw.py` checks the implementation against both IDT's
formula and `Bio.SeqUtils.molecular_weight`, which uses the other convention.

## Substitution matrix

**BLOSUM62** (Henikoff & Henikoff 1992). Checked cell by cell against
`Bio.Align.substitution_matrices`.

## Gel ladders

**Band sizes** for the NEB 1 kb, 1 kb Plus and 100 bp ladders, the Thermo
GeneRuler 1 kb and 1 kb Plus, and lambda DNA cut with HindIII, from the
manufacturers' published band lists. Product names are the manufacturers'
trademarks and are used to identify the products.

## Feature library

**The 27 built-in feature sequences** are short, widely published regulatory
elements (T7 and lac promoters, common terminators and tags) written from
their published sequences.

The project does **not** ship a large feature library, deliberately — section
10 of the evaluation argues that curated open libraries already exist and
should be imported rather than rebuilt. If you import one, its licence
travels with it and is recorded per entry:

- **ApE**'s default feature file — ApE's own terms.
- **FPbase** — open API, CC-BY for most content; check per protein.
- **pLannotate** is GPL-3.0, and its primary feature set traces through
  GenoLIB to SnapGene's curation. Using its *approach* carries no obligation;
  incorporating its code would impose copyleft, and bundling its data
  warrants reading the chain of terms.
- **GenoLIB** (Adames et al., *NAR* 2015) — the paper is CC BY 4.0.

An imported library is stored in the browser, not in this repository, and
every entry keeps the source and licence it was imported under.

## Test corpus

The validation suite compares against **Biopython 1.88** (Cock et al. 2009,
*Bioinformatics*; BSD-3-Clause / Biopython License) and reads real and
deliberately malformed files from **Biopython's own `Tests/` directory**,
which is fetched at validation time and is not redistributed here. See
`tests/README.md` for the one command that fetches it.

## AI assistance

Parts of this codebase were written with the help of Anthropic's Claude,
including engine modules, the validation suites and this documentation. The
scientific choices, the validation strategy, and the decisions about what the
tool should and should not claim are the author's. Section 5.4 of the
evaluation asks for this to be recorded, and it is recorded here and in
`AGENTS.md`.

Every algorithm in the engine is explained in a comment beside it, and every
numerical claim has an outside oracle in `tests/oracles/` — which is the only
thing that makes AI-assisted scientific code reviewable.
