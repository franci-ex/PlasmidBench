# Fixtures

Small files the suites read. Real ones, not synthetic, because the bugs are
in real files.

| file | what it is for |
|---|---|
| `pGro7.gb` | a chaperone vector whose CmR gene is reverse-strand, spliced **and** crosses the origin — all three at once, which is the case the read-interpretation module's own documentation names as the hard one |
| `pRSFDuet-1.dna` | a SnapGene binary file, for the `.dna` reader |
| `Primer_stock.tsv` | an oligo inventory in the shape a lab's own spreadsheet comes in |
| `Plasmid_stock.tsv` | the same for a plasmid collection |
| `Competent_cells.tsv` | the same for a strain list |

The larger corpus — real and deliberately malformed GenBank, EMBL, SnapGene
and ABI files — is **not** here. It is Biopython's `Tests/` directory, fetched
at validation time; see `tests/README.md`. Keeping it out makes its
provenance obvious and keeps this repository small.
