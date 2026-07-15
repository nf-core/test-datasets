# bioBakery test data

Test datasets for [bioBakery](https://huttenhower.sph.harvard.edu/biobakery-workflows/) tools.

## MetaPhlAn toy database

A minimal CHOCOPhlAn database for running MetaPhlAn without downloading the full 60 GB index.
Built from `mpa_vJan25_CHOCOPhlAnSGB_202503` (produced 2026-04-30) by randomly sampling 10 SGBs
per kingdom (Bacteria, Archaea, Eukaryota) and keeping up to 200 markers per SGB, then building
a Bowtie2 large index from the resulting FASTA.

| Metric              | Value                        |
| ------------------- | ---------------------------- |
| Kingdoms            | Bacteria, Archaea, Eukaryota |
| SGBs per kingdom    | 10 (30 total)                |
| Max markers per SGB | 200 (6,000 total)            |
| Archive size        | ~14 MB                       |

Files:

- `mpa_vJan25_TOY_CHOCOPhlAnSGB_202503.tar.gz` — archive containing the marker database pickle,
  marker sequences (FASTA), and Bowtie2 large index files.
- `test_genefamilies.tsv` – mock file for gene family abundance data, to test the HUMAnN3 module.
