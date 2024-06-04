# test-datasets: `reportho`
This branch contains data used for automated testing of [nf-core/reportho](https://github.com/nf-core/reportho).

## Content of this repository
- `databases/` – minified databases for testing local search
  - `1_members-mini.tsv.gz` – minified [EggNOG](http://eggnog5.embl.de/#/app/home); contains all identifiers necessary for BicD2, plus 20 random lines; gzipped
  - `AllOrthologs.txt` – minified [PANTHER](https://www.pantherdb.org/); contains all identifiers necessary for BicD2, plus 50 random lines
  - `latest.Eukaryota-mini.tsv.gz` –  minified Uniprot-EggNOG ID map; contains all identifiers necessary for BicD2, plus 20 random lines; gzipped
  - `oma-ensembl-mini.txt.gz` – minified OMA-Ensembl ID map; contains 100 random lines; gzipped
  - `oma-mini.txt.gz` – minified OMA groups database; contains all identifiers necessary for BicD2, plus 10 random lines; gzipped
  - `oma-refseq-mini.txt.gz` – minified OMA-Refseq ID map; contains 100 random lines; gzipped
  - `oma-uniprot-mini.txt.gz` – minified OMA-Uniprot ID map; contains all data necessary for BicD2, plus 20 random lines; gzipped
- `samplesheet/` – various test samplesheets
  - `samplesheet.csv` – samplesheet for running with Uniprot IDs as query
  - `samplesheet_fasta.csv` – samplesheet for running with FASTA files as query
- `sequences/` – sequences used in the FASTA test
  - `ste2.fa`, `ste3.fa` – FASTA sequences
 
## ID and sequence information

### Uniprot samplesheet

| Sample ID | Description                   | Uniprot                                                   |
|-----------|-------------------------------|-----------------------------------------------------------|
| BicD2     | Human bicaudal protein D2     | [Uniprot](https://www.uniprot.org/uniprotkb/Q8TD16/entry) |
| HBB       | Human hemoglobin subunit beta | [Uniprot](https://www.uniprot.org/uniprotkb/P68871/entry) |

### FASTA samplesheet

Sequences are sourced from the Uniprot entry.

| Sample ID | Description                   | Uniprot                                                   |
|-----------|-------------------------------|-----------------------------------------------------------|
| STE2      | Yeast pheromone receptor STE2 | [Uniprot](https://www.uniprot.org/uniprotkb/D6VTK4/entry) |
| STE3      | Yeast pheromone receptor STE3 | [Uniprot](https://www.uniprot.org/uniprotkb/P06783/entry) |
