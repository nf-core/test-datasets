# test-datasets: `plasmodiumdrugres`

This branch contains test data to be used for automated testing with the [nf-core/plasmodiumdrugres](https://github.com/nf-core/plasmodiumdrugres) pipeline.

## Content of this repository


The pipeline analyzes *Plasmodium* drug resistance markers from allele tables or [Portable Microhaplotype Object (PMO)](https://plasmogenepi.github.io/PMO_Docs/) files, translating loci of interest to amino acid changes and estimating single-locus and multi-locus allele frequencies and prevalences.

All files live under `testdata/`.

### Pipeline inputs

- `testdata/example_PMO.json`: Example Portable Microhaplotype Object (PMO) file for the PMO entry point.
- `testdata/allele_table.tsv`: Microhaplotype allele table (specimen, target, allele, reads, sequence) for the allele-table entry point extracted from `example_PMO.json`.
- `testdata/panel_info.bed`: Panel information BED including insert reference sequences extracted from `example_PMO.json`.
- `testdata/panel_info_no_ref.bed`: Panel information BED without reference sequences extracted from `example_PMO.json`.
- `testdata/dummy_panel_info_fake_chroms.bed`: Panel information BED with synthetic chromosome names to test adding ref from whole genome using `insert_refseqs.fasta`.
- `testdata/insert_refseqs.fasta`: Targeted insert/reference sequences matching panel target names.
- `testdata/loci_of_interest.bed`: Drug-resistance loci of interest (amino acid positions) used for translation.
- `testdata/loci_groups.tsv`: Multi-locus groups (e.g. `pfdhfr_pfdhps`, `crt`) for multi-locus frequency estimation.
- `testdata/population_map.tsv`: Specimen-to-population assignment table.

### Intermediate and module test inputs

- `testdata/amino_acid_calls.tsv`: Amino acid call table produced from translation of loci of interest.
- `testdata/loci_of_interest_mhaps.bed`: Loci-of-interest annotations linked to microhaplotype targets (for microhaplotype-based SLAF tests).
- `testdata/mhaps_slaf.tsv`: Microhaplotype single-locus allele frequencies.
- `testdata/aa_mlaf.tsv`: Amino acid multi-locus allele frequencies (MLAF).
- `testdata/mlaf_pop1.tsv`: Multi-locus allele frequencies for population `pop1`.
- `testdata/mlaf_pop2.tsv`: Multi-locus allele frequencies for population `pop2`.
- `testdata/population_map_indexed.tsv`: Specimen-to-population assignment including population index IDs.
- `testdata/population_map_with_spaces.tsv`: Population assignment table with spaces in population names.
- `testdata/population_index_lookup.tsv`: Lookup table mapping population index IDs to population names.
- `testdata/empty_population_index_lookup.tsv`: Empty population index lookup file for edge-case tests.
- `testdata/allele_prev.tsv`: Expected allele prevalence estimates.
- `testdata/aa_slaf.tsv`: Expected amino acid single-locus allele frequencies (SLAF).
- `testdata/amino_acid_calls.aa_sl_from_ml.tsv`: Expected single-locus amino acid frequencies derived from multi-locus frequencies.
- `testdata/sl_from_ml.tsv`: Expected single-locus frequencies derived from multi-locus allele frequencies.
- `testdata/sl_from_mlaf_pop1.tsv`: Expected single-locus frequencies derived from MLAF for population `pop1`.
- `testdata/sl_from_mlaf_pop2.tsv`: Expected single-locus frequencies derived from MLAF for population `pop2`.
- `testdata/sl_pop1.tsv`: Expected single-locus allele frequencies and prevalences for population `pop1`.
- `testdata/sl_pop2.tsv`: Expected single-locus allele frequencies and prevalences for population `pop2`.
