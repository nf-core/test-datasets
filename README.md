# test-datasets: `eupathpopgen`

This branch contains test data for automated testing with the [nf-core/eupathpopgen](https://github.com/nf-core/eupathpopgen) pipeline.

## Content of this repository

All files are located under `testdata/`. The allele table and PMO examples represent the same dataset in different formats.

The data were simulated using [recombuddy](https://github.com/PlasmoGenEpi/recombuddy) from a background population derived from publicly available whole-genome sequencing data from Southeast Asia (SEA). The simulated parasite genomes were then processed in silico using the MAD4HaTteR targeted sequencing panel to generate realistic targeted-sequencing FASTQ files. These FASTQ files were subsequently processed using the SeekDeep pipeline.

The PMO file also contains synthetic sample metadata.

### Pipeline inputs

- `testdata/allele_table.tsv`: Microhaplotype allele table (`specimen_name`, `target_name`, `reads`, `seq`, plus optional `bioinformatics_run_name` / `allele` columns from the drugres extract).
- `testdata/population_assignment.tsv`: Specimen-to-population map (`specimen_name`, `population`) for optional `--population_map` runs.
- `testdata/example_PMO.json`: Example Portable Microhaplotype Object for `--pmo` / `-profile test_pmo`.
- `testdata/population_assignment_unassigned.tsv`: Specimen-to-population assignment table missing some specimens (`Vietnam2018-23`, `Vietnam2018-24`) to test handling of unassigned specimens.

### Module test inputs

Intermediate files produced by running the pipeline on `allele_table.tsv` with `population_assignment.tsv` (naive COI). Used as inputs to the module nf-tests.

- `testdata/population_map_indexed.tsv`: Indexed population map (`specimen_name`, `population`, `population_index`) from `INDEX_POPULATION_ASSIGNMENT`.
- `testdata/population_index_lookup.tsv`: Population index to population label lookup from `INDEX_POPULATION_ASSIGNMENT`.
- `testdata/popidx_001.allele_table.tsv.gz`: Allele table for one population from `SPLIT_ALLELE_TABLE_BY_POP`.
- `testdata/popidx_001.coi_table.tsv`: Naive COI table from `ESTIMATE_COI_NAIVE`.
- `testdata/popidx_001.naive.coi_table.tsv`: Standardised COI table from `STANDARDIZE_COI_TABLE`.
- `testdata/popidx_001.naive.coi_table.labeled.tsv`, `testdata/popidx_002.naive.coi_table.labeled.tsv`: COI tables with a population column from `ADD_POPULATION_COLUMN`.
- `testdata/popidx_001.allele_summary_by_target.tsv`: Per-locus allele summary from `ALLELE_PER_LOCUS_SUMMARY`.
- `testdata/popidx_001.per_locus_popgen_summary.tsv`: Per-locus popgen summary from `PER_LOCUS_POPGEN_SUMMARY`.
