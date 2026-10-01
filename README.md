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
