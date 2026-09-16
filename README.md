# test-datasets: `eupathpopgen`

This branch contains test data for automated testing with the [nf-core/eupathpopgen](https://github.com/nf-core/eupathpopgen) pipeline.

## Content of this repository

All files live under `testdata/`.

### Pipeline inputs

- `testdata/allele_table.tsv`: Microhaplotype allele table (`specimen_name`, `target_name`, `reads`, `seq`, plus optional `bioinformatics_run_name` / `allele` columns from the drugres extract).
- `testdata/population_assignment.tsv`: Specimen-to-population map (`specimen_name`, `population`) for optional `--population_map` runs.
- `testdata/example_PMO.json`: Example Portable Microhaplotype Object for `--pmo` / `-profile test_pmo`.
