# test-datasets: `eupathpopgen`

This branch contains test data for automated testing with the [nf-core/eupathpopgen](https://github.com/nf-core/eupathpopgen) pipeline (currently developed at [PlasmoGenEpi/eupathpopgen](https://github.com/PlasmoGenEpi/eupathpopgen)).

## Content of this repository

All files live under `testdata/`. These are the same simulated microhaplotype inputs used by [nf-core/plasmodiumdrugres](https://github.com/nf-core/plasmodiumdrugres) (shared from the `plasmodiumdrugres` test-datasets branch) so related PlasmoGenEpi pipelines stay consistent.

### Pipeline inputs

- `testdata/allele_table.tsv`: Microhaplotype allele table (`specimen_name`, `target_name`, `reads`, `seq`, plus optional `bioinformatics_run_name` / `allele` columns from the drugres extract).
- `testdata/population_assignment.tsv`: Specimen-to-population map (`specimen_name`, `population`) for optional `--population_map` runs.

## Usage

```bash
nextflow run nf-core/eupathpopgen \
  -profile test,docker \
  --outdir results
```

Raw URLs (branch `eupathpopgen`):

- `https://raw.githubusercontent.com/nf-core/test-datasets/eupathpopgen/testdata/allele_table.tsv`
- `https://raw.githubusercontent.com/nf-core/test-datasets/eupathpopgen/testdata/population_assignment.tsv`
