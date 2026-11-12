# test-datasets: tcrpmhcbinding

Test data to be used for automated testing with the nf-core pipelines.

This branch contains test data for the nf-core/tcrpmhcbinding pipeline.

## Content

`testdata-default` contains the test data to run the pipeline for a set of predictors.

- `samplesheet_test.csv` samplesheet for test data
- `samplesheet_test_airr.csv` samplesheet for `tcrpmhc_sample_airr.tsv`
- `tcrpmhc_sample_airr.tsv` the rows of `tcrpmhc_sample.tsv` with [AIRR](https://docs.airr-community.org/en/stable/datarep/rearrangements.html) column names (`species`, `mhc_allele_1`, `peptide_sequence_aa`, `v_call_tra`, `j_call_tra`, `junction_aa_tra`, `v_call_trb`, `j_call_trb`, `junction_aa_trb`) and canonical values (`mhc_allele_1` as `HLA-A*02:01`, `mhc_class` as `MHC-I`). Used by nf-core/tcrpmhcbinding from the next release on.
- `sample_gil.tsv` test dataset derived from this [publication](https://www.biorxiv.org/content/10.64898/2026.02.27.708137v2.abstract). 


