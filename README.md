# ![nfcore/test-datasets](docs/images/test-datasets_logo.png)

Test data to be used for automated testing with the nf-core pipelines

## Test data for nf-core/mag

This branch contains test data for the [nf-core/mag](https://github.com/nf-core/mag) pipeline.

## Full-size test data

The `manifext.full.tsv` links to gut metagenome data of antibiotic-treated patients originating from [Bertrand et al. _Nature Biotechnology_ (2019)](https://doi.org/10.1038/s41587-019-0191-2).

| SAMPLE    | ILLUMINA READS: ENA ID | ONT READs: ENA ID |
| --------- | ---------------------- | ----------------- |
| CAPES S11 | ERR3201918             | ERR3201942        |
| CAPES S21 | ERR3201928             | ERR3201952        |
| CAPES S7  | ERR3201914             | ERR3201938        |

## Downloading test data

Due the large number of large files in this repository for each pipeline, we highly recommend cloning only the branches you would use.

```bash
git clone <url> --single-branch --branch <pipeline/modules/branch_name>
```

To subsequently clone other branches[^1]

```bash
git remote set-branches --add origin [remote-branch]
git fetch
```

## Broken samplesheets

For testing input validation, the `samplesheets` directory contains the `broken/` subdirectory containing samplesheets with errors that should be caught by the pipeline.

-`samplesheets/broken/assembly_samplesheet_invalid_assembler.csv`: has a mistyped and missing (mandatory) assembler entry -`samplesheets/broken/assembly_samplesheet_invalid_id.csv`: incorrectly has a space in the sample ID name -`samplesheets/broken/assembly_samplesheet_missing_group.csv`: missing the mandatory group column -`samplesheets/broken/assembly_samplesheet_missing_id.csv`: missing the mandatory id column -`samplesheets/broken/assembly_samplesheet_nonunique_fasta.csv`: has duplicate FASTA files in two rows -`samplesheets/broken/samplesheet_empty_group.csv`: group column has a header but empty entries (must have a group) -`samplesheets/broken/samplesheet_empty_run.csv`: group column has a header but empty entries (if column is present, can't be empty) -`samplesheets/broken/samplesheet_lr_without_sr2.csv`: long-reads given, but missing read 2 required for hybrid assembly -`samplesheets/broken/samplesheet_missing_r1.csv`: missing mandatory read 1 column -`samplesheets/broken/samplesheet_missing_sample.csv`: missing mandatory sample column -`samplesheets/broken/samplesheet_nonunique_sample_run_combination.csv`: has invalid duplicate sample-run combinations -`samplesheets/broken/samplesheet_spaces_in_name.csv`: incorrect sample name with spaces

## Support

For further information or help, don't hesitate to get in touch on our [Slack organisation](https://nf-co.re/join/slack) (a tool for instant messaging).

[^1]: From [stackoverflow](https://stackoverflow.com/a/60846265/11502856)
