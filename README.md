# ![nfcore/test-datasets](docs/images/test-datasets_logo.png)
Test data to be used for automated testing with the nf-core pipelines

## Introduction

This branch contains test data for [nf-core/detaxizer](https://nf-co.re/detaxizer/) pipeline.

## Full-size test data 
The `samplesheet.full.csv` links to gut metagenome data of antibiotic-treated patients originating from [Bertrand et al. *Nature Biotechnology* (2019)](https://doi.org/10.1038/s41587-019-0191-2). This dataset is used by [nf-core/mag](https://nf-co.re/mag/) as a full test and is now also used for nf-core/detaxizer.

| SAMPLE    | ILLUMINA READS: ENA ID | ONT READs: ENA ID |
|-----------|------------------------|-------------------|
| CAPES S11 | ERR3201918             | ERR3201942        |
| CAPES S21 | ERR3201928             | ERR3201952        |
| CAPES S7  | ERR3201914             | ERR3201938        |

## Test data
For a functionality test of nf-core/detaxizer the small test data sets of nf-core/mag (paired-end short reads: `test_minigut_*`) and [nf-core/bacass](https://nf-co.re/bacass) (long reads: `subset350.fq.gz`) are used together with the `minigut_kraken.tgz` of nf-core/mag as `kraken2` database and the host reference (`genome.hg38.chr21_10000bp_region.fa`) also from nf-core/mag test data set.

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

## Support

For further information or help, don't hesitate to get in touch on our [Slack organisation](https://nf-co.re/join/slack) (a tool for instant messaging).

[^1]: From [stackoverflow](https://stackoverflow.com/a/60846265/11502856)
