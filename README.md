# ![nfcore/test-datasets](docs/images/test-datasets_logo.png)
Test data to be used for automated testing with the nf-core pipelines

> ⚠️ **Do not merge your test data to `master`! Each pipeline has a dedicated branch (and a special one for modules)**

## Introduction

nf-core is a collection of high quality Nextflow pipelines. This repository contains various files for CI and unit testing of nf-core pipelines and infrastructure.

The principle for nf-core test data is as small as possible, as large as necessary. Please see the [guidelines](https://nf-co.re/docs/contributing/test_data_guidelines) for more detailed information. Always ask for guidance on the [nf-core slack](https://nf-co.re/join) before adding new test data.

## Documentation

nf-core/test-datasets comes with documentation in the `docs/` directory:

01. [Add a new  test dataset](https://github.com/nf-core/test-datasets/blob/master/docs/ADD_NEW_DATA.md)
02. [Use an existing test dataset](https://github.com/nf-core/test-datasets/blob/master/docs/USE_EXISTING_DATA.md)

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

## createtaxdb CI test specific information

### FASTA files

FASTA reference files used for building databases are copies of the nf-core/modules test dataset files (`sarscov2` and `haemophilus_influenzae` files) as of December 2023. 

- [sarscov2.fasta](https://github.com/nf-core/test-datasets/blob/0d5006780e17a3b11a36437d220c372c2e6e4ed0/data/genomics/sarscov2/genome/genome.fasta)
- [sarscov2.faa](https://github.com/nf-core/test-datasets/blob/89f6476aa0006451c1e9ea789ce4e4173c892319/data/genomics/sarscov2/genome/proteome.fasta)
- [haemophilus_influenzae.fna.gz](https://github.com/nf-core/test-datasets/blob/575e27aa850e186d4bcf85afc5572648aa35f2f4/data/genomics/prokaryotes/haemophilus_influenzae/genome/genome.fna.gz)


### taxonomy files

These are NCBI taxdump re-constructed files, where the entries only include those of the two FASTA files above (rather than the entire tax dump).

- Prot taxdump: as of December 2023