# ![nfcore/test-datasets](docs/images/test-datasets_logo.png)
Test data to be used for automated testing with the nf-core pipelines

## Introduction

nf-core is a collection of high quality Nextflow pipelines. This repository contains various files for CI and unit testing of nf-core pipelines and infrastructure.

The principle for nf-core test data is as small as possible, as large as necessary. Please see the [guidelines](https://nf-co.re/docs/contributing/test_data_guidelines) for more detailed information. Always ask for guidance on the [nf-core slack](https://nf-co.re/join) before adding new test data.

## Documentation

nf-core/test-datasets comes with documentation in the `docs/` directory:

01. [Add a new  test dataset](https://github.com/nf-core/test-datasets/blob/master/docs/ADD_NEW_DATA.md)
02. [Use an existing test dataset](https://github.com/nf-core/test-datasets/blob/master/docs/USE_EXISTING_DATA.md)

## Data generation

We use `SB26` as a test dataset from [_Keren-Shaul et al., Nat Prot. 2019_](https://tanaylab.github.io/old_resources/pages/672.html). We subset it only for Amplification batch `AB339` as shown below. The helper script is part of the [nf-core/marsseq](https://github.com/nf-core/marsseq) pipeline.

```bash
filter_reads.py --input ../SB26-orig --output SB26/ --batches AB339 && \
  mv SB26/Undetermined_S0_L001_R1_001.fastq.gz SB26/AB339_R1.fastq.gz && \
  mv SB26/Undetermined_S0_L001_R2_001.fastq.gz SB26/AB339_R2.fastq.gz
```

## Citations

Keren-Shaul, H., Kenigsberg, E., Jaitin, D.A. et al. MARS-seq2.0: an experimental and analytical pipeline for indexed sorting combined with single-cell RNA sequencing. Nat Protoc 14, 1841–1862 (2019). [https://doi.org/10.1038/s41596-019-0164-4](https://doi.org/10.1038/s41596-019-0164-4)

## Support

For further information or help, don't hesitate to get in touch on our [Slack organisation](https://nf-co.re/join/slack) (a tool for instant messaging).

[^1]: From [stackoverflow](https://stackoverflow.com/a/60846265/11502856)
