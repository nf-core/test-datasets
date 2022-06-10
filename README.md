# ![nfcore/test-datasets](docs/images/test-datasets_logo.png)
Test data to be used for automated testing with the nf-core pipeline **hgtseq**

## Introduction

nf-core is a collection of high quality Nextflow pipelines. This repository contains various files for CI and unit testing of nf-core pipelines and infrastructure.

The principle for nf-core test data is as small as possible, as large as necessary. Always ask for guidance on the [nf-core slack](https://nf-co.re/join) before adding new test data.

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

## Content

### CSV files

The `testdata/csv` folder contains test input files, and in particular:

- `input_bams.csv` contains example data to test the pipeline with provided bam files (i.e. aligned data)
- `input_fastqs.csv` containes example data to test the pipeline with provided fastq files (i.e. raw data): these are paired end.

### FASTQ files

The folder `testdata/fastq` contains raw reads to be used as primary input to the pipeline and in particular 3 samples:

- `testsample01` is extracted from *SRR13106578* 
- `testsample02` is extracted from *SRR13106582*
- `testsample03` is extracted from *SRR17085829*

All reads have been first processed with the pipeline and only those classified with *kraken2* according to the pipeline criteria have been used to create the sample fastq files. All reads map to chromosome 21 of the Human Genome.

### BAM files

The folder `testdata/bam`containes aligned reads to be used as alternative input to the pipeline. The bam files have been obtained by aligning the sample fastq reads in `testadata/fastq` to chromosome 21 of the Human Genome:

- `testsample01.bam` corresponds to *SRR13106578* sample reads
- `testsample02.bam` corresponds to *SRR13106582* sample reads

### Classified files

The folder `testdata/classified` contains reads extracted from the bam files in `testdata/bam`, using samtools flag and are categorised in pairs where both reads are unmapped (in `testdata/classified/both_unmapped`) and in reads where only one of the pair is not mapped (in `testdata/classified/single_unmapped`).

Like for the above sample data, they represent reads extracted from *SRR13106578* (testsample01) and *SRR13106582* (testsample02). These files are example inputs for the downstream subworkflows of this pipeline.

### Reference

The folder `testdata/reference` only contains the dictionary file of the Human reference in GRCh38, used to plot circular plots in R with *ggbio*

## Support

For further information or help, don't hesitate to get in touch on our [Slack organisation](https://nf-co.re/join/slack) (a tool for instant messaging).

[^1]: From [stackoverflow](https://stackoverflow.com/a/60846265/11502856)
