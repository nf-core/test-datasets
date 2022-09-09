# ![nfcore/test-datasets](docs/images/test-datasets_logo.png)

Test data to be used for automated testing with the nf-core pipelines.

This branch contains test data for the [nf-core/crisprseq](https://github.com/nf-core/crisprseq) pipeline.

## Content

### Analysis of CRISPR editing

`testdata-edition` contains the a minimal test dataset needed to run the pipeline for the analysis of CRISPR editing.

- `samplesheet_test.csv` contains the input samplesheet file to run the pipeline.
- The `fastq`files are samples from a simulated run (`chr6`) and from a MiSeq amplicon sequencing run of the target sites `TRAC` and `AAVS1`.
