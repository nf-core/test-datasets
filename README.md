# ![nfcore/test-datasets](docs/images/test-datasets_logo.png)

Test data to be used for automated testing with the nf-core pipelines.

This branch contains test data for the [nf-core/crisprseq](https://github.com/nf-core/crisprseq) pipeline.

## Content

### Analysis of CRISPR editing

`testdata-edition` contains the a minimal test dataset needed to run the pipeline for the analysis of CRISPR editing.

- `samplesheet_test.csv` contains the input samplesheet file to run the pipeline.
- The `fastq`files are samples from a simulated run (`chr6`) and from a MiSeq amplicon sequencing run of the target sites `TRAC` and `AAVS1`.

- `samplesheet_test_full.csv` contains the input samplesheet file to run the pipeline with full size data. Data is obtained from ENA project [PRJNA326019](https://www.ebi.ac.uk/ena/browser/view/PRJNA326019) obtained by [VanOverbeek et al. 2016](https://doi.org/10.1016/j.molcel.2016.06.037)
