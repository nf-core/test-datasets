# ![nfcore/test-datasets](docs/images/test-datasets_logo.png)

Test data to be used for automated testing with the nf-core pipelines.

This branch contains test data for the [nf-core/crisprseq](https://github.com/nf-core/crisprseq) pipeline.

## Content

### Analysis of CRISPR editing

`testdata-edition` contains the a minimal test dataset needed to run the pipeline for the analysis of CRISPR editing.

- `samplesheet_test.csv` contains the input samplesheet file to run the pipeline.
- The `fastq`files are samples from a simulated run (`chr6`) and from a MiSeq amplicon sequencing run of the target sites `TRAC` and `AAVS1`.

- `samplesheet_test_full.csv` contains the input samplesheet file to run the pipeline with full size data. Data is obtained from ENA project [PRJNA326019](https://www.ebi.ac.uk/ena/browser/view/PRJNA326019) obtained by [VanOverbeek et al. 2016](https://doi.org/10.1016/j.molcel.2016.06.037)

- `samplesheet_test_umis.csv` contains the input samplesheet file to run the pipeline with the option of UMIs. Samples were obtained from ENA project [ERR9897751](https://www.ebi.ac.uk/ena/browser/view/ERR9897751) and [ERR9897765](https://www.ebi.ac.uk/ena/browser/view/ERR9897765) and subsampled to with seqtk and a seed of 100.
  ```
    seqtk sample -s100 ERR9897751_1.fastq 632577 > ERR9897751_1_subsample.fastq
    seqtk sample -s100 ERR9897751_2.fastq 632577 > ERR9897751_2_subsample.fastq
    seqtk sample -s100 ERR9897765_1.fastq 606960 > ERR9897765_1_subsample.fastq
    seqtk sample -s100 ERR9897765_2.fastq 606960 > ERR9897765_2_subsample.fastq
  ```