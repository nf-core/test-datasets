# ![nfcore/test-datasets](docs/images/test-datasets_logo.png)
Test data to be used for automated testing with the nf-core pipelines

## nf-core/circrna
This branch contains test-data for `nf-core/circrna`.

### Contents of branch:
* `fastq/` 9 FASTQ read pairs.
* `reference/` Reference annotation files
* `phenotype.csv` metadata file for DESeq2.
* `samples.csv` input test-dataset csv file

### Test-dataset generation strategy:
Gencode GRCh38 (v34) GTF file was subsampled to chromosome 1 (protein coding only) and mock datasets were made for each phenotype in the experimental design:

1. `control`: chromosome 1 arm 2
2. `lung`: chromosome 1 arm 1
3. `melanoma`: chromosome 1 arm 1 + 2

3 replicates of each phenotype were generated, with 2 of the 3 replicates being identical to ensure differentially expressed circRNAs can be detected by `DESeq2`.
