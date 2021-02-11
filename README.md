# ![nfcore/test-datasets](docs/images/test-datasets_logo.png)
Test data to be used for automated testing with the nf-core pipelines

## nf-core/circrna
This branch contains test-data for `nf-core/circrna`.

### Contents of branch:
* `fastq/` 6 FASTQ read pairs.
* `reference/` Reference annotation files
* `source/` Scripts used to simulate circRNA FASTQ reads.
* `phenotype.txt` metadata file for DESeq2.
* `samples.csv` input test-dataset csv file

### Test-dataset generation strategy:
Gencode GRCh38 (v34) GTF file was sub-sampled to chromosome 1, and split into 2 GTF files per chromosome arm, including only protein coding regions.

`CIRI_simulator.pl` was used on each sub-sampled chromosome 1 arm to generate reads mimicking differential circRNA expression.
