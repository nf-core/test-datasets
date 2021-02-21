# test-datasets: `bcellmagic`

Test data to be used for automated testing with the nf-core pipeline.

This branch contains test data for the [nf-core/bcellmagic](https://github.com/nf-core/bcellmagic) pipeline.

## Content

### BCR test data

`testdata-bcr` contains the test data needed to run the pipeline for BCRseq data.

* `metadata.tsv` contains the paths to the BCR test data.
* `C_primers.fasta` and `V_primers.fasta` contain fake primers, do not use them!
* The `fastq` files are random samples of a human BCRseq MiSeq amplicon sequencing data, for 4 different B-cell populations and two time points. The `*I1.fastq.gz` files are the index reads with illumina and UMI barcodes.

### TCR test data

`testdata-tcr` contains the test data needed to run the pipeline for BCRseq data.

* `metadata.tsv` contains the paths to the TCR test data.
* `C_primers.fasta` and `linker.fasta` are the primer sequences, and linker sequences for the 5' RACE amplification of the TCR.
* The `fastq` files are random samples of a human TCRseq 5'RACE sequencing data produced with the TAKARA protocol, for 2 different samples.
