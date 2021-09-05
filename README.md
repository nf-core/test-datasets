# test-datasets: `airrflow`

Test data to be used for automated testing with the nf-core pipeline.

This branch contains test data for the [nf-core/airrflow](https://github.com/nf-core/airrflow) pipeline.

## Content

### BCR test data

`testdata-bcr` contains the test data needed to run the pipeline for BCRseq data.

* `Metadata_test.tsv` contains the paths to the BCR test data.
* `C_primers.fasta` and `V_primers.fasta` contain fake primers, do not use them!
* The `fastq` files are random samples of a human BCRseq MiSeq amplicon sequencing data, for 4 different B-cell populations and two time points. The `*I1.fastq.gz` files are the index reads with illumina and UMI barcodes.

`testdata-no-umi` contains BCR test data to run the pipeline without UMIs. This dataset is a minimal example to test pipeline function based on the Illumina MiSeq 2x250 BCR mRNA workflow from the presto docs.

* `metadata_test-no-umi.tsv` contains the paths to the associated test data.
* `Greiff2014_Cprimers.fasta` and `Greiff2014_Vprimers.fasta` contain the C and
  V primers, respectively
* The `fastq` files are 2k read subsamples of a human BCRseq MiSeq amplicon
  sequencing data. The I1 / UMI barcode field in the metadata file is left empty
  in this case.

### TCR test data

`testdata-tcr` contains the test data needed to run the pipeline for TCRseq data.

* `TCR_metadata.tsv` contains the paths to the TCR test data.
* `C_primers.fasta` and `linker.fasta` are the primer sequences, and linker sequences for the 5' RACE amplification of the TCR.
* The `fastq` files are random samples of a human TCRseq 5'RACE sequencing data produced with the TAKARA protocol, for 2 different samples.

### Reveal test data

`testdata-reveal` contains test data needed to run the AIRR-seq workflow Reveal. Further details can be found in the README.md file inside the folder.
