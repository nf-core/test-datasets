# test-datasets: `airrflow`

Test data to be used for automated testing with the nf-core pipeline.

This branch contains test data for the [nf-core/airrflow](https://github.com/nf-core/airrflow) pipeline.

## Content

### BCR test data

`testdata-bcr` contains the test data needed to run the pipeline for BCRseq data.

- `Metadata_test.tsv` samplesheet for BCR test data for nf-core/airrflow version < 3.0.
- `Metadata_test_airr.tsv` samplesheet for BCR test data for nf-core/airrflow version >= 3.0.
- `C_primers.fasta` and `V_primers.fasta` contain fake primers, do not use them!
- The `fastq` files are random samples of a human BCRseq MiSeq amplicon sequencing data, for 4 different B-cell populations and two time points. The `*I1.fastq.gz` files are the index reads with illumina and UMI barcodes.

- `metadata_pcr_umi_airr.tsv` contains the samplesheet for the pipeline full tests for nf-core/airrflow version < 3.0.
- `metadata_pcr_umi_airr_300.tsv` contains the samplesheet for the pipeline full tests for nf-core/airrflow version >3.0.

`testdata-no-umi` contains BCR test data to run the pipeline without UMIs. This dataset is a minimal example to test pipeline function based on the Illumina MiSeq 2x250 BCR mRNA workflow from the presto docs.

- `metadata_test-no-umi.tsv` contains the paths to the associated test data for nf-core/airrflow < 3.0.
- `metadata_test-no-umi_airr.tsv` conatins the paths to the associated test data for nf-core/airrflow > 3.0.
- `Greiff2014_Cprimers.fasta` and `Greiff2014_Vprimers.fasta` contain the C and
  V primers, respectively
- The `fastq` files are 2k read subsamples of a human BCRseq MiSeq amplicon
  sequencing data. The I1 / UMI barcode field in the metadata file is left empty
  in this case.

`testdata-clontech` contains the test data needed to run the pipeline on data generated with the clontech-bcr-umi protocol. The first 25K reads of the SRR23055142 sample were downloaded.
`testdata-neb` contains the test data needed to run the pipeline on data generated with the neb-bcr-umi protocol. The first 25K reads of the SRR4026043 sample were downloaded.

### TCR test data

`testdata-tcr` contains the test data needed to run the pipeline for TCRseq data.

- `TCR_metadata.tsv` contains the paths to the TCR test data for the nf-core/airrflow < 3.0.
- `TCR_metadata_airr.tsv` contains the paths to the TCR test data for nf-core/airrflow > 3.0.
- `C_primers.fasta` and `linker.fasta` are the primer sequences, and linker sequences for the 5' RACE amplification of the TCR.
- The `fastq` files are random samples of a human TCRseq 5'RACE sequencing data produced with the TAKARA protocol, for 2 different samples.

### single-cell test data

`testdata-sc` contains the test data required to run the pipeline for 10xGenomics derived single-cell AIRR-seq data, currently only test data wit TCR sequences are available.

- `10x_sc_raw.tsv` contains the paths to the sc test data.
- `refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0.tar.gz` is the 10xGenomics reference dataset.
- The `fastq` files are subsampled samples from a healthy donor made publically available by 10xGenomics ([dataset](https://www.10xgenomics.com/datasets/human-t-cells-from-a-healthy-donor-1-k-cells-multi-v-2-2-standard-5-0-0)).

### RNA-seq test data

`testdata-rnaseq` contains the test data required to run the pipeline from raw RNA-seq data, both for bulk and single-cell derived FASTQ files.

- `rnaseq_metadata.tsv` contains the paths to the RNA-seq test data
- `IMGT+C.fa` is the reference file downloaded from the IMGT website using TRUST4's `BuildImgtAnnot.pl` script ([docs](https://github.com/liulab-dfci/TRUST4?tab=readme-ov-file#build-custom-vjc-gene-database-files-for--f-and---ref))
- The bulk `fastq` files are derived from [TRUST4's example data](https://github.com/liulab-dfci/TRUST4/tree/master/example)

## Database cache

`database cache` contains the IMGT and Igblast db human and mouse database caches for running the pipeline tests.

### Reveal test data

`testdata-reveal` contains test data needed to run the AIRR-seq workflow Reveal. Further details can be found in the README.md file inside the folder.
