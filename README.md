# test-datasets: `bcellmagic`

Test data to be used for automated testing with the nf-core pipeline.

This branch contains test data for the [nf-core/bcellmagic](https://github.com/nf-core/bcellmagic) pipeline.

## Content

### Metadata file and primer fasta files

`Metadata_test.tsv` contains the metadata sheet needed to run the pipeline.
`C_primers.fasta` and `V_primers.fasta` contain fake primers, do not use them!

### Paired-end test data

For each test Sample 1-8
`*R1.fastq.gz` and `*R2.fastq.gz`: paired end MiSeq data (~10K reads).
`*I1.fastq.gz` : containing illumina barcodes and UMI codes.
