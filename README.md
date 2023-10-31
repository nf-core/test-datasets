# test-datasets: `callingcards`

This branch contains test data to be used for automated testing with the [nf-core/callingcards](https://github.com/nf-core/callingcards) pipeline.

## Content of this repository

### Mammals data

`mammals/`: Subsampled test data and other resources meant to test the mammal specific functionality

`mammals/*.fastq.gz`: Significantly downsampled raw hops reads
`mammals/barcode_details.json`: A configuration file describing the read barcode structure suitable for human and mouse data
`mammals/chr1.fa`: A short subset of GRCh38 chr1
`mammals/chr1.gtf`: A subset of the gencode v38 human chr1 gtf format annotations
`mammals/samplesheet.csv`: A samplesheet for the minimal mammals test data
`mammals/samplesheet_full.csv`: A samplesheet for the full mammals test data

### Yeast data

`yeast/`: Subsampled test data and other resources meant to test the yeast specific functionality

`yeast/*.fastq.gz`: Significantly downsampled raw hops reads
`yeast/run_6177_barcode_details.json`: A configuration file describing the read barcode structure suitable for yeast data
`yeast/samplesheet.csv`: input samplesheet for the yeast test data


`mammals/samplesheet.csv`: Experiment design file for minimal test dataset

## Minimal test dataset origin

- mammals: a significant subsampling of unpublished data meant for process testing only
- yeast: a significant subsampling of unpublished data meant for process testing only

### Sampling information

### Sampling procedure

