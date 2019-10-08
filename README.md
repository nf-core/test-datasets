# test-datasets: `kmermaid`
Test data to be used for automated testing with the nf-core pipelines

This branch contains test data for the [nf-core/kmermaid](https://github.com/nf-core/kmermaid) pipeline.

## View `.bam` file obtained from 10x datasets

The `.bam` file can be viewed using samtools, and is processed in this pipeline using `pysam` in `sourmash`:

```samtools view 10x-example/possorted_genome_bam.bam | less -S```

## `.bai` file obtained from 10x datasets

	A bai file isn't an indexed form of a bam - it's a companion to your bam that contains the index.

## Barcodes are cell barcodes in `.tsv` file

Cell barcodes provided in `barcodes.tsv` could be unfiltered for "good barcodes" (high quality, not repeated, etc) so it could be the whole set of them which for the 10x V2 chemistry is 700k and for the 10x V3 chemistry is 2 million possible barcodes.
