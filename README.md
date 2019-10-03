# test-datasets: `kmermaid`
Test data to be used for automated testing with the nf-core pipelines

This branch contains test data for the [nf-core/kmermaid](https://github.com/nf-core/kmermaid) pipeline.

## View .bam file obtained from 10x datasets

	A .bam file can be viewed using samtools, and is processed using pysam in sourmash (wrapper around samtools in kmermaid):

	```samtools view 10x-example/possorted_genome_bam.bam | less -S```

## Barcodes are cell barcodes in .tsv file

	 	Cell barcodes provided in barcodes.tsv could be unfiltered for "good barcodes" (high quality, not repeated, etc) so it could be the whole set of them which for the 10x v2 chemistry is 700k and for the 10x v3 chemistry is 2 million possible barcodes

## Renamer barcodes in .tsv file
	
	For bam/10x files, Use this to specify the location of your tsv (tab separated file) containing map of cell barcodes and their corresponding new names(e.g row in the tsv file: AAATGCCCAAACTGCT-1    lung_epithelial_cell|AAATGCCCAAACTGCT-1). 
