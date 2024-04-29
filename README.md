# test-datasets: `tfactivity`

This branch contains test data to be used for automated testing with the [nf-core/tfactivity](https://github.com/nf-core/tfactivity) pipeline.

## Content of this repository

`reference/`: Sub-sampled genome-specific files. Based on mouse genome mm10, but only chromosome 1 is included.

`peaks`: Peak files for testing the pipeline. Based on HM ChIP-seq data, but only chromosome 1 is included.
`rna-seq`: RNA-seq count files for testing the pipeline.
`bams`: BAM files that are only used in the `test_full` pipeline test.

`samplesheet`: Sample sheets for each of `peaks`, `rna-seq`, and `bams` folders.

## Test dataset origin

*M. musculus* dataset was obtained from:

> Lee HK, Willi M, Kuhns T, Liu C, Hennighausen L. Redundant and non-redundant cytokine-activated enhancers control Csn1s2b expression in the lactating mouse mammary gland. Nat Commun. 2021 Apr 14;12(1):2239. doi: 10.1038/s41467-021-22500-w. [Pubmed](https://pubmed.ncbi.nlm.nih.gov/33854063/) [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE161620)

### Test dataset pre-processing
The RNA-Seq data was processed using the `nf-core/rnaseq` pipeline. 
The ChIP-Seq data was processed using the `nf-core/chipseq` pipeline.