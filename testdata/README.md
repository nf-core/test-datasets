# Riboseq test data

FASTQ data in this folder are down-sampled read files from GSE182201, which is a mixture of RNA-seq and Ribo-seq data. Down-sampling has been done by taking alignment files, subsetting to chromosome 20, and cross-referencing the aligned reads with raw FASTQ files. See [make_test_data.sh](make_test_data.sh) for details.

hg38-mature-tRNAs-dna.fasta is human tRNAs derived from [GtRNAdb](http://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi38/Hsapi38-seq.html) in March 2023. 
