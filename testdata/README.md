# Riboseq test data

Data in this folder are down-sampled read files from GSE182201, which is a mixture of RNA-seq and Ribo-seq data. Down-sampling has been done by taking alignment files, subsetting to chromosome 20, and cross-referencing the aligned reads with raw FASTQ files. See [make_test_data.sh](make_test_data.sh) for details. 
