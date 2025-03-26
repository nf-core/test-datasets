This folder holds test datasets for modules (Vireo, Popscle, ScSplit, and Souporcell) that deconvolve single-cell multiplexed datasets. 
The dataset, provided by Popscle (https://github.com/statgen/popscle/), contains genetic data from 2 donors.
For compatibility and testing within nf-core pipelines, the BAM and VCF files have been downsampled to include only data from chromosome 21. `barcode.tsv` file was obtained by `samtools view chr21.bam | awk '{for(i=12;i<=NF;i++) if($i ~ /^CB:Z:/) print substr($i,6)}' | sort | uniq > barcodes.tsv`
