# test-datasets `porepatrol`
Test data to be used for automated testing with the [nf-core/porepatrol](https://github.com/nf-core/porepatrol) pipeline.

## Content of this repository

`nanopore_reads.fastq.gz`: Basecalled ONT reads, gzipped. 

## Dataset origin

This dataset was obtained from:

https://melbourne.figshare.com/articles/Basecalled_ONT_reads/5170843, by Ryan Wick.

The sample is from a bacteria, sequenced with MinION, and basecalled with Albacore. 

I downloaded file `barcode04.fastq.gz`

I subsampled to 0.5% with seqtk:

```bash
gunzip barcode04.fastq.gz
seqtk sample barcode04.fastq 0.005 > nanopore_reads.fastq
gzip nanopore_reads.fastq
```


## Expected output

Running the pipeline locally with the default parameters:

```bash
nextflow run porepatrol --reads nanopore_reads.fastq.gz 
```

Expected file sizes:

* Gzipped input file size: 634 KB
* Unzipped file size: 1,264 KB
* After adapters are chopped: 1,259 KB
* After filtering with defaults (-q 12): 935 KB



