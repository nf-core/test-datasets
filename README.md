# test-datasets: `circdna`

This branch contains test data to be used for automated testing with the [nf-core/circdna](https://github.com/nf-core/circdna) pipeline.

## Content of this repository

`design.csv`: Experiment design file for minimal test dataset  
`reference/`: Homo Sapiens Assembly 38 (GRCh38) - Chr20:1Mbp-21Mbp (20Mb subset of Chr20 of GRCh38)   
`testdata/` : 220,000 FastQ paired-end reads

## Minimal test dataset origin
The data set was generated using Circle-Map Simulate (see [Circle-Map](https://github.com/iprada/Circle-Map) and InSilicoSeq (see [InSilicoSeq](https://github.com/HadrienG/InSilicoSeq). Circle-Map simulated 200,000 paired-end reads originated from circle-seq data and InSilicoSeq simulated random sequencing data from the reference genome.

### Data Generation

The example below was used to generate the raw paired-end FastQ files.

``` bash
Circle-Map Simulate -c 25 -g reference/Homo_sapiens_assembly38_chr20_20Mb.fasta -N 200000 -o simulated -r 150
iss generate --genomes ../../reference/Homo_sapiens_assembly38_chr20_20Mb.fasta -o iss_simulated --seed 1 --n_reads 20000 --model hiseq
cat iss_simulated_R1.fastq simulated_1.fastq > circdna_R1.fastq
cat iss_simulated_R2.fastq simulated_2.fastq > circdna_R2.fastq
```

### Expected output

To track and test the reproducibility of the pipeline with default parameters below are some of the expected outputs.

### Number of `Circle-Map Realign` circles

| sample	              | circles	|
|-----------------------|-------|
| Test	| 9953	  |

### Number of `Circexplorer2` circles

| sample	              | circles	|
|-----------------------|-------|
| TEST	    | 10257	|

### Number of `circle_finder` circles

| sample	              | circles	|
|-----------------------|-------|
| TEST	    | 8401	|

These are just guidelines and will change with the use of different software, and with any restructuring of the pipeline away from the current defaults.

