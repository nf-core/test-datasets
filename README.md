# test-datasets: `circdna`

This branch contains test data to be used for automated testing with the [nf-core/circdna](https://github.com/nf-core/circdna) pipeline.

## Content of this repository

`reference/`: Genome reference files (iGenomes R64-1-1 Ensembl release)

`testdata/` : 200,000 FastQ paired-end reads

## Minimal test dataset origin
The data set was generated using Circle-Map Simulate (see [Circle-Map](https://github.com/iprada/Circle-Map) and InSilicoSeq (see [InSilicoSeq](https://github.com/HadrienG/InSilicoSeq). Circle-Map simulated 120,000 paired-end reads originated from circle-seq data and InSilicoSeq simulated 80,000 random reads from the reference genome.

### Data Generation

The example below was used to generate the raw paired-end FastQ files.

``` bash
Circle-Map Simulate -c 200 -g genome.fa -N 120000 -r 150 -b cm_1 -p 10
Circle-Map Simulate -c 200 -g genome.fa -N 120000 -r 150 -b cm_2 -p 10
Circle-Map Simulate -c 200 -g genome.fa -N 120000 -r 150 -b cm_3 -p 10
wgsim -1 150 -2 150 -N 80000 genome.fa wgsim_1_R1.fastq wgsim_1_R2.fastq -S 1
wgsim -1 150 -2 150 -N 80000 genome.fa wgsim_2_R1.fastq wgsim_2_R2.fastq -S 1
wgsim -1 150 -2 150 -N 80000 genome.fa wgsim_3_R1.fastq wgsim_3_R2.fastq -S 1
cat cm_1_2.fastq wgsim_1_R2.fastq | gzip --no-name > ../testdata/circdna_1_R2.fastq.gz
cat cm_2_2.fastq wgsim_2_R2.fastq | gzip --no-name > ../testdata/circdna_2_R2.fastq.gz
cat cm_3_2.fastq wgsim_3_R2.fastq | gzip --no-name > ../testdata/circdna_3_R2.fastq.gz

cat cm_1_1.fastq wgsim_1_R1.fastq | gzip --no-name > ../testdata/circdna_1_R1.fastq.gz
cat cm_2_1.fastq wgsim_2_R1.fastq | gzip --no-name > ../testdata/circdna_2_R1.fastq.gz
cat cm_3_1.fastq wgsim_3_R1.fastq | gzip --no-name > ../testdata/circdna_3_R1.fastq.gz
```

### Expected output

To track and test the reproducibility of the pipeline with default parameters below are some of the expected outputs.

### Number of `Circle-Map Realign` circles

| sample	              | circles	|
|-----------------------|-------|
| circdna_1	| 926	  |
| circdna_2	| 924	  |
| circdna_3	| 931	  |

### Number of `Circexplorer2` circles

| sample	              | circles	|
|-----------------------|-------|
| circdna_1	| 1042	  |
| circdna_2	| 1070	  |
| circdna_3	| 1059	  |

### Number of `circle_finder` circles

| sample	              | circles	|
|-----------------------|-------|
| circdna_1	| 848	  |
| circdna_2	| 841	  |
| circdna_3	| 848	  |

### Number of `unicycler` circles

| sample	              | circles	|
|-----------------------|-------|
| circdna_1	| 1	  |
| circdna_2	| 0	  |
| circdna_3	| 0	  |

These are just guidelines and will change with the use of different software, and with any restructuring of the pipeline away from the current defaults.

