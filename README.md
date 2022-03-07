<<<<<<< HEAD
# test-datasets: `circdna`

This branch contains test data to be used for automated testing with the [nf-core/circdna](https://github.com/nf-core/circdna) pipeline.
=======
# test-datasets: `atacseq`

This branch contains test data to be used for automated testing with the [nf-core/atacseq](https://github.com/nf-core/atacseq) pipeline.
>>>>>>> 0c58a9f36205cc5f8c6bbb2ca03c401c61cb849d

## Content of this repository

`design.csv`: Experiment design file for minimal test dataset  
<<<<<<< HEAD
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

=======
`design_full.csv`: Experiment design file for full test dataset  
`reference/`: Genome reference files (iGenomes R64-1-1 Ensembl release)   
`testdata/` : FastQ files sub-sampled to 100,000 paired-end reads   

## Minimal test dataset origin

*S. cerevisiae* paired-end ATAC-seq dataset was obtained from:

Schep AN, Buenrostro JD, Denny SK, Schwartz K, Sherlock G, Greenleaf WJ. Structured nucleosome fingerprints enable high-resolution mapping of chromatin architecture within regulatory regions. Genome Res 2015 Nov;25(11):1757-70. [Pubmed](https://www.ncbi.nlm.nih.gov/pubmed/26314830) [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE66386)

### Sample information

| GEO_ID	    | SRA_ID	    | SAMPLE_NAME	                  |
|-------------|-------------|-------------------------------|
| GSM1621339	| SRR1822153	| Osmotic Stress Time 0 A rep1	|
| GSM1621340	| SRR1822154	| Osmotic Stress Time 0 A rep2	|
| GSM1621343	| SRR1822157	| Osmotic Stress Time 15 C rep1	|
| GSM1621344	| SRR1822158	| Osmotic Stress Time 15 C rep2	|

### Sampling procedure

The example command below was used to sub-sample the raw paired-end FastQ files to 100,000 reads (see [seqtk](https://github.com/lh3/seqtk)).

```bash
mkdir -p sample
seqtk sample -s100 SRR1822153_1.fastq.gz 100000 | gzip > ./sample/SRR1822153_1.fastq.gz
seqtk sample -s100 SRR1822153_2.fastq.gz 100000 | gzip > ./sample/SRR1822153_2.fastq.gz
```

### Expected output

To track and test the reproducibility of the pipeline with default parameters below are some of the expected outputs.

### Number of `mergedLibrary` broadPeaks

| sample	              | peaks	|
|-----------------------|-------|
| OSMOTIC_STRESS_T0_R1	| 890	  |
| OSMOTIC_STRESS_T0_R2	| 627	  |
| OSMOTIC_STRESS_T15_R1	| 1133	|
| OSMOTIC_STRESS_T15_R2	| 1135  |

### Number of `mergedReplicate` broadPeaks

| sample	              | peaks	|
|-----------------------|-------|
| OSMOTIC_STRESS_T0	    | 1130	|
| OSMOTIC_STRESS_T15	  | 1395	|

These are just guidelines and will change with the use of different software, and with any restructuring of the pipeline away from the current defaults.

## Full test dataset origin

*H. sapiens* paired-end ATAC-seq dataset was obtained from:

M Ryan Corces *et al*. An Improved ATAC-seq Protocol Reduces Background and Enables Interrogation of Frozen Tissues. Nat Methods. 2017 Oct;14(10):959-962. doi: 10.1038/nmeth.4396.
[Pubmed](https://pubmed.ncbi.nlm.nih.gov/28846090/) [SRA](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA380283)


### Sample information

| SRA ID     | SAMPLE NAME                                                         |
|------------|---------------------------------------------------------------------|
| SRR5427884 | ATAC-seq of in vitro culture GM12878 using the Standard ATAC method |
| SRR5427885 | ATAC-seq of in vitro culture GM12878 using the Standard ATAC method |
| SRR5427886 | ATAC-seq of in vitro culture GM12878 using the Omni-ATAC method	   |
| SRR5427887 | ATAC-seq of in vitro culture GM12878 using the Omni-ATAC method	   |
| SRR5427888 | ATAC-seq of in vitro culture GM12878 using the Fast-ATAC method	   |
| SRR5427889 | ATAC-seq of in vitro culture GM12878 using the Fast-ATAC method	   |
>>>>>>> 0c58a9f36205cc5f8c6bbb2ca03c401c61cb849d
