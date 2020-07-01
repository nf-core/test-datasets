# test-datasets: `atacseq`

This branch contains test data to be used for automated testing with the [nf-core/atacseq](https://github.com/nf-core/atacseq) pipeline.

## Content of this repository

`design.csv`: Experiment design file for minimal test dataset  
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
