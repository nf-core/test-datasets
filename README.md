# test-datasets: `atacseq`

This branch contains test data to be used for automated testing with the [nf-core/atacseq](https://github.com/nf-core/atacseq) pipeline.

## Content of this repository

`design.csv`: Experiment design file
`reference/`: Genome reference files (iGenomes R64-1-1 Ensembl release)
`testdata/` : FastQ files randomly sub-sampled to 100,000 paired-end reads mapping to chromosome I

## Dataset origin

*S. cerevisiae* paired-end ATAC-seq dataset was obtained from:

Schep AN, Buenrostro JD, Denny SK, Schwartz K, Sherlock G, Greenleaf WJ. Structured nucleosome fingerprints enable high-resolution mapping of chromatin architecture within regulatory regions. Genome Res 2015 Nov;25(11):1757-70.

https://www.ncbi.nlm.nih.gov/pubmed/26314830
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE66386

### Sample information

| GEO_ID	    | SRA_ID	    | SAMPLE_NAME	                  | GROUP               |
|-------------|-------------|-------------------------------|---------------------|
| GSM1621339	| SRR1822153	| Osmotic Stress Time 0 A rep1	| OSMOTIC_STRESS_A_T0 |
| GSM1621340	| SRR1822154	| Osmotic Stress Time 0 A rep2	| OSMOTIC_STRESS_A_T0 |
| GSM1621341	| SRR1822155	| Osmotic Stress Time 0 B rep1	| OSMOTIC_STRESS_B_T0 |
| GSM1621342	| SRR1822156	| Osmotic Stress Time 0 B rep2	| OSMOTIC_STRESS_B_T0 |
| GSM1621343	| SRR1822157	| Osmotic Stress Time 15 C rep1	| OSMOTIC_STRESS_T15  |
| GSM1621344	| SRR1822158	| Osmotic Stress Time 15 C rep2	| OSMOTIC_STRESS_T15  |
| GSM1621349	| SRR1822163	| Osmotic Stress Time 60 F rep1	| OSMOTIC_STRESS_T60  |
| GSM1621350	| SRR1822164	| Osmotic Stress Time 60 F rep2	| OSMOTIC_STRESS_T60  |

## Expected output
