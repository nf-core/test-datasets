# test-datasets: `slamseq`

This branch contains test data to be used for automated testing with the [nf-core/slamseq](https://github.com/nf-core/slamseq) pipeline.

## Content of this repository

`reference/`: Genome reference files (*Homo sapiens* hg38 + RefSeq 3'UTR subset to chr8)   
`testdata/` : FastQ files with SR100 reads from chr8

## Dataset origin

*H. sapiens* single-end 100 bp SLAM-seq dataset was obtained from:

Muhar M, Ebert A, Neumann T, Umkehrer C, Jude J, Wieshofer C, Rescheneder P, Lipp JJ, Herzog VA, Reichholf B, Cisneros DA, Hoffmann T, Schlapansky MF, Bhat P, von Haeseler A, Köcher T, Obenauf AC, Popow J, Ameres SL & Zuber J: SLAM-seq defines direct gene-regulatory functions of the BRD4-MYC axis. Science, 2018. Vol. 360, Issue 6390, pp. 800-805. [Pubmed](https://www.ncbi.nlm.nih.gov/pubmed/29622725) [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111463)

### Sample information

| GEO_ID	    | SRA_ID	    | SAMPLE_NAME	      |
|-------------|-------------|-------------------|
| GSM2691882	| SRR5806795	| MOLM-13_dmso_1	  |
| GSM2691883	| SRR5806796	| MOLM-13_dmso_2	  |
| GSM2691884	| SRR5806797	| MOLM-13_dmso_3	  |
| GSM2691891	| SRR5806804	| MOLM-13_nvp.hi_1	|
| GSM2691892	| SRR5806805	| MOLM-13_nvp.hi_2	|
| GSM2691893	| SRR5806806	| MOLM-13_nvp.hi_3	|

## Expected output

To track and test the reproducibility of the pipeline with default parameters below are some of the expected outputs.
