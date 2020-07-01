# test-datasets: `chipseq`
Test data to be used for automated testing with the nf-core pipelines

This branch contains test data for the [nf-core/chipseq](https://github.com/nf-core/chipseq) pipeline.

## Content of this repository

`design.csv`: Experiment design file for minimal test dataset  
`design_full.csv`: Experiment design file for full test dataset  
`reference/`: Genome reference files  
`testdata/` : Sub-sampled FastQ files sub-sampled

## Full test dataset origin

*H. sapiens* single-end ChIP-seq dataset was obtained from 2 separate studies as described [here](https://academic.oup.com/bib/article/17/6/953/2453197#47712047).

### Transcription factor data

[Pubmed](https://pubmed.ncbi.nlm.nih.gov/25752574/), [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE59530)

#### Sample information

| SRA ID     | GEO ID     | SAMPLE NAME                 |
|------------|------------|-----------------------------|
| SRR1635435 | GSM1534712 | ChIP-seq_Input_Vehicle_rep1 |
| SRR1635436 | GSM1534713 | ChIP-seq_Input_Vehicle_rep2 |
| SRR1635437 | GSM1534714 | ChIP-seq_Input_E2_rep1      |
| SRR1635438 | GSM1534715 | ChIP-seq_Input_E2_rep2      |
| SRR1635459 | GSM1534736 | ChIP-seq_FoxA1_Vehicle_rep1 |
| SRR1635460 | GSM1534737 | ChIP-seq_FoxA1_Vehicle_rep2 |
| SRR1635461 | GSM1534738 | ChIP-seq_FoxA1_E2_rep1      |
| SRR1635462 | GSM1534739 | ChIP-seq_FoxA1_E2_rep2      |

### Broad histone data

[Pubmed](https://pubmed.ncbi.nlm.nih.gov/25188243/), [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE57632)

#### Sample information

| SRA ID     | GEO ID     | SAMPLE NAME          |
|------------|------------|----------------------|
| SRR1285070 | GSM1385748 | NTKO_EZH2_ChIPseq_1  |
| SRR1285071 | GSM1385749 | NTKO_EZH2_ChIPseq_2  |
| SRR1285072 | GSM1385750 | TKO_EZH2_ChIPseq_1   |
| SRR1285073 | GSM1385751 | TKO_EZH2_ChIPseq_2   |
| SRR1285074 | GSM1385752 | NTKO_Input_ChIPseq_1 |
| SRR1285075 | GSM1385753 | NTKO_Input_ChIPseq_2 |
| SRR1285076 | GSM1385754 | TKO_Input_ChIPseq_1  |
| SRR1285077 | GSM1385755 | TKO_Input_ChIPseq_2  |
