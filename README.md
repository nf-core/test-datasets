# test-datasets: `mnaseseq`

This branch contains test data to be used for automated testing with the [nf-core/mnaseseq](https://github.com/nf-core/mnaseseq) pipeline.

## Content of this repository  

`design.csv`: Sample design file  
`fastq/`: Raw FastQ files required to test the pipeline. Data originated from *S. cerevisiae*.

## Dataset origin

Raw fastq files were downloaded from [GSE117881](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117881) and sub-sampled to 100,000 read pairs per sample.

| GEO ID     | Sample description            |
|------------|-------------------------------|
| GSM3314439 | WT biological replicate 1     |
| GSM3314440 | WT biological replicate 2     |
| GSM3314441 | DEGRON biological replicate 1 |
| GSM3314442 | DEGRON biological replicate 2 |

**NOTE**: Genome reference files required to run the pipeline will be obtained from the existing [nf-core/chipseq test-dataset](https://github.com/nf-core/test-datasets/tree/chipseq).
