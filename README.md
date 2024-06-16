# test-datasets: `scdownstream`

This branch contains test data to be used for automated testing with the [nf-core/scdownstream](https://github.com/nf-core/scdownstream) pipeline.

## Content of this repository

`samples/`: The pipeline input in `h5ad` and `rds` formats.
`samplesheet.csv`: The samplesheet used as the pipeline input.

## Test dataset origin

The data used in this test dataset is a subset of the data from the study by [He et al. (2020)](https://doi.org/10.1016/j.jaci.2020.01.042). The data was downloaded from the [Gene Expression Omnibus (GEO)](https://www.ncbi.nlm.nih.gov/geo/) with the accession number [GSE147424](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE147424).

### Test dataset pre-processing
The scRNA-Seq data was processed using the `nf-core/scrnaseq` pipeline. 
