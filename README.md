# test-datasets: `scdownstream`

This branch contains test data to be used for automated testing with the [nf-core/scdownstream](https://github.com/nf-core/scdownstream) pipeline.

## Content of this repository

`samples/`: The pipeline input in `h5ad` and `rds` formats.
`samplesheet.csv`: The samplesheet used as the pipeline input.

## Test dataset origin

The test dataset is the output of the test configuration of the [nf-core/scrnaseq](https://github.com/nf-core/scrnaseq) pipeline. For more information on the test dataset, please refer to the [documentation of the nf-core/scrnaseq test data](https://github.com/nf-core/test-datasets/tree/scrnaseq).

### Test dataset pre-processing
The scRNA-Seq data was processed using the `nf-core/scrnaseq` pipeline. 
